# ==============================================================================
# 03_signals_mono.R - monopharmacy analysis: SSRI reports with no other drug
# entry, against a reference set restricted to monopharmacy reports.
#
# Same steps as 02 (PT signals, group composites, H-S composite, case lists)
# with both the group and the reference set restricted to monopharmacy.
#
# Main outputs (Supplementary "DA", right-hand block):
#   ssri_mono_targeted_signals.csv, ssri_mono_perdrug_pt_signals.csv,
#   varenicline_mono_targeted_signals.csv, group_pooled_mono.csv,
#   mono_homicide_suicide_composite_phv.csv, ssri_mono_perdrug_hs_composite.csv,
#   varenicline_mono_hs_composite_phv.csv, hs_composite_all_mono.csv,
#   mono_cases_homicide_suicide.csv
# ==============================================================================
source("R/load_all.R")
need_faers()
t0 <- Sys.time()
log_msg("=== 03 SIGNALS: MONOPHARMACY ===")

# ── SSRI monopharmacy IDs (from step 02) ──────────────────────────────────────
codrug   <- fread(require_res("ssri_codrug_per_report.csv"))[, primaryid := bit64::as.integer64(primaryid)]
# The same rule the reference set is built with (mono_reference_ids): zero
# co-drugs AND at least one drug named, because n_codrugs is floored at zero
# and so reads 0 for a report naming nothing. The two must agree, or the check
# at the end of this block stops the run.
mono_ids <- codrug[n_codrugs == 0 & n_distinct_drugs >= 1L, primaryid]
codrug_mono <- codrug[n_codrugs == 0 & n_distinct_drugs >= 1L]
n_ssri_all  <- nrow(codrug)

# ── monopharmacy reference set ────────────────────────────────────────────────
master <- apply_blacklist(readRDS(require_res("master_mono.rds")), "monopharmacy reference set")
n_master    <- nrow(faers_get(master, "demo"))
norm_master <- pt_norm_lookup(master)
full_ids    <- pt_ids_table(build_obj_reac_norm(master, norm_master))
group_bg    <- group_background(master, norm_master)
n_hs_master <- hs_master_count(full_ids)
log_msg("Monopharmacy reference set: ", sfmt(n_master), " | H-S co-reports: ", sfmt(n_hs_master))
use_ebgm_prior(master, norm_master, n_master, "mono")

# ── groups ───────────────────────────────────────────────────────────────────
ssri <- drug_filter(master, ssri_pattern)
ssri <- faers_keep(ssri, primaryid = id_intersect(faers_get(ssri, "demo")$primaryid, mono_ids))
var  <- apply_blacklist(drug_filter(master, varenicline_pattern), "varenicline monopharmacy group")
n_ssri <- nrow(faers_get(ssri, "demo"))
n_var  <- nrow(faers_get(var, "demo"))
# One rule decides monopharmacy, so the SSRI reports selected by the group
# filter and those in the reference set must be the same reports. A difference
# means the rule was applied inconsistently somewhere.
if (n_ssri != length(mono_ids))
  stop("SSRI monopharmacy group (", n_ssri, ") != monopharmacy IDs from step 02 (", length(mono_ids),
          "). Check the reference set and data/mono_overrides.csv.")
log_msg("SSRI monopharmacy: ", sfmt(n_ssri), " (", spct(n_ssri, n_ssri_all), " of SSRI group) | varenicline: ", sfmt(n_var))

# ── PT signals ────────────────────────────────────────────────────────────────
pt_signals <- function(obj, n) compute_pt_signals(obj, full_ids, n, n_master, norm_master)

tgt <- cached("ssri_mono_targeted_signals.rds", pt_signals(ssri, n_ssri))
fwrite(tgt[, !"ids"], res_path("ssri_mono_targeted_signals.csv"))

perdrug <- rbindlist(lapply(names(ssri_drugs), function(nm) {
  res <- cached(paste0("mono_perdrug_pt_", nm, ".rds"), {
    sub <- drug_filter(ssri, ssri_drug_pattern[[nm]])
    pt_signals(sub, nrow(faers_get(sub, "demo")))
  })
  res[, drug := nm]
}), fill = TRUE)
fwrite(perdrug[, !"ids"], res_path("ssri_mono_perdrug_pt_signals.csv"))
for (p in targets)   # each monopharmacy report names one SSRI, so these must match
  if (perdrug[pt == p, sum(a)] != tgt[pt == p, a]) warning("per-drug sum != pooled a for ", p)

var_tgt <- cached("varenicline_mono_targeted_signals.rds", pt_signals(var, n_var)[, drug := "varenicline"])
fwrite(var_tgt[, !"ids"], res_path("varenicline_mono_targeted_signals.csv"))

# ── Suicidality / Homicidality group composites ───────────────────────────────
grp <- cached("group_pooled_mono.rds", {
  sets <- c(list(list(obj = ssri, n = n_ssri, label = "SSRI (POOLED)")),
            lapply(names(ssri_drugs), function(nm) {
              sub <- drug_filter(ssri, ssri_drug_pattern[[nm]])
              list(obj = sub, n = nrow(faers_get(sub, "demo")), label = nm) }),
            list(list(obj = var, n = n_var, label = "Varenicline")))
  # fill = TRUE on both levels: a composite with fewer than three cases carries
  # no estimate, and that row has fewer columns than one that does.
  add_signals(rbindlist(lapply(sets, function(s) rbindlist(lapply(names(group_defs), function(g)
    compute_group_signal(s$obj, g, group_bg, s$n, n_master, norm_master, s$label)),
    fill = TRUE)), fill = TRUE))
})
fwrite(grp[, !"ids"], res_path("group_pooled_mono.csv"))

# ── homicide-suicide composite ────────────────────────────────────────────────
hs_ids <- function(sig) id_intersect(get_ids_for_pt(sig, "Homicide"), get_ids_for_pt(sig, "Completed Suicide"))

hs_perdrug <- rbindlist(lapply(names(ssri_drugs), function(nm) {
  res <- perdrug[drug == nm]
  n_d <- (as.numeric(res$a) + as.numeric(res$b))[1]
  hs_composite(length(hs_ids(res)), n_d, n_hs_master, n_master)[, `:=`(drug = nm, pt = HS_PT, n_group = n_d)]
}), fill = TRUE)
fwrite(hs_perdrug, res_path("ssri_mono_perdrug_hs_composite.csv"))

HS_PT_MONO <- "Homicide + Completed suicide (co-reported, mono)"
ssri_hs <- hs_ids(tgt)
# hs_composite already returns an empty-estimate row below the three-case floor,
# so write it either way: the pair must exist under both analyses.
hs_phv <- hs_composite(length(ssri_hs), n_ssri, n_hs_master, n_master)[, pt := HS_PT_MONO]
fwrite(hs_phv, res_path("mono_homicide_suicide_composite_phv.csv"))
hs_var <- hs_composite(length(hs_ids(var_tgt)), n_var, n_hs_master, n_master)[
  , `:=`(pt = HS_PT_MONO, drug_group = "varenicline_mono", n_group = n_var)]
fwrite(hs_var, res_path("varenicline_mono_hs_composite_phv.csv"))
write_res(rbindlist(list(copy(hs_phv)[, `:=`(drug_group = "ssri_mono", n_group = n_ssri)], hs_var),
                    fill = TRUE), "hs_composite_all_mono.csv")
if (sum(hs_perdrug$a) != hs_phv$a) warning("per-drug H-S sum != pooled H-S a")

# ── case list ─────────────────────────────────────────────────────────────────
# Only the homicide-suicide listing is written here. The per-PT monopharmacy
# counts are in the signal tables; the reports behind them are in the
# unrestricted listings, which carry the monopharmacy flag.
if (length(ssri_hs)) {
  log_msg("Case list...")
  write_res(build_case_list(ssri, ssri_hs, codrug_mono), "mono_cases_homicide_suicide.csv")
}
rm(master); gc2()

log_msg("=== 03 COMPLETE | ", elapsed(t0), " ===")
