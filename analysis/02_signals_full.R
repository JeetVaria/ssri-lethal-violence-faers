# ==============================================================================
# 02_signals_full.R - unrestricted analysis: SSRI reports (SSRI = primary
# suspect, any number of co-drugs) against the full FAERS reference set.
#
#   - SSRI group, monopharmacy/polypharmacy split, and the varenicline
#     comparator group (all after the manual blacklist)
#   - disproportionality for each target PT: pooled SSRIs, each SSRI, varenicline
#   - Suicidality / Homicidality group composites
#   - Homicide-suicide (H-S) composite: pooled, per SSRI, varenicline
#   - per-PT case lists (identical report IDs to the signal counts)
#   - SSRI group characteristics tables
#
# Main outputs (Supplementary "DA", "SUMMARY", "COUNTS"):
#   ssri_targeted_signals.csv, ssri_perdrug_pt_signals.csv,
#   group_pooled_full.csv, homicide_suicide_composite_phv.csv,
#   ssri_perdrug_hs_composite.csv, hs_composite_all_pooled.csv,
#   cases_<pt>.csv, cases_homicide_suicide.csv, ssri_codrug_per_report.csv,
#   ssri_codrug_distribution.csv, ssri_drug_breakdown.csv, ssri_by_year.csv,
# ==============================================================================
source("R/load_all.R")
need_faers()
t0 <- Sys.time()
log_msg("=== 02 SIGNALS: UNRESTRICTED ===")

# ── reference set ─────────────────────────────────────────────────────────────
master <- apply_blacklist(readRDS(require_res("master_dedup.rds")), "master")
n_master    <- nrow(faers_get(master, "demo"))
norm_master <- pt_norm_lookup(master, applied = "llt_mappings_applied.csv")
full_ids    <- pt_ids_table(build_obj_reac_norm(master, norm_master))
group_bg    <- group_background(master, norm_master)
n_hs_master <- hs_master_count(full_ids)
log_msg("Reference set: ", sfmt(n_master), " reports | H-S co-reports: ", sfmt(n_hs_master))
use_ebgm_prior(master, norm_master, n_master, "full")

# ── groups ───────────────────────────────────────────────────────────────────
ssri <- apply_blacklist(drug_filter(master, ssri_pattern), "SSRI group")
var  <- apply_blacklist(drug_filter(master, varenicline_pattern), "varenicline group")
saveRDS(ssri, res_path("ssri_group.rds"))
saveRDS(var,  res_path("varenicline_group.rds"))
n_ssri <- nrow(faers_get(ssri, "demo"))
n_var  <- nrow(faers_get(var, "demo"))

codrug  <- build_codrug_table(ssri)
fwrite(codrug, res_path("ssri_codrug_per_report.csv"))
# Monopharmacy is zero co-drugs AND at least one drug named. n_codrugs is
# floored at zero, so it also reads 0 for a report naming no drug at all - the
# group filter can match on prod_ai without any usable drugname. Same rule as
# mono_reference_ids(); step 03 stops if the two disagree.
codrug[, is_mono := n_codrugs == 0 & n_distinct_drugs >= 1L]
n_mono   <- codrug[(is_mono), .N]
n_nodrug <- codrug[n_distinct_drugs == 0L, .N]
if (n_nodrug)
  log_msg("  ", sfmt(n_nodrug), " SSRI report(s) name no drug: not counted as monopharmacy")
log_msg("SSRI: ", sfmt(n_ssri), " (mono ", sfmt(n_mono), ", poly ", sfmt(n_ssri - n_mono),
        ") | varenicline: ", sfmt(n_var))

# ── PT signals ────────────────────────────────────────────────────────────────
pt_signals <- function(obj, n) compute_pt_signals(obj, full_ids, n, n_master, norm_master)

tgt <- cached("ssri_targeted_signals.rds", pt_signals(ssri, n_ssri))
fwrite(tgt[, !"ids"], res_path("ssri_targeted_signals.csv"))

perdrug <- rbindlist(lapply(names(ssri_drugs), function(nm) {
  res <- cached(paste0("perdrug_pt_", nm, ".rds"), {
    sub <- drug_filter(ssri, ssri_drug_pattern[[nm]])
    pt_signals(sub, nrow(faers_get(sub, "demo")))
  })
  res[, drug := nm]
}), fill = TRUE)
fwrite(perdrug[, !"ids"], res_path("ssri_perdrug_pt_signals.csv"))

var_tgt <- cached("varenicline_targeted_signals.rds", pt_signals(var, n_var)[, drug := "varenicline"])
fwrite(var_tgt[, !"ids"], res_path("varenicline_targeted_signals.csv"))

# Check: each SSRI report has one SSRI as primary suspect in almost all cases, so
# per-drug counts should add up to (or slightly exceed) the pooled count.
for (p in targets) log_msg("  ", p, ": pooled a=", sfmt(tgt[pt == p, a]), " | sum per-drug a=", sfmt(perdrug[pt == p, sum(a)]))

# ── Suicidality / Homicidality group composites ───────────────────────────────
grp <- cached("group_pooled_full.rds", {
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
fwrite(grp[, !"ids"], res_path("group_pooled_full.csv"))

# ── homicide-suicide composite ────────────────────────────────────────────────
hs_ids <- function(sig) id_intersect(get_ids_for_pt(sig, "Homicide"), get_ids_for_pt(sig, "Completed Suicide"))

hs_perdrug <- rbindlist(lapply(names(ssri_drugs), function(nm) {
  res <- perdrug[drug == nm]
  n_d <- (as.numeric(res$a) + as.numeric(res$b))[1]
  hs_composite(length(hs_ids(res)), n_d, n_hs_master, n_master)[, `:=`(drug = nm, pt = HS_PT, n_group = n_d)]
}), fill = TRUE)
fwrite(hs_perdrug, res_path("ssri_perdrug_hs_composite.csv"))

ssri_hs <- hs_ids(tgt)
# hs_composite already returns an empty-estimate row below the three-case floor,
# so write it either way: the pair must exist under both analyses.
hs_phv <- hs_composite(length(ssri_hs), n_ssri, n_hs_master, n_master)[, pt := HS_PT]
fwrite(hs_phv, res_path("homicide_suicide_composite_phv.csv"))
hs_var <- hs_composite(length(hs_ids(var_tgt)), n_var, n_hs_master, n_master)[, `:=`(pt = HS_PT, drug_group = "varenicline", n_group = n_var)]
fwrite(hs_var, res_path("varenicline_hs_composite_phv.csv"))
write_res(rbindlist(list(copy(hs_phv)[, `:=`(drug_group = "ssri", n_group = n_ssri)], hs_var),
                    fill = TRUE), "hs_composite_all_pooled.csv")
log_msg("H-S composite (SSRI): a=", length(ssri_hs),
        if (!is.na(hs_phv$ror)) paste0(" ROR=", round(hs_phv$ror, 2), " [", round(hs_phv$ror_ci_low, 2), "-", round(hs_phv$ror_ci_high, 2), "]"))

# ── case lists ────────────────────────────────────────────────────────────────
log_msg("Case lists...")
write_case_lists(pt_files, list(
  list(name = "ssri",        obj = ssri, signals = tgt,     codrug = codrug),
  list(name = "varenicline", obj = var,  signals = var_tgt, codrug = build_codrug_table(var))))
if (length(ssri_hs)) write_res(build_case_list(ssri, ssri_hs, codrug), "cases_homicide_suicide.csv")
rm(master); gc2()

# ── SSRI group characteristics ────────────────────────────────────────────────────────
# Split on is_mono as well, so the 0 bucket cannot label a report that names no
# drug as monopharmacy. With no such report the table is exactly as before.
cd_dist <- codrug[, .N, by = .(n_codrugs, is_mono)][order(n_codrugs, -is_mono)]
cd_dist[, `:=`(pct = round(N / sum(N) * 100, 1), cumulative_pct = round(cumsum(N) / sum(N) * 100, 1),
               label = ifelse(is_mono, "0 (monopharmacy)",
                       ifelse(n_codrugs == 0, "0 (no drug named)",
                       ifelse(n_codrugs <= 5, as.character(n_codrugs),
                       ifelse(n_codrugs <= 10, "6-10", "11+")))))]
write_res(cd_dist[, .(n_codrugs, label, N, pct, cumulative_pct)], "ssri_codrug_distribution.csv")

# reports per SSRI (as primary suspect); a report naming two SSRIs counts for both
ps_drug <- faers_get(ssri, "drug")[toupper(role_cod) == "PS"]
breakdown <- rbindlist(lapply(names(ssri_drugs), function(nm) {
  pat <- ssri_drug_pattern[[nm]]
  hit <- !is.na(ps_drug$drugname) & grepl(pat, ps_drug$drugname, ignore.case = TRUE, perl = TRUE)
  if ("prod_ai" %in% names(ps_drug))
    hit <- hit | (!is.na(ps_drug$prod_ai) & grepl(pat, ps_drug$prod_ai, ignore.case = TRUE, perl = TRUE))
  data.table(drug = nm, n = uniqueN(ps_drug$primaryid[hit]))
}))
breakdown[, pct := round(n / n_ssri * 100, 1)]
write_res(breakdown, "ssri_drug_breakdown.csv")

by_year <- faers_get(ssri, "demo")[, .(n = uniqueN(primaryid)), by = year][order(year)]
yr_drug <- rbindlist(lapply(names(ssri_drugs), function(nm)
  faers_get(drug_filter(ssri, ssri_drug_pattern[[nm]]), "demo")[, .(n = uniqueN(primaryid)), by = year][, drug := nm]))
write_res(merge(by_year, dcast(yr_drug, year ~ drug, value.var = "n", fill = 0L), by = "year", all = TRUE), "ssri_by_year.csv")

if (WRITE_FULL_SSRI_CASE_DUMP)
  write_res(build_case_list(ssri, faers_get(ssri, "demo")$primaryid, codrug), "ssri_cases_ps.csv")

log_msg("=== 02 COMPLETE | ", elapsed(t0), " ===")
