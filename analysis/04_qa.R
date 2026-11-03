# ==============================================================================
# 04_qa.R - pipeline quality assurance (Supplementary "QA").
#
# QA pairs with an expected signal: well-established drug-event pairs run
# through exactly the same path as the SSRIs (drug_filter -> compute_pt_signals, PS reports only,
# blacklisted reference set). Recovering them confirms drug matching,
# deduplication and the 2 x 2 construction work.
#   statin -> rhabdomyolysis, myopathy
#   warfarin -> haemorrhage, GI haemorrhage, cerebral haemorrhage
#   ACE inhibitor -> angioedema
# QA pairs with an expected null: vitamin D across every target PT (computed
# here, since vitamin D takes no part in the analysis) and ingrowing nail for
# the pooled SSRIs.
#
# PASS rule (pre-specified for these large, unambiguous pairs): a >= 3, lower
# 95% CI of the ROR > 1 and ROR >= expect_floor in data/qa_pairs.csv.
#
# Outputs: qa_signals.csv, qa_supplementary.csv
# ==============================================================================
source("R/load_all.R")
need_faers()
t0 <- Sys.time()
log_msg("=== 04 QA ===")
qa_def <- fread("data/qa_pairs.csv")

master   <- apply_blacklist(readRDS(require_res("master_dedup.rds")), "master")
n_master <- nrow(faers_get(master, "demo"))
reac     <- faers_get(master, "reac")[, .(primaryid, pt)]
# the same roll-up the study outcomes get, so a QA pair is ascertained the way
# a study outcome is; reac is passed so the reaction table is read once
norm_master <- pt_norm_lookup(reac = reac)
# Only the QA PTs are looked up, so the background table is restricted to them.
reac[, pt_norm := norm_master[pt]]
bg <- reac[pt_norm %in% union(qa_def$pt, targets), .(a = uniqueN(primaryid), ids = list(unique(primaryid))), by = pt_norm]
rm(reac); gc2()
# The QA pairs are compared with the same unrestricted reference set as step 02,
# so they are scored against the same prior: both the weighted points and the
# fitted hyperparameters are cached, so this reuses 02's fit rather than
# producing a second one and overwriting the table that describes it.
use_ebgm_prior(master, norm_master, n_master, "full")

qa <- rbindlist(lapply(unique(qa_def$drug), function(dr) {
  def <- qa_def[drug == dr]
  sub <- drug_filter(master, terms_pattern(qa_drugs[[dr]]))
  n   <- nrow(faers_get(sub, "demo"))
  log_msg("  ", def$label[1], ": ", sfmt(n), " PS reports")
  res <- compute_pt_signals(sub, bg, n, n_master, norm_master, pts = def$pt)[, ids := NULL]
  # phv_signal returns its own numeric "expected" column, which is why ours is
  # called expectation; drop any clash so merge() keeps our labels.
  clash <- intersect(names(res), c("qa_drug", "expect_floor", "expectation", "n_group"))
  if (length(clash)) res[, (clash) := NULL]
  merge(res, def[, .(pt, qa_drug = label, expect_floor, expectation)], by = "pt")[, n_group := n]
}), fill = TRUE)
qa[, group := "QA (expected signal)"]
qa[, qa_flag := fifelse(!is.na(ror) & a >= 3 & ror_ci_low > 1 & ror >= expect_floor, "PASS",
                fifelse(!is.na(ror) & a >= 3 & ror_ci_low > 1, "PASS (weak vs expected)",
                fifelse(is.na(a) | a < 3, "FLAG: <3 cases (check drug pattern / PT name)", "FLAG: no signal (ROR CI includes 1)")))]
write_res(qa[, .(qa_drug, pt, n_group, a, ror, ror_ci_low, ror_ci_high, prr, ebgm, ebgm_ci_low,
                 bcpnn_mcmc_ic, bcpnn_mcmc_ic_ci_low, full_consensus, expect_floor, expectation, qa_flag)],
          "qa_signals.csv")

# ── supplementary table: all QA pairs ─────────────────────────────────────────
# Vitamin D is a quality-assurance drug only and takes no part in the analysis,
# so its group is built here rather than alongside the SSRIs in step 02. It is
# run through exactly the same path: drug_filter -> compute_pt_signals, primary
# suspect only, against the blacklisted reference set.
vd_sub  <- drug_filter(master, terms_pattern(qa_drugs$vitamind))
n_vd_qa <- nrow(faers_get(vd_sub, "demo"))
log_msg("  Vitamin D (expected null): ", sfmt(n_vd_qa), " PS reports")
vd_res  <- compute_pt_signals(vd_sub, bg, n_vd_qa, n_master, norm_master, pts = targets)
if ("ids" %in% names(vd_res)) vd_res[, ids := NULL]
vd_res[, qa_drug := "Vitamin D"]

neg <- rbindlist(list(
  vd_res,
  read_res("ssri_targeted_signals.csv")[pt == "Ingrowing Nail"][, qa_drug := "SSRI (pooled)"]), fill = TRUE)
fmt2 <- function(x) ifelse(is.na(x), "NA", sprintf("%.2f", x))
supp <- rbindlist(list(qa, neg), fill = TRUE)[, .(
  PT = pt, DRUG = qa_drug, a, b, c, d,
  `ROR (95% CI)`  = paste0(fmt2(ror), " (", fmt2(ror_ci_low), " - ", fmt2(ror_ci_high), ")"),
  `PRR (chi2)`      = paste0(fmt2(prr), " (", ifelse(is.na(chisq), "NA", sprintf("%.1f", chisq)), ")"),
  `MGPS (EBGM05)` = paste0(fmt2(ebgm), " (", fmt2(ebgm_ci_low), ")"),
  `BCPNN (IC025)` = paste0(fmt2(bcpnn_mcmc_ic), " (", fmt2(bcpnn_mcmc_ic_ci_low), ")"),
  CONSENSUS = ifelse(full_consensus, "Yes", "No"))]
write_res(supp, "qa_supplementary.csv")
log_msg("=== 04 COMPLETE | ", elapsed(t0), " ===")
