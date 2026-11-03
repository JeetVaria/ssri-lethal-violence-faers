# ==============================================================================
# R/signals.R - disproportionality analysis.
#
# PT counting method. FAERS stores the same MedDRA PT in different cases over
# time ("COMPLETED SUICIDE" before 2012, "Completed suicide" after). Every PT is
# normalized with pt_std() in R/definitions.R and a report counts once per
# normalized PT (union of report IDs across raw spellings). The same report IDs
# are reused to build the case lists, so `a` always equals the number of rows
# in the matching case list.
#
# 2 x 2 table for a drug group (the drug of interest) against a reference set (.full):
#   a = group reports with the PT         b = group reports without it
#   c = other reports with the PT         d = other reports without it
# ==============================================================================

# Lookup raw PT string -> standardized preferred term, built once per reference
# set. Case folding handles the spelling FAERS changed in 2012; it does not
# handle a reaction recorded at Lowest Level Term (LLT), which is a different
# term for the same PT, so llt_map() (R/definitions.R) rolls those up.
#
# `applied` records the mappings that actually fired, with the number of reports
# behind each, so every count the roll-up changed can be traced to its term.
pt_norm_lookup <- function(obj = NULL, applied = NULL, reac = NULL) {
  if (is.null(reac)) reac <- faers_get(obj, "reac")
  raw  <- unique(reac$pt)
  norm <- pt_std(raw)                      # the normalization llt_map() uses too
  m    <- llt_map()
  hit  <- match(norm, names(m))
  if (!is.null(applied) && any(!is.na(hit))) {
    n <- reac[, .(reports = uniqueN(primaryid)), by = pt]
    fwrite(data.table(llt = norm[!is.na(hit)], pt = unname(m[hit[!is.na(hit)]]),
                      raw_term = raw[!is.na(hit)])[
             n, reports := i.reports, on = c(raw_term = "pt")][order(-reports)],
           res_path(applied))
  }
  norm[!is.na(hit)] <- m[hit[!is.na(hit)]]
  setNames(norm, raw)
}

# (primaryid, pt, pt_norm) for one FAERS object
build_obj_reac_norm <- function(obj, norm_lookup) {
  rt <- faers_get(obj, "reac")[, .(primaryid, pt)]
  rt[, pt_norm := norm_lookup[pt]]
  rt[!is.na(pt_norm) & pt_norm != ""]
}

# (pt_norm, a, ids): one row per normalized PT
pt_ids_table <- function(reac_norm) {
  reac_norm[, .(a = uniqueN(primaryid), ids = list(unique(primaryid))), by = pt_norm]
}

NA_SIGNAL_ROW <- data.table(
  ror = NA_real_, ror_ci_low = NA_real_, ror_ci_high = NA_real_,
  prr = NA_real_, prr_ci_low = NA_real_, prr_ci_high = NA_real_,
  chisq = NA_real_, chisq_pvalue = NA_real_,
  bcpnn_mcmc_ic = NA_real_, bcpnn_mcmc_ic_ci_low = NA_real_, bcpnn_mcmc_ic_ci_high = NA_real_,
  ebgm = NA_real_, ebgm_ci_low = NA_real_, ebgm_ci_high = NA_real_)

# ROR, PRR and chi-squared are closed forms on one table, and BCPNN's posterior
# is drawn from that table's own cells, so phv_signal() computes all four row by
# row. MGPS is not of that kind: its prior is a property of the whole reference
# set, so it is scored here against the prior that use_ebgm_prior() fitted for
# the set this row is being compared with. See R/ebgm.R.
run_phv <- function(a, b, c, d) {
  out <- as.data.table(phv_signal(a, b, c, d, methods = PHV_PER_ROW,
                                  n_mcmc = PHV_N_MCMC, BPPARAM = PHV_BPPARAM))
  if (!"ebgm" %chin% PHV_METHODS) return(out)
  if (is.null(EBGM_PRIOR))
    stop("no MGPS prior is in force: call use_ebgm_prior() for this reference ",
         "set before computing any signal against it.", call. = FALSE)
  # phv_signal() returns the expected count for this table, so MGPS is scored on
  # the same E the other algorithms used rather than on a second calculation
  cbind(out, score_ebgm(EBGM_PRIOR, a = a, E = out$expected))
}

# Signal thresholds. A signal is reported as significant only when all four
# algorithms agree (full_consensus):
#   ROR  lower 95% CI > 1
#   PRR  >= 2 with a >= 3
#   BCPNN IC025 > 0
#   MGPS EBGM05 > 2
add_signals <- function(dt) {
  num_cols <- c("a", "ror", "ror_ci_low", "ror_ci_high", "prr", "prr_ci_low", "prr_ci_high",
                "chisq", "ebgm", "ebgm_ci_low", "bcpnn_mcmc_ic", "bcpnn_mcmc_ic_ci_low",
                "bcpnn_norm_ic", "bcpnn_norm_ic_ci_low")
  for (col in num_cols) if (col %in% names(dt)) dt[, (col) := as.numeric(get(col))]
  # Quick-test mode skips the Bayesian methods; keep the columns so the rules still run.
  for (col in c("chisq", "ebgm", "ebgm_ci_low", "ebgm_ci_high")) if (!col %in% names(dt)) dt[, (col) := NA_real_]
  if (!any(c("bcpnn_mcmc_ic_ci_low", "bcpnn_norm_ic_ci_low") %in% names(dt)))
    for (col in c("bcpnn_mcmc_ic", "bcpnn_mcmc_ic_ci_low", "bcpnn_mcmc_ic_ci_high")) dt[, (col) := NA_real_]
  # one of the two is guaranteed by the block above
  IC025 <- if ("bcpnn_mcmc_ic_ci_low" %in% names(dt)) "bcpnn_mcmc_ic_ci_low" else "bcpnn_norm_ic_ci_low"
  dt[, ror_signal     := !is.na(ror_ci_low) & ror_ci_low > 1]
  dt[, prr_signal     := !is.na(prr) & prr >= 2 & !is.na(a) & a >= 3]
  dt[, ic_signal      := !is.na(get(IC025)) & get(IC025) > 0]
  dt[, ebgm_signal    := !is.na(ebgm_ci_low) & ebgm_ci_low > 2]
  dt[, full_consensus := ror_signal & prr_signal & ic_signal & ebgm_signal]
  # With a = 0 the Bayesian priors still return small EBGM/IC values; blank them.
  zero_cols <- c("ror", "ror_ci_low", "ror_ci_high", "prr", "prr_ci_low", "prr_ci_high",
                 "chisq", "ebgm", "ebgm_ci_low", "ebgm_ci_high",
                 "bcpnn_mcmc_ic", "bcpnn_mcmc_ic_ci_low", "bcpnn_mcmc_ic_ci_high",
                 "bcpnn_norm_ic", "bcpnn_norm_ic_ci_low")
  for (col in zero_cols) if (col %in% names(dt)) dt[!is.na(a) & a == 0, (col) := NA_real_]
  sig_cols <- c("ror_signal", "prr_signal", "ic_signal", "ebgm_signal", "full_consensus")
  for (col in sig_cols) if (col %in% names(dt)) dt[!is.na(a) & a == 0, (col) := FALSE]
  dt
}

# Signals for each row of a (pt_norm, a, ids) table. No estimate when a < 3.
signals_from_ids_table <- function(idt, full_ids_table, n_obj, n_full) {
  rbindlist(lapply(seq_len(nrow(idt)), function(i) {
    a <- idt$a[i]; b <- n_obj - a
    full_row <- full_ids_table[pt_norm == idt$pt_norm[i]]
    a_full   <- if (nrow(full_row)) full_row$a[1] else a
    c <- a_full - a
    d <- n_full - n_obj - c
    sig <- if (a < 3) copy(NA_SIGNAL_ROW) else run_phv(a, b, c, d)
    sig[, `:=`(pt = idt$pt_norm[i], a = a, b = b, c = c, d = d)]
    sig[, ids := list(list(idt$ids[[i]]))]
    sig
  }), fill = TRUE)
}

# Signals for the target PTs of one drug group. Targets with no reports get an
# explicit a = 0 row so every group returns the same set of PTs.
compute_pt_signals <- function(obj, full_ids_table, n_obj, n_full, norm_lookup, pts = targets) {
  idt <- pt_ids_table(build_obj_reac_norm(obj, norm_lookup))[pt_norm %in% pts]
  missing_pts <- setdiff(pts, idt$pt_norm)
  if (length(missing_pts))
    idt <- rbindlist(list(idt, data.table(pt_norm = missing_pts, a = 0L,
                                          ids = replicate(length(missing_pts), integer(0), simplify = FALSE))),
                     use.names = TRUE)
  out <- signals_from_ids_table(idt, full_ids_table, n_obj, n_full)
  out[, group := groups[pt]]
  add_signals(out)
}

get_ids_for_pt <- function(signal_dt, target_pt) {
  row <- signal_dt[pt == target_pt]
  if (!nrow(row)) return(integer(0))
  row$ids[[1]]
}

# Background union counts for the Suicidality / Homicidality groups.
group_background <- function(master_obj, norm_lookup) {
  rn <- build_obj_reac_norm(master_obj, norm_lookup)
  rbindlist(lapply(names(group_defs), function(g) {
    ids <- unique(rn[pt_norm %in% group_defs[[g]]]$primaryid)
    data.table(group = g, a = length(ids), ids = list(ids))
  }))
}

# Group composite: a report counts once if it carries ANY PT of the group.
compute_group_signal <- function(obj, group_label, bg, n_obj, n_full, norm_lookup, drug_label) {
  rn  <- build_obj_reac_norm(obj, norm_lookup)
  ids <- unique(rn[pt_norm %in% group_defs[[group_label]]]$primaryid)
  a <- length(ids)
  a_full <- bg[group == group_label]$a[1]
  b <- n_obj - a; c <- a_full - a; d <- n_full - n_obj - c
  sig <- if (a < 3) copy(NA_SIGNAL_ROW) else run_phv(a, b, c, d)
  sig[, `:=`(group = group_label, drug = drug_label, a = a, b = b, c = c, d = d)]
  sig[, ids := list(list(ids))]
  sig
}

# Homicide-suicide (H-S) composite: reports co-reporting Homicide AND Completed
# Suicide. n_hs_master = H-S reports in the whole reference set.
# min_a: groups with fewer co-reports than this get an empty (a, b, c, d) row.
hs_composite <- function(a, n_group, n_hs_master, n_master, min_a = 3L) {
  b <- n_group - a
  c <- max(0, n_hs_master - a)
  d <- n_master - a - b - c
  if (a < min_a) {
    # The columns the fitted branch returns, in the same order. Returning fewer
    # made the column set depend on the data: a group below the floor wrote a
    # file without expected, prr_ci_low, prr_ci_high and chisq_pvalue, and
    # tools/compare_to_reference.R reads a changed column set as a failure to
    # reproduce.
    out <- cbind(data.table(expected = NA_real_), copy(NA_SIGNAL_ROW),
                 data.table(a = as.integer(a), b = b, c = c, d = d))
  } else {
    out <- run_phv(a, b, c, d)
    out[, `:=`(a = a, b = b, c = c, d = d)]
  }
  add_signals(out)
}

hs_master_count <- function(full_ids_table) {
  length(id_intersect(full_ids_table[pt_norm == "Homicide"]$ids[[1]],
                   full_ids_table[pt_norm == "Completed Suicide"]$ids[[1]]))
}

# ── comparing two RORs ────────────────────────────────────────────────────────
# Normalized difference (ND) = (x - y) / (x + y), bounded -1..+1.
nd_index <- function(x, y) ifelse(!is.na(x) & !is.na(y) & (x + y) > 0, (x - y) / (x + y), NA_real_)

# Two-sided z-test on log(ROR); each SE derived from its 95% CI.
# The two RORs are treated as independent; monopharmacy is a subset of the
# unrestricted set, so the test is conservative (stated in the Methods).
ror_z_test <- function(ror1, lo1, hi1, ror2, lo2, hi2) {
  se <- function(lo, hi) (log(hi) - log(lo)) / (2 * 1.96)
  z  <- (log(ror1) - log(ror2)) / sqrt(se(lo1, hi1)^2 + se(lo2, hi2)^2)
  list(z = z, p = 2 * pnorm(abs(z), lower.tail = FALSE))
}
