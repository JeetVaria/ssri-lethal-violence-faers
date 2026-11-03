# ==============================================================================
# tools/refit_ebgm_prior.R - refit the MGPS prior from many starting values and
# report whether the likelihood has one optimum or several.
#
# WHY. The prior is the one quantity in this analysis that is estimated rather
# than counted, and a mixture likelihood can have more than one optimum. The
# pipeline runs every starting value in EBGM_THETA_INIT and warns when they
# disagree, but a start that errors inside hyperEM leaves no comparison at all.
# This tool refits from nine widely spread starting values. On both reference sets
# every start that converged reached the same optimum, agreeing to better than
# 0.001% in every hyperparameter; see ebgm_prior_restarts_full.csv and
# ebgm_prior_restarts_mono.csv.
#
# It is cheap because the weighted points are cached: the matrix is not rebuilt.
#
# HOW IT DECIDES. Comparing fitted parameters is the weaker test, because the
# mixture is invariant under exchanging its two components and because an
# absolute tolerance means nothing when beta1 is of order 0.006. The decisive
# quantity is the log-likelihood: two labelings of one optimum give exactly the
# same value, whereas two genuinely different optima do not. So the converged
# fits are grouped by log-likelihood, and only within a group are the parameters
# compared - relative to their own size, not against a fixed tolerance.
#
# The starting values below are deliberately NOT drawn from the fitted answer.
# Starting next to a known solution and arriving back at it shows only that the
# solution is a local attractor, which was never in doubt. They span several
# orders of magnitude in every parameter, and the first is DuMouchel's published
# prior, which is also the pipeline's first start.
#
# It writes nothing the pipeline reads - not ebgm_prior_<tag>.csv and not the
# cached prior - so it cannot change a published number. Output goes to
# ebgm_prior_restarts_<tag>.csv, including the error message from any start that
# failed, since characterising those failures is half the point.
#
# Usage:  FAERS_RESULTS_DIR=full_run Rscript tools/refit_ebgm_prior.R full
#         FAERS_RESULTS_DIR=full_run Rscript tools/refit_ebgm_prior.R mono
# ==============================================================================
source("R/load_all.R")

tag <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(tag)) tag <- "full"
f <- res_path(paste0("ebgm_points_", tag, ".rds"))
if (!file.exists(f))
  stop("no cached points for the ", tag, " reference set: ", f,
       "\nRun the pipeline for this period first - this tool refits the prior, ",
       "it does not build the matrix.", call. = FALSE)

pts   <- as.data.table(readRDS(f))
stamp <- attr(readRDS(f), "ebgm_stamp")
log_msg("=== REFIT ", toupper(tag), " MGPS PRIOR ===")
log_msg("  ", sfmt(nrow(pts)), " weighted points (", sfmt(sum(pts$weight)), " cells) from ",
        basename(f), ", built from ", if (is.null(stamp)) "an unrecorded run" else stamp)

starts <- data.frame(
  alpha1 = c(0.20, 0.02, 1.00, 0.05, 0.50, 0.01, 2.00, 0.10, 0.30),
  beta1  = c(0.10, 0.005, 1.00, 0.01, 0.20, 0.10, 0.50, 0.05, 2.00),
  alpha2 = c(2.00, 0.50, 5.00, 1.00, 3.00, 0.20, 0.50, 8.00, 1.50),
  beta2  = c(4.00, 0.50, 5.00, 2.00, 3.00, 1.00, 2.00, 8.00, 0.50),
  P      = c(1/3,  0.05, 0.50, 0.10, 0.70, 0.30, 0.20, 0.40, 0.60))

d <- as.data.frame(pts[, .(N, E, weight)])
# the same tolerance the pipeline derived: the minimum over EBGM_THETA_INIT, not
# over the starting values used here, or the comparison would be against a
# different convergence criterion than the published fit used
ll0 <- min(vapply(seq_len(nrow(EBGM_THETA_INIT)), function(i)
  openEBGM::negLLzeroSquash(as.numeric(EBGM_THETA_INIT[i, ]),
                            ni = d$N, ei = d$E, wi = d$weight), numeric(1)))
tol <- max(EBGM_EM_LL_TOL_MIN, EBGM_EM_LL_REL * abs(ll0))
pub <- tryCatch(read_res(paste0("ebgm_prior_", tag, ".csv"), required = FALSE),
                error = function(e) NULL)
if (!is.null(pub) && "ll_tol" %chin% names(pub) && !isTRUE(all.equal(pub$ll_tol[1], tol)))
  log_msg("  WARNING: the published fit used a tolerance of ", sprintf("%.4g", pub$ll_tol[1]),
          " and this refit derives ", sprintf("%.4g", tol), " - the two are not comparable")
log_msg("  convergence when the log-likelihood moves by less than ", sprintf("%.4g", tol),
        " for ", EBGM_EM_CONSECUTIVE, " iterations")

# seeded once before the loop, as the pipeline seeds once before its own, so
# that the first start reproduces the published fit exactly
set.seed(PHV_SEED)
rows <- rbindlist(lapply(seq_len(nrow(starts)), function(i) {
  t0  <- Sys.time()
  fit <- tryCatch(
    openEBGM::hyperEM(data = d, theta_init_vec = as.numeric(starts[i, ]),
                      squashed = TRUE, zeroes = TRUE, N_star = NULL,
                      method = "score", LL_tol = tol,
                      consecutive = EBGM_EM_CONSECUTIVE,
                      param_lower = EBGM_PAR_MIN, param_upper = EBGM_PAR_MAX,
                      print_level = 0L, max_iter = EBGM_EM_MAX_ITER),
    error = function(e) structure(list(msg = conditionMessage(e)), class = "refit_failed"))
  mins  <- round(as.numeric(Sys.time() - t0, units = "mins"), 1)
  begin <- paste(sprintf("%.4g", as.numeric(starts[i, ])), collapse = " ")
  if (inherits(fit, "refit_failed")) {
    log_msg("  start ", i, " (", begin, "): no prior after ", mins, " min - ", fit$msg)
    return(data.table(start = i, from = begin, converged = FALSE, minutes = mins,
                      error = fit$msg))
  }
  th <- canon_theta(as.numeric(fit$estimates))
  log_msg("  start ", i, " (", begin, "): ", fit$iters, " iterations in ", mins,
          " min, logL ", sprintf("%.8g", fit$maximum), " -> ",
          paste(sprintf("%.5f", th), collapse = " "))
  data.table(start = i, from = begin, converged = TRUE, minutes = mins,
             iterations = fit$iters, log_likelihood = fit$maximum,
             score_norm = fit$score_norm,
             alpha1 = th[1], beta1 = th[2], alpha2 = th[3], beta2 = th[4], P = th[5],
             mean_high = th[1] / th[2], mean_low = th[3] / th[4], error = NA_character_)
}), fill = TRUE)

ok <- rows[converged == TRUE]
log_msg("  ", nrow(ok), " of ", nrow(starts), " starting values returned a prior")

if (nrow(ok) >= 2L) {
  # One optimum or several. Two labelings of the same optimum share a
  # log-likelihood exactly; different optima do not. The grouping width is ten
  # times the convergence tolerance, so a fit that merely stopped slightly
  # short is not counted as a different answer.
  setorder(ok, -log_likelihood)
  best <- ok$log_likelihood[1]
  ok[, at_best := (best - log_likelihood) <= 10 * tol]
  log_msg("  log-likelihood range ", sprintf("%.8g", min(ok$log_likelihood)), " to ",
          sprintf("%.8g", best), "; ", sum(ok$at_best), " of ", nrow(ok),
          " reached the best value found")
  if (!all(ok$at_best))
    log_msg("  MORE THAN ONE OPTIMUM: ", sum(!ok$at_best), " fit(s) settled ",
            sprintf("%.4g", best - min(ok$log_likelihood)),
            " below the best - the prior is not determined by the data alone, and the ",
            "published fit is whichever optimum its starting value happened to reach")
  at <- ok[at_best == TRUE]
  if (nrow(at) >= 2L) {
    # relative spread, because an absolute tolerance is meaningless when beta1
    # is of order 0.006 - autoHyper's 0.2 would admit beta1 anywhere in
    # [1e-5, 0.2] and the high component's mean anywhere from 0.4 to 8,000
    rel <- vapply(c("alpha1", "beta1", "alpha2", "beta2"),
                  function(c) diff(range(at[[c]])) / stats::median(at[[c]]), numeric(1))
    pabs <- diff(range(at$P))
    log_msg("  within that optimum, relative spread: ",
            paste(sprintf("%s %.2g%%", names(rel), 100 * rel), collapse = "  "),
            "  P +/-", sprintf("%.4f", pabs))
    log_msg(if (all(rel < 0.02) && pabs < 0.01)
      "  the fits agree to better than 2% in every parameter: the prior is identified" else
      "  THE FITS DISAGREE within one optimum, which should not happen - check score_norm")
  }
} else {
  log_msg("  too few fits to compare: this says nothing about whether the prior is identified")
}
write_res(rows, paste0("ebgm_prior_restarts_", tag, ".csv"))
log_msg("=== REFIT COMPLETE ===")
