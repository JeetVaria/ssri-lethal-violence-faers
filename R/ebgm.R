# ==============================================================================
# R/ebgm.R - MGPS, the multi-item gamma Poisson shrinker.
#
# The gamma-mixture prior belongs to the reference set, not to one drug-event
# pair. DuMouchel (1999) estimates it across the whole drug-event matrix, with
# the cells that hold no report counted as the zeros they are: they are the bulk
# of the matrix and they are what tells the prior how often a drug and an event
# do not occur together. Fitting it to a single 2 x 2 table instead makes each
# pair its own background, which shrinks a large ratio resting on few reports to
# below the null and loses signals the other three algorithms report.
#
# Why this is not done through faers. phv_ebgm(a, b, c, d, alpha, theta_init,
# squashing) takes theta_init as a starting guess for a fit it then performs on
# the table it was handed; it has no argument that accepts an estimated prior,
# so there is no way to call it such that the prior comes from the matrix. Its
# own documentation refers the reader to openEBGM for the hyperparameters, which
# is the one thing openEBGM is used for here. ROR, PRR, chi-squared and BCPNN
# all still come from faers::phv_signal.
#
# The matrix is built here rather than with openEBGM::processRaw(), which needs
# the whole report x drug x PT frame resident at once; at FAERS scale that is
# some hundreds of millions of rows.
#
# Two things the prior is NOT, and which the Methods must say:
#   - its drug axis is one canonical drug name per report (drug_canon(), so a
#     brand and its generic count once), over every drug role, with the
#     placeholders that mean "no co-medication" dropped first (drop_non_drugs(),
#     as report_drug_counts() does) - every placeholder here, NON_DRUG_ANY, not
#     only the excluded ones, because UNKNOWN is not a drug either and would be
#     one of the highest-count rows of the axis the prior is fitted over. The
#     filter runs twice, on the raw name and again on the canonical one. The
#     second pass catches nothing, since names key as written, and is kept
#     only as a guard against a canonical form equal to a placeholder. The study rows
#     are scored as drug GROUPS - pooled SSRIs, each SSRI, varenicline, primary
#     suspect only - so a study row is not itself a cell of this matrix. Its E
#     comes from its own 2 x 2 table, as it does for the other three algorithms.
#   - it is not stratified by age, sex or year.
#
#   EBGM          geometric mean of the posterior of the true relative risk
#   ebgm_ci_low   5th percentile of that posterior - EBGM05 in Table 2
#   ebgm_ci_high  95th percentile. phv_signal() wrote the 97.5th here, so this
#                 column changes meaning as well as value; nothing downstream
#                 reads it.
#
# References: DuMouchel W. Am Stat 1999;53:177-90.
#             Canida T, Ihrie J. openEBGM. R Journal 2017;9(2):436-49.
# ==============================================================================

# ── the observed cells ────────────────────────────────────────────────────────
# One (N, E) point for every drug-event pair that occurs at least once:
#
#   N  = reports naming the drug and carrying the PT
#   E  = (reports naming the drug) x (reports carrying the PT) / n_total
#
# nd and np are carried out of here as well, because the cells holding no report
# are reconstructed from the two margins rather than enumerated.
#
# Built in blocks of reports because the report x drug x PT expansion runs to
# some hundreds of millions of rows, and drug and PT are replaced by integer
# codes first for the same reason. The same remainder is taken on both sides, so
# a report's drugs and its reactions always land in the same block; blocks hold
# disjoint sets of reports, so a pair's counts add across them.
ebgm_pair_counts <- function(obj, norm_lookup, n_total, chunks = EBGM_CHUNKS) {
  # narrowed to two columns before the filter, as report_drug_counts() does:
  # row-subsetting the whole ~22-column table first copies 176 bytes a row over
  # 80 million rows, and this stage already holds the database in memory.
  dr <- drop_non_drugs(faers_get(obj, "drug")[, .(primaryid, drugname)],
                       drop_terms = NON_DRUG_ANY)[
          , .(primaryid, drug = drug_canon(drugname))]
  # Again on the canonical name. Since names key as written this catches
  # nothing - a name reaching here equal to a placeholder was already equal to
  # one before canonicalisation, and the first filter took it. Kept as a guard:
  # it cannot remove a real drug, because no dictionary term canonicalises to a
  # placeholder.
  dr <- dr[!non_drug_key(drug) %chin% NON_DRUG_ANY]
  dr <- unique(dr[!is.na(drug) & nzchar(drug)])
  rc <- unique(build_obj_reac_norm(obj, norm_lookup)[, .(primaryid, pt = pt_norm)])

  # faers holds primaryid as character (unify_ascii), so it is converted before
  # the remainder is taken. An id that will not convert gives a block of NA,
  # which matches no block, and that report's pairs would drop out in silence.
  blk <- function(id) {
    i <- bit64::as.integer64(id)
    if (anyNA(i))
      stop("primaryid does not convert to an integer on ", sum(is.na(i)),
           " rows, so their reports cannot be blocked", call. = FALSE)
    as.integer(i %% chunks)
  }

  dlv <- sort(unique(dr$drug)); plv <- sort(unique(rc$pt))
  dr[, `:=`(di = match(drug, dlv), drug = NULL, ck = blk(primaryid))]
  rc[, `:=`(ti = match(pt,   plv), pt   = NULL, ck = blk(primaryid))]
  log_msg("    ", sfmt(length(dlv)), " drug names x ", sfmt(length(plv)), " PTs = ",
          sfmt(as.numeric(length(dlv)) * length(plv)), " cells")

  # one aggregation at the end rather than one per block: re-grouping a table
  # that grows to the full pair count, once per block, is the dominant cost
  parts <- vector("list", chunks)
  for (k in 0L:(chunks - 1L)) {
    parts[[k + 1L]] <- merge(dr[ck == k, .(primaryid, di)], rc[ck == k, .(primaryid, ti)],
                             by = "primaryid", allow.cartesian = TRUE)[, .(N = .N), by = .(di, ti)]
    log_msg("    block ", k + 1L, "/", chunks, ": ", sfmt(nrow(parts[[k + 1L]])), " pairs")
  }
  acc <- rbindlist(parts)[, .(N = sum(N)), by = .(di, ti)]
  rm(parts); gc2()

  # dr and rc are unique on (primaryid, code) - unique() runs before the code is
  # assigned and match() is injective - so .N by code counts distinct reports.
  # as.numeric() before the product: .N is integer, and an integer product
  # overflows to NA above 2^31 - 1, which here is every pair with E above ~109.
  nd <- dr[, .(nd = as.numeric(.N)), by = di]
  np <- rc[, .(np = as.numeric(.N)), by = ti]
  setkey(acc, di); acc[nd, nd := i.nd, on = "di"]
  setkey(acc, ti); acc[np, np := i.np, on = "ti"]

  list(obs     = acc[, .(N = as.numeric(N), E = nd * np / n_total, nd, np)],
       nd      = nd$nd,
       np      = np$np,
       n_drugs = length(dlv),
       n_pts   = length(plv))
}

# ── the cells holding no report ───────────────────────────────────────────────
# Every cell's E is nd x np / n_total, so all cells sharing a pair of margin
# counts share one E. That is what makes the zero stratum tractable: it is never
# enumerated. The two margins are grouped, the groups are crossed, each class
# knows how many cells it holds (drugs in the class x PTs in the class) and how
# many of them were observed, and the difference is a count of zero cells all
# carrying the same E.
#
# Margin counts are grouped on a log scale when there are more distinct values
# than EBGM_MARGIN_BINS, the weighted mean standing for the group. At a thousand
# groups a side the counts inside one group span about one percent, and the only
# error that survives is in the expected count charged to the observed pairs
# being subtracted - the sum over every cell of a class is exact, because it
# factorises. Measured against a complete build of the matrix, the sum of
# expected counts over the zero cells is then within 1e-6 of exact; the
# threshold where that error reaches 1e-3 is around a hundred groups. The
# equal-width cut on a skewed log density leaves roughly two thirds of the
# groups occupied, so the effective resolution is lower than the setting.
margin_groups <- function(v, bins = EBGM_MARGIN_BINS) {
  # a margin count is a number of reports, so it is finite and at least 1. An NA
  # would be carried through to a class with an NA expected count and the right
  # weight, which every total below would still balance.
  if (!length(v) || anyNA(v) || any(!is.finite(v)) || any(v < 1))
    stop("a margin of the drug-event matrix holds ", sum(is.na(v) | !is.finite(v) | v < 1),
         " counts that are not a positive number of reports", call. = FALSE)
  t <- data.table(v = v)[, .(k = as.numeric(.N)), by = v][order(v)]
  t[, g := if (.N <= bins) seq_len(.N) else as.integer(cut(log(v), breaks = bins, labels = FALSE))]
  list(t = t, tab = t[, .(v = sum(v * k) / sum(k), k = sum(k)), by = g])
}

ebgm_zero_points <- function(mx, n_total, bins = EBGM_ZERO_BINS) {
  # not named gd/gp: data.table evaluates i with the table's own columns in
  # scope, so a name shared with a column of cls would resolve to the column
  md <- margin_groups(mx$nd); mp <- margin_groups(mx$np)
  exact <- nrow(md$tab) == nrow(md$t) && nrow(mp$tab) == nrow(mp$t)
  log_msg("    margins: ", sfmt(nrow(md$t)), " distinct drug counts and ",
          sfmt(nrow(mp$t)), " distinct term counts", if (exact)
            ", neither grouped, so every cell keeps its own expected count" else
            paste0(", grouped to ", sfmt(nrow(md$tab)), " x ", sfmt(nrow(mp$tab)),
                   " - cells in a class then share one expected count"))
  cls <- CJ(gd = md$tab$g, gp = mp$tab$g)
  cls[md$tab, `:=`(nd = i.v, cd = i.k), on = c(gd = "g")]
  cls[mp$tab, `:=`(np = i.v, cp = i.k), on = c(gp = "g")]
  cls[, cells := cd * cp]

  o <- mx$obs[, .(nd, np)]
  o[md$t, gd := i.g, on = c(nd = "v")]
  o[mp$t, gp := i.g, on = c(np = "v")]
  o <- o[, .(obs = as.numeric(.N)), by = .(gd, gp)]
  cls[o, obs := i.obs, on = c("gd", "gp")]
  cls[is.na(obs), obs := 0]

  # A class cannot hold more observed pairs than it holds cells. If it does, the
  # margins and the pair table disagree and the whole matrix is wrong.
  if (any(cls$obs > cls$cells))
    stop("the drug-event matrix is inconsistent: ", sum(cls$obs > cls$cells),
         " margin classes hold more observed pairs than they hold cells", call. = FALSE)
  # Both totals below are compared exactly rather than within a tolerance. Every
  # term is a count: cells is cd x cp, at most the size of the matrix, and the
  # running total stays far below 2^53, so each addition is exact in double
  # precision and sum()'s accumulator is wider still. A tolerance of even one
  # cell would pass the case that matters - an observed pair whose margin count
  # belongs to no class takes gd = NA, is never subtracted, and is published as
  # a cell holding no report, moving the total by exactly one.
  n_cells <- sum(cls$cells)
  if (n_cells != as.numeric(mx$n_drugs) * mx$n_pts)
    stop("the margin classes cover ", sfmt(n_cells), " cells but the matrix has ",
         sfmt(as.numeric(mx$n_drugs) * mx$n_pts), call. = FALSE)

  z <- cls[cells > obs, .(E = nd * np / n_total, w = cells - obs)]
  if (!nrow(z))
    stop("every cell of the drug-event matrix holds a report, so there is no zero ",
         "stratum and the prior cannot be fitted with the zero cells included",
         call. = FALSE)
  z[, b := as.integer(cut(log(E), breaks = bins, labels = FALSE))]
  out <- z[, .(N = 0, E = sum(E * w) / sum(w), weight = sum(w)), by = b][, .(N, E, weight)]
  if (sum(out$weight) != n_cells - nrow(mx$obs))
    stop("the zero stratum totals ", sfmt(sum(out$weight)), " cells but the matrix has ",
         sfmt(n_cells - nrow(mx$obs)), " cells holding no report", call. = FALSE)
  log_msg("    ", sfmt(n_cells - nrow(mx$obs)), " cells hold no report, grouped to ",
          sfmt(nrow(out)), " weighted points")
  out
}

# ── the points the prior is fitted on ────────────────────────────────────────
# Observed cells are binned by autoSquash(), which keeps a count stratum whole
# when it holds fewer points than keep_pts; a FAERS matrix has thousands of
# sparse strata, so expect of the order 1e4 to 1e5 weighted points rather than
# the couple of thousand the vignette's example produces. The zero cells arrive
# already grouped, from the margins.
ebgm_points <- function(mx, n_total) {
  stopifnot(nrow(mx$obs) > 0L, !anyNA(mx$obs$N), !anyNA(mx$obs$E))
  sq <- as.data.table(openEBGM::autoSquash(as.data.frame(mx$obs[, .(N, E)])))
  log_msg("    squashed ", sfmt(nrow(mx$obs)), " observed pairs to ", sfmt(nrow(sq)),
          " weighted points")
  out <- rbindlist(list(sq[, .(N = as.numeric(N), E = as.numeric(E),
                               weight = as.numeric(weight))],
                        ebgm_zero_points(mx, n_total)), use.names = TRUE)
  stopifnot(all(is.finite(out$N)), all(is.finite(out$E)), all(out$E > 0),
            all(is.finite(out$weight)), all(out$weight > 0))
  out
}

# ── the five hyperparameters ─────────────────────────────────────────────────
# hyperEM() is openEBGM's EM (strictly ECM) estimator. It is used rather than
# autoHyper() for three reasons: autoHyper caps itself at 20,000 points and
# declines beyond it; the EM takes the zero cells, which is what identifies the
# low component at all; and with method = "score" it reparameterises to
# log/log/tan internally, so the parameter space is unconstrained and the edges
# that defeat a direct minimization are not edges any more.
#
# N_star is NULL because nothing is truncated: the matrix is complete, so the
# likelihood is the untruncated one (negLLzeroSquash) rather than the one
# conditioned on N >= 1. Excluding the zero cells leaves the low component with
# nothing to fit, and the estimate walks its shape parameter to zero - a point
# mass, not a gamma, and not a prior.
#
# Every starting value is run and they are required to agree, which is the check
# autoHyper performs and the only thing a single run would lose. score_norm is
# logged: it is the Euclidean norm of the score vector at the estimate and
# should be near zero at a maximum, so it is the one number that says whether
# the answer is a maximum rather than wherever the algorithm stopped.
# The mixture likelihood is invariant under exchanging the two components, and
# so are Qn(), ebgm() and quantBisect(), so a fit may return either labeling
# with the same log-likelihood. Left alone that has two consequences: which
# labeling reaches ebgm_prior_<tag>.csv is decided by floating-point noise
# between starting values that agree perfectly, and the spread check below
# compares alpha1 from one labeling against alpha1 from the other and reports
# a disagreement that is not there. Components are therefore ordered by their
# mean, larger first, which is the parameterisation DuMouchel publishes
# (alpha1 = 0.2, beta1 = 0.1, mean 2; alpha2 = 2, beta2 = 4, mean 0.5) and so
# the one a reader will compare the printed prior against.
canon_theta <- function(th) {
  m1 <- th[1] / th[2]; m2 <- th[3] / th[4]
  if (is.finite(m1) && is.finite(m2) && m1 < m2)
    th <- c(th[3], th[4], th[1], th[2], 1 - th[5])
  th
}

fit_ebgm_prior <- function(pts) {
  d <- as.data.frame(pts[, .(N, E, weight)])
  # The convergence test is on an absolute change in the log-likelihood, so the
  # tolerance has to be relative to the size of that likelihood - see
  # EBGM_EM_LL_REL in config.R. It is scaled by the likelihood itself, evaluated
  # at the starting values, and not by the total weight: the zero stratum is
  # over a billion cells but nearly all of them have an expected count near
  # zero and contribute almost nothing, so the total weight runs some eighty
  # times the likelihood - and that ratio grows with the number of one-report
  # drug names, which would leave the criterion several times looser on the
  # 22-year run than on the two-year run it was set from.
  ll0 <- min(vapply(seq_len(nrow(EBGM_THETA_INIT)), function(i)
    openEBGM::negLLzeroSquash(as.numeric(EBGM_THETA_INIT[i, ]),
                              ni = d$N, ei = d$E, wi = d$weight), numeric(1)))
  # not guarded, a single non-finite value here passes through min() and then
  # max(1e-4, NA) is NA, which hyperEM would take as its convergence tolerance
  # and compare against - failing every start with the same opaque message
  if (!is.finite(ll0))
    stop("the log-likelihood is not finite at any of the ", nrow(EBGM_THETA_INIT),
         " starting values, so no convergence tolerance can be derived from it",
         call. = FALSE)
  tol <- max(EBGM_EM_LL_TOL_MIN, EBGM_EM_LL_REL * abs(ll0))
  log_msg("    log-likelihood at the starting values ", sprintf("%.4g", -ll0),
          ": convergence when it moves by less than ", sprintf("%.3g", tol),
          " for ", EBGM_EM_CONSECUTIVE, " iterations")
  fits <- lapply(seq_len(nrow(EBGM_THETA_INIT)), function(i)
    tryCatch(openEBGM::hyperEM(data = d,
                               theta_init_vec = as.numeric(EBGM_THETA_INIT[i, ]),
                               squashed = TRUE, zeroes = TRUE, N_star = NULL,
                               method = "score", LL_tol = tol,
                               consecutive = EBGM_EM_CONSECUTIVE,
                               param_lower = EBGM_PAR_MIN, param_upper = EBGM_PAR_MAX,
                               print_level = 0L, max_iter = EBGM_EM_MAX_ITER),
             error = function(e) list(estimates = rep(NA_real_, 5L), maximum = -Inf,
                                      score_norm = NA_real_, iters = NA_integer_,
                                      message = conditionMessage(e))))
  # the log-likelihood has to be finite too: which.max() over a single NA
  # returns integer(0), which would index the parameter matrix to nothing
  ok <- vapply(fits, function(f) all(is.finite(f$estimates)) && is.finite(f$maximum),
               logical(1))
  # the score vector is summed over the weights, so its norm has no scale of its
  # own: an unnormalized 1e4 looks alarming and is 1e-5 of the total weight.
  # Both are logged, and the second is the one to read.
  w <- sum(pts$weight)
  for (i in which(ok))
    log_msg("    start ", i, ": ", fits[[i]]$iters, " iterations, logL ",
            sprintf("%.4g", fits[[i]]$maximum), ", score norm ",
            sprintf("%.3g", fits[[i]]$score_norm), " (",
            sprintf("%.3g", fits[[i]]$score_norm / w), " per cell)")
  for (i in which(!ok))
    log_msg("    start ", i, " did not return a prior: ", fits[[i]]$message)
  if (!any(ok))
    stop("no starting value produced an MGPS prior: ", fits[[1]]$message, call. = FALSE)

  par <- do.call(rbind, lapply(fits[ok], function(f) canon_theta(as.numeric(f$estimates))))
  # autoHyper's own tolerances for calling two fits the same answer
  spread <- apply(par, 2L, function(x) diff(range(x)))
  if (sum(ok) > 1L && any(spread > c(0.2, 0.2, 0.2, 0.2, 0.05)))
    log_msg("    WARNING: the fits disagree by ",
            paste(sprintf("%.4f", spread), collapse = " "),
            " - the likelihood has more than one optimum and the prior is not ",
            "determined by the data alone")
  best <- which.max(vapply(fits[ok], function(f) f$maximum, numeric(1)))
  out  <- setNames(par[best, ], c("alpha1", "beta1", "alpha2", "beta2", "P"))
  attr(out, "score_norm") <- fits[ok][[best]]$score_norm
  attr(out, "iters")      <- fits[ok][[best]]$iters
  attr(out, "logL")       <- fits[ok][[best]]$maximum
  # Recorded because a prior resting on one surviving start is a different
  # claim from one three starts agreed on, and the five numbers alone do not
  # say which it is. The tolerance goes with them: it is the criterion the fit
  # stopped on, and it is derived at run time rather than set in config.R.
  attr(out, "starts")     <- nrow(EBGM_THETA_INIT)
  attr(out, "converged")  <- sum(ok)
  attr(out, "ll_tol")     <- tol
  if (sum(ok) < 2L)
    log_msg("    NOTE: only ", sum(ok), " of ", nrow(EBGM_THETA_INIT),
            " starting values returned a prior, so the agreement check could ",
            "not be made. tools/refit_ebgm_prior.R refits from more starts ",
            "without rebuilding the matrix.")
  out
}

# A prior the estimator is content with can still be one no reader would accept.
# exploreHypers() rejects a gamma parameter above param_limit = 100; a gamma
# whose shape has reached the floor is a point mass rather than a gamma, and a
# near-point-mass prior at RR = 2 would score every pair at EBGM 2, failing MGPS
# everywhere and reproducing the very result this correction exists to undo -
# while passing every arithmetic check, because the arithmetic would be right.
check_ebgm_prior <- function(theta) {
  bad <- c(
    if (any(theta[1:4] >= EBGM_PAR_MAX * 0.999)) "a gamma parameter is at the upper limit",
    # on the shapes, not on all four: a shape going to zero is the point mass,
    # whereas a small rate only makes a component diffuse, which is legitimate
    if (any(theta[c(1L, 3L)] < EBGM_PAR_DEGENERATE))
      paste0("a gamma shape parameter is below ", EBGM_PAR_DEGENERATE,
             " - that component is a point mass, not a gamma"),
    if (theta[["P"]] <= 0.001 || theta[["P"]] >= 0.999)
      "the mixture fraction has collapsed onto one component")
  if (length(bad))
    stop("the fitted MGPS prior is not usable: ", paste(bad, collapse = "; "), ".\n  ",
         paste(sprintf("%s=%.5f", names(theta), theta), collapse = "  "),
         "\nThe fit is seeded and deterministic, so running again returns this ",
         "same prior: the mixture is not identified on this matrix and the ",
         "reason has to be found rather than retried.", call. = FALSE)
  invisible(theta)
}

# openEBGM rewrites the hyperparameters before scoring with them, in
# .reinInTheta(): each is held inside min_theta/max_theta, and then alpha is
# held within half to twice its own beta. The ratio rule applies however wide
# the bounds are and cannot be switched off, so it is applied here as well -
# otherwise ebgm_prior_<set>.csv would report a prior that was never used.
# (faers::phv_ebgm scored with private copies of these functions that have the
# rule removed, which is why the published run was never subject to it.)
rein_in_theta <- function(theta, min_theta = EBGM_MIN_THETA, max_theta = EBGM_MAX_THETA) {
  th <- pmin(pmax(as.numeric(theta), min_theta), max_theta)
  th[1] <- min(th[1], 2 * th[2]); th[3] <- min(th[3], 2 * th[4])
  th[1] <- max(th[1], 0.5 * th[2]); th[3] <- max(th[3], 0.5 * th[4])
  setNames(th, names(theta))
}

# EBGM and its 90% credible interval for one or more tables, with the prior held
# fixed. digits = 6 so the four decimals the CSVs carry are not pre-rounded.
# E is the expected count phv_signal() already returned for the same table, so
# the two cannot drift apart.
# Whether the installed openEBGM clamps theta before scoring with it. The
# min_theta/max_theta arguments and the .reinInTheta rule behind them arrived
# in a later release; where they are absent the fitted prior is used as it is,
# and where they are present they are widened so they cannot bind, because
# their 0.25 floor sits above what a real fit returns.
ebgm_clamps <- function()
  all(c("min_theta", "max_theta") %chin% names(formals(openEBGM::Qn)))

score_ebgm <- function(theta, a, E) {
  theta <- unname(as.numeric(theta))
  bd <- if (ebgm_clamps()) list(min_theta = EBGM_MIN_THETA, max_theta = EBGM_MAX_THETA) else list()
  qn <- do.call(openEBGM::Qn, c(list(theta_hat = theta, N = a, E = E), bd))
  eb <- do.call(openEBGM::ebgm,
                c(list(theta_hat = theta, N = a, E = E, qn = qn, digits = 6), bd))
  qb <- function(p) do.call(openEBGM::quantBisect,
                            c(list(p, theta_hat = theta, N = a, E = E, qn = qn, digits = 6), bd))
  lo <- qb(EBGM_PCT_LOW); hi <- qb(EBGM_PCT_HIGH)
  # The percentiles are checked for scale in config.R, where they are set. This
  # only catches a bound returned on the wrong side of the point estimate.
  stopifnot(all(lo > 0), all(lo <= eb + 1e-6), all(hi >= eb - 1e-6))
  data.table(ebgm = eb, ebgm_ci_low = lo, ebgm_ci_high = hi)
}

# The prior in force. Set by the step that builds a reference set, before any row
# is scored against that set; every row then uses these five numbers, so a pair's
# shrinkage comes from the background and never from the pair itself.
EBGM_PRIOR <- NULL

set_ebgm_prior <- function(theta) {
  stopifnot(is.numeric(theta), length(theta) == 5L, all(is.finite(theta)))
  EBGM_PRIOR <<- theta
  log_msg("  MGPS prior: ", paste(sprintf("%s=%.5f", names(theta), theta), collapse = "  "))
  invisible(theta)
}

# Fit the prior for one reference set and make it the one in force. `tag` names
# the set ("full", "mono") and its output files.
#
# The weighted points are what is cached, not the prior: they are a few tens of
# thousands of rows where the matrix is over a billion cells, so a fit that has
# to be repeated costs seconds rather than the whole build. The cache records
# the reference set it was built from, because neither the study period nor the
# set size is part of the cache key.
use_ebgm_prior <- function(obj, norm_lookup, n_total, tag) {
  if (!"ebgm" %chin% PHV_METHODS) {
    log_msg("  MGPS is not in PHV_METHODS: prior not fitted")
    return(invisible(NULL))
  }
  need_faers()
  if (!requireNamespace("openEBGM", quietly = TRUE))
    stop("MGPS needs the openEBGM package:\n  install.packages(\"openEBGM\")\n",
         "It is checked here rather than when it is first called, which is after ",
         "the drug-event matrix has been built.", call. = FALSE)
  if (!"hyperEM" %chin% getNamespaceExports("openEBGM"))
    stop("the installed openEBGM (", as.character(utils::packageVersion("openEBGM")),
         ") has no hyperEM(), so the prior cannot be fitted with the zero cells ",
         "included. Upgrade openEBGM.", call. = FALSE)

  stamp <- paste0(START_YEAR, "-", END_YEAR, "/", n_total)
  pts <- cached(paste0("ebgm_points_", tag, ".rds"), {
    log_msg("  MGPS: building the ", tag, " drug-event matrix over ", sfmt(n_total), " reports")
    out <- ebgm_points(ebgm_pair_counts(obj, norm_lookup, n_total), n_total)
    attr(out, "ebgm_stamp") <- stamp
    out
  })
  if (!identical(attr(pts, "ebgm_stamp"), stamp))
    stop("ebgm_points_", tag, ".rds was built from ", attr(pts, "ebgm_stamp"),
         " and this run is ", stamp, ". Delete ", res_path(paste0("ebgm_points_", tag, ".rds")),
         " and run again. Do not set FORCE_REBUILD for this: it would also ",
         "discard every cached MCMC result.", call. = FALSE)

  # The prior is cached as well as the points. Step 04 scores the same reference
  # set, and refitting there would overwrite ebgm_prior_<tag>.csv with a second
  # fit - so the table describing the prior could describe a different fit from
  # the one the published EBGM columns were computed under.
  raw <- cached(paste0("ebgm_theta_", tag, ".rds"), {
    log_msg("  MGPS: fitting the ", tag, " prior on ", sfmt(nrow(pts)),
            " weighted points (", sfmt(sum(pts$weight)), " cells)")
    # hyperEM's maximization step replaces a root that escapes its bounds with a
    # random draw inside them, so the fit is not deterministic unless it is
    # seeded here. Unseeded, the prior - and every EBGM05 and consensus call
    # resting on it - would differ from one run to the next.
    set.seed(PHV_SEED)
    out <- check_ebgm_prior(fit_ebgm_prior(pts))
    attr(out, "ebgm_stamp") <- stamp
    out
  })
  if (!identical(attr(raw, "ebgm_stamp"), stamp))
    stop("ebgm_theta_", tag, ".rds was fitted on ", attr(raw, "ebgm_stamp"),
         " and this run is ", stamp, ". Delete ",
         res_path(paste0("ebgm_theta_", tag, ".rds")), " and run again.", call. = FALSE)
  # only where the installed openEBGM would apply it itself
  theta <- if (ebgm_clamps()) rein_in_theta(raw) else raw
  moved <- as.numeric(raw) != as.numeric(theta)
  # The ratio rule binds by construction on a prior whose high component mean
  # exceeds 2, which every prior this study fits does. Where it fires it pulls
  # EBGM05 down hardest on the low-count rows - the ones where the consensus
  # arm can flip - so a run that hits it is not comparable with the published
  # numbers. openEBGM 0.9.1, the version used, has no such rule; stop rather
  # than publish figures that differ for a reason nobody would look for.
  if (any(moved))
    stop("openEBGM ", as.character(utils::packageVersion("openEBGM")),
         " applies a ratio rule that moved the fitted prior: ",
         paste(sprintf("%s %.5f -> %.5f", names(raw)[moved], as.numeric(raw)[moved],
                       as.numeric(theta)[moved]), collapse = "  "),
         ".\nThe published results were produced with openEBGM 0.9.1, which has no ",
         "such rule, and the two are not comparable on low-count rows.\n",
         '  remotes::install_version("openEBGM", "0.9.1")', call. = FALSE)
  set_ebgm_prior(setNames(as.numeric(theta), names(raw)))
  log_msg("  MGPS: score norm ", sprintf("%.3g", attr(raw, "score_norm")), " after ",
          attr(raw, "iters"), " iterations")
  # both are written: `fitted` is what the likelihood returned, `used` is what
  # scored every row. The check above stops the run if any parameter moved, so
  # in a run that completes they are equal; the two columns exist so that a
  # reader can see that rather than take it on trust.
  write_res(data.table(reference_set = tag, reports = n_total,
                       weighted_points = nrow(pts), cells = sum(pts$weight),
                       log_likelihood = attr(raw, "logL"),
                       score_norm = attr(raw, "score_norm"),
                       iterations = attr(raw, "iters"),
                       ll_tol = attr(raw, "ll_tol"),
                       starts_attempted = attr(raw, "starts"),
                       starts_converged = attr(raw, "converged"),
                       parameter = names(raw), fitted = unname(as.numeric(raw)),
                       used = unname(as.numeric(theta))),
            paste0("ebgm_prior_", tag, ".csv"))
  invisible(EBGM_PRIOR)
}
