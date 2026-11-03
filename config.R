# ==============================================================================
# config.R - the only file you should need to edit.
#
# Paths can also be set through environment variables (useful on a server):
#   FAERS_DATA_DIR, FAERS_RESULTS_DIR, FAERS_CACHE_DIR
# ==============================================================================

# Folder holding the raw inputs:
#   <DATA_DIR>/FAERS_ASCII/         quarterly FAERS ASCII zips (faers_ascii_2004q1.zip, ...)
#   <DATA_DIR>/MedDRA_29_0_English/ MedDRA v29.0 distribution (licensed, not shared)
DATA_DIR <- Sys.getenv("FAERS_DATA_DIR", unset = "~/faers_data")

# Where all tables, intermediate .rds files and figures are written.
# A relative path here is read against the working directory, so running from
# the repository with FAERS_RESULTS_DIR=my_run writes beside the code, not
# where a folder of that name already sits elsewhere. Give an absolute path to
# reuse an existing run: a wrong one is not an error, it is a parse of all 88
# quarters into a new empty folder.
RES <- Sys.getenv("FAERS_RESULTS_DIR", unset = "results")

# Per-year standardized FAERS objects (large; safe to delete after step 01).
CACHE_DIR <- Sys.getenv("FAERS_CACHE_DIR", unset = file.path(RES, "cache"))

FAERS_DIR  <- file.path(path.expand(DATA_DIR), "FAERS_ASCII")
MEDDRA_DIR <- file.path(path.expand(DATA_DIR), "MedDRA_29_0_English")
FIG_DIR    <- file.path(RES, "figures")

# Study period (manuscript: Q1 2004 - Q4 2025). Edit these two lines to run a
# shorter period - for example 2024L and 2025L to watch the pipeline run end to
# end on two years.
START_YEAR <- 2004L
END_YEAR   <- 2025L

# Shell overrides, so a short run needs no edit at all and leaves nothing to
# change back:  FAERS_START_YEAR=2024 FAERS_END_YEAR=2025 Rscript run_all.R
# Nothing to edit below.
if (nzchar(Sys.getenv("FAERS_START_YEAR"))) START_YEAR <- as.integer(Sys.getenv("FAERS_START_YEAR"))
if (nzchar(Sys.getenv("FAERS_END_YEAR")))   END_YEAR   <- as.integer(Sys.getenv("FAERS_END_YEAR"))

# Bayesian signal detection (faers::phv_signal). 10,000 MCMC draws with a fixed
# seed so the BCPNN estimates are reproducible run to run. EBGM is not drawn:
# it is excluded from phv_signal below and scored deterministically in R/ebgm.R,
# so PHV_N_MCMC does not reach it.
PHV_N_MCMC  <- 10000L
PHV_SEED    <- 1L
PHV_METHODS <- c("ror", "prr", "chisq", "bcpnn_mcmc", "ebgm")

# MGPS settings (R/ebgm.R). EBGM05 is the 5th percentile of the posterior, as
# Table 2 states, so the 90% credible interval runs from the 5th to the 95th.
# openEBGM takes these as percentages, not proportions: a proportion passed by
# mistake returns a lower bound that looks entirely plausible, so the scale is
# checked here rather than left to be noticed in the results.
EBGM_PCT_LOW    <- 5
EBGM_PCT_HIGH   <- 95
stopifnot(EBGM_PCT_LOW >= 1, EBGM_PCT_HIGH <= 99, EBGM_PCT_LOW < EBGM_PCT_HIGH)

# Blocks of reports the drug-event matrix is built in, to keep the
# report x drug x PT expansion inside memory. A pair appearing in several blocks
# is a separate row until the counts are added, so fewer blocks means less to
# add up: raise this only if the build runs short of memory.
EBGM_CHUNKS     <- 8L

# The cells of the matrix that hold no report are the bulk of it - over a
# billion - and the prior is fitted with them, because without them the low
# component of the mixture has nothing to fit and its shape parameter goes to
# zero. They are never enumerated: cells sharing a pair of margin counts share
# one expected count, so the stratum is reconstructed from the two margins.
# EBGM_MARGIN_BINS groups each margin (exactly, where there are fewer distinct
# counts than this); EBGM_ZERO_BINS groups the resulting expected counts.
#
# EBGM_MARGIN_BINS is set above the number of distinct margin counts the
# database holds - about 1,100 and 850 over two years, 4,302 and 3,734 over 22 -
# so that the margins are not grouped and every cell keeps its own expected
# count. R/ebgm.R logs which of the two it did, and crossing this limit changes
# the fitted prior rather than merely blurring it, so check that line before
# using the numbers from a larger data cut. The cost is the size of the cross
# product: a few thousand distinct counts a side is a few million rows.
#
# The zero cells are then grouped by expected count into EBGM_ZERO_BINS points
# regardless, so the stratum enters the fit squashed, as the observed cells do.
# What EBGM_MARGIN_BINS controls is whether the expected counts going into that
# squashing are each cell's own or a class average.
EBGM_MARGIN_BINS <- 5000L
EBGM_ZERO_BINS   <- 500L

# Starting values for the hyperparameter fit. Not openEBGM's vignette values,
# which cannot be used here: openEBGM forms the mixture as log(p * NB1 +
# (1 - p) * NB2) rather than in log space, so a point whose density falls below
# the smallest double rounds to zero and its log is -Inf. Each of the vignette's
# three starts puts at least one real cell there - a drug reported almost
# entirely with one term, such as a gadolinium contrast agent with nephrogenic
# systemic fibrosis - and the fit cannot begin. The three below leave every
# squashed point of the unrestricted set representable, that being the larger
# and the sharper of the two, while staying spread across the space so the
# agreement check still has something to check; the guard in R/ebgm.R tests the
# monopharmacy set the same way when it is built, and stops the run rather than
# fitting if it ever fails. They are starting values and nothing more: all
# three converge to the same prior.
EBGM_THETA_INIT <- data.frame(alpha1 = c(0.05, 0.20, 0.30), beta1 = c(0.005, 0.020, 0.030),
                              alpha2 = c(0.5,  2.0,  1.5),  beta2 = c(0.5,   2.0,   1.5),
                              p      = c(0.05, 0.20, 0.30))

# The range a gamma parameter of the prior may take, passed to hyperEM() as
# param_lower and param_upper. openEBGM's exploreHypers refuses anything above
# param_limit = 100, and hyperEM's own floor is 1e-5.
EBGM_PAR_MIN    <- 1e-5
EBGM_PAR_MAX    <- 100

# The separate threshold at which R/ebgm.R calls a fitted component degenerate
# and refuses to publish it. It is NOT the optimizer's bound and must sit well
# above it: hyperEM replaces a root that escapes its bounds with a draw strictly
# inside them, so it never leaves a parameter sitting on the bound and a test
# against the bound itself would never fire. A shape parameter of 1e-3 is a
# point mass for every practical purpose, and a prior like that scores every
# pair at the prior mean - which is the failure this correction exists to undo,
# and it would pass every arithmetic check because the arithmetic would be
# right. DuMouchel's own prior has alpha1 = 0.2, so 0.01 is two orders of
# magnitude clear of a real fit.
EBGM_PAR_DEGENERATE <- 0.01

# hyperEM() stops when the log-likelihood changes by less than LL_tol for
# EBGM_EM_CONSECUTIVE iterations in a row. That test is on an absolute change,
# and this likelihood is weighted by every cell of the matrix - of order 1e9 -
# so a fixed 1e-4 is unreachable and the fit would always run to max_iter. The
# tolerance is therefore set relative to the log-likelihood at the starting
# values, not to the total weight: at 1e-10 of it the hyperparameters are
# settled to about five decimal places. R/ebgm.R sets out why the total weight
# is the wrong scale to use.
EBGM_EM_LL_REL     <- 1e-10
EBGM_EM_LL_TOL_MIN <- 1e-4
EBGM_EM_CONSECUTIVE <- 100L
EBGM_EM_MAX_ITER    <- 5000L

EBGM_MIN_THETA  <- rep(1e-6, 5)
EBGM_MAX_THETA  <- c(1e4, 1e4, 1e4, 1e4, 1 - 1e-6)

# Set TRUE to recompute everything: the quarterly parse and the deduplication
# in step 01 as well as the cached signal results (the MCMC steps are slow).
# Route 3 has no quarterly files to parse, so do not set it there. Delete the
# cached signal files instead - ssri_targeted_signals.rds, perdrug_pt_*.rds,
# group_pooled_*.rds, ebgm_points_*.rds, ebgm_theta_*.rds - and on no account
# master_dedup.rds, which is 3.17 GB of the download and cannot be rebuilt
# without the quarterly files. master_mono.rds is the other 1.18 GB, but step
# 01 does rebuild that one from master_dedup.rds whenever master_mono_stamp.txt
# is absent or does not match, writing to a .tmp file and renaming, so the
# download is never left half-replaced. A rebuilt file will not match the
# SHA-256 the README publishes for the downloaded one.
FORCE_REBUILD <- FALSE

# Quick verification mode: frequentist methods only (ROR, PRR, chi-squared),
# no MCMC, for checking a change end to end. Never use for published results.
QUICK_TEST <- Sys.getenv("FAERS_QUICK_TEST", unset = "FALSE") == "TRUE"
if (QUICK_TEST) PHV_METHODS <- c("ror", "prr", "chisq")

# What phv_signal() is asked for: everything except MGPS, which R/ebgm.R scores
# from a prior fitted once per reference set rather than once per row.
PHV_PER_ROW <- setdiff(PHV_METHODS, "ebgm")

# Cached signal files record this key and are rebuilt when it changes, so a
# cache written by different arithmetic is never silently reused. Bump the
# version whenever a change alters the numbers a cached file holds.
#   2  MGPS prior fitted once per reference set, not once per 2 x 2 table;
#      EBGM05 taken at the 5th percentile of the posterior, not the 2.5th
#   3  MGPS prior fitted with the zero cells of the matrix included, by
#      hyperEM(). Fitted without them, one mixture component is unidentified
#      and its shape parameter goes to zero.
#   4  withdrawn before release; no cache was ever written under it
#   5  the drug axis of the MGPS matrix drops every placeholder reporters write
#      in the drug name field (NON_DRUG_ANY), the "unknown drug" ones as well as
#      the "no co-medication" ones: neither is a drug, and UNKNOWN would be one
#      of the highest-count rows of the axis the prior is fitted over. The
#      monopharmacy count keeps the "unknown" ones on purpose, which is why
#      there are two lists. PHV_DEF_HASH does not catch this: the definitions
#      file is unchanged and the edit is in R/ebgm.R, so without this bump the
#      cached points and prior from an earlier run would be reused and the
#      change would do nothing.
#   6  report_drug_counts() resolves a drug by its active ingredient when the
#      name itself is not in the dictionary, so one drug written two ways is
#      counted once. The edit is in R/groups.R, which PHV_DEF_HASH does not
#      hash, so the bump is what makes the cached counts rebuild.
#   7  report_drug_counts() counts a drug once per report. A line still takes
#      the ingredient recorded against it, but a line with none of its own now
#      takes the one its same-named siblings agree on, which stops one drug
#      written twice with prod_ai filled in only once from counting as two.
#      The sibling grouping is on the name as written, not on the canonical
#      key: drug_key() rewrites a name to its last parenthetical, so
#      "SERTRALINE (TABLETS)" and "FLUOXETINE (TABLETS)" share a key while
#      being different drugs, and grouping on the key would count two SSRIs as
#      one. A prod_ai holding several ingredients - FAERS marks a combination
#      there with a backslash - is not used to resolve a name, since taking its
#      first match would merge two drugs into one. A name itself is read as the
#      dictionary reads it for group membership, with no parsing of separators.
#      The edits are in R/groups.R, which PHV_DEF_HASH does not hash, so the
#      bump is what rebuilds the counts.
#   8  The name-parsing rule tried under 7 (commas, slashes and a list of form
#      words) is withdrawn: it read manufacturers and foreign form text as
#      second drugs. A name is read only as the dictionary reads it. Bumped
#      so that nothing counted under 7 is reused.
#   9  drug_canon() keys a name by its last parenthetical, so two different
#      drugs carrying the same clarifier - "CITALOPRAM (NO PREF. NAME)" and
#      "COCAINE (NO PREF. NAME)", as older reports write them - collapsed into
#      one and made a polypharmacy report look like monopharmacy. The counter
#      now reads the dictionary over the whole name first, as the group filter
#      already did, so the two no longer disagree.
#  10  drug_key() keys a name by its last bracket, which usually holds the
#      ingredient. Older reports put a clarifier there instead - the same one
#      on every line - so every drug on such a report merged into one, and a
#      polypharmacy report counted as monopharmacy. A clarifier is not an
#      ingredient, so a name whose bracket holds one keeps its own form.
#      data/non_drug_terms.csv names them, and gains NO PREF. NAME and NGX.
#  11  Withdrawn in favour of not unwrapping at all. Reporters put anything in
#      that bracket - a form, a strength, a manufacturer, an annotation - and
#      no list keeps up with every strength. A name now keys as it is written,
#      which is what the Methods describes: counted by name, with the study
#      and comparator drugs resolved through the dictionary. For a drug
#      outside the dictionary a brand and its ingredient on one report count
#      as two, which can only move a report out of monopharmacy.
# The study period is part of the key as well: without it, changing the two
# lines above and running into the same results folder would reuse every cached
# signal from the previous period without a word.
# The settings that change a cached number are part of the key as well, so that
# editing one of them in this file and running into the same results folder
# rebuilds what it affects instead of silently reusing it. A config edit is not
# a version bump, and a sensitivity check that reused the cache would return the
# original answer and look like a confirmation.
PHV_CACHE_VERSION <- 11L

# The study definitions live in data/*.csv, not in this file, and editing one
# changes the numbers just as surely as editing a setting here. They are hashed
# into the key so that a definition edit rebuilds what it affects. Without this
# the edit is invisible: every count still agrees with every other count,
# because they all come from the same stale cache.
PHV_DEF_FILES <- c("data/drug_dictionary.csv", "data/pt_definitions.csv",
                   "data/mono_overrides.csv", "data/non_drug_terms.csv",
                   Sys.getenv("FAERS_BLACKLIST", unset = "data/manual_dedup_blacklist.csv"))
PHV_DEF_HASH <- local({
  # sorted by the C locale explicitly: R/setup.R has not run yet, so the
  # ambient collation could otherwise order a custom blacklist path differently
  # on two machines and change the key for no reason.
  f <- PHV_DEF_FILES[file.exists(PHV_DEF_FILES)]
  if (!length(f)) return("nodefs")
  f <- f[order(f, method = "radix")]
  t <- tempfile(); on.exit(unlink(t))
  writeLines(paste(tools::md5sum(f), collapse = ""), t)
  substr(unname(tools::md5sum(t)), 1, 10)
})

# The MGPS prior is cached too, and these settings decide what it converges to.
PHV_EBGM_HASH <- local({
  v <- paste(c(unlist(EBGM_THETA_INIT), EBGM_PAR_MIN, EBGM_PAR_MAX, EBGM_PAR_DEGENERATE,
               EBGM_EM_LL_REL, EBGM_EM_LL_TOL_MIN, EBGM_EM_CONSECUTIVE, EBGM_EM_MAX_ITER,
               EBGM_MIN_THETA, EBGM_MAX_THETA), collapse = "|")
  t <- tempfile(); on.exit(unlink(t)); writeLines(v, t)
  substr(unname(tools::md5sum(t)), 1, 8)
})

PHV_CACHE_KEY     <- c(PHV_METHODS, paste0("v", PHV_CACHE_VERSION),
                       paste0(START_YEAR, "-", END_YEAR),
                       paste0("pct", EBGM_PCT_LOW, "-", EBGM_PCT_HIGH),
                       paste0("bins", EBGM_MARGIN_BINS, "/", EBGM_ZERO_BINS),
                       paste0("chunks", EBGM_CHUNKS),
                       paste0("mcmc", PHV_N_MCMC, "/seed", PHV_SEED),
                       paste0("defs", PHV_DEF_HASH),
                       paste0("ebgm", PHV_EBGM_HASH))

# Optional: write every PS-filtered SSRI report (~160k rows) as an audit file.
WRITE_FULL_SSRI_CASE_DUMP <- FALSE
