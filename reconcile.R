# ==============================================================================
# reconcile.R - cross-checks the numbers the pipeline wrote against each other.
#
# Nothing here recomputes the analysis. It reads the files in results/ and tests
# them against the arithmetic they imply: each contingency table must sum to the
# reference set, each estimate must follow from its four cells, each signal flag
# must follow from its threshold, the monopharmacy rows must be a subset of the
# unrestricted rows, the per-drug rows must partition the pooled row, the case
# listings must hold as many reports as cell a claims, and every report named in
# the removal list as the retained duplicate must itself be absent from it.
#
# Reads  results/disproportionality_all.csv   (run combine_da.R first)
#        results/cases_*.csv, results/mono_cases_*.csv
#        data/manual_dedup_blacklist.csv
# Writes results/reconciliation_report.csv
#
# Usage: Rscript reconcile.R          exit status 1 if any check fails
# ==============================================================================
# bit64 as well as data.table: without it fread reads a long PRIMARYID as a
# double and stringifies it in exponent form, and the check that no removed
# report appears in a case listing would then pass without comparing anything.
suppressPackageStartupMessages({ library(data.table); library(bit64) })

# Folder to check: the first argument if it names one, otherwise
# FAERS_RESULTS_DIR, otherwise results/. Pass "expected_2004_2025" to check the files
# behind the article:  Rscript reconcile.R expected_2004_2025
.args <- commandArgs(trailingOnly = TRUE)
RES  <- Sys.getenv("FAERS_RESULTS_DIR", unset = "results")
if (length(.args) && dir.exists(.args[1])) RES <- .args[1]
DATA <- "data"   # study definitions, always in the repository
TOL  <- 1e-3          # relative tolerance, and an absolute floor of 1e-3: enough
                      # for a CSV written at four decimals or more

logdt <- data.table(check = character(), status = character(), detail = character())
fails <- list()

rec <- function(k, ok, detail = "", rows = NULL) {
  st <- if (length(ok) != 1L || is.na(ok)) "SKIP" else if (isTRUE(ok)) "PASS" else "FAIL"
  logdt <<- rbind(logdt, data.table(check = k, status = st, detail = as.character(detail)[1]))
  if (st == "FAIL" && !is.null(rows) && nrow(rows) > 0L) fails[[k]] <<- head(rows, 15L)
  invisible(NULL)
}

rd <- function(f, dir = RES) {
  p <- file.path(dir, f)
  if (!file.exists(p)) return(NULL)
  x <- fread(p, showProgress = FALSE)
  if ("ids" %in% names(x)) x[, ids := NULL]
  x[]
}

# Equal to three decimals, or to within TOL in relative terms. A missing value
# compares equal: a row below the three-case floor carries no estimate on either
# side, and there are many of them.
#
# NaN does NOT compare equal, though is.na() is true of it. That distinction is
# the difference between a check and a formality: a recomputation that overflows
# or divides by zero returns NaN on every row, and treating that as "missing"
# would have this function report perfect agreement precisely when it had
# computed nothing at all.
agree <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  broken <- is.nan(x) | is.nan(y)
  (is.na(x) | is.na(y) | abs(x - y) <= 1e-3 | abs(x - y) <= TOL * pmax(abs(x), abs(y))) & !broken
}
fmt <- function(x) formatC(x, big.mark = ",", format = "d")

# ── 1. the combined estimate table ───────────────────────────────────────────
da <- rd("disproportionality_all.csv")
if (is.null(da)) stop(file.path(RES, "disproportionality_all.csv"), " not found - run combine_da.R first")
cat("\nread ", file.path(RES, "disproportionality_all.csv"), ": ", nrow(da), " rows\n", sep = "")

da[, `:=`(A = as.numeric(a), B = as.numeric(b), C = as.numeric(c), D = as.numeric(d))]

inc <- da[is.na(A) | is.na(B) | is.na(C) | is.na(D)]
rec("contingency cells present on every row", nrow(inc) == 0L,
    sprintf("%d of %d rows incomplete", nrow(inc), nrow(da)), inc)

da[, N := A + B + C + D]
nn <- da[, .(k = uniqueN(N), vals = paste(fmt(sort(unique(N))), collapse = " / ")), by = analysis]
rec("reference set size identical across files", all(nn$k == 1L),
    paste(paste0(nn$analysis, ": ", nn$vals), collapse = " | "),
    da[, .N, by = .(analysis, N)])

if ("n_group" %in% names(da)) {
  bb <- da[!is.na(n_group) & A + B != as.numeric(n_group)]
  rec("a + b = group size", nrow(bb) == 0L, sprintf("%d rows disagree", nrow(bb)),
      bb[, .(analysis, drug_group, drug, outcome, a, b, n_group)])
} else rec("a + b = group size", NA, "n_group not carried in the file")

# a + c is the number of reports with that PT in the reference set, so it cannot
# depend on which drug the row is about
ptt <- da[, .(k = uniqueN(A + C), vals = paste(fmt(sort(unique(A + C))), collapse = " / ")),
          by = .(analysis, outcome)]
rec("PT total (a + c) does not depend on the drug", all(ptt$k == 1L),
    sprintf("%d of %d outcome blocks disagree", sum(ptt$k != 1L), nrow(ptt)), ptt[k != 1L])

key <- c("analysis", "drug_group", "drug", "outcome")
cf <- da[, .(k = uniqueN(paste(a, b, c, d))), by = key][k > 1L]
rec("rows repeated across source files agree", nrow(cf) == 0L,
    sprintf("%d keys conflict", nrow(cf)), cf)

# ── 2. every estimate recomputed from the four cells ─────────────────────────
e <- da[A > 0 & B > 0 & C > 0 & D > 0]
e[, `:=`(
  ror_x     = (A / B) / (C / D),
  se_ror    = sqrt(1/A + 1/B + 1/C + 1/D),
  prr_x     = (A / (A + B)) / (C / (C + D)),
  se_prr    = sqrt(1/A - 1/(A + B) + 1/C - 1/(C + D)),
  chisq_x   = N * pmax(abs(A * D - B * C) - N / 2, 0)^2 /
              ((A + B) * (C + D) * (A + C) * (B + D)))]
e[, `:=`(ror_lo_x = exp(log(ror_x) - 1.96 * se_ror), ror_hi_x = exp(log(ror_x) + 1.96 * se_ror),
         prr_lo_x = exp(log(prr_x) - 1.96 * se_prr), prr_hi_x = exp(log(prr_x) + 1.96 * se_prr))]
# The BCPNN point estimate follows from the cells alone (Table 2). Its lower
# bound IC025 is the 2.5th percentile of an MCMC posterior, so it cannot be
# rebuilt here and is left to the signal step. The MGPS estimates can be, once
# the hyperparameters are known: see section 9.
e[, ic_x := log2((A + 0.5) / ((A + B) * (A + C) / N + 0.5))]

for (p in list(c("ror", "ror_x"), c("ror_ci_low", "ror_lo_x"), c("ror_ci_high", "ror_hi_x"),
               c("prr", "prr_x"), c("prr_ci_low", "prr_lo_x"), c("prr_ci_high", "prr_hi_x"),
               c("chisq", "chisq_x"), c("ic", "ic_x"))) {
  if (!all(p[1] %in% names(e))) { rec(paste(p[1], "follows from a, b, c, d"), NA, "column absent"); next }
  bad <- e[!agree(get(p[1]), get(p[2]))]
  rec(paste(p[1], "follows from a, b, c, d"), nrow(bad) == 0L,
      sprintf("%d of %d rows differ", nrow(bad), nrow(e)),
      bad[, c(key, "a", "b", "c", "d", p[1], p[2]), with = FALSE])
}

ord <- da[(!is.na(ror)  & !is.na(ror_ci_low)  & (ror  < ror_ci_low  | ror  > ror_ci_high)) |
          (!is.na(prr)  & !is.na(prr_ci_low)  & (prr  < prr_ci_low  | prr  > prr_ci_high)) |
          (!is.na(ic)   & !is.na(ic025)       &  ic025 > ic) |
          (!is.na(ebgm) & !is.na(ebgm05)      &  ebgm05 > ebgm)]
rec("point estimates lie inside their own intervals", nrow(ord) == 0L,
    sprintf("%d rows outside", nrow(ord)),
    ord[, .(analysis, drug_group, drug, outcome, ror, ror_ci_low, ror_ci_high, ic, ic025, ebgm, ebgm05)])

# ── 3. the four signal flags and the consensus ───────────────────────────────
flg <- c("ror_signal", "prr_signal", "ic_signal", "ebgm_signal")
if (all(flg %in% names(da))) {
  da[, cons_x := ror_signal & prr_signal & ic_signal & ebgm_signal]
  bad <- da[!is.na(signal) & !is.na(cons_x) & signal != cons_x]
  rec("consensus is the conjunction of the four flags", nrow(bad) == 0L,
      sprintf("%d rows disagree", nrow(bad)),
      bad[, c(key, flg, "signal"), with = FALSE])

  # the thresholds as add_signals() in R/signals.R states them: the three-case
  # minimum sits on the PRR, so it governs the consensus through the conjunction
  da[, `:=`(ror_x_sig  = !is.na(ror_ci_low) & ror_ci_low > 1,
            prr_x_sig  = !is.na(prr) & prr >= 2 & !is.na(A) & A >= 3,
            ic_x_sig   = !is.na(ic025)  & ic025 > 0,
            ebgm_x_sig = !is.na(ebgm05) & ebgm05 > 2)]
  for (f in flg) { xf <- sub("_signal$", "_x_sig", f); da[A == 0, (xf) := FALSE] }
  for (f in flg) {
    xf  <- sub("_signal$", "_x_sig", f)
    bad <- da[!is.na(get(f)) & !is.na(get(xf)) & get(f) != get(xf)]
    rec(paste(f, "follows from its threshold"), nrow(bad) == 0L,
        sprintf("%d of %d rows disagree", nrow(bad), nrow(da)),
        bad[, c(key, "a", "ror_ci_low", "prr", "chisq", "ic025", "ebgm05", f, xf), with = FALSE])
  }
} else rec("signal flags recomputed", NA, "flag columns absent")

# the three-case minimum must bind on every reported signal
low <- da[signal %in% TRUE & A < 3]
rec("no reported signal rests on fewer than three cases", nrow(low) == 0L,
    sprintf("%d of %d reported signals", nrow(low), sum(da$signal %in% TRUE)),
    low[, .(analysis, drug_group, drug, outcome, a, ror, signal)])

# a PT row with fewer than three cases carries no estimate at all. The composite
# rows are matched on their own label rather than by PT, so they are held out of
# this check; they carry the same floor (R/signals.R, min_a = 3L).
pt_low <- da[!grepl("co-report|^\\(|^\\[", outcome) & A > 0 & A < 3 & !is.na(ror)]
rec("PT rows under three cases carry no estimate", nrow(pt_low) == 0L,
    sprintf("%d rows estimated below the floor", nrow(pt_low)),
    pt_low[, .(analysis, drug_group, drug, outcome, a, ror, ror_ci_low, signal)])

# ── 4. monopharmacy against unrestricted ────────────────────────────────────
kk <- c("drug_group", "drug", "outcome")
u  <- unique(da[analysis == "Unrestricted"], by = kk)[, c(kk, "a"), with = FALSE]
m  <- unique(da[analysis == "Monopharmacy"], by = kk)[, c(kk, "a"), with = FALSE]
setnames(u, "a", "a_unr"); setnames(m, "a", "a_mono")
j <- merge(u, m, by = kk)
rec("monopharmacy cases never exceed unrestricted cases", nrow(j[a_mono > a_unr]) == 0L,
    sprintf("%d of %d outcome-drug pairs exceed", nrow(j[a_mono > a_unr]), nrow(j)),
    j[a_mono > a_unr])

# ── 5. per-drug rows partition the pooled row ───────────────────────────────
pooled <- c("SSRI", "Varenicline")
for (an in c("Unrestricted", "Monopharmacy")) {
  pd <- unique(da[analysis == an & drug_group == "SSRI" & !drug %chin% pooled], by = c("drug", "outcome"))
  po <- unique(da[analysis == an & drug_group == "SSRI" &  drug == "SSRI"], by = "outcome")
  if (!nrow(pd) || !nrow(po)) { rec(paste(an, "per-drug sums to pooled"), NA, "rows absent"); next }
  s <- merge(pd[, .(sum_drug = sum(A), drugs = .N), by = outcome],
             po[, .(outcome, a_pooled = A)], by = "outcome")
  bad <- s[sum_drug != a_pooled]
  rec(paste(an, "per-drug cases sum to the pooled count"), nrow(bad) == 0L,
      sprintf("%d of %d outcomes disagree", nrow(bad), nrow(s)), bad)
}

# A pipeline run that stops early leaves one of the twelve source files unwritten.
# combine_da.R then writes a shorter table, and the Supplementary sheet empties a
# whole block of cells that reads exactly like the n < 3 suppression - so nothing
# downstream notices. Every outcome-drug pair must appear under both analyses.
ap <- unique(da[, .(outcome, drug_group, drug, analysis)])
cv <- dcast(ap, outcome + drug_group + drug ~ analysis, value.var = "analysis",
            fun.aggregate = length)
if (all(c("Unrestricted", "Monopharmacy") %in% names(cv))) {
  bad <- cv[Unrestricted == 0L | Monopharmacy == 0L]
  rec("every outcome-drug pair is present under both analyses", nrow(bad) == 0L,
      sprintf("%d of %d pairs appear under only one analysis", nrow(bad), nrow(cv)), bad)
} else {
  rec("every outcome-drug pair is present under both analyses", NA, "one analysis absent")
}

# the homicide-suicide composite can only be as large as its smaller component
hsn <- grep("co-report|^\\[H-S\\]", unique(da$outcome), value = TRUE, ignore.case = TRUE)
if (length(hsn) == 1L) {
  comp <- dcast(da[outcome %chin% c(hsn, "Homicide", "Completed Suicide")],
                analysis + drug_group + drug ~ outcome, value.var = "A", fun.aggregate = max)
  setnames(comp, hsn, "hs")
  if (all(c("Homicide", "Completed Suicide") %in% names(comp))) {
    comp[, lim := pmin(Homicide, `Completed Suicide`)]
    bad <- comp[is.finite(hs) & is.finite(lim) & hs > lim]
    rec("H-S composite does not exceed either component", nrow(bad) == 0L,
        sprintf("%d rows exceed", nrow(bad)), bad)
  } else rec("H-S composite does not exceed either component", NA, "component rows absent")
} else rec("H-S composite does not exceed either component", NA, "composite row not found")

# ── 6. case listings ────────────────────────────────────────────────────────
# The list actually in force, so the checks below apply to a reviewer's own list
# rather than always to ours.
BLACKLIST_PATH <- Sys.getenv("FAERS_BLACKLIST", unset = file.path(DATA, "manual_dedup_blacklist.csv"))
bl <- if (file.exists(BLACKLIST_PATH)) fread(BLACKLIST_PATH, showProgress = FALSE)[] else NULL
blids <- if (!is.null(bl)) as.character(bl$primaryid) else character()

cls <- list.files(RES, pattern = "^(mono_)?(cases|caselist)_.*\\.csv$|^SI_DATA_.*CASE_LIST\\.csv$", full.names = TRUE)
# the quality-assurance listings belong to step 04, not to a study outcome
cls <- cls[!grepl("ingrowing_nail", basename(cls))]
if (!length(cls)) rec("case listings", NA, "no case listing files yet") else {
  # file name to the outcome it lists, so each listing can be held against cell a
  norm_label <- function(x) gsub(" +", " ", trimws(gsub("[^a-z0-9]+", " ", tolower(x))))
  outcome_of <- function(f) {
    s <- sub("\\.csv$", "", f)
    if (grepl("^SI_DATA_HS_CASE_LIST$", s)) return(if (length(hsn)) hsn[1] else NA_character_)
    if (grepl("^SI_DATA_H_CASE_LIST$", s))  return("Homicide")
    s <- sub("^(mono_)?(cases|caselist)_", "", s)
    s <- sub("_(integrated|with_indications|detail|full)$", "", s)
    if (grepl("^(hs|homicide_suicide)$", s)) return(if (length(hsn)) hsn[1] else NA_character_)
    o <- unique(da$outcome)
    s <- norm_label(s)
    hit <- o[norm_label(o) == s]
    # a listing may be named for the outcome without its trailing noun,
    # as violence_related is for "Violence-Related Symptom"
    if (!length(hit)) hit <- o[startsWith(norm_label(o), paste0(s, " "))]
    if (length(hit) == 1L) hit[1] else NA_character_
  }
  tally <- rbindlist(lapply(cls, function(p) {
    x  <- fread(p, showProgress = FALSE)
    pc <- intersect(c("primaryid", "PRIMARYID", "primary_id"), names(x))[1]
    cc <- intersect(c("drug_group", "group_label"), names(x))[1]
    id <- if (is.na(pc)) NULL else as.character(x[[pc]])
    ssri <- if (is.na(pc)) NA_integer_
            else if (is.na(cc)) uniqueN(id)
            else uniqueN(id[!grepl("varenicline", tolower(as.character(x[[cc]])))])
    data.table(file     = basename(p),
               analysis = if (grepl("^mono_", basename(p))) "Monopharmacy" else "Unrestricted",
               outcome  = outcome_of(basename(p)),
               rows     = nrow(x),
               ids      = if (is.na(pc)) NA_integer_ else uniqueN(id),
               ssri     = ssri,
               listed   = if (is.na(pc)) NA_integer_ else sum(id %chin% blids))
  }))
  rec("every case listing maps to an outcome", !any(is.na(tally$outcome)),
      sprintf("%d of %d listings could not be matched to an outcome",
              sum(is.na(tally$outcome)), nrow(tally)), tally[is.na(outcome)])
  rec("case listings hold one row per report", all(tally$rows == tally$ids, na.rm = TRUE),
      sprintf("%d files repeat a report", sum(tally$rows != tally$ids, na.rm = TRUE)),
      tally[rows != ids])
  rec("no removed report appears in a case listing", all(tally$listed == 0L, na.rm = TRUE),
      sprintf("%d removed reports found", sum(tally$listed, na.rm = TRUE)), tally[listed > 0L])

  pool <- unique(da[drug_group == "SSRI" & drug == "SSRI", .(analysis, outcome, a = A)],
                 by = c("analysis", "outcome"))
  tc <- merge(tally[!is.na(outcome)], pool, by = c("analysis", "outcome"))
  if (nrow(tc)) {
    bad <- tc[ssri != a]
    rec("case listings hold as many reports as cell a", nrow(bad) == 0L,
        sprintf("%d of %d matched listings disagree", nrow(bad), nrow(tc)),
        bad[, .(file, analysis, outcome, ssri, a)])
    cat("\ncase listings held against cell a\n")
    print(tc[, .(file, outcome, listing = ssri, a, ok = ssri == a)])
  } else rec("case listings hold as many reports as cell a", NA, "no listing matched an outcome")
  cat("\nall case listings\n"); print(tally)
}

# ── 7. Table 1 ──────────────────────────────────────────────────────────────
# Table 1 is built from the case listings through its own lookups, so none of
# the checks above touch it. Each block is held against the case count it reports.
for (tf in list.files(RES, pattern = "^table1_.*\\.csv$", full.names = TRUE)) {
  t1 <- fread(tf, showProgress = FALSE, colClasses = "character")
  nm <- basename(tf)
  if (!all(c("Characteristic", "n", "%") %in% names(t1))) {
    rec(paste(nm, "has the expected columns"), FALSE, paste(names(t1), collapse = ", ")); next
  }
  # counts only: rows such as "Mean age, years (SD)" carry a mean and an SD in
  # the same two columns
  t1[, `:=`(ch = trimws(Characteristic),
            N = suppressWarnings(as.integer(ifelse(grepl("^[0-9]+$", trimws(n)), trimws(n), NA))),
            P = suppressWarnings(as.numeric(`%`)))]
  t1[grepl("Median|Mean", ch), `:=`(N = NA_integer_, P = NA_real_)]
  N1 <- t1$N[1]
  rec(paste(nm, "table 1 N matches cell a"),
      !is.na(N1) && N1 %in% da[drug_group == "SSRI" & drug == "SSRI", A],
      sprintf("table says %s", N1))

  bad <- t1[!is.na(N) & !is.na(P) & abs(P - round(100 * N / N1, 1)) > 0.051]
  rec(paste(nm, "every percentage is n / N"), nrow(bad) == 0L,
      sprintf("%d rows disagree", nrow(bad)), bad[, .(ch, N, P, expected = round(100 * N / N1, 1))])

  # The three mutually exclusive labels must account for every case. Each is
  # counted directly in 05_case_tables.R, none derived as the residual, so this
  # can actually fail: a case that falls through every label is a case the
  # monopharmacy rule did not classify.
  mp <- t1[grepl("^\\s*((Mono|Poly)pharmacy|No drug named)", ch), sum(N)]
  rec(paste(nm, "monopharmacy + polypharmacy + no drug = N"), length(mp) == 1L && mp == N1,
      sprintf("%s vs %s", mp, N1))

  # Rows under a heading, stopping at the next heading. A heading is a row with
  # no count of its own; the file's indentation does not survive being read back.
  is_head <- !grepl("^[0-9]+$", trimws(t1$n))
  blk <- function(start) {
    i <- match(start, t1$ch)
    if (is.na(i)) return(t1[0])
    k <- integer(0); r <- i + 1L
    while (r <= nrow(t1) && !is_head[r]) { k <- c(k, r); r <- r + 1L }
    t1[k]
  }
  dr <- blk("SSRI drug")
  rec(paste(nm, "SSRI drug rows sum to N"), nrow(dr) == 0L || sum(dr$N) == N1,
      sprintf("%s vs %s", sum(dr$N), N1), dr[, .(ch, N)])

  ct <- blk("Top reporting countries")
  ct <- ct[!grepl("Other/not recorded", ch)]
  rec(paste(nm, "countries are in descending order"),
      nrow(ct) < 2L || !is.unsorted(rev(ct$N)),
      paste(ct$N, collapse = " "), ct[, .(ch, N)])

  cd <- blk("Top 10 co-drugs")
  spl <- if (nrow(cd) < 2L) character() else {
    u <- toupper(cd$ch)
    unlist(lapply(u, function(a) u[u != a & startsWith(u, paste0(a, " "))]))
  }
  rec(paste(nm, "no co-drug is a salt form of another"), length(spl) == 0L,
      paste(spl, collapse = "; "))
}

# ── 8. the removal list ─────────────────────────────────────────────────────
if (is.null(bl)) rec("removal list reads", FALSE, paste(BLACKLIST_PATH, "not found")) else {
  rec("removal list identifiers are unique", !any(duplicated(blids)),
      sprintf("%d repeated", sum(duplicated(blids))))
  rsn <- bl[, .N, by = reason][order(-N)]
  rec("removal reasons are the expected three",
      identical(sort(unique(as.character(bl$reason))),
                c("duplicate_report", "implausible_date", "serial_reporter")),
      paste(paste0(rsn$reason, " ", rsn$N), collapse = " | "))
  rec("every row carries a description",
      all(nzchar(trimws(as.character(bl$case_description)))),
      sprintf("%d blank", sum(!nzchar(trimws(as.character(bl$case_description))))))

  # the retained report comes from its own column where the list carries one;
  # otherwise from the wording, taking the identifier that follows it rather than
  # the first number in the sentence, because dates are numbers too
  dup <- bl$reason == "duplicate_report"
  if ("duplicate_of" %in% names(bl)) {
    bl[, named := as.character(duplicate_of)]
    bl[!is.na(named) & !nzchar(trimws(named)), named := NA_character_]
  } else {
    bl[, named := NA_character_]
    pt  <- "(?i)(duplicate of|retained|retains|keeps?|kept)[^0-9]{0,20}([0-9]{6,10})"
    got <- rep(NA_character_, sum(dup))
    hit <- regexpr(pt, bl$case_description[dup], perl = TRUE) > 0
    if (any(hit)) got[hit] <- sub(".*?([0-9]{6,10})$", "\\1",
                                  regmatches(bl$case_description[dup][hit],
                                             regexpr(pt, bl$case_description[dup][hit], perl = TRUE)))
    bl[dup, named := got]
  }
  n_named <- sum(!is.na(bl$named) & dup)
  rec("retained report named on duplicate rows", n_named > 0L,
      sprintf("%d of %d duplicate rows, and %d of %d rows overall",
              n_named, sum(dup), sum(!is.na(bl$named)), nrow(bl)))

  dead <- bl[!is.na(named) & named %chin% blids]
  rec("every named retained report is itself retained", nrow(dead) == 0L,
      sprintf("%d notes point at a report that was also removed", nrow(dead)),
      dead[, .(primaryid, reason, named, case_description)])

  rec("a counterpart is recorded only where one was retained",
      all(is.na(bl$named[bl$reason == "implausible_date"])),
      sprintf("%d rows removed on their own grounds name one",
              sum(!is.na(bl$named[bl$reason == "implausible_date"]))))
  cat("\nremoval list\n"); print(rsn)
}

# ── 9. MGPS, recomputed from the published hyperparameters ──────────────────
# EBGM and EBGM05 are determined by the cells and the five hyperparameters of
# the gamma mixture, so once ebgm_prior_<set>.csv is in the folder they can be
# held against the cells like the other three algorithms. Transcribed from
# openEBGM 0.9.1 (f_posteriorCalcs.R); base R only, so this file still needs
# nothing but data.table.
#
# This is the check that was missing when MGPS was fitted one row at a time:
# the only MGPS test then was that the lower bound sat below the estimate, and
# an EBGM of 0.4257 with a bound of 0.10 satisfied it.
ebgm_qn <- function(th, N, E) {
  p1 <- pmin(pmax(th[2] / (th[2] + E), 1e-10), 1 - 1e-10)
  p2 <- pmin(pmax(th[4] / (th[4] + E), 1e-10), 1 - 1e-10)
  num <- log(th[5]) + dnbinom(N, size = th[1], prob = p1, log = TRUE)
  den <- exp(num) + exp(log(1 - th[5]) + dnbinom(N, size = th[3], prob = p2, log = TRUE))
  exp(num - log(den))
}
ebgm_point <- function(th, N, E, qn)
  exp(qn * (digamma(th[1] + N) - log(th[2] + E)) +
      (1 - qn) * (digamma(th[3] + N) - log(th[4] + E)))
# bisection on the posterior CDF, as quantBisect does
ebgm_quant <- function(th, N, E, qn, pct) {
  cdf <- function(x) exp(log(qn) + pgamma(x, th[1] + N, rate = th[2] + E, log.p = TRUE)) +
                     exp(log(1 - qn) + pgamma(x, th[3] + N, rate = th[4] + E, log.p = TRUE))
  lo <- rep(0, length(N)); hi <- rep(1e5, length(N))
  for (i in seq_len(100L)) { m <- 0.5 * (lo + hi); up <- cdf(m) > pct / 100
                             hi[up] <- m[up]; lo[!up] <- m[!up] }
  0.5 * (lo + hi)
}

SETS <- c(Unrestricted = "full", Monopharmacy = "mono")
pri  <- lapply(SETS, function(tag) rd(paste0("ebgm_prior_", tag, ".csv")))
if (any(vapply(pri, is.null, logical(1)))) {
  rec("MGPS hyperparameters published", NA,
      paste("no", paste0("ebgm_prior_", SETS[vapply(pri, is.null, logical(1))], ".csv",
                         collapse = " or "), "in this folder"))
} else {
  PARS <- c("alpha1", "beta1", "alpha2", "beta2", "P")
  th   <- lapply(pri, function(p) setNames(as.numeric(p$used), p$parameter)[PARS])
  rec("MGPS hyperparameters are five finite numbers per reference set",
      all(vapply(th, function(x) length(x) == 5L && all(is.finite(x)) && all(x > 0), logical(1))),
      paste(vapply(names(th), function(k)
        paste0(k, ": ", paste(sprintf("%.5f", th[[k]]), collapse = "/")), character(1)),
        collapse = " | "))
  # the two sets are different tables, so one prior fitted to both would be wrong
  rec("each reference set has its own prior", !isTRUE(all.equal(th[[1]], th[[2]])),
      "unrestricted and monopharmacy hyperparameters differ")
  # openEBGM rewrites a prior it considers extreme before scoring with it; the
  # run records what it fitted and what it used, and those should be the same
  mv <- rbindlist(lapply(names(pri), function(k)
    pri[[k]][abs(fitted - used) > 1e-12][, set := k]), fill = TRUE)
  rec("the prior used is the prior fitted", nrow(mv) == 0L,
      sprintf("%d parameters were moved before scoring", nrow(mv)), mv)

  g <- da[!is.na(ebgm) & A >= 3 & analysis %chin% names(th)]
  if (!nrow(g)) rec("MGPS follows from a, b, c, d and the prior", NA, "no estimated rows") else {
    g[, E := (A + B) / N * (A + C)]
    g[, c("ebgm_x", "ebgm05_x") := {
      t5 <- th[[.BY$analysis]]
      q  <- ebgm_qn(t5, A, E)
      .(ebgm_point(t5, A, E, q), ebgm_quant(t5, A, E, q, 5))
    }, by = analysis]
    for (p in list(c("ebgm", "ebgm_x"), c("ebgm05", "ebgm05_x"))) {
      bad <- g[!agree(get(p[1]), get(p[2]))]
      rec(paste(p[1], "follows from a, b, c, d and the prior"), nrow(bad) == 0L,
          sprintf("%d of %d rows differ", nrow(bad), nrow(g)),
          bad[, c(key, "a", "ebgm", "ebgm05", p[2]), with = FALSE])
    }
    # EBGM05 is the 5th percentile, not the 2.5th that faers::phv_ebgm returned
    half <- g[, .(analysis, A, E)][, {
      t5 <- th[[.BY$analysis]]; q <- ebgm_qn(t5, A, E)
      .(q025 = ebgm_quant(t5, A, E, q, 2.5))
    }, by = analysis]$q025
    rec("EBGM05 is the 5th percentile and not the 2.5th",
        !all(agree(g$ebgm05, half)),
        sprintf("%d of %d rows would also match the 2.5th", sum(agree(g$ebgm05, half)), nrow(g)))
  }
}

# ── ledger ──────────────────────────────────────────────────────────────────
cat("\n== headline numbers ==\n")
led <- da[drug_group == "SSRI" & drug == "SSRI" &
          outcome %chin% c("Homicide", "Completed Suicide", hsn),
          .(analysis, outcome, a, ror, ror_ci_low, ror_ci_high, signal)]
print(unique(led, by = c("analysis", "outcome"))[order(analysis, outcome)])
cat("\nreference set and group sizes\n")
print(unique(da[, .(reference = fmt(max(N, na.rm = TRUE)), group_size = fmt(max(A + B, na.rm = TRUE))),
                by = .(analysis, drug_group)])[order(analysis, drug_group)])

# ── report ──────────────────────────────────────────────────────────────────
cat("\n")
for (i in seq_len(nrow(logdt)))
  cat(sprintf("  %-5s %-55s %s\n", logdt$status[i], substr(logdt$check[i], 1, 55), logdt$detail[i]))

np <- sum(logdt$status == "PASS"); nf <- sum(logdt$status == "FAIL"); ns <- sum(logdt$status == "SKIP")
cat(sprintf("\n%d passed, %d failed, %d skipped\n", np, nf, ns))

if (length(fails)) {
  cat("\n== failing rows ==\n")
  for (k in names(fails)) { cat("\n", k, "\n", sep = ""); print(fails[[k]]) }
}

fwrite(logdt, file.path(RES, "reconciliation_report.csv"))
cat("\nwrote ", file.path(RES, "reconciliation_report.csv"), "\n", sep = "")
# Exit non-zero so a failure stops run_all.R and fails a CI job. Not when the
# session is interactive: sourcing this in RStudio would close the IDE just as
# a check failed, taking the report above with it.
if (nf > 0L && !interactive()) quit(status = 1L)
