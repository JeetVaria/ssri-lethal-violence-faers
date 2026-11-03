# ==============================================================================
# tools/compare_to_reference.R - check a fresh run reproduces earlier results.
#
# Compares every CSV that exists in both folders: same rows (after sorting on
# all shared columns), same values; numbers compared with a small tolerance.
#
#   Rscript tools/compare_to_reference.R <reference_folder> [results_folder] [--ignore-bayes]
#
# --ignore-bayes: skip BCPNN/EBGM columns and the signal flags that depend on
# them (for runs made with FAERS_QUICK_TEST=TRUE).
# ==============================================================================
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
ignore_bayes <- "--ignore-bayes" %in% args
args <- args[args != "--ignore-bayes"]   # not setdiff(): it would also drop a repeated folder
BAYES_COLS <- "^(bcpnn|ebgm|ic_signal|ebgm_signal|full_consensus|MGPS|BCPNN|CONSENSUS)"
if (!length(args)) stop("usage: Rscript tools/compare_to_reference.R <reference_folder> [results_folder]")
ref_dir <- args[1]
new_dir <- if (length(args) > 1) args[2] else "results"
tol <- 1e-6

ref_files <- setdiff(list.files(ref_dir, pattern = "\\.csv$"), "pipeline_log.txt")
new_files <- setdiff(list.files(new_dir, pattern = "\\.csv$"), "pipeline_log.txt")
files <- intersect(ref_files, new_files)

# Comparing only the shared files would pass a run that produced a different
# set of them, so report the difference before anything else.
# Written by a run but deliberately absent from the reference folders.
# removed_by_dedup.csv is far too large for the repository (hundreds of MB over
# the full period). llt_mappings_applied.csv holds LLT-to-preferred-term pairs,
# which are MedDRA content and cannot be redistributed, so it is never shipped.
# removed_by_dedup_summary.csv and qa_table_counts.csv describe the parse and
# deduplication, so step 01 writes them only when it rebuilds: a run that reuses
# master_dedup.rds does not, and route 3 skips step 01 altogether. Their absence
# says nothing about whether the published numbers reproduce.
NOT_SHIPPED <- c("removed_by_dedup.csv", "llt_mappings_applied.csv",
                 "removed_by_dedup_summary.csv", "qa_table_counts.csv")
# Produced by tools/refit_ebgm_prior.R, which run_all.R does not call, so the
# reference carries them and an ordinary run does not.
TOOL_ONLY <- c("ebgm_prior_restarts_full.csv", "ebgm_prior_restarts_mono.csv")
# NOT_SHIPPED is subtracted from both sides: the two-year reference folder
# carries removed_by_dedup_summary.csv and qa_table_counts.csv, and a run that reuses
# master_dedup.rds does not write them, so on the reference side their absence
# from the run says nothing either.
only_ref <- setdiff(setdiff(ref_files, new_files), c(TOOL_ONLY, NOT_SHIPPED))
only_new <- setdiff(setdiff(new_files, ref_files), NOT_SHIPPED)
if (length(only_ref))
  cat(sprintf("\n%d file(s) in %s but NOT produced by this run:\n  %s\n",
              length(only_ref), ref_dir, paste(only_ref, collapse = "\n  ")))
if (length(only_new))
  cat(sprintf("\n%d file(s) produced by this run but NOT in %s:\n  %s\n",
              length(only_new), ref_dir, paste(only_new, collapse = "\n  ")))
report <- rbindlist(lapply(files, function(f) {
  a <- fread(file.path(ref_dir, f), colClasses = "character")
  b <- fread(file.path(new_dir, f), colClasses = "character")
  a_names <- names(a); b_names <- names(b)
  cols <- intersect(names(a), names(b))
  if (ignore_bayes) cols <- cols[!grepl(BAYES_COLS, cols)]
  # A column the run gained or lost is a failure to reproduce just as much as a
  # differing value is. Comparing only the shared columns reported "identical"
  # for a listing that had gained a column or had one renamed, which is how two
  # reference files drifted a version behind the code without this noticing.
  gained <- setdiff(b_names, a_names); lost <- setdiff(a_names, b_names)
  cols_msg <- c(if (length(lost))   paste("COLUMNS MISSING FROM THE RUN:", paste(lost, collapse = ", ")),
                if (length(gained)) paste("COLUMNS NEW IN THE RUN:", paste(gained, collapse = ", ")))
  # The shared columns are compared whether or not the column set differs, so one
  # line names every way the file differs. Reporting the column difference alone
  # hid a changed value or a changed row count in the same file.
  rest_msg <- if (nrow(a) != nrow(b)) {
    sprintf("ROWS DIFFER (%d vs %d)", nrow(a), nrow(b))
  } else if (!length(cols)) NULL else {
    a <- a[, ..cols]; b <- b[, ..cols]
    setorderv(a, cols); setorderv(b, cols)
    bad <- vapply(cols, function(cn) {
      x <- a[[cn]]; y <- b[[cn]]
      nx <- suppressWarnings(as.numeric(x)); ny <- suppressWarnings(as.numeric(y))
      num <- !is.na(nx) & !is.na(ny)
      any(abs(nx[num] - ny[num]) > tol * pmax(1, abs(nx[num]))) || any(x[!num] != y[!num], na.rm = TRUE) ||
        any(is.na(x[!num]) != is.na(y[!num]))
    }, logical(1))
    if (any(bad)) paste("VALUES DIFFER in:", paste(cols[bad], collapse = ", ")) else NULL
  }
  status <- if (length(cols_msg) || length(rest_msg)) paste(c(cols_msg, rest_msg), collapse = "; ") else "identical"
  data.table(file = f, status = status,
             columns_only_in_reference = paste(setdiff(a_names, b_names), collapse = ", "))
}))
print(report, nrows = Inf)
# written beside the repository, not into the results folder, so a later run
# does not compare the report to itself
fwrite(report, "comparison_report.csv")
cat(sprintf("\n%d of %d shared files identical\n", report[status == "identical", .N], nrow(report)))
# Figures are rendered from the tables above and their bytes depend on the
# graphics device, so they are not compared here; check them by eye.
figs <- length(list.files(file.path(new_dir, "figures")))
if (figs) cat(sprintf("%d figure(s) written but not compared - they render from the tables above\n", figs))
differed <- report[status != "identical", .N]
if (length(only_ref) || length(only_new))
  cat(sprintf("%d file(s) only in %s, %d only in %s - the run did not reproduce the file set\n",
              length(only_ref), ref_dir, length(only_new), new_dir))
# A differing value is a failure to reproduce just as much as a missing file
# is, so both set a non-zero exit status rather than only printing. This is the
# script the README names as the check for routes 1 to 3, so anything chaining
# it with && or running it in CI reads that status rather than the output.
if (differed)
  cat(sprintf("%d file(s) differ in their values - see comparison_report.csv\n", differed))
if ((length(only_ref) || length(only_new) || differed) && !interactive())
  quit(status = 1L)
