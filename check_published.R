# ==============================================================================
# check_published.R - check the numbers behind the article.
#
#   RStudio:   open this file and click Source.
#   Terminal:  Rscript check_published.R
#
# Needs the data.table and bit64 packages and nothing else: no FAERS download,
# no MedDRA license, about two minutes on an ordinary laptop.
#
# Rebuilds the estimate table from the per-step signal files, then runs the
# reconciliation checks over it. Reads and writes expected_2004_2025/ only, never results/, so
# your own runs and the article's output cannot be confused.
# ==============================================================================
if (!dir.exists("expected_2004_2025"))
  stop("expected_2004_2025/ is not here, so this is not the repository root.\n",
       "In RStudio: Session > Set Working Directory > To Source File Location,\n",
       "then click Source again. In a terminal: cd to the folder holding\n",
       "run_all.R and re-run.", call. = FALSE)

# Put back however this script ends, including on an error. Without that, an
# RStudio session that sourced it would keep writing every later run into
# expected_2004_2025/ and overwrite the published output that routes 1 and 3
# are compared against. on.exit() is no use here: at the top level of a sourced
# file it binds to the expression that registered it and fires straight away,
# undoing the setting before the body has run.
.prev_res_dir <- Sys.getenv("FAERS_RESULTS_DIR", unset = NA_character_)
.restore_res_dir <- function() {
  if (is.na(.prev_res_dir)) Sys.unsetenv("FAERS_RESULTS_DIR")
  else Sys.setenv(FAERS_RESULTS_DIR = .prev_res_dir)
}
Sys.setenv(FAERS_RESULTS_DIR = "expected_2004_2025")

tryCatch({
  source("combine_da.R")
  source("reconcile.R")
}, finally = .restore_res_dir())

# reconcile.R prints a PASS/FAIL line per check and writes
# expected_2004_2025/reconciliation_report.csv. No FAIL is the whole article
# reconciling against itself; any FAIL names the file and row to look at. The
# number of checks is not fixed - some run once per outcome file and per
# analysis - so read the FAIL count, not the total. A SKIP means an input that
# check needs is absent from the folder.

if (is.na(.prev_res_dir)) Sys.unsetenv("FAERS_RESULTS_DIR") else
  Sys.setenv(FAERS_RESULTS_DIR = .prev_res_dir)
rm(.prev_res_dir)
