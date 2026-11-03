# ==============================================================================
# tools/new_cases_for_review.R - the reports a re-run added to the case
# listings, and the ones it dropped.
#
# The manual removal list was built by reading the SSRI case listings, so a
# report entering a listing it was never in has not been through that review.
# Standardizing reaction terms to their preferred term (R/definitions.R) brings such
# reports in: MURDER counts towards Homicide, SELF MUTILATION towards
# Intentional self-injury. This writes those arrivals out in full, with the
# fields duplicate review needs, so the check covers the arrivals rather than
# every listing again.
#
# Each arrival is matched against every other report in the listing it joins:
#   same_caseid       another report carries the same CASEID, which is how a
#                     FAERS follow-up duplicate usually shows
#   match_to          the report it most resembles, and
#   match_score       the share of fields they agree on, counting only fields
#                     recorded on both - so two reports are not called alike
#                     because the same boxes are empty on each
#   match_on          which fields those were
#
# The second kind is the one that matters here. A case resubmitted under a new
# CASEID shares no identifier with its original, so only the patient and the
# event give it away.
#
# Reads  <new>/cases_*.csv and <new>/mono_cases_*.csv, the same files under
#        <reference>/, and data/manual_dedup_blacklist.csv
# Writes <new>/new_cases_for_review.csv
#
# Usage: Rscript tools/new_cases_for_review.R [new_dir] [reference_dir]
#        defaults to results/ against expected_2004_2025/
# ==============================================================================
suppressPackageStartupMessages({ library(data.table) })

.args <- commandArgs(trailingOnly = TRUE)
NEW <- if (length(.args) >= 1L) .args[1] else Sys.getenv("FAERS_RESULTS_DIR", unset = "results")
REF <- if (length(.args) >= 2L) .args[2] else "expected_2004_2025"
for (d in c(NEW, REF)) if (!dir.exists(d)) stop(d, " not found", call. = FALSE)

PATTERN <- "^(mono_)?cases_.*\\.csv$"
ID      <- function(x) as.character(x$primaryid)
rd      <- function(dir, f) {
  p <- file.path(dir, f)
  if (!file.exists(p)) return(NULL)
  fread(p, showProgress = FALSE, colClasses = list(character = "primaryid"))
}

files <- union(list.files(NEW, PATTERN), list.files(REF, PATTERN))
# the with-indications listings repeat their parent file's reports
files <- sort(files[!grepl("_with_indications\\.csv$", files)])
if (!length(files)) stop("no case listings in either folder", call. = FALSE)

# Fields that can actually tell two reports apart. Sex, country and the SSRI
# take too few values to be evidence of anything - a male in the United States
# on paroxetine describes thousands of reports - so they rule a pair out rather
# than in: a pair is dropped if any of them is recorded on both and disagrees,
# but agreeing on them earns no score.
DISCRIM <- c("age_in_years", "wt", "event_dt", "start_dt")
CONSIST <- c("sex", "reporter_country", "ssri_drug")
PROFILE <- c(DISCRIM, CONSIST)
MIN_SHARED <- 2L   # agreeing discriminating fields before a pair is worth a look

# For each arrival, the report in the same listing it most resembles. Only
# fields recorded on both sides are compared, so a pair of sparse reports is
# not called a match on the strength of what neither of them records.
best_match <- function(add, all) {
  blank <- function(x) is.na(x) | trimws(as.character(x)) %chin% c("", "NA")
  dsc <- intersect(DISCRIM, names(all)); cns <- intersect(CONSIST, names(all))
  gv  <- function(d, k) as.character(d[[k]])
  A <- lapply(c(dsc, cns), function(k) gv(add, k)); names(A) <- c(dsc, cns)
  B <- lapply(c(dsc, cns), function(k) gv(all, k)); names(B) <- c(dsc, cns)
  ia <- as.character(add$primaryid); ib <- as.character(all$primaryid)
  none <- data.table(match_to = NA_character_, match_score = NA_real_,
                     shared = NA_integer_, match_on = NA_character_)
  rbindlist(lapply(seq_len(nrow(add)), function(i) {
    cmp <- function(k) !blank(A[[k]][i]) & !blank(B[[k]])
    agr <- function(k) cmp(k) & A[[k]][i] == B[[k]]
    # a pair is out if anything low-cardinality is recorded on both and differs
    keep <- Reduce(`&`, lapply(cns, function(k) !cmp(k) | agr(k)), rep(TRUE, length(ib)))
    keep <- keep & ib != ia[i]
    ok   <- Reduce(`+`, lapply(dsc, function(k) as.integer(cmp(k))))
    hit  <- Reduce(`+`, lapply(dsc, function(k) as.integer(agr(k))))
    elig <- keep & hit >= MIN_SHARED
    if (!any(elig)) return(none)
    sc <- ifelse(elig, hit + hit / pmax(ok, 1L), -1)   # agreements first, then purity
    j  <- which.max(sc)
    data.table(match_to = ib[j], match_score = round(hit[j] / ok[j], 3), shared = ok[j],
               match_on = paste(c(dsc[vapply(dsc, function(k) agr(k)[j], logical(1))],
                                  cns[vapply(cns, function(k) agr(k)[j], logical(1))]),
                                collapse = ","))
  }))
}

out <- rbindlist(lapply(files, function(f) {
  a <- rd(NEW, f); b <- rd(REF, f)
  if (is.null(a)) {
    cat(sprintf("  %-44s listing is missing from %s\n", f, NEW)); return(NULL)
  }
  new_ids <- if (is.null(b)) ID(a) else setdiff(ID(a), ID(b))
  gone    <- if (is.null(b)) character() else setdiff(ID(b), ID(a))
  cat(sprintf("  %-44s %5d rows  %+5d added  %+5d dropped\n",
              f, nrow(a), length(new_ids), -length(gone)))
  if (!length(new_ids) && !length(gone)) return(NULL)

  add <- a[ID(a) %chin% new_ids]
  if (nrow(add)) {
    add[, same_caseid := if ("caseid" %chin% names(a))
          vapply(as.character(caseid), function(k) sum(as.character(a$caseid) == k) - 1L,
                 integer(1)) else NA_integer_]
    add <- cbind(add, best_match(add, a))
  }
  drop <- if (!length(gone)) NULL else b[ID(b) %chin% gone]
  rbindlist(list(
    if (nrow(add))                    cbind(data.table(listing = f, change = "added"),   add),
    if (!is.null(drop) && nrow(drop)) cbind(data.table(listing = f, change = "dropped"), drop)),
    use.names = TRUE, fill = TRUE)
}), use.names = TRUE, fill = TRUE)

if (is.null(out) || !nrow(out)) {
  cat("\nno listing changed: nothing to review\n"); quit(save = "no")
}

# A report already on the removal list cannot reach a listing, because the list
# is applied to the reference set before any signal is computed. If one does,
# the listings and the list disagree and the run should not be trusted.
BL <- file.path("data", "manual_dedup_blacklist.csv")
if (file.exists(BL)) {
  blids <- as.character(fread(BL, showProgress = FALSE)$primaryid)
  bad   <- out[change == "added" & as.character(primaryid) %chin% blids]
  if (nrow(bad))
    stop(nrow(bad), " arriving reports are on the removal list - the listings and ",
         BL, " disagree:\n  ", paste(head(unique(bad$primaryid), 10), collapse = ", "),
         call. = FALSE)
}

lead <- intersect(c("listing", "change", "match_score", "shared", "match_on", "match_to",
                    "same_caseid", "primaryid", "caseid", "year", "quarter", "event_dt",
                    "start_dt", "fda_dt", "sex", "age_in_years", "wt", "reporter_country",
                    "occp_cod", "i_f_code", "ssri_drug", "all_pts"), names(out))
setcolorder(out, lead)
setorderv(out, c("change", "match_score", "same_caseid", "listing"), c(1L, -1L, -1L, 1L))
fwrite(out, file.path(NEW, "new_cases_for_review.csv"))

cat("\n", nrow(out[change == "added"]), " added, ", nrow(out[change == "dropped"]),
    " dropped, across ", uniqueN(out$listing), " listings\n", sep = "")
cat("\nadded, by listing\n")
print(out[change == "added",
          .(reports = .N, same_caseid = sum(same_caseid > 0L, na.rm = TRUE),
            every_field_agrees = sum(match_score == 1, na.rm = TRUE),
            most_agree = sum(match_score >= 0.75 & match_score < 1, na.rm = TRUE),
            nothing_to_compare = sum(is.na(match_score))),
          by = listing][order(-reports)])

# the pairs worth a person's time, strongest first
top <- out[change == "added" & !is.na(match_score) & match_score >= 0.75]
if (nrow(top)) {
  cat("\n", nrow(top), " arrival(s) agree with an existing report on at least three",
      " quarters of the fields both record:\n", sep = "")
  show <- intersect(c("listing", "match_score", "shared", "match_on", "primaryid", "match_to",
                      "caseid", "year", "sex", "age_in_years", "event_dt", "ssri_drug"),
                    names(top))
  print(head(top[, ..show], 40L))
} else cat("\nno arrival closely resembles an existing report\n")
cat("\nwrote ", file.path(NEW, "new_cases_for_review.csv"), "\n", sep = "")
cat("Review the added reports and add any duplicate to ", BL, ".\n", sep = "")
