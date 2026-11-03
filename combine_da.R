# ==============================================================================
# combine_da.R - every disproportionality estimate in one table.
#
# The pipeline writes the estimates across a dozen files, split by drug group,
# by analysis, and by whether the row is pooled, per drug, a PT-group composite
# or the homicide-suicide composite. This gathers all of them into one tidy CSV with
# four columns identifying the row and the contingency table alongside every
# estimate, so the whole analysis can be read, sorted and checked in one place.
#
# Reads  results/*.csv written by steps 02 and 03
# Writes results/disproportionality_all.csv and results/SI_DATA_DA.csv
#        (or FAERS_RESULTS_DIR, where that is set)
#
# Usage: Rscript combine_da.R
# ==============================================================================
suppressPackageStartupMessages({ library(data.table) })

# Folder to read and write: the first argument if it names one, otherwise
# FAERS_RESULTS_DIR, otherwise results/. Pass "expected_2004_2025" to work on the files
# behind the article:  Rscript combine_da.R expected_2004_2025
.args <- commandArgs(trailingOnly = TRUE)
RES <- Sys.getenv("FAERS_RESULTS_DIR", unset = "results")
if (length(.args) && dir.exists(.args[1])) RES <- .args[1]
rd  <- function(f) {
  p <- file.path(RES, f)
  if (!file.exists(p)) { message("  (missing, skipped) ", f); return(NULL) }
  x <- fread(p)
  if ("ids" %in% names(x)) x[, ids := NULL]
  x
}

# label each source: analysis, drug group, and how the drug column is derived
sources <- list(
  list(f = "ssri_targeted_signals.csv",            an = "Unrestricted",  co = "SSRI",        dr = "All SSRIs"),
  list(f = "ssri_perdrug_pt_signals.csv",          an = "Unrestricted",  co = "SSRI",        dr = NA),
  list(f = "varenicline_targeted_signals.csv",     an = "Unrestricted",  co = "Varenicline", dr = "Varenicline"),
  list(f = "group_pooled_full.csv",                an = "Unrestricted",  co = NA,            dr = NA),
  list(f = "hs_composite_all_pooled.csv",          an = "Unrestricted",  co = NA,            dr = NA),
  list(f = "ssri_perdrug_hs_composite.csv",        an = "Unrestricted",  co = "SSRI",        dr = NA),
  list(f = "ssri_mono_targeted_signals.csv",       an = "Monopharmacy",  co = "SSRI",        dr = "All SSRIs"),
  list(f = "ssri_mono_perdrug_pt_signals.csv",     an = "Monopharmacy",  co = "SSRI",        dr = NA),
  list(f = "varenicline_mono_targeted_signals.csv",an = "Monopharmacy",  co = "Varenicline", dr = "Varenicline"),
  list(f = "group_pooled_mono.csv",                an = "Monopharmacy",  co = NA,            dr = NA),
  list(f = "hs_composite_all_mono.csv",            an = "Monopharmacy",  co = NA,            dr = NA),
  list(f = "ssri_mono_perdrug_hs_composite.csv",   an = "Monopharmacy",  co = "SSRI",        dr = NA)
)

# A published run writes all twelve. A missing one quietly shortens the combined
# table and empties a whole block of the Supplementary sheet, where it reads
# exactly like the n < 3 suppression - so refuse rather than skip. Set
# FAERS_ALLOW_PARTIAL=TRUE to combine a deliberately partial run.
gone <- Filter(function(f) !file.exists(file.path(RES, f)),
               vapply(sources, `[[`, character(1), "f"))
if (length(gone) && !identical(Sys.getenv("FAERS_ALLOW_PARTIAL"), "TRUE"))
  stop("missing ", length(gone), " of ", length(sources), " source files in ", RES, ":\n  ",
       paste(gone, collapse = "\n  "),
       "\nrun the pipeline first, or set FAERS_ALLOW_PARTIAL=TRUE to combine what is there")

parts <- lapply(sources, function(s) {
  x <- rd(s$f); if (is.null(x) || !nrow(x)) return(NULL)
  x[, analysis := s$an]
  # drug group: from the file's own drug_group column where it has one (the H-S
  # composites do); the group column in group_pooled_*.csv is the PT grouping
  if (!is.na(s$co)) x[, group_label := s$co]
  else if ("drug_group" %in% names(x)) x[, group_label := fifelse(grepl("varenicline", drug_group), "Varenicline", "SSRI")]
  else x[, group_label := NA_character_]
  # drug: named in the file, or fixed by the source. The H-S composite files
  # name no drug and hold both groups, so the comparator row must take its name
  # from drug_group or it collides with the SSRI row and is lost.
  if (!is.na(s$dr)) x[, drug_label := s$dr]
  else if ("drug" %in% names(x)) x[, drug_label := as.character(drug)]
  else if ("drug_group" %in% names(x))
    x[, drug_label := fifelse(grepl("varenicline", drug_group), "Varenicline", "All SSRIs")]
  else x[, drug_label := "All SSRIs"]
  x[, source_file := s$f]
  x
})
parts <- Filter(Negate(is.null), parts)
if (!length(parts)) stop("no result files found in ", RES, " - run the pipeline first")

all <- rbindlist(parts, use.names = TRUE, fill = TRUE)

# Ingrowing nail is the expected-null reaction for quality assurance, reported
# in the Supplementary "QA" sheet by step 04. It is not a study outcome, so it
# is kept out of the disproportionality table and the "DA" sheet.
if ("pt" %in% names(all)) all <- all[is.na(pt) | pt != "Ingrowing Nail"]

# the PT-group tables call the pooled SSRI drug "SSRI (POOLED)" and carry no
# drug_group column
all[, drug_label := fifelse(drug_label == "SSRI (POOLED)", "All SSRIs", drug_label)]
all[tolower(drug_label) == "varenicline", drug_label := "Varenicline"]
all[is.na(group_label),
    group_label := fifelse(drug_label == "Varenicline", "Varenicline", "SSRI")]

# the row is identified by its PT, or by its PT group where the table is a
# PT-group composite; the monopharmacy composite carries the same outcome under a
# different name, which would otherwise split the row in the wide sheet
oc  <- as.character(all$pt)
grp <- if ("group" %in% names(all)) as.character(all$group) else rep(NA_character_, nrow(all))
oc  <- fifelse(is.na(oc) | !nzchar(oc),
               fifelse(is.na(grp) | !nzchar(grp), NA_character_, paste0("(", grp, ")")), oc)
oc  <- sub(", mono\\)", ")", oc)
if (anyNA(oc)) stop(sum(is.na(oc)), " rows carry neither a PT nor a PT-group label")

ic   <- intersect(c("bcpnn_mcmc_ic", "bcpnn_norm_ic"), names(all))[1]
ic25 <- intersect(c("bcpnn_mcmc_ic_ci_low", "bcpnn_norm_ic_ci_low"), names(all))[1]

out <- data.table(
  analysis   = all$analysis,
  drug_group = all$group_label,
  drug       = all$drug_label,
  outcome    = oc,
  n_group    = if ("n_group" %in% names(all)) all$n_group else NA_integer_,
  a = all$a, b = all$b, c = all$c, d = all$d,
  ror = all$ror, ror_ci_low = all$ror_ci_low, ror_ci_high = all$ror_ci_high,
  prr = all$prr, prr_ci_low = all$prr_ci_low, prr_ci_high = all$prr_ci_high,
  chisq = all$chisq,
  ic    = if (!is.na(ic))   all[[ic]]   else NA_real_,
  ic025 = if (!is.na(ic25)) all[[ic25]] else NA_real_,
  ebgm  = all$ebgm, ebgm05 = all$ebgm_ci_low,
  ror_signal = all$ror_signal, prr_signal = all$prr_signal,
  ic_signal  = all$ic_signal,  ebgm_signal = all$ebgm_signal,
  signal     = all$full_consensus,
  source_file = all$source_file
)

num <- c("ror","ror_ci_low","ror_ci_high","prr","prr_ci_low","prr_ci_high","chisq","ic","ic025","ebgm","ebgm05")
# disproportionality_all.csv publishes these at 4 dp. The wide table below is
# printed at 2 dp and is formatted from the unrounded values held here: going
# through 4 dp first rounds twice, and a value just under a half crosses it -
# the homicide ROR 29.904986 becomes 29.905 and would print as 29.91, not 29.90.
for (cl in num) out[, (paste0(cl, ".raw")) := as.numeric(get(cl))]
for (cl in num) out[, (cl) := round(as.numeric(get(cl)), 4)]

# labels exactly as the Supplementary "DA" sheet prints them
out[, outcome := fifelse(grepl("^Homicide \\+ Completed [Ss]uicide", outcome), "[H-S]", outcome)]
out[, drug := fifelse(drug == "All SSRIs", "SSRI", tools::toTitleCase(drug))]

setorder(out, analysis, drug_group, outcome, -a, na.last = TRUE)
fwrite(out[, setdiff(names(out), paste0(num, ".raw")), with = FALSE],
       file.path(RES, "disproportionality_all.csv"))

cat("\nwrote ", file.path(RES, "disproportionality_all.csv"), "\n", sep = "")
cat("rows: ", nrow(out), "\n\n", sep = "")
print(out[, .N, by = .(analysis, drug_group)][order(analysis, drug_group)])
cat("\nrows meeting all four criteria: ", sum(out$signal %in% TRUE), "\n", sep = "")
miss <- c("a","b","c","d")[colSums(!is.na(out[, .(a,b,c,d)])) == 0L]
if (length(miss)) cat("WARNING: contingency cells missing for some rows: ", paste(miss, collapse = ", "), "\n")

# ── wide table matching the Supplementary "DA" sheet ─────────────────────────
# One row per outcome and drug, with the unrestricted and monopharmacy blocks
# side by side, formatted as the sheet presents them.

f2 <- function(x) ifelse(is.na(x), "", sprintf("%.2f", x))
f1 <- function(x) ifelse(is.na(x), "", sprintf("%.1f", x))
blk <- function(d) {
  # the n >= 3 rule the PT-level estimators and Fig. 3 already apply. The
  # composite rows carry it too (hs_composite() takes min_a = 3L), so this is a
  # guard rather than a correction: it keeps the rule visible in one place.
  shown <- !is.na(d$a) & d$a >= 3
  data.table(
  a = d$a, b = d$b, c = d$c, d_ = d$d,
  pt_total = ifelse(is.na(d$a) | is.na(d$c), NA_integer_, d$a + d$c),
  `ROR (95% CI)`  = ifelse(!shown | is.na(d$ror), "", paste0(f2(d$ror.raw), " (", f2(d$ror_ci_low.raw), " - ", f2(d$ror_ci_high.raw), ")")),
  `PRR (chi2)`      = ifelse(!shown | is.na(d$prr), "", paste0(f2(d$prr.raw), " (", f1(d$chisq.raw), ")")),
  `BCPNN (IC025)` = ifelse(!shown | is.na(d$ic),  "", paste0(f2(d$ic.raw), " (", f2(d$ic025.raw), ")")),
  `MGPS (EBGM05)` = ifelse(!shown | is.na(d$ebgm),"", paste0(f2(d$ebgm.raw), " (", f2(d$ebgm05.raw), ")")),
  CONSENSUS       = ifelse(is.na(d$signal), "", ifelse(d$signal, "Yes", "No")),
  `SUM (a + b)`   = ifelse(is.na(d$a) | is.na(d$b), NA_integer_, d$a + d$b),
  `SUM (N)`       = ifelse(is.na(d$a) | is.na(d$b) | is.na(d$c) | is.na(d$d),
                           NA_integer_, d$a + d$b + d$c + d$d))
}

key <- c("outcome", "drug")
u <- out[analysis == "Unrestricted"]
m <- out[analysis == "Monopharmacy"]
u <- u[!duplicated(u[, ..key])]
m <- m[!duplicated(m[, ..key])]

spine <- unique(rbind(u[, ..key], m[, ..key]))
drug_rank <- function(x) fifelse(x == "SSRI", 0L,
                        fifelse(x == "Varenicline", 9L, 5L))
spine <- spine[order(outcome, drug_rank(drug), drug)]

U <- merge(spine, u, by = key, all.x = TRUE, sort = FALSE)
M <- merge(spine, m, by = key, all.x = TRUE, sort = FALSE)

wide <- cbind(
  data.table(`PT, (PT group) or [Composite]` = spine$outcome, Drug = spine$drug, Unrestricted = ""),
  blk(U),
  data.table(Monopharmacy = ""),
  blk(M))
# the monopharmacy block repeats the same column names; make them unique for CSV
nm <- names(wide); n2 <- nm
dup <- duplicated(nm); n2[dup] <- paste0(nm[dup], "_mono")
n2[n2 == "d_"] <- "d"; n2[n2 == "d__mono"] <- "d_mono"
setnames(wide, n2)

fwrite(wide, file.path(RES, "SI_DATA_DA.csv"))
cat("wrote ", file.path(RES, "SI_DATA_DA.csv"), "  (", nrow(wide), " rows, ", ncol(wide), " columns)\n", sep = "")
