# =============================================================================
# tools/check_drug_counts.R - how many distinct drugs a report names.
#
#   Rscript tools/check_drug_counts.R
#
# The monopharmacy reference set is every report naming exactly one drug, so
# the count this makes is what decides which reports are in it, and a report
# miscounted here moves between the two analyses without anything else noticing.
# It is also easy to get wrong in both directions, so the cases below are kept
# as a record of what the rule is meant to do.
#
# The block under test is read out of R/groups.R rather than copied here, so
# this cannot drift away from the code it checks.
#
# Two FAERS fields carry the drug: `drugname`, written by the reporter, and
# `prod_ai`, the active ingredient. The dictionary holds the study drugs under
# their generic and brand names, but not every brand in FAERS, so a report may
# name one drug by a brand the dictionary does not hold and again by its
# ingredient. Resolving the name through prod_ai is what stops that counting
# twice - and what the cases below keep honest.
# =============================================================================
suppressPackageStartupMessages(library(data.table))
source("config.R"); source("R/definitions.R"); source("R/groups.R")

known <- toupper(names(c(drug_terms("study"), drug_terms("comparator"))))

src   <- readLines("R/groups.R")
# The whole counting body, from the drop of non-drug entries to the end of the
# prod_ai block. Anything left outside this range is untested: the cases below
# would then pass or fail on code that is not the code that runs, which has
# already happened twice.
start <- grep('^  drug <- drop_non_drugs\\(drug\\)$', src)
ifln  <- grep('^  if \\("prod_ai" %in% dcols\\) \\{', src)
stopifnot(length(start) == 1L, length(ifln) == 1L, ifln > start)
depth <- 0L; end <- NA_integer_
for (i in seq(ifln, length(src))) {
  depth <- depth + lengths(regmatches(src[i], gregexpr("\\{", src[i])))
  depth <- depth - lengths(regmatches(src[i], gregexpr("\\}", src[i])))
  if (depth == 0L) { end <- i; break }
}
block <- parse(text = paste(src[start:end], collapse = "\n"))

count <- function(drugname, prod_ai) {
  e <- new.env(parent = globalenv())
  e$dcols <- c("primaryid", "drugname", "prod_ai")
  e$drug  <- data.table(primaryid = 1L, drugname = drugname, prod_ai = prod_ai)
  eval(block, envir = e)
  if (any(c("aikey", "rawname", "nm", "dict") %in% names(e$drug)))
    stop("the block left a working column behind")
  uniqueN(e$drug$dkey)
}

SYM <- "FLUOXETINE HYDROCHLORIDE\\OLANZAPINE"   # how FAERS writes a combination

cases <- list(
  # ---- one drug written more than once -------------------------------------
  list("same unlisted brand twice, prod_ai on one line only",
       c("BRANDX", "BRANDX"), c("SERTRALINE", NA), 1L),
  list("same, three lines, prod_ai on two of them",
       c("BRANDX", "BRANDX", "BRANDX"), c("SERTRALINE", NA, "SERTRALINE"), 1L),
  list("same, the empty line is \"\" rather than NA",
       c("BRANDX", "BRANDX"), c("SERTRALINE", ""), 1L),
  list("brand once and its ingredient named separately",
       c("BRANDX", "SERTRALINE"), c("SERTRALINE", NA), 1L),

  # ---- names that share a canonical key but are different drugs ------------
  # Two SSRIs whose names share a trailing bracket. Keying on that bracket
  # would make them one drug; they must stay two.
  list("two SSRIs whose names both key as TABLETS",
       c("SERTRALINE (TABLETS)", "FLUOXETINE (TABLETS)"), c("SERTRALINE", "FLUOXETINE"), 2L),
  list("the same, with prod_ai missing on one line",
       c("SERTRALINE (TABLETS)", "FLUOXETINE (TABLETS)"), c("SERTRALINE", NA), 2L),

  # ---- a clarifier in brackets is not an ingredient -----------------------
  # Older reports put the same clarifier on every line. When a name was keyed
  # by its last bracket, every drug on such a report merged into one.
  list("two non-study drugs sharing (UNKNOWN)",
       c("ALPRAZOLAM (UNKNOWN)", "ZOLPIDEM (UNKNOWN)"), c(NA, NA), 2L),
  list("two non-study drugs sharing (NO PREF. NAME)",
       c("ACETAMINOPHEN/OXYCODONE (NO PREF. NAME)", "COCAINE (NO PREF. NAME)"), c(NA, NA), 2L),
  list("a drug, and the same drug with a clarifier",
       c("SERTRALINE", "SERTRALINE (UNKNOWN)"), c(NA, NA), 1L),

  # ---- placeholders in the drug field --------------------------------------
  # "no co-medication was given" is not a drug; "a drug was given but not
  # named" is one, and keeps the report out of monopharmacy (Methods).
  list("a drug and a no-co-medication placeholder",
       c("FLUOXETINE", "NO CONCURRENT MEDICATION"), c(NA, NA), 1L),
  list("a drug and an unidentified co-medication",
       c("FLUOXETINE", "UNKNOWN DRUG"), c(NA, NA), 2L),

  # ---- two drugs that share a clarifier in brackets ------------------------
  # Real FAERS shapes: older reports put the same clarifier on every line, and
  # a name keyed by its last bracket collapsed them all into one drug.
  list("two drugs both ending (NGX)",
       c("SERTRALINE (NGX)", "METHYLENE BLUE (NGX)"), c(NA, NA), 2L),
  list("two drugs both ending (NO PREF. NAME)",
       c("CITALOPRAM (NO PREF. NAME)", "COCAINE (NO PREF. NAME)"), c(NA, NA), 2L),
  list("two SSRIs sharing a clarifier, no prod_ai on either",
       c("SERTRALINE (TABLETS)", "FLUOXETINE (TABLETS)"), c(NA, NA), 2L),

  # ---- combination products -------------------------------------------------
  # FAERS marks a combination in prod_ai with a backslash. Such a prod_ai is
  # never used to resolve a name, so it cannot merge two drugs into one.
  list("combination in prod_ai, plus one ingredient separately",
       c("SYMBYAX", "FLUOXETINE"), c(SYM, NA), 2L),
  # A name is read as the dictionary reads it for group membership, with no
  # parsing of separators: a combination written in the name counts as the
  # study drug it contains. This is the accepted limit of the rule.
  list("combination in the name, plus one ingredient separately",
       c("OLANZAPINE AND FLUOXETINE", "FLUOXETINE"), c(NA, NA), 1L),
  list("combination in a parenthetical, plus one ingredient separately",
       c("SYMBYAX (OLANZAPINE, FLUOXETINE)", "FLUOXETINE"), c(NA, NA), 1L),

  # ---- one drug written with a strength, form, schedule or maker -----------
  list("brand, and the same brand written with a strength",
       c("PAXIL", "PAXIL, 20 MG"), c(NA, NA), 1L),
  list("brand, and the generic written with a dose form",
       c("ZOLOFT", "SERTRALINE HYDROCHLORIDE, TABLETS"), c(NA, NA), 1L),
  list("generic, then a strength with a multiplier",
       c("ZOLOFT", "ZOLOFT, 25MG X 2"), c(NA, NA), 1L),
  list("generic, then a strength and a FAERS frequency",
       c("CELEXA", "CITALOPRAM, 20 MG, 1 IN 1 D"), c(NA, NA), 1L),
  list("generic, then a manufacturer",
       c("PAXIL", "PAROXETINE, MYLAN"), c(NA, NA), 1L),
  list("foreign brand with foreign form text",
       c("PAXIL", "DEROXAT 20 MG, COMPRIME PELLICULE SECABLE"), c(NA, NA), 1L),
  list("brand and generic written with a slash",
       c("SERTRALINE", "ZOLOFT/SERTRALINE"), c(NA, NA), 1L),
  list("generic, a slash and a manufacturer",
       c("FLUOXETINE", "FLUOXETINE HCL/SANDOZ"), c(NA, NA), 1L),
  list("generic, then a strength written as a ratio",
       c("ZOLOFT", "SERTRALINE 50 MG/ML ORAL CONCENTRATE"), c(NA, NA), 1L),
  list("brand, then the manufacturer",
       c("CHANTIX", "CHANTIX, PFIZER"), c(NA, NA), 1L),

  # ---- reports that really do name more than one drug ----------------------
  list("one brand and one unrelated drug",
       c("BRANDX", "IBUPROFEN"), c("SERTRALINE", NA), 2L),
  list("two different SSRIs under plain names",
       c("ZOLOFT", "PROZAC"), c(NA, NA), 2L),
  list("one name recorded against two different ingredients",
       c("BRANDX", "BRANDX"), c("SERTRALINE", "FLUOXETINE"), 2L),
  list("two unrelated names, no prod_ai at all",
       c("BRANDX", "BRANDY"), c(NA, NA), 2L),
  list("all six SSRIs",
       c("FLUOXETINE", "PAROXETINE", "SERTRALINE",
         "CITALOPRAM", "ESCITALOPRAM", "FLUVOXAMINE"), rep(NA_character_, 6L), 6L),

  # ---- single lines ---------------------------------------------------------
  list("one SSRI", "SERTRALINE", NA_character_, 1L),
  list("one drug the dictionary does not hold", "ASPIRIN", NA_character_, 1L))

cat(sprintf("%-62s %6s %6s\n", "report names", "drugs", "counted"))
failed <- 0L
for (cs in cases) {
  got <- count(cs[[2]], cs[[3]]); want <- cs[[4]]
  if (got != want) failed <- failed + 1L
  cat(sprintf("%-62s %6d %6d  %s\n", cs[[1]], want, got,
              if (got == want) "" else "WRONG"))
}
cat(sprintf("\n%d of %d correct\n", length(cases) - failed, length(cases)))
if (failed) quit(status = 1L)
