# ==============================================================================
# R/definitions.R - study definitions, read from data/*.csv so they can be
# reviewed (and cited in the Supplementary Information) without reading code.
# ==============================================================================

.dict <- fread("data/drug_dictionary.csv")

# named list: drug -> character vector of generic/INN/salt/brand terms
drug_terms <- function(role_) {
  d <- .dict[role == role_]
  sapply(unique(d$drug), function(x) d[drug == x, term], simplify = FALSE)
}
# whole-word, case-insensitive regex matching any of the terms
terms_pattern <- function(terms) paste0("\\b(", paste(unique(unlist(terms, use.names = FALSE)), collapse = "|"), ")\\b")

ssri_drugs          <- drug_terms("study")          # 6 SSRIs (Supplementary "DRUG")
ssri_pattern        <- terms_pattern(ssri_drugs)
ssri_drug_pattern   <- lapply(ssri_drugs, terms_pattern)  # one pattern per SSRI
comparators         <- drug_terms("comparator")
varenicline_pattern <- terms_pattern(comparators$varenicline)
qa_drugs            <- drug_terms("qa")

# ── target preferred terms (normalized Title Case: see pt_std() below) ──
.pts     <- fread("data/pt_definitions.csv")
targets  <- .pts$pt
groups   <- setNames(.pts$group, .pts$pt)
pt_files <- setNames(as.list(.pts$case_file), .pts$pt)                 # unrestricted case lists
group_defs <- list(Suicidality  = names(groups[groups == "Suicidality"]),
                   Homicidality = names(groups[groups == "Homicidality"]))
HS_PT <- "Homicide + Completed suicide (co-reported)"

# ── standardizing reaction terms ─────────────────────────────────────────────
# FAERS stores the term as it was coded, and MedDRA is coded at Lowest Level
# Term (LLT): MURDER is an LLT of PT Homicide, SELF MUTILATION of Intentional
# self-injury, VIOLENT IDEATION of Violence-related symptom. Rolling an LLT up
# to its PT is the analyst's step, not the FDA's, so it happens here.
#
# The map is read from the MedDRA distribution the run is configured with, so
# every REACTION term in the database is standardized and not only the study
# outcomes. It is applied to the reaction table, which is what the counts, the
# case lists and the drug-event matrix are built from. Indications are not
# rolled up: Figure 2 and Table 1 report them as recorded, so an indication
# coded at Lowest Level Term appears under that name there.
#
# MedDRA is licensed and none of it is kept in this repository: the map is built
# at run time from the reader's own llt.asc and pt.asc and never written out
# whole. Anyone able to run this step already holds a license, because step 01
# cannot standardize the database without one. What each run does write is
# llt_mappings_applied.csv - the handful of terms that actually moved a count,
# with the reports behind each. That is the audit trail, but it pairs LLT text
# against preferred-term text, so it is MedDRA content and is not shipped in
# this repository either; running step 02 produces your own copy.
#
# Non-current LLTs are kept deliberately: a term retired from MedDRA is exactly
# the kind that sits in the older reports.

# The one way a reaction term is normalized, used on both sides of the lookup.
# FAERS carries leading, trailing and doubled spaces, which is why the faers
# package trims and collapses them too: a term normalized one way on one side
# of a join and another way on the other matches nothing, and the reports
# behind it disappear from the counts and the listings together.
# tools::toTitleCase() lowercases a standalone Roman numeral, so "bipolar i
# disorder" comes back as "Bipolar i Disorder". Restore those. Single letters
# are deliberately left alone: "a" is an article in several MedDRA terms
# ("...due to a general medical condition"), and uppercasing it would break
# more terms than it fixes.
# Standalone Roman numerals, uppercased. Both renderers use it - Title Case for
# the summary tables and sentence case for the case listings - so one MedDRA
# term is spelled the same way wherever it appears.
roman_up <- function(x) {
  gsub("(?i)(?<![[:alpha:]])(i|ii|iii|iv|v|vi|vii|viii|ix|x|xi|xii|xiii)(?![[:alpha:]])",
       "\\U\\1", x, perl = TRUE)
}

title_case <- function(x) roman_up(tools::toTitleCase(x))

pt_std <- function(x) title_case(tolower(trimws(gsub("[[:space:]]+", " ", x))))

.LLT_MAP <- NULL
llt_map <- function() {
  if (!is.null(.LLT_MAP)) return(.LLT_MAP)
  nz <- pt_std
  # distributions put the ASCII files in a MedAscii subfolder, loose in the
  # release folder, or a level deeper again, so the name is searched for
  asc <- function(nm) {
    p <- list.files(MEDDRA_DIR, pattern = paste0("^", nm, "$"), recursive = TRUE,
                    full.names = TRUE, ignore.case = TRUE)
    if (length(p)) p[1] else NA_character_
  }
  f <- c(asc("llt.asc"), asc("pt.asc"))
  if (anyNA(f))
    stop("reaction terms cannot be standardized: llt.asc and pt.asc were not found ",
         "anywhere under ", MEDDRA_DIR, ".\nSet MEDDRA_DIR in config.R to your MedDRA ",
         "distribution. MedDRA is licensed, so no copy of it is kept here.", call. = FALSE)
  log_msg("  MedDRA: reading ", dirname(f[1]))
  l <- fread(f[1], sep = "$", header = FALSE, quote = "", showProgress = FALSE,
             select = 1:3, col.names = c("llt_code", "llt_name", "pt_code"),
             colClasses = "character")
  p <- fread(f[2], sep = "$", header = FALSE, quote = "", showProgress = FALSE,
             select = 1:2, col.names = c("pt_code", "pt_name"), colClasses = "character")
  x <- merge(l, p, by = "pt_code")[, .(llt = nz(llt_name), pt = nz(pt_name))]
  m <- unique(x[llt != pt & nzchar(llt) & nzchar(pt)])       # identity rows do nothing
  # One LLT belongs to one PT. A repeat means the files were misread.
  if (anyDuplicated(m$llt))
    stop("the MedDRA term map has ", sum(duplicated(m$llt)),
         " repeated Lowest Level Terms - check the layout of ", f[1], call. = FALSE)
  .LLT_MAP <<- setNames(m$pt, m$llt)
  log_msg("  MedDRA: ", sfmt(length(.LLT_MAP)), " Lowest Level Terms map to a different ",
          "preferred term")
  .LLT_MAP
}

# A case listing's recorded-terms string ("A; B; C") rolled up to preferred
# terms. The listings themselves keep the terms as FAERS recorded them, because
# that is the traceability back to the report - but anything that counts or
# classifies those terms has to use the preferred terms the counts were built
# from. Otherwise a report counted under Homicide through the Lowest Level Term
# MURDER is classified as though murder were a co-reported symptom of itself:
# the index term is excluded by name, and MURDER is not that name.
pts_to_pt <- function(s, sep = "; ") {
  m <- llt_map()
  vapply(strsplit(as.character(s), sep, fixed = TRUE), function(v) {
    if (!length(v)) return("")
    v <- pt_std(v)
    h <- match(v, names(m))
    v[!is.na(h)] <- unname(m[h[!is.na(h)]])
    paste(sort(unique(v[nzchar(v)])), collapse = sep)
  }, character(1), USE.NAMES = FALSE)
}

# ── non-drug entries in the drug field ───────────────────────────────────────
# FAERS drug names are free text, and reporters use the field to say that there
# was nothing to report. "NO CONCURRENT MEDICATION" is not a medication, and
# counting it as one turns a monotherapy report into polypharmacy.
#
# data/non_drug_terms.csv lists the placeholders and what is done with each:
#
#   exclude  the report states there was no co-medication, so the entry is
#            dropped before drugs are counted
#   keep     a co-medication was given but not named ("UNKNOWN DRUG"). It is
#            counted as a co-drug, which keeps the report out of monopharmacy.
#            That is the conservative reading: it can only shrink the
#            monopharmacy set and weaken the monopharmacy signals.
#
# Matching is exact, on the whole name after upper-casing and stripping trailing
# punctuation (non_drug_key in R/groups.R), so a real drug whose name merely
# contains one of these words is unaffected: CANRENONE and MILRINONE stay.
.nondrug <- fread("data/non_drug_terms.csv")
stopifnot(all(c("term", "treatment") %chin% names(.nondrug)),
          all(.nondrug$treatment %chin% c("exclude", "keep")))
NON_DRUG_DROP <- toupper(trimws(.nondrug[treatment == "exclude", term]))

# Every placeholder, kept and excluded alike. The MGPS matrix uses this rather
# than NON_DRUG_DROP: "keep" is the conservative reading when counting a
# report's drugs, because an unnamed co-medication can only move a report out of
# monopharmacy. On the drug axis of the matrix it does the opposite - UNKNOWN
# becomes a pseudo-drug carrying a large share of reports, and the prior is
# fitted over that axis. A placeholder is not a drug, so it is not a row there.
NON_DRUG_ANY <- toupper(trimws(.nondrug$term))

# ── manual corrections ────────────────────────────────────────────────────────
# Monopharmacy reclassifications (brand/generic duplicates of the same drug).
mono_overrides <- fread("data/mono_overrides.csv", colClasses = c(primaryid = "character"))
mono_overrides[, primaryid := bit64::as.integer64(primaryid)]
# Manual deduplication blacklist (post hoc clinical review): see README.
# The removal list used in the article. Point FAERS_BLACKLIST at your own file
# to run the pipeline against a list you compiled yourself.
BLACKLIST_FILE <- Sys.getenv("FAERS_BLACKLIST", unset = "data/manual_dedup_blacklist.csv")

# ── code labels ───────────────────────────────────────────────────────────────
outc_labels <- c(OT = "Other serious", HO = "Hospitalization", DE = "Death",
                 LT = "Life-threatening", DS = "Disability",
                 RI = "Required intervention", CA = "Congenital anomaly")
rpsr_labels <- c(PHY = "Physician", PHR = "Pharmacist", OTH = "Other health professional",
                 LW = "Lawyer", CN = "Consumer", UF = "User facility", MFR = "Manufacturer",
                 DT = "Distributor", HP = "Health professional",
                 CSM = "Consumer/non-health professional", FGN = "Foreign",
                 LIT = "Literature", CR = "Consumer reporter", SDY = "Study")
country_codes <- c(US = "UNITED STATES", U.S. = "UNITED STATES", USA = "UNITED STATES",
                   U.S.A. = "UNITED STATES", GB = "UNITED KINGDOM", U.K. = "UNITED KINGDOM",
                   UK = "UNITED KINGDOM", CA = "CANADA", JP = "JAPAN", FR = "FRANCE",
                   DE = "GERMANY", AU = "AUSTRALIA", IT = "ITALY", ES = "SPAIN",
                   NL = "NETHERLANDS", SE = "SWEDEN", BE = "BELGIUM", CH = "SWITZERLAND",
                   KR = "SOUTH KOREA", TW = "TAIWAN", IN = "INDIA", BR = "BRAZIL")

# Reporter country as written in FAERS mapped to one spelling.
norm_country <- function(x) { x <- as.character(x); hit <- country_codes[x]; ifelse(is.na(hit), x, hit) }

# FAERS indication placeholders: entries that record the absence of a stated
# indication rather than an indication. Counting them as indications makes a
# report with nothing but "Product used for unknown indication" look as though
# its indication were known. 05_case_tables.R and tools/group_characteristics.R
# both drop them, from one list, so the Table 1 blocks and the SSRI group block
# that sit side by side in SI "GLOBALS" count the same way.
INDI_PLACEHOLDERS <- c("product used for unknown indication", "off label use",
                       "no suspected drug", "indication not provided",
                       "drug use for unknown indication")
