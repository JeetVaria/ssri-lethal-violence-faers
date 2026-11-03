# ==============================================================================
# tools/group_characteristics.R - characteristics of the SSRI group as a whole.
#
# The per-PT case tables in step 05 describe the outcome case series; this
# describes the whole SSRI group behind them (the "GLOBALS" block of the
# Supplementary Information). It reads the saved group object rather than
# rebuilding anything, so it runs in seconds.
#
# Reads  results/ssri_group.rds (step 02), results/ssri_codrug_per_report.csv
# Writes results/ssri_group_characteristics.csv - the "SSRI group" block of SI GLOBALS
#
# Usage: Rscript tools/group_characteristics.R
# ==============================================================================
source("R/load_all.R")
need_faers()

obj  <- readRDS(require_res("ssri_group.rds"))
demo <- faers_get(obj, "demo")
N    <- nrow(demo)
cd   <- fread(require_res("ssri_codrug_per_report.csv"))
# The same three-way split as Table 1 and the reference set, each side counted
# directly. Deriving polypharmacy as N - n_mono would file a report naming no
# drug under polypharmacy, and this block is concatenated with Table 1 and
# summary_master.csv into one Supplementary sheet.
n_mono   <- cd[n_codrugs == 0 & n_distinct_drugs >= 1L, .N]
n_nodrug <- cd[n_distinct_drugs == 0L, .N]
n_poly   <- cd[n_codrugs > 0, .N]

kv <- function(metric, value) data.table(metric, value = as.character(value))
hd <- function(metric) kv(metric, "-")
ind <- function(x) paste0("  ", x)

# per drug, by the same dictionary the group filter uses
drugs <- rbindlist(lapply(names(ssri_drugs), function(nm)
  data.table(drug = tools::toTitleCase(nm),
             n = nrow(faers_get(drug_filter(obj, ssri_drug_pattern[[nm]]), "demo")))))

sex  <- tolower(demo$sex)
age  <- suppressWarnings(as.numeric(demo$age_in_years))
ok   <- !is.na(age) & age >= 0 & age <= 120
ctry <- demo[, .N, by = reporter_country][, reporter_country := tools::toTitleCase(tolower(norm_country(reporter_country)))][
  , .(N = sum(N)), by = reporter_country][
  !is.na(reporter_country) & reporter_country != "" &
  !grepl("NOT SPECIFIED|UNKNOWN|UNK", reporter_country, ignore.case = TRUE)][order(-N)]
outc <- faers_get(obj, "outc")[, .(N = uniqueN(primaryid)), by = outc_cod][order(-N)]
rpsr <- faers_get(obj, "rpsr")[, .(N = uniqueN(primaryid)), by = rpsr_cod][order(-N)]

# indications: one count per report per distinct indication. The placeholder
# entries are dropped on the same list the case tables use, so this block and
# the Table 1 blocks beside it in SI "GLOBALS" count an indication the same way.
it   <- faers_get(obj, "indi")[!is.na(indi_pt) & indi_pt != "" &
                               !tolower(indi_pt) %in% INDI_PLACEHOLDERS]
it[, indication := title_case(tolower(indi_pt))]
n_indi <- uniqueN(it$primaryid)
indi <- unique(it[, .(primaryid, indication)])[, .(N = .N), by = indication][order(-N)][1:10]

out <- rbindlist(list(
  kv("SSRI group", N),
  kv(ind("Monopharmacy (0 co-drugs) Case Count"), n_mono),
  kv(ind("Polypharmacy (>=1 co-drug) Case Count"), n_poly),
  if (n_nodrug > 0L) kv(ind("No drug named Case Count"), n_nodrug) else NULL,
  hd("SSRI drug"), kv(ind(drugs$drug), drugs$n),
  hd("Sex"),
  kv(ind("Female"), sum(sex == "female", na.rm = TRUE)),
  kv(ind("Male"), sum(sex == "male", na.rm = TRUE)),
  kv(ind("Unknown"), N - sum(sex %in% c("female", "male"))),
  hd("Age"),
  kv(ind("Median age, years (IQR)"), siq(age[ok])),
  kv(ind("Mean age, years +/- SD"), sub("\\((.*)\\)", "+/- \\1", smsd(age[ok]))),
  kv("Age recorded", sum(ok)),
  kv(ind("<18 years"), sum(ok & age < 18)),
  kv(ind("18-44 years"), sum(ok & age >= 18 & age < 45)),
  kv(ind("45-64 years"), sum(ok & age >= 45 & age < 65)),
  kv(ind("65-74 years"), sum(ok & age >= 65 & age < 75)),
  kv(ind(">=75 years"), sum(ok & age >= 75)),
  hd("Top 5 reporting countries"), kv(ind(ctry$reporter_country[1:5]), ctry$N[1:5]),
  hd("Reported outcomes"), kv(ind(lab(outc$outc_cod, outc_labels)), outc$N),
  hd("Report source (rpsr_cod; a report may carry several)"), kv(ind(lab(rpsr$rpsr_cod, rpsr_labels)), rpsr$N),
  kv(ind("Not recorded"), N - uniqueN(faers_get(obj, "rpsr")$primaryid)),
  hd("Top 10 indications for use"), kv(ind(indi$indication), indi$N),
  kv(ind("No indication recorded"), N - n_indi)
))
write_res(out, "ssri_group_characteristics.csv")
log_msg("SSRI group characteristics: ", nrow(out), " rows | group ", sfmt(N),
        " | monopharmacy ", sfmt(n_mono))
