# ==============================================================================
# 01_build_master.R - FAERS Q1 2004 - Q4 2025: parse, standardize, deduplicate.
#
# 1. Parse each quarter's ASCII files and standardize terms against MedDRA v29.0
#    (faers::faers_standardize). One cached .rds per year in CACHE_DIR.
# 2. Combine all years and deduplicate (faers::faers_dedup): collapses repeated
#    PRIMARYIDs, keeps the most recent version of each CASEID, and removes
#    cross-source duplicates.
# 3. Build the monopharmacy reference set (reports with a single distinct drug).
#
# Outputs (RES):
#   master_dedup.rds, master_mono.rds          FAERS objects used by later steps
#   SI_DATA_COUNTS.csv                    raw / pre- / post-dedup counts per year
#   qa_table_counts.csv                        row counts per FAERS table
#   removed_by_dedup.csv (+ _summary)          reports removed by deduplication
#     - the last two are written only when the quarterly parse runs, and are too
#       large to ship, so they are absent from expected_*/ by design
#   summary_master.csv                         reference-set characteristics
#                                              (Supplementary "COUNTS" / "GLOBALS")
# ==============================================================================
source("R/load_all.R")
need_faers()
t0 <- Sys.time()
log_msg("=== 01 BUILD MASTER: ", START_YEAR, "-", END_YEAR, " ===")
master_file <- res_path("master_dedup.rds")
rebuild <- !file.exists(master_file) || FORCE_REBUILD

# ── 1. per-year parse + standardize ───────────────────────────────────────────
# Only reached when the deduplicated database has to be built. Starting from an
# existing master_dedup.rds skips the 88-quarter parse altogether, and so needs
# neither the FAERS ASCII files nor a MedDRA license.
if (!rebuild) log_msg("  skipping the quarterly parse: master_dedup.rds is present")
yearly_qa <- if (!rebuild) NULL else rbindlist(lapply(YEARS, function(yr) {
  sp <- std_cache_path(yr)
  if (!file.exists(sp)) {
    log_msg("  parsing ", yr)
    qs <- Filter(Negate(is.null), lapply(paste0("q", 1:4), function(q) parse_quarter(yr, q)))
    if (!length(qs)) { log_msg("  no ASCII files for ", yr); return(data.table(year = yr, ascii_n = NA_integer_)) }
    raw <- if (length(qs) == 1) qs[[1]] else do.call(faers_combine, qs)
    rm(qs); gc2()
    std <- faers_standardize(raw, meddra_path = MEDDRA_DIR)
    rm(raw); gc2()
    saveRDS(std, sp)
  } else std <- readRDS(sp)
  n <- nrow(faers_get(std, "demo"))
  rm(std); gc2()
  log_msg("  ", yr, ": ", sfmt(n), " reports")
  data.table(year = yr, ascii_n = as.integer(n))
}))

if (rebuild) {
  # ── 2. combine + deduplicate ────────────────────────────────────────────────
  have <- YEARS[file.exists(vapply(YEARS, std_cache_path, ""))]
  master <- readRDS(std_cache_path(have[1]))
  for (yr in have[-1]) { master <- faers_combine(master, readRDS(std_cache_path(yr))); gc2() }
  pre_demo <- faers_get(master, "demo")[, .(primaryid, caseid, year, quarter, sex, age_in_years,
                                            reporter_country, event_dt, mfr_sndr, fda_dt)]
  pre_drug <- faers_get(master, "drug")[, .(n_drugs = .N, drugnames = paste(head(unique(drugname), 3), collapse = "; ")),
                                        by = primaryid]
  pre_counts <- count_tables(master, "pre_dedup")
  pre_by_yr  <- pre_demo[, .(pre_dedup_n = .N), by = year]
  fda_deleted_total <- count_fda_deleted(YEARS)

  # remove_deleted_cases is the package default; stated here so the analysis does
  # not change if that default ever does
  master_dedup <- faers_dedup(master, remove_deleted_cases = TRUE); rm(master); gc2()
  saveRDS(master_dedup, master_file)
  log_msg("Dedup: ", sfmt(nrow(pre_demo)), " -> ", sfmt(nrow(faers_get(master_dedup, "demo"))))

  # Reports removed by deduplication (QA list). This selects by PRIMARYID, so
  # it lists the reports that disappear entirely and not the superseded
  # versions of reports whose newest version was kept - those keep their
  # PRIMARYID in post_ids. The file therefore undercounts: in the two-year run
  # it holds 456,556 rows against the 456,682 the per-year flow removes. The
  # published counts do not come from here. summary_master.csv takes
  # dedup_removed_n from pre_dedup_n - post_dedup_n, which is complete.
  post_ids <- as.character(faers_get(master_dedup, "demo")$primaryid)
  removed  <- merge(pre_demo[!as.character(primaryid) %in% post_ids], pre_drug, by = "primaryid", all.x = TRUE)
  removed[, reason := "deduplication - superseded by newer version"]
  del_files <- list.files(FAERS_DIR, pattern = "deleted.*.txt", recursive = TRUE, ignore.case = TRUE, full.names = TRUE)
  del_ids <- unique(unlist(lapply(del_files, function(f)
    tryCatch(as.character(fread(f, header = FALSE)[[1]]), error = function(e) character(0)))))
  if (length(del_ids)) removed[as.character(primaryid) %in% del_ids, reason := "FDA deletion list"]
  write_res(removed, "removed_by_dedup.csv")
  write_res(removed[, .N, by = .(year, reason)][order(year)], "removed_by_dedup_summary.csv")
  rm(pre_demo, pre_drug, removed); gc2()

  # table counts and per-year flow
  post_counts <- count_tables(master_dedup, "post_dedup")
  tc <- merge(pre_counts[, .(table, pre_dedup_n = n)], post_counts[, .(table, post_dedup_n = n)], by = "table")
  tc[, `:=`(removed_n = pre_dedup_n - post_dedup_n, removed_pct = round((pre_dedup_n - post_dedup_n) / pre_dedup_n * 100, 1))]
  write_res(tc, "qa_table_counts.csv")
  yearly_qa <- merge(yearly_qa, pre_by_yr, by = "year", all.x = TRUE)
  yearly_qa <- merge(yearly_qa, faers_get(master_dedup, "demo")[, .(post_dedup_n = .N), by = year], by = "year", all.x = TRUE)
  yearly_qa[, `:=`(removed_dedup_n   = pre_dedup_n - post_dedup_n,
                   removed_dedup_pct = round((pre_dedup_n - post_dedup_n) / pre_dedup_n * 100, 1),
                   retention_pct     = round(post_dedup_n / pre_dedup_n * 100, 1))]
  write_res(yearly_qa[order(year)], "SI_DATA_COUNTS.csv")
  writeLines(as.character(fda_deleted_total), res_path("fda_deleted_total.txt"))
} else {
  log_msg("master_dedup.rds exists (set FORCE_REBUILD = TRUE to redo deduplication)")
  master_dedup <- readRDS(master_file)
  # Counted from the quarterly deletion files so the flow table can report it
  # when the deduplicated database is reused. If the quarterly files are not
  # present the count is unknowable here, and the table must say so rather than
  # record a zero.
  if (!file.exists(res_path("fda_deleted_total.txt"))) {
    n_del <- count_fda_deleted(YEARS)
    if (n_del > 0L) writeLines(as.character(n_del), res_path("fda_deleted_total.txt"))
    else log_msg("FDA deletion list not counted: quarterly files absent (reported as NA)")
  }
}

# ── 3. monopharmacy reference set ─────────────────────────────────────────────
# master_mono.rds is not covered by PHV_CACHE_KEY, so a change to
# mono_reference_ids() or to data/non_drug_terms.csv would otherwise be
# invisible: the next run would reuse the old reference set, and every
# monopharmacy number with it, without a word. The set is therefore written
# with a stamp beside it, and rebuilt whenever the stamp does not match. An
# unstamped file - the one published on Zenodo before this check existed -
# counts as not matching, and rebuilding it from master_dedup.rds needs no
# quarterly parse.
mono_file  <- res_path("master_mono.rds")
mono_sfile <- res_path("master_mono_stamp.txt")
mono_stamp <- paste0(START_YEAR, "-", END_YEAR, "/defs", PHV_DEF_HASH, "/rule8")
nodrug_file <- res_path("master_nodrug_ids.txt")
mono_ok <- !FORCE_REBUILD && file.exists(mono_file) && file.exists(mono_sfile) &&
  file.exists(nodrug_file) &&
  identical(readLines(mono_sfile, warn = FALSE)[1], mono_stamp)
if (!FORCE_REBUILD && file.exists(mono_file) && !mono_ok)
  log_msg("  master_mono.rds does not match this run (", mono_stamp, "): rebuilding")
if (mono_ok) {
  master_mono <- readRDS(mono_file)
  mono_ids <- faers_get(master_mono, "demo")$primaryid
  nodrug_ids <- bit64::as.integer64(readLines(nodrug_file, warn = FALSE))
} else {
  log_msg("  rebuilding the monopharmacy reference set (no quarterly parse needed)")
  # One report_drug_counts() pass over the whole database, the most expensive
  # thing in this script outside the parse. It yields both the reference set
  # and the reports naming no drug, which the summary block below needs; the
  # second list is written out so the cached path never has to repeat the pass.
  split <- mono_reference_ids(master_dedup)
  mono_ids <- split$mono_ids
  nodrug_ids <- split$nodrug_ids
  writeLines(as.character(nodrug_ids), nodrug_file)
  master_mono <- faers_keep(master_dedup, primaryid = mono_ids)
  # Written beside the target and renamed, never over it. On route 3 this file
  # IS the 1.18 GB download, and a save that runs out of memory or is
  # interrupted part way would otherwise leave neither the old file nor a new
  # one. The stamp is written last, so a crash between the two fails safe: the
  # set is rebuilt again next time rather than trusted.
  mono_tmp <- paste0(mono_file, ".tmp")
  saveRDS(master_mono, mono_tmp)
  # The rename is tried first and the old file deleted only if it fails, so on
  # POSIX the verified download is never removed before its replacement is in
  # place. Windows refuses to rename onto an existing file, hence the retry.
  if (!file.rename(mono_tmp, mono_file)) {
    unlink(mono_file)                       # Windows will not rename onto a file
    if (!file.rename(mono_tmp, mono_file))
      stop("could not put the rebuilt reference set at ", mono_file,
           ": it is in ", mono_tmp, ", rename it by hand.", call. = FALSE)
  }
  writeLines(mono_stamp, mono_sfile)
}
log_msg("Monopharmacy reference set: ", sfmt(length(mono_ids)), " (before manual removals)")

# ── reference-set summary (Supplementary "GLOBALS", first column) ─────────────
# Demographics describe the final analysis set: the deduplicated database after
# the manual blacklist removals, which is the denominator every percentage uses.
demo0 <- faers_get(master_dedup, "demo")
n_tot <- nrow(demo0)
bl    <- bit64::as.integer64(fread(BLACKLIST_FILE, colClasses = "character")$primaryid)
in_bl <- demo0$primaryid %in% bl
n_final     <- n_tot - sum(in_bl)
demo     <- demo0[!in_bl]
keep_ids <- demo$primaryid
n_mono_bl   <- sum(faers_get(master_mono, "demo")$primaryid %in% keep_ids)
# Reports naming no drug, which is not the same as reports with no row in the
# drug file: a report whose only entry is NO CONCURRENT MEDICATION has a row
# and names nothing. Counted by the rule that decides monopharmacy, so that the
# three rows below - monopharmacy, polypharmacy, no drug - are a clean split of
# the final count. Counting rows instead would leave every placeholder-only
# report in the polypharmacy residual.
n_nodrug_bl <- sum(nodrug_ids %in% keep_ids)
counts_file <- res_path("SI_DATA_COUNTS.csv")
if (!file.exists(counts_file))
  stop("SI_DATA_COUNTS.csv is not in ", RES, ".\n",
       "It records the per-year report flow, which can only be counted while the\n",
       "quarterly files are being parsed. Starting from an existing master_dedup.rds\n",
       "skips that parse, so copy SI_DATA_COUNTS.csv from the published results/\n",
       "folder in alongside the .rds files, or set FORCE_REBUILD <- TRUE in config.R\n",
       "to rebuild everything from the FAERS ASCII files.", call. = FALSE)
yq <- fread(counts_file)

age <- demo$age_in_years
age_ok <- !is.na(age) & age >= 0 & age <= 120
# as.numeric first: wt arrives as a 64-bit integer, which quantile() and mean()
# do not handle and silently return NA for
wt_raw <- suppressWarnings(as.numeric(demo$wt))
wt_raw <- wt_raw[is.finite(wt_raw) & wt_raw >= 1 & wt_raw <= 1000]
# wt carries no unit flag that can be relied on, so the unit is inferred: a
# median above 100 is read as pounds. Logged because the published median alone
# does not say which branch ran - 165 lb and 74.8 kg are the same number.
wt_in_lb <- length(wt_raw) > 100 && median(wt_raw) > 100
wt <- if (wt_in_lb) wt_raw * 0.453592 else wt_raw[wt_raw <= 500]
log_msg("  weight: ", if (wt_in_lb) "read as pounds and converted" else "read as kilograms",
        " | ", length(wt), " of ", length(wt_raw), " values kept")
ctry <- demo[, .N, by = reporter_country]
ctry[, reporter_country := norm_country(reporter_country)]
ctry <- ctry[, .(N = sum(N)), by = reporter_country][
  !is.na(reporter_country) & reporter_country != "" & !grepl("NOT SPECIFIED|UNKNOWN|UNK", reporter_country, ignore.case = TRUE)][order(-N)]
# counted as reports, not table rows: a report can carry the same code twice
outc <- faers_get(master_dedup, "outc")[primaryid %in% keep_ids][
          , .(N = uniqueN(primaryid)), by = outc_cod][order(-N)]
rpsr_all <- faers_get(master_dedup, "rpsr")[primaryid %in% keep_ids]
rpsr <- rpsr_all[, .(N = uniqueN(primaryid)), by = rpsr_cod][order(-N)]
n_rpsr_any <- uniqueN(rpsr_all$primaryid)
# Rows in each of the seven FAERS tables after deduplication, which is the
# "File / Count" block of Supplementary "COUNTS". They came straight off the
# parse and were never written out, so a reader had no way to check them
# against anything. Counted on master_dedup, before the manual blacklist, so
# DEMO here equals post_dedup_n above.
dedup_rows <- vapply(c("demo", "drug", "indi", "outc", "reac", "rpsr", "ther"),
                     function(tb) tryCatch(as.double(nrow(faers_get(master_dedup, tb))),
                                           error = function(e) NA_real_),
                     numeric(1))

kv   <- function(section, metric, value) data.table(section, metric, value = as.character(value))

summary_master <- rbindlist(list(
  kv("Pipeline", c("years", "raw_case_count", "dedup_removed_n", "post_dedup_n",
                   "manual_blacklist_removed_n", "final_count", "fda_deletion_list_n"),
     c(paste0(START_YEAR, "-", END_YEAR), sfmt(sum(yq$ascii_n, na.rm = TRUE)),
       sfmt(sum(yq$pre_dedup_n, na.rm = TRUE) - sum(yq$post_dedup_n, na.rm = TRUE)),
       sfmt(n_tot), sfmt(sum(in_bl)), sfmt(n_final),
       if (file.exists(res_path("fda_deleted_total.txt"))) readLines(res_path("fda_deleted_total.txt")) else "NA")),
  kv("Deduplicated database (rows per FAERS file)", toupper(names(dedup_rows)), sfmt(dedup_rows)),
  kv("Monopharmacy (final count)", c("n_monopharmacy", "n_polypharmacy", "n_no_drug_entry", "definition"),
     c(sfmt(n_mono_bl), sfmt(n_final - n_mono_bl - n_nodrug_bl), sfmt(n_nodrug_bl),
       "exactly one distinct drug named on the report, counted by canonical name, role ignored")),
  kv("Sex", c("female", "male", "unknown"),
     c(np(sum(tolower(demo$sex) == "female", na.rm = TRUE), n_final), np(sum(tolower(demo$sex) == "male", na.rm = TRUE), n_final),
       np(n_final - sum(tolower(demo$sex) %in% c("female", "male")), n_final))),
  kv("Age", c("recorded", "median (IQR)", "mean (SD)", "<18", "18-44", "45-64", "65-74", ">=75"),
     # bands are half-open and gated on age_ok, so they sum to the recorded count:
     # ages are not whole years, and an age outside 0-120 belongs to no band
     c(np(sum(age_ok), n_final), siq(age[age_ok]), smsd(age[age_ok]),
       sfmt(sum(age_ok & age < 18)), sfmt(sum(age_ok & age >= 18 & age < 45)),
       sfmt(sum(age_ok & age >= 45 & age < 65)), sfmt(sum(age_ok & age >= 65 & age < 75)),
       sfmt(sum(age_ok & age >= 75)))),
  kv("Weight (kg)", c("recorded", "median (IQR)", "mean (SD)"), c(np(length(wt), n_final), siq(wt), smsd(wt))),
  kv("Reporter country (top 5)", ctry$reporter_country[1:5], np(ctry$N[1:5], n_final)),
  kv("Outcomes (non-exclusive)", paste0(lab(outc$outc_cod, outc_labels), " (", outc$outc_cod, ")"),
     np(outc$N, n_final)),
  # most reports name no source at all; without this row the block looks short
  kv("Report source", c(paste0(lab(rpsr$rpsr_cod, rpsr_labels), " (", rpsr$rpsr_cod, ")"), "Not recorded"),
     c(np(rpsr$N, n_final), np(n_final - n_rpsr_any, n_final)))
))
write_res(summary_master, "summary_master.csv")
log_msg("=== 01 COMPLETE | ", elapsed(t0), " ===")
