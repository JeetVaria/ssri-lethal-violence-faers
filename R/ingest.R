# ==============================================================================
# R/ingest.R - reading the FAERS quarterly ASCII files.
# ==============================================================================

# Path to one quarter's file (FAERS "faers_ascii_" or legacy AERS "aers_ascii_",
# zipped or unzipped, either case of the quarter letter); NULL if absent.
quarter_file <- function(yr, q, zipped_only = FALSE) {
  stems <- c(paste0("faers_ascii_", yr, c(toupper(q), q)), paste0("aers_ascii_", yr, c(toupper(q), q)))
  cands <- file.path(FAERS_DIR, c(paste0(stems, ".zip"), if (!zipped_only) stems))
  found <- cands[file.exists(cands)]
  if (length(found)) found[1] else NULL
}

parse_quarter <- function(yr, q) {
  f <- quarter_file(yr, q)
  if (is.null(f)) return(NULL)
  tryCatch(faers_parse(path = f, year = as.integer(yr), quarter = q),
           error = function(e) tryCatch(faers_parse(path = f),
             error = function(e2) { log_msg("  WARNING: cannot parse ", basename(f)); NULL }))
}

# Number of unique report IDs listed in the FDA's quarterly deletion files.
count_fda_deleted <- function(years) {
  ids <- character(0)
  for (yr in years) for (q in paste0("q", 1:4)) {
    zf <- quarter_file(yr, q, zipped_only = TRUE)
    if (is.null(zf)) next
    del <- tryCatch(unzip(zf, list = TRUE)$Name, error = function(e) character(0))
    del <- del[grepl("delet", del, ignore.case = TRUE) & grepl("[.]txt$", del, ignore.case = TRUE)]
    if (!length(del)) next
    tryCatch({
      unzip(zf, files = del[1], exdir = tempdir())
      x <- file.path(tempdir(), del[1])
      v <- trimws(fread(x, header = FALSE, fill = TRUE, colClasses = "character", showProgress = FALSE)[[1]])
      ids <- unique(c(ids, v[nchar(v) > 0]))
      unlink(x)
    }, error = function(e) NULL)
  }
  length(ids)
}

count_tables <- function(obj, label = "") {
  rbindlist(lapply(c("demo", "drug", "reac", "indi", "ther", "rpsr", "outc"), function(t)
    data.table(stage = label, table = t,
               n = as.integer(tryCatch(nrow(faers_get(obj, t)), error = function(e) NA_integer_)))))
}

std_cache_path <- function(yr) file.path(CACHE_DIR, paste0("std_", yr, ".rds"))
