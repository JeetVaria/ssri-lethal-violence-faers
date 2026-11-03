# ==============================================================================
# R/setup.R - packages, folders, logging, small formatting helpers.
# Sourced by R/load_all.R; not meant to be run on its own.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(bit64)
  library(survival)
})

# faers and BiocParallel are needed only to build the database and compute the
# signals (steps 01-05, tools/group_characteristics.R). Everything after that
# reads the CSVs those steps wrote, so the figures, the time-to-onset tables
# and the Supplementary GLOBALS block regenerate from expected_2004_2025/
# with CRAN packages alone. Steps that do need them stop with need_faers().
HAVE_FAERS <- requireNamespace("faers", quietly = TRUE) &&
              requireNamespace("BiocParallel", quietly = TRUE)
if (HAVE_FAERS) suppressPackageStartupMessages({
  library(faers)         # Peng et al., Bioconductor
  library(BiocParallel)
})
need_faers <- function() {
  if (!HAVE_FAERS)
    stop("this step reads the FAERS database and needs the faers and ",
         "BiocParallel packages from Bioconductor:\n",
         '  BiocManager::install(c("faers", "BiocParallel"))\n',
         "Steps 06 and 07, tools/build_si_globals.R, combine_da.R and ",
         "reconcile.R do not need them.", call. = FALSE)
  invisible(TRUE)
}

# Case listings hold sorted lists inside a cell - every PT on a report, every
# co-drug. Base sort() follows the session's collation, so the same data gives
# "lorazepam; QUETIAPINE" under one locale and "QUETIAPINE; lorazepam" under
# another, and a replicator running in a container (where LC_ALL=C is usual)
# would see those columns differ for no reason. Collation is fixed here so the
# published strings are reproducible on any machine.
Sys.setlocale("LC_COLLATE", "C")

for (d in c(RES, CACHE_DIR, FIG_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

PHV_BPPARAM <- if (HAVE_FAERS) SerialParam(RNGseed = PHV_SEED) else NULL
YEARS       <- START_YEAR:END_YEAR

# ── logging ───────────────────────────────────────────────────────────────────
LOG_FILE <- file.path(RES, "pipeline_log.txt")
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "[%H:%M:%S] "), ...)
  message(msg)
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}
elapsed <- function(t0) paste0(round(as.numeric(Sys.time() - t0, units = "mins"), 1), " min")
gc2     <- function() { gc(verbose = FALSE); gc(verbose = FALSE); invisible() }

# ── file helpers ──────────────────────────────────────────────────────────────
res_path  <- function(...) file.path(RES, ...)
write_res <- function(dt, file) { fwrite(dt, res_path(file)); log_msg("  wrote ", file, " (", nrow(dt), " rows)") }
read_res  <- function(file, required = TRUE) {
  p <- res_path(file)
  if (!file.exists(p)) {
    if (required) stop(file, " not found in ", RES, " - run the earlier pipeline steps first.")
    return(NULL)
  }
  fread(p)
}
require_res <- function(file) {
  if (!file.exists(res_path(file))) stop(file, " not found in ", RES, " - run the earlier pipeline steps first.")
  invisible(res_path(file))
}

# Return the cached .rds result if present (and FORCE_REBUILD is off); otherwise
# evaluate `expr`, save it and return it. Used for the slow MCMC signal steps.
# Cached results record the methods and the arithmetic version that produced
# them (PHV_CACHE_KEY). Without that, a run made with FAERS_QUICK_TEST=TRUE
# (frequentist only) leaves .rds files that a later full run would silently
# reuse, dropping IC025 and EBGM05 from results that still look complete - and
# a change to how an estimate is computed would go unnoticed in the same way.
cached <- function(file, expr) {
  p <- res_path(file)
  if (file.exists(p) && !FORCE_REBUILD) {
    out  <- readRDS(p)
    made <- attr(out, "phv_methods")
    if (identical(made, PHV_CACHE_KEY)) { log_msg("  cached: ", file); return(out) }
    log_msg("  rebuilding ", file, ": cached under ",
            if (is.null(made)) "unrecorded" else paste(made, collapse = "/"),
            ", now ", paste(PHV_CACHE_KEY, collapse = "/"))
  }
  out <- expr
  attr(out, "phv_methods") <- PHV_CACHE_KEY
  saveRDS(out, p)
  out
}

# ── number formatting for summary tables ──────────────────────────────────────
# formatted from the double, not through as.integer(), which returns NA above
# 2^31 - 1: the drug-event matrix runs to billions of cells, and the two checks
# on its totals would otherwise report "NA cells" at exactly the scale where
# they are the only thing standing between a mis-built matrix and a published
# prior.
sfmt <- function(x) ifelse(is.na(x), "NA", formatC(as.numeric(x), format = "f", digits = 0, big.mark = ","))
spct <- function(a, b) sprintf("%.1f%%", a / max(b, 1) * 100)
np   <- function(n, d) paste0(sfmt(n), " (", spct(n, d), ")")
# FAERS dates are YYYYMMDD; partial dates are left missing rather than guessed.
pdate <- function(x) suppressWarnings(as.Date(as.character(x), format = "%Y%m%d"))
# a code rendered with its label, or left as the code when unmapped
lab  <- function(code, labels) ifelse(code %in% names(labels), labels[code], code)

snum <- function(x) tryCatch({ x <- as.numeric(x[[1]]); if (is.na(x)) "NA" else as.character(round(x, 1)) },
                             error = function(e) "NA")
siq  <- function(v) {
  if (!length(v) || all(is.na(v))) return("NA")
  paste0(snum(median(v, na.rm = TRUE)), " (", snum(quantile(v, .25, na.rm = TRUE)[[1]]), "-",
         snum(quantile(v, .75, na.rm = TRUE)[[1]]), ")")
}
smsd <- function(v) {
  if (!length(v) || all(is.na(v))) return("NA")
  paste0(snum(mean(v, na.rm = TRUE)), " (", snum(sd(v, na.rm = TRUE)), ")")
}

# Set operations on report IDs. FAERS IDs are 64-bit integers (bit64::integer64);
# base union/intersect/setdiff drop that class on older bit64 versions.
id_intersect <- function(x, y) unique(x[x %in% y])
id_setdiff   <- function(x, y) unique(x[!(x %in% y)])
id_union     <- function(x, y) unique(c(x, y))
