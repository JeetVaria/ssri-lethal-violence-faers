# ==============================================================================
# 08_rxnorm_lookups.R - rebuild the RxNorm co-drug lookups used by step 05.
#
# OPTIONAL. Step 05 reads the frozen lookups in data/ (rxnorm_name_map.csv,
# rxnorm_class_map.csv), which reproduce the published co-drug class counts.
# This step shows how those files were made, and lets anyone rebuild them.
# RxNorm is updated monthly, so a fresh rebuild can differ slightly from the
# frozen files; the differences are listed in rxnorm_lookup_comparison.csv.
#
# Input:   the case lists written by step 05 (Homicide, Completed Suicide,
#          Serotonin Syndrome; H-S cases are a subset of these).
# Method:  1. Clean each co-drug name (role code, trailing clarifier, "/nnnnnn/"
#             suffix, "NOS" removed; ETHANOL -> ALCOHOL).
#          2. Resolve it to an RxNorm ingredient with the rxref package (NLM RxNav
#             APIs): find_ingredients(); if that fails, resolve() the brand name
#             and take its ingredient with ingredients_for_rxcui(). A precise
#             ingredient (PIN, e.g. a salt) with no ATC code takes the ATC of its
#             base ingredient.
#          3. Fetch ATC (get_atc) and VA (get_va) classes for each ingredient.
#          4. Name map: names resolving to the same RxCUI share one canonical
#             name. Low-confidence matches (string similarity < 0.5 to the
#             ingredient, or brand-match score < 50, unless the same RxCUI is
#             also a confident direct match) keep their own name.
#          5. Class map: each canonical name gets one class, in this order:
#             substance list (alcohol, cocaine, cannabis) -> ATC and VA class
#             agree -> keyword list (R/codrugs.R) -> any ATC class -> VA class
#             -> "Other".
#
# Usage:   Rscript analysis/08_rxnorm_lookups.R            (needs internet; about
#          1.5 h for ~4,000 names; resumable - finished names are cached)
#          RXNORM_RESOLUTION=path/to/resolution.csv Rscript analysis/08_rxnorm_lookups.R
#          (offline: skips steps 1-3 and rebuilds the maps from a saved
#          resolution table, e.g. results/rxnorm/rxnorm_resolution.csv)
#
# Outputs (results/rxnorm/): rxnorm_resolution.csv, rxnorm_name_map.csv,
#          rxnorm_class_map.csv, rxnorm_lookup_comparison.csv
# To use the rebuilt maps in step 05, copy them into data/ and re-run 05.
# ==============================================================================
source("R/load_all.R")
t0 <- Sys.time()
log_msg("=== 08 RXNORM LOOKUPS ===")

OUT <- res_path("rxnorm"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
out_file <- function(f) file.path(OUT, f)

CASE_FILES <- c(Homicide             = "cases_homicide_with_indications.csv",
                `Serotonin Syndrome` = "cases_serotonin_syndrome_with_indications.csv",
                `Completed Suicide`  = "cases_completed_suicide_with_indications.csv")

APPROX_MAX_ENTRIES       <- 5
BATCH_SIZE               <- 25
LOW_CONF_THRESHOLD       <- 0.5
BRAND_FALLBACK_MIN_SCORE <- 50
SSRI_NODE_SUFFIX         <- " (Co-reported index SSRI)"

# ---- 1. cleaning ---------------------------------------------------------------
clean_rx_tokens <- function(v) {
  v <- v[!is.na(v) & v != "" & v != "(None Recorded)"]
  if (!length(v)) return(character(0))
  v <- sub("\\[[^]]*\\]\\s*$", "", v)              # trailing [ROLE]
  v <- toupper(trimws(v))
  v <- sub("\\.$", "", v)                          # trailing period
  v <- sub("\\s*\\([^)]*\\)\\s*$", "", v)          # trailing "(clarifier)"
  v <- sub("\\s*/\\s*[0-9]+\\s*/\\s*$", "", v)     # trailing "/00139501/"
  v <- sub("\\s+NOS$", "", v)
  v <- trimws(v)
  v <- ifelse(v == "ETHANOL", "ALCOHOL", v)
  v[nzchar(v)]
}

# ---- 2-3. resolution (rxref) -----------------------------------------------------
pick_col <- function(df, candidates) { hit <- candidates[candidates %in% names(df)]; if (length(hit)) hit[1] else NA_character_ }

resolve_ingredients <- function(nms) {
  empty <- data.table(input = character(0), rxcui = character(0), matched_name = character(0),
                      score = numeric(0), tty = character(0))
  ing <- tryCatch(rxref::find_ingredients(nms, max_entries = APPROX_MAX_ENTRIES, show_progress = FALSE),
                  error = function(e) { log_msg("    find_ingredients() error: ", conditionMessage(e)); NULL })
  if (is.null(ing) || !nrow(ing)) return(empty)
  ing <- as.data.table(ing)
  if (!all(c("input", "rxcui") %in% names(ing))) return(empty)
  ing[, `:=`(rxcui = as.character(rxcui),
             score = suppressWarnings(as.numeric(if ("score" %in% names(ing)) score else NA_real_)))]
  setorder(ing, input, -score, na.last = TRUE)
  best <- ing[, .SD[1], by = input]
  data.table(input = best$input, rxcui = best$rxcui,
             matched_name = if ("name" %in% names(best)) best$name else NA_character_,
             score = best$score,
             tty = if ("tty" %in% names(best)) toupper(as.character(best$tty)) else NA_character_)
}

resolve_brands <- function(nms) {
  empty <- data.table(input = character(0), rxcui = character(0), matched_name = character(0),
                      score = numeric(0), tty = character(0), note = character(0))
  nms <- unique(nms[!is.na(nms) & nzchar(nms)])
  if (!length(nms)) return(empty)
  res <- tryCatch(rxref::resolve(nms, type = "name", max_entries = 1, show_progress = FALSE), error = function(e) NULL)
  if (is.null(res) || !nrow(res)) return(empty)
  res <- as.data.table(res)
  ci <- pick_col(res, "input"); cr <- pick_col(res, "rxcui")
  if (is.na(ci) || is.na(cr)) return(empty)
  cn <- pick_col(res, "name"); ct <- pick_col(res, "tty"); cs <- pick_col(res, "score")
  res <- data.table(input = res[[ci]], rxcui = as.character(res[[cr]]),
                    name  = if (is.na(cn)) NA_character_ else as.character(res[[cn]]),
                    tty   = if (is.na(ct)) NA_character_ else toupper(as.character(res[[ct]])),
                    score = if (is.na(cs)) NA_real_ else suppressWarnings(as.numeric(res[[cs]])))
  res <- res[!is.na(rxcui) & nzchar(rxcui)]
  if (!nrow(res)) return(empty)
  out <- list()
  direct <- res[!is.na(tty) & tty %in% c("IN", "PIN")]
  if (nrow(direct)) out[[1]] <- data.table(input = direct$input, rxcui = direct$rxcui, matched_name = direct$name,
                                           score = direct$score, tty = direct$tty,
                                           note = "resolve() matched an ingredient-level concept directly")
  walk <- res[is.na(tty) | !tty %in% c("IN", "PIN")]
  if (nrow(walk)) {
    ing <- tryCatch(rxref::ingredients_for_rxcui(unique(walk$rxcui), include_pin = TRUE, include_min = FALSE,
                                                 show_progress = FALSE), error = function(e) NULL)
    ing <- if (!is.null(ing)) as.data.table(ing) else NULL
    ok  <- !is.null(ing) && nrow(ing) && all(c("rxcui", "ingredient_rxcui", "ingredient_name") %in% names(ing))
    if (ok) {
      ing[, rxcui := as.character(rxcui)]
      ing <- unique(ing[!is.na(ingredient_rxcui) & nzchar(as.character(ingredient_rxcui))])
      if ("ingredient_tty" %in% names(ing)) ing[, rk := ifelse(toupper(ingredient_tty) == "IN", 1L, 2L)] else ing[, rk := 1L]
      setorder(ing, rxcui, rk, ingredient_name)
      n_ing <- ing[, .(n_ing = uniqueN(ingredient_rxcui)), by = rxcui]
      best  <- merge(ing[, .SD[1], by = rxcui], n_ing, by = "rxcui")
    }
    for (j in seq_len(nrow(walk))) {
      b <- if (ok) best[rxcui == walk$rxcui[j]] else NULL
      if (is.null(b) || !nrow(b)) {
        out[[length(out) + 1]] <- data.table(input = walk$input[j], rxcui = NA_character_, matched_name = NA_character_,
          score = NA_real_, tty = NA_character_,
          note = paste0("resolve() matched '", walk$name[j], "' (tty=", walk$tty[j], ") but no ingredient was found"))
        next
      }
      out[[length(out) + 1]] <- data.table(input = walk$input[j], rxcui = as.character(b$ingredient_rxcui[1]),
        matched_name = as.character(b$ingredient_name[1]), score = walk$score[j],
        tty = if ("ingredient_tty" %in% names(b)) toupper(as.character(b$ingredient_tty[1])) else NA_character_,
        note = paste0("brand fallback: '", walk$name[j], "' (tty=", walk$tty[j], ") -> ",
                      if (b$n_ing[1] > 1) paste0(b$n_ing[1], " ingredients, showing first (", b$ingredient_name[1], ")")
                      else b$ingredient_name[1]))
    }
  }
  if (length(out)) rbindlist(out, fill = TRUE) else empty
}

get_classes <- function(rxcuis, fetch_fn, extra_args = list()) {
  empty <- data.table(rxcui = character(0), class_id = character(0), class_name = character(0))
  rxcuis <- unique(rxcuis[!is.na(rxcuis) & nzchar(rxcuis)])
  if (!length(rxcuis)) return(empty)
  res <- tryCatch(do.call(fetch_fn, c(list(rxcuis, by = "rxcui", keep_input = TRUE, show_progress = FALSE), extra_args)),
                  error = function(e) NULL)
  if (is.null(res) || !nrow(res)) return(empty)
  res <- as.data.table(res)
  cr <- pick_col(res, c("input", "rxcui")); ci <- pick_col(res, c("class_id", "classId", "atc_code", "code"))
  cn <- pick_col(res, c("class_name", "className", "atc_class_name", "class_label", "name"))
  if (is.na(cr) || is.na(ci)) return(empty)
  out <- data.table(rxcui = as.character(res[[cr]]), class_id = as.character(res[[ci]]),
                    class_name = if (is.na(cn)) NA_character_ else as.character(res[[cn]]))
  unique(out[!is.na(class_id) & nzchar(class_id)])
}
get_atc_classes <- function(rx) get_classes(rx, rxref::get_atc)
get_va_classes  <- function(rx) get_classes(rx, rxref::get_va, list(extended = TRUE))

pin_to_ingredient <- function(pins) {
  rows <- lapply(unique(pins[!is.na(pins) & nzchar(pins)]), function(pr) {
    ing <- tryCatch(rxref::ingredients_for_rxcui(pr, include_pin = TRUE, include_min = FALSE, show_progress = FALSE),
                    error = function(e) NULL)
    if (is.null(ing) || !nrow(ing)) return(NULL)
    ing <- as.data.table(ing)
    if (!all(c("ingredient_rxcui", "ingredient_name") %in% names(ing))) return(NULL)
    inr <- if ("ingredient_tty" %in% names(ing) && any(toupper(ing$ingredient_tty) == "IN")) ing[toupper(ingredient_tty) == "IN"] else ing
    inr <- inr[!is.na(ingredient_rxcui) & nzchar(as.character(ingredient_rxcui))]
    if (!nrow(inr)) return(NULL)
    data.table(pin_rxcui = pr, in_rxcui = as.character(inr$ingredient_rxcui[1]), in_name = as.character(inr$ingredient_name[1]))
  })
  rows <- rbindlist(rows)
  if (!nrow(rows)) data.table(pin_rxcui = character(0), in_rxcui = character(0), in_name = character(0)) else rows
}

match_similarity <- function(raw, matched) {
  if (is.na(raw) || is.na(matched) || !nzchar(raw) || !nzchar(matched)) return(NA_real_)
  a <- gsub("[^A-Z0-9]", "", toupper(raw)); b <- gsub("[^A-Z0-9]", "", toupper(matched))
  if (!nzchar(a) || !nzchar(b)) return(NA_real_)
  1 - adist(a, b)[1, 1] / max(nchar(a), nchar(b))
}

resolution_row <- function(nm, n, rxcui = NA_character_, matched = NA_character_, score = NA_real_,
                           method = "unresolved", atc = NULL, va = NULL, note = NA_character_) {
  atc <- if (is.null(atc) || !nrow(atc)) NULL else atc[order(class_id)]
  va  <- if (is.null(va)  || !nrow(va))  NULL else va[order(class_id)]
  sim <- if (identical(method, "direct_ingredient")) match_similarity(nm, matched) else NA_real_
  lc  <- if (!is.na(sim)) sim < LOW_CONF_THRESHOLD else
         if (identical(method, "brand_fallback") && !is.na(score)) score < BRAND_FALLBACK_MIN_SCORE else NA
  data.table(raw_name = nm, occurrence_count = n, rxcui = rxcui, rxnorm_matched_name = matched, match_score = score,
             resolution_method = method, low_confidence_match = lc, match_similarity = sim,
             atc_code_primary   = if (is.null(atc)) NA_character_ else atc$class_id[1],
             atc_code_all       = if (is.null(atc)) NA_character_ else paste(atc$class_id, collapse = "; "),
             atc_class_name_all = if (is.null(atc)) NA_character_ else paste(atc$class_name, collapse = "; "),
             va_class_name_all  = if (is.null(va))  NA_character_ else paste(va$class_name, collapse = "; "),
             resolution_note = note)
}

resolve_names <- function(nms, counts) {
  ing <- resolve_ingredients(nms)
  ii  <- match(nms, ing$input)
  unres <- is.na(ii) | is.na(ing$rxcui[ii])
  fb  <- resolve_brands(unique(nms[unres]))
  fi  <- match(nms, fb$input)
  r <- rbindlist(lapply(seq_along(nms), function(j) {
    if (!unres[j]) return(data.table(name = nms[j], rxcui = ing$rxcui[ii[j]], matched = ing$matched_name[ii[j]],
                                     score = ing$score[ii[j]], tty = ing$tty[ii[j]], method = "direct_ingredient", note = NA_character_))
    if (!is.na(fi[j]) && !is.na(fb$rxcui[fi[j]]))
      return(data.table(name = nms[j], rxcui = fb$rxcui[fi[j]], matched = fb$matched_name[fi[j]],
                        score = fb$score[fi[j]], tty = fb$tty[fi[j]], method = "brand_fallback", note = fb$note[fi[j]]))
    data.table(name = nms[j], rxcui = NA_character_, matched = NA_character_, score = NA_real_, tty = NA_character_,
               method = "unresolved", note = if (!is.na(fi[j])) fb$note[fi[j]] else "no RxNorm ingredient match found")
  }))
  rx  <- unique(na.omit(r$rxcui))
  atc <- get_atc_classes(rx); va <- get_va_classes(rx)
  pins <- unique(r[!is.na(rxcui) & tty %in% "PIN" & !(rxcui %in% atc$rxcui), rxcui])
  walk <- pin_to_ingredient(pins)
  atc_walk <- if (nrow(walk)) get_atc_classes(unique(walk$in_rxcui)) else atc[0]
  rbindlist(lapply(seq_len(nrow(r)), function(j) {
    x <- r[j]; n <- unname(counts[x$name])
    if (is.na(x$rxcui)) return(resolution_row(x$name, n, method = x$method, note = x$note))
    a <- atc[rxcui == x$rxcui]; v <- va[rxcui == x$rxcui]; note <- x$note
    if (!nrow(a) && identical(x$tty, "PIN")) {
      w <- walk[pin_rxcui == x$rxcui]
      if (nrow(w)) {
        a <- atc_walk[rxcui == w$in_rxcui[1]]
        if (nrow(a)) note <- paste0(if (!is.na(note)) paste0(note, "; ") else "",
                                    "ATC from parent ingredient ", w$in_name[1], " (RxCUI ", w$in_rxcui[1], ")")
      }
    }
    if (!nrow(a)) note <- paste0("RxCUI resolved but no ATC class", if (!is.na(note)) paste0(" -- ", note) else "")
    resolution_row(x$name, n, x$rxcui, x$matched, x$score, x$method, a, v, note)
  }), fill = TRUE)
}

read_resolution <- function(path) {
  d <- fread(path, na.strings = c("NA", ""), colClasses = list(character = "rxcui"))
  need <- c("raw_name", "rxcui", "rxnorm_matched_name", "resolution_method", "low_confidence_match",
            "atc_code_primary", "atc_code_all", "va_class_name_all")
  miss <- setdiff(need, names(d))
  if (length(miss)) stop("resolution table is missing column(s): ", paste(miss, collapse = ", "))
  d
}

# ---- run steps 1-3, or load a saved resolution table ---------------------------
res_file <- Sys.getenv("RXNORM_RESOLUTION", unset = "")
if (nzchar(res_file)) {
  log_msg("  using saved resolution table: ", res_file, " (no RxNav calls)")
  lookup <- read_resolution(res_file)
} else {
  if (!requireNamespace("rxref", quietly = TRUE))
    stop("Package 'rxref' is needed to query RxNav: install.packages(\"rxref\")")
  all_names <- unlist(lapply(names(CASE_FILES), function(cl) {
    d <- read_res(CASE_FILES[[cl]])
    unlist(lapply(strsplit(as.character(d$codrug_names), "; ", fixed = TRUE),
                  function(v) unique(clean_rx_tokens(v))), use.names = FALSE)
  }), use.names = FALSE)
  counts <- sort(table(all_names), decreasing = TRUE)
  counts <- setNames(as.integer(counts), names(counts))
  log_msg("  ", length(counts), " distinct co-drug names (", length(all_names), " mentions)")

  cache_path <- out_file("rxnorm_resolution.csv")
  done <- if (file.exists(cache_path)) read_resolution(cache_path) else NULL
  todo <- setdiff(names(counts), done$raw_name)
  log_msg("  ", length(done$raw_name), " cached, ", length(todo), " to resolve")
  results <- list(done)
  for (b in seq_len(ceiling(length(todo) / BATCH_SIZE))) {
    nms <- todo[((b - 1L) * BATCH_SIZE + 1L):min(b * BATCH_SIZE, length(todo))]
    rows <- tryCatch(resolve_names(nms, counts),
                     error = function(e) { log_msg("    batch ", b, " failed (re-run to retry): ", conditionMessage(e)); NULL })
    results[[length(results) + 1]] <- rows
    tmp <- paste0(cache_path, ".tmp")
    fwrite(rbindlist(results, fill = TRUE)[order(-occurrence_count, raw_name)], tmp)
    file.rename(tmp, cache_path)
    if (b %% 10 == 0) log_msg("    ", min(b * BATCH_SIZE, length(todo)), "/", length(todo), " resolved")
  }
  log_msg("  wrote rxnorm_resolution.csv")
  lookup <- read_resolution(cache_path)
}
log_msg("  resolution: ", nrow(lookup), " names; ", sum(!is.na(lookup$rxcui)), " to an RxCUI; ",
        sum(!is.na(lookup$atc_code_primary)), " with an ATC class")

# A brand-fallback match flagged low-confidence is accepted when the same RxCUI is
# also a confident direct ingredient match elsewhere in the table.
confident <- lookup$resolution_method == "direct_ingredient" & lookup$low_confidence_match %in% FALSE & !is.na(lookup$rxcui)
corrob <- lookup$resolution_method == "brand_fallback" & lookup$low_confidence_match %in% TRUE &
          !is.na(lookup$rxcui) & lookup$rxcui %in% unique(lookup$rxcui[confident])
lookup[corrob, low_confidence_match := FALSE]
log_msg("  ", sum(corrob), " low-confidence brand matches confirmed by a direct match")

# ---- 4. name map -----------------------------------------------------------------
low  <- lookup$low_confidence_match %in% TRUE
good <- !is.na(lookup$rxcui) & nzchar(lookup$rxcui) & !low
name_map <- data.table(raw_name = lookup$raw_name, rxcui = lookup$rxcui,
                       canonical_name = ifelse(good, toupper(lookup$rxnorm_matched_name), lookup$raw_name),
                       canonical_key  = ifelse(good, lookup$rxcui,
                                               ifelse(low, paste0("LOW_CONFIDENCE::", lookup$raw_name),
                                                      paste0("UNRESOLVED::", lookup$raw_name))))

# ---- 5. class map ----------------------------------------------------------------
ATC_PREFIX_TO_CLASS <- c(
  "^N05AN" = "Mood Stabilizers", "^N05CF" = "Z-drugs (Sedative-Hypnotics)", "^N05CD" = "Anti-Anxiety Agents",
  "^N06A"  = "Antidepressants",  "^N06B"  = "Stimulants",                   "^N05B"  = "Anti-Anxiety Agents",
  "^N05A"  = "Antipsychotics",   "^N03A"  = "Mood Stabilizers",             "^N02"   = "Analgesics",
  "^M01A"  = "Analgesics",       "^B01"   = "Cardiovascular Drugs",         "^C"     = "Cardiovascular Drugs",
  "^J01"   = "Antibiotics",      "^R06"   = "Antihistamines",               "^A10"   = "Antidiabetics",
  "^H02"   = "Hormonal Therapies", "^H03" = "Hormonal Therapies",           "^G03"   = "Hormonal Therapies",
  "^J05"   = "Antivirals",       "^A02"   = "Gastrointestinal Drugs",       "^A03"   = "Gastrointestinal Drugs",
  "^A04"   = "Gastrointestinal Drugs", "^A07" = "Gastrointestinal Drugs")
VA_NAME_TO_CLASS <- c(                      # substring of the VA class name -> class (first match wins)
  "ANTIDEPRESSANT" = "Antidepressants", "ANXIOLYTIC" = "Anti-Anxiety Agents",
  "BENZODIAZEPINE DERIVATIVE SEDATIVE" = "Anti-Anxiety Agents", "ANTIPSYCHOTIC" = "Antipsychotics",
  "MOOD STABILIZ" = "Mood Stabilizers", "ANTICONVULSANT" = "Mood Stabilizers", "CNS STIMULANT" = "Stimulants",
  "ANALGESIC" = "Analgesics", "OPIOID" = "Analgesics", "NSAID" = "Analgesics",
  "ANTIHYPERTENSIVE" = "Cardiovascular Drugs", "ANTILIPEMIC" = "Cardiovascular Drugs",
  "ANTICOAGULANT" = "Cardiovascular Drugs", "CARDIA" = "Cardiovascular Drugs",
  "ANTIBACTERIAL,SYSTEMIC" = "Antibiotics", "ANTIHISTAMINE" = "Antihistamines",
  "ANTIDIABETIC" = "Antidiabetics", "HYPOGLYCEMIC" = "Antidiabetics", "THYROID" = "Hormonal Therapies",
  "CORTICOSTEROID" = "Hormonal Therapies", "ESTROGEN" = "Hormonal Therapies", "ANDROGEN" = "Hormonal Therapies",
  "ANTIVIRAL" = "Antivirals", "ANTIEMETIC" = "Gastrointestinal Drugs", "GASTROINTESTINAL" = "Gastrointestinal Drugs",
  "PROTON PUMP" = "Gastrointestinal Drugs")

atc_class <- function(code) { for (p in names(ATC_PREFIX_TO_CLASS)) if (grepl(p, code)) return(ATC_PREFIX_TO_CLASS[[p]]); NA_character_ }
va_class  <- function(nm)   { for (p in names(VA_NAME_TO_CLASS)) if (grepl(p, nm, fixed = TRUE)) return(VA_NAME_TO_CLASS[[p]]); NA_character_ }
first_class <- function(all, fn) {
  if (is.na(all) || !nzchar(all)) return(NA_character_)
  for (x in trimws(strsplit(all, "; ", fixed = TRUE)[[1]])) { cl <- fn(x); if (!is.na(cl)) return(cl) }
  NA_character_
}
has_keyword <- function(name_uc, kws) any(vapply(toupper(kws), function(k) grepl(k, name_uc, fixed = TRUE), logical(1)))

# rx_class() is called with the canonical ingredient names, so key these tables
# by the same name. Keyed on raw_name, a brand that resolved to a differently
# spelled ingredient filed its ATC under the brand, the lookup never found it,
# and the ingredient fell through to "Other": MEVALOTIN (pravastatin) was classed
# as unclassified rather than cardiovascular.
key  <- toupper(name_map$canonical_name)
# several recorded names can share one ingredient: keep the row that has a class
ord  <- order(is.na(lookup$atc_code_primary), is.na(lookup$va_class_name_all))
keep <- ord[!duplicated(key[ord])]
ATC_PRIMARY <- setNames(lookup$atc_code_primary[keep],  key[keep])
ATC_ALL     <- setNames(lookup$atc_code_all[keep],      key[keep])
VA_ALL      <- setNames(lookup$va_class_name_all[keep], key[keep])

rx_class <- function(name) {
  nu <- toupper(name)
  for (cl in names(substance_use_classes)) if (has_keyword(nu, substance_use_classes[[cl]])) return(cl)
  ap <- ATC_PRIMARY[nu]
  a1 <- if (!is.na(ap) && nzchar(ap)) atc_class(ap) else NA_character_
  vc <- first_class(VA_ALL[nu], va_class)
  if (!is.na(a1) && !is.na(vc) && identical(a1, vc)) return(a1)
  for (cl in names(codrug_classes)) if (has_keyword(nu, codrug_classes[[cl]])) return(cl)
  ac <- first_class(ATC_ALL[nu], atc_class)
  if (!is.na(ac)) return(ac)
  if (!is.na(vc)) return(vc)
  "Other"
}

ssri_nodes <- paste0(toupper(ssri_levels), SSRI_NODE_SUFFIX)   # co-reported index SSRIs, one entry each
base  <- sort(unique(c(name_map$canonical_name, ssri_nodes)))
cls   <- vapply(base, rx_class, character(1), USE.NAMES = FALSE)
label <- ifelse(base %in% ssri_nodes, base, paste0(base, " (", cls, ")"))
dup   <- duplicated(label) | duplicated(label, fromLast = TRUE)
label[dup] <- base[dup]
class_map <- data.table(canonical_name = label, assigned_class = cls)

fwrite(name_map,  out_file("rxnorm_name_map.csv"));  log_msg("  wrote rxnorm/rxnorm_name_map.csv (", nrow(name_map), " rows)")
fwrite(class_map, out_file("rxnorm_class_map.csv")); log_msg("  wrote rxnorm/rxnorm_class_map.csv (", nrow(class_map), " rows)")

# ---- compare with the frozen lookups in data/ ------------------------------------
cmp <- list()
if (file.exists("data/rxnorm_name_map.csv")) {
  old <- fread("data/rxnorm_name_map.csv", colClasses = list(character = "rxcui"))
  m <- suppressWarnings(merge(old[, .(raw_name, canonical_old = canonical_name)], name_map[, .(raw_name, canonical_new = canonical_name)],
             by = "raw_name", all = TRUE))
  m <- m[is.na(canonical_old) | is.na(canonical_new) | canonical_old != canonical_new]
  if (nrow(m)) cmp[[1]] <- m[, .(map = "name", key = raw_name, frozen = canonical_old, rebuilt = canonical_new)]
  log_msg("  name map vs data/: ", nrow(m), " of ", nrow(old), " raw names differ")
}
if (file.exists("data/rxnorm_class_map.csv")) {
  old <- fread("data/rxnorm_class_map.csv")
  m <- suppressWarnings(merge(old[, .(canonical_name, class_old = assigned_class)], class_map[, .(canonical_name, class_new = assigned_class)],
             by = "canonical_name", all = TRUE))
  m <- m[is.na(class_old) | is.na(class_new) | class_old != class_new]
  if (nrow(m)) cmp[[2]] <- m[, .(map = "class", key = canonical_name, frozen = class_old, rebuilt = class_new)]
  log_msg("  class map vs data/: ", nrow(m), " of ", nrow(old), " entries differ")
}
cmp <- suppressWarnings(rbindlist(cmp))
if (!nrow(cmp)) cmp <- data.table(map = character(0), key = character(0), frozen = character(0), rebuilt = character(0))
fwrite(cmp, out_file("rxnorm_lookup_comparison.csv"))
log_msg("  wrote rxnorm/rxnorm_lookup_comparison.csv (", nrow(cmp), " differences)")
log_msg("=== 08 COMPLETE | ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min ===")
