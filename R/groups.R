# ==============================================================================
# R/groups.R - building drug groups and the monopharmacy classification.
# ==============================================================================

# Reports where a drug matching `pattern` (drug name or active ingredient) is
# recorded; ps_only = TRUE keeps only reports where it is the primary suspect.
drug_filter <- function(data, pattern, ps_only = TRUE) {
  faers_filter(data, .fn = function(x) {
    hit <- !is.na(x$drugname) & grepl(pattern, x$drugname, ignore.case = TRUE, perl = TRUE)
    if ("prod_ai" %in% names(x))
      hit <- hit | (!is.na(x$prod_ai) & grepl(pattern, x$prod_ai, ignore.case = TRUE, perl = TRUE))
    if (ps_only && "role_cod" %in% names(x))
      hit <- hit & !is.na(x$role_cod) & toupper(x$role_cod) == "PS"
    x[hit, unique(primaryid)]
  }, .field = "drug")
}

# Remove the reports identified for exclusion by manual clinical review
# (duplicate submissions of the same event, implausible dates, and repeat
# filings of one case by several manufacturers). Applied to drug groups AND
# reference sets.
apply_blacklist <- function(obj, label = "object") {
  if (!file.exists(BLACKLIST_FILE)) stop(BLACKLIST_FILE, " not found.")
  bl  <- bit64::as.integer64(fread(BLACKLIST_FILE, colClasses = "character")$primaryid)
  if (!length(bl)) stop(BLACKLIST_FILE, " is empty: add the manual-review report IDs before running.")
  ids <- faers_get(obj, "demo")$primaryid
  log_msg("Blacklist: ", sfmt(length(bl)), " listed | ", sfmt(sum(ids %in% bl)), " found in ", label)
  faers_keep(obj, primaryid = id_setdiff(ids, bl))
}

# Label each drug row with the SSRI it matches (first match in dictionary order).
tag_ssri <- function(drug) {
  drug[, ssri := NA_character_]
  has_ai <- "prod_ai" %in% names(drug)
  for (nm in names(ssri_drug_pattern)) {
    pat <- ssri_drug_pattern[[nm]]
    if (has_ai) {
      drug[is.na(ssri) & (grepl(pat, drugname, ignore.case = TRUE, perl = TRUE) |
             (!is.na(prod_ai) & grepl(pat, prod_ai, ignore.case = TRUE, perl = TRUE))), ssri := nm]
    } else {
      drug[is.na(ssri) & grepl(pat, drugname, ignore.case = TRUE, perl = TRUE), ssri := nm]
    }
  }
  drug
}

# ── monopharmacy ───────────────────────────────────────────────────────────────
# Monopharmacy, one definition, used by the study group, the comparator group
# and the whole-FAERS reference set alike.
#
#   A report is monopharmacy when it names exactly one distinct drug.
#
# Distinct drugs are counted by canonical name, irrespective of the role each
# entry carries (PS/SS/C/I). Two things make a report name one drug twice:
#
#   a repeated entry   the same drug recorded as primary suspect and again as a
#                      concomitant. Counting drug entries would make this
#                      polypharmacy; counting distinct names does not.
#   a brand and its generic   "PROZAC" on one line and "FLUOXETINE" on another.
#                      The published drug dictionary (data/drug_dictionary.csv,
#                      Supplementary "DRUG") lists every brand, INN variant and
#                      salt for each study drug, so both resolve to fluoxetine.
#
# The same dictionary decides group membership and which SSRI a report is
# attributed to, so counting uses it too.
#
# Two different study drugs remain two drugs: sertraline plus fluvoxamine is
# polypharmacy. Two spellings of one drug do not: Chantix plus varenicline is
# monotherapy.
#
# Manual reclassifications in data/mono_overrides.csv are applied last.

# Normalized drug name: upper case, trailing punctuation and repeated spaces
# removed. The name is NOT unwrapped to its last bracket.
#
# A bracket often holds the ingredient - "ZOLOFT (SERTRALINE HYDROCHLORIDE)" -
# but reporters put anything there: a form (TABLETS), a strength (50
# MILLIGRAM), a manufacturer (WATSON LABORATORIES), an annotation (CON.), a
# clarifier (UNKNOWN). Keying on the bracket merges every drug on such a
# report into one, so a polypharmacy report counts as monopharmacy, and no
# word list keeps up with every strength.
#
# A name therefore keys as it is written, which is what the Methods
# describes: drugs counted by name, with the study and comparator drugs
# resolved through the dictionary by drug_canon(). For a drug outside the
# dictionary a brand and its ingredient on one report - TYLENOL
# (ACETAMINOPHEN) and ACETAMINOPHEN - count as two. That error can
# only move a report OUT of monopharmacy, never into it.
drug_key <- function(x) {
  x <- as.character(x)
  u <- unique(x)
  v <- toupper(trimws(u))
  v <- sub("[[:punct:]]+$", "", v)
  v <- trimws(gsub("[[:space:]]+", " ", v))
  v[match(x, u)]
}

# Canonical drug name: any term in the study or comparator dictionary becomes
# that drug; every other name keeps its normalized form. Computed over the
# distinct names and mapped back, which is the same result element by element.

drug_canon <- function(x) {
  k <- drug_key(x)
  u <- unique(k)
  canon <- u
  done  <- rep(FALSE, length(u))
  pats  <- lapply(c(drug_terms("study"), drug_terms("comparator")), terms_pattern)
  for (nm in names(pats)) {
    hit <- !done & !is.na(u) & grepl(pats[[nm]], u, ignore.case = TRUE, perl = TRUE)
    canon[hit] <- toupper(nm)
    done <- done | hit
  }
  canon[match(k, u)]
}

# Entries that are not drugs: the placeholders reporters put in the drug field
# to say there was no co-medication. Dropped from the raw name, before
# drug_canon() runs, so the same rule reaches n_codrugs, the mono/poly flag,
# the monopharmacy reference set and the printed co-drug lists alike.
# See data/non_drug_terms.csv.
#
# The normalization here matches drug_key() but is kept separate: this one is
# about recognising a placeholder written as the whole name, not about keying
# and drop it. Upper case, trailing punctuation and repeated spaces only, so
# "NONE." and "NO  CONCURRENT  MEDICATION" match and no real drug can.
non_drug_key <- function(x) {
  v <- toupper(trimws(as.character(x)))
  v <- sub("[[:punct:]]+$", "", v)
  trimws(gsub("[[:space:]]+", " ", v))
}
drop_non_drugs <- function(dt, col = "drugname", drop_terms = NON_DRUG_DROP) {
  dt[!non_drug_key(get(col)) %chin% drop_terms]
}

# Distinct drugs named on each report, zero where no drug is recorded.
report_drug_counts <- function(obj) {
  dcols <- intersect(c("primaryid", "drugname", "prod_ai"), names(faers_get(obj, "drug")))
  drug <- faers_get(obj, "drug")[, ..dcols]
  drug <- drug[!is.na(drugname) & nzchar(trimws(drugname))]
  drug <- drop_non_drugs(drug)
  drug[, dkey := drug_canon(drugname)]
  # The group filter matches drugname OR prod_ai, so the counter must too: a
  # report naming one drug by a brand the dictionary does not hold, and again
  # by its ingredient, is one drug, not two.
  if ("prod_ai" %in% dcols) {
    known <- toupper(names(c(drug_terms("study"), drug_terms("comparator"))))
    # rawname keeps the name as written, upper case and trimmed, so two lines
    # share a rawname only when they really do name the same thing.
    drug[, rawname := non_drug_key(drugname)]
    drug[, aikey := drug_canon(prod_ai)]
    # A name is the drug the dictionary finds in it, the same rule that puts a
    # report in a drug group, so a combination written in the name counts as
    # the study drug it contains. prod_ai is used only to resolve a name the
    # dictionary does not hold, and only where it records a single ingredient:
    # FAERS marks a combination there with a backslash, and resolving that to
    # its first match would merge two drugs into one.
    drug[grepl("\\", prod_ai, fixed = TRUE), aikey := NA_character_]
    # One name on one report is one drug, so where any line of a name resolves
    # to a known ingredient, and the lines agree on which, every line of that
    # name takes it. Resolving line by line instead would split a drug written
    # twice with prod_ai filled in on only one of the lines.
    #
    # The grouping is on the name, not on dkey: two names that resolve to the
    # same dictionary drug share a dkey, and grouping on that would hand them
    # one ingredient between them.
    drug[, resolvable := !(dkey %chin% known) & !is.na(aikey) & aikey %chin% known]
    res <- drug[resolvable == TRUE,
                .(newkey = aikey[1L], n_ai = uniqueN(aikey)),
                by = .(primaryid, rawname)]
    # Each line first takes the ingredient recorded against it, so one name
    # written twice against two different ingredients stays two drugs.
    drug[resolvable == TRUE, dkey := aikey]
    # Then a line carrying no ingredient of its own takes the one its
    # same-named siblings agree on, which is the split this is here to prevent.
    if (nrow(res[n_ai == 1L]))
      drug[res[n_ai == 1L], on = .(primaryid, rawname),
           dkey := fifelse(dkey %chin% known, dkey, i.newkey)]
    drug[, c("aikey", "rawname", "resolvable") := NULL]
  }
  n <- drug[, .(n_distinct_drugs = uniqueN(dkey)), by = primaryid]
  out <- merge(data.table(primaryid = faers_get(obj, "demo")$primaryid), n,
               by = "primaryid", all.x = TRUE)
  out[is.na(n_distinct_drugs), n_distinct_drugs := 0L]
  out[, primaryid := bit64::as.integer64(primaryid)]
  out[]
}

# The one place a report is labelled mono or poly, so that the case lists, the
# Table 1 blocks, the TTO strata and Figure 5b cannot drift from the rule the
# reference set is built with (mono_reference_ids). Three states, not two: a
# report naming no drug is neither, and calling it "poly" would be as wrong as
# calling it "mono". "no drug" should be vanishingly rare on a drug group's
# case list - the group filter also matches prod_ai, so it is possible - and if
# it ever appears it must be visible rather than folded into a side.
# A missing drug count returns NA rather than a guess: it can only come from a
# report that fell out of the merge, and defaulting it to "mono" would hide
# that behind the very rule this function replaces. reconcile.R adds the three
# labels and compares the total to N, so an NA shows up as a FAIL.
mono_poly_label <- function(n_codrugs, n_distinct_drugs) {
  fifelse(is.na(n_distinct_drugs), NA_character_,
          fifelse(n_distinct_drugs == 0L, "no drug",
                  fifelse(!is.na(n_codrugs) & n_codrugs > 0, "poly", "mono")))
}

# Co-drugs on each report: every distinct drug beyond the first.
build_codrug_table <- function(group_obj) {
  cd <- report_drug_counts(group_obj)
  cd[, n_codrugs := pmax(0L, n_distinct_drugs - 1L)]
  cd[primaryid %in% mono_overrides$primaryid, n_codrugs := 0L]
  cd[]
}

# The same rule over every report in the database: the denominator for the
# monopharmacy disproportionality analysis.
# Returns both halves of the split, because the caller needs the reports that
# name no drug as well and report_drug_counts() over the whole database is the
# most expensive thing in this file: one pass, not two.
mono_reference_ids <- function(master_obj) {
  cd  <- build_codrug_table(master_obj)
  # n_codrugs is n_distinct_drugs - 1 floored at zero, so it is 0 both for a
  # report naming one drug and for one naming none - a report with no drug row,
  # or whose only entries were placeholders. Monopharmacy is "exactly one
  # distinct drug named on the report", as summary_master.csv states, so the
  # second kind is not in the reference set. The manual overrides are reports
  # whose duplicate brand and generic the dictionary did not merge, so they name
  # at least one drug and are unaffected.
  # The overrides are re-added unconditionally, so the new condition can never
  # drop one whatever its drug count.
  ids <- cd[n_codrugs == 0 & n_distinct_drugs >= 1L, primaryid]
  list(mono_ids = id_union(ids, mono_overrides$primaryid[
         mono_overrides$primaryid %in% faers_get(master_obj, "demo")$primaryid]),
       nodrug_ids = cd[n_distinct_drugs == 0L, primaryid])
}

# ── per-report case lists ─────────────────────────────────────────────────────
# One row per report: demographics, all PTs, outcomes, the SSRI(s), co-drugs with
# roles, earliest therapy start date and mono/poly status.
case_list_parts <- function(obj, ids) {
  i64  <- bit64::as.integer64
  demo <- faers_get(obj, "demo")[primaryid %in% ids]
  demo <- demo[, intersect(c("primaryid", "caseid", "year", "quarter", "sex", "age_in_years", "wt",
                             "reporter_country", "event_dt", "fda_dt", "mfr_sndr", "occp_cod", "i_f_code"),
                           names(demo)), with = FALSE]
  reac <- faers_get(obj, "reac")[primaryid %in% ids,
                                 .(all_pts = paste(sort(unique(pt)), collapse = "; ")), by = primaryid]
  outc <- tryCatch(faers_get(obj, "outc")[primaryid %in% ids,
                     .(all_outcomes = paste(sort(unique(outc_cod)), collapse = "; ")), by = primaryid],
                   error = function(e) data.table(primaryid = i64(NA_integer_), all_outcomes = NA_character_))
  # The earliest start date that was recorded in full. FAERS stores some THER
  # lines as a year ("2001") or a year and month ("200106"). The column is text,
  # so min() compares lexicographically, and a partial date sorts below a fuller
  # one from the same era - "2001" and "200106" both below "20010619" - so a
  # plain min() over a report's lines returned the partial in nearly every case
  # where one was present (99 of the 123 such reports). That date
  # does not parse, and the report then drops out of the time-to-onset analysis
  # even though another of its lines carried a usable date. Prefer the earliest
  # eight-digit date; fall back to the earliest partial one only when a report
  # has nothing else, so the listing still shows what was recorded.
  ther <- tryCatch(faers_get(obj, "ther")[primaryid %in% ids,
                     .(start_dt = {
                        s <- trimws(as.character(start_dt))
                        s <- s[!is.na(s) & nzchar(s)]
                        f <- s[nchar(s) == 8L]
                        if (length(f)) min(f) else if (length(s)) min(s) else NA_character_
                      }), by = primaryid],
                   error = function(e) data.table(primaryid = i64(NA_integer_), start_dt = NA_character_))
  drug <- faers_get(obj, "drug")[primaryid %in% ids]
  demo[, primaryid := i64(primaryid)]; reac[, primaryid := i64(primaryid)]
  outc[, primaryid := i64(primaryid)]; ther[, primaryid := i64(primaryid)]
  drug[, primaryid := i64(primaryid)]
  list(demo = demo, reac = reac, outc = outc, ther = ther, drug = drug)
}

# Case listings are written in a fixed order so a re-run reproduces them exactly.
order_cases <- function(out) {
  setorder(out, year, primaryid)
  out
}

merge_all <- function(parts) Reduce(function(a, b) merge(a, b, by = "primaryid", all.x = TRUE), parts)

build_case_list <- function(obj, ids, codrug_dt) {
  if (!length(ids)) return(data.table())
  p    <- case_list_parts(obj, ids)
  drug <- tag_ssri(p$drug)
  # The report is attributed to the SSRI(s) recorded as primary suspect, which is
  # the rule the groups and the disproportionality tables use. SSRIs carried in
  # any other role are listed separately, so the full drug record stays visible.
  # Names are sorted so that two reports naming the same SSRIs in a different
  # order produce the same string.
  # Each non-primary SSRI carries its own role, the same way a co-drug does
  # below. Without it the column says only that another SSRI was on the report,
  # and role_cod cannot supply the answer: it holds the set of roles on the
  # report, so "C; PS; SS" over two other SSRIs does not say which held which.
  # A different SSRI, with the role it was recorded in, the same way a co-drug
  # carries its role below. The report's own SSRI is excluded even where it is
  # entered a second time in another role: role_cod already records that, and
  # listing it here made homicide rows name the index drug back.
  sw   <- drug[!is.na(ssri)][, {
          is_ps <- !is.na(role_cod) & toupper(role_cod) == "PS"
          ps    <- unique(ssri[is_ps])
          other <- !is_ps & !(ssri %in% ps)
          .(ssri_drug        = paste(sort(ps), collapse = "; "),
            # the index drug with every role it was recorded in, so a report
            # naming it twice still shows both. Roles are joined with a comma
            # inside the brackets because "; " separates one drug from the next.
            ssri_drug_roles  = paste(vapply(sort(ps), function(dn) {
              rl <- unique(toupper(role_cod[ssri == dn]))
              rl <- rl[order(match(rl, c("PS", "SS", "C", "I")), rl)]   # PS first, as FAERS orders them
              paste0(dn, " [", paste(rl, collapse = ", "), "]") }, character(1)), collapse = "; "),
            ssri_other_roles = if (any(other))
              paste(sort(unique(paste0(ssri[other], " [", toupper(role_cod[other]), "]"))), collapse = "; ") else "",
            role_cod         = paste(sort(unique(toupper(role_cod))), collapse = "; "))
          }, by = primaryid]
  cn   <- drop_non_drugs(drug[is.na(ssri) & !is.na(drugname)], "drugname")[
               , .(codrug_names = paste(sort(unique(paste0(drugname, " [", toupper(role_cod), "]"))), collapse = "; ")),
               by = primaryid]
  cd   <- copy(codrug_dt)[, .(primaryid = bit64::as.integer64(primaryid), n_codrugs, n_distinct_drugs)]
  out  <- merge_all(list(p$demo, p$reac, p$outc, sw, cd, cn, p$ther))
  out[, polypharmacy := mono_poly_label(n_codrugs, n_distinct_drugs)]
  out <- order_cases(out)
  out[, intersect(c("primaryid", "caseid", "year", "quarter", "sex", "age_in_years", "wt",
                    "reporter_country", "event_dt", "start_dt", "fda_dt", "occp_cod", "i_f_code",
                    "ssri_drug", "ssri_drug_roles", "ssri_other_roles", "role_cod", "all_pts", "all_outcomes",
                    "n_codrugs", "codrug_names", "polypharmacy"),
                  names(out)), with = FALSE]
}

build_ctrl_case_list <- function(obj, ids, ctrl_codrug_dt, ctrl_name) {
  if (!length(ids)) return(data.table())
  p  <- case_list_parts(obj, ids)
  cn <- drop_non_drugs(p$drug[!is.na(drugname)], "drugname")[
               , .(all_drugs = paste(sort(unique(paste0(drugname, " [", toupper(role_cod), "]"))), collapse = "; ")),
               by = primaryid]
  out <- merge_all(list(p$demo[, setdiff(names(p$demo), "mfr_sndr"), with = FALSE], p$reac, p$outc,
                        ctrl_codrug_dt[, .(primaryid, n_codrugs, n_distinct_drugs)], cn, p$ther))
  out[, comparator_drug := ctrl_name]
  out[, polypharmacy := mono_poly_label(n_codrugs, n_distinct_drugs)]
  out <- order_cases(out)
  # n_distinct_drugs is carried only to label the row, and is dropped again
  # here: without this the comparator rows would add a column that the SSRI
  # rows do not have, and rbindlist(fill = TRUE) would put it into every
  # deposited case list with NA down the SSRI side.
  out[, setdiff(names(out), "n_distinct_drugs"), with = FALSE]
}

# Case lists for every target PT, SSRI + comparator rows in one file (`drug_group` column).
write_case_lists <- function(files, drug_groups) {
  for (pt_term in names(files)) {
    parts <- list()
    for (co in drug_groups) {
      ids <- get_ids_for_pt(co$signals, pt_term)
      if (!length(ids)) next
      cl <- if (co$name == "ssri") build_case_list(co$obj, ids, co$codrug)
            else build_ctrl_case_list(co$obj, ids, co$codrug, co$name)
      cl[, drug_group := co$name]
      parts[[co$name]] <- cl
    }
    if (!length(parts)) { log_msg("  ", pt_term, ": no cases"); next }
    fwrite(rbindlist(parts, fill = TRUE, use.names = TRUE), res_path(files[[pt_term]]))
    log_msg("  ", pt_term, ": ", paste(names(parts), vapply(parts, nrow, 1L), sep = "=", collapse = " "))
  }
}
