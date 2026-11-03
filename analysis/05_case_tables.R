# ==============================================================================
# 05_case_tables.R - case-level description of the SSRI lethal-violence reports.
#
#   - indications and report source added to the case lists
#   - top-20 indications per outcome (Fig 2)
#   - Table 1 for the Homicide and H-S case lists (Supplementary "GLOBALS")
#   - integrated per-case listings (Supplementary "H CASE LIST", "H-S CASE LIST")
#   - co-reported PTs by symptom domain, co-drug classes, circumstances
#
# Reads the case lists and signal tables written by steps 02 and 03.
# ==============================================================================
source("R/load_all.R")
need_faers()
t0 <- Sys.time()
log_msg("=== 05 CASE TABLES ===")

# ── indications + report source per case ──────────────────────────────────────
read_cases <- function(file, drug_group_ = NULL) {
  d <- read_res(file, required = FALSE)
  if (is.null(d)) return(NULL)
  d[, primaryid := bit64::as.integer64(primaryid)]
  if (!is.null(drug_group_) && "drug_group" %in% names(d)) d <- d[drug_group == drug_group_]
  d
}
get_indi <- function(obj) {
  it <- as.data.table(faers_get(obj, "indi"))
  u  <- unique(it$indi_pt)
  it[, indi_pt := setNames(title_case(tolower(u)), u)[indi_pt]]
  it <- it[!tolower(indi_pt) %in% INDI_PLACEHOLDERS & !is.na(indi_pt) & indi_pt != ""]
  it[, primaryid := bit64::as.integer64(primaryid)]
}
get_rpsr <- function(obj) {
  rp <- as.data.table(faers_get(obj, "rpsr"))
  rp <- rp[, .(rpsr_cod = paste(sort(unique(toupper(rpsr_cod))), collapse = "; ")), by = primaryid]
  rp[, primaryid := bit64::as.integer64(primaryid)]
}

ssri_obj  <- readRDS(require_res("ssri_group.rds"))
ssri_indi <- get_indi(ssri_obj)
ssri_rpsr <- get_rpsr(ssri_obj)
rm(ssri_obj); gc2()

with_indications <- list(
  c("cases_homicide_suicide.csv",   "cases_homicide_suicide_with_indications.csv"),
  c("cases_homicide.csv",           "cases_homicide_with_indications.csv"),
  c("cases_completed_suicide.csv",  "cases_completed_suicide_with_indications.csv"),
  c("cases_serotonin_syndrome.csv", "cases_serotonin_syndrome_with_indications.csv"))
for (f in with_indications) {
  d <- read_cases(f[1], "ssri")
  if (is.null(d) || !nrow(d)) next
  pc  <- ssri_indi[primaryid %in% d$primaryid, .(indications = paste(sort(unique(indi_pt)), collapse = "; "),
                                                 n_indications = uniqueN(indi_pt)), by = primaryid]
  out <- merge(d, pc, by = "primaryid", all.x = TRUE)
  out[is.na(indications),   indications := "(none recorded)"]
  out[is.na(n_indications), n_indications := 0L]
  out <- merge(out, ssri_rpsr, by = "primaryid", all.x = TRUE)
  write_res(out, f[2])
}

# ── top-20 indications per outcome (Fig 2) ────────────────────────────────────
top_indi <- function(ids, adr) {
  ids <- unique(ids); n_adr <- length(ids)
  sub <- ssri_indi[primaryid %in% ids]
  if (!nrow(sub)) return(data.table(DRUG_GROUP = "SSRI", ADR = adr, adr_n_cases = n_adr, n_with_indication = 0L,
                                    rank = 1L, indication = "(none recorded)", n = 0L, pct = 0))
  t <- sub[, .(n = uniqueN(primaryid)), by = indi_pt][order(-n)][seq_len(min(20L, .N))]
  t[, .(DRUG_GROUP = "SSRI", ADR = adr, adr_n_cases = n_adr, n_with_indication = uniqueN(sub$primaryid),
        rank = seq_len(.N), indication = indi_pt, n, pct = round(100 * n / n_adr, 1))]
}
ssri_ids <- function(pt) { d <- read_cases(pt_files[[pt]], "ssri"); if (is.null(d)) NULL else d$primaryid }
adr_list <- c("Serotonin Syndrome", "Homicide", "Aggression", "Homicidal Ideation", "Violence-Related Symptom",
              "Completed Suicide", "Suicidal Ideation", "Suicide Attempt", "Intentional Self-Injury")
indi_top <- c(
  lapply(adr_list, function(p) top_indi(ssri_ids(p), p)),
  list(top_indi(read_cases("cases_homicide_suicide.csv")$primaryid, "Homicide + Completed Suicide (co-reported)")),
  lapply(c("Homicidality", "Suicidality"), function(g)        # union of the group's PT case lists
    top_indi(do.call(c, lapply(group_defs[[g]], ssri_ids)), paste(g, "Group"))))
write_res(rbindlist(indi_top, fill = TRUE), "adr_indication_top20.csv")
rm(ssri_indi); gc2()

# ── Table 1: Homicide and H-S case lists ──────────────────────────────────────
occp_map <- c(MD = "Health professional", PH = "Health professional", RN = "Health professional",
              HP = "Health professional",
              OT = "Other health professional", CN = "Consumer/non-health professional", LW = "Lawyer/legal")
country_norm <- function(x) {
  x <- toupper(trimws(as.character(x)))
  map <- c(US = "UNITED STATES", USA = "UNITED STATES", "UNITED STATES OF AMERICA" = "UNITED STATES",
           GB = "UNITED KINGDOM", UK = "UNITED KINGDOM", FR = "FRANCE", DE = "GERMANY", IT = "ITALY",
           NL = "NETHERLANDS", CA = "CANADA", ES = "SPAIN", JP = "JAPAN", AU = "AUSTRALIA", SE = "SWEDEN",
           CH = "SWITZERLAND", BR = "BRAZIL", IE = "IRELAND", BE = "BELGIUM", DK = "DENMARK", NO = "NORWAY",
           FI = "FINLAND", AT = "AUSTRIA", PT = "PORTUGAL", NZ = "NEW ZEALAND",
           TR = "TURKEY", AR = "ARGENTINA", CN = "CHINA", IN = "INDIA", KR = "SOUTH KOREA",
           TW = "TAIWAN", MX = "MEXICO", ZA = "SOUTH AFRICA", PL = "POLAND", GR = "GREECE",
           CZ = "CZECHIA", HU = "HUNGARY", IL = "ISRAEL", SG = "SINGAPORE", RU = "RUSSIA")
  x[x %in% c("COUNTRY NOT SPECIFIED", "UNKNOWN", "NOT SPECIFIED")] <- ""
  out <- unname(map[x]); ifelse(is.na(out), x, out)
}
pct <- function(n, N) if (N > 0) sprintf("%.1f", 100 * n / N) else ""

build_table1 <- function(d, label) {
  N <- nrow(d); rows <- list()
  add <- function(char, num = "", prc = "") rows[[length(rows) + 1]] <<- data.table(Characteristic = char, n = as.character(num), `%` = as.character(prc))
  tab_rows <- function(tab, order_names = names(tab)) for (nm in order_names) add(paste0("  ", nm), as.integer(tab[nm]), pct(as.integer(tab[nm]), N))

  add(sprintf("%s group", label), N)
  # Read the label and nothing else: an OR on n_codrugs == 0 would put a
  # report naming no drug back into monopharmacy, which is what
  # mono_poly_label() exists to prevent. Any such report shows as its own row rather than being
  # absorbed, so the three still sum to N.
  # Each of the three counted directly, none as the residual: a residual always
  # sums to N and so can never disagree with anything. reconcile.R adds them
  # and compares with N, which is a real check only because of this.
  lbl      <- tolower(as.character(d$polypharmacy))
  mono_n   <- sum(lbl %in% c("mono", "monopharmacy"))
  poly_n   <- sum(lbl %in% c("poly", "polypharmacy"))
  nodrug_n <- sum(lbl %in% "no drug")
  add("  Monopharmacy (0 co-drugs)", mono_n, pct(mono_n, N))
  add("  Polypharmacy (>=1 co-drug)", poly_n, pct(poly_n, N))
  if (nodrug_n > 0L) add("  No drug named", nodrug_n, pct(nodrug_n, N))

  add("SSRI drug")
  for (dr in ssri_levels) { n <- sum(tolower(d$ssri_drug) == dr); add(paste0("  ", tools::toTitleCase(dr)), n, pct(n, N)) }

  pc <- per_case_codrugs(d)
  add("Top 10 co-drugs")
  dtab <- sort(table(unlist(pc)), decreasing = TRUE); dtab <- dtab[names(dtab) != ""]
  if (length(dtab)) for (nm in head(names(dtab), 10L)) add(paste0("  ", tools::toTitleCase(tolower(nm))), as.integer(dtab[nm]), pct(as.integer(dtab[nm]), N))
  else add("  (no co-drugs recorded)", 0, "0.0")

  add("Top 10 co-drug classes")      # not mutually exclusive
  ctab <- sort(table(unlist(per_case_classes(pc, as.character(d$ssri_drug)))), decreasing = TRUE)
  if (length(ctab)) tab_rows(ctab, c(head(setdiff(names(ctab), "Other"), 9L), if ("Other" %in% names(ctab)) "Other"))
  else add("  (no co-drugs recorded)", 0, "0.0")

  add("Substance Use")
  stab <- sort(table(unlist(per_case_substances(pc))), decreasing = TRUE)
  if (length(stab)) tab_rows(stab) else add("  (no substance use recorded)", 0, "0.0")

  add("Sex")
  sx <- toupper(substr(as.character(d$sex), 1, 1))
  add("  Female", sum(sx == "F"), pct(sum(sx == "F"), N)); add("  Male", sum(sx == "M"), pct(sum(sx == "M"), N))
  add("  Unknown", sum(!sx %in% c("F", "M")), pct(sum(!sx %in% c("F", "M")), N))

  add("Age")
  ag <- suppressWarnings(as.numeric(d$age_in_years)); ag <- ag[is.finite(ag)]
  if (length(ag)) {
    add("  Median age, years (IQR)", siq(ag))
    add("  Mean age, years +/- SD", sub("\\((.*)\\)", "+/- \\1", smsd(ag)))
    add("  Age recorded", length(ag), pct(length(ag), N))
    # half-open bands: ages are not whole years, so 44.5 must fall in a band
    bn <- c(`<18` = sum(ag < 18), `18-44` = sum(ag >= 18 & ag < 45), `45-64` = sum(ag >= 45 & ag < 65),
            `65-74` = sum(ag >= 65 & ag < 75), `>=75` = sum(ag >= 75))
    for (b in names(bn)) add(paste0("    ", b, " years"), bn[[b]], pct(bn[[b]], N))
  } else add("  Age recorded", 0, "0.0")

  # strictly by frequency: pinning a fixed set ahead of the order hides
  # countries that outrank them and contradicts the case listings
  add("Top reporting countries")
  ct  <- country_norm(d$reporter_country)
  tab <- sort(table(ct[ct != "" & !is.na(ct)]), decreasing = TRUE)
  shown <- head(names(tab), 8L)
  for (c0 in shown) add(paste0("  ", tools::toTitleCase(tolower(c0))), as.integer(tab[c0]), pct(as.integer(tab[c0]), N))
  rest <- N - sum(as.integer(tab[shown]))
  if (rest > 0) add("  Other/not recorded", rest, pct(rest, N))

  add("Report by (reporter occupation)")
  oc  <- toupper(trimws(as.character(d$occp_cod)))
  grp <- ifelse(oc %in% names(occp_map), occp_map[oc], "Unknown/other")
  for (g in c("Health professional", "Other health professional", "Consumer/non-health professional", "Lawyer/legal", "Unknown/other")) {
    n <- sum(grp == g); if (n > 0) add(paste0("  ", g), n, pct(n, N)) }

  if ("rpsr_cod" %in% names(d)) {
    add("Report source (rpsr_cod; a report may carry several)")
    rs <- toupper(as.character(d$rpsr_cod))
    for (code in names(rpsr_labels)) { n <- sum(grepl(paste0("\\b", code, "\\b"), rs)); if (n > 0) add(paste0("  ", rpsr_labels[[code]]), n, pct(n, N)) }
    # a report may carry more than one source code and most carry none, so the
    # rows above do not sum to N without this one
    n_none <- sum(is.na(rs) | !nzchar(trimws(rs)))
    add("  Not recorded", n_none, pct(n_none, N))
  }

  add("Top 10 indications for use")
  ind <- as.character(d$indications)
  unknown_n <- sum(ind == "(none recorded)" | ind == "" | is.na(ind))
  add("  Indication unknown / not recorded", unknown_n, pct(unknown_n, N))
  it <- sort(table(unlist(strsplit(ind[ind != "(none recorded)" & ind != "" & !is.na(ind)], "; ", fixed = TRUE))), decreasing = TRUE)
  for (nm in head(names(it), 10L)) add(paste0("  ", nm), as.integer(it[nm]), pct(as.integer(it[nm]), N))

  add("Indications by category")    # same domain taxonomy as the PTs; not mutually exclusive
  icat <- lapply(strsplit(ind, "; ", fixed = TRUE), function(v) {
    v <- v[!is.na(v) & v != "" & v != "(none recorded)"]
    if (!length(v)) character(0) else unique(vapply(v, pt_cat_broad, character(1))) })
  itab <- sort(table(unlist(icat)), decreasing = TRUE)
  if (length(itab)) tab_rows(itab, c(setdiff(names(itab), "Other"), if ("Other" %in% names(itab)) "Other"))
  else add("  (no indications recorded)", 0, "0.0")
  rbindlist(rows)
}
hom <- read_res("cases_homicide_with_indications.csv")
hs  <- read_res("cases_homicide_suicide_with_indications.csv")
write_res(build_table1(hs,  "SSRI H-S (unrestricted)"),      "table1_hs_ssri.csv")
write_res(build_table1(hom, "SSRI Homicide (unrestricted)"), "table1_homicide_ssri.csv")

# ── integrated per-case listings (Supplementary case lists) ───────────────────
cap1 <- function(s) ifelse(is.na(s) | s == "", s, paste0(toupper(substr(s, 1, 1)), tolower(substr(s, 2, nchar(s)))))
# the same string in sentence case, with Roman numerals kept upper: cap1() alone
# turns "Bipolar I Disorder" back into "Bipolar i disorder".
cap1_term <- function(s) roman_up(cap1(s))
split_apply <- function(x, f) vapply(strsplit(as.character(x), "; ", fixed = TRUE), function(v) {
  v <- v[v != ""]; if (!length(v)) "" else paste(vapply(v, f, character(1)), collapse = "; ") }, character(1))
sc_codrug <- function(s) {
  role <- regmatches(s, regexpr("\\[[^]]*\\]\\s*$", s))
  paste0(cap1(trimws(sub("\\[[^]]*\\]\\s*$", "", s))), if (length(role) && nzchar(role)) paste0(" ", toupper(trimws(role))) else "")
}
expand_codes <- function(x, map) vapply(strsplit(toupper(as.character(x)), "; ", fixed = TRUE), function(v) {
  v <- v[v != ""]; if (!length(v)) "" else paste(ifelse(v %in% names(map), map[v], v), collapse = "; ") }, character(1))
# partial FAERS dates (YYYY or YYYYMM) are left missing rather than guessed

make_case_table <- function(in_file, out_file) {
  d <- read_res(in_file, required = FALSE)
  if (is.null(d) || !nrow(d)) return(invisible())
  # ssri_drug_roles comes from step 02/03. Run 05 alone against listings written
  # before that column existed and data.table recycles a zero-length DRUG to NA
  # with a warning, shipping the sheet with the column empty. Stop instead.
  if (!"ssri_drug_roles" %in% names(d))
    stop(in_file, " has no ssri_drug_roles column: rerun steps 02 and 03 before 05.",
         call. = FALSE)
  tto <- suppressWarnings(as.integer(pdate(d$event_dt) - pdate(d$start_dt)))
  tto[!is.finite(tto) | tto < 0] <- NA
  out <- data.table(
    PRIMARYID = d$primaryid, CASEID = d$caseid, YEAR = d$year, QUARTER = d$quarter,
    SEX = cap1(d$sex), AGE = d$age_in_years, WEIGHT = d$wt,
    COUNTRY = tools::toTitleCase(tolower(country_norm(d$reporter_country))),
    `OCCP COD` = d$occp_cod,
    `REPORT SOURCE` = expand_codes(d$rpsr_cod, rpsr_labels),
    # with the role(s) the index drug was recorded in: a report naming it as
    # primary suspect and again as secondary shows [PS, SS] here, which is the
    # only place that fact appears, since the other-SSRI column holds only a
    # different drug. This sheet carries no role column of its own.
    DRUG = split_apply(d$ssri_drug_roles, sc_codrug),
    # all_pts holds the reaction terms as the reporter recorded them, which is
    # not the same as the preferred terms the signals were counted on: a term
    # recorded at lowest level term appears here under that name and in the
    # analysis under its PT ("Schizophrenia, paranoid type" under
    # "Schizophrenia"). Heading the column PTs invited the reader to match it
    # against SI "SYMPTOM DOMAINS", where two terms then appear to be missing.
    `REACTION TERMS (AS RECORDED)` = split_apply(d$all_pts, cap1_term),
    `ALL OUTCOMES` = expand_codes(d$all_outcomes, outc_labels),
    # A second SSRI counts as SSRI exposure, not as co-medication, so it is
    # excluded from codrug_names - but it is still counted in n_codrugs and it
    # is still what makes the report polypharmacy. The other SSRI therefore
    # has a column of its own and the count is always shown. Showing
    # codrug_names alone, falling back to "(None Recorded)" when it is empty,
    # would read 23 of the 334 homicide reports and 4 of the 47 H-S reports as
    # monopharmacy when they are not: the homicide split would come to 161/173
    # against the true 138/196.
    `OTHER SSRI(S) [ROLE]` = ifelse(is.na(d$ssri_other_roles) | d$ssri_other_roles == "",
      "(None)", split_apply(d$ssri_other_roles, sc_codrug)),
    `CODRUG(S) [ROLE] (NO. OF OTHER DRUGS ON REPORT)` = paste0(
      ifelse(is.na(d$codrug_names) | d$codrug_names == "", "(None Recorded)",
             split_apply(d$codrug_names, sc_codrug)),
      " (", ifelse(is.na(d$n_codrugs), 0, d$n_codrugs), ")"),
    # The co-drug names above are as the reporter wrote them, but classes were
    # assigned from RxNorm-resolved ingredients, so Xanax counts under
    # alprazolam and Ambien under zolpidem. Without this column the class
    # counts in SI "GLOBALS" (anti-anxiety 72, Z-drugs 40, alcohol 32) cannot
    # be reached from the listing: 89 of the 177 distinct recorded names are
    # brand or salt forms that appear in no class list. Counting the rows that
    # name a class reproduces each figure directly.
    `CO-DRUG CLASS(ES)` = vapply(
      per_case_classes(per_case_codrugs(d), as.character(d$ssri_drug)),
      function(v) if (!length(v)) "(None)" else paste(sort(v), collapse = "; "),
      character(1)),
    `TIME TO ONSET (days)` = tto,
    INDICATION = ifelse(is.na(d$indications) | d$indications == "(none recorded)", "(None Recorded)", split_apply(d$indications, cap1_term)))
  # ordered on the drug name, not the DRUG column: that column carries the
  # role in brackets, and sorting on it would file "Citalopram [C, PS]" before
  # "Citalopram [PS]" and move most rows of the sheet for no reason.
  write_res(out[order(tools::toTitleCase(as.character(d$ssri_drug)), PRIMARYID)], out_file)
}
make_case_table("cases_homicide_with_indications.csv",           "SI_DATA_H_CASE_LIST.csv")
make_case_table("cases_homicide_suicide_with_indications.csv",    "SI_DATA_HS_CASE_LIST.csv")
make_case_table("cases_completed_suicide_with_indications.csv",   "caselist_completed_suicide_integrated.csv")
make_case_table("cases_serotonin_syndrome_with_indications.csv",  "caselist_serotonin_syndrome_integrated.csv")

# ── co-reported PTs by symptom domain (mono / poly) ───────────────────────────
# TOTAL rows count DISTINCT cases with >= 1 PT in a domain (not the sum of PT rows).
co_pt_matrix <- function(d, exclude_pts, out_file) {
  d  <- copy(d)[, rid := .I]
  # the label itself, not "anything that is not poly": negating it would put a
  # report naming no drug back under mono, which is not the rule the rest of
  # this file uses, and the two would contradict each other in one SI.
  d[, strat := tolower(as.character(polypharmacy))]
  d <- d[strat %chin% c("mono", "poly")]
  # pts_to_pt first: all_pts holds the terms as recorded, so a term coded at
  # Lowest Level Term would otherwise escape the index-term exclusion below and
  # be counted as a co-reported symptom of the outcome it IS
  ex <- d[, .(pt = toupper(trimws(unlist(strsplit(pts_to_pt(all_pts), "; ", fixed = TRUE))))), by = .(rid, strat)]
  ex <- ex[!is.na(pt) & pt != "" & !pt %in% toupper(exclude_pts)]
  ex[, `:=`(cat = vapply(pt, pt_cat, character(1)), bcat = vapply(pt, pt_cat_broad, character(1)))]
  mat <- dcast(ex[, .N, by = .(pt, strat)], pt ~ strat, value.var = "N", fill = 0)
  for (col in c("mono", "poly")) if (!col %in% names(mat)) mat[, (col) := 0L]
  mat[, total := mono + poly]
  mat <- mat[, .(PT = title_case(tolower(pt)), category = vapply(pt, pt_cat, character(1)),
                 broad_category = vapply(pt, pt_cat_broad, character(1)), mono, poly, total)][order(-total)]
  union_row <- function(ids, label, category, broad) {
    if (!nrow(ids)) return(NULL)
    w <- ids[, .N, by = strat]
    mono <- w[strat == "mono", sum(N)]; poly <- w[strat == "poly", sum(N)]
    data.table(PT = label, category = category, broad_category = broad, mono = mono, poly = poly, total = mono + poly)
  }
  totals <- rbindlist(c(
    lapply(c("Manic/hypomanic", "Psychosis"), function(g)
      union_row(unique(ex[cat == g, .(rid, strat)]), paste0("TOTAL (>=1 PT, ", g, ")"), g, NA_character_)),
    lapply(setdiff(broad_cat_levels, "Other"), function(g)
      union_row(unique(ex[bcat == g, .(rid, strat)]), paste0("TOTAL BROAD (>=1 PT, ", g, ")"), NA_character_, g))))
  out <- rbind(totals, mat)[, n_cases := nrow(d)]
  write_res(out, out_file)
  out
}
copt <- list(Homicide = co_pt_matrix(hom, "Homicide", "copt_matrix_homicide.csv"),
             `H-S`    = co_pt_matrix(hs, c("Homicide", "Completed Suicide"), "copt_matrix_hs.csv"))

# Supplementary "SYMPTOM DOMAINS": every co-reported PT assigned to a domain, with counts.
dom_pts <- rbindlist(lapply(names(copt), function(cl)
  copt[[cl]][!grepl("^TOTAL", PT) & broad_category != "Other", .(outcome = cl, PT, broad_category, n = total)]))
write_res(dcast(dom_pts, broad_category + PT ~ outcome, value.var = "n", fill = 0)[
  order(match(broad_category, broad_cat_levels), PT)], "symptom_domain_pts.csv")

# Domain prevalence per outcome, plus "any psychiatric domain" (distinct cases).
cases_by_outcome <- list(Homicide = hom, `H-S` = hs)
write_res(rbindlist(lapply(names(copt), function(cl) {
  tot <- copt[[cl]][grepl("^TOTAL BROAD", PT)]
  n   <- nrow(cases_by_outcome[[cl]])
  any_n <- sum(vapply(pts_to_pt(cases_by_outcome[[cl]]$all_pts), function(s)
    any(vapply(strsplit(as.character(s), "; ", fixed = TRUE)[[1]], function(p) pt_cat_broad(p) != "Other", logical(1))), logical(1)))
  rbind(tot[, .(outcome = cl, category = broad_category, n = total, n_cases = n, pct = round(100 * total / n, 1))],
        data.table(outcome = cl, category = "SUM (non-exclusive, categories overlap)", n = sum(tot$total), n_cases = n,
                   pct = round(100 * sum(tot$total) / n, 1)),
        data.table(outcome = cl, category = "ANY category (distinct cases, true union)", n = any_n, n_cases = n,
                   pct = round(100 * any_n / n, 1)))
}))[order(outcome, -pct)], "pt_category_summary.csv")

# ── circumstances of violence (descriptive) ───────────────────────────────────
circ_terms <- c("aggression", "physical assault", "gun shot wound", "imprisonment", "anger", "irritability",
                "child abuse", "sexual abuse", "stab wound", "drug interaction", "overdose", "intentional overdose",
                "toxicity to various agents", "drug withdrawal syndrome", "withdrawal syndrome", "drug ineffective",
                "alcohol use", "condition aggravated")
# Some terms contain others as substrings ("overdose" matches "intentional
# overdose"), so the rows are not additive. Each row therefore reports both the
# reports it matches and the reports it matches that no narrower term does, and
# names the terms it subsumes.
circ_narrower <- lapply(circ_terms, function(t) circ_terms[circ_terms != t & grepl(t, circ_terms, fixed = TRUE)])
names(circ_narrower) <- circ_terms

circ <- rbindlist(lapply(names(cases_by_outcome), function(cl) {
  pts <- tolower(pts_to_pt(cases_by_outcome[[cl]]$all_pts)); n <- length(pts)
  rbindlist(lapply(circ_terms, function(t) {
    hit  <- grepl(t, pts, fixed = TRUE)
    nar  <- circ_narrower[[t]]
    only <- hit & !Reduce(`|`, lapply(nar, function(x) grepl(x, pts, fixed = TRUE)), init = rep(FALSE, n))
    data.table(outcome = cl, term = tools::toTitleCase(t),
               n = sum(hit), n_exclusive = sum(only),
               subsumes = paste(tools::toTitleCase(nar), collapse = "; "),
               n_cases = n, pct = round(100 * sum(hit) / n, 1)) }))
}))
write_res(circ[n > 0][order(outcome, -pct)], "circumstances_of_violence_summary.csv")

# ── co-drug classes and individual co-drugs ───────────────────────────────────
write_res(rbindlist(lapply(names(cases_by_outcome), function(cl) {
  d <- cases_by_outcome[[cl]]; n <- nrow(d); pc <- per_case_codrugs(d)
  ctab <- table(unlist(per_case_classes(pc, as.character(d$ssri_drug))))
  stab <- table(unlist(per_case_substances(pc)))
  rbind(data.table(outcome = cl, class = names(ctab), type = "Pharmaceutical", n = as.integer(ctab), n_cases = n),
        data.table(outcome = cl, class = names(stab), type = "Substance Use",  n = as.integer(stab), n_cases = n))
}))[, pct := round(100 * n / n_cases, 1)][order(outcome, type, -pct)], "drug_class_summary.csv")

# Supplementary "DRUG CLASSES": every co-reported ingredient with its class and how the
# class was assigned (Other SSRIs rule, RxNorm-resolved class, keyword list, substance list).
class_source <- function(nm, idx) {
  cn <- tolower(nm)
  idxs <- if (is.na(idx) || idx == "") character(0) else tolower(trimws(strsplit(idx, ";", fixed = TRUE)[[1]]))
  if (cn %in% ssri_levels && !(cn %in% idxs)) return("Other SSRI (not the report's index SSRI)")
  if (length(RX_CLASS)) {
    rc <- RX_CLASS[rx_canonical(toupper(trimws(nm)))]
    if (!is.na(rc) && nzchar(rc) && !identical(unname(rc), "Other")) return("RxNorm-resolved class")
  }
  for (cls in names(codrug_classes))
    if (any(vapply(codrug_classes[[cls]], function(k) grepl(k, cn, fixed = TRUE), logical(1)))) return("Keyword list")
  "Unclassified"
}
ingr <- rbindlist(lapply(names(cases_by_outcome), function(cl) {
  d <- cases_by_outcome[[cl]]; pc <- per_case_codrugs(d); idx <- as.character(d$ssri_drug)
  rbindlist(Map(function(v, ix) if (length(v)) data.table(outcome = cl, ingredient = v,
    class = vapply(v, classify_codrug, "", index_drug = ix), assigned_by = vapply(v, class_source, "", idx = ix),
    substance = vapply(v, classify_substance_use, "")), pc, idx))
}))
write_res(dcast(ingr[, .N, by = .(class, assigned_by, substance, ingredient, outcome)],
                class + ingredient + assigned_by + substance ~ outcome, value.var = "N", fill = 0)[order(class, ingredient)],
          "codrug_class_ingredients.csv")

drugs_full <- rbindlist(lapply(names(cases_by_outcome), function(cl) {
  d <- cases_by_outcome[[cl]]
  dtab <- sort(table(unlist(per_case_codrugs(d))), decreasing = TRUE)
  dtab <- dtab[names(dtab) != ""]
  data.table(outcome = cl, drug = tools::toTitleCase(tolower(names(dtab))), n = as.integer(dtab),
             n_cases = nrow(d), pct = round(100 * as.integer(dtab) / nrow(d), 1))
}))
write_res(drugs_full, "drug_summary_full.csv")
# The 15 most co-reported drugs across the outcomes TOGETHER, then split by
# outcome - not a per-outcome top 15. So an outcome's rows here are not its own
# most frequent co-drugs: a drug common in homicide appears under H-S on one
# case, while a drug seen twice in H-S and never in homicide can miss the cut.
# drug_summary_full.csv is the per-outcome ranking, and is what the article and
# the Supplementary Information are built from.
top15 <- drugs_full[, .(tot = sum(n)), by = drug][order(-tot)][seq_len(min(15, .N)), drug]
write_res(drugs_full[drug %in% top15][order(outcome, -pct)], "drug_summary.csv")

log_msg("=== 05 COMPLETE | ", elapsed(t0), " ===")
