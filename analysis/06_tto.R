# ==============================================================================
# 06_tto.R - time to onset (TTO) and Fig 5 (Supplementary "TTO").
#
# TTO = EVENT_DT (DEMO) - earliest START_DT recorded in full (THER), in days. Reports with missing
# or invalid dates, or an event before therapy start, are excluded.
# Kaplan-Meier curves are capped at 1 year (onset within 365 days = event).
#
# Outcomes: Homicide, H-S, Completed Suicide, Serotonin Syndrome (SSRI cases).
# Varenicline cases appear only in the Homicide per-drug view, as a comparison drug.
# In the per-drug view, reports naming two SSRIs are excluded, not split.
#
# Outputs: tto_homicide_vs_hs.csv (outcome comparison), tto_homicide_hs_mono_vs_poly.csv,
#          tto_homicide_hs_perdrug_full.csv, SI_DATA_TTO.csv,
#          figures/fig5_tto_km.png (panels a and b)
# ==============================================================================
source("R/load_all.R")
suppressPackageStartupMessages({ library(ggplot2); library(scales); library(patchwork) })
t0 <- Sys.time()
log_msg("=== 06 TIME TO ONSET ===")

# partial FAERS dates (YYYY or YYYYMM) are left missing rather than guessed

tto_cases <- function(file, drug_group_) {
  d <- read_res(file, required = FALSE)
  if (is.null(d)) return(NULL)
  if ("drug_group" %in% names(d)) d <- d[drug_group == drug_group_] else if (drug_group_ != "ssri") return(NULL)
  if (!nrow(d)) return(NULL)
  tto <- suppressWarnings(as.integer(pdate(d$event_dt) - pdate(d$start_dt)))
  tto[!is.finite(tto) | tto < 0] <- NA
  out <- data.table(primaryid = d$primaryid, tto_days = tto,
                    drug   = if (drug_group_ == "ssri") tools::toTitleCase(as.character(d$ssri_drug)) else "Varenicline",
                    # the label itself, not "anything that is not poly": a
                    # report naming no drug is neither side of this comparison
                    # and maps to NA here rather than silently joining Mono
                    group  = c(poly = "Poly", mono = "Mono")[tolower(as.character(d$polypharmacy))],
                    source = if (drug_group_ == "ssri") "SSRI" else "Varenicline")
  n_nogroup <- sum(!is.na(out$tto_days) & is.na(out$group))
  if (n_nogroup)
    log_msg("  ", file, ": ", n_nogroup, " report(s) name no drug, left out of the mono/poly comparison")
  out[!is.na(tto_days) & !is.na(group)]
}
outcome_files <- c(Homicide = "cases_homicide.csv", `H-S` = "cases_homicide_suicide.csv",
                  `Completed Suicide` = "cases_completed_suicide.csv",
                  `Serotonin Syndrome` = "cases_serotonin_syndrome.csv")
ssri_tto <- lapply(outcome_files, tto_cases, drug_group_ = "ssri")
VAR_COMPARATOR_OUTCOMES <- "Homicide"     # varenicline shown as comparison drug for these outcomes only
perdrug_tto <- lapply(names(outcome_files), function(cl) {
  d <- rbindlist(list(ssri_tto[[cl]], if (cl %in% VAR_COMPARATOR_OUTCOMES) tto_cases(outcome_files[[cl]], "varenicline")))
  d[!grepl("; ", drug, fixed = TRUE)]     # exclude reports naming two SSRIs
})
names(perdrug_tto) <- names(outcome_files)

# ── summary tables ────────────────────────────────────────────────────────────
km_prep  <- function(d) d[, `:=`(tto_cap = pmin(tto_days, 365L), event_cap = as.integer(tto_days <= 365))]
logrank  <- function(d, var) {
  if (uniqueN(d[[var]]) < 2) return(NA_real_)
  lr <- tryCatch(survdiff(as.formula(paste("Surv(tto_cap, event_cap) ~", var)), data = km_prep(copy(d))), error = function(e) NULL)
  if (is.null(lr)) NA_real_ else signif(1 - pchisq(lr$chisq, length(lr$n) - 1), 4)
}
# Kept as the number they are. as.integer() here truncated towards zero, so a
# median of 19.5 days over an even number of cases was published as 19 and a
# quartile of 307.75 as 307 - always downwards, never rounded. Days are whole
# numbers but a median between two of them need not be, and a reader
# recomputing from the case list gets the halves.
describe <- function(d, iqr_min = 1L) data.table(
  n = nrow(d), median_days = as.numeric(median(d$tto_days)),
  q1_days = if (nrow(d) >= iqr_min) as.numeric(quantile(d$tto_days, 0.25)) else NA_real_,
  q3_days = if (nrow(d) >= iqr_min) as.numeric(quantile(d$tto_days, 0.75)) else NA_real_,
  pct_within_30d = round(100 * mean(d$tto_days <= 30), 1),
  pct_within_60d = round(100 * mean(d$tto_days <= 60), 1),
  pct_within_1yr = round(100 * mean(d$tto_days <= 365), 1))

all_outcomes <- rbindlist(lapply(names(ssri_tto), function(cl) copy(ssri_tto[[cl]])[, outcome := cl]))
p_outcome <- logrank(all_outcomes, "outcome")
# A report co-reporting homicide and completed suicide is in the homicide
# cohort and in the H-S cohort, so stacking the four cohorts holds it twice.
# The strata comparison above wants it on both sides, which is what it
# compares; the Combined row is one row per report, so it is taken once here.
# Counted twice it inflates n and drags the median and the percentages toward
# whichever outcomes overlap.
combined <- unique(all_outcomes, by = "primaryid")
summ0 <- rbind(
  rbindlist(lapply(names(ssri_tto), function(cl) cbind(data.table(outcome = cl), describe(ssri_tto[[cl]]),
                                                        logrank_p_outcome_comparison = p_outcome))),
  cbind(data.table(outcome = "Combined"), describe(combined), logrank_p_outcome_comparison = NA_real_))
write_res(summ0, "tto_homicide_vs_hs.csv")

summ1 <- rbindlist(lapply(names(ssri_tto), function(cl) {
  d <- ssri_tto[[cl]]; p <- if (nrow(d) > 4) logrank(d, "group") else NA_real_
  rbindlist(lapply(c("Mono", "Poly"), function(g) if (nrow(d[group == g]))
    cbind(data.table(outcome = cl, group = g), describe(d[group == g]), logrank_p_mono_vs_poly = p)))
}))
write_res(summ1, "tto_homicide_hs_mono_vs_poly.csv")

summ2 <- rbindlist(lapply(names(perdrug_tto), function(cl) {
  d <- perdrug_tto[[cl]]
  rbindlist(lapply(unique(d$drug), function(dr) { s <- d[drug == dr]
    cbind(data.table(outcome = cl, drug = dr), describe(s, iqr_min = 4L)[, `:=`(n_mono = sum(s$group == "Mono"), n_poly = sum(s$group == "Poly"))]) }))
}))[order(outcome, -n)]
write_res(summ2, "tto_homicide_hs_perdrug_full.csv")

# Publication table: "median (Q1-Q3)", each log-rank p shown once per comparison.
#
# Days are reported as whole days. quantile() type 7 over whole days lands on a
# quarter at worst, but a quartile printed as 307.75 claims six-hour precision
# that the source does not carry: a FAERS time to onset is the gap between two
# recorded dates, day-resolution at best and frequently month-only and imputed.
# The exact value is kept, unrounded, in the numeric columns of the tto_*.csv
# files, which are the reproducible record. Everything a reader sees - this
# table, the Supplementary Information sheet built from it, and the figures
# quoted in the article - comes through here, so the three cannot drift apart.
#
# Halves round up rather than to even. round() in R would send 27.5 to 28 and
# 14.5 to 14, which is hard to defend to a reader comparing two rows.
fmt_day <- function(x) ifelse(is.na(x), "", sprintf("%.0f", floor(x + 0.5)))
fmt_iqr <- function(m, q1, q3) ifelse(is.na(m), "", ifelse(is.na(q1) | is.na(q3), fmt_day(m),
                                sprintf("%s (%s - %s)", fmt_day(m), fmt_day(q1), fmt_day(q3))))
fmt_p   <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3)))
first_only <- function(x) { out <- rep("", length(x)); if (length(x)) out[1] <- x[1]; out }
perdrug_p <- sapply(names(perdrug_tto), function(cl) logrank(perdrug_tto[[cl]], "drug"))
s0 <- summ0[order(match(outcome, c(names(outcome_files), "Combined")))]
s1 <- summ1[order(outcome, factor(group, levels = c("Mono", "Poly")))][, p_fmt := first_only(fmt_p(logrank_p_mono_vs_poly)), by = outcome]
s2 <- summ2[order(outcome, median_days)][, p_fmt := first_only(fmt_p(rep(perdrug_p[[outcome[1]]], .N))), by = outcome]
pub <- rbindlist(list(
  s0[, .(section = "(a) Outcome", stratum = outcome, n, median_iqr_days = fmt_iqr(median_days, q1_days, q3_days),
         pct_within_30d, pct_within_1yr, logrank_p = first_only(fmt_p(logrank_p_outcome_comparison)))],
  s1[, .(section = paste0("(b) ", outcome, ": Mono vs Poly"), stratum = group, n, median_iqr_days = fmt_iqr(median_days, q1_days, q3_days),
         pct_within_30d, pct_within_1yr, logrank_p = p_fmt)],
  s2[, .(section = paste0("(c) ", outcome, ": per drug"), stratum = drug, n, median_iqr_days = fmt_iqr(median_days, q1_days, q3_days),
         pct_within_30d, pct_within_1yr, logrank_p = p_fmt)]))
write_res(pub, "SI_DATA_TTO.csv")

# ── Fig 5: Kaplan-Meier, 2 panels ─────────────────────────────────────────────
#  (a) Homicide vs H-S vs Completed Suicide vs Serotonin Syndrome (SSRI cases)
#  (b) Homicide: monopharmacy vs polypharmacy
#  The per-drug comparison (with varenicline) is reported as a table only:
#  tto_homicide_hs_perdrug_full.csv / SI_DATA_TTO.csv section (c).
km_curve <- function(d, var) {
  d <- km_prep(copy(d)); d[, strat := get(var)]
  if (uniqueN(d$strat) < 2) return(NULL)
  fit <- survfit(Surv(tto_cap, event_cap) ~ strat, data = d)
  sf  <- as.data.table(summary(fit, times = seq(0, 365, 1))[c("time", "surv", "lower", "upper", "strata")])
  sf[, strat := sub("strat=", "", strata, fixed = TRUE)]
  # start every curve at 100% so onsets on day 0 show as a drop at day 0
  s0 <- unique(sf[, .(strata, strat)])[, `:=`(time = 0, surv = 1, lower = 1, upper = 1, .ord = 0L)]
  sf <- rbind(s0, sf[, .ord := 1L], fill = TRUE)[order(strat, time, .ord)][, .ord := NULL]
  list(sf = sf, p = logrank(d, var))
}
lr_lab <- function(p) if (!is.na(p) && p < 0.001) "Log-rank p<0.001" else paste0("Log-rank p=", sprintf("%.3f", p))
km_panel <- function(res, cols, ltys, title, ribbon = TRUE, ylab = NULL, key_w = 1.2) {
  res$sf[, strat := factor(strat, levels = intersect(names(cols), unique(strat)))]
  p <- ggplot(res$sf, aes(time, surv, colour = strat, linetype = strat)) +
    geom_step(linewidth = if (ribbon) 1.1 else 1.0)
  if (ribbon) p <- p + geom_ribbon(aes(ymin = lower, ymax = upper, fill = strat), alpha = 0.10, colour = NA) +
    scale_fill_manual(values = cols, name = NULL)
  p + annotate("text", x = 360, y = 0.97, label = lr_lab(res$p), size = 3.2, hjust = 1,
               colour = "grey30", fontface = "italic", family = FONT) +
    scale_colour_manual(values = cols, name = NULL) + scale_linetype_manual(values = ltys, name = NULL) +
    # Bare numbers: the axis is titled in days, so a "d" on every tick repeats
    # it, and "1yr" set days and years side by side on one scale while 0 was
    # left bare - three ways of writing the same axis. 365 is where the data
    # are censored, which the caption says.
    scale_x_continuous(breaks = c(0, 30, 90, 180, 270, 365), limits = c(0, 365)) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
    theme_pub(10, FONT) + theme(legend.key.width = unit(key_w, "cm"),
      axis.title.y = if (is.null(ylab)) element_blank() else element_text(size = 10, face = "italic", angle = 90)) +
    labs(title = title, x = "Days from initiation", y = ylab)
}
res_a <- km_curve(all_outcomes, "outcome")
res_b <- km_curve(ssri_tto$Homicide, "group")

# km_curve returns NULL when a panel has fewer than two strata to compare, which
# a short study period can produce. The tables above are unaffected; only the
# figure needs both panels.
if (is.null(res_a) || is.null(res_b)) {
  log_msg("  Figure 5 not drawn: fewer than two strata in ",
          paste(c("panel (a)", "panel (b)")[c(is.null(res_a), is.null(res_b))], collapse = " and "))
} else {
  pa <- km_panel(res_a, c(Homicide = "black", `H-S` = "grey50", `Completed Suicide` = "#1f77b4", `Serotonin Syndrome` = "#2ca02c"),
                 c(Homicide = "solid", `H-S` = "dashed", `Completed Suicide` = "dotted", `Serotonin Syndrome` = "dotdash"),
                 "(a)", ylab = "Proportion without event") +
    guides(colour = guide_legend(nrow = 2), fill = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))
  pb <- km_panel(res_b, c(Mono = "black", Poly = "grey50"), c(Mono = "solid", Poly = "dashed"), "(b)")
  # no title or caption: the article's figure legend carries them
  fig5 <- (pa | pb) + plot_annotation(theme = theme(text = element_text(family = FONT)))
  save_plot(fig5, "fig5_tto_km", w = 174 / 25.4, h = 4.0)
}
log_msg("=== 06 COMPLETE | ", elapsed(t0), " ===")
