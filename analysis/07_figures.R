# ==============================================================================
# 07_figures.R - manuscript Figures 2, 3 and 4 (Fig 5 is drawn by 06_tto.R;
# Fig 1, the study flow diagram, by tools/make_fig1.py, which carries the counts
# from summary_master.csv and SI_DATA_COUNTS.csv as literals rather than reading
# them - so they need checking by hand after a re-run).
#
#   Fig 2  indications for SSRI reports, by outcome (pie charts)
#   Fig 3  ROR heatmaps: (a) unrestricted, (b) monopharmacy, (c) normalized
#          difference mono vs full with z-test
#   Fig 4  externalizing vs internalizing normalized difference by drug
# ==============================================================================
source("R/load_all.R")
suppressPackageStartupMessages({ library(ggplot2); library(scales); library(patchwork) })
t0 <- Sys.time()
log_msg("=== 07 FIGURES ===")
SSRI6 <- tools::toTitleCase(names(ssri_drugs))
HS    <- "Homicide + Completed Suicide"

# ── sizing for the journal page ───────────────────────────────────────────────
# Figures are drawn at the width they will be printed at, 174 mm, and the type
# is set so it lands inside the journal's 8-12 pt band at that width. Drawing
# large and letting production shrink it does not work: type scales with the
# canvas, so a 24-inch Figure 3 would put its confidence intervals at about
# 2 pt on the page. Eight columns of three lines each is what fits; that is why
# Figure 3 is three panels rather than one.
#
# ggplot's geom_text size is in millimetres, hence the conversion. Everything
# else - theme sizes, element_text - is already in points. Defined here, before
# the first figure, because Figure 2 is saved long before Figure 3 is built.
FIG_W   <- 174 / 25.4          # inches - the journal's full page width
FIG_H   <- 9.0                 # 229 mm, just inside its 234 mm limit
PT      <- 1 / 2.845           # points -> ggplot geom_text size
PT_AXIS <- 8.5

# ── Fig 2: indications ────────────────────────────────────────────────────────
ind <- read_res("adr_indication_top20.csv")[DRUG_GROUP == "SSRI"]
adr_show <- c("Serotonin Syndrome" = "Serotonin Syndrome", "Homicide" = "Homicide", "Completed Suicide" = "Completed Suicide",
              "Homicidality Group" = "Pooled Homicidality", "Suicidality Group" = "Pooled Suicidality",
              "Homicide + Completed Suicide (co-reported)" = "H-S composite")
ind <- ind[ADR %in% names(adr_show)]
adr_ord <- names(adr_show)[names(adr_show) %in% ind$ADR]
first   <- ind[match(adr_ord, ADR)]
lab_map <- setNames(sprintf("%s\n(n = %s / %s)", adr_show[adr_ord], sfmt(first$n_with_indication), sfmt(first$adr_n_cases)), adr_ord)
# The grey wedge is what is left of the panel after the top seven are drawn.
# It was labelled "Other / unreported", which names it as a class of cases, and
# it is not one: a report can record several indications, so the seven overlap
# and the wedge is the remainder after that double counting, not a count of
# anything. For homicide it came to 106 while 119 cases record no indication at
# all. The caption says the same, and the legend agrees with it.
RESIDUAL <- "Balance to 100%"
n_slices <- 7L   # top 7 indications per pie; the rest folds into RESIDUAL
# pct is published at 1 dp in adr_indication_top20.csv, so a label formatted
# from it rounds twice: Anxiety in the H-S panel is 4/47 = 8.5106%, stored as
# 8.5, and sprintf("%.0f") takes 8.5 to 8 rather than 9. The share is recomputed
# from n and the panel denominator so the label and the slice both come off the
# unrounded value.
pie <- rbindlist(lapply(adr_ord, function(a) {
  d <- ind[ADR == a][order(-pct)]; top <- head(d, n_slices); den <- d$adr_n_cases[1]
  rbind(top[, .(ADR = a, indication, n)],
        data.table(ADR = a, indication = RESIDUAL, n = max(den - sum(top$n), 0)))[
          , pct := 100 * n / den][]
}))
pie[, ADR := factor(lab_map[ADR], levels = lab_map[adr_ord])]
slice_ord <- c(head(ind[, .(tn = sum(n)), by = indication][order(-tn)]$indication, n_slices), RESIDUAL)
pie[, indication := factor(indication, levels = rev(unique(c(slice_ord, as.character(pie$indication)))))]
real <- setdiff(levels(pie$indication), RESIDUAL)
slice_cols <- c(setNames(hue_pal()(length(real)), real), setNames("grey80", RESIDUAL))
# Type is set for the printed width, as everywhere else here. This figure used
# to be drawn at 15 in and shrunk, which is why its type was oversized; drawn
# at 174 mm that same type overran the panel titles and ran the 7% and 5%
# labels of the homicide panel into each other.
#
# The key moves under the panels rather than beside them. Nineteen indications
# in a single column need 80 mm of height, which clipped. Three columns across
# the page hold them in seven rows: four columns leave 43 mm a column, and
# "Post-Traumatic Stress Disorder" needs that much for the words alone, before
# the swatch - it ran off the right edge. The key's own title sits above it
# rather than to its left, where it was clipped in turn.
fig2 <- ggplot(pie, aes(x = "", y = pct, fill = indication)) +
  geom_col(colour = "white", linewidth = 0.5, width = 1) +
  # A slice under 6% is too narrow to label without the label meeting its
  # neighbour's: the 7% and 5% indications of the homicide panel sit side by
  # side. The key names every slice, and the Supplementary Information gives
  # all of them with their counts.
  geom_text(aes(label = ifelse(pct >= 6, sprintf("%.0f%%", pct), "")), position = position_stack(vjust = 0.5),
            size = 8 * PT, family = FONT, fontface = "bold") +
  coord_polar(theta = "y") + facet_wrap(~ADR, nrow = 2) +
  scale_fill_manual(values = slice_cols, name = "Indication") +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_legend(ncol = 3, byrow = TRUE, title.position = "top")) +
  theme_void(base_size = 9) +
  theme(text = element_text(family = FONT),
        # the panel titles sit on the top edge without this
        strip.text = element_text(face = "bold", size = PT_AXIS, margin = margin(b = 5)),
        plot.margin = margin(t = 14, r = 6, b = 6, l = 6),
        legend.position = "bottom", legend.title = element_text(size = PT_AXIS, face = "bold"),
        legend.key.size = unit(0.38, "cm"),
        legend.text = element_text(size = 8))
# Drawn at the printed width so the type lands where the journal wants it;
save_plot(fig2, "fig2_indications", w = FIG_W, h = 5.9)

# ── Fig 3: ROR heatmaps, full vs monopharmacy ─────────────────────────────────
sig_cols <- c("drug", "pt", "analysis", "ror", "a", "ror_ci_low", "ror_ci_high", "full_consensus")
take <- function(file, analysis_lab, drug_lab = NULL, pt_lab = NULL, pt_from = NULL) {
  d <- read_res(file, required = FALSE); if (is.null(d)) return(NULL)
  d <- copy(d)[, analysis := analysis_lab]
  if (!is.null(drug_lab)) d[, drug := drug_lab]
  if (!is.null(pt_from))  d[, pt := paste0(get(pt_from), " Group")]
  if (!is.null(pt_lab))   d[, pt := pt_lab]
  d[, intersect(sig_cols, names(d)), with = FALSE]
}
full_mono <- c(full = "Unrestricted", mono = "Monopharmacy")
hm <- rbindlist(fill = TRUE, unlist(recursive = FALSE, lapply(names(full_mono), function(k) {
  co <- full_mono[[k]]; m <- if (k == "mono") "mono_" else ""
  hs_pooled <- if (k == "full") "homicide_suicide_composite_phv.csv" else "mono_homicide_suicide_composite_phv.csv"
  list(take(paste0("ssri_", m, "perdrug_pt_signals.csv"), co),
       take(paste0("ssri_", m, "targeted_signals.csv"), co, drug_lab = "SSRI (POOLED)"),
       take(paste0("varenicline_", m, "targeted_signals.csv"), co, drug_lab = "Varenicline"),
       take(paste0("group_pooled_", k, ".csv"), co, pt_from = "group"),
       take(paste0("ssri_", m, "perdrug_hs_composite.csv"), co, pt_lab = HS),
       take(hs_pooled, co, drug_lab = "SSRI (POOLED)", pt_lab = HS),
       take(paste0("varenicline_", m, "hs_composite_phv.csv"), co, drug_lab = "Varenicline", pt_lab = HS))
})))
hm[drug %in% names(ssri_drugs), drug := tools::toTitleCase(drug)]

drug_ord <- c("SSRI (POOLED)", SSRI6, "Varenicline")
pt_ord   <- c("Serotonin Syndrome", "Homicide", "Aggression", "Homicidal Ideation", "Violence-Related Symptom",
              "Completed Suicide", "Suicide Attempt", "Suicidal Ideation", "Intentional Self-Injury",
              "Homicidality Group", "Suicidality Group", HS)
hm <- merge(CJ(drug = drug_ord, pt = pt_ord, analysis = unname(full_mono)), hm,
            by = c("drug", "pt", "analysis"), all.x = TRUE)
hm[is.na(a), a := 0L]
hm[a < 3, ror := NA_real_]                      # no estimate shown below 3 cases
# "NS" read as "not significant", which is a different and much stronger claim
# than the truth: below three cases nothing was estimated at all, so there is no
# result to call significant or otherwise. The cell says why it is empty.
hm[, label_ror := ifelse(a == 0, "0 cases", ifelse(a < 3, "n < 3", sprintf("%.1f", ror)))]
hm[, cons_ok := (full_consensus == TRUE)]
hm[, drug_f := factor(drug, levels = drug_ord)]
hm[, analysis_f := factor(analysis, levels = unname(full_mono))]
yc <- heatmap_ycoords(hm, pt_ord); hm <- yc$dt
# Twelve rows to eight columns. Fixing the panel to that ratio is what makes
# a cell square, and makes 3a, 3b and 3c one matrix rather than three
# differently proportioned grids: each panel has its own furniture below it -
# 3c's colour bar is captioned over two lines where the others' is one - so
# matching the canvas alone would still leave the cells slightly different.
ar <- length(yc$labels) / length(drug_ord)
nocons <- hm[!is.na(ror) & a >= 3 & cons_ok == FALSE]     # outlined: no full consensus

txt_col <- function(r) ifelse(!is.na(r) & r > 40, "white", "black")

# A bound of the interval, as wide as it needs to be and no wider. Below 100
# a tenth still means something; at 920 it does not, and the two characters it
# costs are the difference between the interval sitting inside its cell and
# running into the borders. The widest interval here, 93-920 on three cases,
# is the one that needed it.
ci_num <- function(x) ifelse(abs(x) >= 100, sprintf("%.0f", x), sprintf("%.1f", x))

heat_theme <- function(base = 10, ratio = ar) theme_pub(base) + theme(
  text = element_text(family = FONT), aspect.ratio = ratio,
  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12, face = "bold", family = FONT),
  axis.text.y = element_text(size = 13, face = "bold", family = FONT),
  axis.title.x = element_blank(), axis.title.y = element_blank(),
  strip.background = element_blank(), plot.margin = margin(4, 4, 4, 4))

# One analysis, eight drugs, twelve outcomes. Built once and called twice, so
# 3a and 3b cannot drift apart.
#
# Three lines a cell: estimate, interval, case count.
#
# At 174 mm a cell is about 15.5 mm across. Set in 8 pt - the smallest the
# journal allows - "(93.2-920.3)" needs about 19 mm of that and collided with
# its neighbours, so the brackets go and the interval is written bare. That
# saves two characters, which is the difference between overflowing the cell
# and sitting inside it, and the row above makes it obvious what it is. The
# count loses its "n = " for the same reason. The three lines are spread
# relatively wide, because at this canvas size the type is large against the
# cell.
heat_panel <- function(which_analysis) {
  d  <- hm[analysis_f == which_analysis]
  nc <- nocons[analysis_f == which_analysis]
  ggplot(d, aes(drug_f, yc, fill = ror)) +
    geom_tile(colour = "grey80", linewidth = 0.5, width = 1, height = 1) +
    geom_tile(data = nc, colour = "#A32D2D", linewidth = 1.0, fill = NA, width = 0.9, height = 0.9) +
    # A cell holds three lines where there is an estimate, two where the count
    # is under the floor of three and there is none, and one where the drug and
    # outcome never occurred together. Each is centred on what it actually
    # holds, so a lone dash sits in the middle of its cell rather than hanging
    # from the top of it.
    geom_text(aes(y = yc + ifelse(a == 0, 0, ifelse(a < 3, 0.13, 0.24)),
                  label = label_ror, colour = txt_col(ror),
                  fontface = ifelse(!is.na(ror) & ror > 10, "bold", "plain")),
              size = 10 * PT, family = FONT) +
    geom_text(aes(y = yc - 0.02, colour = txt_col(ror),
                  label = ifelse(a >= 3 & !is.na(ror_ci_low),
                                 paste0(ci_num(ror_ci_low), "-", ci_num(ror_ci_high)), "")),
              size = 8 * PT, family = FONT) +
    geom_text(aes(y = yc - ifelse(a < 3, 0.18, 0.28), colour = txt_col(ror),
                  label = ifelse(a >= 1, format(a, big.mark = ",", trim = TRUE), "")),
              size = 8 * PT, family = FONT, fontface = "bold") +
    scale_fill_gradientn(colours = c("white", "#ffcccc", "#ff4d4d", "#e60000", "#990000"),
                         values = rescale(c(0, 15, 40, 70, 120)), name = "ROR", na.value = "grey97",
                         limits = c(0, 120), oob = squish) +
    scale_colour_identity() +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    scale_y_continuous(breaks = yc$breaks, labels = disp_lab(yc$labels), expand = expansion(add = 0)) +
    geom_vline(xintercept = c(1.5, 7.5), linetype = "dashed", colour = "grey60", linewidth = 0.6) +
    heat_theme() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = PT_AXIS,
                                     face = "bold", family = FONT),
          axis.text.y = element_text(size = PT_AXIS, face = "bold", family = FONT),
          legend.position = "bottom", legend.key.height = unit(0.35, "cm"),
          legend.key.width = unit(1.6, "cm"),
          legend.title = element_text(size = PT_AXIS, face = "bold"),
          legend.text = element_text(size = PT_AXIS))
}

# (c) normalized difference mono vs full, z-test where both reached full consensus
amp <- dcast(hm[, .(drug, pt, yc, analysis = ifelse(analysis == "Monopharmacy", "mono", "full"), ror, a, ror_ci_low, ror_ci_high, cons_ok)],
             drug + pt + yc ~ analysis, value.var = c("ror", "a", "ror_ci_low", "ror_ci_high", "cons_ok"))
amp[, amp := nd_index(ror_mono, ror_full)]
amp[!is.na(ror_full) & ror_full < 1, amp := NA]
amp[, aflag := (a_full < 3) | (a_mono < 3)]
amp[, both_consensus := (!is.na(cons_ok_full) & cons_ok_full) & (!is.na(cons_ok_mono) & cons_ok_mono)]
amp[, c("z", "p_value") := ror_z_test(ror_mono, ror_ci_low_mono, ror_ci_high_mono, ror_full, ror_ci_low_full, ror_ci_high_full)]
amp[!(both_consensus), c("z", "p_value") := NA_real_]
amp[, drug_f := factor(drug, levels = drug_ord)]
amp[, val := ifelse(!is.na(amp), ifelse(aflag, paste0(sprintf("%+.2f", amp), " (n<3)"), sprintf("%+.2f", amp)),
                    ifelse(!is.na(ror_full) & ror_full < 1, "ROR<1", ""))]
amp[, lab := ifelse(is.na(p_value), val, paste0(val, "\n", ifelse(p_value < 0.001, "p<0.001", paste0("p=", sprintf("%.3f", p_value)))))]
write_res(amp[, .(drug, pt, ror_full, ror_mono, a_full, a_mono, amp, both_consensus, z, p_value)], "fig3c_data.csv")

pC <- ggplot(amp, aes(drug_f, yc, fill = amp)) +
  geom_tile(colour = "grey70", linewidth = 0.3, width = 1, height = 1) +
  geom_text(aes(label = lab, colour = ifelse(!is.na(amp) & (amp > 0.45 | amp < -0.30),
                                             "white", "black"),
                fontface = ifelse(!is.na(amp) & abs(amp) >= 0.2, "bold", "plain")),
            size = 9 * PT, lineheight = 0.9, family = FONT) +
  # The normalized difference is signed: above zero the signal strengthens under
  # monopharmacy, below it the signal weakens. A ramp that runs from white to
  # red cannot show that - it puts -0.19 and +0.02 in nearly the same pale pink,
  # and leaves zero sitting at no particular place on the scale. So the fill
  # diverges from white at exactly zero: red for the strengthening side, keeping
  # faith with (a) and (b) where redder is a larger ROR, and blue for the other.
  # The palette is ColorBrewer's RdBu, which stays legible under the common
  # forms of colour blindness, where a red ramp separated only by lightness
  # does not.
  #
  # The limits are not symmetric, because the data are not: the differences run
  # from about -0.2 to +0.5. Zero is pinned by the colour positions rather than
  # by the midpoint of the range, so the asymmetry costs nothing but a little
  # unused blue.
  scale_fill_gradientn(colours = c("#2166AC", "#67A9CF", "#D1E5F0", "#F7F7F7",
                                   "#FDDBC7", "#EF8A62", "#B2182B"),
                       values = rescale(c(-0.4, -0.27, -0.13, 0, 0.2, 0.4, 0.6)),
                       name = "Normalized\ndifference",
                       na.value = "grey92", limits = c(-0.4, 0.6), oob = squish,
                       breaks = c(-0.3, 0, 0.3, 0.6)) +
  scale_colour_identity() +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  # 3c stands on its own rather than beside the heatmaps, so it carries its
  # own outcome labels rather than borrowing theirs.
  scale_y_continuous(breaks = yc$breaks, labels = disp_lab(yc$labels), expand = expansion(add = 0)) +
  geom_vline(xintercept = c(1.5, 7.5), linetype = "dashed", colour = "grey60", linewidth = 0.5) +
  labs(title = NULL, caption = NULL) +
  heat_theme(9) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6), axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = PT_AXIS,
                                   face = "bold", family = FONT),
        axis.text.y = element_text(size = PT_AXIS, face = "bold", family = FONT),
        legend.position = "bottom", legend.key.height = unit(0.35, "cm"), legend.key.width = unit(1.4, "cm"),
        legend.title = element_text(size = PT_AXIS, face = "bold"),
        legend.text = element_text(size = PT_AXIS))

# Figure 3 goes to the journal as three separate EPS files and its caption
# refers to (a), (b) and (c) throughout, so each panel has to carry its own
# letter. A label added in the manuscript document does not travel: production
# rebuilds the figure from the artwork, where it would not exist, leaving three
# unlabelled panels and a caption pointing at nothing.
panel_tag <- function(p, letter) {
  p + labs(tag = letter) +
    theme(plot.tag = element_text(size = 11, face = "bold", family = FONT, hjust = 0),
          plot.tag.position = c(0.012, 0.992))
}

# Three panels at the printed width, and tall enough that the cells come out
# square. Eight columns across 174 mm leave a cell 15.5 mm wide, so twelve rows
# of the same want a panel about 186 mm deep; with the rotated drug names and
# the colour bar below it that comes to 229 mm, inside the journal's 234 mm
# limit with a little to spare. Left shorter the cells are half again as wide
# as they are tall, which reads as squashed.
save_plot(panel_tag(heat_panel(full_mono[["full"]]), "(a)"),
          "fig3a_disproportionality_unrestricted", w = FIG_W, h = FIG_H)
save_plot(panel_tag(heat_panel(full_mono[["mono"]]), "(b)"),
          "fig3b_disproportionality_monopharmacy", w = FIG_W, h = FIG_H)
save_plot(panel_tag(pC, "(c)"), "fig3c_normalized_difference", w = FIG_W, h = FIG_H)

# ── Fig 4: externalizing vs internalizing ─────────────────────────────────────
# ND = (ROR homicide - ROR completed suicide) / sum, at PT and group level.
ei_src <- hm[analysis == "Unrestricted", .(drug, pt, ror, a, ror_ci_low, ror_ci_high, cons_ok)]
pick <- function(dr, p) ei_src[drug == dr & pt == p]
ei <- rbindlist(lapply(drug_ord, function(dr) rbindlist(lapply(list(
  c("Homicide / Completed Suicide (PT)", "Homicide", "Completed Suicide"),
  c("Homicidality / Suicidality (group)", "Homicidality Group", "Suicidality Group")), function(r) {
    h <- pick(dr, r[2]); s <- pick(dr, r[3])
    data.table(drug = dr, row = r[1], idx = nd_index(h$ror, s$ror),
               p_value = ror_z_test(h$ror, h$ror_ci_low, h$ror_ci_high, s$ror, s$ror_ci_low, s$ror_ci_high)$p,
               both_consensus = isTRUE(h$cons_ok) && isTRUE(s$cons_ok),
               flag = h$a < 3 | s$a < 3) }))))
# A p-value is shown only where both estimates behind it met signal criteria,
# which is the rule the Methods states and Fig 3c already applies. Without it
# this panel printed a p-value for varenicline / completed suicide, which fails
# two of the four criteria, and contradicted SI "SUMMARY".
ei[!(both_consensus), p_value := NA_real_]
ei[, drug := factor(drug, levels = drug_ord)]
ei[, row  := factor(row, levels = c("Homicidality / Suicidality (group)", "Homicide / Completed Suicide (PT)"))]
ei[, p_lab := fifelse(is.na(p_value), "", fifelse(p_value < 0.001, "p<0.001", paste0("p=", sprintf("%.3f", p_value))))]
ei[, lab := fifelse(is.na(idx), "", fifelse(flag, paste0(sprintf("%+.2f", idx), " (n<3)"),
       fifelse(p_lab == "", sprintf("%+.2f", idx), paste0(sprintf("%+.2f", idx), "\n", p_lab))))]
write_res(ei[, .(drug, row, idx, p_value, p_lab, both_consensus, flag_low_n = flag)], "fig4_data.csv")
lim <- ceiling(max(abs(ei$idx), na.rm = TRUE) * 20) / 20
fig4 <- ggplot(ei, aes(drug, row, fill = idx)) +
  geom_tile(colour = "grey85", linewidth = 0.5, width = 1, height = 1) +
  geom_text(aes(label = lab, colour = ifelse(!is.na(idx) & abs(idx) > lim * 0.65, "white", "black")),
            size = 8 * PT, lineheight = 0.95, fontface = "bold", family = FONT) +
  scale_fill_gradientn(colours = c("#5a9bd4", "#a8c6e5", "#e3ecf5", "#f5f5f5", "#f7c6c2", "#ee8f83", "#dc3e2f"),
                       limits = c(-lim, lim), oob = squish, na.value = "grey92", name = "Normalized\ndifference",
                       breaks = pretty_breaks(5)(c(-lim, lim))) +
  scale_colour_identity() + coord_fixed(ratio = 1) +
  labs(x = NULL, y = NULL) +
  theme_pub(11) +
  # Two rows of eight on a 174 mm page leaves a cell about 16 mm wide, so the
  # type comes down to the same 8 pt as Figure 3 and the colour bar moves under
  # the panel: on the right it took 35 mm of the width the cells needed, and
  # the labels ran into each other.
  theme(text = element_text(family = FONT),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = PT_AXIS),
        axis.text.y = element_text(size = PT_AXIS),
        panel.grid = element_blank(),
        legend.position = "bottom", legend.key.height = unit(0.35, "cm"),
        legend.key.width = unit(1.6, "cm"),
        # theme_pub() blanks the legend title, which left this colour bar unlabelled.
        legend.title = element_text(size = PT_AXIS, face = "bold"),
        legend.text = element_text(size = PT_AXIS))
save_plot(fig4, "fig4_externalizing_internalizing", w = FIG_W, h = 3.4)
writeLines(capture.output(sessionInfo()), res_path("sessionInfo.txt"))
log_msg("=== 07 COMPLETE | ", elapsed(t0), " ===")
