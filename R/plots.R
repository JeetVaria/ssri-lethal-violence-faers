# ==============================================================================
# R/plots.R — shared figure styling.
# ==============================================================================

# The journal asks for Helvetica or Arial. Naming a font the machine does not
# have produces "font family not found" and the figure silently falls back to
# the device default, which for cairo is DejaVu Sans - a Bitstream Vera
# derivative with visibly wider letterforms that will not pass as either.
#
# So the candidates are tried in order and the first one installed wins.
# Liberation Sans and Nimbus Sans are metric-compatible clones of Arial and
# Helvetica, drawn to the same widths, and are what a Linux box normally has
# in place of the originals. An empty string, the device default, is the last
# resort and is reported, because a figure set in it needs a font installed
# rather than a re-render.
FONT <- local({
  fams <- tryCatch(if (requireNamespace("systemfonts", quietly = TRUE))
                     systemfonts::system_fonts()$family else character(),
                   error = function(e) character())
  hit <- intersect(c("Arial", "Helvetica", "Liberation Sans", "Nimbus Sans"), fams)
  if (length(hit)) hit[1] else ""
})
if (!nzchar(FONT))
  message("No Arial, Helvetica or metric-compatible clone found. Figures will ",
          "use the device default, which the journal does not accept; install ",
          "fonts-liberation and re-run step 07.")

# Every figure is written twice: a PNG to place in the manuscript and read on
# screen, and an EPS for production. The journal takes vector art as EPS and
# keeps TIFF for photographs, and a chart written as a vector has no resolution
# to fall short of.
#
# EPS carries no alpha channel. Rather than drop a translucent layer, cairo
# rasterises just that layer at fallback_resolution - which is why it is set
# well above the journal's floor, for the sake of the one figure that has any
# (the Kaplan-Meier confidence bands). Everything else stays vector.
save_plot <- function(p, name, w = 14, h = 11, dpi = 300) {
  path <- file.path(FIG_DIR, paste0(name, ".png"))
  ggsave(path, plot = p, width = w, height = h, dpi = dpi, bg = "white")
  log_msg("Saved: ", basename(path))
  if (isTRUE(capabilities("cairo"))) {
    eps <- file.path(FIG_DIR, paste0(name, ".eps"))
    ggsave(eps, plot = p, width = w, height = h, bg = "white",
           device = grDevices::cairo_ps, fallback_resolution = 1200)
    log_msg("Saved: ", basename(eps))
  } else {
    log_msg("NOTE: this R has no cairo device, so ", name, ".eps was not written. ",
            "The PNG is unaffected.")
  }
}

theme_pub <- function(base = 12, family = "") {
  theme_bw(base_size = base, base_family = family) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.4),
      axis.text          = element_text(colour = "black", size = base),
      axis.title.x       = element_text(face = "italic", size = base),
      axis.title.y       = element_blank(),
      strip.background   = element_rect(fill = "grey92", colour = "grey70"),
      strip.text.y       = element_text(face = "bold", size = base - 1, angle = 0, hjust = 0),
      legend.position    = "bottom",
      legend.title       = element_blank(),
      legend.text        = element_text(size = base - 1),
      legend.key.width   = unit(2, "cm"),
      plot.title         = element_text(face = "bold", size = base + 2),
      plot.subtitle      = element_text(colour = "grey40", size = base - 1),
      plot.caption       = element_text(colour = "grey55", size = base - 2, hjust = 0))
}

# Upper-case the aggregate rows so they stand out from individual PTs.
disp_lab <- function(x) {
  x[x == "Homicidality Group"]           <- "HOMICIDALITY GROUP"
  x[x == "Suicidality Group"]            <- "SUICIDALITY GROUP"
  x[x == "Homicide + Completed Suicide"] <- "H-S COMPOSITE"
  x
}

# Numeric y positions (top row = largest y) and three-line cell labels
# (ROR / 95% CI / n) for the ROR heatmaps.
heatmap_ycoords <- function(dt, pt_ord_top) {
  ord  <- pt_ord_top[pt_ord_top %in% unique(dt$pt)]
  pos  <- seq_along(ord) - 1
  ymap <- data.table(pt = ord, yc = max(pos) - pos)
  dt   <- merge(dt, ymap, by = "pt", all.x = TRUE, sort = FALSE)
  # Returns the row order and its breaks only. Cell labels are composed inline
  # in 07_figures.R, with their own formatting, so they are not duplicated here
  # where they would be the wrong place to look.
  list(dt = dt, breaks = ymap$yc, labels = ord)
}
