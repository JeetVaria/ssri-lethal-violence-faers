# ==============================================================================
# tools/build_si_globals.R - assemble the "GLOBALS" sheet of the Supplementary
# Information from the four blocks the pipeline already produces.
#
# Reads  results/summary_master.csv              (step 01)
#        results/ssri_group_characteristics.csv  (tools/group_characteristics.R)
#        results/table1_homicide_ssri.csv        (step 05)
#        results/table1_hs_ssri.csv              (step 05)
# Writes results/SI_DATA_GLOBALS.csv
#
# strip.white = FALSE keeps the two-space indentation that marks sub-rows;
# na.strings = NULL keeps the literal "NA" in fda_deletion_list_n.
#
# Usage: Rscript tools/build_si_globals.R
# ==============================================================================
source("R/load_all.R")

rd <- function(f) fread(require_res(f), colClasses = "character",
                        strip.white = FALSE, na.strings = NULL)

sm <- rd("summary_master.csv")
b1 <- data.table(block  = "All FAERS reports",
                 metric = paste0(sm$section, ": ", sm$metric),
                 value  = sm$value)

gch <- rd("ssri_group_characteristics.csv")
b2  <- data.table(block = "SSRI group", metric = gch$metric, value = gch$value)

tbl1 <- function(file, label) {
  d <- rd(file)
  data.table(block = paste(label, "group"), metric = d$Characteristic, value = d$n)
}
b3 <- tbl1("table1_homicide_ssri.csv", "SSRI Homicide (unrestricted)")
b4 <- tbl1("table1_hs_ssri.csv",       "SSRI H-S (unrestricted)")

out <- rbindlist(list(b1, b2, b3, b4))
write_res(out, "SI_DATA_GLOBALS.csv")
log_msg("SI GLOBALS: ", nrow(out), " rows | all FAERS ", nrow(b1),
        " | SSRI group ", nrow(b2), " | homicide ", nrow(b3), " | H-S ", nrow(b4))
