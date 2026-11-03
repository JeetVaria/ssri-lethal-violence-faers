# Runs the whole pipeline. From the repository root:  Rscript run_all.R
# Each step can also be run on its own (they share files in the results folder).
steps <- c("01_build_master", "02_signals_full", "03_signals_mono", "04_qa",
           "05_case_tables", "06_tto", "07_figures")
for (s in steps) {
  message("\n>>> ", s)
  source(file.path("analysis", paste0(s, ".R")), local = new.env())
}

# Supplementary "GLOBALS" sheet: the SSRI group characteristics, then the
# four-block file assembled from them and the step 01 / 05 outputs.
for (s in c("group_characteristics.R", "build_si_globals.R")) {
  message("\n>>> tools/", s)
  source(file.path("tools", s), local = new.env())
}

# Gather every estimate into one table and the Supplementary sheet, then check
# the numbers against each other. reconcile.R exits non-zero if a check fails.
for (s in c("combine_da.R", "reconcile.R")) {
  message("\n>>> ", s)
  source(s, local = new.env())
}
