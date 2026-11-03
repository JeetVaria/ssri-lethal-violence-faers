# Loads configuration, packages and every helper. Each analysis script starts
# with source("R/load_all.R"); run scripts from the repository root.
source("config.R")
for (f in c("setup", "definitions", "ingest", "signals", "ebgm", "groups", "domains", "codrugs", "plots"))
  source(file.path("R", paste0(f, ".R")))
