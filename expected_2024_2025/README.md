# expected_2024_2025

Expected pipeline output for the 2024–2025 period, against which route 2 is compared:

```bash
FAERS_START_YEAR=2024 FAERS_END_YEAR=2025 FAERS_RESULTS_DIR=two_year Rscript run_all.R
Rscript tools/compare_to_reference.R expected_2024_2025 two_year
```

Generated from the same frozen code as `expected_2004_2025/`, so that a reviewer with two
years of FAERS quarters can test the whole pipeline end to end in about half an hour instead
of three and a half hours.

CSV files only. `compare_to_reference.R` compares the CSVs and merely counts the figures the
run produced, because a PNG is not comparable byte for byte across machines - font
rasterisation and library versions differ - and the run draws its own figures anyway. The
published figures are the ones in `expected_2004_2025/figures/`.
