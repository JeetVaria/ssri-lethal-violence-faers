# ==============================================================================
# R/codrugs.R - classifying co-reported drugs on the Homicide / H-S case lists.
# ==============================================================================

codrug_classes <- list(
  "Antidepressants" = c(
    "fluoxetine", "prozac", "paroxetine", "paxil", "sertraline", "zoloft", "fluvoxamine", "luvox",
    "citalopram", "celexa", "escitalopram", "lexapro",                                         # SSRIs
    "venlafaxine", "effexor", "desvenlafaxine", "pristiq", "duloxetine", "cymbalta",
    "milnacipran", "levomilnacipran", "fetzima",                                               # SNRIs
    "amitriptyline", "nortriptyline", "imipramine", "clomipramine", "doxepin", "desipramine", "trimipramine",
    "phenelzine", "tranylcypromine", "selegiline", "isocarboxazid",                            # MAOIs
    "mirtazapine", "remeron", "bupropion", "wellbutrin", "trazodone", "vortioxetine", "trintellix",
    "vilazodone", "viibryd", "nefazodone"),
  "Anti-Anxiety Agents" = c("clonazepam", "klonopin", "alprazolam", "xanax", "lorazepam", "ativan",
    "diazepam", "valium", "temazepam", "restoril", "oxazepam", "chlordiazepoxide", "clorazepate",
    "midazolam", "triazolam", "flurazepam", "clobazam", "librium", "buspirone", "buspar",
    "hydroxyzine", "vistaril", "atarax"),
  "Z-drugs (Sedative-Hypnotics)" = c("zolpidem", "ambien", "eszopiclone", "lunesta", "zaleplon", "sonata",
    "zopiclone", "stilnox"),
  "Stimulants" = c("amphetamine", "adderall", "dextroamphetamine", "lisdexamfetamine", "vyvanse",
    "methamphetamine", "dexedrine", "methylphenidate", "ritalin", "concerta", "dexmethylphenidate",
    "focalin", "modafinil", "provigil", "armodafinil", "nuvigil"),
  "Antipsychotics" = c("quetiapine", "seroquel", "risperidone", "risperdal", "olanzapine", "zyprexa",
    "aripiprazole", "abilify", "haloperidol", "haldol", "clozapine", "clozaril", "ziprasidone", "geodon",
    "paliperidone", "invega", "lurasidone", "latuda", "chlorpromazine", "perphenazine", "fluphenazine"),
  "Mood Stabilizers" = c("lithium", "lamotrigine", "lamictal", "valproate", "divalproex", "depakote",
    "valproic acid", "carbamazepine", "tegretol", "oxcarbazepine", "trileptal"),
  "Analgesics" = c("acetaminophen", "paracetamol", "tylenol", "ibuprofen", "advil", "motrin", "naproxen",
    "aleve", "aspirin", "celecoxib", "celebrex", "diclofenac", "morphine", "oxycodone", "hydrocodone",
    "tramadol", "fentanyl", "codeine", "hydromorphone", "meperidine"),
  "Cardiovascular Drugs" = c("metoprolol", "atenolol", "propranolol", "carvedilol", "bisoprolol",
    "atorvastatin", "lipitor", "simvastatin", "zocor", "rosuvastatin", "crestor", "pravastatin",
    "lisinopril", "enalapril", "ramipril", "losartan", "valsartan", "furosemide", "lasix",
    "hydrochlorothiazide", "spironolactone", "amlodipine", "diltiazem", "verapamil", "warfarin",
    "apixaban", "eliquis", "rivaroxaban", "xarelto", "clopidogrel", "plavix"),
  "Antibiotics" = c("penicillin", "amoxicillin", "azithromycin", "zithromax", "ciprofloxacin", "cipro",
    "doxycycline", "cephalexin", "clindamycin", "metronidazole", "flagyl", "levofloxacin",
    "trimethoprim", "sulfamethoxazole", "bactrim", "vancomycin"),
  "Antihistamines" = c("cetirizine", "zyrtec", "diphenhydramine", "benadryl", "loratadine", "claritin",
    "fexofenadine", "allegra", "chlorpheniramine", "promethazine", "phenergan"),
  "Antidiabetics" = c("metformin", "insulin", "glipizide", "glyburide", "glimepiride", "sitagliptin",
    "januvia", "empagliflozin", "jardiance", "liraglutide", "victoza", "pioglitazone"),
  "Hormonal Therapies" = c("levothyroxine", "synthroid", "thyroid", "estrogen", "estradiol", "progesterone",
    "testosterone", "ethinyl estradiol", "norethindrone", "prednisone", "prednisolone",
    "hydrocortisone", "dexamethasone"),
  "Antivirals" = c("oseltamivir", "tamiflu", "acyclovir", "valacyclovir", "valtrex", "ritonavir",
    "tenofovir", "entecavir"),
  "Gastrointestinal Drugs" = c("omeprazole", "prilosec", "pantoprazole", "protonix", "esomeprazole",
    "nexium", "ranitidine", "famotidine", "pepcid", "ondansetron", "zofran", "metoclopramide",
    "reglan", "loperamide"))

# Recreational / illicit substances, reported separately from the drug classes.
substance_use_classes <- list(
  "Alcohol"  = c("alcohol", "ethanol"),
  "Cocaine"  = c("cocaine"),
  "Cannabis" = c("dronabinol", "cannabis", "marijuana", "tetrahydrocannabinol", "cannabidiol", "nabilone", "marinol"))

ssri_levels <- names(ssri_drugs)

# Optional RxNorm normalization (brand/salt -> ingredient, ingredient -> class).
# If data/rxnorm_name_map.csv and data/rxnorm_class_map.csv are present they
# take precedence over the keyword lists; otherwise the keyword lists are used.
RX_NAME  <- character(0)
RX_CLASS <- character(0)
if (file.exists("data/rxnorm_name_map.csv")) {
  nm <- fread("data/rxnorm_name_map.csv")
  lk <- setNames(toupper(trimws(nm$canonical_name)), toupper(trimws(nm$raw_name)))
  lk <- lk[!is.na(names(lk)) & nzchar(names(lk)) & !is.na(lk) & nzchar(lk)]
  RX_NAME <- lk[!duplicated(names(lk))]
}
if (file.exists("data/rxnorm_class_map.csv")) {
  cw   <- fread("data/rxnorm_class_map.csv")
  bare <- toupper(trimws(sub("\\s*\\(.*\\)$", "", as.character(cw$canonical_name))))
  cls  <- trimws(as.character(cw$assigned_class))
  keep <- !is.na(bare) & nzchar(bare) & !is.na(cls) & nzchar(cls)
  lk2  <- setNames(cls[keep], bare[keep])
  RX_CLASS <- lk2[!duplicated(names(lk2))]
}
rx_canonical <- function(v) {
  if (!length(v) || !length(RX_NAME)) return(v)
  hit <- RX_NAME[toupper(trimws(v))]
  ifelse(is.na(hit), v, unname(hit))
}

# Raw co-drug strings for one report -> cleaned, de-duplicated upper-case names.
# Salt, ester, hydrate and dosage-form words that follow an ingredient name.
# Without stripping these, "ZOLPIDEM" and "ZOLPIDEM TARTRATE" count as two drugs.
DRUG_SUFFIXES <- c(
  "HYDROCHLORIDE", "HCL", "TARTRATE", "BITARTRATE", "CARBONATE", "BICARBONATE",
  "FUMARATE", "MALEATE", "SUCCINATE", "SULFATE", "SULPHATE", "CITRATE",
  "PHOSPHATE", "ACETATE", "MESYLATE", "MESILATE", "BESYLATE", "TOSYLATE",
  "OXALATE", "NITRATE", "BROMIDE", "CHLORIDE", "HYDROBROMIDE", "EMBONATE",
  "PAMOATE", "DECANOATE", "VALERATE", "PROPIONATE", "SODIUM", "POTASSIUM",
  "CALCIUM", "MAGNESIUM", "HYDRATE", "MONOHYDRATE", "DIHYDRATE", "TRIHYDRATE",
  "ANHYDROUS", "TABLET", "TABLETS", "CAPSULE", "CAPSULES", "INJECTION",
  "SOLUTION", "ORAL", "TOPICAL")

# Misspellings and slash-joined names that the suffix rule cannot reach.
DRUG_ALIASES <- c(
  "CLONAZEPAN" = "CLONAZEPAM",
  "VENLAFAXINE/VENLAFAXINE HYDROCHLORIDE" = "VENLAFAXINE",
  "QUETIAPINE/QUETIAPINE FUMARATE" = "QUETIAPINE")

# Ingredient name: upper case, role tag, parentheticals, salt/ester/form words
# and punctuation removed. Applied after RxNorm resolution so the lookup still
# sees the name as recorded.
strip_drug_suffix <- function(v) {
  pat <- paste0("\\s+(", paste(DRUG_SUFFIXES, collapse = "|"), ")$")
  repeat {
    stripped <- sub(pat, "", v, perl = TRUE)
    if (identical(stripped, v)) break
    v <- stripped
  }
  trimws(v)
}

clean_drug_names <- function(v) {
  v <- v[!is.na(v) & v != "" & v != "(None Recorded)"]
  if (!length(v)) return(character(0))
  v <- sub("\\[[^]]*\\]\\s*$", "", v)         # trailing [ROLE]
  v <- toupper(trimws(v))
  v <- sub("\\.$", "", v)                     # trailing period
  v <- gsub("\\s*\\([^)]*\\)", " ", v)        # any "(clarifier)", not only trailing
  v <- trimws(gsub("\\s{2,}", " ", v))
  v <- ifelse(v == "ETHANOL", "ALCOHOL", v)
  v <- rx_canonical(v)
  v <- toupper(trimws(v))
  hit <- DRUG_ALIASES[v]; v <- ifelse(is.na(hit), v, unname(hit))
  v <- strip_drug_suffix(v)
  unique(v[nzchar(v)])
}

# Class of one co-drug. A co-reported SSRI other than the report's own index SSRI
# is "Other SSRIs". Keyword classes are checked in list order; first match wins.
classify_codrug <- function(codrug_name, index_drug = NA_character_) {
  cn <- tolower(codrug_name)
  idx <- if (is.na(index_drug) || index_drug == "") character(0) else tolower(trimws(strsplit(index_drug, ";", fixed = TRUE)[[1]]))
  if (cn %in% ssri_levels && !(cn %in% idx)) return("Other SSRIs")
  if (length(RX_CLASS)) {
    rc <- RX_CLASS[rx_canonical(toupper(trimws(codrug_name)))]
    if (!is.na(rc) && nzchar(rc) && !identical(unname(rc), "Other")) return(unname(rc))
  }
  for (cls in names(codrug_classes))
    if (any(vapply(codrug_classes[[cls]], function(k) grepl(k, cn, fixed = TRUE), logical(1)))) return(cls)
  "Other"
}

classify_substance_use <- function(codrug_name) {
  cn <- tolower(codrug_name)
  for (cls in names(substance_use_classes))
    if (any(vapply(substance_use_classes[[cls]], function(k) grepl(k, cn, fixed = TRUE), logical(1)))) return(cls)
  NA_character_
}

per_case_codrugs <- function(d) lapply(strsplit(as.character(d$codrug_names), "; ", fixed = TRUE), clean_drug_names)
per_case_classes <- function(per_case, idx_drug) Map(function(v, idx) {
  if (!length(v)) character(0) else unique(vapply(v, classify_codrug, character(1), index_drug = idx))
}, per_case, idx_drug)
per_case_substances <- function(per_case) lapply(per_case, function(v) {
  if (!length(v)) character(0) else unique(na.omit(vapply(v, classify_substance_use, character(1))))
})
