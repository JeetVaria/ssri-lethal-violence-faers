# ==============================================================================
# R/domains.R - psychiatric symptom-domain keywords for the case-level review of
# co-reported PTs (Homicide and H-S case lists).
#
# Taxonomy aligned to DSM-5 chapters; matching is case-insensitive substring/regex
# on the PT text. Bipolar and Depressive are kept apart (separate DSM-5
# chapters). The Akathisia/EPS group is matched by keyword like every other
# group here - it is not the MedDRA "Extrapyramidal syndrome (SMQ)" member
# list, and it catches Restlessness and Tremor as well as Akathisia itself.
# Every term it actually matched is listed in SI "SYMPTOM DOMAINS".
# ==============================================================================

bipolar_kw <- c("\\bmania\\b", "\\bmanic\\b", "\\bhypomania\\b", "\\bhypomanic\\b", "\\bbipolar\\b",
                "mood elevat", "elevated mood", "grandios", "euphoric mood", "mood swing")
depressive_kw <- c("depression", "depressed mood", "dysphoria", "affective disorder",
                   "affect lability", "emotional disorder", "emotional distress", "anhedonia", "crying",
                   "activation syndrome", "decreased interest", "feeling of despair")
psychotic_broad_kw <- c("psychot", "psychosis", "hallucinat", "delusion", "paranoi", "schizophren",
                        "persecut", "disorgani", "thought disorder", "formication",
                        "blunted affect", "flat affect", "emotional poverty", "tangentiality")
dissociative_kw <- c("derealisation", "depersonalis", "dissociat")
anxiety_kw <- c("anxiety", "panic", "phobia", "post-traumatic stress", "posttraumatic stress",
                "ptsd", "acute stress", "obsessive-compulsive", "generalised anxiety",
                "generalized anxiety", "social phobia", "agoraphobia", "adjustment disorder")
personality_kw <- c("personality disorder", "borderline personality", "histrionic",
                    "antisocial personality", "narcissistic", "avoidant personality",
                    "dependent personality", "paranoid personality", "schizoid", "schizotypal")
neurocog_kw <- c("dementia", "delirium", "attention deficit", "adhd", "autis", "intellectual disab",
                 "cognitive impair", "confusional state", "amnesia", "memory impair",
                 "developmental disorder", "alzheimer")
akathisia_kw <- c("akathisia", "restlessness", "\\btremor\\b", "psychomotor hyperactivity",
                  "dyskinesia", "gait disturbance", "movement disorder")

broad_cat_levels <- c("Bipolar and Related Symptoms", "Depressive Symptoms",
                      "Schizophrenia Spectrum & Psychotic Symptoms",
                      "Anxiety, OCD, Trauma-Related and Dissociative Symptoms",
                      "Personality Symptoms", "Neurocognitive and Neurodevelopmental Symptoms",
                      "Akathisia/Extrapyramidal Symptoms", "Other")

# First matching domain wins; Personality is checked first so e.g. "Paranoid
# Personality Disorder" is not classed as psychotic.
pt_cat_broad <- function(pt) {
  p <- tolower(pt)
  hit <- function(kws) any(vapply(kws, function(k) grepl(k, p, perl = TRUE), logical(1)))
  if (hit(personality_kw))                return("Personality Symptoms")
  if (hit(bipolar_kw))                    return("Bipolar and Related Symptoms")
  if (hit(depressive_kw))                 return("Depressive Symptoms")
  if (hit(psychotic_broad_kw))            return("Schizophrenia Spectrum & Psychotic Symptoms")
  if (hit(c(anxiety_kw, dissociative_kw))) return("Anxiety, OCD, Trauma-Related and Dissociative Symptoms")
  if (hit(neurocog_kw))                   return("Neurocognitive and Neurodevelopmental Symptoms")
  if (hit(akathisia_kw))                  return("Akathisia/Extrapyramidal Symptoms")
  "Other"
}

# Narrower two-way tag (manic vs psychotic) reported alongside the broad domains.
manic_kw <- c("\\bmania\\b", "\\bmanic\\b", "\\bhypomania\\b", "\\bhypomanic\\b", "\\bbipolar\\b",
              "mood elevat", "elevated mood", "grandios", "euphoric mood")
psych_kw <- c("psychot", "psychosis", "hallucinat", "delusion", "paranoi", "schizophren",
              "delirium", "persecut", "disorgani", "thought disorder", "derealisation",
              "depersonalis", "formication")
pt_cat <- function(pt) {
  p <- tolower(pt)
  if (any(vapply(manic_kw, function(k) grepl(k, p, perl = TRUE), logical(1)))) return("Manic/hypomanic")
  if (any(vapply(psych_kw, function(k) grepl(k, p, fixed = TRUE), logical(1)))) return("Psychosis")
  "Other"
}
