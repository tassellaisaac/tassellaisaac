
# Isaac T, Simpson G, Smith L, Hasan M, Dambha-Miller H.
# Rate and timing of long-term condition accumulation among adults with multiple long-term conditions in England
# Aging and Health Research (2026)

# Main aims:
# 1) derive transition dates from 2nd to 10th LTC
# 2) summarise timing between successive LTC diagnoses
# 3) produce overall and stage-specific transition summaries
# 4) stratify by North/South and key sociodemographics
# 5) build random forest models for 10-year LTC accumulation
# 6) compare RF against linear regression baseline

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(lubridate)
  library(ranger)
  library(broom)
  library(fastDummies)
  library(forcats)
  library(scales)
})

set.seed(42)
options(datatable.print.nrows = 100)

# =========================================================
# 1. PATHS
# =========================================================

path_in  <- "data/cprd_mltc_person_level.csv"   
path_out <- "outputs_mltc_timing_rate"
dir.create(path_out, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# 2. EXPECTED INPUT STRUCTURE
# =========================================================
# One row per person with:
# - patid
# - age_baseline, sex, ethnicity, imd_quintile, region
# - 56 LTC date columns OR a long diagnosis table converted to wide dates


id_var <- "patid"

ltc_date_cols <- c(
  "hypertension_date", "depression_date", "osteoarthritis_date", "anxiety_date",
  "arrhythmia_date", "drug_alcohol_misuse_date", "cancer_date", "asthma_date",
  "coronary_heart_disease_date", "diabetes_date", "copd_date", "ckd_date",
  "heart_failure_date", "stroke_tia_date", "epilepsy_date", "dementia_date",
  "schizophrenia_bipolar_date", "learning_disability_date", "parkinsons_date",
  "multiple_sclerosis_date", "rheumatoid_arthritis_date", "osteoporosis_date",
  "chronic_pain_date", "psoriasis_eczema_date", "hearing_loss_date",
  "visual_impairment_date", "thyroid_disorder_date", "liver_disease_date",
  "ibs_date", "constipation_date", "migraine_date", "peripheral_vascular_disease_date",
  "chronic_sinusitis_date", "bronchiectasis_date", "inflammatory_bowel_disease_date",
  "diverticular_disease_date", "chronic_fatigue_date", "addison_date", "autism_date",
  "adhd_date", "eating_disorder_date", "personality_disorder_date",
  "substance_misuse_date", "alcohol_problem_date", "sleep_disorder_date",
  "obesity_date", "anaemia_date", "urinary_incontinence_date",
  "prostate_disorder_date", "endometriosis_date", "pcos_date", "gout_date",
  "chronic_pancreatitis_date", "congenital_condition_date", "infection_date",
  "haematological_condition_date"
)

core_covariates <- c(
  "age_baseline", "sex", "ethnicity", "imd_quintile", "region"
)

body_system_flags <- c(
  "cardiovascular", "metabolic_endocrine", "respiratory", "neurological",
  "cancers", "mental_behavioural", "musculoskeletal", "digestive",
  "urogenital", "haematological", "eye", "ear", "infections", "congenital"
)

# =========================================================
# 3. HELPERS
# =========================================================

assert_columns <- function(df, cols, df_name = deparse(substitute(df))) {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) {
    stop(sprintf("%s is missing columns: %s", df_name, paste(miss, collapse = ", ")))
  }
}

as_date_safe <- function(x) {
  if (inherits(x, "Date")) return(x)
  suppressWarnings(as.Date(x))
}

months_between <- function(date1, date2) {
  as.numeric(date2 - date1) / 30.4
}

days_between <- function(date1, date2) {
  as.numeric(date2 - date1)
}

safe_median <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
safe_mean   <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
safe_q25    <- function(x) if (all(is.na(x))) NA_real_ else quantile(x, 0.25, na.rm = TRUE, names = FALSE)
safe_q75    <- function(x) if (all(is.na(x))) NA_real_ else quantile(x, 0.75, na.rm = TRUE, names = FALSE)
safe_min    <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
safe_max    <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)

rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
mae  <- function(obs, pred) mean(abs(obs - pred), na.rm = TRUE)
r2_score <- function(obs, pred) {
  ss_res <- sum((obs - pred)^2, na.rm = TRUE)
  ss_tot <- sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)
  1 - ss_res / ss_tot
}

north_south_map <- function(region) {
  north_regions <- c(
    "North East", "North West", "Yorkshire and The Humber",
    "East Midlands", "West Midlands"
  )
  south_regions <- c(
    "East of England", "London", "South East", "South West"
  )
  dplyr::case_when(
    region %in% north_regions ~ "North",
    region %in% south_regions ~ "South",
    TRUE ~ NA_character_
  )
}

get_sorted_ltc_dates <- function(row_dates) {
  x <- as.Date(unlist(row_dates), origin = "1970-01-01")
  x <- x[!is.na(x)]
  sort(unique(x))
}

extract_transition_dates <- function(date_vec, max_k = 10) {
  out <- rep(as.Date(NA), max_k)
  n <- min(length(date_vec), max_k)
  if (n > 0) out[1:n] <- date_vec[1:n]
  out
}

# =========================================================
# 4. LOAD DATA
# =========================================================

cohort <- fread(path_in)
setDT(cohort)

assert_columns(cohort, c(id_var, core_covariates))

# convert date columns
for (cl in intersect(ltc_date_cols, names(cohort))) {
  cohort[, (cl) := as_date_safe(get(cl))]
}

# LTC date columns 
ltc_date_cols <- intersect(ltc_date_cols, names(cohort))
if (length(ltc_date_cols) < 2) {
  stop("Need at least 2 LTC diagnosis date columns to reconstruct transitions.")
}

# =========================================================
# 5. DERIVE ORDERED LTC TRANSITION DATES PER PERSON
# =========================================================

ordered_dates_list <- apply(cohort[, ..ltc_date_cols], 1, get_sorted_ltc_dates)
transition_mat <- t(sapply(ordered_dates_list, extract_transition_dates, max_k = 10))
transition_dt <- as.data.table(transition_mat)
setnames(transition_dt, paste0("T", 1:10, "_date"))

for (cl in names(transition_dt)) {
  transition_dt[, (cl) := as.Date(get(cl), origin = "1970-01-01")]
}

cohort <- cbind(cohort, transition_dt)

# number of observed LTCs 
cohort[, n_ltc_total := rowSums(!is.na(.SD)), .SDcols = ltc_date_cols]

# include adults with >=2 LTCs and complete 2nd diagnosis date
cohort <- cohort[!is.na(T2_date)]

cohort <- cohort[age_baseline >= 18 | is.na(age_baseline)]

# =========================================================
# 6. QUALITY FILTERS 
# =========================================================
# - exclude implausible time windows >99 years from second LTC diagnosis
# - focus on transitions up to 10th condition

cohort[, years_from_T2_to_T10 := as.numeric(T10_date - T2_date) / 365.25]
cohort[years_from_T2_to_T10 > 99, years_from_T2_to_T10 := NA_real_]

# =========================================================
# 7. BASELINE DERIVATION AT MLTC ONSET (2ND LTC)
# =========================================================

cohort[, baseline_date := T2_date]
cohort[, followup_to_T10_years := as.numeric(T10_date - T2_date) / 365.25]

if (!"age_baseline" %in% names(cohort) && "dob" %in% names(cohort)) {
  cohort[, dob := as_date_safe(dob)]
  cohort[, age_baseline := floor(as.numeric(baseline_date - dob) / 365.25)]
}

# Ethnicity handling 
if ("ethnicity" %in% names(cohort)) {
  cohort[, ethnicity := fct_explicit_na(as.factor(ethnicity), na_level = "Unknown")]
}
if ("sex" %in% names(cohort)) cohort[, sex := as.factor(sex)]
if ("imd_quintile" %in% names(cohort)) cohort[, imd_quintile := as.factor(imd_quintile)]
if ("region" %in% names(cohort)) cohort[, region := as.factor(region)]
cohort[, north_south := north_south_map(as.character(region))]
cohort[, north_south := as.factor(north_south)]

# =========================================================
# 8. DERIVE TRANSITION INTERVALS
# =========================================================

for (k in 2:9) {
  next_k <- k + 1
  dname <- paste0("gap_", k, "_to_", next_k, "_days")
  mname <- paste0("gap_", k, "_to_", next_k, "_months")
  cohort[, (dname) := days_between(get(paste0("T", k, "_date")), get(paste0("T", next_k, "_date")))]
  cohort[, (mname) := months_between(get(paste0("T", k, "_date")), get(paste0("T", next_k, "_date")))]
}

# all successive transition intervals pooled
all_gap_days <- unlist(cohort[, lapply(.SD, identity), .SDcols = patterns("^gap_[2-9]_to_[3-10]_days$")], use.names = FALSE)
all_gap_days <- all_gap_days[!is.na(all_gap_days) & all_gap_days >= 0]
all_gap_months <- all_gap_days / 30.4

summary_overall_transitions <- tibble(
  statistic = c("Minimum", "25th percentile", "Median", "Mean", "75th percentile", "Maximum"),
  value_days = c(
    safe_min(all_gap_days),
    safe_q25(all_gap_days),
    safe_median(all_gap_days),
    safe_mean(all_gap_days),
    safe_q75(all_gap_days),
    safe_max(all_gap_days)
  ),
  value_months = c(
    safe_min(all_gap_months),
    safe_q25(all_gap_months),
    safe_median(all_gap_months),
    safe_mean(all_gap_months),
    safe_q75(all_gap_months),
    safe_max(all_gap_months)
  )
)
write.csv(summary_overall_transitions, file.path(path_out, "table2_overall_transition_summary.csv"), row.names = FALSE)

# =========================================================
# 9. STAGE-SPECIFIC TRANSITION TABLE (2->3 ... 9->10)
# =========================================================

stage_summary <- map_dfr(2:9, function(k) {
  next_k <- k + 1
  v <- cohort[[paste0("gap_", k, "_to_", next_k, "_months")]]
  v <- v[!is.na(v) & v >= 0]
  tibble(
    transition_stage = paste0(k, " -> ", next_k, " conditions"),
    p25_months = safe_q25(v),
    median_months = safe_median(v),
    mean_months = safe_mean(v),
    p75_months = safe_q75(v),
    max_months = safe_max(v),
    n_adults = length(v)
  )
})
write.csv(stage_summary, file.path(path_out, "table3_stage_specific_transition_summary.csv"), row.names = FALSE)

# =========================================================
# 10. BASELINE TABLE
# =========================================================

cohort[, followup_years := as.numeric(pmax(T10_date, T2_date, na.rm = TRUE) - T2_date) / 365.25]
cohort[, reached_5plus := as.integer(!is.na(T5_date))]
cohort[, reached_10 := as.integer(!is.na(T10_date))]
cohort[, baseline_ltc_count := 2L]
cohort[, median_ltc_total_proxy := n_ltc_total]

baseline_table <- list(
  N = nrow(cohort),
  age_mean = mean(cohort$age_baseline, na.rm = TRUE),
  age_sd = sd(cohort$age_baseline, na.rm = TRUE),
  female_n = sum(as.character(cohort$sex) == "Female", na.rm = TRUE),
  female_pct = mean(as.character(cohort$sex) == "Female", na.rm = TRUE) * 100,
  male_n = sum(as.character(cohort$sex) == "Male", na.rm = TRUE),
  male_pct = mean(as.character(cohort$sex) == "Male", na.rm = TRUE) * 100,
  median_ltc_total = median(cohort$n_ltc_total, na.rm = TRUE),
  q1_ltc_total = quantile(cohort$n_ltc_total, 0.25, na.rm = TRUE),
  q3_ltc_total = quantile(cohort$n_ltc_total, 0.75, na.rm = TRUE),
  median_followup_years = median(cohort$followup_years, na.rm = TRUE),
  q1_followup_years = quantile(cohort$followup_years, 0.25, na.rm = TRUE),
  q3_followup_years = quantile(cohort$followup_years, 0.75, na.rm = TRUE),
  reached_5plus_n = sum(cohort$reached_5plus == 1, na.rm = TRUE),
  reached_5plus_pct = mean(cohort$reached_5plus == 1, na.rm = TRUE) * 100,
  reached_10_n = sum(cohort$reached_10 == 1, na.rm = TRUE),
  reached_10_pct = mean(cohort$reached_10 == 1, na.rm = TRUE) * 100
)
write.csv(as.data.frame(baseline_table), file.path(path_out, "table1_overall_baseline.csv"), row.names = FALSE)

ethnicity_tab <- cohort %>%
  as.data.frame() %>%
  count(ethnicity) %>%
  mutate(percent = 100 * n / sum(n))
write.csv(ethnicity_tab, file.path(path_out, "table1_ethnicity_breakdown.csv"), row.names = FALSE)

imd_tab <- cohort %>%
  as.data.frame() %>%
  count(imd_quintile) %>%
  mutate(percent = 100 * n / sum(n))
write.csv(imd_tab, file.path(path_out, "table1_imd_breakdown.csv"), row.names = FALSE)

# =========================================================
# 11. FIGURE 2: DISTRIBUTION OF TIME TO NEXT TRANSITION
# =========================================================

fig2_df <- tibble(
  gap_days = all_gap_days,
  gap_months = all_gap_months,
  gap_years = all_gap_days / 365.25
)

p2a <- ggplot(fig2_df %>% filter(gap_months <= 36), aes(x = gap_months)) +
  geom_histogram(bins = 80) +
  geom_vline(xintercept = c(6, 12, 24), linetype = 2) +
  labs(
    x = "Time to next transition (months)",
    y = "Count",
    title = "Distribution of time to next LTC transition (0-36 months)"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(path_out, "figure2a_gap_distribution_raw.png"), p2a, width = 8, height = 5, dpi = 300)

p2b <- ggplot(fig2_df %>% filter(gap_years > 0), aes(x = gap_years)) +
  geom_histogram(bins = 80) +
  scale_x_log10(labels = label_number()) +
  geom_vline(xintercept = c(3, 5, 10), linetype = 2) +
  labs(
    x = "Time to next transition (years, log scale)",
    y = "Count",
    title = "Distribution of time to next LTC transition (log scale)"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(path_out, "figure2b_gap_distribution_log.png"), p2b, width = 8, height = 5, dpi = 300)

# =========================================================
# 12. FIGURE 3: TIME BETWEEN SUCCESSIVE LTC DIAGNOSES
# =========================================================

fig3_df <- stage_summary %>%
  mutate(stage_num = 2:9)

p3 <- ggplot(fig3_df, aes(x = stage_num, y = median_months)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = p25_months, ymax = p75_months), alpha = 0.2) +
  scale_x_continuous(
    breaks = 2:9,
    labels = c("2→3", "3→4", "4→5", "5→6", "6→7", "7→8", "8→9", "9→10")
  ) +
  labs(
    x = "Transition stage",
    y = "Time to next LTC (months)",
    title = "Time between successive long-term condition diagnoses"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(path_out, "figure3_stage_transition_times.png"), p3, width = 8, height = 5, dpi = 300)

# =========================================================
# 13. NORTH-SOUTH STRATIFIED ANALYSES
# =========================================================

north_south_stage <- map_dfr(levels(na.omit(cohort$north_south)), function(ns) {
  sub <- cohort[north_south == ns]
  map_dfr(2:9, function(k) {
    next_k <- k + 1
    v <- sub[[paste0("gap_", k, "_to_", next_k, "_months")]]
    v <- v[!is.na(v) & v >= 0]
    tibble(
      north_south = ns,
      transition_stage = paste0(k, "->", next_k),
      median_months = safe_median(v),
      p25_months = safe_q25(v),
      p75_months = safe_q75(v),
      n = length(v)
    )
  })
})
write.csv(north_south_stage, file.path(path_out, "supp_north_south_stage_summaries.csv"), row.names = FALSE)

# =========================================================
# 14. SOCIODEMOGRAPHIC STRATIFIED TRANSITION SUMMARIES
# =========================================================

make_stratified_transition_table <- function(df, strat_var, out_file) {
  levs <- unique(df[[strat_var]])
  levs <- levs[!is.na(levs)]

  out <- map_dfr(levs, function(g) {
    sub <- df[get(strat_var) == g]
    map_dfr(2:9, function(k) {
      next_k <- k + 1
      v <- sub[[paste0("gap_", k, "_to_", next_k, "_months")]]
      v <- v[!is.na(v) & v >= 0]
      tibble(
        group = as.character(g),
        strat_var = strat_var,
        transition_stage = paste0(k, "->", next_k),
        median_months = safe_median(v),
        p25_months = safe_q25(v),
        p75_months = safe_q75(v),
        n = length(v)
      )
    })
  })

  write.csv(out, file.path(path_out, out_file), row.names = FALSE)
}

if ("sex" %in% names(cohort)) make_stratified_transition_table(cohort, "sex", "supp_transition_by_sex.csv")
if ("ethnicity" %in% names(cohort)) make_stratified_transition_table(cohort, "ethnicity", "supp_transition_by_ethnicity.csv")
if ("imd_quintile" %in% names(cohort)) make_stratified_transition_table(cohort, "imd_quintile", "supp_transition_by_imd.csv")

# age groups for supplement
cohort[, age_group := cut(age_baseline,
                          breaks = c(18, 40, 50, 60, 70, 80, Inf),
                          right = FALSE,
                          labels = c("18-39", "40-49", "50-59", "60-69", "70-79", "80+"))]
make_stratified_transition_table(cohort, "age_group", "supp_transition_by_age_group.csv")

# =========================================================
# 15. BODY SYSTEM FEATURE ENGINEERING
# =========================================================
if (!all(body_system_flags %in% names(cohort))) {
  condition_to_system <- list(
    cardiovascular = c("hypertension_date", "arrhythmia_date", "coronary_heart_disease_date", "heart_failure_date", "stroke_tia_date", "peripheral_vascular_disease_date"),
    metabolic_endocrine = c("diabetes_date", "thyroid_disorder_date", "obesity_date", "addison_date", "pcos_date"),
    respiratory = c("asthma_date", "copd_date", "bronchiectasis_date", "chronic_sinusitis_date"),
    neurological = c("epilepsy_date", "dementia_date", "parkinsons_date", "multiple_sclerosis_date", "migraine_date"),
    cancers = c("cancer_date"),
    mental_behavioural = c("depression_date", "anxiety_date", "schizophrenia_bipolar_date", "personality_disorder_date", "eating_disorder_date", "drug_alcohol_misuse_date", "substance_misuse_date", "alcohol_problem_date", "sleep_disorder_date", "adhd_date", "autism_date"),
    musculoskeletal = c("osteoarthritis_date", "rheumatoid_arthritis_date", "osteoporosis_date", "chronic_pain_date", "gout_date"),
    digestive = c("ibs_date", "liver_disease_date", "inflammatory_bowel_disease_date", "diverticular_disease_date", "constipation_date", "chronic_pancreatitis_date"),
    urogenital = c("ckd_date", "urinary_incontinence_date", "prostate_disorder_date", "endometriosis_date"),
    haematological = c("anaemia_date", "haematological_condition_date"),
    eye = c("visual_impairment_date"),
    ear = c("hearing_loss_date"),
    infections = c("infection_date"),
    congenital = c("congenital_condition_date", "learning_disability_date")
  )

  for (nm in names(condition_to_system)) {
    cols <- intersect(condition_to_system[[nm]], names(cohort))
    if (length(cols) > 0) {
      cohort[, (nm) := as.integer(rowSums(!is.na(.SD) & (.SD <= baseline_date), na.rm = TRUE) > 0), .SDcols = cols]
    } else {
      cohort[, (nm) := 0L]
    }
  }
}

# =========================================================
# 16. FEATURE ENGINEERING FOR PREDICTIVE MODEL
# =========================================================

# Number of new LTCs within 10 years after T2 baseline
cohort[, T2_plus_10y := baseline_date + years(10)]

# count LTCs whose first date lies within (T2_date, T2_date + 10y]
cohort[, new_ltcs_10y_after_T2 := apply(.SD, 1, function(x) {
  x <- as.Date(x, origin = "1970-01-01")
  base <- x[match("T2_date", c(names(.SD), "T2_date"))]
  NA_integer_
}), .SDcols = ltc_date_cols]

# safer rowwise implementation
new_counts <- vector("numeric", nrow(cohort))
gap_prev_vec <- vector("numeric", nrow(cohort))
for (i in seq_len(nrow(cohort))) {
  ltc_dates <- get_sorted_ltc_dates(cohort[i, ..ltc_date_cols])
  t2 <- cohort$T2_date[i]
  upper <- t2 + years(10)
  new_counts[i] <- sum(ltc_dates > t2 & ltc_dates <= upper, na.rm = TRUE)
  gap_prev_vec[i] <- if (!is.na(cohort$T1_date[i]) && !is.na(cohort$T2_date[i])) months_between(cohort$T1_date[i], cohort$T2_date[i]) else NA_real_
}
cohort[, new_ltcs_10y_after_T2 := new_counts]
cohort[, gap_prev := gap_prev_vec]

# optional stage-specific outcomes: number of new LTCs within 10 years after Tk
for (k in 2:5) {
  out_var <- paste0("new_ltcs_10y_after_T", k)
  tmp <- vector("numeric", nrow(cohort))
  for (i in seq_len(nrow(cohort))) {
    ltc_dates <- get_sorted_ltc_dates(cohort[i, ..ltc_date_cols])
    tk <- cohort[[paste0("T", k, "_date")]][i]
    if (is.na(tk)) {
      tmp[i] <- NA_real_
    } else {
      tmp[i] <- sum(ltc_dates > tk & ltc_dates <= (tk + years(10)), na.rm = TRUE)
    }
  }
  cohort[, (out_var) := tmp]
}

# =========================================================
# 17. REPRESENTATIVE SUBSAMPLE FOR MACHINE LEARNING
# =========================================================

model_df <- cohort %>%
  as.data.frame() %>%
  select(any_of(c(
    id_var, "new_ltcs_10y_after_T2", "age_baseline", "sex", "ethnicity", "imd_quintile",
    "north_south", "gap_prev", body_system_flags
  ))) %>%
  filter(!is.na(new_ltcs_10y_after_T2), !is.na(age_baseline), !is.na(sex), !is.na(imd_quintile))

if (nrow(model_df) > 200000) {
  set.seed(42)
  model_df <- model_df[sample(seq_len(nrow(model_df)), 200000), ]
}
write.csv(model_df[1:min(nrow(model_df), 1000), ], file.path(path_out, "sample_model_input_preview.csv"), row.names = FALSE)

# compare sample vs full cohort for representativeness if full cohort larger
repr_table <- tibble(
  metric = c("mean_age", "female_pct", "mean_new_ltcs_10y"),
  full = c(
    mean(cohort$age_baseline, na.rm = TRUE),
    mean(as.character(cohort$sex) == "Female", na.rm = TRUE) * 100,
    mean(cohort$new_ltcs_10y_after_T2, na.rm = TRUE)
  ),
  sample = c(
    mean(model_df$age_baseline, na.rm = TRUE),
    mean(as.character(model_df$sex) == "Female", na.rm = TRUE) * 100,
    mean(model_df$new_ltcs_10y_after_T2, na.rm = TRUE)
  )
)
write.csv(repr_table, file.path(path_out, "supp_sample_representativeness.csv"), row.names = FALSE)

# =========================================================
# 18. ONE-HOT ENCODING
# =========================================================

cat_vars <- c("sex", "ethnicity", "imd_quintile", "north_south")
cat_vars <- intersect(cat_vars, names(model_df))

model_df_enc <- fastDummies::dummy_cols(
  model_df,
  select_columns = cat_vars,
  remove_first_dummy = FALSE,
  remove_selected_columns = TRUE
)

# =========================================================
# 19. TRAIN / TEST SPLIT
# =========================================================

set.seed(42)
idx <- sample(seq_len(nrow(model_df_enc)), size = floor(0.8 * nrow(model_df_enc)))
train_df <- model_df_enc[idx, ]
test_df  <- model_df_enc[-idx, ]

outcome_var <- "new_ltcs_10y_after_T2"
predictor_vars <- setdiff(names(model_df_enc), c(id_var, outcome_var))

train_x <- train_df[, predictor_vars, drop = FALSE]
test_x  <- test_df[, predictor_vars, drop = FALSE]
train_y <- train_df[[outcome_var]]
test_y  <- test_df[[outcome_var]]

# =========================================================
# 20. RANDOM FOREST REGRESSION
# =========================================================

rf_fit <- ranger(
  formula = as.formula(paste(outcome_var, "~ .")),
  data = train_df[, c(outcome_var, predictor_vars), drop = FALSE],
  num.trees = 500,
  mtry = 5,
  importance = "impurity",
  seed = 42
)

rf_pred <- predict(rf_fit, data = test_df[, predictor_vars, drop = FALSE])$predictions
rf_perf <- tibble(
  model = "random_forest",
  rmse = rmse(test_y, rf_pred),
  mae = mae(test_y, rf_pred),
  r2 = r2_score(test_y, rf_pred)
)
write.csv(rf_perf, file.path(path_out, "rf_performance.csv"), row.names = FALSE)

rf_importance <- tibble(
  feature = names(rf_fit$variable.importance),
  importance = as.numeric(rf_fit$variable.importance)
) %>%
  arrange(desc(importance))
write.csv(rf_importance, file.path(path_out, "rf_variable_importance.csv"), row.names = FALSE)

# calibration by deciles of predicted risk
calib_df <- tibble(obs = test_y, pred = rf_pred) %>%
  mutate(decile = ntile(pred, 10)) %>%
  group_by(decile) %>%
  summarise(observed = mean(obs, na.rm = TRUE), predicted = mean(pred, na.rm = TRUE), .groups = "drop")
write.csv(calib_df, file.path(path_out, "rf_calibration_deciles.csv"), row.names = FALSE)

# =========================================================
# 21. LINEAR REGRESSION BASELINE
# =========================================================

lm_fit <- lm(as.formula(paste(outcome_var, "~ .")), data = train_df[, c(outcome_var, predictor_vars), drop = FALSE])
lm_pred <- predict(lm_fit, newdata = test_df[, predictor_vars, drop = FALSE])

lm_perf <- tibble(
  model = "linear_regression",
  rmse = rmse(test_y, lm_pred),
  mae = mae(test_y, lm_pred),
  r2 = r2_score(test_y, lm_pred)
)
write.csv(lm_perf, file.path(path_out, "lm_performance.csv"), row.names = FALSE)

perf_compare <- bind_rows(rf_perf, lm_perf)
write.csv(perf_compare, file.path(path_out, "model_performance_comparison.csv"), row.names = FALSE)

# =========================================================
# 22. STAGE-SPECIFIC RANDOM FOREST MODELS (T2, T3, T4, T5)
# =========================================================

stage_model_results <- list()
stage_importance_results <- list()

for (k in 2:5) {
  outcome_k <- paste0("new_ltcs_10y_after_T", k)
  base_k <- paste0("T", k, "_date")
  prev_gap_k <- if (k >= 2) paste0("gap_", k - 1, "_to_", k, "_months") else NA_character_

  if (!outcome_k %in% names(cohort)) next

  stage_df <- cohort %>%
    as.data.frame() %>%
    filter(!is.na(.data[[base_k]]), !is.na(.data[[outcome_k]])) %>%
    select(any_of(c(
      id_var, outcome_k, "age_baseline", "sex", "ethnicity", "imd_quintile",
      "north_south", prev_gap_k, body_system_flags
    )))

  if (nrow(stage_df) < 1000) next
  if (nrow(stage_df) > 200000) {
    set.seed(42)
    stage_df <- stage_df[sample(seq_len(nrow(stage_df)), 200000), ]
  }

  cat_vars_stage <- intersect(c("sex", "ethnicity", "imd_quintile", "north_south"), names(stage_df))
  stage_enc <- fastDummies::dummy_cols(
    stage_df,
    select_columns = cat_vars_stage,
    remove_first_dummy = FALSE,
    remove_selected_columns = TRUE
  )

  set.seed(42)
  idx_stage <- sample(seq_len(nrow(stage_enc)), size = floor(0.8 * nrow(stage_enc)))
  train_stage <- stage_enc[idx_stage, ]
  test_stage  <- stage_enc[-idx_stage, ]

  pred_vars_stage <- setdiff(names(stage_enc), c(id_var, outcome_k))

  rf_stage <- ranger(
    formula = as.formula(paste(outcome_k, "~ .")),
    data = train_stage[, c(outcome_k, pred_vars_stage), drop = FALSE],
    num.trees = 500,
    mtry = 5,
    importance = "impurity",
    seed = 42
  )

  pred_stage <- predict(rf_stage, data = test_stage[, pred_vars_stage, drop = FALSE])$predictions

  stage_model_results[[paste0("T", k)]] <- tibble(
    stage = paste0("T", k),
    rmse = rmse(test_stage[[outcome_k]], pred_stage),
    mae = mae(test_stage[[outcome_k]], pred_stage),
    r2 = r2_score(test_stage[[outcome_k]], pred_stage)
  )

  stage_importance_results[[paste0("T", k)]] <- tibble(
    stage = paste0("T", k),
    feature = names(rf_stage$variable.importance),
    importance = as.numeric(rf_stage$variable.importance)
  ) %>% arrange(desc(importance))
}

stage_perf_df <- bind_rows(stage_model_results)
if (nrow(stage_perf_df) > 0) write.csv(stage_perf_df, file.path(path_out, "supp_stage_specific_rf_performance.csv"), row.names = FALSE)

stage_importance_df <- bind_rows(stage_importance_results)
if (nrow(stage_importance_df) > 0) write.csv(stage_importance_df, file.path(path_out, "supp_stage_specific_rf_importance.csv"), row.names = FALSE)

# =========================================================
# 23. FIGURE 4: RANKED VARIABLE IMPORTANCE
# =========================================================

fig4_df <- rf_importance %>% slice_max(order_by = importance, n = 20)

p4 <- ggplot(fig4_df, aes(x = reorder(feature, importance), y = importance)) +
  geom_col() +
  coord_flip() +
  labs(
    x = NULL,
    y = "Variable importance",
    title = "Ranked variable importance for predicting 10-year LTC accumulation"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(path_out, "figure4_rf_variable_importance.png"), p4, width = 8, height = 6, dpi = 300)

# =========================================================
# 24. MOST FREQUENT LTCs AT BASELINE / OVERALL
# =========================================================

ltc_prevalence <- map_dfr(ltc_date_cols, function(cl) {
  tibble(
    ltc = cl,
    n = sum(!is.na(cohort[[cl]])),
    prevalence = mean(!is.na(cohort[[cl]]))
  )
}) %>% arrange(desc(prevalence))
write.csv(ltc_prevalence, file.path(path_out, "supp_ltc_prevalence_ranked.csv"), row.names = FALSE)

# =========================================================
# 25. FLOWCHART COUNTS
# =========================================================

flowchart_counts <- tibble(
  step = c(
    "Initial input rows",
    "Adults with >=2 LTCs",
    "Complete second LTC date",
    "Final analytic cohort"
  ),
  n = c(
    nrow(fread(path_in)),
    sum(rowSums(!is.na(fread(path_in)[, ..ltc_date_cols])) >= 2, na.rm = TRUE),
    nrow(cohort),
    nrow(cohort)
  )
)
write.csv(flowchart_counts, file.path(path_out, "figure1_flowchart_counts.csv"), row.names = FALSE)

# =========================================================
# 26. SESSION INFO
# =========================================================

capture.output(sessionInfo(), file = file.path(path_out, "sessionInfo.txt"))
message("MLTC timing/rate analysis complete. Outputs saved to: ", path_out)
