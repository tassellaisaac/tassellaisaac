****
**# "Variation in multimorbidity estimates due to clinical coding: a cross-sectional study of 7.2 million adults in England"

**
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(binom)
  library(xgboost)
  library(Matrix)
  library(fastDummies)
})

options(datatable.print.nrows = 50)


path_gold  <- "data/cprd_gold_patient_level.csv"
path_aurum <- "data/cprd_aurum_patient_level.csv"
path_hes   <- "data/hes_patient_level.csv"
path_demo  <- "data/demographics_master.csv"
path_out   <- "outputs"

dir.create(path_out, showWarnings = FALSE, recursive = TRUE)


condition_cols <- c(
  
  "hypertension", "depression", "asthma", "anxiety", "diabetes",
  "copd", "ckd", "cancer", "osteoarthritis", "coronary_heart_disease",
  "stroke_tia", "heart_failure", "atrial_fibrillation", "epilepsy",
  "dementia", "schizophrenia_bipolar", "learning_disability",
  "parkinsons", "multiple_sclerosis", "rheumatoid_arthritis",
  "osteoporosis", "chronic_pain", "psoriasis_eczema", "hearing_loss",
  "visual_impairment", "thyroid_disorder", "liver_disease", "ibs",
  "constipation", "migraine", "peripheral_vascular_disease",
  "chronic_sinusitis", "bronchiectasis", "inflammatory_bowel_disease",
  "diverticular_disease", "chronic_fatigue", "addison", "autism",
  "adhd", "eating_disorder", "personality_disorder", "substance_misuse",
  "alcohol_problem", "sleep_disorder", "obesity", "anaemia",
  "chronic_kidney_disease_advanced", "urinary_incontinence",
  "prostate_disorder", "endometriosis", "polycystic_ovary_syndrome",
  "gout", "chronic_pancreatitis", "permanent_arrhythmia"
)


required_covariates <- c(
  "patid", "age", "sex", "ethnicity", "imd_quintile", "region",
  "ae_any", "outpatient_any", "admission_any", "length_of_stay_days",
  "palliative_care", "hes_match"
)

**# =========================================================
# 2. HELPERS
# =========================================================**

assert_columns <- function(df, cols, df_name = deparse(substitute(df))) {
  missing_cols <- setdiff(cols, names(df))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "%s is missing required columns: %s",
      df_name,
      paste(missing_cols, collapse = ", ")
    ))
  }
}

safe01 <- function(x) {
  x <- fifelse(is.na(x), 0L, as.integer(x > 0))
  as.integer(x)
}

wilson_ci <- function(x, n, conf.level = 0.95) {
  ci <- binom::binom.wilson(x, n, conf.level = conf.level)
  data.frame(
    n_success = x,
    n_total = n,
    prevalence = ci$mean,
    lower = ci$lower,
    upper = ci$upper
  )
}

prop_test_2sample <- function(x1, n1, x2, n2) {
  out <- prop.test(x = c(x1, x2), n = c(n1, n2), correct = FALSE)
  tibble(
    estimate1 = unname(out$estimate[1]),
    estimate2 = unname(out$estimate[2]),
    p_value = out$p.value
  )
}

make_mm_flag <- function(df, cols, threshold = 2) {
  mm_count <- rowSums(df[, ..cols], na.rm = TRUE)
  as.integer(mm_count >= threshold)
}

rank_conditions_by_prevalence <- function(df, cols) {
  prev <- sapply(cols, function(cl) mean(df[[cl]] == 1, na.rm = TRUE))
  tibble(condition = cols, prevalence = as.numeric(prev)) %>%
    arrange(desc(prevalence))
}

stepwise_prevalence <- function(df, ranked_cols, threshold = 2, source_name = "GP") {
  n_total <- nrow(df)

  out <- map_dfr(seq_along(ranked_cols), function(k) {
    current_cols <- ranked_cols[1:k]
    mm_flag <- make_mm_flag(df, current_cols, threshold = threshold)
    x <- sum(mm_flag, na.rm = TRUE)
    ci <- wilson_ci(x, n_total)

    tibble(
      source = source_name,
      n_conditions = k,
      n_multimorbid = x,
      n_total = n_total,
      prevalence = ci$prevalence,
      lower = ci$lower,
      upper = ci$upper
    )
  })

  out
}

stepwise_prevalence_selected <- function(step_df, points = c(2, 5, 10, 20, 30, 54)) {
  step_df %>%
    filter(n_conditions %in% points)
}

subgroup_prevalence <- function(df, subgroup_var, condition_cols, threshold = 2, source_name = "GP") {
  subgroup_sym <- rlang::sym(subgroup_var)

  df %>%
    mutate(mm_flag = make_mm_flag(as.data.table(.), condition_cols, threshold)) %>%
    group_by(!!subgroup_sym) %>%
    summarise(
      n_total = n(),
      n_multimorbid = sum(mm_flag, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      prevalence = n_multimorbid / n_total,
      lower = binom::binom.wilson(n_multimorbid, n_total)$lower,
      upper = binom::binom.wilson(n_multimorbid, n_total)$upper,
      source = source_name,
      subgroup = subgroup_var
    ) %>%
    ungroup() %>%
    rename(level = !!subgroup_sym)
}

compare_gp_vs_hospital_by_subgroup <- function(gp_df, hes_df, subgroup_var) {
  gp_tab <- subgroup_prevalence(gp_df, subgroup_var, condition_cols, source_name = "GP")
  hes_tab <- subgroup_prevalence(hes_df, subgroup_var, condition_cols, source_name = "Hospital")

  merged <- gp_tab %>%
    select(level, gp_n_total = n_total, gp_n_mm = n_multimorbid, gp_prev = prevalence) %>%
    inner_join(
      hes_tab %>%
        select(level, hes_n_total = n_total, hes_n_mm = n_multimorbid, hes_prev = prevalence),
      by = "level"
    )

  tests <- pmap_dfr(
    list(merged$gp_n_mm, merged$gp_n_total, merged$hes_n_mm, merged$hes_n_total),
    ~ prop_test_2sample(..1, ..2, ..3, ..4)
  )

  bind_cols(tibble(subgroup = subgroup_var), merged, tests)
}

mean_abs_shap <- function(shap_matrix) {
  tibble(
    feature = colnames(shap_matrix),
    mean_abs_shap = colMeans(abs(shap_matrix), na.rm = TRUE)
  ) %>%
    arrange(desc(mean_abs_shap))
}

pearson_with_shap <- function(feature_df, shap_df) {
  common <- intersect(names(feature_df), names(shap_df))

  map_dfr(common, function(v) {
    x <- feature_df[[v]]
    y <- shap_df[[v]]
    ok <- complete.cases(x, y)

    if (sum(ok) < 10) {
      return(tibble(feature = v, correlation = NA_real_, p_value = NA_real_))
    }

    if (length(unique(x[ok])) < 2) {
      return(tibble(feature = v, correlation = NA_real_, p_value = NA_real_))
    }

    ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "pearson"))
    tibble(feature = v, correlation = unname(ct$estimate), p_value = ct$p.value)
  })
}

**# =========================================================
# 3. LOAD DATA
# =========================================================

gold  <- fread(path_gold)
aurum <- fread(path_aurum)
hes   <- fread(path_hes)
demo  <- fread(path_demo)


setnames(gold, old = intersect(names(gold), c("patient_id", "eid")), new = "patid", skip_absent = TRUE)
setnames(aurum, old = intersect(names(aurum), c("patient_id", "eid")), new = "patid", skip_absent = TRUE)
setnames(hes, old = intersect(names(hes), c("patient_id", "eid")), new = "patid", skip_absent = TRUE)
setnames(demo, old = intersect(names(demo), c("patient_id", "eid")), new = "patid", skip_absent = TRUE)

assert_columns(gold, c("patid", condition_cols))
assert_columns(aurum, c("patid", condition_cols))
assert_columns(hes, c("patid", condition_cols))
assert_columns(demo, required_covariates)

# =========================================================
# 4. COMBINE GOLD + AURUM GP DATA
# =========================================================

# Stack GP sources and de-duplicate at patient level using max() across condition flags
# Assumes same person may appear in one source only; if both, keep union of flags.

gp_all <- rbindlist(list(gold, aurum), use.names = TRUE, fill = TRUE)

for (cl in condition_cols) {
  gp_all[, (cl) := safe01(get(cl))]
  hes[, (cl) := safe01(get(cl))]
}

gp_all <- gp_all[
  , lapply(.SD, max, na.rm = TRUE),
  by = patid,
  .SDcols = condition_cols
]

hes_all <- hes[
  , lapply(.SD, max, na.rm = TRUE),
  by = patid,
  .SDcols = condition_cols
]

# =========================================================
# 5. BUILD MASTER GP / HOSPITAL / COMBINED DATASETS
# =========================================================

master <- as.data.table(demo)

master_gp <- merge(master, gp_all, by = "patid", all.x = TRUE)
master_hes <- merge(master, hes_all, by = "patid", all.x = TRUE)

for (cl in condition_cols) {
  master_gp[, (cl) := safe01(get(cl))]
  master_hes[, (cl) := safe01(get(cl))]
}

combined_flags <- merge(
  gp_all, hes_all,
  by = "patid", all = TRUE,
  suffixes = c("_gp", "_hes")
)

for (cl in condition_cols) {
  gp_col  <- paste0(cl, "_gp")
  hes_col <- paste0(cl, "_hes")
  combined_flags[, (cl) := pmax(
    fifelse(is.na(get(gp_col)), 0L, get(gp_col)),
    fifelse(is.na(get(hes_col)), 0L, get(hes_col))
  )]
}

combined_flags <- combined_flags[, c("patid", condition_cols), with = FALSE]
master_combined <- merge(master, combined_flags, by = "patid", all.x = TRUE)

for (cl in condition_cols) {
  master_combined[, (cl) := safe01(get(cl))]
}

# =========================================================
# 6. PRIMARY PREVALENCE ESTIMATES (>=2 CONDITIONS)
# =========================================================

master_gp[, mm_flag := make_mm_flag(.SD, condition_cols, threshold = 2), .SDcols = condition_cols]
master_hes[, mm_flag := make_mm_flag(.SD, condition_cols, threshold = 2), .SDcols = condition_cols]
master_combined[, mm_flag := make_mm_flag(.SD, condition_cols, threshold = 2), .SDcols = condition_cols]

primary_prev <- bind_rows(
  wilson_ci(sum(master_gp$mm_flag), nrow(master_gp)) %>% mutate(source = "GP"),
  wilson_ci(sum(master_hes$mm_flag), nrow(master_hes)) %>% mutate(source = "Hospital"),
  wilson_ci(sum(master_combined$mm_flag), nrow(master_combined)) %>% mutate(source = "Combined")
) %>%
  select(source, everything())

fwrite(primary_prev, file.path(path_out, "primary_multimorbidity_prevalence.csv"))

# =========================================================
# 7. STEPWISE INCLUSION ANALYSIS
# =========================================================

rank_gp <- rank_conditions_by_prevalence(master_gp, condition_cols)
rank_hes <- rank_conditions_by_prevalence(master_hes, condition_cols)
rank_combined <- rank_conditions_by_prevalence(master_combined, condition_cols)

step_gp <- stepwise_prevalence(master_gp, rank_gp$condition, source_name = "GP")
step_hes <- stepwise_prevalence(master_hes, rank_hes$condition, source_name = "Hospital")
step_combined <- stepwise_prevalence(master_combined, rank_combined$condition, source_name = "Combined")

step_all <- bind_rows(step_gp, step_hes, step_combined)
fwrite(step_all, file.path(path_out, "stepwise_prevalence_all_thresholds.csv"))

step_selected <- bind_rows(
  stepwise_prevalence_selected(step_gp),
  stepwise_prevalence_selected(step_hes),
  stepwise_prevalence_selected(step_combined)
)
fwrite(step_selected, file.path(path_out, "stepwise_prevalence_selected_points.csv"))

# Figure 1 style plot
p_step <- ggplot(step_all, aes(x = n_conditions, y = prevalence, colour = source)) +
  geom_line(linewidth = 1.1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = source), alpha = 0.12, colour = NA) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(
    x = "Number of conditions included",
    y = "Prevalence of multimorbidity (>=2 conditions)",
    title = "Stepwise inclusion of conditions by prevalence",
    colour = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

ggsave(file.path(path_out, "figure1_stepwise_prevalence.png"), p_step, width = 9, height = 6, dpi = 300)

# =========================================================
# 8. GP VS HOSPITAL PROPORTION TEST
# =========================================================

gp_vs_hospital <- prop_test_2sample(
  x1 = sum(master_gp$mm_flag), n1 = nrow(master_gp),
  x2 = sum(master_hes$mm_flag), n2 = nrow(master_hes)
)
write.csv(gp_vs_hospital, file.path(path_out, "gp_vs_hospital_prop_test.csv"), row.names = FALSE)

# =========================================================
# 9. SUBGROUP ANALYSES
# =========================================================

subgroups <- c("age_group", "sex", "ethnicity", "imd_quintile")

# create age groups if absent
if (!"age_group" %in% names(master_gp)) {
  master_gp[, age_group := cut(
    age,
    breaks = c(18, 40, 50, 60, 70, 80, Inf),
    right = FALSE,
    labels = c("18-39", "40-49", "50-59", "60-69", "70-79", "80+")
  )]
  master_hes[, age_group := cut(
    age,
    breaks = c(18, 40, 50, 60, 70, 80, Inf),
    right = FALSE,
    labels = c("18-39", "40-49", "50-59", "60-69", "70-79", "80+")
  )]
  master_combined[, age_group := cut(
    age,
    breaks = c(18, 40, 50, 60, 70, 80, Inf),
    right = FALSE,
    labels = c("18-39", "40-49", "50-59", "60-69", "70-79", "80+")
  )]
}

subgroup_prev_all <- bind_rows(
  map_dfr(subgroups, ~ subgroup_prevalence(master_gp, .x, condition_cols, source_name = "GP")),
  map_dfr(subgroups, ~ subgroup_prevalence(master_hes, .x, condition_cols, source_name = "Hospital")),
  map_dfr(subgroups, ~ subgroup_prevalence(master_combined, .x, condition_cols, source_name = "Combined"))
)

fwrite(subgroup_prev_all, file.path(path_out, "subgroup_prevalence.csv"))

subgroup_comparisons <- map_dfr(subgroups, ~ compare_gp_vs_hospital_by_subgroup(master_gp, master_hes, .x))
fwrite(subgroup_comparisons, file.path(path_out, "gp_vs_hospital_subgroup_comparisons.csv"))

# =========================================================
# 10. CHI-SQUARE TESTS ACROSS SUBGROUPS
# =========================================================

chisq_results <- map_dfr(subgroups, function(v) {
  gp_tab <- table(master_gp[[v]], master_gp$mm_flag)
  hes_tab <- table(master_hes[[v]], master_hes$mm_flag)
  comb_tab <- table(master_combined[[v]], master_combined$mm_flag)

  tibble(
    subgroup = v,
    source = c("GP", "Hospital", "Combined"),
    p_value = c(
      suppressWarnings(chisq.test(gp_tab)$p.value),
      suppressWarnings(chisq.test(hes_tab)$p.value),
      suppressWarnings(chisq.test(comb_tab)$p.value)
    )
  )
})

fwrite(chisq_results, file.path(path_out, "chisq_subgroup_differences.csv"))

# =========================================================
# 11. GP-ONLY OUTCOME FOR XGBOOST
# =========================================================

# Outcome definition used here:
# gp_only = 1 if multimorbid in GP but not multimorbid in hospital

ml_df <- merge(
  master[, ..required_covariates],
  master_gp[, .(patid, gp_mm_flag = mm_flag)],
  by = "patid",
  all.x = TRUE
) %>%
  left_join(master_hes[, .(patid, hes_mm_flag = mm_flag)] %>% as.data.frame(), by = "patid") %>%
  as.data.table()

ml_df[, gp_mm_flag := safe01(gp_mm_flag)]
ml_df[, hes_mm_flag := safe01(hes_mm_flag)]
ml_df[, gp_only := as.integer(gp_mm_flag == 1 & hes_mm_flag == 0)]

model_vars <- c(
  "age", "sex", "ethnicity", "imd_quintile", "region",
  "admission_any", "outpatient_any", "ae_any", "length_of_stay_days",
  "palliative_care", "hes_match"
)

ml_complete <- ml_df[, c("gp_only", model****
