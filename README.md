
**# Isaac T, Jat A, Smith L, Simpson G, Dambha-Miller H.
# Social care need and accumulation in the number of long-term conditions within multimorbidity:
# a cohort study of 7.2 million adults in England.
#**


# Aims:
# 1) quantify whether baseline social care need predicts faster accumulation in number of LTCs
# 2) estimate incident acquisition rate after index date
# 3) model time to next LTC and/or transition to more complex MLTC states
# 4) describe heterogeneity by social care need domains and sociodemographic factors

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(survival)
  library(splines)
  library(MASS)
  library(broom)
  library(broom.helpers)
  library(lubridate)
  library(sandwich)
  library(lmtest)
  library(emmeans)
  library(geepack)
  library(xgboost)
  library(fastDummies)
  library(Matrix)
})

options(datatable.print.nrows = 100)
set.seed(123)

# =========================================================
# 1. PATHS
# =========================================================

path_in  <- "data/df_mltc_yes.csv"   
path_out <- "outputs_social_care_need_mltc"
dir.create(path_out, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# 2. EXPECTED INPUT STRUCTURE
# =========================================================

# Core fields 
# - patid
# - index_date
# - d_crd / d_uts / d_tod
# - age, sex, ethnicity, imd_quintile, region
# - condition_count_baseline     (number of LTCs at index)
# - condition_count_end          (number of LTCs by end of follow-up)
# - next_ltc_date                (date of first newly acquired LTC after index)
# - death_date                   (optional)
# - plus social care variables/domains

# =========================================================
# 3. VARIABLE LISTS
# =========================================================

id_var           <- "patid"
index_date_var   <- "index_date"
end_date_var     <- "d_tod"
uts_date_var     <- "d_uts"
reg_date_var     <- "d_crd"
death_date_var   <- "death_date"

# demographic covariates
core_covariates <- c(
  "age", "sex", "ethnicity", "imd_quintile", "region"
)

# social care need domains
social_need_domains <- c(
  "hearing_disorders",
  "financial_challenges",
  "social_isolation",
  "mobility_impairments",
  "bereavement",
  "learning_disabilities",
  "residential_needs",
  "adl_limitations"
)

# individual social care variables
social_need_vars_optional <- c(
  "hearing_impairment", "visual_impairment", "housebound", "lives_alone",
  "carer_support", "financial_hardship", "wheelchair_need", "falls_need",
  "adl_bathing", "adl_dressing", "adl_feeding", "adl_toileting",
  "supported_housing", "care_home_need", "bereavement_need", "learning_disability_need"
)


count_baseline_var <- "condition_count_baseline"
count_end_var      <- "condition_count_end"
next_ltc_date_var  <- "next_ltc_date"


extra_covariates <- c(
  "consultation_rate", "hospital_admission_any", "ae_any", "outpatient_any"
)

# =========================================================
# 4. HELPERS
# =========================================================

assert_columns <- function(df, cols, df_name = deparse(substitute(df))) {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) {
    stop(sprintf("%s is missing columns: %s", df_name, paste(miss, collapse = ", ")))
  }
}

safe01 <- function(x) as.integer(fifelse(is.na(x), 0L, as.integer(x > 0)))

as_date_safe <- function(x) {
  if (inherits(x, "Date")) return(x)
  as.Date(x)
}

mode_value <- function(x) {
  ux <- unique(na.omit(x))
  ux[which.max(tabulate(match(x, ux)))]
}

calc_followup_years <- function(start_date, end_date) {
  as.numeric(difftime(end_date, start_date, units = "days")) / 365.25
}

make_binary_domain <- function(df, vars, out_name) {
  present <- intersect(vars, names(df))
  if (length(present) == 0) {
    warning(sprintf("No source variables found for domain %s", out_name))
    df[[out_name]] <- NA_integer_
  } else {
    df[[out_name]] <- as.integer(rowSums(df[, ..present], na.rm = TRUE) > 0)
  }
  df
}

robust_poisson <- function(formula, data) {
  fit <- glm(formula, family = poisson(link = "log"), data = data)
  vc  <- sandwich::vcovHC(fit, type = "HC0")
  out <- lmtest::coeftest(fit, vcov. = vc)
  list(fit = fit, robust = out, vcov = vc)
}

fit_nb_or_poisson <- function(formula, data) {
  pois <- glm(formula, family = poisson(), data = data)
  disp <- sum(residuals(pois, type = "pearson")^2, na.rm = TRUE) / pois$df.residual

  if (is.finite(disp) && disp > 1.5) {
    nb <- MASS::glm.nb(formula, data = data)
    return(list(model = nb, family = "negative binomial", dispersion = disp))
  } else {
    return(list(model = pois, family = "poisson", dispersion = disp))
  }
}

extract_rr <- function(model, conf.level = 0.95) {
  broom::tidy(model, conf.int = TRUE, conf.level = conf.level, exponentiate = TRUE)
}

extract_hr <- function(model, conf.level = 0.95) {
  broom::tidy(model, conf.int = TRUE, conf.level = conf.level, exponentiate = TRUE)
}

# =========================================================
# 5. LOAD COHORT
# =========================================================

cohort <- fread(path_in)
setDT(cohort)

# standardise date fields 
for (v in intersect(c(index_date_var, end_date_var, uts_date_var, reg_date_var, death_date_var, next_ltc_date_var), names(cohort))) {
  cohort[, (v) := as_date_safe(get(v))]
}

assert_columns(cohort, c(id_var, index_date_var, core_covariates))

# =========================================================
# 6. DOMAIN DERIVATION FROM INDIVIDUAL VARIABLES
# =========================================================


if (!all(social_need_domains %in% names(cohort))) {
  domain_map <- list(
    hearing_disorders   = c("hearing_impairment", "deafness_need", "hearing_aid_need"),
    financial_challenges = c("financial_hardship", "benefits_need", "debt_need"),
    social_isolation    = c("lives_alone", "social_isolation_need", "loneliness_need"),
    mobility_impairments = c("mobility_need", "wheelchair_need", "housebound", "falls_need"),
    bereavement         = c("bereavement_need"),
    learning_disabilities = c("learning_disability_need", "learning_disability"),
    residential_needs   = c("supported_housing", "care_home_need", "housing_need"),
    adl_limitations     = c("adl_bathing", "adl_dressing", "adl_feeding", "adl_toileting")
  )

  for (nm in names(domain_map)) {
    if (!nm %in% names(cohort)) {
      cohort <- make_binary_domain(cohort, domain_map[[nm]], nm)
    }
  }
}


for (v in intersect(social_need_domains, names(cohort))) {
  cohort[, (v) := safe01(get(v))]
}

# social need summary measures
cohort[, social_need_count := rowSums(.SD, na.rm = TRUE), .SDcols = intersect(social_need_domains, names(cohort))]
cohort[, any_social_need := as.integer(social_need_count >= 1)]
cohort[, social_need_cat := cut(
  social_need_count,
  breaks = c(-Inf, 0, 1, 2, Inf),
  labels = c("0", "1", "2", "3+")
)]

# =========================================================
# 7. OUTCOME DERIVATION
# =========================================================

if (!count_baseline_var %in% names(cohort)) {
  stop("condition_count_baseline not found. Add the precomputed baseline LTC count or adapt derivation section.")
}
if (!count_end_var %in% names(cohort)) {
  stop("condition_count_end not found. Add the precomputed follow-up LTC count or adapt derivation section.")
}

# end of observation date
if (!end_date_var %in% names(cohort)) {
  cohort[, study_end := as.Date("2020-12-31")]
  end_date_var <- "study_end"
}

cohort[, followup_end := get(end_date_var)]
cohort[!is.na(get(death_date_var)) & get(death_date_var) < followup_end, followup_end := get(death_date_var)]
cohort[, followup_years := calc_followup_years(get(index_date_var), followup_end)]
cohort <- cohort[followup_years > 0]

# accumulation outcome
cohort[, ltc_gain := pmax(get(count_end_var) - get(count_baseline_var), 0)]
cohort[, annualised_gain := ltc_gain / followup_years]
cohort[, gained_any_ltc := as.integer(ltc_gain >= 1)]
cohort[, gained_2plus_ltc := as.integer(ltc_gain >= 2)]

# time-to-next-LTC outcome if available
if (next_ltc_date_var %in% names(cohort)) {
  cohort[, event_next_ltc := as.integer(!is.na(get(next_ltc_date_var)) & get(next_ltc_date_var) <= followup_end)]
  cohort[, time_to_next_ltc_days := fifelse(
    event_next_ltc == 1,
    as.numeric(get(next_ltc_date_var) - get(index_date_var)),
    as.numeric(followup_end - get(index_date_var))
  )]
}

cohort[, progressed_to_5plus := as.integer(get(count_baseline_var) < 5 & get(count_end_var) >= 5)]
cohort[, progressed_to_10plus := as.integer(get(count_baseline_var) < 10 & get(count_end_var) >= 10)]

# =========================================================
# 8. CLEAN ANALYSIS DATASET
# =========================================================

analysis_vars <- unique(c(
  id_var, index_date_var, "followup_end", "followup_years",
  core_covariates,
  intersect(extra_covariates, names(cohort)),
  intersect(social_need_domains, names(cohort)),
  "social_need_count", "social_need_cat", "any_social_need",
  count_baseline_var, count_end_var,
  "ltc_gain", "annualised_gain", "gained_any_ltc", "gained_2plus_ltc",
  intersect(c("event_next_ltc", "time_to_next_ltc_days"), names(cohort)),
  "progressed_to_5plus", "progressed_to_10plus"
))

ana <- cohort[, ..analysis_vars]

# cleaning
ana[, sex := as.factor(sex)]
ana[, ethnicity := as.factor(ethnicity)]
ana[, imd_quintile := as.factor(imd_quintile)]
if ("region" %in% names(ana)) ana[, region := as.factor(region)]
ana[, social_need_cat := as.factor(social_need_cat)]

# complete-case dataset for main adjusted models
main_model_vars <- unique(c(
  "ltc_gain", "followup_years", count_baseline_var,
  core_covariates, "social_need_cat", "social_need_count", "any_social_need"
))
main_model_vars <- intersect(main_model_vars, names(ana))
ana_cc <- na.omit(ana[, ..main_model_vars])

# =========================================================
# 9. DESCRIPTIVE TABLES
# =========================================================

baseline_summary <- ana %>%
  as.data.frame() %>%
  summarise(
    N = n(),
    mean_age = mean(age, na.rm = TRUE),
    sd_age = sd(age, na.rm = TRUE),
    female_n = sum(sex == levels(sex)[grep("female", tolower(levels(sex)))[1]], na.rm = TRUE),
    median_baseline_conditions = median(.data[[count_baseline_var]], na.rm = TRUE),
    mean_followup_years = mean(followup_years, na.rm = TRUE),
    any_social_need_n = sum(any_social_need == 1, na.rm = TRUE),
    any_social_need_pct = mean(any_social_need == 1, na.rm = TRUE) * 100
  )
write.csv(baseline_summary, file.path(path_out, "table1_overall_summary.csv"), row.names = FALSE)

social_need_dist <- ana %>%
  as.data.frame() %>%
  count(social_need_cat) %>%
  mutate(percent = 100 * n / sum(n))
write.csv(social_need_dist, file.path(path_out, "table1_social_need_distribution.csv"), row.names = FALSE)

baseline_by_social_need <- ana %>%
  as.data.frame() %>%
  group_by(social_need_cat) %>%
  summarise(
    n = n(),
    mean_age = mean(age, na.rm = TRUE),
    female_pct = mean(sex == "Female", na.rm = TRUE) * 100,
    baseline_median = median(.data[[count_baseline_var]], na.rm = TRUE),
    followup_mean = mean(followup_years, na.rm = TRUE),
    mean_gain = mean(ltc_gain, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(baseline_by_social_need, file.path(path_out, "table1_baseline_by_social_need.csv"), row.names = FALSE)

# =========================================================
# 10. MAIN MODEL: LTC ACCUMULATION COUNT
# =========================================================
# Outcome: number of additional LTCs gained over follow-up
# Offset: log(follow-up years)
# Main exposure: social_need_cat and social_need_count

form_cat <- as.formula(paste0(
  "ltc_gain ~ social_need_cat + ", count_baseline_var, " + age + sex + ethnicity + imd_quintile",
  if ("region" %in% names(ana_cc)) " + region" else "",
  " + offset(log(followup_years))"
))

count_model_cat <- fit_nb_or_poisson(form_cat, ana_cc)
count_model_cat_tidy <- extract_rr(count_model_cat$model)
write.csv(count_model_cat_tidy, file.path(path_out, "model1_social_need_category_count_gain.csv"), row.names = FALSE)

form_count <- as.formula(paste0(
  "ltc_gain ~ social_need_count + ", count_baseline_var, " + age + sex + ethnicity + imd_quintile",
  if ("region" %in% names(ana_cc)) " + region" else "",
  " + offset(log(followup_years))"
))

count_model_count <- fit_nb_or_poisson(form_count, ana_cc)
count_model_count_tidy <- extract_rr(count_model_count$model)
write.csv(count_model_count_tidy, file.path(path_out, "model2_social_need_count_count_gain.csv"), row.names = FALSE)

# =========================================================
# 11. DOMAIN-SPECIFIC MODELS
# =========================================================

domain_terms <- intersect(social_need_domains, names(ana))
ana_domain <- na.omit(ana[, ..unique(c("ltc_gain", "followup_years", count_baseline_var, core_covariates, domain_terms))])

domain_formula <- as.formula(paste0(
  "ltc_gain ~ ", paste(domain_terms, collapse = " + "), " + ", count_baseline_var,
  " + age + sex + ethnicity + imd_quintile",
  if ("region" %in% names(ana_domain)) " + region" else "",
  " + offset(log(followup_years))"
))

domain_model <- fit_nb_or_poisson(domain_formula, ana_domain)
domain_model_tidy <- extract_rr(domain_model$model)
write.csv(domain_model_tidy, file.path(path_out, "model3_domain_specific_count_gain.csv"), row.names = FALSE)

# forest-plot ready table
forest_domains <- domain_model_tidy %>%
  filter(term %in% domain_terms)
write.csv(forest_domains, file.path(path_out, "figure_domain_forest_data.csv"), row.names = FALSE)

# =========================================================
# 12. TIME TO NEXT LTC (COX MODEL)
# =========================================================

if (all(c("event_next_ltc", "time_to_next_ltc_days") %in% names(ana))) {
  ana_surv <- na.omit(ana[, ..unique(c(
    "time_to_next_ltc_days", "event_next_ltc", "social_need_cat", "social_need_count",
    count_baseline_var, core_covariates, domain_terms
  ))])

  cox_cat <- coxph(
    as.formula(paste0(
      "Surv(time_to_next_ltc_days, event_next_ltc) ~ social_need_cat + ", count_baseline_var,
      " + age + sex + ethnicity + imd_quintile",
      if ("region" %in% names(ana_surv)) " + region" else ""
    )),
    data = ana_surv
  )
  write.csv(extract_hr(cox_cat), file.path(path_out, "cox_social_need_category_time_to_next_ltc.csv"), row.names = FALSE)

  cox_count <- coxph(
    as.formula(paste0(
      "Surv(time_to_next_ltc_days, event_next_ltc) ~ social_need_count + ", count_baseline_var,
      " + age + sex + ethnicity + imd_quintile",
      if ("region" %in% names(ana_surv)) " + region" else ""
    )),
    data = ana_surv
  )
  write.csv(extract_hr(cox_count), file.path(path_out, "cox_social_need_count_time_to_next_ltc.csv"), row.names = FALSE)

  # Kaplan-Meier style descriptive curves
  sf <- survfit(Surv(time_to_next_ltc_days, event_next_ltc) ~ social_need_cat, data = ana_surv)
  km_df <- broom::tidy(sf)
  write.csv(km_df, file.path(path_out, "km_social_need_category.csv"), row.names = FALSE)
}

# =========================================================
# 13. TRANSITION TO HIGHER COMPLEXITY STATES
# =========================================================

ana_prog5 <- na.omit(ana[, ..unique(c("progressed_to_5plus", "social_need_cat", "social_need_count", count_baseline_var, core_covariates))])
mod_prog5 <- glm(
  as.formula(paste0(
    "progressed_to_5plus ~ social_need_cat + ", count_baseline_var,
    " + age + sex + ethnicity + imd_quintile",
    if ("region" %in% names(ana_prog5)) " + region" else ""
  )),
  family = binomial(),
  data = ana_prog5
)
write.csv(broom::tidy(mod_prog5, conf.int = TRUE, exponentiate = TRUE), file.path(path_out, "logit_progression_to_5plus.csv"), row.names = FALSE)

ana_prog10 <- na.omit(ana[, ..unique(c("progressed_to_10plus", "social_need_cat", "social_need_count", count_baseline_var, core_covariates))])
mod_prog10 <- glm(
  as.formula(paste0(
    "progressed_to_10plus ~ social_need_cat + ", count_baseline_var,
    " + age + sex + ethnicity + imd_quintile",
    if ("region" %in% names(ana_prog10)) " + region" else ""
  )),
  family = binomial(),
  data = ana_prog10
)
write.csv(broom::tidy(mod_prog10, conf.int = TRUE, exponentiate = TRUE), file.path(path_out, "logit_progression_to_10plus.csv"), row.names = FALSE)

# =========================================================
# 14. INTERACTION ANALYSES
# =========================================================

interaction_specs <- list(
  age = "social_need_count * ns(age, df = 4)",
  sex = "social_need_count * sex",
  imd = "social_need_count * imd_quintile"
)

for (nm in names(interaction_specs)) {
  f <- as.formula(paste0(
    "ltc_gain ~ ", interaction_specs[[nm]], " + ", count_baseline_var,
    if (nm != "age") " + age" else "",
    if (nm != "sex") " + sex" else "",
    " + ethnicity + imd_quintile",
    if (nm != "imd") "" else "",
    if ("region" %in% names(ana_cc)) " + region" else "",
    " + offset(log(followup_years))"
  ))

  fit_i <- fit_nb_or_poisson(f, ana_cc)
  write.csv(
    extract_rr(fit_i$model),
    file.path(path_out, paste0("interaction_", nm, "_count_gain.csv")),
    row.names = FALSE
  )
}

# =========================================================
# 15. SENSITIVITY ANALYSES
# =========================================================

# Sensitivity 1: complete follow-up >= 1 year
ana_1y <- ana[followup_years >= 1]
ana_1y_cc <- na.omit(ana_1y[, ..main_model_vars])
fit_1y <- fit_nb_or_poisson(form_count, ana_1y_cc)
write.csv(extract_rr(fit_1y$model), file.path(path_out, "sensitivity_followup_ge_1y.csv"), row.names = FALSE)

# Sensitivity 2: baseline MLTC count 2-4 only
ana_24 <- ana[get(count_baseline_var) >= 2 & get(count_baseline_var) <= 4]
ana_24_cc <- na.omit(ana_24[, ..main_model_vars])
fit_24 <- fit_nb_or_poisson(form_count, ana_24_cc)
write.csv(extract_rr(fit_24$model), file.path(path_out, "sensitivity_baseline_2_to_4.csv"), row.names = FALSE)

# Sensitivity 3: binary outcome gained >=2 LTCs using robust Poisson
ana_bin <- na.omit(ana[, ..unique(c("gained_2plus_ltc", "social_need_count", count_baseline_var, core_covariates))])
rob_fit <- robust_poisson(
  as.formula(paste0(
    "gained_2plus_ltc ~ social_need_count + ", count_baseline_var,
    " + age + sex + ethnicity + imd_quintile",
    if ("region" %in% names(ana_bin)) " + region" else ""
  )),
  data = ana_bin
)
write.csv(
  broom::tidy(rob_fit$fit, conf.int = TRUE, exponentiate = TRUE),
  file.path(path_out, "sensitivity_robust_poisson_gained_2plus.csv"),
  row.names = FALSE
)

# =========================================================
# 16. LONG FORMAT OPTIONAL: REPEATED ANNUAL COUNTS / GEE
# =========================================================

path_long <- "data/df_mltc_long_annual.csv"
if (file.exists(path_long)) {
  long_df <- fread(path_long)
  setDT(long_df)

  if (all(c(id_var, "year_since_index", "condition_count", "social_need_count") %in% names(long_df))) {
    gee_df <- na.omit(long_df[, ..unique(c(
      id_var, "year_since_index", "condition_count", "social_need_count", core_covariates
    ))])

    gee_fit <- geepack::geeglm(
      as.formula(paste0(
        "condition_count ~ year_since_index * social_need_count + age + sex + ethnicity + imd_quintile",
        if ("region" %in% names(gee_df)) " + region" else ""
      )),
      id = gee_df[[id_var]],
      data = gee_df,
      family = poisson(link = "log"),
      corstr = "exchangeable"
    )

    write.csv(
      broom::tidy(gee_fit, conf.int = TRUE, exponentiate = TRUE),
      file.path(path_out, "gee_repeated_condition_counts.csv"),
      row.names = FALSE
    )
  }
}

# =========================================================
# 17. MACHINE LEARNING EXPLORATORY MODEL
# =========================================================
# Predict faster accumulation: top quartile annualised gain

ml_vars <- unique(c(
  "annualised_gain", core_covariates,
  intersect(extra_covariates, names(ana)),
  intersect(social_need_domains, names(ana)),
  "social_need_count", count_baseline_var
))
ml_df <- na.omit(ana[, ..ml_vars]) %>% as.data.frame()

if (nrow(ml_df) > 1000) {
  thr <- quantile(ml_df$annualised_gain, probs = 0.75, na.rm = TRUE)
  ml_df$fast_accumulation <- as.integer(ml_df$annualised_gain >= thr)

  cat_vars <- names(ml_df)[sapply(ml_df, function(x) is.factor(x) || is.character(x))]
  ml_enc <- fastDummies::dummy_cols(
    ml_df[, setdiff(names(ml_df), "annualised_gain"), drop = FALSE],
    select_columns = cat_vars,
    remove_first_dummy = FALSE,
    remove_selected_columns = TRUE
  )

  idx <- sample(seq_len(nrow(ml_enc)), size = floor(0.8 * nrow(ml_enc)))
  train_df <- ml_enc[idx, ]
  test_df  <- ml_enc[-idx, ]

  x_train <- as.matrix(train_df[, setdiff(names(train_df), "fast_accumulation")])
  y_train <- train_df$fast_accumulation
  x_test  <- as.matrix(test_df[, setdiff(names(test_df), "fast_accumulation")])
  y_test  <- test_df$fast_accumulation

  dtrain <- xgb.DMatrix(x_train, label = y_train)
  dtest  <- xgb.DMatrix(x_test, label = y_test)

  xgb_fit <- xgb.train(
    params = list(objective = "binary:logistic", eval_metric = "auc", eta = 0.1, max_depth = 4),
    data = dtrain,
    nrounds = 100,
    watchlist = list(train = dtrain, test = dtest),
    verbose = 0
  )

  pred <- predict(xgb_fit, dtest)
  perf <- data.frame(
    auc_proxy = NA_real_,
    accuracy = mean(ifelse(pred >= 0.5, 1, 0) == y_test)
  )
  write.csv(perf, file.path(path_out, "xgboost_fast_accumulation_performance.csv"), row.names = FALSE)

  shap_vals <- as.data.frame(predict(xgb_fit, dtest, predcontrib = TRUE))
  if ("BIAS" %in% names(shap_vals)) shap_vals <- shap_vals[, setdiff(names(shap_vals), "BIAS"), drop = FALSE]
  shap_imp <- data.frame(
    feature = names(shap_vals),
    mean_abs_shap = colMeans(abs(as.matrix(shap_vals)), na.rm = TRUE)
  ) %>% arrange(desc(mean_abs_shap))
  write.csv(shap_imp, file.path(path_out, "xgboost_fast_accumulation_shap.csv"), row.names = FALSE)
}

# =========================================================
# 18. FIGURES
# =========================================================

# Figure 1: mean annualised gain by social need category
fig1_df <- ana %>%
  as.data.frame() %>%
  group_by(social_need_cat) %>%
  summarise(
    mean_annualised_gain = mean(annualised_gain, na.rm = TRUE),
    se = sd(annualised_gain, na.rm = TRUE) / sqrt(sum(!is.na(annualised_gain))),
    .groups = "drop"
  )

p1 <- ggplot(fig1_df, aes(x = social_need_cat, y = mean_annualised_gain)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_annualised_gain - 1.96 * se, ymax = mean_annualised_gain + 1.96 * se), width = 0.15) +
  labs(x = "Social care need count category", y = "Mean annualised LTC gain", title = "Accumulation in LTC count by baseline social care need") +
  theme_minimal(base_size = 12)

ggsave(file.path(path_out, "figure1_mean_annualised_gain.png"), p1, width = 7, height = 5, dpi = 300)

# Figure 2: domain-specific IRRs
if (nrow(forest_domains) > 0) {
  p2 <- ggplot(forest_domains, aes(x = reorder(term, estimate), y = estimate)) +
    geom_point() +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.15) +
    geom_hline(yintercept = 1, linetype = 2) +
    coord_flip() +
    labs(x = NULL, y = "Incidence rate ratio", title = "Association of social care need domains with LTC accumulation") +
    theme_minimal(base_size = 12)

  ggsave(file.path(path_out, "figure2_domain_forest.png"), p2, width = 7, height = 5.5, dpi = 300)
}

# =========================================================
# 19. RESULT TABLES
# =========================================================

main_results <- bind_rows(
  count_model_cat_tidy %>% mutate(model = "Count gain ~ social need category"),
  count_model_count_tidy %>% mutate(model = "Count gain ~ social need count")
)
write.csv(main_results, file.path(path_out, "main_results_combined.csv"), row.names = FALSE)

# =========================================================
# 20. SESSION INFO
# =========================================================

capture.output(sessionInfo(), file = file.path(path_out, "sessionInfo.txt"))

message("social-care/MLTC accumulation analysis complete. Outputs saved to: ", path_out)
