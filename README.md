
**# Association between Extreme Temperature Events and Hospitalisation/Mortality
**# in Individuals with Multimorbidity in England**
#**

# Main components:
# 1) derive multimorbidity onset (2nd LTC date)
# 2) restrict follow-up to age >=65 years
# 3) prepare weekly region-specific heatwave/cold-spell exposures from HadUK-Grid
# 4) create person-week time-to-event datasets for unplanned hospitalisation and mortality
# 5) fit Cox models with time-varying weekly exposures
# 6) fit DLNM-style exposure-lag models using weekly temperature summaries
# 7) run subgroup and sensitivity analyses
# 8) generate outputs for manuscript figures/tables

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(lubridate)
  library(survival)
  library(splines)
  library(dlnm)
  library(ggplot2)
  library(broom)
  library(purrr)
  library(forcats)
  library(zoo)
})

options(datatable.print.nrows = 100)
set.seed(42)

# =========================================================
# 1. PATHS
# =========================================================

path_cprd    <- "data/cprd_mltc_person_level.csv"
path_hes     <- "data/hes_admissions.csv"
path_ons     <- "data/ons_deaths.csv"
path_haduk   <- "data/haduk_grid_daily_region_temp.csv"
path_out     <- "outputs_extreme_temperature_mltc"
dir.create(path_out, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 2. USER-EDITABLE VARIABLE NAMES
# =========================================================

id_var              <- "patid"
practice_var        <- "practice_id"
region_var          <- "cprd_region"
reg_start_var       <- "reg_start"
reg_end_var         <- "reg_end"
dob_var             <- "dob"
sex_var             <- "sex"
ethnicity_var       <- "ethnicity"
imd_var             <- "imd_quintile"

# 54 LTC first-diagnosis date columns
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
  "chronic_pancreatitis_date", "infection_date"
)

# HES
hes_admidate_var    <- "admission_date"
hes_disdate_var     <- "discharge_date"
hes_admi_method2    <- "admi_method2"
hes_admi_method     <- "admi_method"

# ONS
ons_deathdate_var   <- "death_date"

# HadUK-Grid expected columns
haduk_region_var    <- "met_region"
haduk_date_var      <- "date"
haduk_tasmax_var    <- "tasmax"
haduk_tasmin_var    <- "tasmin"

# =========================================================
# 3. HELPERS
# =========================================================

assert_cols <- function(df, cols, nm = deparse(substitute(df))) {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) {
    stop(sprintf("%s missing columns: %s", nm, paste(miss, collapse = ", ")))
  }
}

as_date_safe <- function(x) {
  if (inherits(x, "Date")) return(x)
  suppressWarnings(as.Date(x))
}

calc_age_date <- function(dob, ref_date) {
  floor(as.numeric(ref_date - dob) / 365.25)
}

week_monday <- function(x) {
  x - lubridate::wday(x, week_start = 1) + 1
}

mode_value <- function(x) {
  ux <- unique(na.omit(x))
  if (length(ux) == 0) return(NA)
  ux[which.max(tabulate(match(x, ux)))]
}

safe_first <- function(x) {
  y <- x[!is.na(x)]
  if (length(y) == 0) NA else y[1]
}

make_region_map <- function() {
  tibble(
    cprd_region = c(
      "North East", "Yorkshire and The Humber",
      "North West",
      "West Midlands", "East Midlands",
      "East of England",
      "South Central", "South East Coast", "London",
      "South West",
      "Wales",
      "Scotland"
    ),
    met_region = c(
      "England E and NE", "England E and NE",
      "England NW and N Wales",
      "Midlands", "Midlands",
      "East Anglia",
      "England SE and Central South", "England SE and Central South", "England SE and Central South",
      "England SW and S Wales",
      "England SW and S Wales",
      "Scotland"
    )
  )
}

get_sorted_unique_dates <- function(row_dates) {
  x <- as.Date(unlist(row_dates), origin = "1970-01-01")
  x <- x[!is.na(x)]
  sort(unique(x))
}

extract_kth_date <- function(date_vec, k) {
  if (length(date_vec) < k) return(as.Date(NA))
  date_vec[k]
}

# =========================================================
# 4. LOAD DATA
# =========================================================

cprd <- fread(path_cprd)
hes  <- fread(path_hes)
ons  <- fread(path_ons)
haduk <- fread(path_haduk)

setDT(cprd); setDT(hes); setDT(ons); setDT(haduk)

assert_cols(cprd, c(id_var, practice_var, region_var, reg_start_var, reg_end_var,
                    dob_var, sex_var, ethnicity_var, imd_var))
assert_cols(hes, c(id_var, hes_admidate_var))
assert_cols(ons, c(id_var, ons_deathdate_var))
assert_cols(haduk, c(haduk_region_var, haduk_date_var, haduk_tasmax_var, haduk_tasmin_var))

for (cl in intersect(ltc_date_cols, names(cprd))) cprd[, (cl) := as_date_safe(get(cl))]
cprd[, (reg_start_var) := as_date_safe(get(reg_start_var))]
cprd[, (reg_end_var)   := as_date_safe(get(reg_end_var))]
cprd[, (dob_var)       := as_date_safe(get(dob_var))]
hes[,  (hes_admidate_var) := as_date_safe(get(hes_admidate_var))]
hes[,  (hes_disdate_var)  := as_date_safe(get(hes_disdate_var))]
ons[,  (ons_deathdate_var) := as_date_safe(get(ons_deathdate_var))]
haduk[, (haduk_date_var) := as_date_safe(get(haduk_date_var))]

ltc_date_cols <- intersect(ltc_date_cols, names(cprd))
if (length(ltc_date_cols) < 2) stop("Need at least two LTC date columns in CPRD input.")

# =========================================================
# 5. DERIVE MULTIMORBIDITY ONSET (SECOND LTC DATE)
# =========================================================

sorted_ltc_dates <- apply(cprd[, ..ltc_date_cols], 1, get_sorted_unique_dates)
cprd[, first_ltc_date  := as.Date(sapply(sorted_ltc_dates, extract_kth_date, k = 1), origin = "1970-01-01")]
cprd[, second_ltc_date := as.Date(sapply(sorted_ltc_dates, extract_kth_date, k = 2), origin = "1970-01-01")]
cprd[, ltc_count_total := rowSums(!is.na(.SD)), .SDcols = ltc_date_cols]

# body-system or type indicators 
cprd[, ltc_group_cat := cut(ltc_count_total,
                            breaks = c(1, 5, 10, Inf),
                            labels = c("2-5", "6-10", "11+"),
                            right = TRUE)]

# =========================================================
# 6. LINK REGION TO MET OFFICE REGION
# =========================================================

region_map <- make_region_map()
cprd <- merge(cprd, region_map, by.x = region_var, by.y = "cprd_region", all.x = TRUE)

# =========================================================
# 7. HES: DEFINE FIRST UNPLANNED ADMISSION
# =========================================================

hes[, unplanned_flag := fifelse(!is.na(get(hes_admi_method2)),
                                as.integer(get(hes_admi_method2) == 2),
                                as.integer(get(hes_admi_method) %in% c(2, 3, 4)))]

hes_unplanned <- hes[unplanned_flag == 1]
setorder(hes_unplanned, get(id_var), get(hes_admidate_var))
hes_first_unplanned <- hes_unplanned[, .(first_unplanned_admission = safe_first(get(hes_admidate_var))), by = id_var]

# =========================================================
# 8. ONS: DEATH DATE
# =========================================================

ons_first <- ons[, .(death_date = safe_first(get(ons_deathdate_var))), by = id_var]

# =========================================================
# 9. BUILD MAIN COHORT AND AGE>=65 LEFT TRUNCATION
# =========================================================

cohort <- merge(cprd, hes_first_unplanned, by = id_var, all.x = TRUE)
cohort <- merge(cohort, ons_first, by = id_var, all.x = TRUE)

cohort <- cohort[!is.na(second_ltc_date)]
cohort[, age_at_index := calc_age_date(get(dob_var), second_ltc_date)]
cohort <- cohort[age_at_index >= 18]

# age 65 date
cohort[, age65_date := get(dob_var) + years(65)]

study_start <- as.Date("2001-01-01")
study_end   <- as.Date("2019-12-29")

# analysis window start = max(index, study_start, age65_date, registration start)
cohort[, entry_date := pmax(second_ltc_date, study_start, age65_date, get(reg_start_var), na.rm = TRUE)]
# exit = min(reg end, study_end)
cohort[, exit_date := pmin(get(reg_end_var), study_end, na.rm = TRUE)]
cohort <- cohort[entry_date <= exit_date]

# linkage-eligible at index assumed already upstream; retain if linked outcome data structure present
cohort[, followup_days := as.numeric(exit_date - entry_date)]

# =========================================================
# 10. WEEKLY HADUK EXPOSURES
# =========================================================
# Thresholds pooled over 2001-2019 by region and ISO week.

haduk <- haduk[get(haduk_date_var) >= study_start & get(haduk_date_var) <= study_end]
haduk[, iso_week := isoweek(get(haduk_date_var))]
haduk[, year := year(get(haduk_date_var))]
haduk[, week_start := week_monday(get(haduk_date_var))]

thresholds <- haduk[
  , .(
    tasmax_p95 = quantile(get(haduk_tasmax_var), probs = 0.95, na.rm = TRUE),
    tasmin_p05 = quantile(get(haduk_tasmin_var), probs = 0.05, na.rm = TRUE)
  ),
  by = .(get(haduk_region_var), iso_week)
]
setnames(thresholds, "get", haduk_region_var)

haduk <- merge(haduk, thresholds, by = c(haduk_region_var, "iso_week"), all.x = TRUE)

# daily extreme indicators
haduk[, heat_day := as.integer(get(haduk_tasmax_var) >= tasmax_p95)]
haduk[, cold_day := as.integer(get(haduk_tasmin_var) <= tasmin_p05)]

# consecutive-day logic within region ordered by date
setorder(haduk, get(haduk_region_var), get(haduk_date_var))
haduk[, heat_run := sequence(rle(heat_day)$lengths), by = get(haduk_region_var)]
haduk[heat_day == 0, heat_run := 0]
haduk[, cold_run := sequence(rle(cold_day)$lengths), by = get(haduk_region_var)]
haduk[cold_day == 0, cold_run := 0]

haduk[, heatwave_day := as.integer(heat_day == 1 & heat_run >= 3)]
haduk[, coldspell_day := as.integer(cold_day == 1 & cold_run >= 2)]

# intensity = exceedance beyond threshold (degree-days)
haduk[, heat_intensity_dd := pmax(get(haduk_tasmax_var) - tasmax_p95, 0)]
haduk[, cold_intensity_dd := pmax(tasmin_p05 - get(haduk_tasmin_var), 0)]

# weekly aggregation Monday-start
weekly_exp <- haduk[
  , .(
    weekly_heatwave = as.integer(any(heatwave_day == 1, na.rm = TRUE)),
    weekly_coldspell = as.integer(any(coldspell_day == 1, na.rm = TRUE)),
    heat_intensity_week = sum(heat_intensity_dd, na.rm = TRUE),
    cold_intensity_week = sum(cold_intensity_dd, na.rm = TRUE),
    mean_tasmax = mean(get(haduk_tasmax_var), na.rm = TRUE),
    mean_tasmin = mean(get(haduk_tasmin_var), na.rm = TRUE),
    seasonal_avg_temp = mean((get(haduk_tasmax_var) + get(haduk_tasmin_var))/2, na.rm = TRUE)
  ),
  by = .(get(haduk_region_var), week_start)
]
setnames(weekly_exp, "get", haduk_region_var)

weekly_exp[, weekly_exposure_cat := fifelse(weekly_heatwave == 0 & weekly_coldspell == 0, 0L,
                                     fifelse(weekly_heatwave == 1 & weekly_coldspell == 0, 1L,
                                     fifelse(weekly_heatwave == 0 & weekly_coldspell == 1, 2L, 3L)))]

# lags 0-21 days approximated on weekly scale as lag0-lag3 weeks for main person-week dataset;
# daily DLNM section later uses continuous temp summaries.
setorder(weekly_exp, get(haduk_region_var), week_start)
for (lagk in 0:3) {
  weekly_exp[, paste0("heatwave_lagw", lagk) := dplyr::lag(weekly_heatwave, n = lagk, default = 0), by = get(haduk_region_var)]
  weekly_exp[, paste0("coldspell_lagw", lagk) := dplyr::lag(weekly_coldspell, n = lagk, default = 0), by = get(haduk_region_var)]
}

# =========================================================
# 11. PERSON-WEEK DATASET CREATION
# =========================================================

create_person_week_dataset <- function(df_cohort) {
  out_list <- vector("list", nrow(df_cohort))

  for (i in seq_len(nrow(df_cohort))) {
    rowi <- df_cohort[i]
    weeks <- seq(week_monday(rowi$entry_date), week_monday(rowi$exit_date), by = "week")
    tmp <- data.table(
      patid = rowi[[id_var]],
      week_start = weeks,
      tstart = as.numeric(weeks - rowi$entry_date),
      tstop = as.numeric(pmin(weeks + 6, rowi$exit_date) - rowi$entry_date) + 1,
      entry_date = rowi$entry_date,
      exit_date = rowi$exit_date,
      age65_date = rowi$age65_date,
      second_ltc_date = rowi$second_ltc_date,
      sex = rowi[[sex_var]],
      ethnicity = rowi[[ethnicity_var]],
      imd_quintile = rowi[[imd_var]],
      cprd_region = rowi[[region_var]],
      met_region = rowi$met_region,
      ltc_count_total = rowi$ltc_count_total,
      ltc_group_cat = rowi$ltc_group_cat,
      dob = rowi[[dob_var]],
      first_unplanned_admission = rowi$first_unplanned_admission,
      death_date = rowi$death_date
    )

    tmp[, age_this_week := calc_age_date(dob, week_start)]
    tmp[, age_band_tv := cut(age_this_week,
                             breaks = c(65, 75, 85, Inf),
                             labels = c("65-74", "75-84", "85+"),
                             right = FALSE)]
    out_list[[i]] <- tmp
  }

  rbindlist(out_list, fill = TRUE)
}

pw <- create_person_week_dataset(cohort)

# events: first unplanned hospital admission and death within that week
pw[, hosp_event := as.integer(!is.na(first_unplanned_admission) & first_unplanned_admission >= week_start & first_unplanned_admission <= (week_start + 6))]
pw[, death_event := as.integer(!is.na(death_date) & death_date >= week_start & death_date <= (week_start + 6))]

# keep person-time until first event for event-specific analyses
pw_hosp <- pw[, .SD[week_start <= ifelse(any(hosp_event == 1), min(week_start[hosp_event == 1]), max(week_start))], by = patid]
pw_death <- pw[, .SD[week_start <= ifelse(any(death_event == 1), min(week_start[death_event == 1]), max(week_start))], by = patid]

# merge weekly exposures
pw_hosp  <- merge(pw_hosp, weekly_exp, by.x = c("met_region", "week_start"), by.y = c(haduk_region_var, "week_start"), all.x = TRUE)
pw_death <- merge(pw_death, weekly_exp, by.x = c("met_region", "week_start"), by.y = c(haduk_region_var, "week_start"), all.x = TRUE)

# temporal covariates
for (obj in c("pw_hosp", "pw_death")) {
  X <- get(obj)
  X[, cal_year := year(week_start)]
  X[, month := month(week_start)]
  X[, season := factor(case_when(
    month %in% c(12, 1, 2) ~ "Winter",
    month %in% c(3, 4, 5) ~ "Spring",
    month %in% c(6, 7, 8) ~ "Summer",
    TRUE ~ "Autumn"
  ))]
  assign(obj, X)
}

# =========================================================
# 12. COX MODELS WITH TIME-VARYING EXPOSURES
# =========================================================

fit_cox_model <- function(df, event_var, exposure_var, sensitivity = FALSE) {
  fml <- as.formula(paste0(
    "Surv(tstart, tstop, ", event_var, ") ~ ", exposure_var,
    " + age_band_tv + sex + ethnicity + imd_quintile + ltc_group_cat + cprd_region",
    if (sensitivity) " + seasonal_avg_temp" else "",
    " + strata(cal_year, season)"
  ))
  coxph(fml, data = as.data.frame(df), ties = "efron")
}

cox_hosp_heat <- fit_cox_model(pw_hosp, "hosp_event", "weekly_heatwave")
cox_hosp_cold <- fit_cox_model(pw_hosp, "hosp_event", "weekly_coldspell")
cox_death_heat <- fit_cox_model(pw_death, "death_event", "weekly_heatwave")
cox_death_cold <- fit_cox_model(pw_death, "death_event", "weekly_coldspell")

write.csv(broom::tidy(cox_hosp_heat, exponentiate = TRUE, conf.int = TRUE), file.path(path_out, "cox_hosp_heat.csv"), row.names = FALSE)
write.csv(broom::tidy(cox_hosp_cold, exponentiate = TRUE, conf.int = TRUE), file.path(path_out, "cox_hosp_cold.csv"), row.names = FALSE)
write.csv(broom::tidy(cox_death_heat, exponentiate = TRUE, conf.int = TRUE), file.path(path_out, "cox_death_heat.csv"), row.names = FALSE)
write.csv(broom::tidy(cox_death_cold, exponentiate = TRUE, conf.int = TRUE), file.path(path_out, "cox_death_cold.csv"), row.names = FALSE)

# proportional hazards checks
ph_hosp_heat  <- cox.zph(cox_hosp_heat)
ph_hosp_cold  <- cox.zph(cox_hosp_cold)
ph_death_heat <- cox.zph(cox_death_heat)
ph_death_cold <- cox.zph(cox_death_cold)

capture.output(ph_hosp_heat,  file = file.path(path_out, "ph_hosp_heat.txt"))
capture.output(ph_hosp_cold,  file = file.path(path_out, "ph_hosp_cold.txt"))
capture.output(ph_death_heat, file = file.path(path_out, "ph_death_heat.txt"))
capture.output(ph_death_cold, file = file.path(path_out, "ph_death_cold.txt"))

# sensitivity adjustment for seasonal average temp
cox_hosp_heat_sens  <- fit_cox_model(pw_hosp, "hosp_event", "weekly_heatwave", sensitivity = TRUE)
cox_death_heat_sens <- fit_cox_model(pw_death, "death_event", "weekly_heatwave", sensitivity = TRUE)
write.csv(broom::tidy(cox_hosp_heat_sens, exponentiate = TRUE, conf.int = TRUE), file.path(path_out, "cox_hosp_heat_sensitivity.csv"), row.names = FALSE)
write.csv(broom::tidy(cox_death_heat_sens, exponentiate = TRUE, conf.int = TRUE), file.path(path_out, "cox_death_heat_sensitivity.csv"), row.names = FALSE)

# =========================================================
# 13. SUBGROUP / EFFECT MODIFICATION
# =========================================================

fit_subgroup_models <- function(df, event_var, exposure_var, subgroup_var, out_prefix) {
  levs <- unique(df[[subgroup_var]])
  levs <- levs[!is.na(levs)]
  out <- map_dfr(levs, function(g) {
    dsub <- df[get(subgroup_var) == g]
    if (sum(dsub[[event_var]], na.rm = TRUE) < 20) return(NULL)
    fit <- coxph(as.formula(paste0(
      "Surv(tstart, tstop, ", event_var, ") ~ ", exposure_var,
      " + age_band_tv + sex + ethnicity + imd_quintile + ltc_group_cat + cprd_region + strata(cal_year, season)"
    )), data = as.data.frame(dsub), ties = "efron")
    broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
      mutate(subgroup = subgroup_var, level = as.character(g))
  })
  write.csv(out, file.path(path_out, paste0(out_prefix, "_", subgroup_var, ".csv")), row.names = FALSE)
}

subgroups <- c("age_band_tv", "sex", "ethnicity", "imd_quintile", "cprd_region", "ltc_group_cat")
for (sg in subgroups) {
  fit_subgroup_models(pw_hosp,  "hosp_event",  "weekly_heatwave", sg, "subgroup_hosp_heat")
  fit_subgroup_models(pw_hosp,  "hosp_event",  "weekly_coldspell", sg, "subgroup_hosp_cold")
  fit_subgroup_models(pw_death, "death_event", "weekly_heatwave", sg, "subgroup_death_heat")
  fit_subgroup_models(pw_death, "death_event", "weekly_coldspell", sg, "subgroup_death_cold")
}

# interaction example
cox_hosp_heat_int_age <- coxph(
  Surv(tstart, tstop, hosp_event) ~ weekly_heatwave * age_band_tv + sex + ethnicity + imd_quintile + ltc_group_cat + cprd_region + strata(cal_year, season),
  data = as.data.frame(pw_hosp), ties = "efron"
)
write.csv(broom::tidy(cox_hosp_heat_int_age, exponentiate = TRUE, conf.int = TRUE), file.path(path_out, "interaction_hosp_heat_age.csv"), row.names = FALSE)

# =========================================================
# 14. DAILY DATA FOR DLNM-TYPE ANALYSIS
# =========================================================

create_person_day_sample <- function(df_cohort, max_n = 50000) {
  set.seed(42)
  samp <- df_cohort[sample(.N, min(.N, max_n))]
  out <- vector("list", nrow(samp))

  for (i in seq_len(nrow(samp))) {
    rowi <- samp[i]
    days <- seq(rowi$entry_date, rowi$exit_date, by = "day")
    tmp <- data.table(
      patid = rowi[[id_var]],
      date = days,
      tstart = as.numeric(days - rowi$entry_date),
      tstop = as.numeric(days - rowi$entry_date) + 1,
      met_region = rowi$met_region,
      sex = rowi[[sex_var]],
      ethnicity = rowi[[ethnicity_var]],
      imd_quintile = rowi[[imd_var]],
      cprd_region = rowi[[region_var]],
      ltc_group_cat = rowi$ltc_group_cat,
      dob = rowi[[dob_var]],
      age_band_tv = cut(calc_age_date(rowi[[dob_var]], days), breaks = c(65,75,85,Inf), labels = c("65-74","75-84","85+"), right = FALSE),
      hosp_event = as.integer(!is.na(rowi$first_unplanned_admission) & rowi$first_unplanned_admission == days),
      death_event = as.integer(!is.na(rowi$death_date) & rowi$death_date == days)
    )
    out[[i]] <- tmp
  }
  rbindlist(out, fill = TRUE)
}

pd_sample <- create_person_day_sample(cohort)

haduk_daily <- haduk[, .(met_region = get(haduk_region_var), date = get(haduk_date_var), tasmax = get(haduk_tasmax_var), tasmin = get(haduk_tasmin_var))]
pd_sample <- merge(pd_sample, haduk_daily, by = c("met_region", "date"), all.x = TRUE)
pd_sample[, cal_year := year(date)]
pd_sample[, month := month(date)]
pd_sample[, season := factor(case_when(
  month %in% c(12,1,2) ~ "Winter",
  month %in% c(3,4,5) ~ "Spring",
  month %in% c(6,7,8) ~ "Summer",
  TRUE ~ "Autumn"
))]

# crossbasis on tasmax and tasmin with lag 0-21 days
cb_heat <- crossbasis(pd_sample$tasmax, lag = 21,
                      argvar = list(fun = "ns", df = 4),
                      arglag = list(fun = "ns", df = 4))
cb_cold <- crossbasis(pd_sample$tasmin, lag = 21,
                      argvar = list(fun = "ns", df = 4),
                      arglag = list(fun = "ns", df = 4))

# DLNM-style Cox models on sample
cox_dlnm_hosp_heat <- coxph(
  Surv(tstart, tstop, hosp_event) ~ cb_heat + age_band_tv + sex + ethnicity + imd_quintile + cprd_region + ltc_group_cat + strata(cal_year, season),
  data = as.data.frame(pd_sample), ties = "efron"
)
cox_dlnm_death_heat <- coxph(
  Surv(tstart, tstop, death_event) ~ cb_heat + age_band_tv + sex + ethnicity + imd_quintile + cprd_region + ltc_group_cat + strata(cal_year, season),
  data = as.data.frame(pd_sample), ties = "efron"
)

capture.output(summary(cox_dlnm_hosp_heat), file = file.path(path_out, "cox_dlnm_hosp_heat.txt"))
capture.output(summary(cox_dlnm_death_heat), file = file.path(path_out, "cox_dlnm_death_heat.txt"))

# prediction surfaces for heat
pred_heat_hosp <- crosspred(cb_heat, cox_dlnm_hosp_heat, cen = median(pd_sample$tasmax, na.rm = TRUE), by = 1)
png(file.path(path_out, "dlnm_heat_hosp_surface.png"), width = 1200, height = 900, res = 150)
plot(pred_heat_hosp, xlab = "Tasmax", zlab = "HR", ylab = "Lag days", theta = 210, phi = 35, ltheta = 170)
dev.off()

pred_heat_death <- crosspred(cb_heat, cox_dlnm_death_heat, cen = median(pd_sample$tasmax, na.rm = TRUE), by = 1)
png(file.path(path_out, "dlnm_heat_death_surface.png"), width = 1200, height = 900, res = 150)
plot(pred_heat_death, xlab = "Tasmax", zlab = "HR", ylab = "Lag days", theta = 210, phi = 35, ltheta = 170)
dev.off()

# =========================================================
# 15. DESCRIPTIVE TABLES
# =========================================================

baseline_table <- cohort %>%
  as.data.frame() %>%
  summarise(
    N = n(),
    mean_age_at_entry = mean(calc_age_date(dob, entry_date), na.rm = TRUE),
    sd_age_at_entry = sd(calc_age_date(dob, entry_date), na.rm = TRUE),
    female_pct = mean(sex == "Female", na.rm = TRUE) * 100,
    mean_ltc = mean(ltc_count_total, na.rm = TRUE),
    median_ltc = median(ltc_count_total, na.rm = TRUE),
    deaths = sum(!is.na(death_date) & death_date >= entry_date & death_date <= exit_date, na.rm = TRUE),
    first_unplanned = sum(!is.na(first_unplanned_admission) & first_unplanned_admission >= entry_date & first_unplanned_admission <= exit_date, na.rm = TRUE)
  )
write.csv(baseline_table, file.path(path_out, "table1_baseline.csv"), row.names = FALSE)

baseline_by_region <- cohort %>%
  as.data.frame() %>%
  count(cprd_region, ltc_group_cat)
write.csv(baseline_by_region, file.path(path_out, "baseline_by_region_ltcgroup.csv"), row.names = FALSE)

exposure_counts <- weekly_exp %>%
  as.data.frame() %>%
  group_by(met_region, year = year(week_start)) %>%
  summarise(
    n_heatwave_weeks = sum(weekly_heatwave, na.rm = TRUE),
    n_coldspell_weeks = sum(weekly_coldspell, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(exposure_counts, file.path(path_out, "supp_exposure_counts_by_region_year.csv"), row.names = FALSE)

# =========================================================
# 16. FIGURES
# =========================================================

# Forest plot data
forest_main <- bind_rows(
  broom::tidy(cox_hosp_heat, exponentiate = TRUE, conf.int = TRUE) %>% mutate(model = "Hospitalisation", exposure = "Heatwave"),
  broom::tidy(cox_hosp_cold, exponentiate = TRUE, conf.int = TRUE) %>% mutate(model = "Hospitalisation", exposure = "Cold spell"),
  broom::tidy(cox_death_heat, exponentiate = TRUE, conf.int = TRUE) %>% mutate(model = "Mortality", exposure = "Heatwave"),
  broom::tidy(cox_death_cold, exponentiate = TRUE, conf.int = TRUE) %>% mutate(model = "Mortality", exposure = "Cold spell")
) %>% filter(term %in% c("weekly_heatwave", "weekly_coldspell"))
write.csv(forest_main, file.path(path_out, "forest_main_results.csv"), row.names = FALSE)

p_forest <- ggplot(forest_main, aes(x = exposure, y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_pointrange() +
  geom_hline(yintercept = 1, linetype = 2) +
  facet_wrap(~model) +
  coord_flip() +
  labs(x = NULL, y = "Hazard ratio", title = "Extreme temperature events and outcomes") +
  theme_minimal(base_size = 12)

ggsave(file.path(path_out, "figure_main_forest.png"), p_forest, width = 8, height = 5, dpi = 300)

# KM curves by weekly exposure at entry week (descriptive only)
cohort_entry <- pw_death %>%
  as.data.frame() %>%
  group_by(patid) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(entry_heat = factor(weekly_heatwave, labels = c("No heatwave", "Heatwave")))

sf <- survfit(Surv(tstop, death_event) ~ entry_heat, data = cohort_entry)
png(file.path(path_out, "km_death_by_entry_heat.png"), width = 1200, height = 900, res = 150)
plot(sf, xlab = "Days since entry", ylab = "Survival probability", col = 1:2)
legend("bottomleft", legend = levels(cohort_entry$entry_heat), col = 1:2, lty = 1)
dev.off()

# Regional heatmap of exposure weeks
heatmap_df <- exposure_counts %>%
  group_by(met_region) %>%
  summarise(total_heatwave_weeks = sum(n_heatwave_weeks), total_coldspell_weeks = sum(n_coldspell_weeks), .groups = "drop")

p_heatmap <- ggplot(heatmap_df, aes(x = "Heatwave weeks", y = fct_reorder(met_region, total_heatwave_weeks), fill = total_heatwave_weeks)) +
  geom_tile() +
  labs(x = NULL, y = NULL, title = "Regional heatwave exposure weeks") +
  theme_minimal(base_size = 12)

ggsave(file.path(path_out, "regional_heatmap_heatwave.png"), p_heatmap, width = 6, height = 4, dpi = 300)

# =========================================================
# 17. FLOWCHART COUNTS
# =========================================================

flowchart_counts <- tibble(
  step = c(
    "Initial CPRD input",
    "With >=2 LTCs",
    "Age >=18 at MM onset",
    "Eligible at age >=65 during study period",
    "Final cohort"
  ),
  n = c(
    nrow(cprd),
    sum(!is.na(cprd$second_ltc_date), na.rm = TRUE),
    sum(cprd$age_at_index >= 18, na.rm = TRUE),
    sum(cohort$entry_date <= cohort$exit_date, na.rm = TRUE),
    nrow(cohort)
  )
)
write.csv(flowchart_counts, file.path(path_out, "figure1_flowchart_counts.csv"), row.names = FALSE)

# =========================================================
# 18. SESSION INFO
# =========================================================

capture.output(sessionInfo(), file = file.path(path_out, "sessionInfo.txt"))
message("Extreme temperature / multimorbidity reconstructed analysis complete.")
