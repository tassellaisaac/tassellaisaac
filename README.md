# --------------------------------------------------
# SOCIAL NEED FRAMEWORK - EHR MAPPING ANALYSIS
# --------------------------------------------------

library(data.table)
library(dplyr)
library(stringr)
library(broom)

# --------------------------------------------------
# LOAD DATA
# --------------------------------------------------
df <- fread("cprd_social_need_data.csv")

# --------------------------------------------------
# DEFINE SOCIAL NEED DOMAINS 
# --------------------------------------------------

adl_codes <- c("adl1","adl2")
mobility_codes <- c("mob1","mob2")
financial_codes <- c("fin1","fin2")
disability_codes <- c("dis1","dis2")
community_codes <- c("com1","com2")
housing_codes <- c("house1","house2")
social_codes <- c("soc1","soc2")
bereavement_codes <- c("bereav1","bereav2")

# --------------------------------------------------
# MAP DOMAINS (binary flags)
# --------------------------------------------------

df <- df %>%
  mutate(
    adl = as.integer(code %in% adl_codes),
    mobility = as.integer(code %in% mobility_codes),
    financial = as.integer(code %in% financial_codes),
    disability = as.integer(code %in% disability_codes),
    community = as.integer(code %in% community_codes),
    housing = as.integer(code %in% housing_codes),
    social = as.integer(code %in% social_codes),
    bereavement = as.integer(code %in% bereavement_codes)
  )

# --------------------------------------------------
# COLLAPSE TO PATIENT LEVEL
# --------------------------------------------------

df_patient <- df %>%
  group_by(patient_id) %>%
  summarise(
    across(adl:bereavement, max, na.rm = TRUE),
    age = first(age),
    sex = first(sex),
    ethnicity = first(ethnicity),
    imd = first(imd),
    region = first(region),
    ltc_count = first(ltc_count)
  ) %>%
  ungroup()

# --------------------------------------------------
# DERIVE SOCIAL NEED MEASURES
# --------------------------------------------------

df_patient <- df_patient %>%
  mutate(
    any_social_need = as.integer(rowSums(across(adl:bereavement)) > 0),
    social_need_count = rowSums(across(adl:bereavement))
  )

# --------------------------------------------------
# DESCRIPTIVE STATISTICS
# --------------------------------------------------

summary_table <- df_patient %>%
  summarise(
    total_n = n(),
    any_need = mean(any_social_need),
    mean_social_need = mean(social_need_count),
    mean_ltc = mean(ltc_count)
  )

print(summary_table)

# --------------------------------------------------
# LOGISTIC REGRESSION: ANY SOCIAL NEED
# --------------------------------------------------

logit_model <- glm(
  any_social_need ~ ltc_count + age + sex + ethnicity + imd + region,
  data = df_patient,
  family = binomial()
)

summary(logit_model)
exp(cbind(OR = coef(logit_model), confint(logit_model)))

# --------------------------------------------------
# DOMAIN-SPECIFIC MODELS
# --------------------------------------------------

domains <- c("adl","mobility","financial","disability",
             "community","housing","social","bereavement")

results <- lapply(domains, function(dom){
  model <- glm(
    as.formula(paste(dom, "~ ltc_count + age + sex + ethnicity + imd + region")),
    data = df_patient,
    family = binomial()
  )
  tidy(model, exponentiate = TRUE, conf.int = TRUE)
})

names(results) <- domains

# --------------------------------------------------
# LINEAR MODEL: SOCIAL NEED COUNT
# --------------------------------------------------

lm_model <- lm(
  social_need_count ~ ltc_count + age + sex + ethnicity + imd + region,
  data = df_patient
)

summary(lm_model)

# --------------------------------------------------
# VISUALISATION
# --------------------------------------------------

library(ggplot2)

ggplot(df_patient, aes(x = ltc_count, y = social_need_count)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "lm", color = "blue") +
  labs(title = "LTC vs Social Need Burden")

# --------------------------------------------------
# STRATIFIED ANALYSIS
# --------------------------------------------------

age_groups <- df_patient %>%
  mutate(age_group = cut(age, breaks = c(18,40,60,80,100))) %>%
  group_by(age_group) %>%
  summarise(mean_need = mean(any_social_need))

print(age_groups)

# --------------------------------------------------
# END
# --------------------------------------------------
