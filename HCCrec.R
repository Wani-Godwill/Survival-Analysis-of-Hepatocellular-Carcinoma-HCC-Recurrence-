## --------------------------------------------------------------------------------------------------------------------------
# Installing the necessary packges
packages <- c(
  "readxl", "dplyr", "tidyr", "ggplot2", "survival", "survminer",
  "cmprsk", "gtsummary", "broom", "purrr", "stringr",
  "riskRegression", "pec", "rms", "dcurves"
)

to_install <- packages[!packages %in% installed.packages()[, "Package"]]
if(length(to_install) > 0) install.packages(to_install)



## --------------------------------------------------------------------------------------------------------------------------
# Import Libraries
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(survminer)
library(cmprsk)
library(gtsummary)
library(broom)
library(purrr)
library(stringr)
library(riskRegression)
library(pec)
library(rms)
library(dcurves)


## --------------------------------------------------------------------------------------------------------------------------
# Assign the arranged HCC recurrence dataset to the working object.
library(readxl)

hcc <- read_excel("HCCrec_arranged.xlsx")


## --------------------------------------------------------------------------------------------------------------------------
# Inspect
glimpse(hcc)
head(hcc)


## --------------------------------------------------------------------------------------------------------------------------
# Clean and record variables
hcc <- hcc %>%
  mutate(
    idpat = as.integer(idpat),
    Age = as.numeric(Age),
    Gender = factor(gsub('"', '', Gender), levels = c("F", "M")),
    
    RecMultinodular = factor(RecMultinodular,
                             levels = c(0, 1),
                             labels = c("Single", "Multinodular")),
    
    RecNoduleLargeSize = factor(RecNoduleLargeSize,
                                levels = c(0, 1),
                                labels = c("<=5cm", ">5cm")),
    
    RecExtrahepatic = factor(RecExtrahepatic,
                             levels = c(0, 1),
                             labels = c("Liver only", "Extrahepatic")),
    TimeToFirstRecMonths = as.numeric(TimeToFirstRecMonths),
    FupAfterFirstRecMonths = as.numeric(FupAfterFirstRecMonths),
    SecondRecOrDeath = as.integer(SecondRecOrDeath),
    RecTreat = factor(gsub('"', '', RecTreat), levels = c("CUR", "PAL")),
  )

hcc


## --------------------------------------------------------------------------------------------------------------------------
# Create Events Indicators
hcc <- hcc %>%
  mutate(
    event_secondrec = ifelse(SecondRecOrDeath == 1, 1, 0),
    event_death = ifelse(SecondRecOrDeath == 2, 1, 0),
    event_composite = ifelse(SecondRecOrDeath %in% c(1, 2), 1, 0)
  )
head(hcc)


## --------------------------------------------------------------------------------------------------------------------------
# Check Missing Values
colSums(is.na(hcc))


## --------------------------------------------------------------------------------------------------------------------------
# Check Event Distribution
table(hcc$SecondRecOrDeath)
prop.table(table(hcc$SecondRecOrDeath))


## --------------------------------------------------------------------------------------------------------------------------
# Descriptive Analyses
tbl_overall <- hcc %>%
  select(
    Age, Gender, RecMultinodular, RecNoduleLargeSize, RecExtrahepatic,
    TimeToFirstRecMonths, FupAfterFirstRecMonths, SecondRecOrDeath, RecTreat
  ) %>%
  tbl_summary(
    statistic = list(
      all_continuous() ~ "{mean} ({sd}) | median={median} [Q1={p25}, Q3={p75}]",
      all_categorical() ~ "{n} / {N} ({p}%)"
    ),
    digits = all_continuous() ~ 2
  ) %>%
  bold_labels()

tbl_overall


## --------------------------------------------------------------------------------------------------------------------------
# 1. Age distribution
ggplot(hcc, aes(x = Age)) +
  geom_histogram(bins = 20) +
  theme_minimal() +
  labs(title = "Distribution of Age", x = "Age", y = "Frequency")

# 2. Time to first recurrence distribution
ggplot(hcc, aes(x = TimeToFirstRecMonths)) +
  geom_histogram(bins = 20) +
  theme_minimal() +
  labs(title = "Distribution of Time to First Recurrence", x = "Months", y = "Frequency")

# 3. Follow-up time distribution
ggplot(hcc, aes(x = FupAfterFirstRecMonths)) +
  geom_histogram(bins = 20) +
  theme_minimal() +
  labs(title = "Distribution of Follow-up Time", x = "Months", y = "Frequency")

# 4. Treatment distribution
ggplot(hcc, aes(x = RecTreat)) +
  geom_bar() +
  theme_minimal() +
  labs(title = "Treatment Distribution", x = "Treatment", y = "Count")

# 6. Age by treatment
ggplot(hcc, aes(x = RecTreat, y = Age)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Age by Treatment", x = "Treatment", y = "Age")

# 7. Time to first recurrence by treatment
ggplot(hcc, aes(x = RecTreat, y = TimeToFirstRecMonths)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Time to First Recurrence by Treatment", x = "Treatment", y = "Months")



## --------------------------------------------------------------------------------------------------------------------------
# Estimate cumulative incidence curves by treatment group.
library(cmprsk)

ci_treat <- cuminc(
  ftime = hcc$FupAfterFirstRecMonths,
  fstatus = hcc$SecondRecOrDeath,
  group = hcc$RecTreat,
  cencode=0
)

plot(ci_treat, lty = 1, col = c("blue", "red", "black", "green"),
     xlab = "Months after first recurrence",
     ylab = "Cumulative incidence",
     main = "Aalen-Johansen cumulative incidence curves by treatment")



## --------------------------------------------------------------------------------------------------------------------------
# Show Gray's test results for differences between treatment groups.
ci_treat$Tests


## --------------------------------------------------------------------------------------------------------------------------
# Define covariates
covariates <- c(
  "Age", "Gender", "RecMultinodular", "RecNoduleLargeSize",
  "RecExtrahepatic", "TimeToFirstRecMonths", "RecTreat"
)


## --------------------------------------------------------------------------------------------------------------------------
# Univariate Cox for second recurrence
# Patients with death are censored.

run_uni_cox <- function(var, data, time, event) {
  fml <- as.formula(paste0("Surv(", time, ", ", event, ") ~ ", var))
  model <- coxph(fml, data = data)
  broom::tidy(model, exponentiate = TRUE, conf.int = TRUE)
}

uni_secondrec <- map_dfr(
  covariates,
  ~ run_uni_cox(.x, hcc, "FupAfterFirstRecMonths", "event_secondrec") %>%
    mutate(variable = .x, outcome = "Second recurrence")
)

uni_secondrec


## --------------------------------------------------------------------------------------------------------------------------
# Univariate Cox for death without second recurrence
# Patients with second recurrence are censored.

uni_death <- map_dfr(
  covariates,
  ~ run_uni_cox(.x, hcc, "FupAfterFirstRecMonths", "event_death") %>%
    mutate(variable = .x, outcome = "Death without second recurrence")
)

uni_death


## --------------------------------------------------------------------------------------------------------------------------
# Univariate Cox for composite endpoint
uni_composite <- map_dfr(
  covariates,
  ~ run_uni_cox(.x, hcc, "FupAfterFirstRecMonths", "event_composite") %>%
    mutate(variable = .x, outcome = "Composite endpoint")
)

uni_composite


## --------------------------------------------------------------------------------------------------------------------------
# Combine all results
uni_all <- bind_rows(uni_secondrec, uni_death, uni_composite) %>%
  select(outcome, variable, term, estimate, conf.low, conf.high, p.value) %>%
  rename(
    HR = estimate,
    CI_low = conf.low,
    CI_high = conf.high,
    p_value = p.value
  )

uni_all


## --------------------------------------------------------------------------------------------------------------------------
# Multivariable Cox model for the composite endpoint
cox_comp <- coxph(
  Surv(FupAfterFirstRecMonths, event_composite) ~
    Age + Gender + RecMultinodular + RecNoduleLargeSize +
    RecExtrahepatic + TimeToFirstRecMonths + RecTreat,
  data = hcc,
  x = TRUE,
  y = TRUE
)

summary(cox_comp)


## --------------------------------------------------------------------------------------------------------------------------
 #Present the multivariable Cox model as a formatted regression table.
tbl_regression(cox_comp, exponentiate = TRUE)


## --------------------------------------------------------------------------------------------------------------------------
####Set up rms environment
dd <- datadist(hcc)
options(datadist = "dd")


## --------------------------------------------------------------------------------------------------------------------------
cox_age_linear <- cph(
  Surv(FupAfterFirstRecMonths, event_composite) ~
    Age + Gender + RecMultinodular + RecNoduleLargeSize +
    RecExtrahepatic + TimeToFirstRecMonths + RecTreat,
  data = hcc, x = TRUE, y = TRUE, surv = TRUE
)

cox_age_spline <- cph(
  Surv(FupAfterFirstRecMonths, event_composite) ~
    rcs(Age, 4) + Gender + RecMultinodular + RecNoduleLargeSize +
    RecExtrahepatic + TimeToFirstRecMonths + RecTreat,
  data = hcc, x = TRUE, y = TRUE, surv = TRUE
)

anova(cox_age_linear)
anova(cox_age_spline)


## --------------------------------------------------------------------------------------------------------------------------
library(rms)
library(survival)

dd <- datadist(hcc)
options(datadist = "dd")

cox_age_spline <- cph(
  Surv(FupAfterFirstRecMonths, event_composite) ~
    rcs(Age, 4) + Gender + RecMultinodular + RecNoduleLargeSize +
    RecExtrahepatic + TimeToFirstRecMonths + RecTreat,
  data = hcc,
  x = TRUE, y = TRUE, surv = TRUE
)

anova(cox_age_spline)
# Plot the hazard ratio pattern for Age from the spline model.
plot(Predict(cox_age_spline, Age, fun = exp),
     xlab = "Age",
     ylab = "Hazard ratio",
     main = "Functional form of Age")


## --------------------------------------------------------------------------------------------------------------------------
cox_ttf_linear <- cph(
  Surv(FupAfterFirstRecMonths, event_composite) ~
    Age + Gender + RecMultinodular + RecNoduleLargeSize +
    RecExtrahepatic + TimeToFirstRecMonths + RecTreat,
  data = hcc, x = TRUE, y = TRUE, surv = TRUE
)

cox_ttf_spline <- cph(
  Surv(FupAfterFirstRecMonths, event_composite) ~
    Age + Gender + RecMultinodular + RecNoduleLargeSize +
    RecExtrahepatic + rcs(TimeToFirstRecMonths, 4) + RecTreat,
  data = hcc, x = TRUE, y = TRUE, surv = TRUE
)

anova(cox_ttf_linear)
anova(cox_ttf_spline)

plot(Predict(cox_ttf_spline, TimeToFirstRecMonths, fun = exp),
     xlab = "Time to first recurrence (months)",
     ylab = "Hazard ratio",
     main = "Functional form of TimeToFirstRecMonths")


## --------------------------------------------------------------------------------------------------------------------------
# Test the proportional hazards assumption for the Cox model.
ph_test <- cox.zph(cox_comp)
ph_test


## --------------------------------------------------------------------------------------------------------------------------
# Plot Schoenfeld residual diagnostics for proportional hazards.
plot(ph_test)


## --------------------------------------------------------------------------------------------------------------------------
# Evaluate 12-month discrimination and calibration for the Cox model.
score_obj <- Score(
  object = list("Cox model" = cox_comp),
  formula = Surv(FupAfterFirstRecMonths, event_composite) ~ 1,
  data = hcc,
  times = 12,
  metrics = c("auc", "brier"),
  plots = c("calibration", "roc"),
  conf.int = TRUE
)

score_obj


## --------------------------------------------------------------------------------------------------------------------------
# Extract the time-dependent AUC at 12 months.
score_obj$AUC$score


## --------------------------------------------------------------------------------------------------------------------------
# Extract the Brier score at 12 months.
score_obj$Brier$score


## --------------------------------------------------------------------------------------------------------------------------
# Plot calibration at 12 months for predicted versus observed risk.
plotCalibration(score_obj, times = 12)


## --------------------------------------------------------------------------------------------------------------------------
# Plot the ROC curve at 12 months.
plotROC(score_obj, times = 12)


## --------------------------------------------------------------------------------------------------------------------------
library(survival)
library(dcurves)
library(dplyr)
# Attempt decision curve analysis directly from the Cox model formula.
# 1. Fit Cox model
cox_comp <- coxph(
  Surv(FupAfterFirstRecMonths, event_composite) ~
    Age + Gender + RecMultinodular + RecNoduleLargeSize +
    RecExtrahepatic + TimeToFirstRecMonths + RecTreat,
  data = hcc,
  x = TRUE
)

# 2. Get baseline cumulative hazard
bh <- basehaz(cox_comp, centered = FALSE)

# cumulative baseline hazard at 12 months
H0_12 <- bh$hazard[max(which(bh$time <= 12))]

# 3. Linear predictor for each subject
lp <- predict(cox_comp, newdata = hcc, type = "lp")

# 4. Predicted event risk at 12 months
risk_12 <- 1 - exp(-H0_12 * exp(lp))

# 5. Create DCA dataset explicitly
hcc_dca <- hcc %>%
  mutate(risk_12 = risk_12) %>%
  select(FupAfterFirstRecMonths, event_composite, risk_12)

# Check
names(hcc_dca)
summary(hcc_dca$risk_12)
range(hcc_dca$risk_12, na.rm = TRUE)

# 6. Run DCA
dca_model <- dca(
  Surv(FupAfterFirstRecMonths, event_composite) ~ risk_12,
  data = hcc_dca,
  time = 12,
  thresholds = seq(0.01, 0.50, by = 0.01),
  label = list(risk_12 = "Cox model")
)

# 7. Plot
plot(dca_model, smooth = TRUE)


## --------------------------------------------------------------------------------------------------------------------------
# Select three random patients from the dataset for prediction.
new_patients <- hcc %>%
  slice_sample(n = 3) %>%
  select(
    idpat, Age, Gender, RecMultinodular, RecNoduleLargeSize,
    RecExtrahepatic, TimeToFirstRecMonths, RecTreat
  )

new_patients


## --------------------------------------------------------------------------------------------------------------------------
# Estimate 12-month survival for the selected patients.
sf <- survfit(cox_comp, newdata = new_patients)
summary(sf, times = 12)


## --------------------------------------------------------------------------------------------------------------------------

# Fit survival curves for the 3 selected patients
sf <- survfit(cox_comp, newdata = new_patients)

# Extract only the 12-month survival probabilities as a plain numeric vector
pred12 <- as.numeric(summary(sf, times = 12)$surv)

# Build a clean table: one row per patient
risk_12 <- data.frame(
  idpat = new_patients$idpat,
  surv_12 = pred12,
  risk_12 = 1 - pred12
)

print(risk_12)

