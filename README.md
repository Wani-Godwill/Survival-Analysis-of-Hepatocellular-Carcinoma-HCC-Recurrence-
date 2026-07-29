# Big Data in Health Care: HCC Recurrence Survival Analysis

## Survival Analysis of Hepatocellular Carcinoma Recurrence

![R](https://img.shields.io/badge/R-Survival%20Analysis-blue)
![Healthcare](https://img.shields.io/badge/Domain-Healthcare-green)
![Cox Model](https://img.shields.io/badge/Model-Cox%20Regression-orange)
![License](https://img.shields.io/badge/License-Academic-lightgrey)

---

## Project Overview

This project investigates clinical outcomes after the first recurrence of **hepatocellular carcinoma (HCC)** using survival analysis methods.

The analysis focuses on three outcomes:

1. Second recurrence
2. Death without second recurrence
3. A composite endpoint combining second recurrence and death

The project applies both competing-risk and Cox proportional hazards methods to identify prognostic factors, assess model assumptions, evaluate prediction performance at 12 months, and estimate individual patient risk.

---

## Clinical Background

Hepatocellular carcinoma is a major form of primary liver cancer. Even after curative treatment, many patients experience recurrence.

After the first recurrence, patients may experience:

- A second recurrence
- Death before a second recurrence
- Censoring if neither event is observed during follow-up

Because second recurrence and death are competing clinical outcomes, this project uses competing-risk methods alongside standard survival analysis.

---

## Dataset

The dataset contains **317 patient records** with no missing values.

The main variables are:

| Variable | Description |
|---|---|
| `idpat` | Patient identifier |
| `Age` | Age at first recurrence |
| `Gender` | Female or male |
| `RecMultinodular` | Single or multinodular recurrence |
| `RecNoduleLargeSize` | Nodule size ≤5 cm or >5 cm |
| `RecExtrahepatic` | Liver-only or extrahepatic recurrence |
| `TimeToFirstRecMonths` | Time from primary treatment to first recurrence |
| `FupAfterFirstRecMonths` | Follow-up time after first recurrence |
| `SecondRecOrDeath` | 0 = censored, 1 = second recurrence, 2 = death |
| `RecTreat` | Curative or palliative treatment |

Three event indicators are derived:

- `event_secondrec`
- `event_death`
- `event_composite`

The R script explicitly recodes the variables and creates these three event indicators before analysis. :contentReference[oaicite:1]{index=1}

---

## Cohort Summary

The cohort includes:

- 317 patients
- Mean age: 68.88 years
- Median age: 71 years
- 77% male
- 56% with multinodular recurrence
- 14% with nodules larger than 5 cm
- 19% with extrahepatic spread
- 63% receiving palliative treatment
- 103 second recurrences
- 30 deaths without second recurrence
- 184 censored observations

The median time to first recurrence is 17.28 months, while the median follow-up after recurrence is 12 months. :contentReference[oaicite:2]{index=2} :contentReference[oaicite:3]{index=3}

---

## Project Objectives

The main objectives are to:

- Describe the clinical characteristics of the cohort
- Estimate cumulative incidence for competing outcomes
- Compare curative and palliative treatment groups
- Identify prognostic factors using Cox regression
- Fit a multivariable prediction model
- Assess linearity of continuous predictors
- Test the proportional hazards assumption
- Evaluate discrimination and calibration at 12 months
- Perform decision-curve analysis
- Generate patient-specific survival predictions

---



---

## Statistical Methods

### 1. Descriptive Analysis

The analysis includes:

- Summary statistics
- Histograms
- Treatment distribution
- Boxplots by treatment group
- Event-frequency tables

The exploratory analysis shows that age is concentrated among older adults, while time to first recurrence and follow-up duration are right-skewed. Palliative treatment is more frequent than curative treatment. :contentReference[oaicite:4]{index=4}

---

### 2. Competing-Risk Analysis

The project uses the **Aalen–Johansen estimator** to estimate cumulative incidence curves for:

- Second recurrence
- Death without second recurrence

The curves are stratified by treatment group.

Gray’s test is used to compare cumulative incidence between curative and palliative treatment groups.

Results:

- Second recurrence: statistically significant difference between treatment groups
- Death without second recurrence: no statistically significant difference

Gray’s test gives:

- Second recurrence: p = 0.0013
- Death without recurrence: p = 0.637

:contentReference[oaicite:5]{index=5} :contentReference[oaicite:6]{index=6}

---

### 3. Univariate Cox Regression

Separate univariate Cox models are fitted for:

- Second recurrence
- Death without second recurrence
- Composite endpoint

The covariates include:

- Age
- Gender
- Multinodular recurrence
- Large nodule size
- Extrahepatic spread
- Time to first recurrence
- Treatment type

For second recurrence, the main adverse factors are:

- Multinodular recurrence
- Nodule size greater than 5 cm
- Palliative treatment

A longer time to first recurrence is associated with a lower hazard. :contentReference[oaicite:7]{index=7}

---

### 4. Multivariable Cox Regression

A multivariable Cox proportional hazards model is fitted for the composite endpoint.

The model includes:

```r
Age +
Gender +
RecMultinodular +
RecNoduleLargeSize +
RecExtrahepatic +
TimeToFirstRecMonths +
RecTreat
```

The R implementation fits the model using `coxph()` and reports exponentiated coefficients as hazard ratios. :contentReference[oaicite:8]{index=8}

---

## Main Multivariable Results

Significant adverse prognostic factors are:

| Predictor | Adjusted HR | Interpretation |
|---|---:|---|
| Multinodular recurrence | 1.59 | Higher hazard than single recurrence |
| Nodule size >5 cm | 1.86 | Higher hazard than nodules ≤5 cm |
| Palliative treatment | 1.56 | Higher hazard than curative treatment |

A longer time to first recurrence is protective:

| Predictor | Adjusted HR |
|---|---:|
| Time to first recurrence | 0.98 per month |

Age, gender, and extrahepatic recurrence are not statistically significant after adjustment.

The model has moderate discrimination with a C-index of approximately 0.67. :contentReference[oaicite:9]{index=9}

---

## Functional-Form Assessment

Restricted cubic splines with four knots are used to assess whether the continuous predictors should be modelled linearly.

The variables evaluated are:

- Age
- Time to first recurrence

The analysis finds no strong evidence of non-linearity.

The nonlinear component for time to first recurrence is not significant:

- Chi-square = 0.82
- p = 0.66

Therefore, both continuous predictors are retained in approximately linear form. :contentReference[oaicite:10]{index=10}

---

## Proportional Hazards Assumption

The proportional hazards assumption is assessed using:

- Schoenfeld residual tests
- Schoenfeld residual plots
- Global proportional hazards test

The global test gives:

- p = 0.280

All individual covariates also have p-values above 0.05.

Therefore, there is no strong evidence that the proportional hazards assumption is violated. :contentReference[oaicite:11]{index=11}

---

## Model Performance at 12 Months

The Cox model is evaluated at a 12-month prediction horizon using:

- Time-dependent AUC
- Brier score
- Calibration curve
- ROC curve

The analysis reports:

- C-index: approximately 0.67
- 12-month AUC: approximately 0.69

This indicates moderate discrimination.

The R script evaluates the model using the `Score()` function from the `riskRegression` package and produces calibration and ROC plots at 12 months. :contentReference[oaicite:12]{index=12}

---

## Decision-Curve Analysis

Decision-curve analysis is used to assess whether the Cox model provides clinical value compared with:

- Treating all patients as high risk
- Treating no patients as high risk

The model’s predicted 12-month risk is calculated from the baseline cumulative hazard and individual linear predictors.

The analysis evaluates threshold probabilities between 1% and 50%.

The report indicates that the model offers clinical utility mainly at moderate risk thresholds, approximately 25% to 40%. :contentReference[oaicite:13]{index=13}

---

## Individual Risk Prediction

The final model is used to estimate 12-month survival probabilities for selected patients.

The process includes:

1. Selecting patient profiles
2. Applying the fitted Cox model
3. Estimating patient-specific survival curves
4. Extracting survival probability at 12 months
5. Converting survival probability into event risk when required

The R script uses `survfit()` to generate these patient-specific predictions. :contentReference[oaicite:14]{index=14}

---

## Key Findings

- Multinodular recurrence is associated with increased post-recurrence risk.
- Nodules larger than 5 cm are an important adverse prognostic factor.
- Palliative treatment is associated with a higher hazard than curative treatment.
- A longer interval before first recurrence is associated with better prognosis.
- Age, gender, and extrahepatic spread are not independently significant in the adjusted model.
- The Cox model satisfies the proportional hazards assumption.
- Predictive performance is moderate rather than strong.
- The model may have clinical value at selected risk thresholds.

---

## Technologies and R Packages

The project is implemented in **R**.

Main packages include:

- `readxl`
- `dplyr`
- `tidyr`
- `ggplot2`
- `survival`
- `survminer`
- `cmprsk`
- `gtsummary`
- `broom`
- `purrr`
- `riskRegression`
- `pec`
- `rms`
- `dcurves`

These packages are installed and loaded directly in the R script. :contentReference[oaicite:15]{index=15}

---



---

## How to Run

### 1. Clone the repository

```bash
git clone https://github.com/your-username/HCC-Recurrence-Survival-Analysis.git
cd HCC-Recurrence-Survival-Analysis
```

### 2. Place the dataset in the project folder

The R script expects:

```text
HCCrec_arranged.xlsx
```

### 3. Open the R script

Open:

```text
HCCrec.R
```

using RStudio or another R environment.

### 4. Run the script

The script will:

- Install missing packages
- Load the dataset
- Clean and recode the variables
- Create event indicators
- Perform descriptive and survival analyses
- Fit Cox models
- Generate diagnostic plots
- Evaluate the model at 12 months
- Produce individual predictions

---


---

## Limitations

- The cohort contains only 317 patients.
- The number of death events is relatively small.
- The analysis is based on a single dataset.
- Internal performance assessment does not replace external validation.
- Treatment effects should not be interpreted causally because treatment was not randomly assigned.
- The model shows moderate predictive performance and is not ready for direct clinical deployment.

---

## Future Improvements

Possible extensions include:

- External validation on an independent HCC cohort
- Fine–Gray regression for cause-specific competing-risk prediction
- Bootstrap validation and optimism correction
- Penalized Cox regression
- Time-varying covariates
- Dynamic prediction using longitudinal biomarkers
- Comparison with machine learning survival models
- Development of an interactive clinical risk calculator

---

## Disclaimer

This project was developed for academic and research purposes.

The model is not intended for direct clinical diagnosis, treatment selection, or patient management without further validation.

---

## License

This repository is intended for academic and educational use.
