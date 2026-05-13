# Survival Prediction Models for Liver Disease Patients

This project compares several survival prediction models for patients with liver disease using clinical, laboratory, bile acid, and Raman spectroscopy predictors.

The goal is to evaluate whether machine learning and hybrid survival modelling methods can improve prediction performance compared to traditional Cox-based models.

This project was developed as part of my Bachelor thesis at TU Dortmund University.

---

## Project Overview

Survival prediction is important in clinical decision-making because it helps estimate the risk of an event over time while accounting for censored observations.

In this project, the event of interest is **all-cause mortality**, with transplant cases treated as right-censored observations. The analysis compares different modelling strategies using time-to-event data from liver disease patients.

The full thesis compares 15 predefined predictor combinations across four predictor domains: laboratory measurements, bile acid profiles, Raman spectroscopy features, and clinical characteristics.

---

## Objectives

The main objectives of this project are:

- Compare traditional and machine learning survival prediction models
- Evaluate the contribution of different predictor groups
- Assess model performance using nested cross-validation
- Compare predictive accuracy using Brier Score and Index of Prediction Accuracy (IPA)
- Identify which modelling approach performs best across different predictor combinations

---

## Methods

The following survival modelling approaches were compared:

### 1. LASSO-Cox Model
A Cox Proportional Hazards model with LASSO regularization for variable selection and coefficient shrinkage.

### 2. Random Survival Forests
A tree-based ensemble method for right-censored survival data that can capture nonlinear effects and interactions.

### 3. Tree-Cluster Cox Model
A hybrid approach where survival trees are used to create patient clusters, followed by Cox models fitted within each cluster.

### 4. Hybrid RSF-Cox Model
Random Survival Forests are used for variable ranking and selection, followed by Cox modelling on the selected predictors.

---

## Dataset

The original dataset contains clinical information from patients with liver disease, including:

- Survival time and event status
- Clinical characteristics
- Laboratory measurements
- Bile acid profile variables
- Raman spectroscopy features

Due to privacy and data protection reasons, the raw patient-level dataset is **not included** in this repository.

Instead, this repository focuses on:

- The report
- Analysis workflow
- Model comparison methodology
- Selected figures and result summaries

---

## Evaluation Metrics

Model performance was evaluated using:

### Brier Score
Measures prediction error for survival probability estimates. Lower values indicate better predictive performance.

### Index of Prediction Accuracy (IPA)
A rescaled version of the Brier Score relative to a null Kaplan-Meier model. Higher values indicate better predictive performance.

### Nested Cross-Validation
A nested cross-validation framework was used to reduce bias during model tuning and performance evaluation.

---

## Key Findings

The main findings from the thesis are:

- Random Survival Forests showed the strongest overall performance.
- RSF achieved the best mean IPA in 9 out of 15 predictor combinations.
- The Hybrid RSF-Cox method achieved the single highest mean IPA for the laboratory + bile acid + Raman predictor combination.
- LASSO-Cox was generally weaker, especially for high-dimensional predictor groups.
- Laboratory and clinical variables were the most consistently useful predictor domains.
- Raman features showed potential when combined with other predictor groups, but were less stable when used alone.

---

## Project Structure

```text
survival-prediction-liver-disease/
│
├── README.md
├── LICENSE
├── .gitignore
│
├── report/
│   └── survival_prediction_report.pdf
│
├── scripts/
│   └── analysis_pipeline.R
│
├── figures/
│   ├── km_curve.png
│   ├── model_comparison.png
│   └── variable_importance.png
│
└── results/
    └── summary_results.csv
```

---

## Tools and Technologies

The analysis was conducted in R.

Main R packages used:

- `survival`
- `glmnet`
- `randomForestSRC`
- `rpart`
- `riskRegression`
- `pec`
- `tidyverse`
- `openxlsx`

---

## How to Run

> Note: The raw clinical dataset is not included in this repository due to privacy restrictions.

If access to a compatible dataset is available, the general workflow is:

```r
# Install required packages
install.packages(c(
  "survival",
  "glmnet",
  "randomForestSRC",
  "rpart",
  "riskRegression",
  "pec",
  "tidyverse",
  "openxlsx"
))
```

Then run the analysis script:

```r
source("scripts/LiverAnalysis.R")
```

The expected workflow includes:

1. Load and preprocess survival data
2. Define predictor groups
3. Perform missing value handling
4. Train survival models
5. Run nested cross-validation
6. Compute Brier Score and IPA
7. Export model comparison results

---

## Report

The full thesis report is available in the `report/` folder:

```text
report/survival_prediction_report.pdf
```

---

## Skills Demonstrated

This project demonstrates experience in:

- Survival analysis
- Statistical modelling
- Machine learning for healthcare data
- Cross-validation and model evaluation
- High-dimensional predictor analysis
- R programming
- Scientific report writing
- Data interpretation and result communication

---

## Author

**Sumed Seeyakmani Kuson**  
B.Sc. Data Science  
TU Dortmund University

---

## Disclaimer

The dataset used in this thesis contains clinical patient information and is not publicly shared in this repository. The repository is intended to showcase the analysis methodology, modelling workflow, and academic report.
