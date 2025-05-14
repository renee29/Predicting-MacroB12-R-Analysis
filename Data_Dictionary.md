# Data Dictionary: Predicting-MacroB12-R-Analysis

This document describes the variables found in the processed datasets (`df_pre_imputation_data_30pct.csv`, `df_pre_imputation_data_40pct.csv`, and implicitly in the `mice_imputation_object.rds` files) used for the study "Predicting MacroB12: Data-Driven Assessment of Pre-PEG B12 Utility and Clinical Caveats."

**Note:** The original raw data contained direct patient identifiers and detailed free-text clinical notes which have been removed or categorized for privacy and ethical considerations. The variables listed below are a result of this preprocessing.

---

## Demographic Variables

| Variable Name   | Type    | Description                                     | Units / Levels                                  | Notes                                                                      |
| :-------------- | :------ | :---------------------------------------------- | :---------------------------------------------- | :------------------------------------------------------------------------- |
| `Age`           | Numeric | Patient's age at the time of B12 assessment.    | Years                                           | Calculated from date of birth and analysis date.                           |
| `Sex`           | Factor  | Patient's biological sex.                       | `Female`, `Male`                                | `Female` is typically used as the reference level in models.             |
| `PatientSource` | Factor  | Origin of the patient sample/referral.          | `Ambulatory`, `Hospital`, `External Center`     | `Ambulatory` is typically used as the reference level.                   |

---

## Vitamin B12 Specific Variables

| Variable Name            | Type    | Description                                                                | Units    | Notes                                                                                                   |
| :----------------------- | :------ | :------------------------------------------------------------------------- | :------- | :------------------------------------------------------------------------------------------------------ |
| `VB12_PrePEG`            | Numeric | Serum Vitamin B12 concentration before Polyethylene Glycol (PEG) precipitation. | pg/mL    | Values originally reported as ">X" (e.g., >2000) were capped at X for pre-imputation data. Imputed in MICE. |
| `PEG_Recovery_Percent`   | Numeric | Percentage of B12 recovered in supernatant after PEG precipitation.        | %        | Calculated as `(VB12_PostPEG / VB12_PrePEG) * 100`. Original VB12_PostPEG values also processed.      |
| `MacroB12_Status`        | Factor  | MacroB12 status based on PEG recovery (<30% = Positive).                   | `Negative`, `Positive` | Primary outcome variable for 30% threshold analysis.                                            |
| `MacroB12_Status_40pct`  | Factor  | MacroB12 status based on PEG recovery (<40% = Positive).                   | `Negative`, `Positive` | Outcome variable for 40% threshold sensitivity analysis.                                       |

---

## Concomitant Laboratory Markers

| Variable Name | Type    | Description                      | Units   | Notes                                                                                                |
| :------------ | :------ | :------------------------------- | :------ | :--------------------------------------------------------------------------------------------------- |
| `CRP`         | Numeric | C-Reactive Protein level.        | mg/L    | Imputed in MICE if missing. Original "no tiene" or "<X" values were handled.                        |
| `Folate`      | Numeric | Serum Folate level.              | ng/mL   | Imputed in MICE if missing. Original "no tiene" or ">X" values were handled.                         |
| `Hemoglobin`  | Numeric | Hemoglobin level.                | g/dL    | Imputed in MICE if missing.                                                                          |
| `MCV`         | Numeric | Mean Corpuscular Volume.         | fL      | Imputed in MICE if missing.                                                                          |
| `RDW`         | Numeric | Red Cell Distribution Width.     | %       | Imputed in MICE if missing.                                                                          |

---

## Diagnostic Category Flags

These variables are binary flags (0 = Absent, 1 = Present) indicating the presence of a diagnosis within a broad category. These were derived from original coded classifications and/or processed free-text diagnoses. A patient may have flags in multiple categories.

| Variable Name                 | Type   | Description                                       | Levels | Notes                                                                      |
| :---------------------------- | :----- | :------------------------------------------------ | :----- | :------------------------------------------------------------------------- |
| `Diag_Hematologic`            | Factor | Presence of a hematologic condition/diagnosis.    | `0`, `1` |                                                                            |
| `Diag_Hepatic`                | Factor | Presence of a hepatic condition/diagnosis.        | `0`, `1` |                                                                            |
| `Diag_Oncologic`              | Factor | Presence of an oncologic condition/diagnosis.     | `0`, `1` |                                                                            |
| `Diag_Rheum_Autoimmune`       | Factor | Presence of a rheumatologic or autoimmune disease. | `0`, `1` |                                                                            |
| `Diag_Infectious`             | Factor | Presence of an infectious disease diagnosis.      | `0`, `1` |                                                                            |
| `Diag_Renal`                  | Factor | Presence of a renal condition/diagnosis.          | `0`, `1` |                                                                            |
| `Diag_GI_Suspicion`           | Factor | Presence of a GI disorder or suspicion of neoplasm. | `0`, `1` |                                                                            |
| `Diag_Hypertension`           | Factor | Presence of a hypertension diagnosis.             | `0`, `1` |                                                                            |
| `Diag_Diabetes`               | Factor | Presence of a diabetes mellitus diagnosis.        | `0`, `1` |                                                                            |
| `Diag_Healthy_Unspecified`  | Factor | Patient considered healthy or diagnosis unspecified, and no other major diagnostic flags present. | `0`, `1` | This category is mutually exclusive with the other Diag_* flags being 1. |

---

**Notes on Data Processing:**

-   **Missing Data:** Missing values for `CRP`, `Folate`, `Hemoglobin`, `MCV`, and `RDW` in the original dataset were handled using Multiple Imputation by Chained Equations (MICE). The provided `.rds` files contain the MICE objects, and the analysis scripts typically use the first imputed dataset. The `.csv` files (`df_pre_imputation_data_...csv`) represent the data *before* MICE imputation.
-   **Laboratory Value Cleaning:** Original laboratory values reported with qualifiers (e.g., "<0.05", ">2000", "no tiene") were converted to numerical representations or NA prior to imputation. For instance, "no tiene" was treated as NA. Values like ">2000" were often capped at 2000 for the pre-imputation dataset, and the imputation process handled these.
-   **Categorization:** Original free-text diagnoses and more granular coded diagnoses were mapped into the broader `Diag_*` categories listed above. The specific keywords and logic used for this mapping can be inferred from the (optional) `01_Data_Ingestion_Cleaning_...R` scripts if those are also shared for full transparency on pre-processing.