# Predicting-MacroB12-R-Analysis

Replication package for the study: 'Predicting MacroB12: Data-Driven Assessment of Pre-PEG B12 Utility and Clinical Caveats.'

## Description
This repository contains the R scripts and processed data necessary to reproduce the analyses and figures presented in the aforementioned study on MacroB12 prediction. The research identifies independent clinical and laboratory predictors of MacroB12 positivity and rigorously evaluates the discriminatory performance of pre-Polyethylene Glycol (PEG) B12 concentrations.

## Repository Structure

The repository is organized as follows:

-   `/R_Outputs_Data_Prep_New_30pct/`: Contains processed data files for analyses where MacroB12 positivity is defined by a <30% PEG recovery threshold.
    -   `mice_imputation_object.rds`: The multiple imputation (MICE) object generated for this threshold.
    -   `df_pre_imputation_data_30pct.csv`: The dataset before MICE imputation, used primarily for obtaining original valid N counts for descriptive statistics.
-   `/R_Outputs_Data_Prep_New_40pct/`: Contains processed data files for sensitivity analyses where MacroB12 positivity is defined by a <40% PEG recovery threshold.
    -   `mice_imputation_object.rds`: The MICE imputation object for this threshold.
    -   `df_pre_imputation_data_40pct.csv`: The dataset before MICE imputation for this threshold.
-   `01_DescriptiveStats_Tables_Figures_30pct.R`: R script to generate descriptive statistics (Table 1 from the main manuscript) for the 30% threshold data.
-   `01_DescriptiveStats_Tables_Figures_40pct.R`: R script to generate descriptive statistics (Table S3 from Supplementary Information) for the 40% threshold data.
-   `02_DescriptiveStats_BoxViolinPlots_30pct.R`: R script to generate box-and-violin plots (Figure 1 from the main manuscript) for the 30% threshold data.
-   `03_Main_Regression_MacroB12Status_30pct.R`: R script for the main multivariable regression models (Figure 2 from the main manuscript) for the 30% threshold.
-   `03_Main_Regression_MacroB12Status_40pct.R`: R script for multivariable regression models for the 40% threshold (results for Table S3).
-   `04_ROC_Analysis_VB12PrePEG_30pct.R`: R script for Receiver Operating Characteristic (ROC) analysis of pre-PEG B12 and other biomarkers (Figure 3 & Supplementary Figure S1) for the 30% threshold.
-   `04_ROC_Analysis_VB12PrePEG_40pct.R`: R script for ROC analysis for the 40% threshold.
-   `05_Supplementary_BetaRegression_RandomForest_30pct.R`: R script for supplementary analyses including Beta regression and Random Forest variable importance (Supplementary Figure S2) for the 30% threshold.
    <!-- Add similar entries for 40pct supplementary analyses if they exist -->
-   `.gitignore`: Specifies intentionally untracked files that Git should ignore (standard R template).
-   `LICENSE`: Contains the MIT License governing the use of the code in this repository.
- `README.md`: This file provides an overview and instructions.
-   `Data_Dictionary.md`: (Strongly Recommended) A file describing the variables in the provided CSV datasets. *(You will need to create this file).*

## Requirements

To execute the analysis scripts, the following are required:

-   **R**: Version 4.3.x or later is recommended.
-   **R Packages**:
    -   `tidyverse` (includes `ggplot2`, `dplyr`, `forcats`, `readr`, etc.)
    -   `data.table`
    -   `here`
    -   `mice`
    -   `sandwich`
    -   `lmtest`
    -   `logistf`
    -   `ggthemes`
    -   `RColorBrewer`
    -   `ggrepel`
    -   `ggdist`
    -   `knitr`
    -   `kableExtra`
    -   `skimr`
    -   `ggpubr`
    -   `paletteer`
    -   `patchwork`
    -   `pROC`
    -   `plotROC`
    -   `svglite`
    -   `betareg`
    -   `performance`
    -   `car`

    The scripts include a basic function at the beginning to attempt installation of any missing required packages from CRAN.

## Instructions for Use

1.  **Clone the Repository:**
    Open your terminal or Git client and clone this repository to your local machine:
    ```bash
    git clone https://github.com/renee29/Predicting-MacroB12-R-Analysis.git
    cd Predicting-MacroB12-R-Analysis
    ```

2.  **Set Up R Environment:**
    Open the R scripts in your preferred R environment (e.g., RStudio, VS Code with R extensions). The scripts use the `here` package to manage file paths, which relies on the cloned repository's root directory (identified by the `.git` folder). No manual setting of the working directory (`setwd()`) is required.

3.  **Install Packages:**
    Ensure all R packages listed under "Requirements" are installed. Running any of the analysis scripts will attempt to install missing packages automatically.

4.  **Run Analysis Scripts:**
    The processed data required for the analysis scripts are located in the `/R_Outputs_Data_Prep_New_30pct/` and `/R_Outputs_Data_Prep_New_40pct/` folders.
    Execute the R scripts (e.g., `01_DescriptiveStats_Tables_Figures_30pct.R`, `03_Main_Regression_MacroB12Status_30pct.R`, etc.) to reproduce the specific analyses and generate outputs. Outputs (tables, figures) will typically be saved into subfolders within an `R_Outputs_[AnalysisType]` directory (e.g., `R_Outputs_Table1/Threshold_30pct/`).

## Data

The datasets provided in the `/R_Outputs_Data_Prep_New_.../` folders are processed versions derived from an original clinical dataset. The original raw dataset is not publicly shared due to patient privacy restrictions and ethical considerations.

A **Data Dictionary** describing the variables in the `df_pre_imputation_data_XXpct.csv` files and implicitly in the MICE objects is crucial for understanding the data. *(Please create and add a `Data_Dictionary.md` file to this repository).*

## License

The R code in this repository is licensed under the MIT License - see the [LICENSE](LICENSE) file for complete details.
The processed data shared in this repository is permitted for research and replication purposes in conjunction with the accompanying code, subject to the terms of the study's ethical approvals.

## Citation

If you use this code or data in your research, please cite our study:

[Placeholder for your full paper citation once available. Example: Frías-Ruiz, C., Gálvez-Navas, J.M., et al. (Year). Predicting MacroB12: Data-Driven Assessment of Pre-PEG B12 Utility and Clinical Caveats. *Journal Name*, Volume(Issue), Pages. DOI: XXXXXX]

This repository is archived on Zenodo: https://doi.org/10.5281/zenodo.15420222

## Contact

For any questions regarding the code or data, please contact:
[Rene Fabregas] - [rfabregas@ugr.es]
