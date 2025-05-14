# ===========================================================
#        Script R para Forest Plot Combinado (Unadj & Adj) OR/PR
#        (Estilo Publicación con Datos MICE Pooled y Orden Tabla 1)
# ===========================================================

# 0. CONFIGURACIÓN E IMPORTACIÓN
# ===========================================================
# --- Cargar Paquetes ---
required_packages <- c("data.table", "tidyverse", "ggplot2", "here", "rstudioapi",
                       "mice", "sandwich", "lmtest", "logistf",
                       "ggthemes", "RColorBrewer", "ggrepel", "ggdist",
                       "knitr", "kableExtra", "skimr", "ggpubr", "paletteer",
                       "dplyr", "stats", "forcats", "patchwork")

install_if_missing <- function(pkg) { if (!requireNamespace(pkg, quietly = TRUE)) { cat(sprintf("Installing package: %s\n", pkg)); tryCatch(install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org/"), error = function(e){ cat(sprintf("  Failed to install %s: %s\n", pkg, e$message))}) } }
sapply(required_packages, install_if_missing)

library(data.table); library(tidyverse); library(ggplot2); library(here); library(rstudioapi)
library(mice); library(sandwich); library(lmtest); # library(logistf)
library(ggthemes); library(RColorBrewer); library(ggrepel); library(ggdist)
library(knitr); library(kableExtra); library(skimr); library(ggpubr); library(paletteer)
library(dplyr); library(stats); library(forcats); library(patchwork)

# --- Directorios ---
script_dir <- tryCatch({ dirname(rstudioapi::getSourceEditorContext()$path) }, error = function(e) { getwd() })
input_data_prep_dir <- file.path(script_dir, "R_Outputs_Data_Prep_New_40pct")
# Directorio de salida para los plots OR/PR combinados
output_analysis_dir <- file.path(script_dir, "Forest_Plots_OR_PR_40pct") # <<--- NUEVO NOMBRE DIRECTORIO
if (!dir.exists(output_analysis_dir)) { dir.create(output_analysis_dir) } else { cat("Output directory exists:", normalizePath(output_analysis_dir), "\n") }

# --- Constantes ---
RANDOM_STATE <- 42; set.seed(RANDOM_STATE)
PR_THRESHOLD_PREVALENCE <- 0.20

# ===========================================================
#        Mapeo de Etiquetas COMPLETAS (Revisado)
#        Claves = Nombres EXACTOS de los TÉRMINOS del Modelo (Variable_Raw)
#        Valores = Etiquetas Descriptivas para Publicación
# ===========================================================
LABEL_MAP <- list(
  # --- Términos del Modelo (Claves EXACTAS de Variable_Raw) ---
  'EDAD' = 'Age (years)',
  'VB12.PRE.PEG' = 'VB12 Pre-PEG (pg/mL)', # Clave con punto como en Variable_Raw
  'PCR' = 'CRP (mg/L)',
  'FOL' = 'Serum Folate (ng/mL)',
  'HB' = 'Hemoglobin (g/dL)',
  'VCM' = 'MCV (fL)',
  'ADE' = 'RDW (%)',
  'SEXOMale' = 'Sex: Male vs Female',         # Clave para el contraste
  'ProcedenciaExternal Center' = 'Source: External Center vs Ambulatory', # Clave para el contraste
  'ProcedenciaHospital' = 'Source: Hospital vs Ambulatory', # Clave para el contraste
  'Diag_Hematologic1' = 'Diagnosis: Hematologic: Present vs Absent',
  'Diag_Hepatic1' = 'Diagnosis: Hepatic: Present vs Absent',
  'Diag_Oncologic1' = 'Diagnosis: Oncologic: Present vs Absent',
  'Diag_Rheum_Autoimmune1' = 'Diagnosis: Rheum./Autoimmune: Present vs Absent',
  'Diag_Infectious1' = 'Diagnosis: Infectious: Present vs Absent',
  'Diag_Renal1' = 'Diagnosis: Renal: Present vs Absent',
  'Diag_GI_Suspicion1' = 'Diagnosis: GI/Neo. Suspicion: Present vs Absent',
  'Diag_Hypertension1' = 'Diagnosis: Hypertension: Present vs Absent',
  'Diag_Diabetes1' = 'Diagnosis: Diabetes: Present vs Absent',
  'Diag_Healthy_Unspecified1' = 'Diagnosis: Healthy/Unspecified: Present vs Absent',

  # --- Nombres Originales/Base (Opcional pero útil si se usan en otros lugares) ---
  'edad' = 'Age (years)',               # Nombre original (si se usa)
  'sexo' = 'Sex',                       # Nombre original (si se usa)
  'procedencia.' = 'Source',              # Nombre original con punto (si se usa)
  'vb12_pre_peg' = 'VB12 Pre-PEG (pg/mL)', # Nombre original (si se usa)
  'pcr' = 'CRP (mg/L)',
  'fol' = 'Serum Folate (ng/mL)',
  'hb' = 'Hemoglobin (g/dL)',
  'vcm' = 'MCV (fL)',
  'ade' = 'RDW (%)',
  'Diag_Hematologic' = 'Diagnosis: Hematologic', # Nombre base del flag
  'Diag_Hepatic' = 'Diagnosis: Hepatic',
  'Diag_Oncologic' = 'Diagnosis: Oncologic',
  'Diag_Rheum_Autoimmune' = 'Diagnosis: Rheum./Autoimmune',
  'Diag_Infectious' = 'Diagnosis: Infectious',
  'Diag_Renal' = 'Diagnosis: Renal',
  'Diag_GI_Suspicion' = 'Diagnosis: GI/Neo. Suspicion',
  'Diag_Hypertension' = 'Diagnosis: Hypertension',
  'Diag_Diabetes' = 'Diagnosis: Diabetes',
  'Diag_Healthy_Unspecified' = 'Diagnosis: Healthy/Unspecified',

  # --- Otros posibles términos o niveles (Añadir si es necesario) ---
  'MacroB12_Positive' = 'MacroB12 Status', # La variable dependiente
  'Male' = 'Male',                         # Nivel de factor
  'Female' = 'Female',                     # Nivel de factor (Referencia)
  'Ambulatory' = 'Ambulatory',               # Nivel de factor (Referencia)
  'Hospital' = 'Hospital',                 # Nivel de factor
  'External Center' = 'External Center'    # Nivel de factor
  # ... añadir cualquier otro término o nivel que necesite etiqueta ...
)


# --- Mapeo de Etiquetas COMPACTAS (Keys = Nombres INTERNOS/Términos Modelo) ---
# !! ASEGÚRATE QUE LAS CLAVES COINCIDAN EXACTAMENTE CON LOS TÉRMINOS DE TU MODELO !!
# !! AJUSTA LAS CLAVES PARA QUE COINCIDAN EXACTAMENTE CON LOS ERRORES/SALIDA !!
COMPACT_LABEL_MAP <- list(
  # Variables (Claves en MAYÚSCULAS como en los errores)
  'EDAD' = 'Age',
  'VB12.PRE.PEG' = 'VB12 Pre',
  'PCR' = 'CRP',
  'FOL' = 'Folate',
  'HB' = 'Hb',
  'VCM' = 'MCV',
  'ADE' = 'RDW',
  # Contrastes (Claves EXACTAS de los errores/salida del modelo)
  'SEXOMale'='Sex (M vs F)',
  'ProcedenciaHospital'='Source (Hosp vs Amb)', # Asumiendo sin punto basado en error
  'ProcedenciaExternal Center'='Source (Ext vs Amb)', # Asumiendo sin punto basado en error
  # Flags Diag (Clave con '1' al final)
  'Diag_Hematologic1'='Diag: Hematologic',
  'Diag_Hepatic1'='Diag: Hepatic',
  'Diag_Oncologic1'='Diag: Oncologic',
  'Diag_Rheum_Autoimmune1'='Diag: Rheum/AI',
  'Diag_Infectious1'='Diag: Infectious',
  'Diag_Renal1'='Diag: Renal',
  'Diag_GI_Suspicion1'='Diag: GI/Neo Susp.',
  'Diag_Hypertension1'='Diag: Hypertension',
  'Diag_Diabetes1'='Diag: Diabetes',
  'Diag_Healthy_Unspecified1'='Diag: Healthy/Unspec.'
)


PLOT_ORDER_VARIABLES <- c(
    "EDAD",                     # Corregido
    "SEXOMale",                 # Corregido
    "ProcedenciaHospital",      # Corregido
    "ProcedenciaExternal Center", # Corregido
    "Diag_Healthy_Unspecified1",
    "Diag_Hematologic1",
    "Diag_Hepatic1",
    "Diag_Oncologic1",
    "Diag_Rheum_Autoimmune1",
    "Diag_Infectious1",
    "Diag_Renal1",
    "Diag_GI_Suspicion1",
    "Diag_Hypertension1",
    "Diag_Diabetes1",
    "VB12.PRE.PEG",             # Corregido
    "PCR",                      # Corregido
    "FOL",                      # Corregido
    "HB",                       # Corregido
    "VCM",                      # Corregido
    "ADE"                       # Corregido
)

# --- Helper Functions ---
# (Sin cambios)
get_label <- function(name, default_fmt=TRUE) { name_str <- as.character(name); label <- LABEL_MAP[[name_str]]; if (!is.null(label)) return(label); clean_name <- gsub("`", "", name_str); if (default_fmt) return(tools::toTitleCase(gsub("[._]", " ", clean_name))) else return(clean_name) }
get_compact_label <- function(name) { name_str <- as.character(name); compact_label <- COMPACT_LABEL_MAP[[name_str]]; if (!is.null(compact_label)) { return(compact_label) } else { warning(paste("No compact label for:", name_str,". Using Variable_Raw.")); return(name_str) } }
format_p_value <- function(p_vector) { sapply(p_vector, function(p) { if (is.na(p)) { return("NA") } else if (p < 0.0001) { return(sprintf("%.1e", p)) } else if (p < 0.001) { return("<0.001") } else if (p < 0.01) { return(sprintf("%.3f", p)) } else { return(sprintf("%.2f", p)) } }) }
add_stars <- function(p_vector) { sapply(p_vector, function(p) { if (is.na(p) || !is.numeric(p)) return('') ; if (p < 0.001) return('***') ; if (p < 0.01) return('**') ; if (p < 0.05) return('*') ; return('') }) }
sanitize_filename <- function(name) { name_sanitized <- gsub('[<>:"/\\|?*%[:cntrl:]]', '_', as.character(name)); name_sanitized <- gsub('[[:punct:]]', '_', name_sanitized); name_sanitized <- gsub('\\s+', '_', name_sanitized); name_sanitized <- gsub('_+', '_', name_sanitized); name_sanitized <- gsub('^_|_$', '', name_sanitized); name_sanitized <- substr(name_sanitized, 1, 100); return(name_sanitized) }


# =======================================================================
# 1. CARGAR DATOS MICE Y PREPARAR PARA GLM
# =======================================================================
# (Código sin cambios)
cat("\n--- 1. Loading MICE Data Object & Preparing for GLM ---\n"); mice_rds_file <- file.path(input_data_prep_dir, "mice_imputation_object_internalNames.rds"); mids_data_original <- NULL; if (file.exists(mice_rds_file)) { tryCatch({ mids_data_original <- readRDS(mice_rds_file); cat("Successfully loaded MICE object (.rds) with", mids_data_original$m, "imputations.\n") }, error = function(e) { cat("Error loading MICE object:", e$message, "\n") }) }; if (is.null(mids_data_original)) { stop("Could not load MICE object.") }
mids_data_glm <- mids_data_original; target_col_internal <- make.names('MacroB12_Positive'); if (!target_col_internal %in% names(mids_data_glm$data)) { stop(paste("Target column '", target_col_internal, "' not found.")) }
cat(sprintf("Converting target variable '%s' to numeric 0/1...\n", target_col_internal)); for (i in 1:mids_data_glm$m) { imp_data <- complete(mids_data_glm, action = i); imp_data[[target_col_internal]] <- as.numeric(imp_data[[target_col_internal]] == "Positive"); for(colname in names(imp_data)){ if(colname %in% names(mids_data_glm$imp)) { mids_data_glm$imp[[colname]][,i] <- imp_data[[colname]][is.na(mids_data_glm$data[[colname]])] } } }; mids_data_glm$data[[target_col_internal]] <- as.numeric(mids_data_glm$data[[target_col_internal]] == "Positive"); cat(sprintf("Target variable '%s' converted.\n", target_col_internal))
cat("Verifying/setting factor types...\n"); sex_col_internal <- make.names("sexo"); proc_col_internal <- make.names("procedencia."); if(sex_col_internal %in% names(mids_data_glm$data)) { mids_data_glm <- mice::complete(mids_data_glm, action="long", include=TRUE) %>% mutate(!!sym(sex_col_internal) := factor(!!sym(sex_col_internal))) %>% as.mids(); if("Female" %in% levels(mids_data_glm$data[[sex_col_internal]])) mids_data_glm$data[[sex_col_internal]] <- relevel(mids_data_glm$data[[sex_col_internal]], ref="Female") }
if(proc_col_internal %in% names(mids_data_glm$data)) { mids_data_glm <- mice::complete(mids_data_glm, action="long", include=TRUE) %>% mutate(!!sym(proc_col_internal) := factor(!!sym(proc_col_internal), levels = c("Ambulatory", "External Center", "Hospital"))) %>% as.mids(); if("Ambulatory" %in% levels(mids_data_glm$data[[proc_col_internal]])) mids_data_glm$data[[proc_col_internal]] <- relevel(mids_data_glm$data[[proc_col_internal]], ref="Ambulatory") }
diag_cols_internal <- names(mids_data_glm$data)[startsWith(names(mids_data_glm$data), "Diag_")]; for(col in diag_cols_internal){ if(col %in% names(mids_data_glm$data)){ mids_data_glm <- mice::complete(mids_data_glm, action="long", include=TRUE) %>% mutate(!!sym(col) := factor(!!sym(col), levels=c(0,1))) %>% as.mids(); mids_data_glm$data[[col]] <- factor(mids_data_glm$data[[col]], levels=c(0,1)) } }; cat("Factor types verified/set.\n")

# =======================================================================
# 2. PREPARAR VARIABLES PARA ANÁLISIS (Usando mids_data_glm)
# =======================================================================
# (Código sin cambios)
cat("\n--- 2. Preparing Variables for Analysis ---\n"); target_col_model <- target_col_internal; predictor_vars_internal <- names(mids_data_glm$data); predictor_vars_internal <- setdiff(predictor_vars_internal, target_col_model); predictor_vars_internal <- setdiff(predictor_vars_internal, make.names(c("vb12_post_peg", "Recovery_pct"))); cat("Excluding 'Recovery_pct' from predictor list.\n")
constant_vars <- names(which(sapply(mids_data_glm$data[, predictor_vars_internal], function(x) n_distinct(x, na.rm = TRUE)) <= 1)); if(length(constant_vars) > 0) { cat("Removing constant predictors:", paste(constant_vars, collapse=", "), "\n"); predictor_vars_internal <- setdiff(predictor_vars_internal, constant_vars) }; predictor_vars_model <- predictor_vars_internal; cat("Target variable:", target_col_model, "\nPredictors:\n", paste(predictor_vars_model, collapse=", "), "\n")

# =======================================================================
# 3. CÁLCULO DE PREVALENCIA
# =======================================================================
# (Código sin cambios)
cat("\n--- 3. Calculating Average Outcome Prevalence ---\n"); 
prevalence <- mean(sapply(complete(mids_data_glm, "all"), function(df) mean(df[[target_col_model]], na.rm = TRUE))); 
if (!is.na(prevalence)) { cat(sprintf("Average Outcome Prevalence (%s == 1): %.1f%%\n", target_col_model, prevalence * 100)); 
calculate_pr <- prevalence > PR_THRESHOLD_PREVALENCE; 
if (calculate_pr) { cat("Prevalence > threshold. PRs will be calculated.\n") } else { cat("Prevalence <= threshold. PRs calculation skipped.\n") } } else { cat("Could not calculate prevalence.\n"); 
calculate_pr <- FALSE }


# =======================================================================
# 4. EJECUTAR MODELOS Y OBTENER RESULTADOS POOLED
# =======================================================================
# (Código sin cambios - Usa la función run_and_pool actualizada que incluye Compact_Var_Label)
cat("\n--- 4. Running Univariate & Multivariate Models (ORs and PRs) ---\n")
run_and_pool <- function(formula_str, family_type, data_mids, exponentiate_est = TRUE) { results_list <- list() ; fit_pooled <- NULL; summary_pooled <- NULL; tryCatch({ fit_models <- with(data = data_mids, expr = glm(formula = as.formula(formula_str), family = family_type)); if(is.null(fit_models) || !inherits(fit_models, "mira")) { stop("'with' did not return 'mira' object.") }; if(length(fit_models$analyses) == 0) { stop("No analyses in 'mira' object.") }; errors_in_analyses <- sapply(fit_models$analyses, function(x) inherits(x, "try-error")); if(any(errors_in_analyses)) { cat("  WARNING: Errors in individual GLM fits.\n") }; fit_pooled <- pool(fit_models); summary_pooled <- summary(fit_pooled, conf.int = TRUE, exponentiate = exponentiate_est); if (!is.null(summary_pooled) && nrow(summary_pooled) > 0) { for (i in 1:nrow(summary_pooled)) { term_row <- summary_pooled[i, ]; term_name <- term_row$term; if (term_name == "(Intercept)") next; est <- term_row$estimate; ci_l <- term_row$conf.low; ci_h <- term_row$conf.high; p_val <- term_row$p.value; if(all(is.finite(c(est, ci_l, ci_h, p_val)))){ results_list[[length(results_list) + 1]] <- list(Variable_Raw=term_name, Variable=get_label(term_name), Compact_Var_Label=get_compact_label(term_name), Estimate=est, CI_Low=ci_l, CI_High=ci_h, P_Value=p_val) } else { cat(sprintf("    Skipping term '%s' due to non-finite results.\n", term_name)) } } } else { cat(sprintf("    Pooling resulted in NULL/empty summary: %s\n", formula_str)) } }, error = function(e) { cat(sprintf("  ERROR pooling for %s: %s\n", formula_str, e$message)) }, warning = function(w) { cat(sprintf("  WARNING pooling for %s: %s\n", formula_str, w$message)) }); if (length(results_list) > 0) { return(bind_rows(results_list)) } else { return(data.frame(Variable_Raw=character(), Variable=character(), Compact_Var_Label=character(), Estimate=numeric(), CI_Low=numeric(), CI_High=numeric(), P_Value=numeric(), stringsAsFactors = FALSE)) } }
cat("Running Univariate Models...\n"); family_logistic <- binomial(link="logit"); univariate_or_results_list <- lapply(predictor_vars_model, function(pred) { formula_uni_str <- paste0("`", target_col_model, "` ~ `", pred, "`"); run_and_pool(formula_uni_str, family_logistic, mids_data_glm, exponentiate_est = TRUE) }); univariate_or_df_pooled <- bind_rows(univariate_or_results_list) %>% rename(OR = Estimate) %>% filter(!is.na(OR))
univariate_pr_df_pooled <- data.frame(); if(calculate_pr){ cat("Running Univariate Poisson Models for PRs...\n"); family_poisson <- poisson(link="log"); univariate_pr_results_list <- lapply(predictor_vars_model, function(pred) { formula_uni_str <- paste0("`", target_col_model, "` ~ `", pred, "`"); run_and_pool(formula_uni_str, family_poisson, mids_data_glm, exponentiate_est = TRUE) }); univariate_pr_df_pooled <- bind_rows(univariate_pr_results_list) %>% rename(PR = Estimate) %>% filter(!is.na(PR)) }
cat("Running Multivariate Models...\n"); multivariate_or_df_pooled <- data.frame(); multivariate_pr_df_pooled <- data.frame()
if (length(predictor_vars_model) > 1) { formula_multi_str <- paste0("`", target_col_model, "` ~ ", paste(paste0("`", predictor_vars_model, "`"), collapse = " + ")); cat("Formula:", formula_multi_str, "\n"); cat("  Multivariate Logistic (ORs)...\n"); multivariate_or_df_pooled <- run_and_pool(formula_multi_str, family_logistic, mids_data_glm, exponentiate_est = TRUE) %>% rename(OR = Estimate) %>% filter(!is.na(OR)); if(calculate_pr){ cat("  Multivariate Poisson (PRs)...\n"); multivariate_pr_df_pooled <- run_and_pool(formula_multi_str, family_poisson, mids_data_glm, exponentiate_est = TRUE) %>% rename(PR = Estimate) %>% filter(!is.na(PR)) } } else { cat("Not enough predictors.\n")}
cat("\n--- Inspecting Pooled Results Before Plotting ---\n"); cat("== Univariate ORs ==\n"); print(univariate_or_df_pooled); cat("\n== Multivariate ORs ==\n"); print(multivariate_or_df_pooled); if(calculate_pr) { cat("\n== Univariate PRs ==\n"); print(univariate_pr_df_pooled); cat("\n== Multivariate PRs ==\n"); print(multivariate_pr_df_pooled) }; cat("--------------------------------------------------\n")
if(nrow(univariate_or_df_pooled) > 0) write.csv(univariate_or_df_pooled, file.path(output_analysis_dir, "Table_ORs_Univariate_Pooled.csv"), row.names=F, na="")
if(nrow(multivariate_or_df_pooled) > 0) write.csv(multivariate_or_df_pooled, file.path(output_analysis_dir, "Table_ORs_Multivariate_Pooled.csv"), row.names=F, na="")
if(calculate_pr && nrow(univariate_pr_df_pooled) > 0) write.csv(univariate_pr_df_pooled, file.path(output_analysis_dir, "Table_PRs_Univariate_Pooled.csv"), row.names=F, na="")
if(calculate_pr && nrow(multivariate_pr_df_pooled) > 0) write.csv(multivariate_pr_df_pooled, file.path(output_analysis_dir, "Table_PRs_Multivariate_Pooled.csv"), row.names=F, na="")


# =======================================================================
# 5. COMBINAR RESULTADOS Y GENERAR FOREST PLOT COMBINADO (MEJORADO)
# =======================================================================
cat("\n--- 5. Combining Results and Generating Enhanced Forest Plots ---\n")

# --- Paleta de colores 'ggthemes::calc' ---
# (Sin cambios desde tu script original)
#palette_calc <- paletteer_d("ggthemes::calc", n = 3)[1:2]
#color_unadjusted <- palette_calc[1] # Primer color (ej. coral/rojo)
#color_adjusted <- palette_calc[2]   # Segundo color (ej. azul)

palette_econ <- paletteer_d("ggthemes::stata_s2color", n = 3)[1:2]
color_unadjusted <- palette_econ[1] # Primer color (ej. coral/rojo)
color_adjusted <- palette_econ[2]   # Segundo color (ej. azul)

# --- FUENTE (Define la familia de fuente deseada - Asegúrate de que esté disponible) ---
font_family <- "Arial" # O "Helvetica", "Calibri", etc.

# --- Función Forest Plot Mejorada (Estilo Publicación con Tabla Integrada) ---
create_enhanced_combined_forest_plot <- function(unadj_df, adj_df, estimate_col_name, plot_title, filename_base) {

    estimate_sym <- sym(estimate_col_name)
    # Asegurarse que la columna compacta existe o crearla como fallback
    if (!"Compact_Var_Label" %in% names(unadj_df)) { unadj_df <- unadj_df %>% mutate(Compact_Var_Label = Variable_Raw) }
    if (!"Compact_Var_Label" %in% names(adj_df)) { adj_df <- adj_df %>% mutate(Compact_Var_Label = Variable_Raw) }

    unadj_df_mod <- unadj_df %>% rename(Estimate = !!estimate_sym) %>% mutate(Model = "Unadjusted")
    adj_df_mod <- adj_df %>% rename(Estimate = !!estimate_sym) %>% mutate(Model = "Adjusted")

    combined_df <- bind_rows(unadj_df_mod, adj_df_mod)

    # Asegurar que Compact_Var_Label no sea NA ANTES de filtrar
    combined_df <- combined_df %>% mutate(Compact_Var_Label = ifelse(is.na(Compact_Var_Label) | Compact_Var_Label == "", Variable_Raw, Compact_Var_Label))

    # --- PREPARACIÓN DE DATOS PARA EL PLOT Y LA TABLA ---
    combined_df <- combined_df %>%
        mutate(
            # Formato P-valor mejorado
            P_Value_Formatted = case_when(
                 is.na(P_Value)   ~ "NA",
                 P_Value < 0.001 ~ "<0.001",
                 P_Value < 0.01  ~ sprintf("%.3f", P_Value),
                 TRUE            ~ sprintf("%.2f", P_Value)
            ),
            # Formato Estimado [IC 95%] para tabla (modelo ajustado)
            Estimate_CI_Formatted = ifelse(Model == "Adjusted",
                                           sprintf("%.2f (%.2f–%.2f)", Estimate, CI_Low, CI_High), # Usar guion largo (–)
                                           NA_character_),
             # Añadir estrellas para tabla (modelo ajustado)
            Stars = ifelse(Model == "Adjusted", add_stars(P_Value), NA_character_),
            P_Value_Formatted_With_Stars = ifelse(Model == "Adjusted" & !is.na(P_Value),
                                                  paste0(P_Value_Formatted, Stars),
                                                  NA_character_),
            # Columna para poner en negrita resultados significativos (modelo ajustado)
            Is_Significant = ifelse(Model == "Adjusted" & !is.na(P_Value) & P_Value < 0.05, "bold", "plain")
        ) %>%
        # Filtro de validez general
        filter(!is.na(Estimate) & is.finite(Estimate) &
               !is.na(CI_Low) & is.finite(CI_Low) &
               !is.na(CI_High) & is.finite(CI_High))
               # No filtrar por P-valor aquí, la tabla mostrará "NA" si aplica

    # --- Comprobar escala logarítmica ---
    use_log_scale <- TRUE
    # Permitir CIs <= 0 solo si NO usamos escala log
    if(any(combined_df$CI_Low <= 0 | combined_df$CI_High <= 0, na.rm = TRUE)) {
        cat("Warning: Non-positive CIs found for", plot_title,". Switching to linear scale.\n");
        use_log_scale <- FALSE
    }
    if(use_log_scale) {
        # Filtrar datos con CI no positivos SOLO si se usa escala log
        combined_df <- combined_df %>% filter(CI_Low > 0 & CI_High > 0)
    }

    if(nrow(combined_df) == 0) { cat(sprintf("  Skipping plot '%s': No valid combined data after filtering.\n", plot_title)); return(NULL) }

    # --- ORDENACIÓN según PLOT_ORDER_VARIABLES ---
    # Asegurar que todas las Variable_Raw existan en los datos antes de convertirlas a factor
    available_vars_in_order <- intersect(PLOT_ORDER_VARIABLES, unique(combined_df$Variable_Raw))
    if(length(available_vars_in_order) == 0) {
         cat(sprintf("  Skipping plot '%s': No variables found matching PLOT_ORDER_VARIABLES.\n", plot_title));
         return(NULL)
    }
    plot_data_ordered <- combined_df %>%
        filter(Variable_Raw %in% available_vars_in_order) %>%
        mutate(
            # Crear factor ordenado para el eje Y usando etiquetas compactas
            compact_var_ordered = factor(Variable_Raw,
                                         levels = rev(available_vars_in_order), # Usar solo variables presentes, en orden inverso
                                         labels = rev(sapply(available_vars_in_order, get_compact_label)) # Obtener etiquetas compactas
                                         ),
            Model = factor(Model, levels = c("Unadjusted", "Adjusted"))
        ) %>%
         # Eliminar filas donde la etiqueta compacta no se pudo crear (aunque get_compact_label tiene fallback)
        filter(!is.na(compact_var_ordered))

    if(nrow(plot_data_ordered) == 0 || nlevels(plot_data_ordered$compact_var_ordered) == 0) {
         cat(sprintf("  Skipping plot '%s': No data left after ordering/labeling.\n", plot_title));
         return(NULL)
    }

   # --- Cálculo de Límites y Escalas (Eje X) ---

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# MODIFICACIÓN PARA LÍMITES FIJOS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# 1. Forzar el uso de escala logarítmica porque el rango 0.001-10 lo requiere.
use_log_scale <- TRUE # Asegúrate de que esta variable se establezca a TRUE antes de este bloque.
                     # En tu script original, 'use_log_scale' se determina antes.
                     # Si quieres que SIEMPRE sea logarítmica para estos límites, puedes ponerlo aquí.

# 2. Definir tus límites fijos
fixed_min_x <- 0.1
fixed_max_x <- 100

# 3. Generar breaks para la escala logarítmica dentro de los límites fijos
#    Puedes ajustar esta lógica para obtener los breaks que desees.
major_breaks <- 10^(floor(log10(fixed_min_x)):ceiling(log10(fixed_max_x))) # e.g., 0.001, 0.01, 0.1, 1, 10

# Opcional: Añadir minor breaks si lo deseas (ejemplo)
minor_breaks_list <- list()
if (length(major_breaks) > 0) {
    for(p_major in major_breaks){
        if (p_major < fixed_max_x) { # Solo añadir si p_major no es el último break
             # Añadir algunos puntos entre p_major y el siguiente orden de magnitud
             minor_breaks_list[[length(minor_breaks_list) + 1]] <- p_major * c(2, 5) # Ejemplo: 0.002, 0.005, 0.02, 0.05, etc.
        }
    }
}
all_potential_breaks <- sort(unique(c(major_breaks, unlist(minor_breaks_list))))
# Filtrar breaks para que estén estrictamente dentro de los límites o sean los límites mismos
valid_breaks <- all_potential_breaks[all_potential_breaks >= fixed_min_x & all_potential_breaks <= fixed_max_x]

# Asegurar que 1 esté en los breaks si está dentro del rango
if (1 >= fixed_min_x && 1 <= fixed_max_x && !(1 %in% valid_breaks)) {
    valid_breaks <- sort(unique(c(1, valid_breaks)))
}
# Si después de filtrar no quedan breaks (poco probable con 0.001-10), al menos usar los límites
if(length(valid_breaks) == 0) {
    valid_breaks <- sort(unique(c(fixed_min_x, fixed_max_x, if(1 >= fixed_min_x && 1 <= fixed_max_x) 1 else NULL)))
}
# Limitar el número de breaks para evitar solapamiento, si es necesario
if(length(valid_breaks) > 10) valid_breaks <- pretty(valid_breaks, n=7)


# 4. Crear la escala del eje X con los límites y breaks fijos
x_scale <- scale_x_log10(
    limits = c(fixed_min_x, fixed_max_x), # Aplicar límites fijos
    breaks = valid_breaks,
    labels = scales::label_number(accuracy = 0.1, trim = TRUE) # Ajustar accuracy para mostrar 0.001
)
vline_intercept <- 1 # Línea vertical en 1 para OR/PR
plot_xlab <- paste(estimate_col_name, "(log scale)") # Etiqueta del eje X

# La parte 'else' de la condición if(use_log_scale) ya no sería necesaria si siempre usas estos límites fijos
# y por lo tanto siempre usas escala logarítmica. Si quieres mantener la opción de escala lineal
# para otros casos, tendrías que envolver esta modificación en otra condición.

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# FIN DE LA MODIFICACIÓN
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    # --- Altura y parámetros del Plot ---
    n_vars <- nlevels(plot_data_ordered$compact_var_ordered)
    plot_height <- max(4, n_vars * 0.35 + 2.0) # Ajustar multiplicador y base según necesidad
    dodge_width <- 0.6
    base_size <- 9 # Tamaño base para fuentes (ajustable)
nx=1
ny=0 #Forzar límites Y para alineación con p_table
           # --- GRÁFICO PRINCIPAL (Puntos y Líneas) ---
    p_main <- ggplot(plot_data_ordered, aes(x = Estimate, y = compact_var_ordered, color = Model)) +
        # Geoms: línea vertical, barras de error, puntos
        geom_vline(xintercept = vline_intercept, linetype = "dashed", color = "black", linewidth = 0.4) +
        geom_errorbarh(aes(xmin = CI_Low, xmax = CI_High),
                       height = 0, # Sin terminaciones
                       linewidth = 0.75, # Grosor línea error
                       position = position_dodge(width = dodge_width)) +
        geom_point(shape = 16, # Círculo relleno
                   size = 3.5,  # Tamaño punto
                   position = position_dodge(width = dodge_width)) +  
                         
        # Escalas: Eje X y Color
        x_scale +
        scale_color_manual(values = c("Unadjusted" = color_unadjusted, "Adjusted" = color_adjusted), name = "Model:") +
        # Labs: Título Eje X, sin título Eje Y
        labs(x = plot_xlab, y = NULL ) +
        # Forzar límites Y para alineación con p_table
        coord_cartesian(
            ylim = c(nx, n_vars + ny),
            clip = "off"
        ) +
        # --- TEMA (Estilo Publicación Limpio con Leyenda Horizontal Abajo) ---
        #if (FALSE){ 
        theme_classic(base_size = base_size, base_family = font_family) +
        theme(
               # Panel y Rejilla
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               # Definir explícitamente las líneas de eje deseadas
               # Ejes: theme_classic da líneas X inferior e Y izquierda por defecto.
               # No necesitamos definir axis.line.x.top aquí porque la dibujamos con annotate.
               #axis.line.x.bottom = element_line(colour = "black", linewidth = 0.4), # Puedes sobreescribir si quieres cambiar estilo
               #axis.line.y.left   = element_line(colour = "black", linewidth = 0.4), # Puedes sobreescribir si quieres cambiar estilo
               # Ejes
               axis.line = element_line(colour = "black", linewidth = 0.4),
               # Leyenda Horizontal Abajo
               axis.title.x = element_text(size = rel(1.5), margin = margin(t = 8, r = 0, b = 0, l = 0)), # Más espacio sobre eje X para leyenda
               axis.text.y = element_text(hjust = 1, size = rel(1.5), colour = "black"),
               axis.text.x = element_text(size = rel(1.5), colour = "black"),
               axis.ticks.y = element_line(linewidth = 0.5, colour = "black"), # Ocultar ticks eje Y
               axis.ticks.x = element_line(linewidth = 0.5, colour = "black"),
               #legend.position = "bottom",           # Posición debajo del gráfico
               #legend.direction = "horizontal",     # Items en fila
               #legend.title.align = 0.5,            # Centrar título sobre items
               # Leyenda Dentro del Panel (Arriba-Derecha)
               legend.position = c(1.1, -0.05),       # Coordenadas relativas (X=0.95, Y=0.95)
               #legend.justification = c("right", "top"), # Anclar por la esquina superior derecha de la leyenda
               #legend.direction = "vertical",         # Disposición vertical
               legend.direction = "horizontal",
               legend.background = element_blank(), # Sin fondo
               legend.box.background = element_blank(), # Sin caja
               # Añadir espacio *encima* de la leyenda (entre eje X y leyenda)
               legend.box.margin = margin(t = 0, r= 0, b = 0, l = 0, unit="pt"),
               # Ajustar espaciado entre símbolo y texto si es necesario
               # legend.key.spacing.x = unit(0.1, 'cm'),
               legend.key.size = unit(0.4, "cm"), # Tamaño de los símbolos en leyenda
               legend.title = element_text(size=rel(1.25), face="bold"), # Estilo título leyenda
               legend.text = element_text(size=rel(1.25)), # Estilo texto leyenda

               # Margen del plot p_main
               plot.margin = ggplot2::margin(5, 5, 5, 5, unit = "pt")
               )
               #}
        #--- TEMA (Estilo Publicación Limpio con Leyenda Interior Arriba-Derecha) ---
        if (FALSE) {theme_classic(base_size = base_size, base_family = font_family) +
        theme(
               # Panel y Rejilla
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),

               # Ejes
               axis.line = element_line(colour = "black", linewidth = 0.4),
               # Margen superior del título del eje X puede ser más pequeño ya que la leyenda no está debajo
               axis.title.x = element_text(size = rel(1.1), margin = margin(t = 5, r = 0, b = 0, l = 0)),
               axis.text.y = element_text(hjust = 1, size = rel(2), colour = "black"),
               axis.text.x = element_text(size = rel(1.1), colour = "black"),
               axis.ticks.y = element_blank(linewidth = 0.5, colour = "black"), # Ocultar ticks eje Y
               axis.ticks.x = element_line(linewidth = 0.5, colour = "black"),

               # Leyenda Dentro del Panel (Arriba-Derecha)
               legend.position = c(0.95, 0.95),       # Coordenadas relativas (X=0.95, Y=0.95)
               legend.justification = c("right", "top"), # Anclar por la esquina superior derecha de la leyenda
               legend.direction = "vertical",         # Disposición vertical
               # legend.title.align = 0,              # Alinear título a la izquierda (opcional, default)
               legend.background = element_rect(fill = alpha("white", 0.7), colour = NA), # Fondo blanco semi-transparente, sin borde
               legend.box.background = element_blank(), # Sin caja exterior
               legend.box.margin = margin(t = 0, r= 0, b = 0, l = 0), # Sin márgenes extra alrededor de la leyenda
               legend.key.size = unit(0.4, "cm"),     # Tamaño de los símbolos
               legend.title = element_text(size=rel(0.9), face="bold"), # Estilo título leyenda
               legend.text = element_text(size=rel(0.9)),   # Estilo texto leyenda

               # Margen del plot p_main (puede mantenerse o ajustarse si es necesario)
               plot.margin = ggplot2::margin(5, 5, 5, 5, unit = "pt")
        )
        }
    # --- TABLA INTEGRADA (Usando ggplot) ---
    table_data <- plot_data_ordered %>%
                    filter(Model == "Adjusted") %>%
                    # Seleccionar datos relevantes para la tabla
                    select(compact_var_ordered, Estimate_CI_Formatted, P_Value_Formatted_With_Stars, Is_Significant) %>%
                    # Asegurar que el factor está ordenado igual que el plot
                    mutate(compact_var_ordered = factor(compact_var_ordered, levels = levels(plot_data_ordered$compact_var_ordered))) %>%
                    distinct() # Asegurar filas únicas por variable

            # --- TABLA INTEGRADA (Usando ggplot) ---
    table_data <- plot_data_ordered %>%
                    filter(Model == "Adjusted") %>%
                    # Seleccionar datos relevantes para la tabla
                    select(compact_var_ordered, Estimate_CI_Formatted, P_Value_Formatted_With_Stars, Is_Significant) %>%
                    # Asegurar que el factor está ordenado igual que el plot
                    mutate(compact_var_ordered = factor(compact_var_ordered, levels = levels(plot_data_ordered$compact_var_ordered))) %>%
                    distinct() # Asegurar filas únicas por variable

    # Plot de la tabla
    p_table <- ggplot(table_data, aes(y = compact_var_ordered)) +
                  # Columna Estimado [IC 95%]
                  geom_text(aes(x = -0.1, label = Estimate_CI_Formatted, fontface = Is_Significant),
                            hjust = 0, size = base_size * 0.4, family = font_family) + # Ajustar tamaño geom_text
                  # Columna P-valor (con estrellas)
                  geom_text(aes(x = 0.75, label = P_Value_Formatted_With_Stars, fontface = Is_Significant), # Ajustar posición X
                            hjust = 0, size = base_size * 0.4, family = font_family) +
                  # Títulos de columna (anotaciones)
                  # La posición Y (n_vars + 0.7) las coloca por encima de la última fila de datos (que está en n_vars)
                  annotate("text", x=-0.1, y=n_vars+0.7, label= paste0("Adj.", estimate_col_name," [95% CI]"), hjust=0, size=base_size*0.4, fontface="bold", family=font_family) +
                  annotate("text", x=0.75, y=n_vars+0.7, label= "P-value", hjust=0, size=base_size*0.4, fontface="bold", family=font_family) +
                  # Tema vacío
                  theme_void(base_family = font_family) +
                  theme(
                      # Margen izquierdo cero para alinear con plot principal
                      plot.margin = ggplot2::margin(t=5, r=5, b=5, l=0, unit = "pt")
                  ) +
                  # <<<--- CONTROL EXPLÍCITO DE LÍMITES Y ---<<<
                  # Eje X para dar espacio a las columnas de texto
                  # Eje Y definido explícitamente para asegurar espacio arriba de las anotaciones
                  coord_cartesian(
                      xlim = c(0, 1.1),  # Ajustar según el ancho del texto más largo
                      # Límites Y: desde 0.5 (debajo de la 1ra fila) hasta n_vars + 1.5 (bastante espacio sobre las anotaciones)
                      ylim = c(nx, n_vars + ny),
                      clip = "off" # Permitir dibujar fuera del panel si es necesario (seguridad)
                  )
                  # Nota: No se necesita scale_y_discrete aquí, coord_cartesian y el factor ordenado manejan el eje Y.

    # --- COMBINAR PLOT Y TABLA CON PATCHWORK ---
    # library(patchwork) # Debe estar cargado
    final_plot <- p_main + p_table +
                  # Ajustar anchos relativos (plot principal vs tabla)
                  plot_layout(widths = c(2.5, 1)) +
                  # Título general y subtítulo
                  plot_annotation(
                      title = plot_title,
                      subtitle = "Comparison of Unadjusted and Adjusted Estimates",
                      theme = theme(
                          plot.title = element_text(hjust = 0.5, size=rel(.15)*base_size, face="bold", family=font_family, margin = margin(b=2)),
                          plot.subtitle = element_text(hjust = 0.5, size=rel(.15)*base_size, colour = "grey30", family=font_family, margin = margin(t=0, r=0, b=10, l=0, unit = "pt")),
                          # Margen general para dar espacio al título
                          plot.margin = margin(t = 5, r = 10, b = -2, l = 5)
                      )
                  )

    # --- GUARDAR EL PLOT COMBINADO ---
    safe_filename <- sanitize_filename(filename_base);
    output_path_jpg <- file.path(output_analysis_dir, paste0(safe_filename, '.jpg'));
    output_path_pdf <- file.path(output_analysis_dir, paste0(safe_filename, '.pdf'))
    output_path_tiff <- file.path(output_analysis_dir, paste0(safe_filename, '.tif'))
    #output_path_svg <- file.path(output_analysis_dir, paste0(safe_filename, '.svg'))
    #output_path_emf <- file.path(output_analysis_dir, paste0(safe_filename, '.emf'))
    # Ancho total ajustado para acomodar tabla (pulgadas)
    total_plot_width = 7.5 # Ajustar según sea necesario

    # Guardar con alta calidad y dispositivo Cairo para PDF
    ggsave(output_path_jpg, plot=final_plot, width=total_plot_width, height=plot_height, dpi=400, device="jpeg", quality=100, limitsize=FALSE);
    cat(sprintf("  Saved Combined Forest Plot (JPG): %s\n", basename(output_path_jpg)))
    ggsave(output_path_pdf, plot=final_plot, width=total_plot_width, height=plot_height, device=cairo_pdf, limitsize=FALSE);
    cat(sprintf("  Saved Combined Forest Plot (PDF): %s\n", basename(output_path_pdf)))
    # <<<--- GUARDAR TIFF (CON ALTA RESOLUCIÓN Y COMPRESIÓN) ---<<<
    ggsave(
        output_path_tiff,
        plot = final_plot,
        width = total_plot_width,
        height = plot_height,
        device = "tiff",           # Especificar dispositivo TIFF
        dpi = 600,                 # ¡MUY IMPORTANTE! Alta resolución
        compression = "lzw",       # Compresión sin pérdidas (recomendado)
        limitsize = FALSE
    )
    cat(sprintf("  Saved Combined Forest Plot (TIFF @ %d DPI): %s\n", 600, basename(output_path_tiff)))
    #ggsave(output_path_svg, plot=final_plot, width=total_plot_width, height=plot_height, device=svglite::svglite, limitsize=FALSE)
    #cat(sprintf("  Saved Combined Forest Plot (SVG): %s\n", basename(output_path_svg)))
    #ggsave(output_path_emf, plot=final_plot, width=total_plot_width, height=plot_height, device = "emf", limitsize=FALSE)
    #cat(sprintf("  Saved Combined Forest Plot (EMF): %s\n", basename(output_path_emf)))
    return(final_plot) # Devolver el objeto patchwork combinado
}


# --- Generar Gráficos Forest Combinados Mejorados ---
# (Llamadas a la función - sin cambios respecto a tu versión funcional)

# ORs
if (nrow(univariate_or_df_pooled) > 0 && nrow(multivariate_or_df_pooled) > 0) {
    create_enhanced_combined_forest_plot(
        unadj_df = univariate_or_df_pooled,
        adj_df = multivariate_or_df_pooled, # Usar multivariado principal
        estimate_col_name = "OR",
        plot_title = "Association with MacroB12 Status (Odds Ratios)",
        filename_base = "Fig_ForestPlot_Combined_ORs_PubStyle") # Nuevo nombre archivo
} else { cat("Skipping Combined OR Forest Plot (missing results).\n") }

# PRs (si se calcularon)
if(calculate_pr) {
    if (nrow(univariate_pr_df_pooled) > 0 && nrow(multivariate_pr_df_pooled) > 0) {
        create_enhanced_combined_forest_plot(
            unadj_df = univariate_pr_df_pooled,
            adj_df = multivariate_pr_df_pooled, # Usar multivariado principal
            estimate_col_name = "PR",
            plot_title = "Association with MacroB12 Status (Prevalence Ratios)",
            filename_base = "Fig_ForestPlot_Combined_PRs_PubStyle") # Nuevo nombre archivo
    } else { cat("Skipping Combined PR Forest Plot (missing results).\n") }
}


# =======================================================================
# 6. GENERAR LEYENDA PARA ETIQUETAS COMPACTAS
# =======================================================================
# (Código sin cambios)
cat("\n--- 6. Generating Legend for Compact Labels ---\n")
plot_labels_df <- bind_rows( multivariate_or_df_pooled %>% select(Variable_Raw, Variable, Compact_Var_Label), multivariate_pr_df_pooled %>% select(Variable_Raw, Variable, Compact_Var_Label) ) %>% distinct(Variable_Raw, .keep_all = TRUE) %>% filter(!is.na(Compact_Var_Label) & !is.na(Variable))
if(nrow(plot_labels_df) > 0) { plot_labels_df <- plot_labels_df %>% arrange(Compact_Var_Label); legend_title_text <- "Variable Abbreviation Legend"; legend_separator <- paste(rep("-", nchar(legend_title_text)), collapse=""); legend_items_text <- paste(sprintf("%-18s = %s", plot_labels_df$Compact_Var_Label, plot_labels_df$Variable), collapse = "\n"); legend_text <- paste(legend_title_text, legend_separator, legend_items_text, sep="\n"); legend_file_path <- file.path(output_analysis_dir, "Plot_Compact_Label_Legend.txt"); tryCatch({ writeLines(legend_text, legend_file_path); cat(sprintf("Compact label legend saved to: %s\n", normalizePath(legend_file_path))) }, error = function(e){ cat(sprintf("ERROR saving label legend: %s\n", e$message)) }); cat("\n"); cat(legend_text); cat("\n") } else { cat("Could not generate compact label legend.\n") }


# =======================================================================
# 7. MENSAJE FINAL
# =======================================================================
# (Código sin cambios)
cat("\n--- Combined Forest Plot Script Completed ---")
cat(sprintf("\nAnalysis performed using pooled results from MICE object (m=%d imputations).", mids_data_original$m))
cat(sprintf("\nAll outputs saved in '%s'.\n", normalizePath(output_analysis_dir)))
cat(sprintf("Script finished at: %s\n", Sys.time()))