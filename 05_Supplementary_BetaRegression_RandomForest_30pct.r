# ====================================================================
#        Script R for Beta Regression on Recovery_pct
#        (MICE Pooled, Pub-Style Plot, Diagnostics, Boxed Axes)
# ====================================================================
# 0. CONFIGURACIÓN E IMPORTACIÓN
# ====================================================================
# --- Cargar Paquetes ---
required_packages <- c("data.table", "tidyverse", "ggplot2", "here", "rstudioapi",
                       "mice",      # Para manejar objeto mids y pool
                       "betareg",   # <<<--- Para Beta Regression
                       "performance", # <<<--- Para Pseudo R-squared
                       "car",       # <<<--- Para VIF
                       "lmtest",    # Para coeftest (si se hiciera manual)
                       "ggthemes", "RColorBrewer", "ggrepel",
                       "knitr", "kableExtra", "skimr", "ggpubr", "paletteer",
                       "dplyr", "stats", "forcats", "patchwork", # Para combinar plots
                       "svglite")   # Para guardar SVG

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing package: %s\n", pkg))
    tryCatch(install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org/"),
             error = function(e){ cat(sprintf("Failed to install %s: %s\n", pkg, e$message))})
  }
}
sapply(required_packages, install_if_missing)

library(data.table); library(tidyverse); library(ggplot2); library(here); library(rstudioapi)
library(mice); library(betareg); library(performance); library(car); library(lmtest);
library(ggthemes); library(RColorBrewer); library(ggrepel);
library(knitr); library(kableExtra); library(skimr); library(ggpubr); library(paletteer)
library(dplyr); library(stats); library(forcats); library(patchwork); library(svglite)

# --- !! ESTABLECER DIRECTORIO DE TRABAJO MANUALMENTE !! ---
# Asegúrate de que esta ruta sea EXACTAMENTE la carpeta donde está este script
tryCatch({
  setwd("F:/WorkInProgre/Collaborators/Chema/Data/R/DataClean")
  cat("Manually set working directory to:", getwd(), "\n")
}, error = function(e) {
  warning("Could not manually set working directory. Paths might be incorrect.", call. = FALSE)
})
# -----------------------------------------------------------

# --- Directorios ---
script_dir <- tryCatch({ dirname(rstudioapi::getSourceEditorContext()$path) }, error = function(e) { getwd() })
input_data_prep_dir <- file.path(script_dir, "R_Outputs_Data_Prep_New_30pct")
output_analysis_dir <- file.path(script_dir, "Supplementary_BetaRegression_30pct") # <<<--- Directorio actualizado
if (!dir.exists(output_analysis_dir)) { dir.create(output_analysis_dir) } else { cat("Output directory exists:", normalizePath(output_analysis_dir), "\n") }


# --- Constantes y Estilos Globales ---
RANDOM_STATE <- 42; set.seed(RANDOM_STATE)
FONT_FAMILY <- "Arial"
BASE_SIZE <- 11 # Tamaño base aumentado
HIGH_RES_DPI <- 600

# Colores consistentes
palette_calc <- paletteer_d("ggthemes::calc", n = 5)[1:4]
color_unadjusted <- palette_calc[3]
color_adjusted <- palette_calc[4]

# --- Cargar Mapeo de Etiquetas COMPLETAS (Keys = Nombres INTERNOS/make.names) ---
# (Usar el LABEL_MAP completo y verificado de scripts anteriores)
LABEL_MAP <- list(
  'EDAD' = 'Age (years)', 'sexo' = 'Sex', 'procedencia.' = 'Source',
  'VB12.PRE.PEG' = 'VB12 Pre-PEG (pg/mL)', 'vb12_post_peg' = 'VB12 Post-PEG (pg/mL)',
  'Recovery_pct' = 'PEG Recovery (%)', 'MacroB12_Positive' = 'MacroB12 Status',
  'pcr' = 'CRP (mg/L)', 'fol' = 'Serum Folate (ng/mL)', 'hb' = 'Hemoglobin (g/dL)',
  'vcm' = 'MCV (fL)', 'ade' = 'RDW (%)', 'Diag_Hematologic' = 'Diagnosis: Hematologic',
  'Diag_Hepatic' = 'Diagnosis: Hepatic', 'Diag_Oncologic' = 'Diagnosis: Oncologic',
  'Diag_Rheum_Autoimmune' = 'Diagnosis: Rheum./Autoimmune', 'Diag_Infectious' = 'Diagnosis: Infectious',
  'Diag_Renal' = 'Diagnosis: Renal', 'Diag_GI_Suspicion' = 'Diagnosis: GI/Neo. Suspicion',
  'Diag_Hypertension' = 'Diagnosis: Hypertension', 'Diag_Diabetes' = 'Diagnosis: Diabetes',
  'Diag_Healthy_Unspecified' = 'Diagnosis: Healthy/Unspecified', 'Male' = 'Male', 'Female' = 'Female',
  'Ambulatory' = 'Ambulatory', 'Hospital' = 'Hospital', 'External Center' = 'External Center',
  'Negative'='Negative', 'Positive'='Positive',
  # Términos de Modelo (Claves EXACTAS)
  'SEXOMale'='Sex: Male vs Female',
  'ProcedenciaHospital'='Source: Hospital vs Ambulatory',
  'ProcedenciaExternal Center'='Source: Ext. Ctr. vs Ambulatory',
  'Diag_Hematologic1'='Diagnosis: Hematologic', # No añadir "Present vs Absent" aquí si es para etiqueta de plot
  'Diag_Hepatic1'='Diagnosis: Hepatic',
  'Diag_Oncologic1'='Diagnosis: Oncologic',
  'Diag_Rheum_Autoimmune1'='Diagnosis: Rheum./Autoimmune',
  'Diag_Infectious1'='Diagnosis: Infectious',
  'Diag_Renal1'='Diagnosis: Renal',
  'Diag_GI_Suspicion1'='Diagnosis: GI/Neo. Suspicion',
  'Diag_Hypertension1'='Diagnosis: Hypertension',
  'Diag_Diabetes1'='Diagnosis: Diabetes',
  'Diag_Healthy_Unspecified1'='Diagnosis: Healthy/Unspecified'
)

# --- Mapeo de Etiquetas COMPACTAS (Revisado) ---
COMPACT_LABEL_MAP <- list(
  'EDAD' = 'Age', 'VB12.PRE.PEG' = 'VB12 Pre', 'PCR' = 'CRP', 'FOL' = 'Folate',
  'HB' = 'Hb', 'VCM' = 'MCV', 'ADE' = 'RDW',
  'SEXOMale'='Sex (M vs F)',
  'ProcedenciaHospital'='Source (Hosp)', # Más corto
  'ProcedenciaExternal Center'='Source (Ext Ctr)', # Más corto
  'Diag_Hematologic1'='Diag: Hematol.', # Más corto
  'Diag_Hepatic1'='Diag: Hepatic',
  'Diag_Oncologic1'='Diag: Oncologic',
  'Diag_Rheum_Autoimmune1'='Diag: Rheum/AI',
  'Diag_Infectious1'='Diag: Infect.', # Más corto
  'Diag_Renal1'='Diag: Renal',
  'Diag_GI_Suspicion1'='Diag: GI/Neo',
  'Diag_Hypertension1'='Diag: Hypertens.', # Más corto
  'Diag_Diabetes1'='Diag: Diabetes',
  'Diag_Healthy_Unspecified1'='Diag: Healthy'
)

# --- Orden de Variables para el Plot (Basado en Forest Plot anterior) ---
# Usar los NOMBRES DE TÉRMINO DEL MODELO
PLOT_ORDER_VARIABLES <- c(
    "EDAD", "SEXOMale", "ProcedenciaHospital", "ProcedenciaExternal Center",
    "Diag_Healthy_Unspecified1", "Diag_Hematologic1", "Diag_Hepatic1",
    "Diag_Oncologic1", "Diag_Rheum_Autoimmune1", "Diag_Infectious1",
    "Diag_Renal1", "Diag_GI_Suspicion1", "Diag_Hypertension1", "Diag_Diabetes1",
    "VB12.PRE.PEG", "PCR", "FOL", "HB", "VCM", "ADE"
)

# --- Cargar mapeo de nombres finales (si existe) ---
# (Se mantiene igual)
final_name_map_file <- file.path(input_data_prep_dir, "column_rename_map_mk_to_eng.rds")
# ... (código para cargar mapa se mantiene) ...

# --- Helper Functions ---
# (get_label, get_compact_label, format_p_value, add_stars, sanitize_filename se mantienen igual)
# ... (definiciones de helper functions) ...
get_label <- function(name, default_fmt=TRUE) { name_str <- as.character(name); internal_name <- name_str ; if(exists("col_names_type") && col_names_type == "final"){ internal_lookup <- final_to_internal_map[[name_str]]; if(!is.null(internal_lookup)) internal_name <- internal_lookup }; label <- LABEL_MAP[[internal_name]]; if (!is.null(label)) return(label); clean_name <- gsub("`", "", name_str); if (default_fmt) return(tools::toTitleCase(gsub("[._]", " ", clean_name))) else return(clean_name) }
get_compact_label <- function(name) { name_str <- as.character(name); internal_name <- name_str; if(exists("col_names_type") && col_names_type == "final"){ internal_lookup <- final_to_internal_map[[name_str]]; if(!is.null(internal_lookup)) internal_name <- internal_lookup }; compact_label <- COMPACT_LABEL_MAP[[internal_name]]; if (!is.null(compact_label)) { return(compact_label) } else { warning(paste("No compact label found for internal name:", internal_name, "(original term:", name_str, "). Using Variable_Raw instead.")); return(name_str) } }
format_p_value <- function(p_vector) { sapply(p_vector, function(p) { if (is.na(p)) { return("NA") } else if (p < 0.0001) { return(sprintf("%.1e", p)) } else if (p < 0.001) { return("<0.001") } else if (p < 0.01) { return(sprintf("%.3f", p)) } else { return(sprintf("%.2f", p)) } }) }
add_stars <- function(p_vector) { sapply(p_vector, function(p) { if (is.na(p) || !is.numeric(p)) return('') ; if (p < 0.001) return('***') ; if (p < 0.01) return('**') ; if (p < 0.05) return('*') ; return('') }) }
sanitize_filename <- function(name) { name_sanitized <- gsub('[<>:"/\\|?*%[:cntrl:]]', '_', as.character(name)); name_sanitized <- gsub('[[:punct:]]', '_', name_sanitized); name_sanitized <- gsub('\\s+', '_', name_sanitized); name_sanitized <- gsub('_+', '_', name_sanitized); name_sanitized <- gsub('^_|_$', '', name_sanitized); name_sanitized <- substr(name_sanitized, 1, 100); return(name_sanitized) }

# =======================================================================
# 1. CARGAR DATOS MICE Y PREPARAR PARA BETA REGRESIÓN
# =======================================================================
cat("\n--- 1. Loading MICE Data Object & Preparing for Beta Regression ---\n")
mice_rds_file <- file.path(input_data_prep_dir, "mice_imputation_object_internalNames.rds")
expected_path <- normalizePath(mice_rds_file, mustWork = FALSE)
cat("DEBUG: Script is looking for the MICE RDS file at this exact path:\n", expected_path, "\n")
mids_data_original <- NULL
if (!file.exists(expected_path)) {
    stop(paste("DEBUG: The file was NOT found at the expected path:", expected_path))
} else {
    cat("DEBUG: File DOES exist at the expected path. Proceeding to load...\n")
}
if (file.exists(mice_rds_file)) { tryCatch({ mids_data_original <- readRDS(mice_rds_file); cat("Successfully loaded MICE object (.rds) with", mids_data_original$m, "imputations.\n") }, error = function(e) { cat("Error loading MICE object:", e$message, "\n") }) }; if (is.null(mids_data_original)) { stop("Could not load MICE object.") }

mids_data_reg <- mids_data_original
outcome_var <- make.names("Recovery_pct") # Nombre original 0-100
outcome_var_beta <- "Recovery_beta"      # Nuevo nombre para variable transformada (0,1)

if (!outcome_var %in% names(mids_data_reg$data)) { stop(paste("Outcome variable '", outcome_var, "' not found.")) }

# Transformar outcome a (0,1) DENTRO de cada imputación
# Usaremos la transformación (y*(n-1)+0.5)/n donde y es la proporción
n_total <- nrow(mids_data_reg$data) # Tamaño muestral original
cat(sprintf("Transforming %s to (0,1) scale for Beta Regression using n=%d...\n", outcome_var, n_total))

# Añadir la variable transformada al objeto mids
# Nota: Esto modifica el objeto mids_data_reg in-place
tryCatch({
  mids_data_reg <- mice::complete(mids_data_reg, action = "long", include = TRUE) %>%
    mutate(
      # Primero, limitar Recovery_pct a un máximo de 100
      Recovery_pct_capped = ifelse(!!sym(outcome_var) > 100, 100, !!sym(outcome_var)),
      # Luego, transformar el valor limitado/truncado
      !!sym(outcome_var_beta) := (Recovery_pct_capped / 100 * (n_total - 1) + 0.5) / n_total
    ) %>%
    # Opcional: eliminar la columna intermedia si no la necesitas
    select(-Recovery_pct_capped) %>%
    as.mids()
   cat("Successfully added transformed outcome variable (capped at 100%):", outcome_var_beta, "\n")
}, error = function(e){
   stop(paste("Error transforming outcome variable within MIDS object:", e$message))
})

# <<<--- NUEVO BLOQUE DE VERIFICACIÓN ---<<<
cat("Verifying transformed outcome range across all imputations...\n")
all_valid <- TRUE  
for (i in 1:mids_data_reg$m) {
    imp_data <- mice::complete(mids_data_reg, action = i)
    current_outcome_vals <- imp_data[[outcome_var_beta]]
    range_vals <- range(current_outcome_vals, na.rm = TRUE)
    num_invalid <- sum(current_outcome_vals <= 0 | current_outcome_vals >= 1 | is.na(current_outcome_vals), na.rm = TRUE)

    cat(sprintf("  Imputation %d: Range = [%.5f, %.5f], N Invalid (<=0 or >=1 or NA) = %d\n",
                i, range_vals[1], range_vals[2], num_invalid))

    if (num_invalid > 0) {
        all_valid <- FALSE
        cat(sprintf("    Problematic values in imp %d:\n", i))
        print(summary(current_outcome_vals[current_outcome_vals <= 0 | current_outcome_vals >= 1 | is.na(current_outcome_vals)]))
    }
}

if (!all_valid) {
    stop("Transformation failed: Found values <= 0 or >= 1 or NA in the Recovery_beta variable. Beta regression cannot proceed.")
} else {
    cat("Transformation successful: All Recovery_beta values are strictly between 0 and 1.\n")
}
# <<<--- FIN NUEVO BLOQUE DE VERIFICACIÓN ---<<<
# Verificar transformación en datos originales
summary(mids_data_reg$data[[outcome_var_beta]])

# --- Preparar Predictores ---
excluded_vars <- c(make.names("MacroB12_Positive"), make.names("vb12_post_peg"),
                   outcome_var, outcome_var_beta); # Excluir original y transformada, y otras no relevantes
predictor_vars_internal <- names(mids_data_reg$data); predictor_vars_internal <- setdiff(predictor_vars_internal, excluded_vars)
cat("Verifying/setting factor types...\n"); # Asegurar factores después de la transformación
sex_col_internal <- make.names("sexo"); proc_col_internal <- make.names("procedencia.")
if(sex_col_internal %in% names(mids_data_reg$data)) { if(!is.factor(mids_data_reg$data[[sex_col_internal]])) {mids_data_reg <- mice::complete(mids_data_reg, action="long", include=TRUE) %>% mutate(!!sym(sex_col_internal) := factor(!!sym(sex_col_internal))) %>% as.mids()}; if("Female" %in% levels(mids_data_reg$data[[sex_col_internal]])) { mids_data_reg <- mice::complete(mids_data_reg, action="long", include=TRUE) %>% mutate(!!sym(sex_col_internal) := relevel(factor(!!sym(sex_col_internal)), ref="Female")) %>% as.mids() } else {warning("Ref 'Female' not found for sexo")} }
if(proc_col_internal %in% names(mids_data_reg$data)) { if(!is.factor(mids_data_reg$data[[proc_col_internal]])) {mids_data_reg <- mice::complete(mids_data_reg, action="long", include=TRUE) %>% mutate(!!sym(proc_col_internal) := factor(!!sym(proc_col_internal), levels = c("Ambulatory", "External Center", "Hospital"))) %>% as.mids()}; if("Ambulatory" %in% levels(mids_data_reg$data[[proc_col_internal]])) { mids_data_reg <- mice::complete(mids_data_reg, action="long", include=TRUE) %>% mutate(!!sym(proc_col_internal) := relevel(factor(!!sym(proc_col_internal), levels = c("Ambulatory", "External Center", "Hospital")), ref="Ambulatory")) %>% as.mids() } else {warning("Ref 'Ambulatory' not found for procedencia.")} }
diag_cols_internal <- names(mids_data_reg$data)[startsWith(names(mids_data_reg$data), "Diag_") & !grepl("1$", names(mids_data_reg$data))]; for(col in diag_cols_internal){ if(col %in% names(mids_data_reg$data)){ if(!is.factor(mids_data_reg$data[[col]])) {mids_data_reg <- mice::complete(mids_data_reg, action="long", include=TRUE) %>% mutate(!!sym(col) := factor(!!sym(col), levels=c(0,1))) %>% as.mids()} } }; cat("Factor types verified/set.\n")
constant_vars <- names(which(sapply(complete(mids_data_reg, 1)[, predictor_vars_internal], function(x) n_distinct(x, na.rm = TRUE)) <= 1)); if(length(constant_vars) > 0) { cat("Removing constant predictors:", paste(constant_vars, collapse=", "), "\n"); predictor_vars_internal <- setdiff(predictor_vars_internal, constant_vars) }; predictor_vars_model <- predictor_vars_internal
cat("Outcome for Beta Regression:", outcome_var_beta, "\nPredictors:\n", paste(predictor_vars_model, collapse=", "), "\n")
col_names_type <- "internal" # Para helper functions

# =======================================================================
# 2. EJECUTAR BETA REGRESIÓN Y OBTENER RESULTADOS POOLED
# =======================================================================
cat("\n--- 2. Running Pooled Beta Regression (Univariate & Multivariate) ---\n")

# --- Función run_and_pool_betareg (Intenta usar pool estándar) ---
run_and_pool_betareg <- function(formula_str, data_mids) {
    results_list <- list() ; fit_pooled <- NULL; summary_pooled <- NULL
    model_fits <- NULL # Guardar los modelos ajustados
    tryCatch({
        # Ajustar el modelo en cada imputación
        fit_models <- with(data = data_mids, expr = betareg::betareg(formula = as.formula(formula_str), link = "logit")) # Logit link por defecto
        if(is.null(fit_models) || !inherits(fit_models, "mira")) { stop("'with' did not return 'mira' object.") }
        if(length(fit_models$analyses) == 0) { stop("No analyses in 'mira' object.") }
        # Verificar errores individuales
        errors_in_analyses <- sapply(fit_models$analyses, function(x) inherits(x, "try-error") || is.null(coef(x)))
        if(any(errors_in_analyses)) { cat("  WARNING: Errors or NULL coefs in individual BetaReg fits. N errors:", sum(errors_in_analyses), "\n") }
        valid_analyses <- fit_models$analyses[!errors_in_analyses]
        if (length(valid_analyses) < 2) { stop("Fewer than 2 valid model fits after checking errors.") }
        # Intentar pooling estándar
        fit_pooled <- suppressWarnings(pool(fit_models)) # Suprimir warnings si pool se queja de la clase
        summary_pooled <- summary(fit_pooled, conf.int = TRUE)
        # Extraer resultados
        if (!is.null(summary_pooled) && nrow(summary_pooled) > 0) {
            for (i in 1:nrow(summary_pooled)) {
                term_row <- summary_pooled[i, ]
                term_name <- term_row$term
                # Excluir intercepto de media y precisión (phi)
                if (grepl("Intercept", term_name) || grepl("phi", term_name)) next
                est <- term_row$estimate # Coeficiente en escala logit
                ci_l <- term_row$conf.low
                ci_h <- term_row$conf.high
                p_val <- term_row$p.value
                if(all(is.finite(c(est, ci_l, ci_h, p_val)))){
                    results_list[[length(results_list) + 1]] <- list(
                        Variable_Raw=term_name, Variable=get_label(term_name),
                        Compact_Var_Label=get_compact_label(term_name),
                        Coefficient=est, CI_Low=ci_l, CI_High=ci_h, P_Value=p_val
                    )
                } else { cat(sprintf("    Skipping term '%s' due to non-finite results.\n", term_name)) }
            }
            # Guardar modelos ajustados válidos para diagnósticos/R2
             attr(results_list, "valid_model_fits") <- valid_analyses
        } else { cat(sprintf("    Pooling resulted in NULL/empty summary: %s\n", formula_str)) }
    }, error = function(e) {
        cat(sprintf("  ERROR pooling BetaReg for %s: %s\n", formula_str, e$message))
        attr(results_list, "valid_model_fits") <- NULL
    }, warning = function(w) {
        cat(sprintf("  WARNING pooling BetaReg for %s: %s\n", formula_str, w.message))
        # Podríamos intentar continuar si es solo un warning
    })
    results_df <- if (length(results_list) > 0) bind_rows(results_list) else data.frame(Variable_Raw=character(), Variable=character(), Compact_Var_Label=character(), Coefficient=numeric(), CI_Low=numeric(), CI_High=numeric(), P_Value=numeric(), stringsAsFactors = FALSE)
    # Adjuntar modelos válidos al dataframe resultante
    attr(results_df, "valid_model_fits") <- attr(results_list, "valid_model_fits")
    return(results_df)
}

# --- Ejecutar Modelos Univariados ---
cat("Running Univariate Beta Regression Models...\n")
formula_lhs_beta <- paste0("`", outcome_var_beta, "`")
univariate_beta_results_list <- lapply(predictor_vars_model, function(pred) {
    formula_uni_str <- paste0(formula_lhs_beta, " ~ `", pred, "`")
    run_and_pool_betareg(formula_uni_str, mids_data_reg)
})
univariate_beta_df_pooled <- bind_rows(univariate_beta_results_list) %>% filter(!is.na(Coefficient))

# --- Ejecutar Modelo Multivariado ---
cat("Running Multivariate Beta Regression Model...\n")
multivariate_beta_df_pooled <- data.frame()
valid_multi_fits <- NULL # Para guardar los modelos válidos

if (length(predictor_vars_model) > 1) {
    formula_beta_multi_str <- paste0(formula_lhs_beta, " ~ ", paste(paste0("`", predictor_vars_model, "`"), collapse = " + "))
    cat("Formula:", formula_beta_multi_str, "\n")
    multi_results <- run_and_pool_betareg(formula_beta_multi_str, mids_data_reg)
    multivariate_beta_df_pooled <- multi_results %>% filter(!is.na(Coefficient))
    valid_multi_fits <- attr(multi_results, "valid_model_fits") # Recuperar los modelos
} else {
    cat("Not enough predictors for multivariate model.\n")
}

# --- Calcular Pseudo R-squared Aproximado ---
if (!is.null(valid_multi_fits) && length(valid_multi_fits) > 0) {
    cat("Calculating approximate pooled Pseudo R-squared (Nagelkerke)...\n")
    r2_values <- sapply(valid_multi_fits, function(fit) {
        tryCatch({ performance::r2_nagelkerke(fit) }, error = function(e) NA)
    })
    r2_values <- r2_values[!is.na(r2_values)]
    if (length(r2_values) > 0) {
        pooled_r2_approx <- mean(r2_values)
        cat(sprintf("  Approximate Pooled Pseudo R-squared (Nagelkerke): %.3f (based on %d valid models)\n",
                    pooled_r2_approx, length(r2_values)))
    } else {
        cat("  Could not calculate Pseudo R-squared.\n")
    }
}

# --- Mostrar y Guardar Resultados ---
cat("\n--- Pooled Multivariate Beta Regression Results ---\n")
if(nrow(multivariate_beta_df_pooled) > 0) {
    beta_results_display <- multivariate_beta_df_pooled %>%
        mutate(Estimate_CI = sprintf("%.3f (%.3f, %.3f)", Coefficient, CI_Low, CI_High),
               P_Value_Formatted = format_p_value(P_Value),
               Signif = add_stars(P_Value)) %>%
        select(Variable = Compact_Var_Label, Estimate_CI, P_Value_Formatted, Signif)
    print(kable(beta_results_display, caption = "Pooled Multivariate Beta Regression Coefficients (Logit Link)") %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE))
} else { cat("No valid results obtained from the multivariate beta regression model.\n") }

# Guardar tablas
if(nrow(univariate_beta_df_pooled) > 0) { beta_table_uni_path <- file.path(output_analysis_dir, "Table_BetaReg_Univariate_Pooled.csv"); write.csv(univariate_beta_df_pooled, beta_table_uni_path, row.names=F, na=""); cat(sprintf("\nUnivariate results saved to: %s\n", normalizePath(beta_table_uni_path))) }
if(nrow(multivariate_beta_df_pooled) > 0) { beta_table_multi_path <- file.path(output_analysis_dir, "Table_BetaReg_Multivariate_Pooled.csv"); write.csv(multivariate_beta_df_pooled, beta_table_multi_path, row.names=F, na=""); cat(sprintf("Multivariate results saved to: %s\n", normalizePath(beta_table_multi_path))) }


# =======================================================================
# 3. GENERAR COEFFICIENT PLOT COMBINADO (Beta Regression)
# =======================================================================
cat("\n--- 3. Generating Combined Coefficient Plot (Beta Regression) ---\n")

# --- Tema ggplot base para publicación con CAJA CERRADA ---
theme_pub_boxed <- theme_bw(base_size = BASE_SIZE, base_family = FONT_FAMILY) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1.3), margin = margin(b=4)),
    plot.subtitle = element_text(hjust = 0.5, size = rel(1.3), colour = "grey30", margin = margin(b=15)),
    axis.title = element_text(face = "plain", size = rel(1.3)),
    axis.text = element_text(size = rel(1.3), colour = "black"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
    axis.line.x.bottom = element_line(colour = "black", linewidth = 0.5),
    axis.line.x.top    = element_line(colour = "black", linewidth = 0.5),
    axis.line.y.left   = element_line(colour = "black", linewidth = 0.5),
    axis.line.y.right  = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.x.top = element_blank(), axis.ticks.y.right = element_blank(),
    axis.ticks.x.bottom = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.y.left = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    legend.background = element_blank(), legend.box.background = element_blank(), legend.key = element_blank(),
    legend.title = element_text(size=rel(1.2), face="bold"),
    legend.text = element_text(size=rel(1.2)),
    plot.margin = margin(10, 15, 10, 10)
  )

# --- Función Coefficient Plot Combinado con Tabla Integrada ---
create_combined_beta_coefficient_plot <- function(unadj_df, adj_df, plot_title, filename_base, order_vars) {

    required_cols <- c("Variable_Raw", "Variable", "Compact_Var_Label", "Coefficient", "CI_Low", "CI_High", "P_Value")
    if (!all(required_cols %in% names(unadj_df)) || !all(required_cols %in% names(adj_df))) { stop("Input DFs missing required columns for plot.") }

    unadj_df_mod <- unadj_df %>% mutate(Model = "Unadjusted")
    adj_df_mod <- adj_df %>% mutate(Model = "Adjusted")
    combined_df <- bind_rows(unadj_df_mod, adj_df_mod)

    # Fallback etiqueta compacta
    if (!"Compact_Var_Label" %in% names(combined_df)) { combined_df <- combined_df %>% mutate(Compact_Var_Label = Variable_Raw)}
    combined_df <- combined_df %>% mutate(Compact_Var_Label = ifelse(is.na(Compact_Var_Label) | Compact_Var_Label == "", Variable_Raw, Compact_Var_Label))

    # Añadir texto y filtrar
    combined_df <- combined_df %>%
        mutate(
            p_formatted = format_p_value(P_Value),
            stars = add_stars(P_Value),
            # Columnas para tabla
            Estimate_Formatted = sprintf("%.2f", Coefficient),
            #CI_Formatted = sprintf("%.2f, %.2f", CI_Low, CI_High),
            #CI_Formatted = sprintf("%.2f – %.2f", CI_Low, CI_High), # Opción 2: En dash (U+2013) - verificar codificación
            CI_Formatted = sprintf("%.2f - %.2f", CI_Low, CI_High), # Opción 3: Guion normal (más seguro)
            P_Value_Formatted_With_Stars = paste0(p_formatted, stars),
            Is_Significant = ifelse(Model == "Adjusted" & !is.na(P_Value) & P_Value < 0.05, "bold", "plain")
         ) %>%
        filter(!is.na(Coefficient) & is.finite(Coefficient) & !is.na(CI_Low) & is.finite(CI_Low) & !is.na(CI_High) & is.finite(CI_High))

    if(nrow(combined_df) == 0) { cat(sprintf(" Skipping plot '%s': No valid data after filtering.\n", plot_title)); return(NULL) }

    # --- Ordenar según PLOT_ORDER_VARIABLES ---
    available_vars_in_order <- intersect(order_vars, unique(combined_df$Variable_Raw))
    if(length(available_vars_in_order) == 0) { cat(sprintf(" Skipping plot '%s': No variables found matching order_vars.\n", plot_title)); return(NULL) }

    combined_df_ordered <- combined_df %>%
        filter(Variable_Raw %in% available_vars_in_order) %>%
        mutate(
            compact_var_ordered = factor(Variable_Raw,
                                         levels = rev(available_vars_in_order),
                                         labels = rev(sapply(available_vars_in_order, get_compact_label))),
            Model = factor(Model, levels = c("Unadjusted", "Adjusted"))
        ) %>% filter(!is.na(compact_var_ordered))

    if(nrow(combined_df_ordered) == 0 || nlevels(combined_df_ordered$compact_var_ordered) == 0) { cat(sprintf(" Skipping plot '%s': No data left after ordering.\n", plot_title)); return(NULL) }

    n_vars <- nlevels(combined_df_ordered$compact_var_ordered)
    plot_height <- max(6, n_vars * 0.4 + 3); dodge_width <- 0.6
nx=1
ny=0 #Forzar límites Y para alineación con p_table
    # --- Plot Principal (Puntos y Líneas) ---
    p_main <- ggplot(combined_df_ordered, aes(x = Coefficient, y = compact_var_ordered, color = Model)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
        geom_errorbarh(aes(xmin = CI_Low, xmax = CI_High),
                       height = 0, linewidth = 0.725, position = position_dodge(width = dodge_width)) +
        geom_point(shape = 16, size = 4, position = position_dodge(width = dodge_width)) +
        scale_color_manual(values = c("Unadjusted" = color_unadjusted, "Adjusted" = color_adjusted), , name = "Model:") +
        labs(x = "Coefficient Estimate (Logit Link)", y = NULL) +
        theme_pub_boxed + # Aplicar tema caja cerrada
        theme(legend.position = c(1.25, -0.0425), # Mover leyenda arriba
              #legend.justification = "right",
              legend.direction = "horizontal",
              plot.margin = margin(5, 5, 5, 5)) +
         # Forzar límites Y para alineación con p_table
        coord_cartesian(
            ylim = c(nx, n_vars + ny),
            clip = "off"
        ) 

    # --- Tabla Integrada ---
    table_data <- combined_df_ordered %>%
                    filter(Model == "Adjusted") %>%
                    select(compact_var_ordered, Estimate_Formatted, CI_Formatted, P_Value_Formatted_With_Stars, Is_Significant) %>%
                    distinct()

    text_size_table <- BASE_SIZE * 0.33 # Tamaño texto tabla

    p_table <- ggplot(table_data, aes(y = compact_var_ordered)) +
                  geom_text(aes(x = -.075, label = Estimate_Formatted, fontface = Is_Significant), hjust = 0.5, size = text_size_table*1.2, family = FONT_FAMILY) +
                  geom_text(aes(x = 0.35, label = CI_Formatted, fontface = Is_Significant), hjust = 0.5, size = text_size_table*1.2, family = FONT_FAMILY) +
                  geom_text(aes(x = 0.95, label = P_Value_Formatted_With_Stars, fontface = Is_Significant), hjust = 1, size = text_size_table*1.2, family = FONT_FAMILY) +
                  # Títulos de columna
                  annotate("text", x=-.075, y=n_vars+0.7, label= "Adj. Coeff.", hjust=0.5, size=text_size_table*1.2, fontface="bold", family=FONT_FAMILY) +
                  annotate("text", x=0.35, y=n_vars+0.7, label= "95% CI", hjust=0.5, size=text_size_table*1.2, fontface="bold", family=FONT_FAMILY) +
                  annotate("text", x=0.95, y=n_vars+0.7, label= "P-value", hjust=1, size=text_size_table*1.2, fontface="bold", family=FONT_FAMILY) +
                  theme_void(base_family = FONT_FAMILY) +
                  theme(plot.margin = margin(t=5, r=0, b=5, l=0)) +
                  coord_cartesian(xlim = c(-0.1, 1.0), ylim = c(nx, n_vars + ny), clip = "off") # Ajustar límites y espacio arriba

    # --- Combinar y Añadir Título ---
    final_plot <- p_main + p_table +
                  plot_layout(widths = c(2.5, 1.5)) + # Ajustar proporción anchos
                  plot_annotation(
                      title = plot_title,
                      subtitle = "Comparison of Unadjusted and Adjusted Beta Regression Coefficients (95% CI)",
                      theme = theme_pub_boxed + theme(plot.margin = margin(t=10,b=5)) # Reaplicar tema base sin margen inferior
                  )

    # --- Guardar Plot Combinado ---
    safe_filename <- sanitize_filename(filename_base);
    output_path_pdf <- file.path(output_analysis_dir, paste0(safe_filename, '.pdf'))
    output_path_png <- file.path(output_analysis_dir, paste0(safe_filename, '.png'))
    output_path_tiff <- file.path(output_analysis_dir, paste0(safe_filename, '.tif'))
    total_plot_width = 9.0 # Ancho total ajustado

    ggsave(output_path_pdf, plot=final_plot, width=total_plot_width, height=plot_height, device=cairo_pdf, limitsize=FALSE); cat(sprintf("  Saved Plot (PDF): %s\n", basename(output_path_pdf)))
    ggsave(output_path_png, plot=final_plot, width=total_plot_width, height=plot_height, dpi = HIGH_RES_DPI, device="png", limitsize=FALSE); cat(sprintf("  Saved Plot (PNG @ %d DPI): %s\n", HIGH_RES_DPI, basename(output_path_png)))
    ggsave(output_path_tiff, plot=final_plot, width=total_plot_width, height=plot_height, dpi = HIGH_RES_DPI, device="tiff", compression="lzw", limitsize=FALSE); cat(sprintf("  Saved Plot (TIFF @ %d DPI): %s\n", HIGH_RES_DPI, basename(output_path_tiff)))

    return(final_plot)
}


# --- Generar el plot combinado con resultados Beta Regression ---
if(nrow(univariate_beta_df_pooled) > 0 && nrow(multivariate_beta_df_pooled) > 0) {
    create_combined_beta_coefficient_plot(
        unadj_df = univariate_beta_df_pooled,
        adj_df = multivariate_beta_df_pooled,
        plot_title = paste("Factors Associated with PEG Recovery (%) - Beta Regression"),
        filename_base = "Fig_CoefficientPlot_Combined_BetaReg_Recovery_PubStyle",
        order_vars = PLOT_ORDER_VARIABLES # Usar orden predefinido
    )
} else {
    cat("Skipping combined beta coefficient plot (missing univariate or multivariate results).\n")
}


# =======================================================================
# 4. DIAGNÓSTICOS DEL MODELO BETA REGRESSION
# =======================================================================
cat("\n--- 4. Beta Regression Model Diagnostics ---\n")

# --- VIF Check (on first imputation) ---
if (!is.null(valid_multi_fits) && length(valid_multi_fits) > 0) {
    cat("Checking Variance Inflation Factors (VIF) using first valid model...\n")
    first_valid_fit <- valid_multi_fits[[1]]
    tryCatch({
        vif_values <- car::vif(first_valid_fit)
        print(vif_values)
        if(any(vif_values > 10)) {
             warning("High VIF values detected (> 10). Check for multicollinearity.", call. = FALSE)
        } else if (any(vif_values > 5)) {
             warning("Moderate VIF values detected (> 5). Consider checking multicollinearity.", call. = FALSE)
        } else {
             cat("  VIF values appear acceptable (all <= 5).\n")
        }
    }, error = function(e){ cat("Could not calculate VIF:", e$message, "\n") })
} else { cat("Skipping VIF check (no valid multivariate models found).\n") }

# --- Residual Plot (on first imputation) ---
if (!is.null(valid_multi_fits) && length(valid_multi_fits) > 0) {
    cat("Generating Residual vs Fitted plot (Quantile Residuals) using first valid model...\n")
    diag_plot_filename <- file.path(output_analysis_dir, "Fig_BetaReg_Residuals_Imp1.pdf")
    tryCatch({
      first_valid_fit <- valid_multi_fits[[1]]
      # Calcular cuantiles residuales (mejor para GLM/Beta) y fitted en la escala del link
      qres <- residuals(first_valid_fit, type = "quantile")
      link_fitted <- predict(first_valid_fit, type = "link")
      diag_df <- data.frame(Fitted = link_fitted, Residuals = qres)

    #   p_diag <- ggplot(diag_df, aes(x = Fitted, y = Residuals)) +
    #       geom_point(alpha=0.6, shape=1, color="orange") +
    #       geom_smooth(method="loess", color="blue", se=FALSE, linewidth=0.8) +
    #       geom_hline(yintercept = 0, linetype="dashed", color="black") +
    #       labs(title = "Quantile Residuals vs Linear Predictor",
    #            subtitle = "Based on first valid model fit",
    #            x = "Linear Predictor (Logit Scale)",
    #            y = "Quantile Residuals") +
    #       theme_bw(base_size = BASE_SIZE*1.3, base_family = FONT_FAMILY) +
    #       theme(plot.title=element_text(hjust=0.5), plot.subtitle=element_text(hjust=0.5))

    p_diag <- ggplot(diag_df, aes(x = Fitted, y = Residuals)) +
        geom_point(alpha = 0.6, shape = 1, color = "red2") +
        geom_smooth(method = "loess", color = "blue", se = FALSE, linewidth = 0.8) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        labs(title = "Quantile Residuals vs Linear Predictor",
             subtitle = "Based on first valid model fit",
             x = "Linear Predictor (Logit Scale)",
             y = "Quantile Residuals") +
        theme_bw(base_size = BASE_SIZE * 1.3, base_family = FONT_FAMILY) +
        theme(
    # Centrar títulos
              plot.title    = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5),

    # Eliminar rejilla mayor y menor
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),

    # Fondo transparente (opcional, theme_bw ya es blanco)
              panel.background = element_blank(),

    # Dibujar ejes en negro
              axis.line = element_line(color = "black"),

    # Asegurar que las líneas de los ejes tengan grosor suficiente
              axis.line.x = element_line(size = 0.5),
              axis.line.y = element_line(size = 0.5)
              )

      ggsave(diag_plot_filename, plot=p_diag, width=6, height=5, device=cairo_pdf)
      cat(sprintf("  Residual plot saved to: %s\n", normalizePath(diag_plot_filename)))
      cat("  Review residual plot for patterns or deviations from normality/homoscedasticity assumptions.\n")

    }, error = function(e){ cat("Could not generate residual plot:", e$message, "\n") })

} else { cat("Skipping residual plot (no valid multivariate models found).\n") }


# =======================================================================
# 5. GENERAR LEYENDA PARA ETIQUETAS COMPACTAS
# =======================================================================
cat("\n--- 5. Generating Legend for Compact Labels ---\n")
# Usar los resultados multivariados para obtener la lista de variables
if(nrow(multivariate_beta_df_pooled) > 0) {
    plot_labels_df <- multivariate_beta_df_pooled %>%
        select(Variable_Raw, Variable, Compact_Var_Label) %>%
        distinct() %>%
        filter(!is.na(Compact_Var_Label) & !is.na(Variable))

    # Filtrar para que coincida con el orden del plot si es posible
    ordered_labels <- data.frame(Variable_Raw = PLOT_ORDER_VARIABLES) %>%
                      left_join(plot_labels_df, by="Variable_Raw") %>%
                      filter(!is.na(Compact_Var_Label))

    if(nrow(ordered_labels) > 0) {
        legend_title_text <- "Variable Abbreviation Legend"; legend_separator <- paste(rep("-", nchar(legend_title_text)), collapse="");
        legend_items_text <- paste(sprintf("%-18s = %s", ordered_labels$Compact_Var_Label, ordered_labels$Variable), collapse = "\n");
        legend_text <- paste(legend_title_text, legend_separator, legend_items_text, sep="\n");
        legend_file_path <- file.path(output_analysis_dir, "Plot_Compact_Label_Legend_BetaReg.txt");
        tryCatch({ writeLines(legend_text, legend_file_path); cat(sprintf("Compact label legend saved to: %s\n", normalizePath(legend_file_path))) }, error = function(e){ cat(sprintf("ERROR saving label legend: %s\n", e$message)) });
        cat("\n"); cat(legend_text); cat("\n")
    } else { cat("Could not generate compact label legend (no matching variables found).\n") }
} else { cat("Could not generate compact label legend (no multivariate results).\n") }


# =======================================================================
# 6. MENSAJE FINAL
# =======================================================================
cat("\n--- Beta Regression Script Completed ---")
cat(sprintf("\nAnalysis performed using pooled results from MICE object (m=%d imputations).", mids_data_original$m))
cat(sprintf("\nAll outputs saved in '%s'.\n", normalizePath(output_analysis_dir)))
cat(sprintf("Script finished at: %s\n", Sys.time()))

