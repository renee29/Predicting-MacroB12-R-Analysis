# ===========================================================
#        Script R para Calcular Estadísticas Descriptivas
#        (Tabla 1 - Basado en Imputación #1 y Datos Pre-MICE)
# ===========================================================

# 0. CONFIGURACIÓN E IMPORTACIÓN
# ===========================================================

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

install_if_missing <- function(pkg) { if (!requireNamespace(pkg, quietly = TRUE)) { cat(sprintf("Installing package: %s\n", pkg)); tryCatch(install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org/"), error = function(e){ cat(sprintf("  Failed to install %s: %s\n", pkg, e$message))}) } }
sapply(required_packages, install_if_missing)

library(data.table); library(tidyverse); library(here); library(rstudioapi)
library(mice); library(knitr); library(kableExtra); library(dplyr); library(stats)

# --- !! ESTABLECER DIRECTORIO DE TRABAJO MANUALMENTE !! ---
# Asegúrate de que esta ruta sea EXACTAMENTE la carpeta donde está este script
#tryCatch({
#  setwd("F:/WorkInProgre/Collaborators/Chema/Data/R/DataCleanPG30_ToGitHub/ToSubmit")
#  cat("Manually set working directory to:", getwd(), "\n")
#}, error = function(e) {
#  warning("Could not manually set working directory. Paths might be incorrect.", call. = FALSE)
#})
# -----------------------------------------------------------
# --- Directorios ---
script_dir <- tryCatch({ dirname(rstudioapi::getSourceEditorContext()$path) }, error = function(e) { getwd() })
# Directorio donde se guardaron los resultados del script MICE
input_data_prep_dir <- file.path(script_dir, "R_Outputs_Data_Prep_New_30pct")
# Nuevo directorio para los resultados de este script de visualización
output_plot_dir <- file.path(script_dir, "R_Outputs_Table1_30pct")
if (!dir.exists(output_plot_dir)) { dir.create(output_plot_dir) } else { cat("Plot output directory exists:", normalizePath(output_plot_dir), "\n") }

# --- Constantes ---
RANDOM_STATE <- 42; set.seed(RANDOM_STATE)
IMPUTATION_TO_USE <- 1
INPUT_FILE_NAME_ORIGINAL <- "analysis_data_final_EngNames_1.csv" # Ajustar si es necesario
INPUT_FILE_PATH_ORIGINAL <- file.path(script_dir, INPUT_FILE_NAME_ORIGINAL)
MACROB12_THRESHOLD <- 30

# --- Cargar Mapeo de Etiquetas (Keys = Nombres REALES en analysis_df) ---
LABEL_MAP <- list(
  'EDAD' = 'Age (years)',
  'SEXO' = 'Sex',
  'Procedencia' = 'Source', # Sin punto
  'VB12.PRE.PEG' = 'VB12 Pre-PEG (pg/mL)',
  # 'VB12.POST.PEG' = 'VB12 Post-PEG (pg/mL)', # Nombre interno si existe
  'Recovery_pct' = 'PEG Recovery (%)',
  'MacroB12_Positive' = 'MacroB12 Status',
  'PCR' = 'CRP (mg/L)',
  'FOL' = 'Serum Folate (ng/mL)',
  'HB' = 'Hemoglobin (g/dL)',
  'VCM' = 'MCV (fL)',
  'ADE' = 'RDW (%)',
  # Nombres base de Flags (como existen en analysis_df)
  'Diag_Hematologic' = 'Diagnosis: Hematologic',
  'Diag_Hepatic' = 'Diagnosis: Hepatic',
  'Diag_Oncologic' = 'Diagnosis: Oncologic',
  'Diag_Rheum_Autoimmune' = 'Diagnosis: Rheum./Autoimmune',
  'Diag_Infectious' = 'Diagnosis: Infectious',
  'Diag_Renal' = 'Diagnosis: Renal',
  'Diag_GI_Suspicion' = 'Diagnosis: GI/Neo. Suspicion',
  'Diag_Hypertension' = 'Diagnosis: Hypertension',
  'Diag_Diabetes' = 'Diagnosis: Diabetes',
  'Diag_Healthy_Unspecified' = 'Diagnosis: Healthy/Unspecified',
  # Niveles (sin cambios)
  'Male' = 'Male', 'Female' = 'Female',
  'Ambulatory' = 'Ambulatory', 'Hospital' = 'Hospital', 'External Center' = 'External Center',
  'Negative'='Negative', 'Positive'='Positive'
  # No se necesitan términos de contraste como 'SEXOMale' aquí
)


# --- Helper Functions Definidas en Sección 0 ---
# get_label ahora usa LABEL_MAP con claves corregidas
get_label <- function(name, default_fmt=TRUE) {
    name_str <- as.character(name)
    label <- LABEL_MAP[[name_str]] # Busca directamente la clave (ej. 'EDAD', 'SEXO')
    if (!is.null(label)) return(label)
    # Fallback simple
    clean_name <- gsub("`", "", name_str)
    if (default_fmt) return(tools::toTitleCase(gsub("[._]", " ", clean_name))) else return(clean_name)
}
median_iqr <- function(x, digits = 1, na.rm = TRUE) { if (all(is.na(x))) return("-"); q <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = na.rm); sprintf(paste0("%.", digits, "f [%.", digits, "f – %.", digits, "f]"), q[2], q[1], q[3]) }
n_pct <- function(x, level=NULL, total_n, digits = 1, na.rm = TRUE) { if (is.factor(x) && is.null(level)) { level <- levels(x)[2] } else if (is.null(level)) { return("-") }; count <- sum(x == level, na.rm = na.rm); if(total_n == 0) return("0 (NA)"); pct <- (count / total_n) * 100; sprintf("%d (%.*f)", count, digits, pct) }
# get_n_valid se definirá después de cargar df_prelim
format_p_value <- function(p_vector) { sapply(p_vector, function(p) { if (is.na(p)) { return("NA") } else if (p < 0.0001) { return(sprintf("%.1e", p)) } else if (p < 0.001) { return("<0.001") } else if (p < 0.01) { return(sprintf("%.3f", p)) } else { return(sprintf("%.2f", p)) } }) }
add_stars <- function(p_vector) { sapply(p_vector, function(p) { if (is.na(p) || !is.numeric(p)) return('') ; if (p < 0.001) return('***') ; if (p < 0.01) return('**') ; if (p < 0.05) return('*') ; return('') }) }
sanitize_filename <- function(name) { name_sanitized <- gsub('[<>:"/\\|?*%[:cntrl:]]', '_', as.character(name)); name_sanitized <- gsub('[[:punct:]]', '_', name_sanitized); name_sanitized <- gsub('\\s+', '_', name_sanitized); name_sanitized <- gsub('_+', '_', name_sanitized); name_sanitized <- gsub('^_|_$', '', name_sanitized); name_sanitized <- substr(name_sanitized, 1, 100); return(name_sanitized) }

# =======================================================================
# 1. CARGAR DATOS (Imputados y Pre-Imputación)
# =======================================================================
cat("\n--- 1. Loading Processed and Pre-Imputation Data ---\n")
mice_rds_file <- file.path(input_data_prep_dir, "mice_imputation_object_internalNames.rds"); mids_data_original <- NULL; analysis_df <- NULL
if (file.exists(mice_rds_file)) { tryCatch({ mids_data_original <- readRDS(mice_rds_file); cat("Successfully loaded MICE object (.rds).\n"); if (IMPUTATION_TO_USE > 0 && IMPUTATION_TO_USE <= mids_data_original$m) { analysis_df <- mice::complete(mids_data_original, action = IMPUTATION_TO_USE); cat(sprintf("Extracted imputed dataset #%d for analysis.\n", IMPUTATION_TO_USE)) } else { stop(paste("Invalid imputation number:", IMPUTATION_TO_USE)) } }, error = function(e) { stop("Error loading MICE object:", e$message) }) } else { stop("MICE object not found:", mice_rds_file)}
df_prelim <- NULL
if(file.exists(INPUT_FILE_PATH_ORIGINAL)){ tryCatch({ dt_prelim <- data.table::fread(INPUT_FILE_PATH_ORIGINAL, na.strings = c("", "NA", "N/A", "#N/A", "#REF!", "#VALUE!", "Unknown", "no tiene", "..."), encoding = "UTF-8", check.names=FALSE, fill=TRUE, quote="\""); df_prelim <- as.data.frame(dt_prelim); original_names_prelim <- names(df_prelim); names(df_prelim) <- make.names(original_names_prelim, unique = TRUE); cat("Successfully loaded pre-imputation data from:", basename(INPUT_FILE_PATH_ORIGINAL), "\n") }, error = function(e){ warning("Could not load or process pre-imputation data:", e$message,"."); df_prelim <- NULL }) } else { warning("Pre-imputation file not found:", INPUT_FILE_PATH_ORIGINAL,".") }

# --- Definir get_n_valid AHORA ---
# Espera el nombre como existe en df_prelim (probablemente minúsculas/punto)
get_n_valid <- function(var_name_prelim) {
    if (!is.null(df_prelim) && var_name_prelim %in% names(df_prelim)) {
        return(sum(!is.na(df_prelim[[var_name_prelim]])))
    } else { return(NA) }
}

# --- Preparar analysis_df (tipos y Ns) ---
# !! Usar los nombres de columna REALES de analysis_df !!
target_col <- "MacroB12_Positive"
num_vars <- c('EDAD', 'VB12.PRE.PEG', 'PCR', 'FOL', 'HB', 'VCM', 'ADE', 'Recovery_pct')
cat_vars <- c('SEXO', 'Procedencia') # Sin punto!
diag_vars <- names(analysis_df)[startsWith(names(analysis_df), "Diag_")]

if (!target_col %in% names(analysis_df)) stop("Target column not found in analysis_df.")
analysis_df[[target_col]] <- factor(analysis_df[[target_col]], levels=c("Negative", "Positive"))
analysis_df <- analysis_df %>% filter(!is.na(!!sym(target_col)))

# Asegurar tipos correctos en analysis_df
cat("Ensuring correct types in analysis_df...\n")
for(col in intersect(num_vars, names(analysis_df))) { if(!is.numeric(analysis_df[[col]])) analysis_df[[col]] <- suppressWarnings(as.numeric(as.character(analysis_df[[col]]))) }
for(col in intersect(cat_vars, names(analysis_df))) { if(!is.factor(analysis_df[[col]])) analysis_df[[col]] <- factor(analysis_df[[col]]) }
for(col in intersect(diag_vars, names(analysis_df))) { if(!is.factor(analysis_df[[col]])) analysis_df[[col]] <- factor(analysis_df[[col]], levels=c(0,1)) }
# Relevel usando nombres correctos
if("SEXO" %in% names(analysis_df) && "Female" %in% levels(analysis_df[["SEXO"]])) analysis_df[["SEXO"]] <- relevel(analysis_df[["SEXO"]], ref="Female")
if("Procedencia" %in% names(analysis_df) && "Ambulatory" %in% levels(analysis_df[["Procedencia"]])) analysis_df[["Procedencia"]] <- relevel(analysis_df[["Procedencia"]], ref="Ambulatory")

n_total <- nrow(analysis_df); n_neg <- sum(analysis_df[[target_col]] == "Negative"); n_pos <- sum(analysis_df[[target_col]] == "Positive"); cat(sprintf("Final N for Table 1: %d (Negative: %d, Positive: %d)\n", n_total, n_neg, n_pos))

# =======================================================================
# 2. CALCULAR ESTADÍSTICAS DESCRIPTIVAS
# =======================================================================
cat("\n--- 2. Calculating Descriptive Statistics ---\n")
table1_results <- list()

# --- Calcular para cada variable (usando nombres REALES como claves y para acceso) ---
# Continuas
for(var in intersect(num_vars, names(analysis_df))) {
    digits <- ifelse(var %in% c('PCR','FOL','HB','VCM','ADE','Recovery_pct'), 1, 0)
    # --- Para N valid, mapear nombre REAL a nombre prelim (minúscula/punto) ---
    var_prelim_name <- case_when(
        var == "EDAD" ~ "edad",
        var == "VB12.PRE.PEG" ~ "vb12_pre_peg",
        var == "PCR" ~ "pcr",
        var == "FOL" ~ "fol",
        var == "HB" ~ "hb",
        var == "VCM" ~ "vcm",
        var == "ADE" ~ "ade",
        var == "Recovery_pct" ~ "Recovery_pct", # Asume que este no estaba en prelim
        TRUE ~ NA_character_ # No mapear otros
    )
    n_valid <- get_n_valid(var_prelim_name)
    # ---------------------------------------------------------------------
    table1_results[[var]] <- list( # CLAVE es el nombre REAL
         N_Valid = ifelse(is.na(n_valid), n_total, n_valid),
         Overall = median_iqr(analysis_df[[var]], digits), # Acceder con nombre REAL
         Negative = median_iqr(analysis_df[[var]][analysis_df[[target_col]]=="Negative"], digits),
         Positive = median_iqr(analysis_df[[var]][analysis_df[[target_col]]=="Positive"], digits)
    )
}

# Categóricas (Sexo, Procedencia)
for(var in intersect(cat_vars, names(analysis_df))) {
     if (is.factor(analysis_df[[var]])) {
         lvls <- levels(analysis_df[[var]])
         for(lvl in lvls) {
             var_level_name <- paste(var, lvl, sep="_") # Clave para resultados usa nombre REAL
             table1_results[[var_level_name]] <- list(
                 Overall = n_pct(analysis_df[[var]], lvl, n_total), # Acceder con nombre REAL
                 Negative = n_pct(analysis_df[[var]][analysis_df[[target_col]]=="Negative"], lvl, n_neg),
                 Positive = n_pct(analysis_df[[var]][analysis_df[[target_col]]=="Positive"], lvl, n_pos)
             )
         }
     } else {warning(paste("Variable", var, "is not a factor."))}
}

# Diagnósticos (Flags 0/1) - Nombres ya son correctos
for (dv in intersect(diag_vars, names(analysis_df))) {
     if (is.factor(analysis_df[[dv]])) {
        table1_results[[paste0(dv, "_Present")]] <- list(
            Overall = n_pct(analysis_df[[dv]], "1", n_total),
            Negative = n_pct(analysis_df[[dv]][analysis_df[[target_col]]=="Negative"], "1", n_neg),
            Positive = n_pct(analysis_df[[dv]][analysis_df[[target_col]]=="Positive"], "1", n_pos)
        )
     } else {warning(paste("Variable", dv, "is not a factor."))}
}

# Placeholder para límites B12
n_pre_gt <- "?"; n_post_lt <- "?"

# =======================================================================
# 3. MOSTRAR RESULTADOS EN FORMATO PARA LATEX
# =======================================================================
cat("\n--- 3. Results Formatted for LaTeX Table ---\n")

# --- Definir print_table_row aquí ---
# Ahora espera el nombre REAL como var_key
print_table_row <- function(var_key, label_override=NULL, type="cont", level=NULL, n_valid_col=NULL) {
    # Obtener etiqueta completa usando el mapeo LABEL_MAP (cuyas claves ahora son nombres reales)
    label <- ifelse(!is.null(label_override), label_override, get_label(var_key))
    result_key <- var_key
    prefix = "\\quad "
    if(type=="cat_header") {result_key <- NULL; prefix=""}
    if(type=="cat_level") {result_key <- paste(var_key, level, sep="_"); label=level} # Clave de resultado para nivel
    if(type=="diag") {result_key <- paste(var_key, "Present", sep="_"); prefix="\\quad "; label=get_label(var_key)} # Usar get_label para el nombre base diag

    if(!is.null(result_key) && result_key %in% names(table1_results)) {
        res <- table1_results[[result_key]]
        if(!is.null(res)) {
            n_valid_text <- ""
            # Usar N_Valid ya calculado y guardado en la lista
            if(!is.null(n_valid_col) && !is.null(res$N_Valid) && !is.na(res$N_Valid) && res$N_Valid < n_total) {
                 n_valid_text <- sprintf(" (N valid=%d)", res$N_Valid)
            }
            overall_val <- res$Overall; neg_val <- res$Negative; pos_val <- res$Positive
            # Añadir \texttt{} para nombres originales si se desea
            # label_latex = gsub("_", "\\\\_", label, fixed=TRUE) # Escapar guiones bajos para LaTeX
            label_latex = label # Usar etiqueta directamente
            cat(sprintf("%s%s%s & %s & %s & %s \\\\\n", prefix, label_latex, n_valid_text, overall_val, neg_val, pos_val))
        } else { cat(sprintf("%s%s & - & - & - \\\\ %% Result key '%s' not found\n", prefix, label, result_key)) }
    } else if (type=="cat_header"){ cat(sprintf("%s%s, N (\\%%) & & & \\\\\n", prefix, label))
    } else { cat(sprintf("%s%s & - & - & - \\\\\n", prefix, label)) }
}

# --- Imprimir Tabla ---
cat(sprintf("\\caption{Baseline Characteristics of the Study Cohort (N=%d)} %% Ajusta caption \n", n_total))
cat(sprintf("\\label{tab:baseline_char} \n"))
cat("\\resizebox{\\textwidth}{!}{% \n")
cat("\\begin{tabular}{l c c c} \n")
cat("\\toprule \n")
cat(sprintf("Characteristic & Overall Cohort & Macro B12 Negative & Macro B12 Positive \\\\ \n"))
cat(sprintf(" & (N = %d) & (N = %d) & (N = %d) \\\\ \n", n_total, n_neg, n_pos))
cat("\\midrule \n")

# Demographics
cat("\\textbf{Demographics} & & & \\\\\n")
print_table_row("EDAD") # Usar clave REAL
print_table_row("SEXO", type="cat_header") # Usar clave REAL
print_table_row("SEXO", type="cat_level", level="Male")
print_table_row("SEXO", type="cat_level", level="Female")
cat("\\addlinespace\n")

# Clinical Context
cat("\\textbf{Clinical Context} & & & \\\\\n")
print_table_row("Procedencia", type="cat_header") # Usar clave REAL
print_table_row("Procedencia", type="cat_level", level="Ambulatory")
print_table_row("Procedencia", type="cat_level", level="Hospital")
print_table_row("Procedencia", type="cat_level", level="External Center")

# Diagnósticos
cat("\\textbf{Diagnosis Present (N (\\%))} & & & \\\\\n")
diag_vars_ordered <- c("Diag_Healthy_Unspecified", "Diag_Hematologic", "Diag_Hepatic", "Diag_Oncologic", "Diag_Rheum_Autoimmune", "Diag_Infectious", "Diag_Renal", "Diag_GI_Suspicion", "Diag_Hypertension", "Diag_Diabetes") # Nombres base internos
for(dv in intersect(diag_vars_ordered, names(analysis_df))){
    print_table_row(dv, type="diag") # Pasa el nombre interno/real de la flag
}
cat("\\addlinespace\n")

# Vitamin B12 Parameters
cat("\\textbf{Vitamin B12 Parameters} & & & \\\\\n")
# Usar clave REAL, y la etiqueta usa get_label con esa clave
print_table_row("VB12.PRE.PEG", label_override = paste0(get_label("VB12.PRE.PEG"), "**"))
print_table_row("Recovery_pct")
cat("\\addlinespace\n")

# Concomitant Laboratory Markers
cat("\\textbf{Concomitant Laboratory Markers} & & & \\\\\n")
lab_vars_ordered <- c("PCR", "FOL", "HB", "VCM", "ADE") # Nombres REALES
for(lv in intersect(lab_vars_ordered, names(analysis_df))){
    print_table_row(lv, n_valid_col=TRUE) # Pasar nombre REAL
}

# Final de tabla
cat("\\bottomrule \n")
cat("\\end{tabular} \n")
cat("} % End resizebox \n")
cat(sprintf("\\caption*{\\small Data presented as Median [Interquartile Range] for continuous variables or N (\\%%) for categorical variables, based on imputed dataset \\#%d. MacroB12 Positive defined as PEG recovery < %.0f\\%%. P-values omitted. IQR: Interquartile Range. CRP: C-reactive protein. MCV: Mean Corpuscular Volume. RDW: Red Cell Distribution Width. \\\\\n", IMPUTATION_TO_USE, MACROB12_THRESHOLD))
cat(sprintf("    ** Includes %s values reported as > limit in original data. Values shown are from imputed data.} \n", n_pre_gt)) # Simplificado nota al pie
cat("\\end{table} \n")

# =======================================================================
# 4. GUARDAR RESULTADOS NUMÉRICOS DE TABLA 1 EN CSV  <--- NUEVA SECCIÓN
# =======================================================================
cat("\n--- 4. Saving Table 1 Numerical Results to CSV ---\n")

# Convertir la lista table1_results a un data frame
# Cada elemento de table1_results es una lista con Overall, Negative, Positive (y N_Valid para continuas)
# Necesitamos reestructurar esto.

table1_df_rows <- list()
row_counter <- 1

# Helper para añadir filas al data frame
add_table_row_to_df <- function(variable_label, characteristic_level, overall, negative, positive, n_valid = NA_integer_) {
  table1_df_rows[[row_counter]] <<- data.frame(
    Variable_Label = variable_label,
    Characteristic_Level = characteristic_level,
    N_Valid_Original = n_valid, # N válidos del archivo pre-imputación
    Overall_Stat = overall,
    Negative_Group_Stat = negative,
    Positive_Group_Stat = positive,
    stringsAsFactors = FALSE
  )
  row_counter <<- row_counter + 1
}

# Procesar Variables Continuas
for(var_key in intersect(num_vars, names(analysis_df))) { # num_vars son los NOMBRES FINALES
    res <- table1_results[[var_key]]
    if (!is.null(res)) {
        add_table_row_to_df(
            variable_label = get_label(var_key), # Usa el nombre final para la etiqueta
            characteristic_level = "Median [IQR]",
            overall = res$Overall,
            negative = res$Negative,
            positive = res$Positive,
            n_valid = res$N_Valid # N_Valid ya está en res
        )
    }
}

# Procesar Variables Categóricas (Sexo, Procedencia)
for(var_key in intersect(cat_vars, names(analysis_df))) { # cat_vars son los NOMBRES FINALES
    if (is.factor(analysis_df[[var_key]])) {
        lvls <- levels(analysis_df[[var_key]])
        # Añadir una fila para el nombre de la variable categórica como cabecera de sección
        add_table_row_to_df(
            variable_label = get_label(var_key),
            characteristic_level = "N (%)", # Tipo de estadística para toda la categoría
            overall = "", negative = "", positive = "", n_valid = n_total # N para toda la categoría
        )
        for(lvl in lvls) {
            result_key_for_level <- paste(var_key, lvl, sep="_")
            res <- table1_results[[result_key_for_level]]
            if (!is.null(res)) {
                add_table_row_to_df(
                    variable_label = "", # Dejar en blanco o indentar para niveles
                    characteristic_level = get_label(lvl), # Etiqueta del nivel
                    overall = res$Overall,
                    negative = res$Negative,
                    positive = res$Positive,
                    n_valid = NA # N_valid no aplica directamente a cada nivel así
                )
            }
        }
    }
}

# Procesar Variables de Diagnóstico (Flags)
# Añadir una fila para la cabecera de la sección de Diagnósticos
add_table_row_to_df(
    variable_label = "Diagnosis Present",
    characteristic_level = "N (%)",
    overall = "", negative = "", positive = "", n_valid = n_total
)
diag_vars_ordered <- c("Diag_Healthy_Unspecified", "Diag_Hematologic", "Diag_Hepatic", "Diag_Oncologic", "Diag_Rheum_Autoimmune", "Diag_Infectious", "Diag_Renal", "Diag_GI_Suspicion", "Diag_Hypertension", "Diag_Diabetes")
for (dv_key in intersect(diag_vars_ordered, names(analysis_df))) { # dv_key es el NOMBRE FINAL
    result_key_for_diag <- paste0(dv_key, "_Present")
    res <- table1_results[[result_key_for_diag]]
    if (!is.null(res)) {
        add_table_row_to_df(
            variable_label = "", # Dejar en blanco o indentar
            characteristic_level = get_label(dv_key), # Etiqueta del diagnóstico
            overall = res$Overall,
            negative = res$Negative,
            positive = res$Positive,
            n_valid = NA # N_valid no aplica directamente aquí
        )
    }
}

# Combinar todas las filas en un solo data frame
table1_summary_df <- do.call(rbind, table1_df_rows)

# Guardar el data frame como CSV
table1_csv_path <- file.path(output_plot_dir, "Table1_Descriptive_Statistics_Numerical.csv")
tryCatch({
  write_csv(table1_summary_df, table1_csv_path, na = "-") # Escribir NAs como "-"
  cat(sprintf("Table 1 numerical results successfully saved to: %s\n", normalizePath(table1_csv_path)))
}, error = function(e) {
  cat(sprintf("Error saving Table 1 numerical results to CSV: %s\n", e$message))
})

# Aquí iría tu sección 3. MOSTRAR RESULTADOS EN FORMATO PARA LATEX
# O puedes poner esta sección de guardado de CSV al final de todo el script.

# =======================================================================
# 5. MENSAJE FINAL
# =======================================================================
cat("\n--- Table 1 Calculation Script Completed ---")
cat(sprintf("\nCalculations based on Imputation #%d.", IMPUTATION_TO_USE))
cat(sprintf("\nCopy the LaTeX output above and fill in any remaining placeholders ('?').\n"))
cat(sprintf("Output saved in '%s'.\n", normalizePath(output_plot_dir)))
cat(sprintf("Script finished at: %s\n", Sys.time()))