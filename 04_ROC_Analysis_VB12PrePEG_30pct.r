# ===========================================================
#        Script R para Análisis ROC y Umbral Óptimo
#        (Estilo Publicación - Caja Cerrada Manual, Fuente Grande, Umbral Natural)
# ===========================================================

# 0. CONFIGURACIÓN E IMPORTACIÓN
# ===========================================================
# --- Cargar Paquetes ---
required_packages <- c("data.table", "tidyverse", "ggplot2", "here", "rstudioapi",
                       "mice", "pROC", "plotROC", "ggthemes", "RColorBrewer",
                       "ggrepel", "ggdist", "knitr", "kableExtra", "skimr",
                       "ggpubr", "paletteer", "dplyr", "stats", "forcats",
                       "svglite") # Para guardar SVG

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing package: %s\n", pkg))
    tryCatch(install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org/"),
             error = function(e){ cat(sprintf("  Failed to install %s: %s\n", pkg, e$message))})
  }
}
sapply(required_packages, install_if_missing)

library(data.table); library(tidyverse); library(ggplot2); library(here); library(rstudioapi)
library(mice); library(pROC); library(plotROC)
library(ggthemes); library(RColorBrewer); library(ggrepel); library(ggdist)
library(knitr); library(kableExtra); library(skimr); library(ggpubr); library(paletteer)
library(dplyr); library(stats); library(forcats)
library(svglite)

# --- Directorios ---
script_dir <- tryCatch({ dirname(rstudioapi::getSourceEditorContext()$path) }, error = function(e) { getwd() })
input_data_prep_dir <- file.path(script_dir, "R_Outputs_Data_Prep_New_30pct")
# Directorio de salida con estilo actualizado
output_roc_dir <- file.path(script_dir, "ROC_Analysis_VB12PrePEG_30pct")
if (!dir.exists(output_roc_dir)) { dir.create(output_roc_dir) } else { cat("ROC Output directory exists:", normalizePath(output_roc_dir), "\n") }

# --- Constantes y Estilos Globales ---
RANDOM_STATE <- 42; set.seed(RANDOM_STATE)
IMPUTATION_NUMBER_TO_ANALYZE <- 1
FONT_FAMILY <- "Arial"
# T3: Aumentar tamaño base de fuente
BASE_SIZE <- 11 # Incrementado de 9 a 11 (ajustar según necesidad)
HIGH_RES_DPI <- 600

# Colores (sin cambios)
color_roc_curve <- "dodgerblue3" # Slightly brighter blue maybe
color_optimal_point <- "red3"
color_sensitivity <- color_roc_curve # Match ROC curve
color_specificity <- "coral3"

# ===========================================================
#        Mapeo de Etiquetas y Funciones Helper
# ===========================================================
# --- MAPEO DE ETIQUETAS COMPLETAS ---
LABEL_MAP <- list(
  'EDAD' = 'Age (years)', 'VB12.PRE.PEG' = 'VB12 Pre-PEG (pg/mL)', 'PCR' = 'CRP (mg/L)',
  'FOL' = 'Serum Folate (ng/mL)', 'HB' = 'Hemoglobin (g/dL)', 'VCM' = 'MCV (fL)', 'ADE' = 'RDW (%)',
  'SEXOMale' = 'Sex: Male vs Female', 'ProcedenciaExternal Center' = 'Source: External Center vs Ambulatory',
  'ProcedenciaHospital' = 'Source: Hospital vs Ambulatory',
  'Diag_Hematologic1' = 'Diagnosis: Hematologic: Present vs Absent', 'Diag_Hepatic1' = 'Diagnosis: Hepatic: Present vs Absent',
  'Diag_Oncologic1' = 'Diagnosis: Oncologic: Present vs Absent', 'Diag_Rheum_Autoimmune1' = 'Diagnosis: Rheum./Autoimmune: Present vs Absent',
  'Diag_Infectious1' = 'Diagnosis: Infectious: Present vs Absent', 'Diag_Renal1' = 'Diagnosis: Renal: Present vs Absent',
  'Diag_GI_Suspicion1' = 'Diagnosis: GI/Neo. Suspicion: Present vs Absent', 'Diag_Hypertension1' = 'Diagnosis: Hypertension: Present vs Absent',
  'Diag_Diabetes1' = 'Diagnosis: Diabetes: Present vs Absent', 'Diag_Healthy_Unspecified1' = 'Diagnosis: Healthy/Unspecified: Present vs Absent',
  # Nombres Originales/Base
  'edad' = 'Age (years)', 'sexo' = 'Sex', 'procedencia.' = 'Source', 'vb12_pre_peg' = 'VB12 Pre-PEG (pg/mL)',
  'pcr' = 'CRP (mg/L)', 'fol' = 'Serum Folate (ng/mL)', 'hb' = 'Hemoglobin (g/dL)', 'vcm' = 'MCV (fL)',
  'ade' = 'RDW (%)', 'Diag_Hematologic' = 'Diagnosis: Hematologic', 'Diag_Hepatic' = 'Diagnosis: Hepatic',
  'Diag_Oncologic' = 'Diagnosis: Oncologic', 'Diag_Rheum_Autoimmune' = 'Diagnosis: Rheum./Autoimmune',
  'Diag_Infectious' = 'Diagnosis: Infectious', 'Diag_Renal' = 'Diagnosis: Renal',
  'Diag_GI_Suspicion' = 'Diagnosis: GI/Neo. Suspicion', 'Diag_Hypertension' = 'Diagnosis: Hypertension',
  'Diag_Diabetes' = 'Diagnosis: Diabetes', 'Diag_Healthy_Unspecified' = 'Diagnosis: Healthy/Unspecified',
  # Otros
  'MacroB12_Positive' = 'MacroB12 Status', 'Male' = 'Male', 'Female' = 'Female',
  'Ambulatory' = 'Ambulatory', 'Hospital' = 'Hospital', 'External Center' = 'External Center'
)

# --- Cargar mapeo de nombres finales (si existe) ---
final_name_map_file <- file.path(input_data_prep_dir, "column_rename_map_mk_to_eng.rds")
internal_to_final_map <- NULL; final_to_internal_map <- list()
if(file.exists(final_name_map_file)){ tryCatch({ internal_to_final_map <- readRDS(final_name_map_file); final_to_internal_map <- setNames(names(internal_to_final_map), unlist(internal_to_final_map)) }, error=function(e){}) }

# --- Helper Functions ---
get_label <- function(name, default_fmt=TRUE) {
     name_str <- as.character(name); internal_name <- name_str ; if(exists("col_names_type") && col_names_type == "final"){ internal_lookup <- final_to_internal_map[[name_str]]; if(!is.null(internal_lookup)) internal_name <- internal_lookup }; label <- LABEL_MAP[[internal_name]]; if (!is.null(label)) return(label); clean_name <- gsub("`", "", name_str); if (default_fmt) return(tools::toTitleCase(gsub("[._]", " ", clean_name))) else return(clean_name)
}
sanitize_filename <- function(name) { name_sanitized <- gsub('[<>:"/\\|?*%[:cntrl:]]', '_', as.character(name)); name_sanitized <- gsub('[[:punct:]]', '_', name_sanitized); name_sanitized <- gsub('\\s+', '_', name_sanitized); name_sanitized <- gsub('_+', '_', name_sanitized); name_sanitized <- gsub('^_|_$', '', name_sanitized); name_sanitized <- substr(name_sanitized, 1, 100); return(name_sanitized) }

# =======================================================================
# 1. CARGAR DATASET SELECCIONADO
# =======================================================================
cat("\n--- 1. Loading Processed Data for ROC Analysis ---\n")
mice_rds_file <- file.path(input_data_prep_dir, "mice_imputation_object_internalNames.rds")
csv_file_path <- file.path(input_data_prep_dir, paste0("analysis_data_final_EngNames_", IMPUTATION_NUMBER_TO_ANALYZE, ".csv"))
roc_df <- NULL; data_source_msg <- NULL; col_names_type <- NULL
if (file.exists(mice_rds_file)) {
    tryCatch({
        mids_data <- readRDS(mice_rds_file); cat("Successfully loaded MICE object (.rds).\n")
        imp_num <- if (IMPUTATION_NUMBER_TO_ANALYZE > 0 && IMPUTATION_NUMBER_TO_ANALYZE <= mids_data$m) IMPUTATION_NUMBER_TO_ANALYZE else 1
        if(imp_num != IMPUTATION_NUMBER_TO_ANALYZE) cat(sprintf("Warning: Imputation number %d not valid (m=%d). Using imputation #%d.\n", IMPUTATION_NUMBER_TO_ANALYZE, mids_data$m, imp_num))
        roc_df <- mice::complete(mids_data, action = imp_num)
        cat(sprintf("Extracted imputed dataset #%d for ROC analysis.\n", imp_num)); data_source_msg <- sprintf("Based on Imputed Dataset #%d", imp_num); col_names_type <- "internal"
    }, error = function(e) { cat("Error loading MICE object or extracting dataset:", e$message, "\n"); roc_df <<- NULL })
}
if (is.null(roc_df)) { # Fallback
    cat("MICE object not used. Attempting to load final name CSV.\n")
    if (file.exists(csv_file_path)) {
        tryCatch({ roc_df <- as.data.frame(fread(csv_file_path, na.strings = c("", "NA"))); cat("Successfully loaded data with final names from CSV:", basename(csv_file_path), "\n"); data_source_msg <- "Based on Single Dataset"; col_names_type <- "final"
        }, error = function(e) { cat("Error loading final name CSV:", e$message, "\n"); roc_df <<- NULL })
    } else { cat("Final name CSV file not found:", csv_file_path, "\n") }
}
if (is.null(roc_df)) stop("Could not load data for ROC analysis.")

# --- Asegurar tipos de columna correctos ---
cat("Ensuring correct column types for ROC analysis...\n")
if(col_names_type == "internal"){
    target_col_roc <- make.names("MacroB12_Positive")
    num_features_expected <- make.names(c('EDAD', 'VB12.PRE.PEG', 'PCR', 'FOL', 'HB', 'VCM', 'ADE'))
} else { # final names
    target_col_roc <- "MacroB12_Status"
    num_features_expected <- c('Age', 'VB12_PrePEG', 'CRP', 'Folate', 'Hemoglobin', 'MCV', 'RDW')
}
if (!target_col_roc %in% names(roc_df)) stop("Target column not found in loaded data.")
if (!is.factor(roc_df[[target_col_roc]])) { roc_df[[target_col_roc]] <- factor(roc_df[[target_col_roc]]) }
roc_df[[target_col_roc]] <- factor(roc_df[[target_col_roc]], levels = c("Negative", "Positive"))
roc_df <- roc_df %>% filter(!is.na(!!sym(target_col_roc)))
cat("Target variable levels:", paste(levels(roc_df[[target_col_roc]]), collapse=", "), "\n")
for(col in intersect(num_features_expected, names(roc_df))){ if(!is.numeric(roc_df[[col]])) roc_df[[col]] <- suppressWarnings(as.numeric(as.character(roc_df[[col]]))) }
cat("Data ready for ROC analysis. Dimensions:", dim(roc_df), "\n")

# =======================================================================
# 2. GENERAR PLOTS ROC Y UMBRAL ÓPTIMO (CAJA MANUAL, UMBRAL NATURAL, FUENTE GRANDE)
# =======================================================================
cat("\n--- 2. Generating Boxed ROC and Threshold Plots (Natural Threshold) ---\n")

# --- Variables numéricas para las que generar plots ROC ---
if (col_names_type == "internal") {
    roc_vars_to_plot <- intersect(make.names(c('VB12.PRE.PEG', 'EDAD', 'PCR', 'FOL', 'HB', 'VCM', 'ADE')), names(roc_df))
} else {
    roc_vars_to_plot <- intersect(c('VB12_PrePEG', 'Age', 'CRP', 'Folate', 'Hemoglobin', 'MCV', 'RDW'), names(roc_df))
}
roc_vars_to_plot <- roc_vars_to_plot[sapply(roc_df[, roc_vars_to_plot, drop=FALSE], is.numeric)]
cat("Generating ROC plots for:", paste(roc_vars_to_plot, collapse=", "), "\n")

# --- TEMA ggplot base (Sin líneas de eje, las añadiremos manualmente) ---
theme_pub_roc_manual_box <- theme_bw(base_size = BASE_SIZE, base_family = FONT_FAMILY) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1.2), margin = margin(b=4)),
    plot.subtitle = element_text(hjust = 0.5, size = rel(1.1), colour = "grey30", margin = margin(b=12)),
    axis.title = element_text(face = "plain", size = rel(1.05)),
    axis.text = element_text(size = rel(1.0), colour = "black"),
    # Eliminar rejilla, borde de panel, Y TODAS las líneas de eje del tema
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    # Ticks solo en bottom y left (se mantienen)
    axis.ticks.x.top = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.ticks.x.bottom = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.y.left = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    # Leyenda (estilos generales)
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_text(size=rel(1.0), face="bold"),
    legend.text = element_text(size=rel(1.0)),
    plot.margin = margin(10, 15, 10, 10)
  )

# --- Bucle para generar plots ---
for (predictor_var in roc_vars_to_plot) {

    predictor_label <- get_label(predictor_var)
    cat(sprintf("  Processing ROC for: %s (%s)...\n", predictor_var, predictor_label))

    roc_data_current <- roc_df %>% filter(!is.na(!!sym(predictor_var)))
    if(nrow(roc_data_current) < 10 || n_distinct(roc_data_current[[target_col_roc]]) < 2) {
        cat("    Skipping ROC: Insufficient data or groups.\n"); next
    }

    # --- T2: Calcular Dirección Correcta (Revisada según AUCs observados) ---
    # Identificar nombre interno para comparación consistente
    predictor_internal_name <- predictor_var
    if(col_names_type == "final"){
        internal_lookup <- final_to_internal_map[[predictor_var]]
        if(!is.null(internal_lookup)) predictor_internal_name <- internal_lookup
    }

    # Definir dirección basada en la relación observada o esperada con el estado "Positive"
    if (predictor_internal_name %in% make.names(c("VB12.PRE.PEG", "HB", "MCV", "EDAD"))) {
        # Para estas, valores BAJOS se asocian con "Positive" (basado en AUCs < 0.5 para Hb, MCV, Age y >0.5 para VB12)
        roc_direction <- "<"
    } else if (predictor_internal_name %in% make.names(c("RDW", "PCR"))) {
        # Para estas, valores ALTOS se asocian con "Positive" (basado en AUCs < 0.5 para ellas)
        roc_direction <- ">"
    } else {
        # Para Folate (AUC~0.5) u otras no listadas, dejar que pROC decida
        roc_direction <- "auto"
        cat(sprintf("    ROC direction for %s set to 'auto'.\n", predictor_internal_name))
    }
    # Imprimir la dirección final usada
    cat(sprintf("    Using FINAL ROC direction: '%s'\n", roc_direction))

    # --- Calcular ROC y Umbral Óptimo (Youden) ---
    roc_obj <- tryCatch({
        pROC::roc(response = roc_data_current[[target_col_roc]],
                  predictor = roc_data_current[[predictor_var]],
                  levels = c("Negative", "Positive"), direction = roc_direction, # Usar dirección corregida
                  ci = TRUE, quiet = TRUE)
    }, error = function(e) { cat("    ERROR calculating ROC:", e$message, "\n"); NULL })

    if(is.null(roc_obj)) next

    # Recalcular AUC con la dirección correcta
    auc_val <- pROC::auc(roc_obj); auc_ci <- pROC::ci.auc(roc_obj);
    # Ahora todos los AUC deberían ser >= 0.5 si hay alguna discriminación
    if (auc_val < 0.5 && roc_direction != "auto") {
      # Esta advertencia ahora indicaría un problema genuino si aparece
      cat(sprintf("    PERSISTENT WARNING: AUC is %.3f (< 0.5) even with direction '%s' for %s.\n",
                  auc_val, roc_direction, predictor_var))
    }
    auc_text <- sprintf("AUC = %.3f (95%% CI %.3f–%.3f)", auc_val, auc_ci[1], auc_ci[3])

# Calcular coords óptimas (Youden) - esto no cambia
    opt_coords <- tryCatch({
      pROC::coords(roc_obj, "best", best.method="youden", ret=c("threshold", "sensitivity", "specificity"), transpose = FALSE)
    }, error = function(e) {cat("    ERROR getting coords for 'best' threshold:", e$message, "\n"); NULL})

    point_to_plot_df <- data.frame(specificity = numeric(0), sensitivity = numeric(0))
    optimal_vals_text <- "Optimal point could not be determined."
    opt_threshold <- NA; opt_sens <- NA; opt_spec <- NA

    if (!is.null(opt_coords) && nrow(opt_coords) > 0) {
        opt_threshold <- opt_coords$threshold[1]
        opt_sens <- opt_coords$sensitivity[1]
        opt_spec <- opt_coords$specificity[1]
        point_to_plot_df <- data.frame(specificity = opt_spec, sensitivity = opt_sens)
        optimal_vals_text <- sprintf("Optimal threshold (Youden index) = %.1f (Sensitivity = %.1f%%, Specificity = %.1f%%).", opt_threshold, opt_sens*100, opt_spec*100)
        if (predictor_internal_name == make.names("VB12.PRE.PEG")) {
             cat(sprintf("    >>> Calculated Youden Threshold for %s: %.1f\n", predictor_var, opt_threshold))
        }
    } else {
        cat("    Could not calculate optimal coordinates. Skipping point display.\n")
    }
    cat("   ", optimal_vals_text, "\n")

    # --- Plot 1: Curva ROC (Estilo Publicación - Caja Manual) ---
    roc_plot_title <- paste("ROC Curve:", predictor_label)
    roc_filename_base <- paste0("Fig_ROC_", sanitize_filename(predictor_var), "_BoxedNatThresh")

    # Definir límites para annotate (basados en escala 0-1)
    x_min_ann <- 0; x_max_ann <- 1; y_min_ann <- 0; y_max_ann <- 1

    p_roc <- ggplot(roc_data_current, aes(d = !!sym(target_col_roc), m = !!sym(predictor_var))) +
        # Línea diagonal
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
        # Curva ROC
        geom_roc(n.cuts = 0, colour = color_roc_curve, size = 1.0) +  # n.cuts=NULL para escalones
        # Punto óptimo
        geom_point(data = point_to_plot_df, aes(x = 1 - specificity, y = sensitivity),
                   color = color_optimal_point, size = 3.5, shape = 16, inherit.aes = FALSE) +

        # Anotación de Texto Umbral (Reactivada)
        { if (!is.na(opt_sens) && !is.na(opt_spec) && !is.na(opt_threshold)) {
             annotate("label", x = 1 - opt_spec, y = opt_sens, parse = FALSE,
                      label = sprintf("Threshold: %.1f\nSens: %.1f%%\nSpec: %.1f%%", opt_threshold, opt_sens*100, opt_spec*100),
                      hjust = ifelse( (1-opt_spec) > 0.6 | opt_sens < 0.4 , 1.1, -0.1),
                      vjust = ifelse(opt_sens < 0.2, -0.2, 1.2),
                      size = BASE_SIZE * 0.30, # Tamaño ajustado
                      family = FONT_FAMILY,
                      lineheight = 0.9, fill = "white", label.size = 0.15,
                      label.padding = unit(0.15, "lines"), colour = "black")
          }
        } +

        # Escalas y Límites
        scale_x_continuous(name="1 – Specificity (FPR)", limits=c(-0.01,1.01), expand=c(0, 0), breaks = seq(0, 1, 0.25)) +
        scale_y_continuous(name="Sensitivity (TPR)", limits=c(-0.01,1.01), expand=c(0, 0), breaks = seq(0, 1, 0.25)) +
        # Título y Subtítulo
        labs(title = roc_plot_title, subtitle = auc_text) +
        # Añadir líneas de eje manualmente usando annotate y límites 0-1
        annotate("segment", x=x_min_ann, xend=x_max_ann, y=y_min_ann, yend=y_min_ann, colour = "black", linewidth=0.5) + # Bottom
        annotate("segment", x=x_min_ann, xend=x_max_ann, y=y_max_ann, yend=y_max_ann, colour = "black", linewidth=0.5) + # Top
        annotate("segment", x=x_min_ann, xend=x_min_ann, y=y_min_ann, yend=y_max_ann, colour = "black", linewidth=0.5) + # Left
        annotate("segment", x=x_max_ann, xend=x_max_ann, y=y_min_ann, yend=y_max_ann, colour = "black", linewidth=0.5) + # Right
        # Forzar límites y permitir anotaciones fuera
        coord_cartesian(xlim = c(0,1), ylim = c(0,1), clip = "off", expand = FALSE) + # expand=FALSE para pegar a los bordes
        # Aplicar tema base (sin líneas de eje propias)
        theme_pub_roc_manual_box +
        theme(legend.position = "none")

    # Guardar en múltiples formatos
    ggsave(file.path(output_roc_dir, paste0(roc_filename_base, ".pdf")), plot = p_roc, width = 5.5, height = 5, device=cairo_pdf)
    ggsave(file.path(output_roc_dir, paste0(roc_filename_base, ".png")), plot = p_roc, width = 5.5, height = 5, dpi = HIGH_RES_DPI, device="png")
    ggsave(file.path(output_roc_dir, paste0(roc_filename_base, ".tif")), plot = p_roc, width = 5.5, height = 5, dpi = HIGH_RES_DPI, device="tiff", compression="lzw")
    cat(sprintf("    Saved ROC plot (PDF, PNG, TIFF): %s\n", roc_filename_base))


    # --- Plot 2: Sens/Spec vs. Umbral (Estilo Publicación - Caja Manual) ---
    roc_data_plot <- data.frame(thresholds = roc_obj$thresholds, sensitivity = roc_obj$sensitivities, specificity = roc_obj$specificities) %>%
      filter(is.finite(thresholds)) %>% arrange(thresholds)
    sens_spec_plot_title <- paste("Performance Metrics vs. Threshold:", predictor_label)
    sens_spec_filename_base <- paste0("Fig_SensSpec_", sanitize_filename(predictor_var), "_BoxedNatThresh")
    plot_data_long <- roc_data_plot %>%
      select(thresholds, sensitivity, specificity) %>%
      pivot_longer(cols = c("sensitivity", "specificity"), names_to = "Metric", values_to = "Value") %>%
      mutate(Metric = factor(Metric, levels = c("sensitivity", "specificity"), labels = c("Sensitivity", "Specificity")))
    finite_threshold_range <- range(roc_data_plot$thresholds, na.rm = TRUE)
    x_min_sens <- finite_threshold_range[1]; x_max_sens <- finite_threshold_range[2];
    y_min_sens <- 0; y_max_sens <- 1;

    p_sens_spec <- ggplot(plot_data_long, aes(x = thresholds, y = Value, color = Metric, group = Metric)) +
        geom_line(linewidth = 1.0) +
        geom_vline(xintercept = opt_threshold, linetype = "dotted", color = "grey50", linewidth=0.6) +
        scale_color_manual(values = c("Sensitivity" = color_sensitivity, "Specificity" = color_specificity), name=NULL) +
        scale_y_continuous(name = "Metric Value", limits = c(-0.01, 1.01), expand = c(0,0), breaks=seq(0,1,0.25)) +
        scale_x_continuous(limits = finite_threshold_range, expand=c(0.01, 0.01)) +
        labs(title = sens_spec_plot_title, x = paste(predictor_label, "Threshold")) +
        # Añadir líneas de eje manualmente usando annotate y límites calculados/fijos
        annotate("segment", x=-Inf, xend=Inf, y=y_min_sens, yend=y_min_sens, colour = "black", linewidth=0.5) + # Bottom
        annotate("segment", x=-Inf, xend=Inf, y=y_max_sens, yend=y_max_sens, colour = "black", linewidth=0.5) + # Top
        annotate("segment", x=x_min_sens, xend=x_min_sens, y=-Inf, yend=Inf, colour = "black", linewidth=0.5) + # Left
        annotate("segment", x=x_max_sens, xend=x_max_sens, y=-Inf, yend=Inf, colour = "black", linewidth=0.5) + # Right
        # Forzar límites y permitir anotaciones fuera
        coord_cartesian(xlim=finite_threshold_range, ylim = c(y_min_sens,y_max_sens), clip = "off", expand=FALSE) +
        # Aplicar tema base (sin líneas de eje propias)
        theme_pub_roc_manual_box +
        theme(legend.position = "bottom")

    # Guardar en múltiples formatos
    ggsave(file.path(output_roc_dir, paste0(sens_spec_filename_base, ".pdf")), plot = p_sens_spec, width = 6.5, height = 5, device=cairo_pdf)
    ggsave(file.path(output_roc_dir, paste0(sens_spec_filename_base, ".png")), plot = p_sens_spec, width = 6.5, height = 5, dpi = HIGH_RES_DPI, device="png")
    ggsave(file.path(output_roc_dir, paste0(sens_spec_filename_base, ".tif")), plot = p_sens_spec, width = 6.5, height = 5, dpi = HIGH_RES_DPI, device="tiff", compression="lzw")
    cat(sprintf("    Saved Sens/Spec plot (PDF, PNG, TIFF): %s\n", sens_spec_filename_base))
    cat("    (Note: Consider if Sens/Spec plot is needed for main text or supplementary)\n")

} # Fin del bucle for

# =======================================================================
# 3. MENSAJE FINAL
# =======================================================================
cat("\n--- ROC Analysis Script (Manual Box, Natural Threshold, Larger Font) Completed ---")
cat(sprintf("\nAnalysis performed on: %s", data_source_msg))
cat(sprintf("\nAll ROC outputs saved in '%s'.\n", normalizePath(output_roc_dir)))
cat(sprintf("Script finished at: %s\n", Sys.time()))