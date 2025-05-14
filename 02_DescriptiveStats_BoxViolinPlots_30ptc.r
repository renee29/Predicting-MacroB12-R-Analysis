# ===========================================================
#        Script R para Visualización Comparativa Post-Imputación
#        (Estilo BoxViolin con Datos MICE - Basado en 1 Imputación)
# ===========================================================

# 0. CONFIGURACIÓN E IMPORTACIÓN
# ===========================================================
# --- Cargar Paquetes ---
required_packages <- c("data.table", "tidyverse", "ggplot2", "here", "rstudioapi",
                       "mice",      # Para manejar objeto mids
                       "ggthemes", "RColorBrewer", "ggrepel", "ggdist",
                       "knitr", "kableExtra", "skimr", "ggpubr", "paletteer",
                       "dplyr", "stats")

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing package: %s\n", pkg))
    tryCatch(install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org/"),
             error = function(e){ cat(sprintf("  Failed to install %s: %s\n", pkg, e$message))})
  }
}
sapply(required_packages, install_if_missing)

library(data.table)
library(tidyverse)
library(ggplot2)
library(here)
library(rstudioapi)
library(mice)
library(ggthemes)
library(RColorBrewer)
library(ggrepel)
library(ggdist)
library(knitr)
library(kableExtra)
library(skimr)
library(ggpubr)
library(paletteer)
library(dplyr)
library(stats)

# -----------------------------------------------------------
# --- Directorios ---
script_dir <- tryCatch({ dirname(rstudioapi::getSourceEditorContext()$path) }, error = function(e) { getwd() })
# Directorio donde se guardaron los resultados del script MICE
input_data_prep_dir <- file.path(script_dir, "R_Outputs_Data_Prep_New_30pct")
# Nuevo directorio para los resultados de este script de visualización
output_plot_dir <- file.path(script_dir, "BoxViolinPlots_30ptc")
if (!dir.exists(output_plot_dir)) { dir.create(output_plot_dir) } else { cat("Plot output directory exists:", normalizePath(output_plot_dir), "\n") }

# --- Constantes ---
RANDOM_STATE <- 42; set.seed(RANDOM_STATE)
IMPUTATION_TO_PLOT <- 1 # Qué dataset imputado usar para los gráficos

# --- Cargar Mapeo de Etiquetas (Keys deben ser los nombres INTERNOS/make.names) ---
# Asegúrate que este mapa esté actualizado desde el script MICE/Análisis
LABEL_MAP <- list(
  'edad' = 'Age (years)',
  'sexo' = 'Sex',
  'procedencia.' = 'Source',
  'vb12_pre_peg' = 'VB12 Pre-PEG (pg/mL)',
  'vb12_post_peg' = 'VB12 Post-PEG (pg/mL)',
  'Recovery_pct' = 'PEG Recovery (%)',
  'MacroB12_Positive' = 'MacroB12 Status',
  'pcr' = 'CRP (mg/L)',
  'fol' = 'Serum Folate (ng/mL)',
  'hb' = 'Hemoglobin (g/dL)',
  'vcm' = 'MCV (fL)',
  'ade' = 'RDW (%)',
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
  'Male' = 'Male', 'Female' = 'Female',
  'Ambulatorio' = 'Ambulatory', 'Hospital' = 'Hospital', 'Centro_Externo' = 'External Center',
  'Negative'='Negative', 'Positive'='Positive'
)
# Cargar mapeo de nombres finales si existe (para reverse lookup si es necesario)
final_name_map_file <- file.path(input_data_prep_dir, "column_rename_map_mk_to_eng.rds")
internal_to_final_map <- NULL
if(file.exists(final_name_map_file)){
    internal_to_final_map <- readRDS(final_name_map_file)
    # Invertir para tener final -> internal
    final_to_internal_map <- setNames(names(internal_to_final_map), unlist(internal_to_final_map))
} else {
    cat("Warning: Final name map file not found. Label lookup might be limited if using final CSV.\n")
    final_to_internal_map <- list() # Crear lista vacía
}


# --- Helper Functions ---
# Ajustado para usar LABEL_MAP con nombres internos/make.names
# y manejar posible reverse lookup si se usan nombres finales
get_label <- function(name, default_fmt=TRUE) {
    name_str <- as.character(name);

    # Intentar buscar el nombre interno correspondiente si estamos usando nombres finales
    internal_name <- name_str # Asumir que ya es el nombre interno
    if(exists("col_names_type") && col_names_type == "final"){
       internal_lookup <- final_to_internal_map[[name_str]]
       if(!is.null(internal_lookup)) internal_name <- internal_lookup
       # else: no se encontró mapeo, se usará el nombre final para buscar en LABEL_MAP
    }

    # Buscar etiqueta usando el nombre (potencialmente interno)
    label <- LABEL_MAP[[internal_name]];
    if (!is.null(label)) return(label);

    # Fallback simple si no se encuentra etiqueta
    clean_name <- gsub("`", "", name_str) # Usar el nombre original pasado si falla el lookup
    if (default_fmt) return(tools::toTitleCase(gsub("[._]", " ", clean_name))) else return(clean_name)
}
format_p_value <- function(p) { if (is.na(p)) { return("NA") } else if (p < 0.0001) { return(sprintf("%.1e", p)) } else if (p < 0.001) { return("<0.001") } else if (p < 0.01) { return(sprintf("%.3f", p)) } else { return(sprintf("%.2f", p)) } }
add_stars <- function(p) { if (is.na(p) || !is.numeric(p)) return('') ; if (p < 0.001) return('***') ; if (p < 0.01) return('**') ; if (p < 0.05) return('*') ; return('') }
sanitize_filename <- function(name) { name_sanitized <- gsub('[<>:"/\\|?*%[:cntrl:]]', '_', as.character(name)); name_sanitized <- gsub('[[:punct:]]', '_', name_sanitized); name_sanitized <- gsub('\\s+', '_', name_sanitized); name_sanitized <- gsub('_+', '_', name_sanitized); name_sanitized <- gsub('^_|_$', '', name_sanitized); name_sanitized <- substr(name_sanitized, 1, 100); return(name_sanitized) }

# =======================================================================
# 1. CARGAR DATOS PROCESADOS (Objeto MICE o CSV Final)
# =======================================================================
cat("\n--- 1. Loading Processed Data for Visualization ---\n")

mice_rds_file <- file.path(input_data_prep_dir, "mice_imputation_object_internalNames.rds")
complete_case_csv_file <- file.path(input_data_prep_dir, paste0("analysis_data_final_EngNames_", IMPUTATION_TO_PLOT, ".csv"))

mids_data <- NULL
plot_df <- NULL
data_source <- NULL
col_names_type <- NULL

# Priorizar MICE object
if (file.exists(mice_rds_file)) {
    tryCatch({
        mids_data <- readRDS(mice_rds_file)
        cat("Successfully loaded MICE object (.rds).\n")
        if (IMPUTATION_TO_PLOT > 0 && IMPUTATION_TO_PLOT <= mids_data$m) {
             plot_df <- mice::complete(mids_data, action = IMPUTATION_TO_PLOT)
             cat(sprintf("Extracted imputed dataset #%d for plotting.\n", IMPUTATION_TO_PLOT))
        } else {
             cat(sprintf("Warning: Imputation number %d not valid (m=%d). Using imputation #1.\n", IMPUTATION_TO_PLOT, mids_data$m))
             plot_df <- mice::complete(mids_data, action = 1)
             IMPUTATION_TO_PLOT <- 1 # Resetear a 1 si no era válido
        }
        data_source <- "mice_imputed"
        col_names_type <- "internal"

        # Asegurar tipos correctos post-extracción de MICE
        target_col_internal <- make.names("MacroB12_Positive") # Nombre interno esperado
        if(target_col_internal %in% names(plot_df)){
             if(!is.factor(plot_df[[target_col_internal]])) plot_df[[target_col_internal]] <- factor(plot_df[[target_col_internal]], levels=c(0,1), labels=c("Negative", "Positive")) # Asume 0/1 numérico antes
             plot_df[[target_col_internal]] <- factor(plot_df[[target_col_internal]], levels=c("Negative", "Positive")) # Asegurar niveles
        } else { stop("Target column (internal name) not found in extracted MICE data.") }

        internal_diag_cols <- names(plot_df)[startsWith(names(plot_df), "Diag_")]
        for(col in internal_diag_cols){
             if(is.numeric(plot_df[[col]])) { plot_df[[col]] <- factor(round(plot_df[[col]]), levels=c(0,1)) }
             else if (!is.factor(plot_df[[col]])) { plot_df[[col]] <- factor(plot_df[[col]], levels=c(0,1)) } # Si era character 0/1
        }
         # Asegurar otros factores categóricos
         sex_col_internal <- make.names("sexo"); proc_col_internal <- make.names("procedencia.")
         if(sex_col_internal %in% names(plot_df) && !is.factor(plot_df[[sex_col_internal]])) plot_df[[sex_col_internal]] <- factor(plot_df[[sex_col_internal]])
         if(proc_col_internal %in% names(plot_df) && !is.factor(plot_df[[proc_col_internal]])) plot_df[[proc_col_internal]] <- factor(plot_df[[proc_col_internal]])
         # Restablecer ref levels
         if(sex_col_internal %in% names(plot_df) && "Female" %in% levels(plot_df[[sex_col_internal]])) plot_df[[sex_col_internal]] <- relevel(plot_df[[sex_col_internal]], ref="Female")
         if(proc_col_internal %in% names(plot_df) && "Ambulatory" %in% levels(plot_df[[proc_col_internal]])) plot_df[[proc_col_internal]] <- relevel(plot_df[[proc_col_internal]], ref="Ambulatory")


    }, error = function(e) {
        cat("Error loading MICE object or extracting dataset:", e$message, "\n")
        mids_data <<- NULL; plot_df <<- NULL
    })
}

# Fallback al CSV si falla MICE
if (is.null(plot_df)) {
    cat("Attempting to load final name CSV as fallback.\n")
    if (file.exists(complete_case_csv_file)) {
        tryCatch({
            plot_df <- fread(complete_case_csv_file, na.strings = c("", "NA"))
            plot_df <- as.data.frame(plot_df)
            # Asegurar tipos desde CSV
            target_col_final <- "MacroB12_Status"; sex_col_final <- "Sex"; source_col_final <- "PatientSource"; diag_cols_final <- names(plot_df)[startsWith(names(plot_df), "Diag_")]
            if (target_col_final %in% names(plot_df)) plot_df[[target_col_final]] <- factor(plot_df[[target_col_final]], levels=c("Negative", "Positive"))
            if (sex_col_final %in% names(plot_df)) plot_df[[sex_col_final]] <- factor(plot_df[[sex_col_final]], levels = c("Female", "Male"))
            if (source_col_final %in% names(plot_df)) plot_df[[source_col_final]] <- factor(plot_df[[source_col_final]], levels = c("Ambulatory", "External Center", "Hospital"))
            for(col in diag_cols_final){ if(!is.factor(plot_df[[col]])) plot_df[[col]] <- factor(plot_df[[col]], levels=c(0,1)) }
            # Ref levels
            if (sex_col_final %in% names(plot_df)) plot_df[[sex_col_final]] <- relevel(plot_df[[sex_col_final]], ref="Female")
            if (source_col_final %in% names(plot_df)) plot_df[[source_col_final]] <- relevel(plot_df[[source_col_final]], ref="Ambulatory")
            # Numéricos
            num_cols_final <- c('Age', 'VB12_PrePEG', 'CRP', 'Folate', 'Hemoglobin', 'MCV', 'RDW', 'PEG_Recovery_Percent')
            for(col in intersect(num_cols_final, names(plot_df))){ if(!is.numeric(plot_df[[col]])) plot_df[[col]] <- suppressWarnings(as.numeric(plot_df[[col]])) }

            cat("Successfully loaded data with final names from CSV:", basename(complete_case_csv_file), "\n")
            data_source <- "csv_final_names" # Cambiado nombre fuente
            col_names_type <- "final"
        }, error = function(e) {
            cat("Error loading final name CSV:", e$message, "\n"); plot_df <<- NULL
        })
    } else {
        cat("Final name CSV file not found:", complete_case_csv_file, "\n")
    }
}

if (is.null(plot_df)) {
    stop("Could not load data for plotting. Check paths and output from previous script.")
}

cat("Data loaded for plotting. Dimensions:", dim(plot_df), "\n")
cat("Using column names type:", col_names_type, "\n")
cat("Column types:\n"); print(sapply(plot_df, function(x) paste(class(x), collapse=", ")))


# =======================================================================
# 2. DEFINIR VARIABLES Y GENERAR GRÁFICOS COMPARATIVOS
# =======================================================================
cat(sprintf("\n--- 2. Generating Comparative Plots (Based on %s) ---\n", data_source))

# Determinar nombres de columnas a usar
if (col_names_type == "internal") {
    target_col_plot <- make.names("MacroB12_Positive") # Nombre interno
    all_num_features <- names(plot_df)[sapply(plot_df, is.numeric)]
    # Excluir columnas numéricas que no son variables clínicas de interés (si las hubiera)
    vars_to_exclude <- c(make.names("vb12_post_peg")) # Añadir otros si es necesario
    num_features_plot <- setdiff(all_num_features, vars_to_exclude)
} else { # col_names_type == "final"
    target_col_plot <- "MacroB12_Status" # Nombre final
    all_num_features <- names(plot_df)[sapply(plot_df, is.numeric)]
    vars_to_exclude <- c("VB12_PostPEG") # Usar nombre final
    num_features_plot <- setdiff(all_num_features, vars_to_exclude)
}

vars_to_plot <- num_features_plot

cat("Numeric variables selected for plotting:", paste(vars_to_plot, collapse=", "), "\n")
cat("Target variable for grouping:", target_col_plot, "\n")

# --- Paleta de Colores 'ggthemes::calc' ---
if (!target_col_plot %in% names(plot_df) || !is.factor(plot_df[[target_col_plot]])) {
    stop(paste("Target column '", target_col_plot, "' not found or is not a factor in the plotting data."))
}
target_levels <- levels(plot_df[[target_col_plot]])
num_colors_needed <- length(target_levels)
if(num_colors_needed < 2) stop("Target variable needs at least 2 levels for color mapping.")
plot_palette_calc <- paletteer_d("ggthemes::calc", n = max(num_colors_needed, 3))[1:num_colors_needed]
names(plot_palette_calc) <- target_levels
cat("Using color palette 'ggthemes::calc':\n"); print(plot_palette_calc)

# --- Tema ggplot base ---
theme_pub_final <- theme_classic(base_size = 11, base_family = "sans") +
  theme( plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1.5)),
         plot.subtitle = element_text(hjust = 0.5, size = rel(1.2), color="black"),
         axis.title = element_text(face = "bold", size = rel(1.5)),
         axis.text = element_text(size = rel(1.5), color = "black"),
         axis.line = element_line(colour = "black", linewidth = 0.5),
         axis.ticks = element_line(colour = "black", linewidth = 0.5),
         panel.grid.major.y = element_line(colour = "grey90", linetype = "dashed", linewidth = 0.3),
         panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_rect(fill = "white", colour = NA),
         plot.background = element_rect(fill = "white", colour = NA),
         plot.margin = ggplot2::margin(10, 15, 10, 10, unit = "pt"),
         legend.position = "none" ,
         panel.border = element_rect(colour = "black", fill=NA, linewidth=0.7)
         )

# --- Bucle de Gráficos ---
for (var_y in vars_to_plot) {

  # Obtener etiqueta bonita (intenta buscar nombre interno si es necesario)
  var_y_label <- get_label(var_y) # get_label ahora maneja la posible búsqueda inversa

  target_label <- get_label(target_col_plot)
  plot_title <- paste(var_y_label, "by", target_label)
  subtitle_info <- ifelse(data_source == "mice_imputed", sprintf("(Based on Imputation #%d)", IMPUTATION_TO_PLOT), "(Based on Single Dataset)")
  safe_filename <- paste0("Fig_BoxViolin_", sanitize_filename(var_y_label), "_by_Target_", data_source,".jpg") # Nombre de archivo más descriptivo
  output_path <- file.path(output_plot_dir, safe_filename)

  cat(sprintf("  Generating plot for %s...\n", var_y_label))

  # Filtrar NAs para la variable Y actual
  plot_data_current <- plot_df %>% filter(!is.na(!!sym(var_y)))
  if(nrow(plot_data_current) < 2 || n_distinct(plot_data_current[[target_col_plot]], na.rm=TRUE) < 2) {
      cat(sprintf("    Skipping %s: Insufficient data or groups.\n", var_y_label)); next
  }

  # Calcular estadísticas para etiquetas
  plot_summary <- plot_data_current %>%
    group_by(!!sym(target_col_plot)) %>%
    summarise( Median = median(!!sym(var_y), na.rm = TRUE), Q1 = quantile(!!sym(var_y), 0.25, na.rm = TRUE), Q3 = quantile(!!sym(var_y), 0.75, na.rm = TRUE), N = n(), Mean = mean(!!sym(var_y), na.rm=TRUE), .groups = 'drop') %>%
    mutate(LabelIQR = sprintf("%.1f [%.1f-%.1f]", Median, Q1, Q3), LabelN = sprintf("(n = %d)", N))
  y_range <- range(plot_data_current[[var_y]], na.rm = TRUE)
  y_padding <- diff(y_range) * 0.08
  y_min_label_n <- y_range[1] - y_padding
  plot_summary <- plot_summary %>% mutate(LabelY_IQR = Median, LabelY_N = y_min_label_n)

  # Calcular p-valor (Wilcoxon)
  p_val_mw <- NA_real_; test_statistic_mw <- NA_real_
  test_result <- tryCatch({ wilcox.test(as.formula(paste0("`",var_y,"` ~ `", target_col_plot, "`")), data = plot_data_current)}, error = function(e) NULL)
  if(!is.null(test_result)) { p_val_mw <- test_result$p.value; test_statistic_mw <- test_result$statistic }
  stats_subtitle_text <- sprintf("Wilcoxon W=%.1f, p=%s %s, n=%d", ifelse(!is.na(test_statistic_mw), test_statistic_mw, NA), format_p_value(p_val_mw), add_stars(p_val_mw), nrow(plot_data_current))
  #stats_subtitle_text <- sprintf("Wilcoxon W=%.1f, p=%s %s, n=%d\n%s", ifelse(!is.na(test_statistic_mw), test_statistic_mw, NA), format_p_value(p_val_mw), add_stars(p_val_mw), nrow(plot_data_current), subtitle_info)
  #stats_subtitle_text <- sprintf("Wilcoxon W=%.1f, p=%s %s, n=%d\n%s", ifelse(!is.na(test_statistic_mw), test_statistic_mw, NA), format_p_value(p_val_mw), add_stars(p_val_mw), nrow(plot_data_current), subtitle_info)

  # Calcular p-valor (t-test) - opcional
  p_val_tt <- NA_real_; test_statistic_tt <- NA_real_
  test_result_tt <- tryCatch({ t.test(as.formula(paste0("`",var_y,"` ~ `", target_col_plot, "`")), data = plot_data_current)}, error = function(e) NULL)
  if(!is.null(test_result_tt)) { p_val_tt <- test_result_tt$p.value; test_statistic_tt <- test_result_tt$statistic }
  # Podrías añadirlo al subtítulo si quieres: sprintf(" | t=%.2f, p=%s", test_statistic_tt, format_p_value(p_val_tt))

  tryCatch({
    # Crear gráfico base con ggplot2
    p <- ggplot(plot_data_current, aes(x = !!sym(target_col_plot), y = !!sym(var_y))) +

      # Capa 1: Violín
      geom_violin(aes(color = !!sym(target_col_plot)), trim = TRUE, alpha = 0.25, linewidth = 0.8, scale = "width", draw_quantiles = NULL) +

      # Capa 2: Puntos Jittered
       geom_jitter(aes(color = !!sym(target_col_plot)), width = 0.20, height = 0, alpha = 0.4, size = 1.75, shape = 16) +

      # Capa 3: Boxplot
      geom_boxplot(width = 0.12, outlier.shape = NA, coef = 1.5, fill = "white", color = "black", alpha = 0.7, linewidth=0.7) +

      # Capa 4: Punto para la media
      stat_summary(aes(color = !!sym(target_col_plot)), fun = mean, geom = "point", shape = 19, size = 3.5) +

     # Capa 5: Etiquetas Mediana [IQR]
      ggrepel::geom_label_repel(
                data = plot_summary,
                aes(label = LabelIQR, y = LabelY_IQR, color = !!sym(target_col_plot), group = !!sym(target_col_plot)),
                size = 4.5, fontface="italic", fill = "white", alpha = 0.85,
                label.padding = unit(0.15, "lines"), label.r = unit(0.1, "lines"), label.size = 0.4,
                box.padding = 0.3, nudge_x = -0.35, direction = "y",
                segment.size = 0.3, segment.color = "grey50", min.segment.length = 0.1, point.padding = 0.1,
                max.overlaps = Inf, seed = RANDOM_STATE, show.legend = FALSE
                ) +

      # Capa 6: Etiqueta N
      geom_text( data = plot_summary, aes(label = LabelN, y = LabelY_N, color = !!sym(target_col_plot)),
                 size = 4.5, vjust = 1.1, show.legend = FALSE) +

      # Aplicar escalas, tema y etiquetas
      scale_color_manual(values = plot_palette_calc) +
      labs(title = plot_title, subtitle = stats_subtitle_text, x = target_label, y = var_y_label) +
      theme_pub_final +
      guides(color = "none", fill = "none") +
      coord_cartesian(clip="off") # Permitir que las etiquetas N salgan un poco

    # Guardar
    ggsave(filename = output_path, plot = p, width = 6, height = 5.5, dpi = 500, device = 'jpeg', quality = 95)
    cat(sprintf("  Saved: %s\n", basename(output_path)))

  }, error = function(e) {
    cat(sprintf("    ERROR generating plot for %s: %s\n", var_y_label, e$message))
  })
}

cat("\nComparative plot generation complete.\n")
cat(sprintf("Plots saved in '%s'.\n", normalizePath(output_plot_dir)))
cat(sprintf("Script finished at: %s\n", Sys.time()))