#!/usr/bin/env Rscript

# =============================================================================
# TF Activity Analysis Script
# 
# This script performs transcription factor activity analysis using:
# 1. DecoupleR with VIPER method
# 2. Differential TF activity analysis using limma
# 3. Visualization with volcano plots (static and interactive)
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(decoupleR)
  library(dorothea)
  library(limma)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(plotly)
  library(ggrepel)
})

# ---- DEBUG / DEV PARAMETERS (uncomment for testing) ----
if (TRUE) {
  opt <- list(
    tpm_file = "raw/all_samples.gene_tpm.tsv",
    samplesheet = "raw/samplesheet.csv",
    contrasts_file = "raw/contrasts.csv",
    output_dir = "results",
    logFC_threshold = 1,
    adj_pval_threshold = 1,
    min_genes = 5,
    dorothea_confidence = "A,B,C",
    viper_method = "scale",
    viper_cores = 1,
    plot_volcano = TRUE,
    interactive_volcano = TRUE,
    seed = 42
  )
}

# ---- define CLI options ----
if (!exists("opt")) {
  option_list <- list(
    make_option("--tpm_file", type="character", 
                help="Path to TPM file (gene_tpm.tsv)"),
    make_option("--samplesheet", type="character", 
                help="Path to samplesheet CSV file with phenotype information"),
    make_option("--contrasts_file", type="character", 
                help="Path to contrasts CSV file with id,variable,reference,target columns"),
    make_option("--output_dir", type="character", default="results",
                help="Output directory for results"),
    make_option("--logFC_threshold", type="numeric", default=1,
                help="Log2 fold change threshold for significance"),
    make_option("--adj_pval_threshold", type="numeric", default=0.05,
                help="Adjusted p-value threshold for significance"),
    make_option("--min_genes", type="integer", default=5,
                help="Minimum number of target genes for TF regulon"),
    make_option("--dorothea_confidence", type="character", default="A,B,C",
                help="Dorothea confidence levels (comma-separated: A,B,C,D)"),
    make_option("--viper_method", type="character", default="scale",
                help="VIPER method: scale, rank, or none"),
    make_option("--viper_cores", type="integer", default=1,
                help="Number of cores for VIPER"),
    make_option("--plot_volcano", type="logical", default=TRUE,
                help="Generate static volcano plot"),
    make_option("--interactive_volcano", type="logical", default=TRUE,
                help="Generate interactive volcano plot (HTML)"),
    make_option("--seed", type="integer", default=42,
                help="Random seed for reproducibility")
  )
  
  opt_parser <- OptionParser(option_list=option_list)
  opt <- parse_args(opt_parser)
}

# ---- set seed ----
set.seed(opt$seed)

# ---- log inputs ----
cat("\n========================================\n")
cat("TF Activity Analysis Pipeline\n")
cat("========================================\n\n")
cat("Input parameters:\n")
print(opt)
cat("\n")

# ---- check required arguments ----
required_args <- c("tpm_file", "samplesheet", "contrasts_file")
missing_args <- required_args[!required_args %in% names(opt) | 
                                sapply(opt[required_args], is.null) | 
                                opt[required_args] == ""]

if (length(missing_args) > 0) {
  if (!exists("opt_parser")) opt_parser <- NULL
  print_help(opt_parser)
  stop("Missing required arguments: ", paste(missing_args, collapse = ", "), "\n")
}

# ---- check input files exist ----
if (!file.exists(opt$tpm_file)) {
  stop("TPM file not found: ", opt$tpm_file)
}
if (!file.exists(opt$samplesheet)) {
  stop("Samplesheet not found: ", opt$samplesheet)
}
if (!file.exists(opt$contrasts_file)) {
  stop("Contrasts file not found: ", opt$contrasts_file)
}

# ---- create output directory ----
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
}

# ---- load data ----
cat("Loading TPM data...\n")
tpm <- fread(opt$tpm_file)

# ---- read samplesheet ----
sample_info <- read_csv(opt$samplesheet, show_col_types = FALSE)

# ---- read contrasts ----
contrasts <- read_csv(opt$contrasts_file, show_col_types = FALSE)

# Validate contrasts file structure
required_contrast_cols <- c("id", "variable", "reference", "target")
if (!all(required_contrast_cols %in% colnames(contrasts))) {
  stop("Contrasts file must contain columns: ", 
       paste(required_contrast_cols, collapse = ", "))
}

cat("Loaded contrasts:\n")
print(contrasts)

# Process all contrasts
for (i in 1:nrow(contrasts)) {
  contrast <- contrasts[i, ]
  variable_col <- contrast$variable
  reference_group <- contrast$reference
  case_group <- contrast$target
  contrast_id <- contrast$id
  
  cat("\n========================================\n")
  cat("Processing contrast:", contrast_id, "\n")
  cat("========================================\n")
  cat("  Variable:", variable_col, "\n")
  cat("  Reference:", reference_group, "\n")
  cat("  Target:", case_group, "\n\n")
  
  # Validate phenotype column exists
  if (!variable_col %in% colnames(sample_info)) {
    warning("Variable column '", variable_col, "' not found in samplesheet. Skipping contrast.")
    next
  }
  
  # ---- prepare expression matrix ----
  cat("Preparing expression matrix...\n")
  expr <- tpm |>
    select(gene_name, starts_with("R_")) |>
    filter(!is.na(gene_name), gene_name != "") |>
    distinct(gene_name, .keep_all = TRUE) |>
    column_to_rownames("gene_name") |>
    as.matrix() + 1 |> 
    log2()
  
  # ---- get sample names from expression matrix ----
  expr_samples <- colnames(expr)
  
  # ---- match samples to phenotype ----
  sample_info_filtered <- sample_info |>
    filter(sample %in% expr_samples) |>
    mutate(
      group = !!sym(variable_col),
      group = factor(group, levels = c(reference_group, case_group))
    )
  
  # Check if we have both groups
  if (length(unique(sample_info_filtered$group)) < 2) {
    warning("Need at least two groups for contrast ", contrast_id, 
            ". Found: ", paste(unique(sample_info_filtered$group), collapse = ", "))
    next
  }
  
  cat("Sample grouping:\n")
  print(table(sample_info_filtered$group))
  
  # ---- filter expression matrix to match samples ----
  expr_filtered <- expr[, sample_info_filtered$sample, drop = FALSE]
  
  # ---- load Dorothea regulons (only once) ----
  if (!exists("regulons")) {
    cat("Loading Dorothea regulons...\n")
    data(dorothea_hs)
    
    # Parse confidence levels
    confidence_levels <- str_split(opt$dorothea_confidence, ",")[[1]] |>
      str_trim()
    
    regulons <- dorothea_hs |>
      filter(confidence %in% confidence_levels)
    
    cat("Loaded ", nrow(regulons), " regulon entries with confidence levels: ", 
        paste(confidence_levels, collapse = ", "), "\n")
  }
  
  # ---- VIPER ----
  cat("Running VIPER TF activity analysis...\n")
  tf_activities <- run_viper(
    input = expr_filtered,
    regulons = regulons,
    options = list(
      method = opt$viper_method,
      minsize = opt$min_genes,
      eset.filter = FALSE,
      cores = opt$viper_cores,
      verbose = FALSE
    )
  )
  
  cat("VIPER TF activity matrix: ", nrow(tf_activities), " TFs x ", 
      ncol(tf_activities), " samples\n")
  
  # ---- Differential TF Activity ----
  cat("Running differential TF activity analysis...\n")
  
  # Ensure groups are in correct order
  sample_info_filtered$group <- factor(
    sample_info_filtered$group,
    levels = c(reference_group, case_group)
  )
  
  # Create design matrix
  design <- model.matrix(~ group, data = sample_info_filtered)
  
  # Fit linear model
  fit <- lmFit(tf_activities, design)
  fit <- eBayes(fit)
  
  # Extract results for case vs reference
  coef_name <- paste0("group", case_group)
  if (!coef_name %in% colnames(design)) {
    warning("Coefficient '", coef_name, "' not found in design matrix for contrast ", contrast_id)
    next
  }
  
  tf_stats <- topTable(
    fit,
    coef = coef_name,
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )
  
  # ---- create results ----
  tf_results <- tf_stats |>
    rownames_to_column("TF") |>
    select(
      TF,
      logFC,
      AveExpr,
      t,
      P.Value,
      adj.P.Val
    )
  
  # ---- filter significant TFs ----
  sig_tf <- tf_results |>
    filter(adj.P.Val < opt$adj_pval_threshold & abs(logFC) > opt$logFC_threshold) |>
    arrange(desc(abs(logFC)))
  
  cat("\nSignificant TFs (adj.P.Val < ", opt$adj_pval_threshold, 
      " and |logFC| > ", opt$logFC_threshold, "): ", nrow(sig_tf), "\n")
  
  if (nrow(sig_tf) > 0) {
    cat("Top 10 significant TFs:\n")
    print(head(sig_tf, 10))
  }
  
  # ---- save results ----
  cat("Saving results...\n")
  
  # Save full results
  write_tsv(
    tf_results,
    file.path(opt$output_dir, paste0(contrast_id, "_tf_results.tsv"))
  )
  
  # Save significant TFs separately
  if (nrow(sig_tf) > 0) {
    write_tsv(
      sig_tf,
      file.path(opt$output_dir, paste0(contrast_id, "_significant_tfs.tsv"))
    )
  }
  
  # ---- static volcano plot ----
  if (opt$plot_volcano && nrow(tf_results) > 0) {
    cat("Generating static volcano plot...\n")
    
    # Add significance labels
    tf_plot <- tf_results |>
      mutate(
        significance = case_when(
          adj.P.Val < opt$adj_pval_threshold & abs(logFC) > opt$logFC_threshold ~ "Significant",
          adj.P.Val < opt$adj_pval_threshold & abs(logFC) <= opt$logFC_threshold ~ "Significant (small FC)",
          TRUE ~ "Not significant"
        ),
        label = if_else(significance == "Significant", TF, "")
      )
    
    # Get top TFs for labeling (top 10 by adjusted p-value)
    top_tfs <- tf_plot |>
      filter(significance == "Significant") |>
      arrange(adj.P.Val) |>
      head(10) |>
      pull(TF)
    
    tf_plot <- tf_plot |>
      mutate(label = if_else(TF %in% top_tfs, TF, ""))
    
    volcano_plot <- ggplot(tf_plot, aes(x = logFC, y = -log10(adj.P.Val))) +
      geom_point(aes(color = significance), alpha = 0.7, size = 2) +
      scale_color_manual(
        values = c(
          "Significant" = "red",
          "Significant (small FC)" = "orange",
          "Not significant" = "gray70"
        )
      ) +
      geom_vline(xintercept = c(-opt$logFC_threshold, opt$logFC_threshold), 
                 linetype = "dashed", alpha = 0.5, color = "darkred") +
      geom_hline(yintercept = -log10(opt$adj_pval_threshold), 
                 linetype = "dashed", alpha = 0.5, color = "darkred") +
      labs(
        title = paste("TF Activity Volcano Plot -", contrast_id),
        x = "Log2 Fold Change",
        y = "-Log10 Adjusted P-value",
        color = "Significance"
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12)
      )
    
    # Add labels for top TFs
    if (length(top_tfs) > 0) {
      volcano_plot <- volcano_plot +
        geom_text_repel(
          data = tf_plot |> filter(TF %in% top_tfs),
          aes(label = TF),
          size = 3,
          max.overlaps = 15,
          box.padding = 0.5,
          point.padding = 0.2
        )
    }
    
    # Save static plot
    ggsave(
      file.path(opt$output_dir, paste0(contrast_id, "_volcano.pdf")),
      volcano_plot,
      width = 10,
      height = 8
    )
    
    ggsave(
      file.path(opt$output_dir, paste0(contrast_id, "_volcano.png")),
      volcano_plot,
      width = 10,
      height = 8,
      dpi = 300
    )
    
    cat("Static volcano plot saved to: ", 
        file.path(opt$output_dir, paste0(contrast_id, "_volcano.pdf")), "\n")
  }
  
  # ---- interactive volcano plot ----
  if (opt$interactive_volcano && nrow(tf_results) > 0) {
    cat("Generating interactive volcano plot...\n")
    
    # Prepare data for interactive plot
    tf_plotly <- tf_results |>
      mutate(
        significance = case_when(
          adj.P.Val < opt$adj_pval_threshold & abs(logFC) > opt$logFC_threshold ~ "Significant",
          adj.P.Val < opt$adj_pval_threshold & abs(logFC) <= opt$logFC_threshold ~ "Significant (small FC)",
          TRUE ~ "Not significant"
        ),
        # Add hover text
        hover_text = paste(
          "TF:", TF,
          "<br>logFC:", round(logFC, 3),
          "<br>adj.P.Val:", format(adj.P.Val, scientific = TRUE, digits = 3),
          "<br>AveExpr:", round(AveExpr, 2),
          "<br>t-statistic:", round(t, 2)
        )
      )
    
    # Create interactive plot
    interactive_plot <- plot_ly(
      data = tf_plotly,
      x = ~logFC,
      y = ~-log10(adj.P.Val),
      type = 'scatter',
      mode = 'markers',
      color = ~significance,
      colors = c(
        "Significant" = "red",
        "Significant (small FC)" = "orange",
        "Not significant" = "gray70"
      ),
      text = ~hover_text,
      hoverinfo = 'text',
      marker = list(
        size = 8,
        opacity = 0.7
      ),
      showlegend = TRUE
    ) |>
      layout(
        title = list(
          text = paste("TF Activity Volcano Plot -", contrast_id),
          x = 0.5
        ),
        xaxis = list(
          title = "Log2 Fold Change",
          zeroline = TRUE,
          zerolinecolor = 'lightgray',
          gridcolor = 'lightgray'
        ),
        yaxis = list(
          title = "-Log10 Adjusted P-value",
          zeroline = TRUE,
          zerolinecolor = 'lightgray',
          gridcolor = 'lightgray'
        ),
        shapes = list(
          # Vertical lines at logFC thresholds
          list(
            type = "line",
            x0 = -opt$logFC_threshold,
            x1 = -opt$logFC_threshold,
            y0 = 0,
            y1 = max(-log10(tf_plotly$adj.P.Val)) * 1.1,
            line = list(color = "darkred", dash = "dash", width = 1)
          ),
          list(
            type = "line",
            x0 = opt$logFC_threshold,
            x1 = opt$logFC_threshold,
            y0 = 0,
            y1 = max(-log10(tf_plotly$adj.P.Val)) * 1.1,
            line = list(color = "darkred", dash = "dash", width = 1)
          ),
          # Horizontal line at p-value threshold
          list(
            type = "line",
            x0 = min(tf_plotly$logFC) * 1.1,
            x1 = max(tf_plotly$logFC) * 1.1,
            y0 = -log10(opt$adj_pval_threshold),
            y1 = -log10(opt$adj_pval_threshold),
            line = list(color = "darkred", dash = "dash", width = 1)
          )
        ),
        legend = list(
          y = 0.5,
          x = 0.02,
          bgcolor = 'rgba(255, 255, 255, 0.8)'
        ),
        hoverlabel = list(
          bgcolor = "white",
          font = list(size = 12)
        )
      )
    
    # Save interactive plot as HTML
    htmlwidgets::saveWidget(
      interactive_plot,
      file.path(opt$output_dir, paste0(contrast_id, "_volcano_interactive.html")),
      selfcontained = TRUE
    )
    
    cat("Interactive volcano plot saved to: ", 
        file.path(opt$output_dir, paste0(contrast_id, "_volcano_interactive.html")), "\n")
  }
  
  # Save VIPER activity matrix (only once, but keep for each contrast)
  write_tsv(
    tf_activities |> as.data.frame() |> rownames_to_column("TF"),
    file.path(opt$output_dir, paste0(contrast_id, "_viper_activity_matrix.tsv"))
  )
}

# ---- generate summary report ----
cat("\nGenerating summary report...\n")

summary_report <- file.path(opt$output_dir, "analysis_summary.txt")
sink(summary_report)

cat("========================================\n")
cat("TF Activity Analysis Summary\n")
cat("========================================\n\n")
cat("Analysis date:", Sys.time(), "\n\n")
cat("Input files:\n")
cat("  TPM:", opt$tpm_file, "\n")
cat("  Samplesheet:", opt$samplesheet, "\n")
cat("  Contrasts:", opt$contrasts_file, "\n\n")
cat("Parameters:\n")
cat("  logFC threshold:", opt$logFC_threshold, "\n")
cat("  Adjusted P-value threshold:", opt$adj_pval_threshold, "\n")
cat("  Dorothea confidence levels:", opt$dorothea_confidence, "\n")
cat("  Min genes per regulon:", opt$min_genes, "\n\n")

cat("Processed contrasts:\n")
for (i in 1:nrow(contrasts)) {
  contrast <- contrasts[i, ]
  cat("  -", contrast$id, ":", contrast$reference, "vs", contrast$target, "\n")
}
cat("\n")

cat("Output files:\n")
for (i in 1:nrow(contrasts)) {
  contrast_id <- contrasts[i, ]$id
  cat("  Contrast:", contrast_id, "\n")
  cat("    -", file.path(opt$output_dir, paste0(contrast_id, "_tf_results.tsv")), "\n")
  cat("    -", file.path(opt$output_dir, paste0(contrast_id, "_significant_tfs.tsv")), "\n")
  cat("    -", file.path(opt$output_dir, paste0(contrast_id, "_viper_activity_matrix.tsv")), "\n")
  if (opt$plot_volcano) {
    cat("    -", file.path(opt$output_dir, paste0(contrast_id, "_volcano.pdf")), "\n")
    cat("    -", file.path(opt$output_dir, paste0(contrast_id, "_volcano.png")), "\n")
  }
  if (opt$interactive_volcano) {
    cat("    -", file.path(opt$output_dir, paste0(contrast_id, "_volcano_interactive.html")), "\n")
  }
}

sink()

cat("========================================\n")
cat("Analysis complete!\n")
cat("Summary report saved to:", summary_report, "\n")
cat("========================================\n")