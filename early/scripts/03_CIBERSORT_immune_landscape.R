# =============================================================================
# 03_CIBERSORT_immune_landscape.R
# PD–OSCC transcriptomic study | GitHub reproducibility CANDIDATE
#
# RUN FROM THE PROJECT ROOT:
#   First set the working directory to this project folder in RStudio.
#   source("scripts/03_CIBERSORT_immune_landscape.R")
#
# PURPOSE:
# Audit archived CIBERSORTx LM22 outputs (Job29: GSE10334 / PD;
# Job30: GSE30784 / OSCC), align GSM sample identifiers, perform historical
# unpaired Wilcoxon comparisons + BH correction, and run a separately labeled
# PD within-patient paired sensitivity analysis.
#
# SCIENTIFIC / REPRODUCIBILITY NOTES:
# - This script DOES NOT run CIBERSORTx or regenerate its immune fractions.
#   It starts from ORIGINAL archived Job29/Job30 results.
# - Both historical jobs: LM22.update-gene-symbols.txt, QN=TRUE, permutations=1
#   (verified against archived runtime logs; 1 permutation is insufficient
#   for reliable permutation-based CIBERSORTx P-value estimation).
# - The archived mixture TXT inputs are NOT identical to the 20,848-row
#   representative-gene matrices created by revised Script 01:
#   PD input: 28,655 rows; OSCC input: 28,704 rows; 21,755 distinct symbols
#   in each, with duplicated symbols. Do NOT silently substitute Script 01's
#   alternative matrices as if they reproduced the historical mixture files.
# - Fractions are estimated relative cell fractions, not direct cell counts.
# - Historical unpaired comparisons treat biopsies as independent. PD paired
#   sensitivity uses 64 patients with BOTH PD and Control; 26 PD-only patients
#   are not included in that paired analysis.
# - These analyses do not prove macrophage phenotype transitions or causality.
# - Outputs are regenerated under results/03_CIBERSORT; originals under
#   data/cibersortx_historical are NEVER modified by this script.
# =============================================================================

# --- Setup ---------------------------------------------------------------
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source_dir <- file.path(root, "data", "cibersortx_historical")
output_dir <- file.path(root, "results", "03_CIBERSORT")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metadata_file <- file.path(root, "results", "01_preprocessing",
                           "01_preprocessed_esets.rds")
pd_result_file <- file.path(source_dir, "CIBERSORTx_Job29_Results.csv")
oscc_result_file <- file.path(source_dir, "CIBERSORTx_Job30_Results.csv")
pd_input_file <- file.path(source_dir, "CIBERSORT_input_PD_linear.txt")
oscc_input_file <- file.path(source_dir, "CIBERSORT_input_OSCC_linear.txt")
pd_log_file <- file.path(source_dir, "CIBERSORTx_Job29_runtime_log.txt")
oscc_log_file <- file.path(source_dir, "CIBERSORTx_Job30_runtime_log.txt")

required_files <- c(metadata_file, pd_result_file, oscc_result_file,
                    pd_input_file, oscc_input_file,
                    pd_log_file, oscc_log_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Missing prerequisite file(s):\n", paste(missing_files, collapse = "\n"))
}

# These are the canonical LM22 cell names expected in the original outputs.
immune_cols <- c(
  "B cells naive", "B cells memory", "Plasma cells", "T cells CD8",
  "T cells CD4 naive", "T cells CD4 memory resting",
  "T cells CD4 memory activated", "T cells follicular helper",
  "T cells regulatory (Tregs)", "T cells gamma delta",
  "NK cells resting", "NK cells activated", "Monocytes",
  "Macrophages M0", "Macrophages M1", "Macrophages M2",
  "Dendritic cells resting", "Dendritic cells activated",
  "Mast cells resting", "Mast cells activated", "Eosinophils",
  "Neutrophils"
)
stopifnot(length(immune_cols) == 22L, !anyDuplicated(immune_cols))

ZERO_THRESHOLD <- 50        # exclude a cell type when >50% values are zero
FDR_THRESHOLD <- 0.05
LOG2FC_PSEUDOCOUNT <- 1e-4  # matches documented historical comparisons

# --- Provenance: archive integrity and historical settings ----------------
original_files <- c(pd_result_file, oscc_result_file,
                    pd_input_file, oscc_input_file,
                    pd_log_file, oscc_log_file)
manifest <- data.frame(
  File = basename(original_files),
  MD5 = unname(tools::md5sum(original_files)),
  Size_bytes = unname(file.info(original_files)$size),
  stringsAsFactors = FALSE
)
utils::write.csv(manifest, file.path(output_dir, "03_input_manifest.csv"),
                 row.names = FALSE)

check_job_log <- function(path, label, input_name) {
  lines <- readLines(path, warn = FALSE)
  required <- c(
    "LM22.update-gene-symbols.txt",
    "QN: TRUE",
    "Number of permutations set to: 1",
    input_name
  )
  ok <- vapply(required, function(item) {
    any(grepl(item, lines, fixed = TRUE))
  }, logical(1))
  if (!all(ok)) {
    stop(label, " log settings differ from documented settings: ",
         paste(required[!ok], collapse = " | "))
  }
  cat(label, ": historical LM22 / QN=TRUE / permutations=1 confirmed\n")
  invisible(TRUE)
}
check_job_log(pd_log_file, "Job29 PD", "CIBERSORT_input_PD_linear.txt")
check_job_log(oscc_log_file, "Job30 OSCC", "CIBERSORT_input_OSCC_linear.txt")

# Verify the original mixture file layout for archival traceability.
# Reading this table is for QC only. It is never submitted to CIBERSORTx.
audit_mixture <- function(path, expected_gsm, expected_rows, label) {
  tab <- utils::read.delim(path, check.names = FALSE,
                           stringsAsFactors = FALSE)
  if (ncol(tab) != length(expected_gsm) + 1L) {
    stop(label, ": unexpected mixture column count")
  }
  if (!identical(names(tab)[-1], as.character(expected_gsm))) {
    stop(label, ": original mixture GSM columns do not match metadata/order")
  }
  genes <- as.character(tab[[1L]])
  output <- data.frame(
    Dataset = label,
    Mixture_rows = nrow(tab),
    Unique_gene_symbols = length(unique(genes)),
    Duplicate_gene_symbol_rows = sum(duplicated(genes)),
    Samples = ncol(tab) - 1L,
    stringsAsFactors = FALSE
  )
  if (output$Mixture_rows != expected_rows ||
      output$Unique_gene_symbols != 21755L) {
    stop(label, ": original mixture row/gene counts changed; audit inputs")
  }
  rm(tab)
  output
}

# --- Step 01 metadata, and archived CIBERSORTx fractions -------------------
pre <- readRDS(metadata_file)
if (!all(c("sample_map_pd", "sample_map_oscc") %in% names(pre))) {
  stop("Script 01 RDS is missing sample maps")
}
map_pd <- pre$sample_map_pd
map_oscc <- pre$sample_map_oscc

validate_map <- function(map, label, disease, n_disease, n_control) {
  required <- c("GSM", "Group")
  if (!all(required %in% names(map)) ||
      anyNA(map$GSM) || anyDuplicated(map$GSM) ||
      anyNA(map$Group)) {
    stop(label, ": invalid GSM or Group metadata")
  }
  groups <- as.character(map$Group)
  if (!all(groups %in% c(disease, "Control")) ||
      sum(groups == disease) != n_disease ||
      sum(groups == "Control") != n_control) {
    stop(label, ": group counts differ from original design")
  }
  invisible(TRUE)
}
validate_map(map_pd, "PD", "PD", 183L, 64L)
validate_map(map_oscc, "OSCC", "OSCC", 167L, 45L)

mixture_audit <- rbind(
  audit_mixture(pd_input_file, map_pd$GSM, 28655L, "GSE10334"),
  audit_mixture(oscc_input_file, map_oscc$GSM, 28704L, "GSE30784")
)
utils::write.csv(mixture_audit,
                 file.path(output_dir, "03_original_mixture_QC.csv"),
                 row.names = FALSE)
print(mixture_audit, row.names = FALSE)

read_and_align_fractions <- function(path, map, label) {
  df <- utils::read.csv(path, check.names = FALSE,
                        stringsAsFactors = FALSE)
  expected_cols <- c("Mixture", immune_cols,
                     "P-value", "Correlation", "RMSE")
  if (!setequal(names(df), expected_cols) ||
      ncol(df) != length(expected_cols)) {
    stop(label, ": expected 22 LM22 fractions + Mixture/P-value/Correlation/RMSE")
  }
  ids <- as.character(df$Mixture)
  gsm <- as.character(map$GSM)
  if (anyNA(ids) || anyDuplicated(ids) || !setequal(ids, gsm) ||
      nrow(df) != nrow(map)) {
    stop(label, ": missing/duplicate/mismatched GSM identifiers")
  }
  same_order <- identical(ids, gsm)
  df <- df[match(gsm, ids), , drop = FALSE]
  rownames(df) <- NULL
  if (!identical(as.character(df$Mixture), gsm)) {
    stop(label, ": sample alignment failed")
  }
  for (cell in immune_cols) {
    v <- df[[cell]]
    if (!is.numeric(v) || anyNA(v) || any(!is.finite(v)) ||
        any(v < 0 | v > 1)) {
      stop(label, ": invalid fraction values in ", cell)
    }
  }
  fraction_sums <- rowSums(df[immune_cols])
  if (any(abs(fraction_sums - 1) > 0.01)) {
    warning(label, ": LM22 fractions do not sum approximately to 1",
            "; inspect the original file")
  }
  list(data = df, same_order = same_order)
}
aligned_pd <- read_and_align_fractions(pd_result_file, map_pd, "Job29 PD")
aligned_oscc <- read_and_align_fractions(oscc_result_file, map_oscc,
                                        "Job30 OSCC")
pd <- aligned_pd$data
oscc <- aligned_oscc$data

alignment_qc <- data.frame(
  Dataset = c("GSE10334", "GSE30784"),
  Samples = c(nrow(pd), nrow(oscc)),
  Same_GSM_set = c(TRUE, TRUE),
  Same_GSM_order_in_original = c(aligned_pd$same_order,
                                  aligned_oscc$same_order),
  Duplicated_GSM = c(anyDuplicated(pd$Mixture),
                     anyDuplicated(oscc$Mixture)),
  stringsAsFactors = FALSE
)
utils::write.csv(alignment_qc, file.path(output_dir, "03_GSM_alignment_QC.csv"),
                 row.names = FALSE)

# --- Zero-value QC (full cohort, prior to statistical comparisons) --------
zero_pd <- vapply(pd[immune_cols], function(v) mean(v == 0) * 100, numeric(1))
zero_oscc <- vapply(oscc[immune_cols],
                    function(v) mean(v == 0) * 100, numeric(1))
reliable_pd <- names(zero_pd)[zero_pd <= ZERO_THRESHOLD]
reliable_oscc <- names(zero_oscc)[zero_oscc <= ZERO_THRESHOLD]
qc_cells <- data.frame(
  Cell = immune_cols,
  PD_zeros = round(unname(zero_pd), 2),
  OSCC_zeros = round(unname(zero_oscc), 2),
  PD_pass = unname(zero_pd <= ZERO_THRESHOLD),
  OSCC_pass = unname(zero_oscc <= ZERO_THRESHOLD),
  stringsAsFactors = FALSE
)
utils::write.csv(qc_cells, file.path(output_dir, "03_cell_zero_QC.csv"),
                 row.names = FALSE)
cat("LM22 QC: PD ", length(reliable_pd), "/22; OSCC ",
    length(reliable_oscc), "/22\n", sep = "")

# --- Historical unpaired Wilcoxon + BH comparisons -------------------------
# Recreates the validated historical table, NOT patient-independent evidence.
analyze_immune <- function(df, map, disease, cells) {
  groups <- as.character(map$Group)
  if (!identical(as.character(df$Mixture), as.character(map$GSM))) {
    stop("GSM alignment was lost")
  }
  out <- do.call(rbind, lapply(cells, function(cell) {
    disease_values <- df[[cell]][groups == disease]
    control_values <- df[[cell]][groups == "Control"]
    p <- stats::wilcox.test(disease_values, control_values,
                            exact = FALSE)$p.value
    data.frame(
      Cell = cell,
      Mean_disease = mean(disease_values),
      Mean_control = mean(control_values),
      log2FC = log2((mean(disease_values) + LOG2FC_PSEUDOCOUNT) /
                   (mean(control_values) + LOG2FC_PSEUDOCOUNT)),
      P_value = p,
      stringsAsFactors = FALSE
    )
  }))
  out$FDR <- stats::p.adjust(out$P_value, method = "BH")
  out[order(out$FDR), , drop = FALSE]
}
diff_pd <- analyze_immune(pd, map_pd, "PD", reliable_pd)
diff_oscc <- analyze_immune(oscc, map_oscc, "OSCC", reliable_oscc)
utils::write.csv(diff_pd,
                 file.path(output_dir, "03_PD_immune_comparison.csv"),
                 row.names = FALSE)
utils::write.csv(diff_oscc,
                 file.path(output_dir, "03_OSCC_immune_comparison.csv"),
                 row.names = FALSE)

# --- PD paired within-patient sensitivity (64 matched patients) ------------
if (!"Patient_ID" %in% names(map_pd) ||
    anyNA(map_pd$Patient_ID) ||
    length(unique(map_pd$Patient_ID)) != 90L) {
  stop("PD sample metadata must have 90 valid Patient_ID values")
}
patient_means <- stats::aggregate(
  pd[reliable_pd],
  by = list(Patient_ID = as.character(map_pd$Patient_ID),
            Group = as.character(map_pd$Group)),
  FUN = mean
)
pd_pat <- patient_means[patient_means$Group == "PD", , drop = FALSE]
control_pat <- patient_means[patient_means$Group == "Control", , drop = FALSE]
paired_ids <- intersect(pd_pat$Patient_ID, control_pat$Patient_ID)
if (length(paired_ids) != 64L) {
  stop("Expected 64 patients with both PD and Control; found ",
       length(paired_ids))
}
pd_pat <- pd_pat[match(paired_ids, pd_pat$Patient_ID), , drop = FALSE]
control_pat <- control_pat[match(paired_ids, control_pat$Patient_ID),
                           , drop = FALSE]
if (!identical(pd_pat$Patient_ID, control_pat$Patient_ID)) {
  stop("Patient pairing did not align")
}

paired_results <- do.call(rbind, lapply(reliable_pd, function(cell) {
  disease_values <- pd_pat[[cell]]
  control_values <- control_pat[[cell]]
  delta <- disease_values - control_values
  p <- stats::wilcox.test(disease_values, control_values,
                          paired = TRUE, exact = FALSE)$p.value
  data.frame(
    Cell = cell,
    Median_difference = stats::median(delta),
    PD_higher = sum(delta > 0),
    PD_lower = sum(delta < 0),
    P_value = p,
    stringsAsFactors = FALSE
  )
}))
paired_results$FDR <- stats::p.adjust(paired_results$P_value, method = "BH")
paired_results <- paired_results[order(paired_results$FDR), , drop = FALSE]
utils::write.csv(
  paired_results,
  file.path(output_dir, "03_PD_paired_patient_sensitivity.csv"),
  row.names = FALSE
)

# Compare ALL QC-passed PD populations, not just the significant subset.
pd_comparison_full <- data.frame(
  Cell = reliable_pd,
  Original_FDR = diff_pd$FDR[match(reliable_pd, diff_pd$Cell)],
  Paired_FDR = paired_results$FDR[
    match(reliable_pd, paired_results$Cell)
  ],
  stringsAsFactors = FALSE
)
pd_comparison_full$Original_significant <-
  pd_comparison_full$Original_FDR < FDR_THRESHOLD
pd_comparison_full$Paired_significant <-
  pd_comparison_full$Paired_FDR < FDR_THRESHOLD

utils::write.csv(
  pd_comparison_full,
  file.path(output_dir, "03_PD_all_population_sensitivity.csv"),
  row.names = FALSE
)
pd_comparison_sig <- pd_comparison_full[
  pd_comparison_full$Original_significant |
    pd_comparison_full$Paired_significant,
  c("Cell", "Original_FDR", "Paired_FDR"), drop = FALSE
]
pd_comparison_sig <- pd_comparison_sig[
  order(pd_comparison_sig$Cell), , drop = FALSE
]
utils::write.csv(
  pd_comparison_sig,
  file.path(output_dir, "03_PD_sensitivity_comparison.csv"),
  row.names = FALSE
)

# --- Figures: newly generated descriptive figures, NOT claimed identical ---
# to figures in the submitted manuscript without visual comparison.
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  warning("Package ggplot2 not installed: CSV files saved, figures skipped.")
} else {
  draw_barplot <- function(df, groups, cells, disease, outfile) {
    fractions <- as.matrix(df[cells])
    long <- data.frame(
      Group = rep(as.character(groups), times = length(cells)),
      Cell_type = rep(cells, each = nrow(df)),
      Fraction = as.vector(fractions),
      stringsAsFactors = FALSE
    )
    means <- stats::aggregate(Fraction ~ Group + Cell_type, long, FUN = mean)
    means$Group <- factor(means$Group, levels = c("Control", disease))
    fig <- ggplot2::ggplot(
      means, ggplot2::aes(x = Cell_type, y = Fraction, fill = Group)
    ) +
      ggplot2::geom_col(position = "dodge", width = 0.72) +
      ggplot2::labs(
        title = paste0("Estimated immune cell fractions: ", disease,
                       " vs Control"),
        x = NULL, y = "Mean estimated fraction"
      ) +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 55, hjust = 1),
        legend.position = "top"
      )
    ggplot2::ggsave(outfile, fig, width = 13, height = 7, dpi = 220,
                    limitsize = FALSE)
  }
  draw_m1_m2_boxplot <- function(df, groups, disease, outfile) {
    cells <- c("Macrophages M1", "Macrophages M2")
    long <- data.frame(
      Group = rep(as.character(groups), times = length(cells)),
      Cell = rep(cells, each = nrow(df)),
      Fraction = as.vector(as.matrix(df[cells])),
      stringsAsFactors = FALSE
    )
    long$Group <- factor(long$Group, levels = c("Control", disease))
    fig <- ggplot2::ggplot(
      long, ggplot2::aes(x = Group, y = Fraction, fill = Group)
    ) +
      ggplot2::geom_boxplot(width = 0.58, outlier.alpha = 0.35) +
      ggplot2::facet_wrap(~Cell, scales = "free_y") +
      ggplot2::labs(
        title = paste0("Estimated macrophage fractions: ", disease,
                       " vs Control"),
        x = NULL, y = "Estimated relative fraction"
      ) +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::theme(legend.position = "none")
    ggplot2::ggsave(outfile, fig, width = 9, height = 5, dpi = 220)
  }
  draw_barplot(pd, map_pd$Group, reliable_pd, "PD",
               file.path(output_dir, "03_PD_cell_mean_barplot.png"))
  draw_barplot(oscc, map_oscc$Group, reliable_oscc, "OSCC",
               file.path(output_dir, "03_OSCC_cell_mean_barplot.png"))
  draw_m1_m2_boxplot(pd, map_pd$Group, "PD",
                     file.path(output_dir, "03_PD_M1_M2_boxplots.png"))
  draw_m1_m2_boxplot(oscc, map_oscc$Group, "OSCC",
                     file.path(output_dir, "03_OSCC_M1_M2_boxplots.png"))
}

# --- Results check (descriptive, never changes thresholds to force match) ---
counts <- data.frame(
  Metric = c("PD_QC_pass", "OSCC_QC_pass",
             "PD_unpaired_significant", "OSCC_unpaired_significant",
             "PD_paired_patients", "PD_paired_significant"),
  Historical_observed = c(16L, 14L, 9L, 8L, 64L, 9L),
  Current = c(length(reliable_pd), length(reliable_oscc),
              sum(diff_pd$FDR < FDR_THRESHOLD),
              sum(diff_oscc$FDR < FDR_THRESHOLD),
              length(paired_ids),
              sum(paired_results$FDR < FDR_THRESHOLD))
)
counts$Match <- counts$Historical_observed == counts$Current
utils::write.csv(counts,
                 file.path(output_dir, "03_reproduction_count_QC.csv"),
                 row.names = FALSE)
print(counts, row.names = FALSE)
if (!all(counts$Match)) {
  warning("Some Script 03 counts differ from the previously observed run. ",
          "Review input metadata, code and package versions; ",
          "do NOT tune thresholds to force agreement.")
}

cat("\nPD M1/M2 (unpaired):\n")
print(diff_pd[diff_pd$Cell %in% c("Macrophages M1", "Macrophages M2"),
              , drop = FALSE], row.names = FALSE)
cat("\nPD M1/M2 (paired):\n")
print(paired_results[paired_results$Cell %in%
                       c("Macrophages M1", "Macrophages M2"),
                     , drop = FALSE], row.names = FALSE)
cat("\nOSCC M0/M1/M2 (unpaired):\n")
print(diff_oscc[diff_oscc$Cell %in%
                  c("Macrophages M0", "Macrophages M1", "Macrophages M2"),
                , drop = FALSE], row.names = FALSE)

if (any(diff_pd$P_value == 0) || any(diff_oscc$P_value == 0)) {
  cat("NOTE: Any reported Wilcoxon P_value=0 is a numerical underflow,",
      " NOT a literal probability of zero.\n")
}
utils::write.csv(
  data.frame(
    Name = c("HistoricalJob_PD", "HistoricalJob_OSCC",
             "Signature", "Quantile_normalization",
             "Permutations", "Zero_filter_percent",
             "Wilcoxon_unpaired", "FDR_method",
             "Paired_analysis", "Original_mix_PD_rows",
             "Original_mix_OSCC_rows",
             "Candidate_Script01_gene_rows",
             "Deconvolution_recomputed"),
    Value = c("Job29", "Job30", "LM22.update-gene-symbols.txt",
              "TRUE (both jobs)", "1 (both jobs)",
              "Exclude cell types with >50% zeros in each cohort",
              "wilcox.test(exact=FALSE) on sample fractions",
              "Benjamini-Hochberg within QC-passed cell types",
              "64 paired patients: mean per patient/group, paired Wilcoxon",
              "28655", "28704", "20848",
              "NO: archived CIBERSORTx results analyzed"),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "03_methods_and_provenance.csv"),
  row.names = FALSE
)
writeLines(capture.output(utils::sessionInfo()),
           file.path(output_dir, "03_sessionInfo.txt"))

cat("\nSTEP 03 COMPLETE. Results: ",
    normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n", sep = "")
cat("Original CIBERSORTx files were read but NOT changed.\n")
