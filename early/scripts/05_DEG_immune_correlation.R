# =============================================================================
# 05_DEG_immune_correlation.R
# PD–OSCC transcriptomic study | manuscript-aligned reproducibility CANDIDATE
#
# From repository root:
#   source("scripts/05_DEG_immune_correlation.R")
#
# This script re-analyses the original DEG lists (Script 02) and archived
# LM22 CIBERSORTx Job29/Job30 fractions (Script 03). It does NOT rerun CIBERSORTx.
#
# Manuscript claims to audit:
#  * Spearman's rho; |rho| > 0.3 and BH-adjusted p < 0.05
#  * PD: 57 genes x 16 LM22 cell types = 912 tested pairs; 346 significant
#  * OSCC: top 50 of 1266 DEGs by |log2FC| x 14 types = 700; 180 significant
#
# Historical scripts selected one Affymetrix probe per gene using the highest
# mean expression within each dataset, after taking FIRST symbol before "///".
# The 20,848-row representative-gene matrices of the revised Script 01 choose
# probes by POOLED PD+OSCC expression instead and are NOT interchangeable here.
# The original code adjusted p-values only after filtering |rho| > 0.3; this
# candidate also calculates BH across ALL 912/700 tests as described in Methods.
# The difference is audited explicitly, not silently hidden.
#
# Note: PD unpaired biopsy-level Spearman tests do NOT account for 90 patient
# clusters; this step reproduces the historical Figure 5 analysis, while
# patient-cluster sensitivity is reported separately in other analyses.
# =============================================================================

required_packages <- c("Biobase", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

COR_THRESHOLD <- 0.30
FDR_THRESHOLD <- 0.05
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "results", "05_DEG_immune_correlation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- c(
  preprocessed = file.path(root, "results", "01_preprocessing", "01_preprocessed_esets.rds"),
  PD_DEG = file.path(root, "results", "02_DEG", "02_PD_unique_DEG_entries.csv"),
  OSCC_DEG = file.path(root, "results", "02_DEG", "02_OSCC_unique_DEG_entries.csv"),
  PD_CIBERSORTx = file.path(root, "data", "cibersortx_historical", "CIBERSORTx_Job29_Results.csv"),
  OSCC_CIBERSORTx = file.path(root, "data", "cibersortx_historical", "CIBERSORTx_Job30_Results.csv"),
  QC_03 = file.path(root, "results", "03_CIBERSORT", "03_cell_zero_QC.csv")
)
if (!all(file.exists(input_files))) {
  stop("Missing required input(s):\n",
       paste(input_files[!file.exists(input_files)], collapse = "\n"))
}

# Optional copies of the historical Correlation_Table.xlsx sheets. They are
# only audit references, NEVER used to calculate the new correlations.
references <- c(
  PD = file.path(root, "data", "reference", "05_reference_PD_pairs.csv"),
  OSCC = file.path(root, "data", "reference", "05_reference_OSCC_pairs.csv")
)
manifest_paths <- c(input_files, references[file.exists(references)])
input_manifest <- data.frame(
  Role = names(manifest_paths),
  File = vapply(manifest_paths, function(x) substring(x, nchar(root) + 2L), character(1)),
  MD5 = unname(tools::md5sum(manifest_paths)),
  Size_bytes = unname(file.info(manifest_paths)$size),
  row.names = NULL
)
write.csv(input_manifest, file.path(out_dir, "05_input_manifest.csv"), row.names = FALSE)

pre <- readRDS(input_files[["preprocessed"]])
if (!all(c("gset_pd", "gset_oscc", "sample_map_pd", "sample_map_oscc") %in% names(pre))) {
  stop("Script 01 RDS has unexpected contents")
}

read_deg <- function(path, dataset, expected) {
  tab <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("Gene", "logFC") %in% names(tab))) {
    stop(dataset, ": Script 02 DEG table lacks Gene/logFC")
  }
  tab$Gene <- trimws(as.character(tab$Gene))
  if (nrow(tab) != expected || anyNA(tab$Gene) || anyDuplicated(tab$Gene) ||
      any(!nzchar(tab$Gene)) || anyNA(tab$logFC) || any(!is.finite(tab$logFC))) {
    stop(dataset, ": DEG count/identifiers/logFC inconsistent with audited Script 02")
  }
  tab
}
pd_deg <- read_deg(input_files[["PD_DEG"]], "PD", 57L)
oscc_deg <- read_deg(input_files[["OSCC_DEG"]], "OSCC", 1266L)
# Original Script 05: dplyr::arrange(desc(abs(logFC))) %>% head(50).
# Use original row order to break any |logFC| ties deterministically.
oscc_top <- oscc_deg[order(-abs(oscc_deg$logFC), seq_len(nrow(oscc_deg))), , drop = FALSE]
oscc_top <- utils::head(oscc_top, 50L)
write.csv(oscc_top[, c("Gene", "logFC")],
          file.path(out_dir, "05_OSCC_top50_gene_selection.csv"), row.names = FALSE)

immune_cols <- c(
  "B cells naive", "B cells memory", "Plasma cells", "T cells CD8",
  "T cells CD4 naive", "T cells CD4 memory resting",
  "T cells CD4 memory activated", "T cells follicular helper",
  "T cells regulatory (Tregs)", "T cells gamma delta",
  "NK cells resting", "NK cells activated", "Monocytes",
  "Macrophages M0", "Macrophages M1", "Macrophages M2",
  "Dendritic cells resting", "Dendritic cells activated",
  "Mast cells resting", "Mast cells activated", "Eosinophils", "Neutrophils"
)
stopifnot(length(immune_cols) == 22L, !anyDuplicated(immune_cols))

load_immune <- function(path, map, dataset) {
  df <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  needed <- c("Mixture", immune_cols, "P-value", "Correlation", "RMSE")
  if (!setequal(names(df), needed) || ncol(df) != length(needed)) {
    stop(dataset, ": CIBERSORTx result columns do not match LM22 layout")
  }
  id <- as.character(df$Mixture)
  gsm <- as.character(map$GSM)
  if (anyNA(id) || anyDuplicated(id) || !setequal(id, gsm) || nrow(df) != length(gsm)) {
    stop(dataset, ": CIBERSORTx GSM set/uniqueness does not match Script 01")
  }
  df <- df[match(gsm, id), , drop = FALSE]
  if (!identical(as.character(df$Mixture), gsm)) stop(dataset, ": GSM alignment failed")
  for (cell in immune_cols) {
    v <- df[[cell]]
    if (!is.numeric(v) || anyNA(v) || any(!is.finite(v)) || any(v < 0) || any(v > 1)) {
      stop(dataset, ": invalid CIBERSORTx LM22 fractions in ", cell)
    }
  }
  zeros <- vapply(df[immune_cols], function(v) 100 * mean(v == 0), numeric(1))
  list(matrix = as.matrix(df[immune_cols]), zeros = zeros, gsm = gsm)
}

immune_pd_info <- load_immune(input_files[["PD_CIBERSORTx"]], pre$sample_map_pd, "PD")
immune_oscc_info <- load_immune(input_files[["OSCC_CIBERSORTx"]], pre$sample_map_oscc, "OSCC")

qc03 <- utils::read.csv(input_files[["QC_03"]], stringsAsFactors = FALSE)
if (!all(c("Cell", "PD_pass", "OSCC_pass") %in% names(qc03)) ||
    !setequal(qc03$Cell, immune_cols) || anyDuplicated(qc03$Cell)) {
  stop("Step 03 QC file has unexpected cell names/columns")
}
actual_pd_pass <- immune_pd_info$zeros <= 50
actual_oscc_pass <- immune_oscc_info$zeros <= 50
if (!identical(unname(actual_pd_pass[qc03$Cell]), as.logical(qc03$PD_pass)) ||
    !identical(unname(actual_oscc_pass[qc03$Cell]), as.logical(qc03$OSCC_pass))) {
  stop("Step 03 saved QC differs from independent recomputation of zero fractions")
}
reliable_pd <- immune_cols[actual_pd_pass[immune_cols]]
reliable_oscc <- immune_cols[actual_oscc_pass[immune_cols]]
if (length(reliable_pd) != 16L || length(reliable_oscc) != 14L) {
  stop("Unexpected QC cell counts; check historical Job29/Job30 data")
}

# Historical Script 01 renamed every probe by the *FIRST* symbol before ///,
# and historical Script 05 retained the highest-mean probe per gene WITHIN
# each disease. Reconstruct directly from Script 01 archived ExpressionSets,
# never silently use the newer 20,848-row pooled-representative matrices.
first_symbol <- function(annotation) {
  vapply(strsplit(as.character(annotation), "///", fixed = TRUE), function(x) {
    if (!length(x) || is.na(x[[1L]])) NA_character_ else trimws(x[[1L]])
  }, character(1))
}

extract_historical_expr <- function(eset, map, wanted_genes, label) {
  expr <- Biobase::exprs(eset)
  annotation <- Biobase::fData(eset)
  if (!"Gene.symbol" %in% names(annotation)) {
    stop(label, ": ExpressionSet lacks Gene.symbol annotation")
  }
  if (!identical(colnames(expr), as.character(map$GSM))) {
    stop(label, ": ExpressionSet sample order differs from audited GSM metadata")
  }
  symbol <- first_symbol(annotation[["Gene.symbol"]])
  if (length(symbol) != nrow(expr)) stop(label, ": invalid probe annotation length")
  candidate_rows <- which(!is.na(symbol) & symbol %in% wanted_genes)
  if (!length(candidate_rows)) stop(label, ": no candidate probes found")
  means <- rowMeans(expr[candidate_rows, , drop = FALSE])
  selected_idx <- vapply(wanted_genes, function(g) {
    ix <- which(symbol[candidate_rows] == g)
    if (!length(ix)) return(NA_integer_)
    # which.max returns first index on ties, matching slice_max(..., with_ties=FALSE)
    candidate_rows[ix[[which.max(means[ix])]]]
  }, integer(1))
  if (anyNA(selected_idx)) {
    stop(label, ": missing representative probe(s) for: ",
         paste(wanted_genes[is.na(selected_idx)], collapse = ", "))
  }
  selected_expr <- expr[selected_idx, , drop = FALSE]
  rownames(selected_expr) <- wanted_genes
  if (anyNA(selected_expr) || any(!is.finite(selected_expr))) {
    stop(label, ": invalid representative-probe expression values")
  }
  probe_map <- data.frame(
    Dataset = label,
    Gene = wanted_genes,
    GPL570_probe_ID = rownames(expr)[selected_idx],
    Candidate_probes_for_gene = vapply(wanted_genes,
                                     function(g) sum(symbol[candidate_rows] == g), integer(1)),
    Mean_log2_expression = rowMeans(selected_expr),
    stringsAsFactors = FALSE
  )
  list(expression = selected_expr, probe_map = probe_map)
}

expr_pd_info <- extract_historical_expr(pre$gset_pd, pre$sample_map_pd,
                                        pd_deg$Gene, "PD")
expr_oscc_info <- extract_historical_expr(pre$gset_oscc, pre$sample_map_oscc,
                                          oscc_top$Gene, "OSCC")
write.csv(rbind(expr_pd_info$probe_map, expr_oscc_info$probe_map),
          file.path(out_dir, "05_historical_probe_selection.csv"), row.names = FALSE)

# Both expression matrices and CIBERSORTx fractions are in identical GSM order.
if (!identical(colnames(expr_pd_info$expression), immune_pd_info$gsm) ||
    !identical(colnames(expr_oscc_info$expression), immune_oscc_info$gsm)) {
  stop("Expression and immune sample GSM orders disagree")
}
expr_pd <- expr_pd_info$expression
expr_oscc <- expr_oscc_info$expression
immune_pd <- immune_pd_info$matrix[, reliable_pd, drop = FALSE]
immune_oscc <- immune_oscc_info$matrix[, reliable_oscc, drop = FALSE]
rownames(immune_pd) <- immune_pd_info$gsm
rownames(immune_oscc) <- immune_oscc_info$gsm

# All gene-cell tests, with original cor.test() Spearman call semantics.
# Defaults (including the approximate/exact calculation R chooses for ties)
# are intentionally unchanged from the historical Script 05.
calc_all_pairs <- function(expr, immune, label) {
  stopifnot(ncol(expr) == nrow(immune),
            identical(colnames(expr), rownames(immune)))
  if (any(!is.finite(expr)) || any(!is.finite(immune))) {
    stop(label, ": NA/Inf prevents correlation analysis")
  }
  n_tests <- nrow(expr) * ncol(immune)
  rows <- vector("list", n_tests)
  k <- 0L
  for (i in seq_len(nrow(expr))) {
    x <- as.numeric(expr[i, ])
    if (length(unique(x)) < 2L) stop(label, ": constant gene expression: ", rownames(expr)[i])
    for (j in seq_len(ncol(immune))) {
      y <- as.numeric(immune[, j])
      if (length(unique(y)) < 2L) stop(label, ": constant immune fraction: ", colnames(immune)[j])
      test <- suppressWarnings(stats::cor.test(x, y, method = "spearman"))
      k <- k + 1L
      rows[[k]] <- data.frame(
        Gene = rownames(expr)[i], Cell_type = colnames(immune)[j],
        Correlation = unname(as.numeric(test$estimate)),
        p_value = unname(as.numeric(test$p.value)),
        stringsAsFactors = FALSE
      )
    }
  }
  ans <- do.call(rbind, rows)
  if (nrow(ans) != n_tests || anyNA(ans$Correlation) ||
      anyNA(ans$p_value) || any(!is.finite(ans$p_value))) {
    stop(label, ": invalid/missing Spearman test output")
  }
  # Correct manuscript-defined family-wise BH for each disease, before filter.
  ans$p_adj <- stats::p.adjust(ans$p_value, method = "BH")
  ans$Significant <- abs(ans$Correlation) > COR_THRESHOLD &
    ans$p_adj < FDR_THRESHOLD
  ans$Pass_abs_r <- abs(ans$Correlation) > COR_THRESHOLD
  ans
}

cat("Calculating PD Spearman tests (", nrow(expr_pd), " x ", ncol(immune_pd), ")...\n", sep = "")
all_pd <- calc_all_pairs(expr_pd, immune_pd, "PD")
cat("Calculating OSCC Spearman tests (", nrow(expr_oscc), " x ", ncol(immune_oscc), ")...\n", sep = "")
all_oscc <- calc_all_pairs(expr_oscc, immune_oscc, "OSCC")

# Compare the corrected BH family to the historical filter-before-BH method.
# Historical Correlation_Table.xlsx used BH only after filtering |rho| > 0.3.
historical_filter_BH <- function(all_table) {
  tab <- all_table[all_table$Pass_abs_r, , drop = FALSE]
  tab$p_adj_historical_filtered <- stats::p.adjust(tab$p_value, method = "BH")
  tab[order(tab$p_adj, -abs(tab$Correlation), tab$Gene, tab$Cell_type),,
      drop = FALSE]
}
legacy_pd <- historical_filter_BH(all_pd)
legacy_oscc <- historical_filter_BH(all_oscc)
cor_pd_sig <- legacy_pd[legacy_pd$Significant, , drop = FALSE]
cor_oscc_sig <- legacy_oscc[legacy_oscc$Significant, , drop = FALSE]

for (spec in list(
  list(name = "05_PD_all_912_pairs.csv", data = all_pd),
  list(name = "05_OSCC_all_700_pairs.csv", data = all_oscc),
  list(name = "05_PD_legacy_threshold_pairs.csv", data = legacy_pd),
  list(name = "05_OSCC_legacy_threshold_pairs.csv", data = legacy_oscc),
  list(name = "05_PD_significant_pairs.csv", data = cor_pd_sig),
  list(name = "05_OSCC_significant_pairs.csv", data = cor_oscc_sig)
)) {
  utils::write.csv(spec$data, file.path(out_dir, spec$name), row.names = FALSE)
}

if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(
    list(PD_correlations = cor_pd_sig,
         OSCC_correlations = cor_oscc_sig),
    file.path(out_dir, "05_Correlation_Table_corrected_BH.xlsx")
  )
} else {
  message("writexl not installed: audit CSV tables are available; xlsx export skipped")
}

# Historical workbook-based comparison: per-pair keys & rho, NOT historical
# adjusted p-values (which used a different multiple-testing family).
compare_reference <- function(current_sig, reference_file, label) {
  if (!file.exists(reference_file)) {
    message(label, ": historical reference CSV not installed; coefficient audit skipped")
    return(data.frame(Dataset = label, Reference_available = FALSE,
                      Historical_pairs = NA_integer_, Matched_pairs = NA_integer_,
                      Missing_original_pairs = NA_integer_, Additional_current_pairs = NA_integer_,
                      Max_abs_rho_delta = NA_real_, Coefficients_match_1e_6 = NA))
  }
  ref <- utils::read.csv(reference_file, stringsAsFactors = FALSE)
  if (!all(c("Gene", "Cell_type", "Correlation") %in% names(ref))) {
    stop(label, ": reference CSV lacks required columns")
  }
  key <- function(tab) paste(tab$Gene, tab$Cell_type, sep = "\r")
  ref_keys <- key(ref)
  cur_keys <- key(current_sig)
  if (anyDuplicated(ref_keys) || anyDuplicated(cur_keys)) {
    stop(label, ": duplicate gene-cell keys in correlation table")
  }
  both <- intersect(ref_keys, cur_keys)
  ref_match <- match(both, ref_keys)
  cur_match <- match(both, cur_keys)
  differences <- abs(ref$Correlation[ref_match] -
                     current_sig$Correlation[cur_match])
  detail <- data.frame(
    Gene = ref$Gene[ref_match],
    Cell_type = ref$Cell_type[ref_match],
    Historical_rho = ref$Correlation[ref_match],
    Current_rho = current_sig$Correlation[cur_match],
    Absolute_rho_delta = differences,
    stringsAsFactors = FALSE
  )
  utils::write.csv(detail,
    file.path(out_dir, paste0("05_", label, "_reference_rho_comparison.csv")),
    row.names = FALSE)
  utils::write.csv(ref[!ref_keys %in% cur_keys, c("Gene", "Cell_type"), drop = FALSE],
    file.path(out_dir, paste0("05_", label, "_missing_historical_pairs.csv")),
    row.names = FALSE)
  utils::write.csv(current_sig[!cur_keys %in% ref_keys, c("Gene", "Cell_type"), drop = FALSE],
    file.path(out_dir, paste0("05_", label, "_new_current_pairs.csv")),
    row.names = FALSE)
  data.frame(
    Dataset = label, Reference_available = TRUE,
    Historical_pairs = nrow(ref), Matched_pairs = length(both),
    Missing_original_pairs = sum(!ref_keys %in% cur_keys),
    Additional_current_pairs = sum(!cur_keys %in% ref_keys),
    Max_abs_rho_delta = if (length(differences)) max(differences) else NA_real_,
    Coefficients_match_1e_6 = length(both) == nrow(ref) &&
      length(both) == nrow(current_sig) &&
      all(differences <= 1e-6)
  )
}
reference_qc <- rbind(
  compare_reference(cor_pd_sig, references[["PD"]], "PD"),
  compare_reference(cor_oscc_sig, references[["OSCC"]], "OSCC")
)
utils::write.csv(reference_qc,
                 file.path(out_dir, "05_historical_correlation_reference_QC.csv"),
                 row.names = FALSE)

summary <- data.frame(
  Dataset = c("PD", "OSCC"),
  Genes_tested = c(nrow(expr_pd), nrow(expr_oscc)),
  Immune_cells_tested = c(ncol(immune_pd), ncol(immune_oscc)),
  Pairs_tested = c(nrow(all_pd), nrow(all_oscc)),
  Pairs_above_abs_r_0_3 = c(nrow(legacy_pd), nrow(legacy_oscc)),
  Significant_BH_all_tests = c(nrow(cor_pd_sig), nrow(cor_oscc_sig)),
  Historical_manuscript_significant = c(346L, 180L),
  stringsAsFactors = FALSE
)
summary$Match_manuscript_count <- summary$Significant_BH_all_tests ==
  summary$Historical_manuscript_significant
utils::write.csv(summary, file.path(out_dir, "05_correlation_summary.csv"),
                 row.names = FALSE)

# Diagnostic-only heatmaps. These are NOT claimed to reproduce Figure 5
# without a direct visual comparison to the submitted manuscript figure.
make_heatmap <- function(all_table, label, file_name) {
  tab <- all_table
  tab$Gene <- factor(tab$Gene, levels = rev(unique(tab$Gene)))
  tab$Cell_type <- factor(tab$Cell_type, levels = unique(tab$Cell_type))
  fig <- ggplot2::ggplot(tab, ggplot2::aes(x = Cell_type, y = Gene, fill = Correlation)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.15) +
    ggplot2::geom_point(data = tab[tab$Significant, , drop = FALSE],
      ggplot2::aes(x = Cell_type, y = Gene), inherit.aes = FALSE,
      size = 0.5, shape = 16) +
    ggplot2::scale_fill_gradient2(low = "#4267a6", mid = "white", high = "#bb3f42",
                                  midpoint = 0, limits = c(-1, 1)) +
    ggplot2::labs(title = paste0(label, ": DEG–LM22 correlation audit"),
                  subtitle = "Spearman rho; black dot: |rho| > 0.3 & BH FDR < 0.05",
                  x = "Immune cell", y = "DEG", fill = "rho") +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 55, hjust = 1),
                   panel.grid = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"))
  ggplot2::ggsave(file.path(out_dir, file_name), fig,
                  width = 14, height = 13, dpi = 240, limitsize = FALSE)
}
make_heatmap(all_pd, "PD", "05_PD_diagnostic_heatmap.png")
make_heatmap(all_oscc, "OSCC", "05_OSCC_diagnostic_heatmap.png")

writeLines(c(
  paste0("Date: ", Sys.Date()),
  "Step 05: manuscript-aligned Spearman correlation audit",
  "Historical DEG input: Script 02 (unblocked) 57 PD / 1266 OSCC",
  "OSCC correlation subset: top 50 DEG entries by absolute logFC",
  "Historical gene/probe selection: first GPL570 Gene.symbol before ///; highest within-dataset mean expression",
  "Source CIBERSORTx: archived LM22 Job29 / Job30, QC zero fraction <=50%",
  "Samples: 247 PD / 212 OSCC; ordered by GSM identifiers",
  "Primary tests: Spearman cor.test as historical Script 05 (R default exact handling)",
  "Primary significance: |rho| > 0.3 AND BH-adjusted p < 0.05 across ALL tested gene-cell pairs per disease",
  "Historic workbook BH applied after |rho| filter: saved only as comparison, not the primary FDR",
  "PD biopsies repeated for 90 patients: historic unblocked analysis; dependence is a limitation",
  "Original Figure 5 itself requires separate figure-level visual audit; output heatmaps are diagnostic",
  "Do not alter gene selection, threshold or p-values merely to match manuscript counts"
), file.path(out_dir, "05_analysis_settings.txt"))
writeLines(capture.output(utils::sessionInfo()),
           file.path(out_dir, "05_sessionInfo.txt"))

cat("\n=== STEP 05 COUNTS ===\n")
print(summary, row.names = FALSE)
cat("\n=== ARCHIVED CORRELATION COEFFICIENT QC ===\n")
print(reference_qc, row.names = FALSE)
if (!all(summary$Match_manuscript_count) ||
    (all(reference_qc$Reference_available) &&
     !all(reference_qc$Coefficients_match_1e_6))) {
  warning("Script 05 is not yet fully reconciled with the archived Figure 5 data; ",
          "review mismatch CSVs and probe selection. Do NOT force the counts.")
}
cat("STEP 05 COMPLETE. Outputs: ", out_dir, "\n", sep = "")
