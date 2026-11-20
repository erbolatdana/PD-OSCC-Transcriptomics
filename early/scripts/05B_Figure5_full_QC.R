# =============================================================================
# 05B_Figure5_full_QC.R
# Compare every cell of the original Figure 5 correlation matrices with
# independent Script 05 results (all 912 PD / 700 OSCC tested pairs).
# Run from the repository root:
#   First set the working directory to this project folder in RStudio.
#   source("scripts/05B_Figure5_full_QC.R")
# Base R only. Does not rerun or change scripts 01-05 or figure source matrices.
# =============================================================================

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "results", "05_DEG_immune_correlation")
ref_dir <- file.path(root, "data", "reference")

specs <- list(
  PD = list(figure = file.path(ref_dir, "Figure5_PD_correlation_matrix.csv"),
            calculated = file.path(out_dir, "05_PD_all_912_pairs.csv"),
            expected_genes = 16L, expected_cells = 16L,
            expected_all_pairs = 912L, expected_sig = 346L),
  OSCC = list(figure = file.path(ref_dir, "Figure5_OSCC_correlation_matrix.csv"),
              calculated = file.path(out_dir, "05_OSCC_all_700_pairs.csv"),
              expected_genes = 14L, expected_cells = 14L,
              expected_all_pairs = 700L, expected_sig = 180L)
)

missing <- unlist(lapply(specs, function(x) c(x$figure, x$calculated)))
missing <- missing[!file.exists(missing)]
if (length(missing)) stop("Missing files:\n", paste(missing, collapse="\n"))

reports <- list()
comparisons <- list()
for (disease in names(specs)) {
  x <- specs[[disease]]
  fig <- read.csv(x$figure, check.names = FALSE, stringsAsFactors = FALSE)
  names(fig)[1] <- "Gene"
  if (nrow(fig) != x$expected_genes || (ncol(fig)-1L) != x$expected_cells) {
    stop(disease, ": unexpected Figure 5 matrix dimensions")
  }
  if (anyNA(fig$Gene) || anyDuplicated(fig$Gene) ||
      anyDuplicated(names(fig)[-1L])) {
    stop(disease, ": duplicate/missing gene or immune-cell label")
  }
  fig_long <- data.frame(
    Gene = rep(fig$Gene, times = ncol(fig)-1L),
    Cell_type = rep(names(fig)[-1L], each = nrow(fig)),
    Figure_rho = as.numeric(as.matrix(fig[-1L])),
    stringsAsFactors = FALSE
  )
  if (anyNA(fig_long$Figure_rho) || any(!is.finite(fig_long$Figure_rho))) {
    stop(disease, ": figure matrix has missing or nonfinite correlation")
  }
  all_tested <- read.csv(x$calculated, check.names=FALSE,
                         stringsAsFactors=FALSE)
  needed <- c("Gene", "Cell_type", "Correlation", "p_adj", "Significant")
  if (!all(needed %in% names(all_tested)) ||
      nrow(all_tested) != x$expected_all_pairs ||
      sum(as.logical(all_tested$Significant)) != x$expected_sig) {
    stop(disease, ": Script 05 output missing/has unexpected test counts")
  }
  key <- function(df) paste(df$Gene, df$Cell_type, sep="\r")
  if (anyDuplicated(key(all_tested)) || anyDuplicated(key(fig_long))) {
    stop(disease, ": duplicate pair keys")
  }
  loc <- match(key(fig_long), key(all_tested))
  if (anyNA(loc)) {
    stop(disease, ": Figure 5 pair(s) absent from Script 05:\n",
         paste(key(fig_long)[is.na(loc)], collapse="\n"))
  }
  fig_long$Reproduced_rho <- all_tested$Correlation[loc]
  fig_long$Reproduced_BH_FDR <- all_tested$p_adj[loc]
  fig_long$Reproduced_significant <- as.logical(all_tested$Significant[loc])
  fig_long$Abs_difference <- abs(fig_long$Figure_rho - fig_long$Reproduced_rho)
  fig_long$Dataset <- disease
  max_delta <- max(fig_long$Abs_difference)
  threshold_agreement <- all((abs(fig_long$Figure_rho) > .30) ==
                             fig_long$Reproduced_significant)
  reports[[disease]] <- data.frame(
    Dataset = disease,
    Figure_genes = nrow(fig),
    Figure_cell_types = ncol(fig) - 1L,
    Figure_cells_checked = nrow(fig_long),
    Figure_abs_r_above_0_3 = sum(abs(fig_long$Figure_rho) > .30),
    Corresponding_BH_significant_pairs = sum(fig_long$Reproduced_significant),
    Max_abs_rho_difference = max_delta,
    All_rho_match_1e_10 = max_delta <= 1e-10,
    Significance_labels_agree = threshold_agreement,
    stringsAsFactors = FALSE
  )
  comparisons[[disease]] <- fig_long
}

summary <- do.call(rbind, reports)
detail <- do.call(rbind, comparisons)
write.csv(summary, file.path(out_dir, "05_Figure5_full_QC_summary.csv"), row.names=FALSE)
write.csv(detail, file.path(out_dir, "05_Figure5_full_QC_pairwise.csv"), row.names=FALSE)
cat("\n=== FIGURE 5 FULL-CELL QC ===\n")
print(summary, row.names=FALSE)
cat("\nTotal Figure 5 cells checked:", sum(summary$Figure_cells_checked), "\n")
if (!all(summary$All_rho_match_1e_10) ||
    !all(summary$Significance_labels_agree)) {
  stop("Figure 5 QC FAILED. Inspect 05_Figure5_full_QC_pairwise.csv")
}
cat("FIGURE 5 FULL-CELL QC PASSED. No Script 05 or Figure files modified.\n")
