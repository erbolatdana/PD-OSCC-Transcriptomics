# =============================================================================
# 02_DEG_analysis.R
# PD–OSCC transcriptomic study | reproducibility candidate for GitHub
#
# RUN FROM THE REPOSITORY ROOT, AFTER SCRIPT 01:
#   source("scripts/02_DEG_analysis.R")
#
# DEFAULT reproduces the historical UNBLOCKED probe-level limma analysis
# that returned 119 PD / 1874 OSCC significant probes, subsequently
# deduplicated to 57 PD / 1266 OSCC reported gene-symbol entries.
#
# IMPORTANT STATISTICAL LIMITATION:
# GSE10334 includes multiple biopsies from some patients. The original model
# does not adjust for within-patient correlation. To run a separate, correctly
# labeled exploratory sensitivity analysis, set:
#   Sys.setenv(PD_DEG_MODEL = "patient_blocked")
#   source("scripts/02_DEG_analysis.R")
# This alternative writes to a DIFFERENT folder and is NOT required to match
# the historical results. Manuscript methods must accurately disclose which
# model produced which figures and numbers.
# =============================================================================

required_packages <- c("Biobase", "limma", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

LOGFC_THRESHOLD <- 1.5
FDR_THRESHOLD <- 0.05
PD_MODEL <- Sys.getenv("PD_DEG_MODEL", unset = "historical_unblocked")
if (!PD_MODEL %in% c("historical_unblocked", "patient_blocked")) {
  stop("PD_DEG_MODEL must be 'historical_unblocked' or 'patient_blocked'")
}

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
input_file <- file.path(root, "results", "01_preprocessing", "01_preprocessed_esets.rds")
if (!file.exists(input_file)) {
  stop("Missing Script 01 output: ", input_file,
       "\nRun source('scripts/01_data_preprocessing.R') first.")
}
input <- readRDS(input_file)
stopifnot(all(c("gset_pd", "gset_oscc", "gs_pd", "gs_oscc",
                "sample_map_pd", "sample_map_oscc") %in% names(input)))
gset_pd <- input$gset_pd
gset_oscc <- input$gset_oscc
gs_pd <- input$gs_pd
gs_oscc <- input$gs_oscc
sample_map_pd <- input$sample_map_pd
sample_map_oscc <- input$sample_map_oscc

# This explicitly reloads Step 01 data, not objects from a previously
# populated interactive R workspace (avoids accidental reuse of stale data).
if (!identical(Biobase::sampleNames(gset_pd), sample_map_pd$GSM) ||
    !identical(Biobase::sampleNames(gset_oscc), sample_map_oscc$GSM)) {
  stop("Sample order differs from the Step 01 audited metadata")
}
if (!identical(as.character(gset_pd$group), as.character(gs_pd)) ||
    !identical(as.character(gset_oscc$group), as.character(gs_oscc))) {
  stop("ExpressionSet group labels do not match saved factors")
}

if (PD_MODEL == "historical_unblocked") {
  output_dir <- file.path(root, "results", "02_DEG")
} else {
  output_dir <- file.path(root, "results", "02_DEG_patient_blocked_sensitivity")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Use exactly the historical probe-wise model for the primary replication.
# Under patient_blocked, only the PD model is adjusted; the OSCC model remains
# standard limma and therefore serves as a common comparison.
run_limma <- function(eset, disease, control, group_levels, block = NULL) {
  design <- stats::model.matrix(~group + 0, eset)
  colnames(design) <- group_levels
  contrast <- limma::makeContrasts(
    contrasts = paste0(disease, " - ", control), levels = design
  )

  within_patient_correlation <- NA_real_
  if (is.null(block)) {
    fit <- limma::lmFit(eset, design)
  } else {
    if (length(block) != ncol(eset) || anyNA(block) ||
        any(!nzchar(as.character(block)))) {
      stop("Missing/misaligned PD patient IDs for duplicateCorrelation")
    }
    if (!requireNamespace("statmod", quietly = TRUE)) {
      stop("Install the statmod package to use patient_blocked mode")
    }
    patient <- factor(block)
    cor_fit <- limma::duplicateCorrelation(
      Biobase::exprs(eset), design, block = patient
    )
    within_patient_correlation <- cor_fit$consensus.correlation
    if (!is.finite(within_patient_correlation)) {
      stop("Patient-level consensus correlation was not finite")
    }
    fit <- limma::lmFit(
      eset, design, block = patient,
      correlation = within_patient_correlation
    )
  }

  # eBayes(..., proportion=0.01) is the setting in the historical Script 02.
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, contrast), proportion = 0.01)
  tT <- limma::topTable(fit2, adjust = "fdr", sort.by = "B", number = Inf)
  decision <- limma::decideTests(fit2, adjust.method = "fdr",
                                p.value = FDR_THRESHOLD,
                                lfc = LOGFC_THRESHOLD)
  list(table = tT, decision = decision,
       patient_correlation = within_patient_correlation)
}

pd_block <- NULL
if (PD_MODEL == "patient_blocked") {
  if (!"Patient_ID" %in% names(sample_map_pd)) {
    stop("Step 01 PD sample map has no Patient_ID")
  }
  pd_block <- sample_map_pd$Patient_ID
  if (length(unique(pd_block)) != 90L) {
    stop("PD patient clusters differ from the 90 audited patients")
  }
}

pd_fit <- run_limma(gset_pd, disease = "PD", control = "Control",
                    group_levels = levels(gs_pd), block = pd_block)
oscc_fit <- run_limma(gset_oscc, disease = "OSCC", control = "Control",
                      group_levels = levels(gs_oscc))
tT_pd <- pd_fit$table
tT_oscc <- oscc_fit$table
dT_pd <- pd_fit$decision
dT_oscc <- oscc_fit$decision

cat("=== DEG results: PD (GSE10334), model: ", PD_MODEL, " ===\n", sep = "")
cat("Total significant probes: ", sum(dT_pd != 0), "\n", sep = "")
cat("Up-regulated: ", sum(dT_pd == 1), "\n", sep = "")
cat("Down-regulated: ", sum(dT_pd == -1), "\n\n", sep = "")
cat("=== DEG results: OSCC (GSE30784) ===\n")
cat("Total significant probes: ", sum(dT_oscc != 0), "\n", sep = "")
cat("Up-regulated: ", sum(dT_oscc == 1), "\n", sep = "")
cat("Down-regulated: ", sum(dT_oscc == -1), "\n\n", sep = "")

# Maintain precisely the historical deduplication procedure AFTER testing:
# 1) filter significant PROBES using BH FDR and |log2FC|;
# 2) take FIRST symbol if the Affymetrix annotation is 'A /// B';
# 3) retain the significant probe with the highest |log2FC| per symbol.
# The procedure is NOT "highest mean expression" and MUST NOT be described
# that way in the manuscript or README. It may overestimate effects because
# a representative probe is picked after testing; a pre-specified gene-level
# analysis is scientifically distinct and requires recomputation.
first_annotated_symbol <- function(annotation) {
  ann <- as.character(annotation)
  pieces <- strsplit(ann, "///", fixed = TRUE)
  vapply(pieces, function(x) {
    if (length(x) == 0L || is.na(x[1L])) NA_character_ else trimws(x[1L])
  }, character(1))
}

deduplicate_as_original <- function(tT, gse) {
  if (!"Gene.symbol" %in% names(tT)) {
    stop("Missing Gene.symbol in ", gse, " topTable; check GEO annotation")
  }
  good <- is.finite(tT$adj.P.Val) & tT$adj.P.Val < FDR_THRESHOLD &
    is.finite(tT$logFC) & abs(tT$logFC) >= LOGFC_THRESHOLD
  selected <- tT[good, , drop = FALSE]
  selected$Gene <- first_annotated_symbol(selected$Gene.symbol)
  selected <- selected[!is.na(selected$Gene) & selected$Gene != "",,
                       drop = FALSE]
  selected <- selected[order(-abs(selected$logFC), seq_len(nrow(selected))),,
                       drop = FALSE]
  selected <- selected[!duplicated(selected$Gene), , drop = FALSE]
  if (any(selected$Gene %in% c("---", "NA"))) {
    warning(gse, ": non-gene placeholder(s) are present in the historical DEG list")
  }
  selected
}

deg_pd_clean <- deduplicate_as_original(tT_pd, "GSE10334")
deg_oscc_clean <- deduplicate_as_original(tT_oscc, "GSE30784")
deg_oscc_full <- deg_oscc_clean$Gene  # compatibility with historical Script 06

cat("Deduplicated unique DEG entries — PD: ", nrow(deg_pd_clean),
    " | OSCC: ", nrow(deg_oscc_clean), "\n", sep = "")

make_volcano <- function(tT, title, n_up, n_down, path) {
  df <- tT
  df$Gene <- first_annotated_symbol(df$Gene.symbol)
  df$Status <- "NS"
  df$Status[df$logFC >= LOGFC_THRESHOLD & df$adj.P.Val < FDR_THRESHOLD] <- "Up"
  df$Status[df$logFC <= -LOGFC_THRESHOLD & df$adj.P.Val < FDR_THRESHOLD] <- "Down"
  df$Status <- factor(df$Status, levels = c("Up", "Down", "NS"))
  df$neg_log10_p <- -log10(pmax(df$adj.P.Val, .Machine$double.xmin))
  p <- ggplot2::ggplot(df,
    ggplot2::aes(x = logFC, y = neg_log10_p, color = Status)) +
    ggplot2::geom_point(alpha = 0.5, size = 0.8) +
    ggplot2::scale_color_manual(
      values = c(Up = "#E63946", Down = "#457B9D", NS = "grey70"),
      labels = c(paste0("Up (n=", n_up, ")"),
                 paste0("Down (n=", n_down, ")"), "NS")) +
    ggplot2::geom_vline(xintercept = c(-LOGFC_THRESHOLD, LOGFC_THRESHOLD),
                        linetype = "dashed", color = "grey40") +
    ggplot2::geom_hline(yintercept = -log10(FDR_THRESHOLD),
                        linetype = "dashed", color = "grey40") +
    ggplot2::labs(title = title, x = "log2 (Fold Change)",
                  y = "-log10 (adj. P-value)", color = "") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(
      hjust = 0.5, face = "bold", size = 13),
      legend.position = "top")
  grDevices::png(path, width = 2400, height = 2000, res = 300)
  tryCatch(print(p), finally = grDevices::dev.off())
}

utils::write.csv(tT_pd, file.path(output_dir, "02_PD_all_probes.csv"))
utils::write.csv(tT_oscc, file.path(output_dir, "02_OSCC_all_probes.csv"))
utils::write.csv(deg_pd_clean,
                 file.path(output_dir, "02_PD_unique_DEG_entries.csv"), row.names = TRUE)
utils::write.csv(deg_oscc_clean,
                 file.path(output_dir, "02_OSCC_unique_DEG_entries.csv"), row.names = TRUE)
make_volcano(tT_pd, "Periodontitis (GSE10334)", sum(dT_pd == 1),
             sum(dT_pd == -1), file.path(output_dir, "Volcano_PD_GSE10334.png"))
make_volcano(tT_oscc, "OSCC (GSE30784)", sum(dT_oscc == 1),
             sum(dT_oscc == -1), file.path(output_dir, "Volcano_OSCC_GSE30784.png"))

summary <- data.frame(
  Dataset = c("GSE10334", "GSE30784"),
  Disease = c("PD", "OSCC"),
  Disease_n = c(183L, 167L), Control_n = c(64L, 45L),
  Model = c(PD_MODEL, "unblocked"),
  Significant_probes = c(sum(dT_pd != 0), sum(dT_oscc != 0)),
  Up_probes = c(sum(dT_pd == 1), sum(dT_oscc == 1)),
  Down_probes = c(sum(dT_pd == -1), sum(dT_oscc == -1)),
  Deduplicated_DEG_entries = c(nrow(deg_pd_clean), nrow(deg_oscc_clean)),
  PD_consensus_correlation = c(pd_fit$patient_correlation, NA_real_),
  stringsAsFactors = FALSE
)
utils::write.csv(summary, file.path(output_dir, "02_DEG_summary.csv"),
                 row.names = FALSE)
saveRDS(list(tT_pd = tT_pd, tT_oscc = tT_oscc,
             deg_pd_clean = deg_pd_clean,
             deg_oscc_clean = deg_oscc_clean,
             deg_oscc_full = deg_oscc_full,
             dT_pd = dT_pd, dT_oscc = dT_oscc,
             model = PD_MODEL),
        file.path(output_dir, "02_DEG_objects.rds"))
writeLines(capture.output(utils::sessionInfo()),
           file.path(output_dir, "02_sessionInfo.txt"))
print(summary, row.names = FALSE)

# Verification against the already-observed original Script 02 results.
# These are historical comparison values, not targets to tune the model to.
expected <- c(pd_probes = 119L, pd_up = 113L, pd_down = 6L,
              oscc_probes = 1874L, oscc_up = 895L, oscc_down = 979L,
              pd_unique = 57L, oscc_unique = 1266L)
observed <- c(pd_probes = sum(dT_pd != 0),
              pd_up = sum(dT_pd == 1), pd_down = sum(dT_pd == -1),
              oscc_probes = sum(dT_oscc != 0),
              oscc_up = sum(dT_oscc == 1), oscc_down = sum(dT_oscc == -1),
              pd_unique = nrow(deg_pd_clean),
              oscc_unique = nrow(deg_oscc_clean))
if (PD_MODEL == "historical_unblocked") {
  qc <- data.frame(Metric = names(expected),
                   Expected = unname(expected),
                   Observed = unname(observed),
                   Match = unname(expected == observed))
  utils::write.csv(qc, file.path(output_dir, "02_historical_count_QC.csv"),
                   row.names = FALSE)
  print(qc, row.names = FALSE)
  if (!all(qc$Match)) {
    stop("Historical DEG counts do not match: review the QC CSV and inputs. ",
         "Do not modify thresholds merely to force a match.")
  }
  cat("HISTORICAL COUNT CHECK PASSED (pending independent clean-session run).\n")
} else {
  cat("PATIENT-BLOCKED sensitivity complete. Differences from the historical",
      " PD counts are expected and must be reported transparently.\n")
}
cat("STEP 02 COMPLETE. Output folder: ", output_dir, "\n", sep = "")
