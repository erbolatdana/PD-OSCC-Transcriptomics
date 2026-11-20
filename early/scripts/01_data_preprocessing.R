# =============================================================================
# 01_data_preprocessing.R
# PD–OSCC transcriptomic study | reproducibility candidate for GitHub
# Source datasets: GSE10334 (periodontitis), GSE30784 (OSCC), GPL570
#
# RUN FROM THE REPOSITORY ROOT:
#   source("scripts/01_data_preprocessing.R")
#
# This stage preserves probe-level matrices for the historical DEG analysis,
# and separately prepares non-ambiguous, single-probe-per-gene matrices.
# Do not substitute one for the other without recalculating downstream results.
# =============================================================================

required_packages <- c("GEOquery", "Biobase", "limma")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R/Bioconductor packages: ", paste(missing_packages, collapse = ", "),
       ". Install them before running the analysis.")
}

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
metadata_dir <- file.path(root, "data", "metadata")
output_dir <- file.path(root, "results", "01_preprocessing")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
options(timeout = max(300, getOption("timeout", 60)))

# The audited sample tables retain public GEO GSM IDs and, for PD, anonymized
# patient IDs. Their group assignments were cross-checked against the original
# GEO2R binary strings; no arbitrary ordering assumptions are needed below.
read_sample_map <- function(gse, allowed) {
  filename <- file.path(metadata_dir, paste0(gse, "_sample_groups.csv"))
  if (!file.exists(filename)) stop("Missing sample metadata: ", filename)
  tab <- utils::read.csv(filename, stringsAsFactors = FALSE)
  if (!all(c("GSM", "Group") %in% names(tab))) {
    stop("Expected GSM and Group columns in ", filename)
  }
  if (anyNA(tab$GSM) || anyDuplicated(tab$GSM) ||
      anyNA(tab$Group) || any(!tab$Group %in% allowed)) {
    stop("Invalid or duplicated sample IDs/group labels in ", filename)
  }
  tab
}

get_gpl570 <- function(gse) {
  cat("Downloading GEO series: ", gse, "\n", sep = "")
  sets <- GEOquery::getGEO(gse, GSEMatrix = TRUE, AnnotGPL = TRUE)
  if (!is.list(sets)) sets <- list(sets)
  platforms <- vapply(sets, Biobase::annotation, character(1))
  idx <- which(platforms == "GPL570")
  if (length(idx) != 1L) {
    stop(gse, ": expected exactly one GPL570 ExpressionSet; found: ",
         paste(platforms, collapse = ", "))
  }
  eset <- sets[[idx]]
  eset <- Biobase::`fvarLabels<-`(
    eset, value = make.names(Biobase::fvarLabels(eset))
  )
  if (!"Gene.symbol" %in% names(Biobase::fData(eset))) {
    stop("GPL570 annotation does not include Gene.symbol for ", gse)
  }
  eset
}

# Mirrors the transformation and between-array normalization of the original
# scripts used for the reported historical DEG totals. These are NOT the same
# as a workflow that simply retains author-processed GEO series values.
transform_original_method <- function(eset) {
  x <- Biobase::exprs(eset)
  qx <- as.numeric(stats::quantile(
    x, c(0, 0.25, 0.50, 0.75, 0.99, 1), na.rm = TRUE
  ))
  needs_log2 <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
  if (needs_log2) {
    x[x <= 0] <- NaN
    x <- log2(x)
  }
  # Retained intentionally for numerical comparison to the original analysis.
  x <- limma::normalizeBetweenArrays(x)
  eset <- Biobase::`exprs<-`(eset, value = x)
  eset <- eset[stats::complete.cases(Biobase::exprs(eset)), ]
  cat("  log2 conversion applied: ", needs_log2, "\n", sep = "")
  cat("  retained probes: ", nrow(eset), "\n", sep = "")
  eset
}

align_groups <- function(eset, sample_map, included, group_levels, gse) {
  gsm <- Biobase::sampleNames(eset)
  if (anyDuplicated(gsm) || !setequal(gsm, sample_map$GSM) ||
      length(gsm) != nrow(sample_map)) {
    stop(gse, ": GEO sample IDs do not match the committed sample map")
  }
  sample_map <- sample_map[match(gsm, sample_map$GSM), , drop = FALSE]
  stopifnot(identical(gsm, sample_map$GSM))
  keep <- sample_map$Group %in% included
  eset <- eset[, keep]
  sample_map <- sample_map[keep, , drop = FALSE]
  groups <- factor(sample_map$Group, levels = group_levels)
  if (anyNA(groups)) stop(gse, ": unknown group label")
  eset$group <- groups
  if (!identical(Biobase::sampleNames(eset), sample_map$GSM)) {
    stop(gse, ": sample order changed unexpectedly")
  }
  list(eset = eset, sample_map = sample_map, groups = groups)
}

# --- 1. PD: 183 disease sites + 64 control sites -----------------------------
map_pd_raw <- read_sample_map("GSE10334", c("PD", "Control"))
if (!"Patient_ID" %in% names(map_pd_raw) ||
    anyNA(map_pd_raw$Patient_ID) || any(!nzchar(map_pd_raw$Patient_ID))) {
  stop("GSE10334 must include a non-missing Patient_ID per biopsy")
}
if (nrow(map_pd_raw) != 247L || length(unique(map_pd_raw$Patient_ID)) != 90L) {
  stop("GSE10334 sample/patient total changed; review clinical metadata")
}
raw_pd <- get_gpl570("GSE10334")
if (ncol(raw_pd) != 247L) stop("GSE10334: expected 247 GEO samples")
gset_pd <- transform_original_method(raw_pd)
pd_aligned <- align_groups(gset_pd, map_pd_raw,
                           included = c("PD", "Control"),
                           group_levels = c("PD", "Control"),
                           gse = "GSE10334")
gset_pd <- pd_aligned$eset
gs_pd <- pd_aligned$groups
sample_map_pd <- pd_aligned$sample_map
if (sum(gs_pd == "PD") != 183L || sum(gs_pd == "Control") != 64L) {
  stop("GSE10334: unexpected PD/control group sizes")
}

# --- 2. OSCC: 167 tumors + 45 controls; exclude 17 dysplasias ----------------
map_oscc_raw <- read_sample_map("GSE30784", c("OSCC", "Control", "Dysplasia"))
raw_oscc <- get_gpl570("GSE30784")
if (ncol(raw_oscc) != 229L) stop("GSE30784: expected 229 original GEO samples")
# Important for replication: normalization precedes dysplasia exclusion,
# exactly as in the historical Script 01.
gset_oscc <- transform_original_method(raw_oscc)
oscc_aligned <- align_groups(gset_oscc, map_oscc_raw,
                             included = c("OSCC", "Control"),
                             group_levels = c("Control", "OSCC"),
                             gse = "GSE30784")
gset_oscc <- oscc_aligned$eset
gs_oscc <- oscc_aligned$groups
sample_map_oscc <- oscc_aligned$sample_map
if (sum(gs_oscc == "OSCC") != 167L ||
    sum(gs_oscc == "Control") != 45L ||
    sum(map_oscc_raw$Group == "Dysplasia") != 17L) {
  stop("GSE30784: unexpected OSCC/control/dysplasia group sizes")
}

# --- 3. Preserve probe-level matrices for the historical DEG model ----------
mat_pd <- Biobase::exprs(gset_pd)
mat_oscc <- Biobase::exprs(gset_oscc)
if (!identical(rownames(mat_pd), rownames(mat_oscc)) ||
    !identical(rownames(mat_pd), rownames(Biobase::fData(gset_pd))) ||
    !identical(rownames(mat_oscc), rownames(Biobase::fData(gset_oscc)))) {
  stop("Probe IDs/annotation rows do not align across datasets")
}
if (!identical(as.character(Biobase::fData(gset_pd)[["Gene.symbol"]]),
               as.character(Biobase::fData(gset_oscc)[["Gene.symbol"]]))) {
  stop("GPL570 gene-symbol annotations differ between the two series")
}
if (!identical(dim(mat_pd), c(54675L, 247L)) ||
    !identical(dim(mat_oscc), c(54675L, 212L))) {
  stop("Expression matrices differ from the audited 54675 x 247/212 layout")
}

# --- 4. Separately create a strictly one-gene-one-probe representation ------
# This does NOT alter probe-level differential expression testing in Script 02.
# Ambiguous probe-to-multiple-gene mappings are omitted, rather than choosing
# the first (potentially incorrect) symbol. Probe choice here uses expression
# only and a deterministic tie-break, not clinical labels.
genes <- trimws(as.character(Biobase::fData(gset_pd)[["Gene.symbol"]]))
valid <- !is.na(genes) & nzchar(genes) & genes != "---" &
  !grepl("///", genes, fixed = TRUE)
valid_idx <- which(valid)
combined_mean <- (rowSums(mat_pd) + rowSums(mat_oscc)) /
  (ncol(mat_pd) + ncol(mat_oscc))
ordered_idx <- valid_idx[order(-combined_mean[valid_idx],
                               rownames(mat_pd)[valid_idx])]
selected <- ordered_idx[!duplicated(genes[ordered_idx])]
mat_pd_gene <- mat_pd[selected, , drop = FALSE]
mat_oscc_gene <- mat_oscc[selected, , drop = FALSE]
rownames(mat_pd_gene) <- genes[selected]
rownames(mat_oscc_gene) <- genes[selected]
if (anyDuplicated(rownames(mat_pd_gene)) ||
    !identical(rownames(mat_pd_gene), rownames(mat_oscc_gene)) ||
    nrow(mat_pd_gene) != 20848L) {
  stop("Gene-mapping QC differs from observed 20,848 unique symbols")
}
probe_gene_map <- data.frame(
  Gene = genes[selected], GPL570_Probe_ID = rownames(mat_pd)[selected],
  stringsAsFactors = FALSE
)

# Gene-level linear matrices are an exploratory starting point for rechecking
# CIBERSORTx; they have NOT been proven identical to the uploaded Job29/30 files.
mat_pd_gene_linear <- 2^mat_pd_gene
mat_oscc_gene_linear <- 2^mat_oscc_gene
# Explicitly avoid the ambiguous legacy variable 'mat_pd_linear': downstream
# code must choose *_gene_linear when using LM22 (not probe IDs).

# --- 5. Save reproducible hand-off files ------------------------------------
preprocessed <- list(
  gset_pd = gset_pd, gset_oscc = gset_oscc,
  gs_pd = gs_pd, gs_oscc = gs_oscc,
  sample_map_pd = sample_map_pd, sample_map_oscc = sample_map_oscc
)
saveRDS(preprocessed, file.path(output_dir, "01_preprocessed_esets.rds"))
saveRDS(list(
  mat_pd_gene = mat_pd_gene, mat_oscc_gene = mat_oscc_gene,
  mat_pd_gene_linear = mat_pd_gene_linear,
  mat_oscc_gene_linear = mat_oscc_gene_linear
), file.path(output_dir, "01_gene_level_matrices.rds"))
utils::write.csv(probe_gene_map,
                 file.path(output_dir, "01_representative_probe_map.csv"),
                 row.names = FALSE)
utils::write.csv(sample_map_pd, file.path(output_dir, "01_PD_sample_map.csv"),
                 row.names = FALSE)
utils::write.csv(sample_map_oscc,
                 file.path(output_dir, "01_OSCC_sample_map.csv"), row.names = FALSE)

qc <- data.frame(
  Dataset = c("GSE10334", "GSE30784"),
  Disease = c("PD", "OSCC"),
  Samples = c(ncol(mat_pd), ncol(mat_oscc)),
  Disease_n = c(183L, 167L),
  Control_n = c(64L, 45L),
  Probe_rows = c(nrow(mat_pd), nrow(mat_oscc)),
  Unique_gene_rows = c(nrow(mat_pd_gene), nrow(mat_oscc_gene)),
  stringsAsFactors = FALSE
)
utils::write.csv(qc, file.path(output_dir, "01_QC_summary.csv"), row.names = FALSE)
writeLines(capture.output(utils::sessionInfo()),
           file.path(output_dir, "01_sessionInfo.txt"))
print(qc, row.names = FALSE)
cat("STEP 01 COMPLETE. Output folder: ", output_dir, "\n", sep = "")
