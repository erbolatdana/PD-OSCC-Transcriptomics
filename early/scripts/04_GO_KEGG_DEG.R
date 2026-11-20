# =============================================================================
# 04_GO_KEGG_DEG.R
# PD–OSCC transcriptomic study | GitHub reproducibility CANDIDATE
#
# RUN FROM PROJECT ROOT after successful Scripts 01 and 02:
#   First set the working directory to this project folder in RStudio.
#   source("scripts/04_GO_KEGG_DEG.R")
#
# PURPOSE: independently re-run the ORIGINAL published DEG-list GO Biological
# Process / KEGG over-representation analysis using clusterProfiler. No
# workspace objects are required; read Script 02's historical, unblocked
# DEG CSVs. The patient-blocked PD analysis is a separate sensitivity model,
# and MUST NOT silently replace the historical 57-gene input used in Figure 4.
#
# MANUSCRIPT METHODS: clusterProfiler + org.Hs.eg.db; SYMBOL -> ENTREZID;
# enrichGO(ont="BP") and enrichKEGG(organism="hsa"); BH; adjusted p < 0.05.
# ORIGINAL SCRIPT 04 DEFAULTS are preserved explicitly: pvalueCutoff=0.05,
# qvalueCutoff=0.2, minGSSize=10, maxGSSize=500, and no user-supplied universe.
# The annotation DB / GO / live KEGG pathway definitions may change over time;
# recorded package versions and input MD5s are therefore essential. A changed
# count is a discrepancy to investigate, NOT a target to force by tuning.
#
# HISTORICAL CHECKS FROM REVISED MANUSCRIPT (not fabricated output):
#   Input DEG entries: PD 57; OSCC 1266
#   Significant GO BP terms: PD 16; OSCC 788
#   Shared significant KEGG: hsa04060, hsa04061 and hsa05323
#   Reported unmapped SYMBOL proportions: ~7.0% PD and ~5.9% OSCC.
# The last two percentages are descriptive references, not exact assertions.
# =============================================================================

# --- 0. Environment and inputs --------------------------------------------
required_packages <- c("clusterProfiler", "org.Hs.eg.db", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "),
       ". Install them before running this script.")
}

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "results", "04_GO_KEGG")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pd_deg_file <- file.path(root, "results", "02_DEG",
                         "02_PD_unique_DEG_entries.csv")
oscc_deg_file <- file.path(root, "results", "02_DEG",
                           "02_OSCC_unique_DEG_entries.csv")
required_files <- c(pd_deg_file, oscc_deg_file)
if (!all(file.exists(required_files))) {
  stop("Missing Script 02 historical DEG CSV(s):\n",
       paste(required_files[!file.exists(required_files)], collapse = "\n"),
       "\nRun the historical/unblocked Script 02 first.")
}

FDR_THRESHOLD <- 0.05
GO_TERMS_PUBLISHED <- c(PD = 16L, OSCC = 788L)
DEG_ENTRIES_PUBLISHED <- c(PD = 57L, OSCC = 1266L)

read_deg_symbols <- function(path, label) {
  tab <- utils::read.csv(path, stringsAsFactors = FALSE,
                         check.names = FALSE)
  if (!"Gene" %in% names(tab)) {
    stop(label, ": missing 'Gene' column in Script 02 CSV: ", path)
  }
  if (anyNA(tab$Gene) || any(!nzchar(trimws(as.character(tab$Gene))))) {
    stop(label, ": missing or empty gene-symbol entries: ", path)
  }
  symbols <- trimws(as.character(tab$Gene))
  if (anyDuplicated(symbols)) {
    stop(label, ": duplicate symbols in Script 02 *deduplicated* DEG table")
  }
  symbols
}

pd_symbols <- read_deg_symbols(pd_deg_file, "PD")
oscc_symbols <- read_deg_symbols(oscc_deg_file, "OSCC")
cat("Script 02 historical DEG entries: PD=", length(pd_symbols),
    " | OSCC=", length(oscc_symbols), "\n", sep = "")

# Record provenance before any external database calls.
input_manifest <- data.frame(
  File = c("02_PD_unique_DEG_entries.csv", "02_OSCC_unique_DEG_entries.csv"),
  MD5 = unname(tools::md5sum(required_files)),
  Size_bytes = unname(file.info(required_files)$size),
  stringsAsFactors = FALSE
)
utils::write.csv(input_manifest,
                 file.path(out_dir, "04_input_manifest.csv"), row.names = FALSE)

# --- 1. SYMBOL -> ENTREZID (preserve the original bitr semantics) -----------
# bitr can return multiple mappings per symbol; report coverage by DISTINCT
# input symbols, not by nrow(bitr(...)). Enrichment uses distinct Entrez IDs.
map_symbols <- function(symbols, label) {
  mapped <- suppressWarnings(clusterProfiler::bitr(
    symbols, fromType = "SYMBOL", toType = "ENTREZID",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  ))
  if (is.null(mapped) || !nrow(mapped)) {
    stop(label, ": no mapped Entrez IDs; inspect org.Hs.eg.db installation")
  }
  mapped <- unique(mapped[, c("SYMBOL", "ENTREZID"), drop = FALSE])
  mapped <- mapped[!is.na(mapped$SYMBOL) & !is.na(mapped$ENTREZID) &
                     nzchar(as.character(mapped$ENTREZID)), , drop = FALSE]
  unmapped <- setdiff(symbols, unique(as.character(mapped$SYMBOL)))
  entrez <- unique(as.character(mapped$ENTREZID))
  report <- data.frame(
    Dataset = label,
    Input_unique_symbols = length(symbols),
    Mapped_unique_symbols = length(setdiff(symbols, unmapped)),
    Unmapped_unique_symbols = length(unmapped),
    Unmapped_pct = round(100 * length(unmapped) / length(symbols), 2),
    Mapping_rows = nrow(mapped),
    Unique_Entrez_IDs = length(entrez),
    stringsAsFactors = FALSE
  )
  cat(label, ": ", report$Mapped_unique_symbols, "/",
      report$Input_unique_symbols, " SYMBOLs mapped (",
      report$Unmapped_pct, "% unmapped), ", length(entrez),
      " unique Entrez IDs\n", sep = "")
  utils::write.csv(mapped, file.path(out_dir,
       paste0("04_", label, "_SYMBOL_to_ENTREZ.csv")), row.names = FALSE)
  utils::write.csv(data.frame(Gene = unmapped), file.path(out_dir,
       paste0("04_", label, "_unmapped_SYMBOLs.csv")), row.names = FALSE)
  list(entrez = entrez, report = report, mapped = mapped)
}

map_pd <- map_symbols(pd_symbols, "PD")
map_oscc <- map_symbols(oscc_symbols, "OSCC")
mapping_summary <- rbind(map_pd$report, map_oscc$report)
utils::write.csv(mapping_summary,
                 file.path(out_dir, "04_mapping_QC.csv"), row.names = FALSE)

# --- 2. GO Biological Process enrichment ---------------------------------
# Explicitly retain the unmodified *original* enrichment settings:
# - default annotation universe (NO explicit universe argument)
# - default min/max gene-set sizes of 10/500
# - raw and adjusted p <= 0.05, and q <= 0.2 for returned enrichResult terms.
run_go <- function(entrez, label) {
  cat("Running GO BP enrichment for ", label, "...\n", sep = "")
  result <- clusterProfiler::enrichGO(
    gene = entrez,
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = FDR_THRESHOLD,
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500,
    readable = TRUE
  )
  saveRDS(result, file.path(out_dir, paste0("04_GO_", label, ".rds")))
  result
}

go_pd <- run_go(map_pd$entrez, "PD")
go_oscc <- run_go(map_oscc$entrez, "OSCC")

# --- 3. KEGG enrichment ---------------------------------------------------
# enrichKEGG uses LIVE KEGG annotation by default; Internet access required.
# Do not replace live KEGG results with an unrecorded offline annotation source.
run_kegg <- function(entrez, label) {
  cat("Running KEGG enrichment for ", label,
      " (may require access to the KEGG API)...\n", sep = "")
  result <- tryCatch(
    clusterProfiler::enrichKEGG(
      gene = entrez,
      organism = "hsa",
      keyType = "kegg",
      pAdjustMethod = "BH",
      pvalueCutoff = FDR_THRESHOLD,
      qvalueCutoff = 0.2,
      minGSSize = 10,
      maxGSSize = 500,
      use_internal_data = FALSE
    ),
    error = function(e) {
      stop("KEGG enrichment failed for ", label, ": ", conditionMessage(e),
           "\nGO outputs have been saved. Check the KEGG connection ",
           "and package versions; do not substitute new input lists.",
           call. = FALSE)
    }
  )
  saveRDS(result, file.path(out_dir, paste0("04_KEGG_", label, ".rds")))
  result
}

kegg_pd <- run_kegg(map_pd$entrez, "PD")
kegg_oscc <- run_kegg(map_oscc$entrez, "OSCC")

# --- 4. Export enrichment tables and compare pathway IDs -----------------
result_table <- function(obj, label) {
  if (is.null(obj)) {
    stop(label, ": enrichment returned NULL; inspect the input IDs/database")
  }
  tab <- as.data.frame(obj)
  if (!all(c("ID", "Description", "p.adjust") %in% names(tab))) {
    stop(label, ": unexpected enrichment result columns")
  }
  tab
}

filter_fdr <- function(tab) {
  tab[!is.na(tab$p.adjust) & is.finite(tab$p.adjust) &
        tab$p.adjust < FDR_THRESHOLD, , drop = FALSE]
}

go_pd_tab <- result_table(go_pd, "PD GO")
go_oscc_tab <- result_table(go_oscc, "OSCC GO")
kegg_pd_tab <- result_table(kegg_pd, "PD KEGG")
kegg_oscc_tab <- result_table(kegg_oscc, "OSCC KEGG")

for (spec in list(
  list(name = "04_GO_PD.csv", data = go_pd_tab),
  list(name = "04_GO_OSCC.csv", data = go_oscc_tab),
  list(name = "04_KEGG_PD.csv", data = kegg_pd_tab),
  list(name = "04_KEGG_OSCC.csv", data = kegg_oscc_tab)
)) {
  utils::write.csv(spec$data, file.path(out_dir, spec$name), row.names = FALSE)
}

go_pd_sig <- filter_fdr(go_pd_tab)
go_oscc_sig <- filter_fdr(go_oscc_tab)
kegg_pd_sig <- filter_fdr(kegg_pd_tab)
kegg_oscc_sig <- filter_fdr(kegg_oscc_tab)

# Compare pathway IDs rather than descriptions (wording can change in KEGG).
shared_ids <- sort(intersect(as.character(kegg_pd_sig$ID),
                             as.character(kegg_oscc_sig$ID)))
shared_paths <- data.frame(
  ID = shared_ids,
  Description_PD = as.character(kegg_pd_sig$Description[
    match(shared_ids, kegg_pd_sig$ID)]),
  Description_OSCC = as.character(kegg_oscc_sig$Description[
    match(shared_ids, kegg_oscc_sig$ID)]),
  PD_FDR = kegg_pd_sig$p.adjust[match(shared_ids, kegg_pd_sig$ID)],
  OSCC_FDR = kegg_oscc_sig$p.adjust[match(shared_ids, kegg_oscc_sig$ID)],
  stringsAsFactors = FALSE
)
utils::write.csv(shared_paths,
                 file.path(out_dir, "04_shared_significant_KEGG.csv"),
                 row.names = FALSE)

# These 3 shared KEGG pathways are explicitly named in the manuscript.
reported_shared <- data.frame(
  KEGG_ID = c("hsa04060", "hsa04061", "hsa05323"),
  Manuscript_pathway = c(
    "Cytokine-cytokine receptor interaction",
    "Viral protein interaction with cytokine and cytokine receptor",
    "Rheumatoid arthritis"
  ),
  stringsAsFactors = FALSE
)
reported_shared$PD_significant <- reported_shared$KEGG_ID %in% kegg_pd_sig$ID
reported_shared$OSCC_significant <- reported_shared$KEGG_ID %in% kegg_oscc_sig$ID
reported_shared$Shared_significant <- reported_shared$PD_significant &
                                      reported_shared$OSCC_significant
utils::write.csv(reported_shared,
                 file.path(out_dir, "04_manuscript_KEGG_pathway_QC.csv"),
                 row.names = FALSE)

# --- 5. Visualization (historical two figures are GO-only bubble plots) ---
# Original script named these GO-only plots GO_KEGG_PD.png/GO_KEGG_OSCC.png.
# Preserve these filenames while documenting that KEGG isn't plotted here.
make_bubble <- function(tab, title, output, low, high, n = 10) {
  tab <- filter_fdr(tab)
  if (!nrow(tab)) {
    warning("No significant terms for bubble plot: ", title)
    return(invisible(FALSE))
  }
  tab <- tab[order(tab$p.adjust, tab$ID), , drop = FALSE]
  tab <- utils::head(tab, n)
  tab$LogP <- -log10(pmax(tab$p.adjust, .Machine$double.xmin))
  tab$Description <- factor(tab$Description,
                            levels = rev(unique(tab$Description)))
  p <- ggplot2::ggplot(tab,
       ggplot2::aes(x = LogP, y = Description, size = Count, color = LogP)) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::scale_color_gradient(low = low, high = high,
                                  name = "-log10\n(adj. p)") +
    ggplot2::scale_size_continuous(range = c(4, 14), name = "Gene\ncount") +
    ggplot2::labs(title = title, x = "-log10 (adj. p-value)", y = NULL) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(
      hjust = 0.5, face = "bold", size = 15))
  grDevices::png(file.path(out_dir, output),
                 width = 3600, height = 2400, res = 300)
  tryCatch(print(p), finally = grDevices::dev.off())
  invisible(TRUE)
}

make_bubble(go_pd_tab, "GO Enrichment — Periodontitis DEGs (GSE10334)",
            "GO_KEGG_PD.png", "#A8DADC", "#1D3563")
make_bubble(go_oscc_tab, "GO Enrichment — OSCC DEGs (GSE30784)",
            "GO_KEGG_OSCC.png", "#F4A261", "#E63946")

# Supplemental KEGG pathway visualizations; NOT claimed to be published plots.
make_bubble(kegg_pd_tab, "KEGG Enrichment — Periodontitis DEGs",
            "04_KEGG_PD_diagnostic.png", "#A8DADC", "#1D3563")
make_bubble(kegg_oscc_tab, "KEGG Enrichment — OSCC DEGs",
            "04_KEGG_OSCC_diagnostic.png", "#F4A261", "#E63946")

# --- 6. Audit summary, settings and session --------------------------------
summary_04 <- data.frame(
  Dataset = c("PD", "OSCC"),
  DEG_entries = c(length(pd_symbols), length(oscc_symbols)),
  Mapped_unique_symbols = mapping_summary$Mapped_unique_symbols,
  Unmapped_pct = mapping_summary$Unmapped_pct,
  Distinct_Entrez = mapping_summary$Unique_Entrez_IDs,
  GO_BP_significant = c(nrow(go_pd_sig), nrow(go_oscc_sig)),
  KEGG_significant = c(nrow(kegg_pd_sig), nrow(kegg_oscc_sig)),
  stringsAsFactors = FALSE
)
utils::write.csv(summary_04, file.path(out_dir, "04_enrichment_summary.csv"),
                 row.names = FALSE)

checks <- data.frame(
  Metric = c("PD_DEG_entries", "OSCC_DEG_entries",
             "PD_GO_BP_significant", "OSCC_GO_BP_significant"),
  Expected_manuscript = c(DEG_ENTRIES_PUBLISHED["PD"],
                          DEG_ENTRIES_PUBLISHED["OSCC"],
                          GO_TERMS_PUBLISHED["PD"],
                          GO_TERMS_PUBLISHED["OSCC"]),
  Observed = c(length(pd_symbols), length(oscc_symbols),
               nrow(go_pd_sig), nrow(go_oscc_sig)),
  stringsAsFactors = FALSE
)
checks$Match <- checks$Expected_manuscript == checks$Observed
utils::write.csv(checks,
                 file.path(out_dir, "04_manuscript_count_QC.csv"),
                 row.names = FALSE)

settings <- c(
  paste0("Date: ", Sys.Date()),
  "Study: PD-OSCC transcriptomic analysis, Script 04 historical DEG lists",
  "Input: results/02_DEG/02_*_unique_DEG_entries.csv",
  "DEG model: historical_unblocked (PD 57, OSCC 1266)",
  "GO: BP; KEGG: hsa; SYMBOL->ENTREZID via clusterProfiler::bitr",
  "Adjustment: Benjamini-Hochberg; report adjusted p < 0.05",
  "pvalueCutoff=0.05, qvalueCutoff=0.2, minGSSize=10, maxGSSize=500",
  "No explicit universe: default annotation-dependent background",
  "KEGG: live online data; pathway versions subject to change",
  "Historical GO reported: PD 16; OSCC 788",
  "Historical unmapped proportions reported: PD 7.0%; OSCC 5.9%",
  "Original GO plotting filenames are historical names, despite GO-only plots"
)
writeLines(settings, file.path(out_dir, "04_analysis_settings.txt"))
writeLines(capture.output(utils::sessionInfo()),
           file.path(out_dir, "04_sessionInfo.txt"))

cat("\n========== SCRIPT 04 SUMMARY ==========\n")
print(summary_04, row.names = FALSE)
cat("\n========== MANUSCRIPT COUNT QC ==========\n")
print(checks, row.names = FALSE)
cat("\n========== REPORTED SHARED KEGG PATHWAY QC ==========\n")
print(reported_shared, row.names = FALSE)
cat("\nAll shared significant KEGG pathways (ID-based):\n")
print(shared_paths, row.names = FALSE)

if (!all(checks$Match) || !all(reported_shared$Shared_significant)) {
  warning("SCRIPT 04 REQUIRES REVIEW: GO counts or published KEGG pathways ",
          "differ from current enrichment results. Do NOT change thresholds ",
          "just to match the article. Check gene mapping, annotation/KEGG ",
          "versions and manuscript methods. See saved audit CSVs.")
} else {
  cat("MANUSCRIPT COUNT & SHARED KEGG CHECKS MATCH.",
      " Gene mapping and top-ranked biological terms still require review.\n")
}
cat("STEP 04 COMPLETE. Outputs: ", out_dir, "\n", sep = "")
