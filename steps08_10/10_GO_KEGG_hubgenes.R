# Step 10: exploratory over-representation analysis of the historical five-gene
# panel. This is NOT the DEG enrichment performed in Step 04. The manuscript's
# 16/788 historical GO counts refer to the DEG lists, not these five genes.
#
# Run mode is the default. To check inputs, set PD_OSCC_10_MODE=check_inputs before
# sourcing. GO uses the installed org.Hs.eg.db; KEGG uses the online database.
# Database and package updates can change terms and adjusted p-values. The
# archived KEGG results are references, not thresholds or tuning targets.
local({
 main <- function() {
  check <- function(x, ...) if(!isTRUE(x)) stop(..., call. = FALSE)
  mode <- Sys.getenv("PD_OSCC_10_MODE", "run")
  check(mode %in% c("check_inputs", "run"), "PD_OSCC_10_MODE must be check_inputs or run")
  paths <- unlist(lapply(sys.frames(), function(frame) {
    v <- frame$ofile
    if (is.character(v) && length(v) == 1L && nzchar(v)) v else NULL
  }), use.names = FALSE)
  matches <- paths[basename(paths) == "10_GO_KEGG_hubgenes.R"]
  source_dir <- if(length(matches)) dirname(tail(matches, 1L)) else ""
  candidates <- unique(c(source_dir, Sys.getenv("PD_OSCC_10_PACKAGE_DIR", ""), getwd()))
  candidates <- candidates[nzchar(candidates)]
  req <- c("10_assets/10_archived_five_gene_panel.csv",
           "10_assets/10_archived_SYMBOL_ENTREZ_map.csv",
           "10_assets/10_archived_five_gene_KEGG_result.csv")
  valid <- vapply(candidates, function(x) all(file.exists(file.path(x, req))), logical(1))
  check(any(valid), "Cannot locate the 10_assets folder alongside Script 10.")
  root <- normalizePath(candidates[which(valid)[1L]], winslash = "/", mustWork = TRUE)
  files <- file.path(root, req); names(files) <- c("panel", "archived_mapping", "archived_KEGG")
  panel_table <- utils::read.csv(files[["panel"]], stringsAsFactors = FALSE)
  ref_map <- utils::read.csv(files[["archived_mapping"]], colClasses = "character")
  ref <- utils::read.csv(files[["archived_KEGG"]], stringsAsFactors = FALSE)
  check(identical(names(panel_table), "Gene") && nrow(panel_table) == 5L &&
          !anyNA(panel_table$Gene) && !anyDuplicated(panel_table$Gene) &&
          setequal(panel_table$Gene, c("CXCL1", "CXCL6", "CXCR4", "RGS1", "CHST2")),
        "Unexpected archived panel")
  check(all(c("SYMBOL", "ENTREZID") %in% names(ref_map)) && nrow(ref_map) == 5L &&
          all(c("ID", "Description", "p.adjust", "geneID") %in% names(ref)) &&
          !anyDuplicated(ref$ID), "Unexpected archived evidence files")
  pkgs <- c("clusterProfiler", "org.Hs.eg.db")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  check(!length(missing), "Missing required R/Bioconductor packages: ", paste(missing, collapse = ", "))
  cat("Archived five-gene panel:", paste(panel_table$Gene, collapse = ", "), "\n")
  cat("Archived KEGG object:", nrow(ref), "stored rows;",
      sum(is.finite(ref$p.adjust) & ref$p.adjust < 0.05), "with adjusted p < 0.05.\n")
  cat("Archived GO result is not available in the checked workspace; no GO count is assumed.\n")
  if(mode == "check_inputs") {
    cat("INPUT CHECK COMPLETE. No enrichment calculations or output files created.\n")
    return(invisible(NULL))
  }

  base <- Sys.getenv("PD_OSCC_10_OUTPUT_BASE", file.path(root, "results/10_five_gene_enrichment"))
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(base, paste0(format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"), "Z_pid", Sys.getpid()))
  check(!dir.exists(out) && dir.create(out), "Cannot create new result directory: ", out)
  status <- file.path(out, "10_run_status.txt"); writeLines("STARTED", status)
  on.exit({if(file.exists(status) && identical(readLines(status, warn=FALSE)[1L], "STARTED"))
    writeLines("INCOMPLETE: inspect console and saved partial outputs", status)}, add = TRUE)
  write_csv <- function(x, name) utils::write.csv(x, file.path(out, name), row.names = FALSE, na = "")
  manifest <- data.frame(Input = req, MD5 = unname(tools::md5sum(files)))
  write_csv(manifest, "10_input_manifest.csv")
  writeLines(c(paste("UTC start:", format(Sys.time(), tz="UTC", usetz=TRUE)),
               "Five-gene exploratory analysis; not the disease-specific DEG enrichment.",
               "GO BP: enrichGO, BH, pvalueCutoff 0.05, qvalueCutoff 0.2, readable TRUE.",
               "KEGG: enrichKEGG, organism hsa, BH, pvalueCutoff 0.05; package defaults otherwise.",
               "No explicit universe was supplied: database/package defaults apply.",
               "No GO archive reference identified; KEGG reference is a saved object, not re-analysis.",
               "Online KEGG is time-dependent; inspect differences without tuning to historical results.",
               capture.output(utils::sessionInfo())), file.path(out, "10_methods_sessionInfo.txt"))

  panel <- as.character(panel_table$Gene)
  mapped <- suppressMessages(clusterProfiler::bitr(panel, fromType = "SYMBOL", toType = "ENTREZID",
                                                    OrgDb = org.Hs.eg.db::org.Hs.eg.db))
  check(is.data.frame(mapped) && all(c("SYMBOL", "ENTREZID") %in% names(mapped)), "Mapping failed")
  mapped$SYMBOL <- as.character(mapped$SYMBOL); mapped$ENTREZID <- as.character(mapped$ENTREZID)
  write_csv(mapped, "10_current_symbol_entrez_mapping.csv")
  check(!anyNA(mapped$SYMBOL) && !anyNA(mapped$ENTREZID) &&
          all(panel %in% mapped$SYMBOL) && !anyDuplicated(mapped$SYMBOL) &&
          !anyDuplicated(mapped$ENTREZID),
        "The current database does not map the five genes one-to-one; inspect saved mapping before proceeding")
  map_comparison <- merge(ref_map, mapped, by = "SYMBOL", all = TRUE, suffixes = c("_archived", "_current"))
  map_comparison$identical_mapping <- with(map_comparison, !is.na(ENTREZID_archived) &
                                            !is.na(ENTREZID_current) & ENTREZID_archived == ENTREZID_current)
  write_csv(map_comparison, "10_mapping_comparison.csv")
  ids <- mapped$ENTREZID
  go <- clusterProfiler::enrichGO(gene = ids, OrgDb = org.Hs.eg.db::org.Hs.eg.db, ont = "BP",
    pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
  saveRDS(go, file.path(out, "10_five_gene_GO_BP.rds"))
  if(is.null(go)) {
    cat("No GO result returned.\n"); go_tab <- data.frame()
  } else go_tab <- as.data.frame(go)
  write_csv(go_tab, "10_five_gene_GO_BP.csv")
  cat("Current GO BP rows:", nrow(go_tab), "(no historical GO reference recovered).\n")

  # Online KEGG call: a failed request stops the run, retaining the GO and QC
  # outputs. It never fabricates a substitute or copies the historical table.
  kegg <- tryCatch(clusterProfiler::enrichKEGG(gene = ids, organism = "hsa",
    pAdjustMethod = "BH", pvalueCutoff = 0.05), error = function(e) {
      writeLines(c("KEGG request failed; GO and mapping outputs retained.", conditionMessage(e)),
                 file.path(out, "10_KEGG_error.txt"))
      stop("KEGG request failed; inspect 10_KEGG_error.txt: ", conditionMessage(e), call. = FALSE)
    })
  saveRDS(kegg, file.path(out, "10_five_gene_KEGG.rds"))
  kegg_tab <- if(is.null(kegg)) data.frame() else as.data.frame(kegg)
  write_csv(kegg_tab, "10_five_gene_KEGG.csv")
  if(nrow(kegg_tab)) {
    check(all(c("ID", "p.adjust") %in% names(kegg_tab)) && !anyDuplicated(kegg_tab$ID),
          "Unexpected live KEGG result schema")
    current_sig <- kegg_tab[is.finite(kegg_tab$p.adjust) & kegg_tab$p.adjust < 0.05, , drop = FALSE]
  } else current_sig <- data.frame(ID = character(), p.adjust = numeric())
  historic_sig <- ref[is.finite(ref$p.adjust) & ref$p.adjust < 0.05, , drop = FALSE]
  ids_union <- union(as.character(historic_sig$ID), as.character(current_sig$ID))
  comparison <- data.frame(ID = ids_union,
    Archived_significant = ids_union %in% historic_sig$ID,
    Current_significant = ids_union %in% current_sig$ID,
    Archived_adjusted_p = historic_sig$p.adjust[match(ids_union, historic_sig$ID)],
    Current_adjusted_p = current_sig$p.adjust[match(ids_union, current_sig$ID)])
  write_csv(comparison, "10_KEGG_significant_ID_comparison.csv")
  cat("Historical KEGG significant pathways:", nrow(historic_sig), "\n")
  cat("Current KEGG significant pathways:", nrow(current_sig), "\n")
  cat("KEGG calls and GO annotations depend on database versions; compare results before drawing conclusions.\n")
  writeLines(c("COMPLETED", "Historical and current KEGG results must be compared; GO has no historical reference.",
               "This analysis uses five selected genes and is separate from Step 04 DEG enrichment."), status)
  cat("Results saved to:", out, "\n")
  invisible(out)
 }
 main()
})
