# =============================================================================
# 06_WGCNA_PPI_hubgenes.R -- historical recipe reconstruction for journal audit
#
# Run from the repository root after Scripts 01 and 02.
# Default: check saved inputs only; no WGCNA calculation and no file writing.
# To explicitly rebuild later:
#   Sys.setenv(PD_OSCC_06_MODE = "rebuild")
#   source("scripts/06_WGCNA_PPI_hubgenes.R")
#
# This reconstructs HISTORICAL probe-row inputs from archived feature manifests,
# retaining repeated gene labels. The original subset-selection algorithm was
# not recovered: first-symbol annotation alone gives a different set of rows.
# It does not use the 20,848-gene matrices or a top-5,000 selection.
# beta=6/3 are fixed historical settings, not asserted successful estimates.
# Reported results are comparison references only, never selection targets.
# A previously completed reconstruction matched the reported counts; see
# audit_reference. A saved STRING export was subsequently recovered, but the
# original STRING version and exact query settings remain unverified.
# Computational reproduction does not validate the original independence or
# gene-background assumptions. Read README.md before interpreting results.
# =============================================================================

local({
  main <- function() {
    fail <- function(...) stop(..., call. = FALSE)
    check <- function(ok, ...) if (!isTRUE(ok)) fail(...)
    mode <- Sys.getenv("PD_OSCC_06_MODE", unset = "check_inputs")
    check(mode %in% c("check_inputs", "rebuild"),
          "PD_OSCC_06_MODE must be check_inputs or rebuild.")
    root <- normalizePath(Sys.getenv("PD_OSCC_PROJECT_ROOT", unset = getwd()),
                          winslash = "/", mustWork = TRUE)
    reference_dir <- Sys.getenv("PD_OSCC_06_REFERENCE_DIR",
                                 unset = file.path(root, "data/reference"))
    files <- c(
      preprocessing = file.path(root, "results/01_preprocessing/01_preprocessed_esets.rds"),
      PD_map = file.path(root, "results/01_preprocessing/01_PD_sample_map.csv"),
      OSCC_map = file.path(root, "results/01_preprocessing/01_OSCC_sample_map.csv"),
      DEG = file.path(root, "results/02_DEG/02_DEG_objects.rds"),
      DEG_summary = file.path(root, "results/02_DEG/02_DEG_summary.csv")
    )
    references <- c(
      PD_features = file.path(reference_dir, "06_PD_historical_features.csv"),
      OSCC_features = file.path(reference_dir, "06_OSCC_historical_features.csv"),
      PD_samples = file.path(reference_dir, "06_PD_historical_samples.csv"),
      OSCC_samples = file.path(reference_dir, "06_OSCC_historical_samples.csv"),
      fingerprints = file.path(reference_dir, "06_historical_matrix_fingerprints.csv")
    )
    files <- c(files, references)
    check(all(file.exists(files)), "Missing saved input/reference file(s): ",
          paste(files[!file.exists(files)], collapse = ", "))
    packages <- c("Biobase", "digest", if (mode == "rebuild") "WGCNA")
    missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
    check(!length(missing), "Install required packages separately: ", paste(missing, collapse = ", "))

    pre <- readRDS(files[["preprocessing"]])
    deg <- readRDS(files[["DEG"]])
    check(is.list(pre) && all(c("gset_pd", "gset_oscc", "gs_pd", "gs_oscc",
                                "sample_map_pd", "sample_map_oscc") %in% names(pre)),
          "Unexpected preprocessing RDS schema.")
    check(is.list(deg) && all(c("model", "tT_pd", "tT_oscc", "deg_pd_clean",
                                "deg_oscc_clean", "deg_oscc_full") %in% names(deg)),
          "Unexpected DEG RDS schema.")
    check(identical(deg$model, "historical_unblocked"),
          "This recipe requires the historical_unblocked Script 02 output.")
    maps <- list(PD = utils::read.csv(files[["PD_map"]], stringsAsFactors = FALSE),
                 OSCC = utils::read.csv(files[["OSCC_map"]], stringsAsFactors = FALSE))
    fingerprints <- utils::read.csv(references[["fingerprints"]], stringsAsFactors = FALSE)
    check(all(c("Disease", "Rows", "Samples", "Unique_gene_labels", "Numeric_matrix_sha256") %in%
                names(fingerprints)) && !anyDuplicated(fingerprints$Disease) &&
            setequal(fingerprints$Disease, c("PD", "OSCC")),
          "Historical matrix fingerprint table is invalid.")
    summary <- utils::read.csv(files[["DEG_summary"]], stringsAsFactors = FALSE)
    check(all(c("Disease", "Model", "Deduplicated_DEG_entries") %in% names(summary)) &&
            !anyDuplicated(summary$Disease) && setequal(summary$Disease, c("PD", "OSCC")),
          "Unexpected DEG summary schema.")
    check(identical(as.character(summary$Model[match(c("PD", "OSCC"), summary$Disease)]),
                    c("historical_unblocked", "unblocked")), "DEG summary model mismatch.")

    first_symbol <- function(annotation) {
      pieces <- strsplit(as.character(annotation), "///", fixed = TRUE)
      vapply(pieces, function(x) {
        if (!length(x) || is.na(x[1L])) NA_character_ else trimws(x[1L])
      }, character(1))
    }
    # Verify completeness against saved all-probe statistics; do not refit DEG.
    verify_deg <- function(disease) {
      suffix <- tolower(disease)
      tt <- deg[[paste0("tT_", suffix)]]
      dd <- deg[[paste0("deg_", suffix, "_clean")]]
      check(is.data.frame(tt) && all(c("Gene.symbol", "logFC", "adj.P.Val") %in% names(tt)),
            disease, ": invalid all-probe DEG table.")
      check(is.data.frame(dd) && all(c("Gene", "logFC", "adj.P.Val") %in% names(dd)),
            disease, ": invalid DEG list.")
      keep <- is.finite(tt$logFC) & abs(tt$logFC) >= 1.5 &
        is.finite(tt$adj.P.Val) & tt$adj.P.Val < 0.05
      expected <- tt[keep, , drop = FALSE]
      expected$Gene <- first_symbol(expected$Gene.symbol)
      expected <- expected[!is.na(expected$Gene) & expected$Gene != "", , drop = FALSE]
      expected <- expected[order(-abs(expected$logFC), seq_len(nrow(expected))), , drop = FALSE]
      expected <- expected[!duplicated(expected$Gene), , drop = FALSE]
      check(!anyNA(dd$Gene) && all(nzchar(dd$Gene)) && !anyDuplicated(dd$Gene),
            disease, ": missing/duplicate DEG names.")
      check(identical(as.character(dd$Gene), as.character(expected$Gene)) &&
              isTRUE(all.equal(dd$logFC, expected$logFC, tolerance = 0)) &&
              isTRUE(all.equal(dd$adj.P.Val, expected$adj.P.Val, tolerance = 0)),
            disease, ": DEG list differs from saved full table; do not use a top-50 subset.")
      check(nrow(dd) == summary$Deduplicated_DEG_entries[summary$Disease == disease],
            disease, ": DEG count differs from saved summary.")
      dd
    }
    deg_lists <- list(PD = verify_deg("PD"), OSCC = verify_deg("OSCC"))
    check(identical(as.character(deg$deg_oscc_full), as.character(deg_lists$OSCC$Gene)),
          "deg_oscc_full does not match the complete OSCC DEG table.")

    sample_expected <- c(PD = 247L, OSCC = 212L)
    disease_expected <- c(PD = 183L, OSCC = 167L)
    control_expected <- c(PD = 64L, OSCC = 45L)
    input <- list()
    for (disease in c("PD", "OSCC")) {
      suffix <- tolower(disease)
      eset <- pre[[paste0("gset_", suffix)]]
      x <- Biobase::exprs(eset)
      ann <- Biobase::fData(eset)
      m <- maps[[disease]]
      saved_map <- pre[[paste0("sample_map_", suffix)]]
      groups <- as.character(pre[[paste0("gs_", suffix)]])
      check(is.matrix(x) && is.numeric(x) && all(is.finite(x)),
            disease, ": invalid/non-finite probe expression matrix.")
      check(all(c("GSM", "Group") %in% names(m)) && !anyNA(m$GSM) &&
              all(nzchar(m$GSM)) && !anyDuplicated(m$GSM) && !anyNA(m$Group),
            disease, ": invalid sample map.")
      check(identical(colnames(x), as.character(m$GSM)) &&
              identical(as.character(saved_map$GSM), as.character(m$GSM)) &&
              identical(as.character(saved_map$Group), as.character(m$Group)) &&
              identical(groups, as.character(m$Group)) &&
              identical(as.character(eset$group), as.character(m$Group)),
            disease, ": sample order/group mismatch; no silent reordering.")
      check(ncol(x) == sample_expected[[disease]] &&
              setequal(m$Group, c(disease, "Control")) &&
              sum(m$Group == disease) == disease_expected[[disease]] &&
              sum(m$Group == "Control") == control_expected[[disease]],
            disease, ": sample/group counts differ from the audited cohort.")
      if (disease == "PD") {
        check("Patient_ID" %in% names(m) && "Patient_ID" %in% names(saved_map),
              "PD patient IDs are missing.")
        check(!anyNA(m$Patient_ID) && all(nzchar(trimws(as.character(m$Patient_ID)))) &&
                length(unique(m$Patient_ID)) == 90L &&
                identical(as.character(m$Patient_ID), as.character(saved_map$Patient_ID)),
              "PD patient IDs differ from the saved 90-patient mapping.")
      }
      check("Gene.symbol" %in% names(ann) && identical(rownames(x), rownames(ann)) &&
              !anyDuplicated(rownames(x)), disease, ": probe annotation is misaligned.")
      features <- utils::read.csv(references[[paste0(disease, "_features")]],
                                   stringsAsFactors = FALSE)
      reference_samples <- utils::read.csv(references[[paste0(disease, "_samples")]],
                                            stringsAsFactors = FALSE)
      check(all(c("Network_row", "Probe_ID", "Gene") %in% names(features)) &&
              identical(as.integer(features$Network_row), seq_len(nrow(features))) &&
              !anyNA(features$Probe_ID) && !anyDuplicated(features$Probe_ID) &&
              !anyNA(features$Gene) && all(nzchar(features$Gene)),
            disease, ": invalid archived feature manifest.")
      check(all(c("Sample_order", "GSM", "Group") %in% names(reference_samples)) &&
              identical(as.integer(reference_samples$Sample_order), seq_len(nrow(m))) &&
              identical(as.character(reference_samples$GSM), as.character(m$GSM)) &&
              identical(as.character(reference_samples$Group), as.character(m$Group)),
            disease, ": sample order/groups differ from archived inputs.")
      retained <- match(features$Probe_ID, rownames(x))
      check(!anyNA(retained), disease, ": historical probe IDs are missing from current Script 01 data.")
      symbol <- first_symbol(ann$Gene.symbol)
      check(identical(as.character(features$Gene), symbol[retained]),
            disease, ": annotation differs from archived gene labels.")
      features$Original_probe_row <- retained
      historical_matrix <- x[retained, , drop = FALSE]
      # Keep the historical repeated labels. Do not collapse probes or make.names.
      rownames(historical_matrix) <- features$Gene
      fp <- fingerprints[fingerprints$Disease == disease, , drop = FALSE]
      numeric_hash <- digest::digest(serialize(unname(historical_matrix), NULL, version = 2),
                                      algo = "sha256", serialize = FALSE)
      check(nrow(historical_matrix) == fp$Rows && ncol(historical_matrix) == fp$Samples &&
              length(unique(rownames(historical_matrix))) == fp$Unique_gene_labels &&
              identical(numeric_hash, as.character(fp$Numeric_matrix_sha256)),
            disease, ": reconstructed expression differs from the archived matrix. Stop before WGCNA.")
      input[[disease]] <- list(matrix = historical_matrix, features = features, map = m)
    }
    input_summary <- do.call(rbind, lapply(names(input), function(d) {
      x <- input[[d]]$matrix
      data.frame(Disease = d, Samples = ncol(x), Probe_rows_with_gene_labels = nrow(x),
                 Unique_gene_labels = length(unique(rownames(x))),
                 Repeated_gene_label_rows = sum(duplicated(rownames(x))),
                 Complete_DEG_entries = nrow(deg_lists[[d]]), stringsAsFactors = FALSE)
    }))
    cat("HISTORICAL RECIPE: input check\n")
    print(input_summary, row.names = FALSE)
    cat("Historical expression matrices match archived numeric SHA256 fingerprints.\n")
    cat("Historical expression units are probe rows with gene labels, not unique genes.\n")
    cat("Archived manifests fix the historical subset; its original selection algorithm is unresolved.\n")
    cat("PD module-trait inference retains the original biopsy-level method;\n",
        "the dependence of 247 biopsies from 90 patients remains a limitation.\n", sep = "")
    if (mode == "check_inputs") {
      cat("CHECK_INPUTS COMPLETE. No networks/statistical models were calculated; no files written.\n")
      return(invisible(input_summary))
    }

    # Everything below executes only in explicitly selected rebuild mode.
    out_base <- file.path(root, "results/06_WGCNA/historical_reconstruction")
    dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
    run_id <- paste0(format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"), "Z_pid", Sys.getpid())
    out <- file.path(out_base, run_id)
    check(!dir.exists(out), "Output run directory already exists; refusing to overwrite.")
    check(dir.create(out, recursive = FALSE), "Cannot create output directory: ", out)
    write_csv <- function(x, name) utils::write.csv(x, file.path(out, name), row.names = FALSE)
    log_connection <- file(file.path(out, "06_execution.log"), open = "wt", encoding = "UTF-8")
    sink(log_connection, split = TRUE)
    on.exit({sink(); close(log_connection)}, add = TRUE)
    status_file <- file.path(out, "06_run_status.txt")
    writeLines("STARTED -- completion and agreement have not been established", status_file)
    writeLines(capture.output(utils::sessionInfo()), file.path(out, "06_sessionInfo.txt"))
    write_csv(input_summary, "06_input_summary.csv")
    source_paths <- file.path(root, "scripts", c("01_data_preprocessing.R", "02_DEG_analysis.R",
                                                "06_WGCNA_PPI_hubgenes.R"))
    manifest_paths <- c(unname(files), source_paths[file.exists(source_paths)])
    write_csv(data.frame(File = manifest_paths, MD5 = unname(tools::md5sum(manifest_paths))),
              "06_input_and_code_manifest.csv")
    powers <- c(PD = 6L, OSCC = 3L)
    parameters <- list(
      Purpose = "Historical recipe reconstruction; results are not presumed identical",
      Expression_units = "Archived probe-row manifests, first-symbol labels, duplicates retained",
      Feature_selection_provenance = "Actual archived input subset; original selection algorithm unresolved",
      Fixed_historical_powers = powers, Diagnostic_candidate_powers = 1:20,
      Diagnostic_R2_cut = 0.80, Diagnostic_blockSize = 2000L,
      Network = list(networkType = "unsigned", TOMType = "unsigned", corType = "pearson",
                     power = "Disease-specific fixed historical power", minModuleSize = 30L,
                     reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE,
                     pamRespectsDendro = FALSE, saveTOMs = FALSE, maxBlockSize = 10000L,
                     randomSeed = 54321L, nThreads = 1L),
      Seed_provenance = "Explicit current WGCNA default; original effective seed not independently recorded",
      Top3_rule = "Rank all eigengenes, including grey if present, by historical 3-significant-digit p-value; stable ties",
      PD_inference = "Unadjusted biopsy-level Pearson correlation and Student p-value; patient clustering not accounted for",
      WGCNA_version = as.character(utils::packageVersion("WGCNA"))
    )
    dput(parameters, file.path(out, "06_parameters.R"))
    saveRDS(lapply(c("blockwiseModules", "pickSoftThreshold", "goodSamplesGenes"),
                    function(nm) formals(getExportedValue("WGCNA", nm))),
             file.path(out, "06_effective_function_defaults.rds"))
    # Preserve RNG state in an interactive session; blockwiseModules also restores its seed.
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(list = ".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)

    calculate <- function(disease) {
      cat("\nDATASET:", disease, "\n")
      dat <- t(input[[disease]]$matrix)
      features <- input[[disease]]$features
      map <- input[[disease]]$map
      qc <- WGCNA::goodSamplesGenes(dat, verbose = 1)
      write_csv(data.frame(GSM = rownames(dat), Good_sample = qc$goodSamples),
                paste0("06_", disease, "_sample_QC.csv"))
      write_csv(cbind(features, Good_feature = qc$goodGenes),
                paste0("06_", disease, "_feature_QC.csv"))
      check(all(qc$goodSamples), disease,
            ": sample QC would change the historical cohort. Review saved QC before proceeding.")
      # Feature filtering follows the original goodSamplesGenes step and is fully recorded.
      dat <- dat[, qc$goodGenes, drop = FALSE]
      features <- features[qc$goodGenes, , drop = FALSE]
      check(identical(rownames(dat), as.character(map$GSM)) &&
              identical(colnames(dat), as.character(features$Gene)),
            disease, ": data/metadata alignment changed after QC.")
      write_csv(map, paste0("06_", disease, "_sample_map_used.csv"))
      saveRDS(list(datExpr = dat, features = features, map = map),
               file.path(out, paste0("06_", disease, "_network_input.rds")))

      # Re-evaluate the claimed criterion on the historical input, not 20,848 genes.
      # Diagnostic outcome NEVER changes the fixed historical beta.
      sft <- WGCNA::pickSoftThreshold(dat, powerVector = 1:20, RsquaredCut = 0.80,
                                      networkType = "unsigned", blockSize = 2000,
                                      corFnc = WGCNA::cor, corOptions = list(use = "p"), verbose = 1)
      saveRDS(sft, file.path(out, paste0("06_", disease, "_soft_threshold.rds")))
      write_csv(sft$fitIndices, paste0("06_", disease, "_soft_threshold.csv"))
      fixed <- powers[[disease]]
      row <- sft$fitIndices[sft$fitIndices$Power == fixed, , drop = FALSE]
      check(nrow(row) == 1L, "Missing diagnostic row for fixed historical power.")
      eligible <- sft$fitIndices$Power[is.finite(sft$fitIndices$SFT.R.sq) &
                                        sft$fitIndices$SFT.R.sq >= 0.80]
      first_eligible <- if (length(eligible)) min(eligible) else NA_real_
      diagnostic <- data.frame(Disease = disease, Fixed_historical_beta = fixed,
                                Recalculated_R2_at_fixed_beta = row$SFT.R.sq,
                                Recalculated_slope = row$slope,
                                Smallest_tested_beta_R2_ge_080 = first_eligible,
                                WGCNA_powerEstimate = sft$powerEstimate,
                                Fixed_beta_meets_R2_080 = is.finite(row$SFT.R.sq) && row$SFT.R.sq >= 0.80,
                                Fixed_beta_is_smallest = is.finite(first_eligible) && first_eligible == fixed)
      write_csv(diagnostic, paste0("06_", disease, "_power_claim_check.csv"))
      print(diagnostic, row.names = FALSE)

      set.seed(54321L)
      net <- WGCNA::blockwiseModules(
        dat, power = fixed, networkType = "unsigned", TOMType = "unsigned", corType = "pearson",
        minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,
        numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = FALSE,
        maxBlockSize = 10000, randomSeed = 54321, nThreads = 1, verbose = 1
      )
      saveRDS(net, file.path(out, paste0("06_", disease, "_network.rds")))
      check(length(net$colors) == ncol(dat) && all(net$goodSamples) && all(net$goodGenes),
            disease, ": unexpected internal network filtering; inspect saved network/QC.")
      check(identical(names(net$colors), colnames(dat)),
            disease, ": network assignments are not aligned with input feature rows.")
      colors <- WGCNA::labels2colors(net$colors)
      assignment <- cbind(features, Module_label = as.integer(net$colors), Color = colors)
      write_csv(assignment, paste0("06_", disease, "_probe_to_module.csv"))
      write_csv(as.data.frame(table(Color = colors)), paste0("06_", disease, "_module_sizes.csv"))
      check(nrow(net$MEs) == nrow(map) && identical(rownames(net$MEs), rownames(dat)),
            disease, ": eigengene sample order is not verifiable.")
      trait <- as.numeric(map$Group == disease)
      r <- stats::cor(net$MEs, trait, use = "pairwise.complete.obs")
      p <- WGCNA::corPvalueStudent(r, nrow(dat))
      labels <- suppressWarnings(as.integer(sub("^ME", "", colnames(net$MEs))))
      check(!anyNA(labels) && all(is.finite(r)) && all(is.finite(p)),
            disease, ": invalid eigengene labels or correlations.")
      modules <- data.frame(Module = colnames(net$MEs), Color = WGCNA::labels2colors(labels),
                             Cor = round(r[, 1L], 3), Pvalue = signif(p[, 1L], 3),
                             Cor_unrounded = as.numeric(r[, 1L]), Pvalue_unrounded = as.numeric(p[, 1L]),
                             Historical_ME_order = seq_len(ncol(net$MEs)), stringsAsFactors = FALSE)
      # Preserve the historical rank-on-rounded-p rule; no new FDR/grey filter.
      modules <- modules[order(modules$Pvalue, modules$Historical_ME_order), , drop = FALSE]
      top3 <- head(modules, 3L)
      check(nrow(top3) == 3L, disease, ": fewer than three eigengenes available.")
      write_csv(modules, paste0("06_", disease, "_all_module_trait.csv"))
      write_csv(top3, paste0("06_", disease, "_top3_module_trait.csv"))
      write_csv(data.frame(GSM = rownames(net$MEs), net$MEs, check.names = FALSE),
                paste0("06_", disease, "_module_eigengenes.csv"))
      selected <- unique(assignment$Gene[assignment$Color %in% top3$Color])
      overlap <- intersect(as.character(deg_lists[[disease]]$Gene), selected)
      write_csv(data.frame(Gene = selected), paste0("06_", disease, "_selected_module_genes.csv"))
      write_csv(data.frame(Gene = overlap), paste0("06_", disease, "_DEG_module_intersection.csv"))
      write_csv(data.frame(Gene = setdiff(deg_lists[[disease]]$Gene, features$Gene)),
                paste0("06_", disease, "_DEG_absent_from_network.csv"))
      list(non_grey = length(setdiff(unique(colors), "grey")), top3 = top3,
           overlap = overlap, universe = unique(features$Gene), diagnostic = diagnostic)
    }

    tryCatch(withCallingHandlers({
      result <- list(PD = calculate("PD"), OSCC = calculate("OSCC"))
      final <- intersect(result$PD$overlap, result$OSCC$overlap)
      write_csv(data.frame(Gene = final), "06_shared_candidates.csv")
      original_panel <- c("CXCL1", "CXCL6", "CXCR4", "RGS1", "CHST2")
      comparisons <- data.frame(
        Metric = c("PD_non_grey_modules", "OSCC_non_grey_modules", "PD_DEG_module_intersection",
                   "OSCC_DEG_module_intersection", "Shared_candidate_count", "Exact_original_panel"),
        Reported = c("29", "64", "32", "662", "5", "TRUE"),
        Observed = as.character(c(result$PD$non_grey, result$OSCC$non_grey,
                                  length(result$PD$overlap), length(result$OSCC$overlap),
                                  length(final), as.integer(setequal(final, original_panel)))),
        stringsAsFactors = FALSE
      )
      comparisons$Observed[6L] <- as.character(setequal(final, original_panel))
      comparisons$Match <- comparisons$Reported == comparisons$Observed
      write_csv(comparisons, "06_manuscript_result_comparison.csv")
      write_csv(data.frame(Gene = union(original_panel, final),
                           Reported_panel = union(original_panel, final) %in% original_panel,
                           Reconstructed_panel = union(original_panel, final) %in% final),
                "06_panel_membership_comparison.csv")

      # Reproduce the original calculation transparently, with its unit mismatch.
      # N is NOT presented as a justified count of unique genes. No replacement
      # background or new inferential result is silently substituted here.
      N <- nrow(input$PD$matrix)
      a <- length(result$PD$overlap); b <- length(result$OSCC$overlap); k <- length(final)
      check(N >= max(a, b) && N - b >= a - k, "Invalid legacy hypergeometric counts.")
      legacy <- data.frame(
        Analysis = "Legacy arithmetic only; probe-row background with gene-level lists",
        Background_probe_rows = N, PD_unique_genes = a, OSCC_unique_genes = b, Shared_unique_genes = k,
        Expected_legacy = a * b / N,
        Pvalue_legacy = stats::phyper(k - 1, b, N - b, a, lower.tail = FALSE),
        Independently_justified_gene_background = FALSE, stringsAsFactors = FALSE
      )
      write_csv(legacy, "06_legacy_hypergeometric_arithmetic.csv")
      write_csv(data.frame(Gene = intersect(result$PD$universe, result$OSCC$universe)),
                "06_common_unique_gene_labels_for_background_review.csv")
      # PPI is a separately archived manual step. Do not invent current STRING edges.
      writeLines(final, file.path(out, "06_STRING_candidate_input.txt"))
      writeLines(c(
        "STATUS: THIS SCRIPT DOES NOT QUERY OR INDEPENDENTLY REPRODUCE STRING",
        "A saved STRING interaction export and image are in audit_reference/string.",
        "The exported table contains four undirected edges represented in both directions.",
        "The image includes CHST2 as an isolated node; absence from an edge table does not exclude a node.",
        "The original STRING database version/date and exact query settings remain unverified.",
        "A current query is a new dated result and does not establish historical edge scores.",
        "Shared DEG/module candidates are not established network-centrality hubs by intersection alone."
      ), file.path(out, "06_PPI_audit_status.txt"))
      diagnostics <- rbind(result$PD$diagnostic, result$OSCC$diagnostic)
      write_csv(diagnostics, "06_power_claim_comparison.csv")
      saveRDS(list(parameters = parameters, results = result, final_candidates = final,
                    manuscript_comparison = comparisons, legacy_hypergeometric = legacy),
               file.path(out, "06_run_summary.rds"))
      issues <- c(
        "Historical feature manifests reproduce saved input matrices; the original subset-selection algorithm is unresolved.",
        if (!all(comparisons$Match)) "Reported network/intersection/panel results did not all match.",
        if (!all(diagnostics$Fixed_beta_meets_R2_080 & diagnostics$Fixed_beta_is_smallest))
          "The stated smallest-beta/R2 selection rule did not match both fixed historical powers.",
        if (any(vapply(result, function(x) "grey" %in% x$top3$Color, logical(1))))
          "Grey was selected by the original top-three ranking rule.",
        "Original PD module-trait p-values assume independent biopsies; patient clustering is unaccounted for.",
        "Legacy hypergeometric background counts probe rows rather than unique genes.",
        "A STRING export and image were recovered; the historical query version/settings remain unverified."
      )
      writeLines(issues, file.path(out, "06_interpretation_and_open_issues.txt"))
      writeLines(c("COMPUTATION_COMPLETED -- this is not a declaration of methodological validity",
                    paste("All six manuscript result checks match:", all(comparisons$Match)),
                    "See power checks, interpretation issues and PPI audit status."), status_file)
      print(comparisons, row.names = FALSE)
      cat("\nResults saved to:", out, "\n")
    }, warning = function(w) {
      cat("WARNING:", conditionMessage(w), "\n")
      invokeRestart("muffleWarning")
    }), error = function(e) {
      writeLines(c("FAILED -- partial outputs are retained for audit", conditionMessage(e)), status_file)
      cat("ERROR:", conditionMessage(e), "\n")
      stop(e)
    })
    invisible(out)
  }
  main()
})
