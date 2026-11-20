# Script 07: reproduce the SAVED fixed-panel corrected analyses (07A + 07D).
# The source scripts in ./source use the input root supplied by this wrapper.
# Their output directory assignments are redirected for this separate run.
# No code in this wrapper derives a panel or performs feature selection.
# Default mode is read-only input verification. Run mode is explicitly enabled
# with PD_OSCC_07_MODE=run; results go to a new, timestamped project folder.
# This fixed-panel analysis is not nested CV; prior selection used full data.
# GSE16134 overlaps the PD training cohort and is not independently validated.
local({
  main <- function() {
    fail <- function(...) stop(..., call. = FALSE)
    check <- function(x, ...) if (!isTRUE(x)) fail(...)
    mode <- Sys.getenv("PD_OSCC_07_MODE", "check_inputs")
    check(mode %in% c("check_inputs", "run"), "Unknown PD_OSCC_07_MODE: ", mode)
    project <- normalizePath(Sys.getenv("PD_OSCC_PROJECT_ROOT", getwd()),
                             winslash = "/", mustWork = TRUE)
    # Find this script when sourced in R/RStudio, or executed with Rscript.
    # Environment overrides remain available for nonstandard launch methods.
    supplied_package <- Sys.getenv("PD_OSCC_07_PACKAGE_DIR", "")
    if (nzchar(supplied_package)) {
      package_dir <- supplied_package
    } else {
      source_paths <- vapply(sys.frames(), function(frame) {
        path <- frame$ofile
        if (is.null(path) || !length(path)) "" else as.character(path[[1L]])
      }, character(1))
      source_paths <- source_paths[nzchar(source_paths)]
      script_path <- if (length(source_paths)) tail(source_paths, 1L) else ""
      if (!nzchar(script_path)) {
        file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
        if (length(file_arg)) script_path <- sub("^--file=", "", file_arg[[1L]])
      }
      package_dir <- if (nzchar(script_path)) dirname(script_path) else getwd()
    }
    package_dir <- normalizePath(package_dir, winslash = "/", mustWork = TRUE)
    old_root <- Sys.getenv("PD_OSCC_07_SOURCE_DATA", file.path(package_dir, "data_root"))
    check(dir.exists(old_root), "Input data directory not found: ", old_root,
          ". See README_RU.md for the required files and layout.")
    old_root <- normalizePath(old_root, winslash = "/", mustWork = TRUE)
    source_a <- file.path(package_dir, "source/07A_Major4_grouped_CV_fixed_panel_v2.R")
    source_d <- file.path(package_dir, "source/07D_Major4_OSCC_external_fixed_panel_validation.R")
    files <- c(
      PD_matrix = file.path(old_root, "data/processed/GSE10334_01D_log2_gene_matrix.rds"),
      PD_map = file.path(old_root, "data/processed/GSE10334_01D_sample_map.csv"),
      OSCC_matrix = file.path(old_root, "data/processed/GSE30784_01D_log2_gene_matrix.rds"),
      OSCC_map = file.path(old_root, "data/processed/GSE30784_01D_sample_map.csv"),
      external = file.path(old_root, "data/external_validation/original_GEO/GSE37991_ORIGINAL_GEO_ExpressionSet.rds"),
      external_map = file.path(old_root, "data/external_validation/original_GEO/GSE37991_VERIFIED_sample_map.csv"),
      external_annotation = file.path(old_root, "results/07_Major4_external_validation/02_mapping_overlap_QC/GSE37991_five_gene_probe_annotation.csv"),
      source_07A = source_a, source_07D = source_d
    )
    check(all(file.exists(files)), "Missing input(s): ", paste(files[!file.exists(files)], collapse = ", "))
    packages <- c("Biobase", "pROC", "caret")
    missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
    check(!length(missing), "Required packages unavailable: ", paste(missing, collapse = ", "))
    md5 <- unname(tools::md5sum(files))
    check(!anyNA(md5), "Could not calculate an input hash.")
    # The five external-analysis input hashes come from the saved 07D manifest.
    # The two PD fingerprints record currently saved files; the previous 07A
    # output did not include an input-hash manifest, so they are not proof of
    # which PD files were used in that earlier run.
    expected <- c(
      PD_matrix = "8f2260e8a9a973d2d6a1f49a68cc4895",
      PD_map = "a676088e8c3112c68a7be5a81a3bf684",
      OSCC_matrix = "c21f165e4a84c3d45aed18b24d2488b9",
      OSCC_map = "bbd12655e680657eacfcb9953784d02b",
      external = "95d25fcdb4fdec48f25a3ce14f604eea",
      external_map = "199002fc348d6652e69944c3bae80a42",
      external_annotation = "6f506ee7e718e4b5e75eb699719fa28a",
      source_07A = "50b5c8ef4d3d28023cf235c8c4119d7b",
      source_07D = "688404994365415821ded4e995235d37"
    )
    mismatch <- names(files)[tolower(md5) != expected[names(files)]]
    check(!length(mismatch), "Input/source fingerprint differs; investigate before running: ",
          paste(mismatch, collapse = ", "))
    panel <- c("CXCL1", "CXCL6", "CXCR4", "RGS1", "CHST2")
    for (disease in c("PD", "OSCC")) {
      x <- readRDS(files[[paste0(disease, "_matrix")]])
      m <- utils::read.csv(files[[paste0(disease, "_map")]], stringsAsFactors = FALSE)
      expected_n <- if (disease == "PD") 247L else 212L
      check(is.matrix(x) && is.numeric(x) && ncol(x) == expected_n &&
              !anyDuplicated(rownames(x)) && !anyDuplicated(colnames(x)) &&
              all(panel %in% rownames(x)) && all(is.finite(x[panel, , drop = FALSE])),
            disease, ": invalid expression matrix/panel.")
      check(all(c("GSM", "Group") %in% names(m)) &&
              !anyNA(m$GSM) && !anyDuplicated(m$GSM) &&
              identical(as.character(m$GSM), colnames(x)) && !anyNA(m$Group) &&
              setequal(m$Group, c(disease, "Control")), disease, ": metadata mismatch.")
      if (disease == "PD") {
        check("Patient_ID" %in% names(m) && !anyNA(m$Patient_ID) &&
                all(nzchar(as.character(m$Patient_ID))) && length(unique(m$Patient_ID)) == 90L,
              "PD patient IDs are invalid.")
      }
      cat(disease, ": ", ncol(x), " samples, ", length(panel),
          " fixed genes; GSM order verified.\n", sep = "")
    }
    ext_map <- utils::read.csv(files[["external_map"]], stringsAsFactors = FALSE)
    check(all(c("GSM", "Group", "Patient_ID") %in% names(ext_map)) &&
            nrow(ext_map) == 80L && !anyNA(ext_map$Patient_ID) &&
            !anyDuplicated(ext_map$GSM) && length(unique(ext_map$Patient_ID)) == 40L,
          "External sample map is invalid.")
    pairs <- table(ext_map$Patient_ID, ext_map$Group)
    check(setequal(colnames(pairs), c("Control", "OSCC")) && all(pairs == 1L),
          "External pair structure is invalid.")
    cat("GSE37991: 80 samples, 40 patient pairs.\n")
    cat("Seven input and two source hashes match the expected fingerprints.\n")
    cat("PD input hashes document today's files; earlier 07A inputs have no saved hash manifest.\n")
    if (mode == "check_inputs") {
      cat("INPUT CHECK COMPLETE. No models fitted; no files written.\n")
      return(invisible(NULL))
    }

    base <- file.path(project, "results/07_machine_learning/reproduction")
    dir.create(base, showWarnings = FALSE, recursive = TRUE)
    run_id <- paste0(format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"), "Z_pid", Sys.getpid())
    run <- file.path(base, run_id)
    check(!dir.exists(run) && dir.create(run), "Cannot create new output directory: ", run)
    cv_out <- file.path(run, "01_fixed_panel_CV")
    ext_out <- file.path(run, "02_external_OSCC")
    check(dir.create(cv_out) && dir.create(ext_out), "Cannot create result subdirectories.")
    status <- file.path(run, "07_run_status.txt")
    writeLines("STARTED: completion not established", status)
    on.exit({
      if (file.exists(status) && identical(readLines(status, warn = FALSE)[1L],
                                             "STARTED: completion not established"))
        writeLines("INCOMPLETE: inspect console and partial outputs", status)
    }, add = TRUE)
    utils::write.csv(data.frame(Role = names(files), File = unname(files), MD5 = md5),
                     file.path(run, "07_input_and_source_manifest.csv"), row.names = FALSE)
    writeLines(capture.output(utils::sessionInfo()), file.path(run, "07_sessionInfo.txt"))

    # Patch just output locations in memory, leaving algorithms/seeds
    # identical to the two archived scripts. No source file is changed at runtime.
    read_redirected <- function(file, replacements) {
      lines <- readLines(file, warn = FALSE, encoding = "UTF-8")
      for (entry in replacements) {
        hits <- which(lines == entry[[1L]])
        check(length(hits) == 1L, "Source layout changed; refusing automatic patch: ", entry[[1L]])
        lines[hits] <- entry[[2L]]
      }
      parse(text = lines, keep.source = TRUE)
    }
    parsed_a <- read_redirected(source_a, list(
      c('out_dir <- file.path(project_dir, "results", "07_Major4_fixed_panel_CV")',
        'out_dir <- .script07_cv_out')
    ))
    parsed_d <- read_redirected(source_d, list(
      c('out <- file.path(project_dir, "results", "07_Major4_external_validation", "03_OSCC_fixed_panel_external")',
        'out <- .script07_ext_out')
    ))
    private <- new.env(parent = globalenv())
    private$.script07_data_root <- old_root
    private$.script07_cv_out <- cv_out
    private$.script07_ext_out <- ext_out
    cat("\n=== 07A: FIXED-PANEL CROSS-VALIDATION ===\n")
    eval(parsed_a, envir = private)
    cat("\n=== 07D: EXTERNAL OSCC VALIDATION ===\n")
    eval(parsed_d, envir = private)

    cv <- utils::read.csv(file.path(cv_out, "Major4_fixed_panel_corrected_CV_summary.csv"))
    external <- utils::read.csv(file.path(ext_out, "Major4_OSCC_GSE37991_external_AUC_summary.csv"))
    check(nrow(cv) == 2L && setequal(cv$Dataset, c("PD", "OSCC")) && nrow(external) == 1L,
          "Unexpected generated summaries.")
    # Compare actual numbers; reference numbers are not used to train/tune.
    ref <- c(PD_mean = 0.904568216175359, PD_pooled = 0.899504781420765,
             OSCC_mean = 0.967371323529412, OSCC_pooled = 0.959414504324684,
             OSCC_external = 0.79)
    got <- c(PD_mean = cv$Mean_fold_AUC[cv$Dataset == "PD"],
             PD_pooled = cv$Pooled_OOF_AUC[cv$Dataset == "PD"],
             OSCC_mean = cv$Mean_fold_AUC[cv$Dataset == "OSCC"],
             OSCC_pooled = cv$Pooled_OOF_AUC[cv$Dataset == "OSCC"],
             OSCC_external = external$AUC)
    comparison <- data.frame(Metric = names(ref), Saved_reference = unname(ref),
                             New_result = unname(got), Absolute_difference = abs(got - ref),
                             Match_within_1e_6 = abs(got - ref) <= 1e-6,
                             row.names = NULL)
    utils::write.csv(comparison, file.path(run, "07_result_comparison.csv"), row.names = FALSE)
    writeLines(c("COMPLETED", paste("All AUC checks within tolerance:", all(comparison$Match_within_1e_6)),
                 "Fixed panel was selected on full discovery data; CV is not nested.",
                 "GSE16134 is not counted as independent PD validation."), status)
    cat("\n=== COMPARISON WITH PREVIOUS SAVED RESULTS ===\n")
    print(comparison, row.names = FALSE, digits = 9)
    cat("\nResults saved to:", run, "\n")
    invisible(run)
  }
  main()
})
