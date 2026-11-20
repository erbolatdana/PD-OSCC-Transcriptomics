# Prepare the small, analysis-ready inputs for Step 09 from archived Step 01,
# CIBERSORTx Job29/Job30, and reanalysis 01D files. This is provenance code;
# the already-prepared inputs are included in the Step 09 package.
#
# Configure the three environment variables below before running in a clean
# directory. This script reads existing files; it does not run GEOquery,
# CIBERSORTx, WGCNA, DEG or any statistical model. It saves new RDS files only.
local({
  main <- function() {
    check <- function(x, ...) if (!isTRUE(x)) stop(..., call. = FALSE)
    check(requireNamespace("Biobase", quietly = TRUE), "Biobase is required")
    required_env <- c("PD_OSCC_PROJECT_ROOT", "PD_OSCC_09_REANALYSIS", "PD_OSCC_09_PREP_OUT")
    env <- Sys.getenv(required_env, unset = ""); names(env) <- required_env
    check(all(nzchar(env)), "Set environment variables: ", paste(required_env, collapse = ", "))
    root <- normalizePath(env[["PD_OSCC_PROJECT_ROOT"]], winslash = "/", mustWork = TRUE)
    reanalysis <- normalizePath(env[["PD_OSCC_09_REANALYSIS"]], winslash = "/", mustWork = TRUE)
    output <- env[["PD_OSCC_09_PREP_OUT"]]
    check(!dir.exists(output) && dir.create(output, recursive = TRUE),
          "Output directory exists or could not be created; choose a new directory.")
    files <- c(
      step01 = file.path(root, "results/01_preprocessing/01_preprocessed_esets.rds"),
      PD_map = file.path(root, "results/01_preprocessing/01_PD_sample_map.csv"),
      OSCC_map = file.path(root, "results/01_preprocessing/01_OSCC_sample_map.csv"),
      immune_QC = file.path(root, "results/03_CIBERSORT/03_cell_zero_QC.csv"),
      PD_job = file.path(root, "data/cibersortx_historical/CIBERSORTx_Job29_Results.csv"),
      OSCC_job = file.path(root, "data/cibersortx_historical/CIBERSORTx_Job30_Results.csv"),
      sensitivity_PD = file.path(reanalysis, "data/processed/GSE10334_01D_log2_gene_matrix.rds"),
      sensitivity_OSCC = file.path(reanalysis, "data/processed/GSE30784_01D_log2_gene_matrix.rds"),
      sensitivity_PD_map = file.path(reanalysis, "data/processed/GSE10334_01D_sample_map.csv"),
      sensitivity_OSCC_map = file.path(reanalysis, "data/processed/GSE30784_01D_sample_map.csv")
    )
    check(all(file.exists(files)), "Missing inputs: ", paste(files[!file.exists(files)], collapse = ", "))
    md5 <- unname(tools::md5sum(files)); check(!anyNA(md5), "Cannot fingerprint inputs")
    manifest <- data.frame(Role = names(files), File_basename = basename(files), MD5 = md5)
    utils::write.csv(manifest, file.path(output, "09_source_input_manifest.csv"), row.names = FALSE)
    genes <- c("CXCL1", "CXCL6", "CXCR4", "RGS1", "CHST2")
    qc <- utils::read.csv(files[["immune_QC"]], stringsAsFactors = FALSE)
    check(all(c("Cell", "PD_pass", "OSCC_pass") %in% names(qc)) &&
            !anyDuplicated(qc$Cell), "Invalid immune QC input")
    qc$PD_pass <- as.logical(qc$PD_pass); qc$OSCC_pass <- as.logical(qc$OSCC_pass)
    check(!anyNA(qc$PD_pass) && !anyNA(qc$OSCC_pass), "Missing immune QC flags")
    pre <- readRDS(files[["step01"]])
    first_symbol <- function(x) vapply(strsplit(as.character(x), "///", fixed = TRUE),
      function(z) if (!length(z) || is.na(z[1L])) NA_character_ else trimws(z[1L]), character(1))
    for (disease in c("PD", "OSCC")) {
      suffix <- tolower(disease); is_pd <- disease == "PD"; n <- if (is_pd) 247L else 212L
      eset <- pre[[paste0("gset_", suffix)]]
      x <- Biobase::exprs(eset); ann <- Biobase::fData(eset)
      meta <- utils::read.csv(files[[paste0(disease, "_map")]], stringsAsFactors = FALSE)
      job <- utils::read.csv(files[[paste0(disease, "_job")]], check.names = FALSE,
                             stringsAsFactors = FALSE)
      sens <- readRDS(files[[paste0("sensitivity_", disease)]])
      sens_meta <- utils::read.csv(files[[paste0("sensitivity_", disease, "_map")]], stringsAsFactors = FALSE)
      check(ncol(x) == n && identical(rownames(x), rownames(ann)) &&
              identical(as.character(meta$GSM), colnames(x)) &&
              identical(as.character(sens_meta$GSM), colnames(sens)) &&
              identical(colnames(x), colnames(sens)) &&
              identical(as.character(meta$Group), as.character(sens_meta$Group)) &&
              setequal(as.character(job$Mixture), colnames(x)) && !anyDuplicated(job$Mixture),
            disease, ": source samples/groups do not match")
      if (is_pd) {
        check("Patient_ID" %in% names(meta) && "Patient_ID" %in% names(sens_meta) &&
                !anyNA(meta$Patient_ID) && all(nzchar(as.character(meta$Patient_ID))) &&
                identical(as.character(meta$Patient_ID), as.character(sens_meta$Patient_ID)) &&
                length(unique(meta$Patient_ID)) == 90L, "PD patient groups invalid")
        # Preserve lexicographic group ordering to retain the historical seeded
        # cluster-bootstrap sequence, while not exporting original patient IDs.
        original_ids <- as.character(meta$Patient_ID)
        unit <- sprintf("P%03d", match(original_ids, sort(unique(original_ids))))
      } else unit <- as.character(meta$GSM)
      cells <- qc$Cell[if (is_pd) qc$PD_pass else qc$OSCC_pass]
      check(length(cells) == if (is_pd) 16L else 14L, "Unexpected cell counts")
      check(all(cells %in% names(job)) && all(genes %in% rownames(sens)), "Missing genes or cells")
      job <- job[match(colnames(x), job$Mixture), , drop = FALSE]
      immune <- as.matrix(job[, cells, drop = FALSE]); storage.mode(immune) <- "double"
      sym <- first_symbol(ann$Gene.symbol)
      chosen <- vapply(genes, function(g) {
        ix <- which(!is.na(sym) & sym == g); check(length(ix) > 0, "Missing probe: ", g)
        ix[which.max(rowMeans(x[ix, , drop = FALSE]))]
      }, integer(1))
      historical_expr <- t(x[chosen, , drop = FALSE]); colnames(historical_expr) <- genes
      sensitivity_expr <- t(sens[genes, , drop = FALSE])
      check(all(is.finite(historical_expr)) && all(is.finite(sensitivity_expr)) &&
              all(is.finite(immune)), "Non-finite input values")
      data <- list(historical_expr = historical_expr, sensitivity_expr = sensitivity_expr,
                   immune = immune, unit = unit,
                   map = data.frame(GSM = as.character(meta$GSM), Group = as.character(meta$Group),
                                    Analysis_unit = unit, stringsAsFactors = FALSE),
                   probe_map = data.frame(Disease = disease, Gene = genes,
                                          GPL570_Probe_ID = rownames(x)[chosen]))
      saveRDS(data, file.path(output, paste0("09_", disease, "_analysis_inputs.rds")), version = 2)
      cat(disease, ": saved ", n, " samples; anonymized PD group IDs where applicable.\n", sep = "")
    }
  }
  main()
})
