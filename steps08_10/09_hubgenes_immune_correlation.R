# Step 09: historical Figure 9 / Fisher-BH and cluster-bootstrap sensitivity.
# Analysis-ready input matrices are supplied under 09_assets/data. They are
# derived from archived inputs, with PD patient IDs pseudonymised while preserving
# patient clusters and their order. See 09_README.md and 09_assets/source.
# The two branches use different expression preprocessing and must remain separate.
# This script uses base R only; it does not install packages or access the internet.
# Run from the repository root: source("09_hubgenes_immune_correlation.R")
# To check inputs only, set PD_OSCC_09_MODE=check_inputs before sourcing.
# Historical Fisher p-values are approximate for Spearman and do not account
# for repeated PD biopsies. Bootstrap CIs are nominal, not FDR-adjusted.
local({
  main <- function() {
    check <- function(ok, ...) if (!isTRUE(ok)) stop(..., call. = FALSE)
    mode <- Sys.getenv("PD_OSCC_09_MODE", "run")
    check(mode %in% c("check_inputs", "run"), "PD_OSCC_09_MODE must be 'run' or 'check_inputs'")
    paths <- unlist(lapply(sys.frames(), function(frame) {
      v <- frame$ofile
      if (is.character(v) && length(v) == 1L && nzchar(v)) v else NULL
    }), use.names = FALSE)
    matches <- paths[basename(paths) == "09_hubgenes_immune_correlation.R"]
    source_dir <- if (length(matches)) dirname(tail(matches, 1L)) else ""
    # A directory left over from the older Step 09 script may be set in the
    # RStudio session. Prefer the folder holding the current script and use
    # another candidate only if its complete set of assets actually exists.
    candidates <- unique(c(source_dir, Sys.getenv("PD_OSCC_09_PACKAGE_DIR", ""), getwd()))
    candidates <- candidates[nzchar(candidates)]
    required_relative <- c("09_assets/data/09_PD_analysis_inputs.rds",
      "09_assets/data/09_OSCC_analysis_inputs.rds",
      "09_assets/reference/09_original_figure9_reference.csv",
      "09_assets/reference/09_major3_bootstrap_reference.csv",
      "09_assets/data/09_source_input_manifest.csv")
    valid <- vapply(candidates, function(dir) all(file.exists(file.path(dir, required_relative))), logical(1))
    check(any(valid), "Cannot find the complete 09_assets folder next to the current script. ",
          "Use source('FULL/PATH/09_hubgenes_immune_correlation.R') or set ",
          "PD_OSCC_09_PACKAGE_DIR to the folder containing 09_assets. Checked: ",
          paste(candidates, collapse = ' | '))
    package_dir <- normalizePath(candidates[which(valid)[1L]], winslash = "/", mustWork = TRUE)
    files <- c(
      PD_analysis_inputs = file.path(package_dir, "09_assets/data/09_PD_analysis_inputs.rds"),
      OSCC_analysis_inputs = file.path(package_dir, "09_assets/data/09_OSCC_analysis_inputs.rds"),
      original_reference = file.path(package_dir, "09_assets/reference/09_original_figure9_reference.csv"),
      bootstrap_reference = file.path(package_dir, "09_assets/reference/09_major3_bootstrap_reference.csv"),
      source_manifest = file.path(package_dir, "09_assets/data/09_source_input_manifest.csv")
    )
    check(all(file.exists(files)), "Missing Step 09 package inputs: ",
          paste(files[!file.exists(files)], collapse = ", "))
    hashes <- unname(tools::md5sum(files)); check(!anyNA(hashes), "Cannot fingerprint packaged inputs")
    # Relative paths in the run manifest avoid disclosing local account names.
    manifest <- data.frame(Role = names(files), File = unname(files),
                           MD5 = hashes, stringsAsFactors = FALSE)
    manifest$File <- c("09_assets/data/09_PD_analysis_inputs.rds",
      "09_assets/data/09_OSCC_analysis_inputs.rds",
      "09_assets/reference/09_original_figure9_reference.csv",
      "09_assets/reference/09_major3_bootstrap_reference.csv",
      "09_assets/data/09_source_input_manifest.csv")
    genes <- c("CXCL1", "CXCL6", "CXCR4", "RGS1", "CHST2")
    shared_cells <- c("Dendritic cells activated", "Dendritic cells resting", "Macrophages M0",
      "Macrophages M1", "Macrophages M2", "Mast cells activated", "NK cells resting",
      "Neutrophils", "Plasma cells", "T cells CD4 memory activated",
      "T cells CD4 memory resting", "T cells CD4 naive", "T cells follicular helper")
    ref <- utils::read.csv(files[["original_reference"]], check.names = FALSE)
    boot_ref <- utils::read.csv(files[["bootstrap_reference"]], check.names = FALSE)
    check(nrow(ref) == 65L && nrow(boot_ref) == 65L &&
            all(c("Gene", "Cell_type", "r_PD", "r_OSCC", "p_Fisher", "p_Fisher_adj") %in% names(ref)) &&
            !anyDuplicated(paste(ref$Gene, ref$Cell_type, sep = "\r")) &&
            !anyDuplicated(paste(boot_ref$Gene, boot_ref$Cell_type, sep = "\r")) &&
            setequal(unique(ref$Gene), genes) && setequal(unique(ref$Cell_type), shared_cells),
          "Reference comparison tables are invalid.")
    pd <- readRDS(files[["PD_analysis_inputs"]]); os <- readRDS(files[["OSCC_analysis_inputs"]])
    validate_cohort <- function(dat, disease, expected_n, expected_cells) {
      check(is.list(dat) && all(c("historical_expr", "sensitivity_expr", "immune", "unit", "map", "probe_map") %in% names(dat)),
            disease, ": invalid input schema")
      x <- dat$historical_expr; y <- dat$sensitivity_expr; im <- dat$immune
      check(is.matrix(x) && is.numeric(x) && is.matrix(y) && is.numeric(y) &&
              is.matrix(im) && is.numeric(im) && nrow(x) == expected_n &&
              identical(dim(x), dim(y)) && identical(colnames(x), genes) &&
              identical(colnames(y), genes) && identical(rownames(x), rownames(y)) &&
              nrow(im) == expected_n && ncol(im) == expected_cells &&
              !anyDuplicated(rownames(x)) && all(is.finite(x)) && all(is.finite(y)) &&
              all(is.finite(im)) && all(im >= 0) && all(shared_cells %in% colnames(im)),
            disease, ": invalid expression or immune matrices")
      check(length(dat$unit) == expected_n && !anyNA(dat$unit) && all(nzchar(dat$unit)) &&
              all(c("GSM", "Group", "Analysis_unit") %in% names(dat$map)) &&
              identical(as.character(dat$map$GSM), rownames(x)) &&
              identical(as.character(dat$map$Analysis_unit), dat$unit) &&
              setequal(dat$map$Group, c("Control", disease)), disease, ": invalid sample metadata")
      check(nrow(dat$probe_map) == length(genes) && identical(as.character(dat$probe_map$Gene), genes),
            disease, ": invalid probe mapping")
      if (disease == "PD") check(length(unique(dat$unit)) == 90L, "PD requires 90 patient clusters")
      else check(!anyDuplicated(dat$unit), "OSCC requires one sampling unit per sample")
      check(all(colMeans(im == 0) <= 0.5), disease, ": retained cell has >50% zeros")
      invisible(NULL)
    }
    validate_cohort(pd, "PD", 247L, 16L); validate_cohort(os, "OSCC", 212L, 14L)
    check(setequal(intersect(colnames(pd$immune), colnames(os$immune)), shared_cells),
          "Expected 13 shared cell types")
    cat("Input check: PD 247 biopsies / 90 patients / 16 retained cells; ",
        "OSCC 212 samples / 14 retained cells.\n", sep = "")
    cat("Historical branch and sensitivity branch use separate expression matrices.\n")
    if (mode == "check_inputs") {
      cat("INPUT CHECK COMPLETE. No correlations or bootstrap calculated; no files written.\n")
      return(invisible(NULL))
    }
    base <- Sys.getenv("PD_OSCC_09_OUTPUT_BASE", file.path(package_dir, "results/09_hubgenes_immune_correlation"))
    dir.create(base, recursive = TRUE, showWarnings = FALSE)
    run <- file.path(base, paste0(format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"),
                                  "Z_pid", Sys.getpid()))
    check(!dir.exists(run) && dir.create(run), "Cannot create unique output directory")
    status <- file.path(run, "09_run_status.txt")
    writeLines("STARTED", status)
    on.exit({if (file.exists(status) && identical(readLines(status, warn = FALSE)[1L], "STARTED"))
      writeLines("INCOMPLETE: inspect console and partial outputs", status)}, add = TRUE)
    write_csv <- function(x, name) utils::write.csv(x, file.path(run, name), row.names = FALSE)
    write_csv(manifest, "09_input_manifest.csv")
    write_csv(rbind(pd$probe_map, os$probe_map), "09_historical_probe_selection.csv")
    writeLines(capture.output(utils::sessionInfo()), file.path(run, "09_sessionInfo.txt"))
    # STATISTICAL_CALCULATION_INSERTION_POINT

    # Historical Spearman correlations and cor.test p-values, first over all
    # retained cells in each cohort; BH adjustment is a separate audit column.
    correlation_table <- function(dat, disease) {
      cells <- colnames(dat$immune)
      res <- do.call(rbind, lapply(genes, function(g) do.call(rbind, lapply(cells, function(cell) {
        test <- suppressWarnings(stats::cor.test(dat$historical_expr[, g], dat$immune[, cell],
                                                 method = "spearman"))
        data.frame(Disease = disease, Gene = g, Cell_type = cell,
                   Rho = unname(test$estimate), Pvalue = test$p.value)
      }))))
      res$P_BH_all_pairs <- stats::p.adjust(res$Pvalue, method = "BH")
      res$Abs_rho_gt_030 <- abs(res$Rho) > 0.3
      res$BH_all_and_abs_rho_gt_030 <- res$P_BH_all_pairs < 0.05 & res$Abs_rho_gt_030
      res
    }
    hist_pd <- correlation_table(pd, "PD"); hist_os <- correlation_table(os, "OSCC")
    write_csv(hist_pd, "09_original_PD_80_pairs.csv")
    write_csv(hist_os, "09_original_OSCC_70_pairs.csv")
    key <- function(d) paste(d$Gene, d$Cell_type, sep = "\r")
    grid <- expand.grid(Gene = genes, Cell_type = shared_cells, stringsAsFactors = FALSE)
    idx_pd <- match(key(grid), key(hist_pd)); idx_os <- match(key(grid), key(hist_os))
    check(!anyNA(idx_pd) && !anyNA(idx_os), "Historical shared pairs missing")
    grid$rho_PD <- hist_pd$Rho[idx_pd]; grid$rho_OSCC <- hist_os$Rho[idx_os]
    check(all(abs(grid$rho_PD) < 1) && all(abs(grid$rho_OSCC) < 1), "Fisher z is undefined for |rho|=1")
    grid$delta_OSCC_minus_PD <- grid$rho_OSCC - grid$rho_PD
    z <- (atanh(grid$rho_OSCC) - atanh(grid$rho_PD)) / sqrt(1 / (247 - 3) + 1 / (212 - 3))
    grid$P_Fisher_approx <- 2 * stats::pnorm(-abs(z))
    grid$P_Fisher_BH <- stats::p.adjust(grid$P_Fisher_approx, method = "BH")
    grid$Fisher_BH_lt_005 <- grid$P_Fisher_BH < 0.05
    grid$Opposite_sign_and_Fisher_BH_lt_005 <- grid$rho_PD * grid$rho_OSCC < 0 & grid$Fisher_BH_lt_005
    write_csv(grid, "09_original_Fisher_BH_65_pairs.csv")
    ri <- match(key(grid), key(ref)); check(!anyNA(ri), "Original reference pairs missing")
    hist_diff <- c(rho_PD = max(abs(grid$rho_PD - ref$r_PD[ri])),
                   rho_OSCC = max(abs(grid$rho_OSCC - ref$r_OSCC[ri])),
                   Fisher_p = max(abs(grid$P_Fisher_approx - ref$p_Fisher[ri])),
                   Fisher_BH = max(abs(grid$P_Fisher_BH - ref$p_Fisher_adj[ri])))
    hist_counts <- c(PD_abs_rho_gt_030 = sum(hist_pd$Abs_rho_gt_030),
                     OSCC_abs_rho_gt_030 = sum(hist_os$Abs_rho_gt_030),
                     PD_BH_and_magnitude = sum(hist_pd$BH_all_and_abs_rho_gt_030),
                     OSCC_BH_and_magnitude = sum(hist_os$BH_all_and_abs_rho_gt_030),
                     Fisher_BH_65_pairs = sum(grid$Fisher_BH_lt_005))

    # Replicate Major3's seed, sequence, cluster handling and cell order.
    # PD resamples all biopsies from each drawn patient; OSCC resamples samples.
    # Pair independent cohort bootstrap replicate b to form delta_b.
    nrep <- 1200L
    rho <- function(x, y) suppressWarnings(stats::cor(x, y, method = "spearman", use = "pairwise.complete.obs"))
    calc_boot <- function(dat) {
      indices <- split(seq_along(dat$unit), dat$unit); ncluster <- length(indices)
      out <- array(NA_real_, dim = c(length(genes), length(shared_cells), nrep),
                   dimnames = list(genes, shared_cells, NULL))
      for (b in seq_len(nrep)) {
        sampled <- unlist(indices[sample.int(ncluster, ncluster, replace = TRUE)], use.names = FALSE)
        out[, , b] <- rho(dat$sensitivity_expr[sampled, , drop = FALSE],
                           dat$immune[sampled, shared_cells, drop = FALSE])
      }
      out
    }
    old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (old_seed_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({if (old_seed_exists) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)}, add = TRUE)
    set.seed(20261023L)
    cat("Running 1,200 PD patient-cluster and OSCC sample bootstrap replicates...\n")
    bp <- calc_boot(pd); bo <- calc_boot(os)
    obs_pd <- rho(pd$sensitivity_expr, pd$immune[, shared_cells, drop = FALSE])
    obs_os <- rho(os$sensitivity_expr, os$immune[, shared_cells, drop = FALSE])
    mean_rho <- function(dat) {
      ids <- split(seq_along(dat$unit), dat$unit)
      means <- function(x) do.call(rbind, lapply(ids, function(i) colMeans(x[i, , drop = FALSE])))
      rho(means(dat$sensitivity_expr), means(dat$immune[, shared_cells, drop = FALSE]))
    }
    mean_pd <- mean_rho(pd); mean_os <- mean_rho(os)
    interval <- function(v) {
      v <- v[is.finite(v)]
      if (length(v) < as.integer(nrep * 0.9)) return(c(NA_real_, NA_real_))
      unname(stats::quantile(v, c(.025, .975), na.rm = TRUE))
    }
    boot <- do.call(rbind, lapply(genes, function(g) do.call(rbind, lapply(shared_cells, function(cell) {
      p <- bp[g, cell, ]; o <- bo[g, cell, ]; p_ci <- interval(p); o_ci <- interval(o); d_ci <- interval(o - p)
      data.frame(Gene = g, Cell_type = cell, rho_PD = unname(obs_pd[g, cell]),
        rho_OSCC = unname(obs_os[g, cell]), delta_OSCC_minus_PD = unname(obs_os[g, cell] - obs_pd[g, cell]),
        PD_cluster_CI_low = p_ci[1], PD_cluster_CI_high = p_ci[2],
        OSCC_CI_low = o_ci[1], OSCC_CI_high = o_ci[2], delta_boot_CI_low = d_ci[1],
        delta_boot_CI_high = d_ci[2],
        delta_CI_excludes_zero = if (all(is.finite(d_ci))) d_ci[1] > 0 || d_ci[2] < 0 else NA,
        rho_PD_patient_averaged = unname(mean_pd[g, cell]),
        rho_OSCC_patient_averaged = unname(mean_os[g, cell]),
        n_PD_biopsies = nrow(pd$sensitivity_expr), n_PD_patients = length(unique(pd$unit)),
        n_OSCC_samples = nrow(os$sensitivity_expr), bootstrap_replicates = nrep)
    }))))
    write_csv(boot, "09_sensitivity_65_pairs_patient_bootstrap.csv")
    write_csv(boot[boot$Cell_type == "Macrophages M1", , drop = FALSE], "09_sensitivity_M1_five_genes.csv")
    bi <- match(key(boot), key(boot_ref)); check(!anyNA(bi), "Bootstrap reference pairs missing")
    numeric_cols <- c("rho_PD", "rho_OSCC", "delta_OSCC_minus_PD", "PD_cluster_CI_low",
      "PD_cluster_CI_high", "OSCC_CI_low", "OSCC_CI_high", "delta_boot_CI_low", "delta_boot_CI_high",
      "rho_PD_patient_averaged", "rho_OSCC_patient_averaged")
    check(all(numeric_cols %in% names(boot_ref)), "Bootstrap reference columns missing")
    boot_diff <- vapply(numeric_cols, function(nm) max(abs(boot[[nm]] - boot_ref[[nm]][bi])), numeric(1))
    checks <- data.frame(Check = c(names(hist_diff), names(hist_counts), paste0("bootstrap_", names(boot_diff)),
                                   "bootstrap_CI_excludes_zero"),
                         Observed = c(hist_diff, hist_counts, boot_diff, sum(boot$delta_CI_excludes_zero, na.rm = TRUE)),
                         Reference = c(rep(0, length(hist_diff)), 26, 22, 26, 22, 36,
                                       rep(0, length(boot_diff)), 38),
                         Criterion = c(rep("max difference <= 1e-10", length(hist_diff)), rep("count", length(hist_counts)),
                                       rep("max difference <= 1e-10", length(boot_diff)), "count"),
                         row.names = NULL)
    checks$Match <- ifelse(checks$Criterion == "count", checks$Observed == checks$Reference,
                           checks$Observed <= 1e-10)
    write_csv(checks, "09_reference_comparison.csv")
    issues <- c(
      "Historical rho and Fisher-BH and updated 01D bootstrap must be interpreted as separate analyses.",
      "The Fisher transformation is approximate for Spearman and treats 247 PD biopsies as independent.",
      "Bootstrap confidence intervals are nominal percentile intervals without FDR correction.",
      "OSCC samples were resampled individually because patient identifiers were not available.",
      "LM22 CIBERSORTx values are relative estimates; correlations do not establish immune reprogramming.",
      "Original 26/22 counts were saved by filtering |rho| > 0.3; the BH audit is calculated separately.",
      "Historical five-gene Figure 9 data use normalized Step 01; sensitivity uses GEO-processed 01D data."
    )
    writeLines(issues, file.path(run, "09_interpretation_and_limitations.txt"))
    writeLines(c("COMPLETED", paste("All reference checks match:", all(checks$Match)),
                 "Review 09_reference_comparison.csv and 09_interpretation_and_limitations.txt."), status)
    cat("\n=== REFERENCE COMPARISON ===\n"); print(checks, row.names = FALSE, digits = 5)
    cat("\nM1 bootstrap sensitivity:\n"); print(boot[boot$Cell_type == "Macrophages M1",
      c("Gene", "delta_OSCC_minus_PD", "delta_boot_CI_low", "delta_boot_CI_high")], row.names = FALSE)
    cat("\nResults saved to:", run, "\n")
    invisible(run)
  }
  main()
})
