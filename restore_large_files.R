# Restore the three supporting files from binary parts.
# Run from the repository root: source("restore_large_files.R")
# Required package: openssl (https://cran.r-project.org/package=openssl).
# Existing destination files are never overwritten.

local({
  if (!requireNamespace("openssl", quietly = TRUE)) {
    stop("Install the R package 'openssl' before running this script.")
  }
  manifest_path <- "large_files_manifest.csv"
  index_path <- "largefiles/parts_manifest.csv"
  if (!file.exists(manifest_path) || !file.exists(index_path)) {
    stop("Run this script from the repository root, with largefiles/ present.")
  }
  expected <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
  parts <- read.csv(index_path, stringsAsFactors = FALSE, check.names = FALSE)
  required_paths <- c(
    "step06/results/01_preprocessing/01_preprocessed_esets.rds",
    "early/data/cibersortx_historical/CIBERSORT_input_PD_linear.txt",
    "early/data/cibersortx_historical/CIBERSORT_input_OSCC_linear.txt"
  )
  if (anyDuplicated(expected$relative_path) ||
      !all(required_paths %in% expected$relative_path) ||
      !setequal(parts$relative_path, required_paths)) {
    stop("Unexpected file manifest; no files restored.")
  }
  hash_file <- function(path) {
    input <- file(path, "rb")
    on.exit(close(input))
    unclass(as.character(openssl::sha256(input)))
  }
  copy_stream <- function(input, output) {
    repeat {
      bytes <- readBin(input, what = "raw", n = 1024L * 1024L)
      if (!length(bytes)) break
      writeBin(bytes, output)
    }
  }
  for (dest in required_paths) {
    row <- expected[expected$relative_path == dest, , drop = FALSE]
    entries <- parts[parts$relative_path == dest, , drop = FALSE]
    entries <- entries[order(entries$part_number), , drop = FALSE]
    if (!identical(as.integer(entries$part_number), seq_len(nrow(entries))) ||
        any(!grepl("^largefiles/[a-zA-Z0-9_.-]+$", entries$part_file)) ||
        length(unique(entries$encoding)) != 1L ||
        !unique(entries$encoding) %in% c("raw", "gzip")) {
      stop("Invalid parts manifest for ", dest)
    }
    if (file.exists(dest)) {
      if (!identical(tolower(hash_file(dest)), tolower(row$sha256))) {
        stop("Existing destination differs; refusing to overwrite: ", dest)
      }
      cat("Already present and verified: ", dest, "\n", sep = "")
      next
    }
    if (!all(file.exists(entries$part_file))) stop("Missing parts for ", dest)
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    combined <- tempfile("combined_")
    restored <- tempfile("restore_", tmpdir = dirname(dest))
    tryCatch({
      out <- file(combined, "wb")
      tryCatch({
        for (part in entries$part_file) {
          input <- file(part, "rb")
          tryCatch(copy_stream(input, out), finally = close(input))
        }
      }, finally = close(out))
      input <- if (entries$encoding[1] == "gzip") gzfile(combined, "rb") else file(combined, "rb")
      output <- file(restored, "wb")
      tryCatch(copy_stream(input, output), finally = {close(input); close(output)})
      restored_size <- file.info(restored)$size
      restored_hash <- hash_file(restored)
      if (restored_size != row$size_bytes ||
          !identical(tolower(restored_hash), tolower(row$sha256))) {
        stop("Restored file failed size/SHA-256 verification: ", dest,
             " (size ", restored_size, "; SHA-256 ", restored_hash, ")")
      }
      if (file.exists(dest) || !file.rename(restored, dest)) {
        stop("Could not safely move verified file into place: ", dest)
      }
      cat("Restored and SHA-256 verified: ", dest, "\n", sep = "")
    }, finally = {unlink(combined); unlink(restored)})
  }
  # Scripts 03 and 06 use the same saved RDS under separate working roots.
  source_rds <- required_paths[1]
  early_rds <- "early/results/01_preprocessing/01_preprocessed_esets.rds"
  expected_rds <- expected[expected$relative_path == source_rds, , drop = FALSE]
  if (file.exists(early_rds)) {
    if (file.info(early_rds)$size != expected_rds$size_bytes ||
        !identical(tolower(hash_file(early_rds)), tolower(expected_rds$sha256))) {
      stop("Existing destination differs; refusing to overwrite: ", early_rds)
    }
    cat("Already present and verified: ", early_rds, "\n", sep = "")
  } else {
    dir.create(dirname(early_rds), recursive = TRUE, showWarnings = FALSE)
    copy_temp <- tempfile("copy_", tmpdir = dirname(early_rds))
    tryCatch({
      if (!file.copy(source_rds, copy_temp, overwrite = FALSE) ||
          file.info(copy_temp)$size != expected_rds$size_bytes ||
          !identical(tolower(hash_file(copy_temp)), tolower(expected_rds$sha256))) {
        stop("RDS copy failed verification: ", early_rds)
      }
      if (file.exists(early_rds) || !file.rename(copy_temp, early_rds)) {
        stop("Could not safely move verified RDS copy into place: ", early_rds)
      }
      cat("Copied and SHA-256 verified: ", early_rds, "\n", sep = "")
    }, finally = unlink(copy_temp))
  }
  cat("All three source files and the additional RDS copy are verified.\n")
})
