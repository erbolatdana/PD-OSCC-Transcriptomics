# Restore supporting files

Run this step only after obtaining the complete `largefiles/` directory and
`large_files_manifest.csv`. From the repository root in R, install the CRAN
package `openssl` if needed, then run:

```r
source("restore_large_files.R")
```

The script combines binary parts, decompresses the text files and compares each
restored file's size and SHA-256 with `large_files_manifest.csv`. It places the
preprocessing RDS under both `step06/results/01_preprocessing/` and
`early/results/01_preprocessing/`, and the archived mixtures under
`early/data/cibersortx_historical/`. It verifies both copies of the RDS.

An existing file with a matching hash is kept; a differing file causes an error
and is not overwritten. The preprocessing RDS is a saved intermediate and does
not verify reconstruction of the original analysis from raw GEO data. The
exact identity of archived mixtures with historical submitted files has not been
independently established. Redistribution rights for sample-level data and
third-party materials must be reviewed before release.
