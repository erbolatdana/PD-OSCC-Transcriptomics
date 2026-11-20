# Step 09: five-gene and immune-cell correlations

This folder contains an analysis-ready reproduction of the two distinct Step 09 analyses reported in the manuscript. `09_hubgenes_immune_correlation.R` uses only base R, reads the included files in `09_assets`, and does not download data, rerun CIBERSORTx, perform gene selection, or install packages.

Run from this folder with R:

```r
source("09_hubgenes_immune_correlation.R", encoding = "UTF-8")
```

The run creates a new time-stamped directory in `results/09_hubgenes_immune_correlation/`. Earlier outputs are not overwritten. For input checks only, set `Sys.setenv(PD_OSCC_09_MODE = "check_inputs")` before sourcing. An alternative output parent can be supplied via `PD_OSCC_09_OUTPUT_BASE`.

## Included data and provenance

The two RDS files in `09_assets/data/` contain, for GSE10334 (PD) and GSE30784 (OSCC), five-gene expression for the historical and sensitivity branches, retained CIBERSORTx immune-cell fractions, sample IDs, group labels and sampling units. These are small, analysis-ready subsets; they are not the full source datasets and are not a substitute for the original CIBERSORTx input and output provenance. Original PD patient identifiers were replaced by pseudonyms, retaining the original clusters and their sorted order so that seeded cluster bootstrap draws are reproducible. OSCC sampling units are individual GSM samples because patient IDs were unavailable.

`09_assets/source/09_prepare_inputs_from_archived_sources.R` records how the inputs were extracted and checked from the archived project files, without running CIBERSORTx or statistical analysis. It is provenance code and requires the original Step 01/03 and 01D files; these full source files are not included in this package. `09_source_input_manifest.csv` records their basenames and MD5 fingerprints. The input file `09_original_figure9_reference.csv` reproduces the historical Figure 9 / Table S4 comparison data; `09_major3_bootstrap_reference.csv` holds the later saved sensitivity results. They are comparison references, not independent external data.

Historical branch expression comes from the normalized, probe-level Step 01 data, with a probe chosen per cohort by the highest within-cohort mean for each gene. Sensitivity-branch expression comes from the later 01D gene-level matrices, which retained GEO-processed expression without additional between-array normalization. The two matrices differ and must not be interchanged.

## Statistical methods and interpretation

Historical within-disease Spearman coefficients and tests use 16 retained cell types for 247 PD biopsies and 14 for 212 OSCC samples. The script reports both the historical `|rho| > 0.3` counts and BH-adjusted p-values across all tested gene-cell pairs within each cohort. The 13 shared cell types produce 65 between-cohort comparisons. Fisher's z transformation is applied to Spearman correlations as an approximation, using sample sizes 247 and 212; the 65 p-values are BH-adjusted. A sign reversal alone is not a statistically significant difference.

The historical comparison treats PD biopsies as independent despite 90 patient clusters. The separate sensitivity analysis resamples 90 PD patient clusters and OSCC samples independently for 1,200 replicates, computes `rho(OSCC) - rho(PD)`, and reports percentile intervals. These intervals are nominal and not adjusted for the 65 comparisons. `38` intervals excluding zero must not be described as 38 FDR-significant results. Spearman correlations and estimated relative immune-cell fractions do not establish immune reprogramming or causality.

In the checked run, historical Figure 9 correlation and Fisher-BH tables match the archived reference within 1e-10; the counts are 26 PD and 22 OSCC pairs above the magnitude threshold, and 36 of 65 Fisher-BH comparisons below 0.05. The separately calculated 1,200-replicate sensitivity values match the archived table within 1e-10; 38 of 65 nominal intervals exclude zero. These comparisons verify the computation against saved study results, not biological validity or independent replication.

## Limitations and further audit

The package starts from precomputed, study-specific analysis inputs. Full raw-data-to-CIBERSORTx reproduction requires the previously used GEO datasets, preprocessing, LM22 and archived CIBERSORTx jobs and should be documented alongside the full repository. The exact CIBERSORTx upload matrices have not been independently proven identical to the historical job inputs. The calculations do not by themselves regenerate the final Figure 9 graphical layout. Different R versions, numerical libraries or RNG settings may produce small differences; inspect `09_reference_comparison.csv` and `09_run_status.txt`, rather than changing methods to force a match. Review the included sample-level data and data-use terms before public release.
