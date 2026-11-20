# Step 10: exploratory five-gene GO and KEGG enrichment

`10_GO_KEGG_hubgenes.R` analyses the historical five-gene panel: CXCL1, CXCL6, CXCR4, RGS1 and CHST2. This analysis is distinct from the disease-specific DEG enrichment in Step 04. An earlier manuscript version reported 16 PD and 788 OSCC GO terms for the DEG lists, not this five-gene analysis. The separate Step 04 audit found 18 and 729; the author subsequently showed revised text containing those numbers. Check the saved final manuscript and Figure 4 against the Step 04 results. The reason for the historical difference remains unverified and is not resolved by Step 10.

To run from the folder containing this README and `10_assets/`, install the R/Bioconductor packages `clusterProfiler` and `org.Hs.eg.db`, ensure that the KEGG REST service is reachable, and execute:

```r
source("10_GO_KEGG_hubgenes.R", encoding = "UTF-8")
```

The script does not install packages. It creates a new time-stamped folder under `results/10_five_gene_enrichment/`, records the R session, input hashes and settings, and saves symbol–Entrez mappings, GO and KEGG tables and R objects. If the KEGG request fails, the run is marked incomplete; previously written GO and mapping outputs are retained. To check the package without computing or accessing KEGG, set `Sys.setenv(PD_OSCC_10_MODE = "check_inputs")` before sourcing. Use `Sys.unsetenv("PD_OSCC_10_MODE")` to restore the default run mode.

## Methods and interpretation

The five genes are converted from SYMBOL to Entrez IDs with the installed `org.Hs.eg.db`; the archived five mappings are included for comparison. The script stops for review if the current mapping is incomplete or not one-to-one.

The original Step 10 settings are retained: `enrichGO` for Biological Process with BH adjustment, `pvalueCutoff = 0.05`, `qvalueCutoff = 0.2`, and readable gene symbols; `enrichKEGG` for organism `hsa` with BH adjustment and `pvalueCutoff = 0.05`, using package defaults for other parameters. An explicit gene universe was not supplied in the historical code; the database/function defaults are used, not a prespecified set of genes measured on GPL570. KEGG is queried online and GO annotations depend on the installed package. Both the database version and background affect the results.

The panel was selected from the study data, contains only five genes, and is not an independent validation dataset. Treat enriched terms as exploratory associations, not evidence of a causal pathway or independent confirmation. The current manuscript describes DEG enrichment but does not explicitly report this five-gene GO/KEGG analysis; including additional results in the manuscript would require an explicit Methods/Results update and review.

## Archived and checked results

`10_assets/10_archived_five_gene_KEGG_result.csv` contains 25 rows recovered from the saved historical `kegg_hub` object; seven have adjusted p < 0.05. `10_assets/10_archived_five_gene_panel.csv` and `10_archived_SYMBOL_ENTREZ_map.csv` record its panel and saved mappings. The original five-gene GO result object was not found in the inspected archive, so a historical GO comparison is unavailable.

The saved test run in `10_assets/verified_example/` was performed on 2026-11-19 with R 4.6.1. It returned 38 GO BP rows with adjusted p < 0.05; these are current five-gene results, not a reproduced historical GO count. The seven significant KEGG pathway IDs, names and contributing-gene sets matched the archived result. Their adjusted p-values were close but not identical. The background recorded in the result changed from 9,416 to 9,472 annotated IDs, which is relevant to the comparison; the available evidence does not establish it as the sole cause of every numerical difference. Refer to the saved session information and comparison tables. Do not tune thresholds or inputs to recover prior results.

This package provides code and archived comparison tables for Step 10; it does not freeze KEGG's online database or independently reproduce earlier GO results. Before submission, review the broader repository's data provenance, database versions and the manuscript's Step 04 GO count discrepancy.
