# Steps 01–05: internal reproducibility draft

**Draft for review, not a complete or publication-ready repository.** These are
copies of the audited project scripts and selected supporting files. Scripts 02, 03, 05 and 05B completed in isolated local tests with saved
inputs and selected results matched archived comparisons. Script 04 produced
18/729 GO terms but its KEGG portion did not complete in that test. The full
workflow has not been tested end-to-end or on another computer.
Review all contents, sample-level data, rights to redistribute third-party
materials, and manuscript agreement before publication.

## Files and dependencies

Run scripts from the project root, with `scripts/`, `data/`, and `results/`
at that root. The current directory must be selected by the user. In this draft,
comments in scripts 03, 04 and 05B were changed to remove a computer-specific
working-directory example; their executable code was not changed. Script 01
downloads GEO series GSE10334 and GSE30784 with GPL570
annotation and reads `data/metadata/` sample maps (GSM and, for PD, coded
patient IDs). Scripts 02 and 04 consume outputs from earlier steps. Script 05B
is an optional audit of Figure 5 references, after Script 05.

Core packages cited in these scripts include GEOquery, Biobase, limma,
ggplot2, clusterProfiler, and org.Hs.eg.db. Other package use can be conditional
(including dplyr and writexl); package versions and a complete installation
guide still need verification. Network access is needed for GEO and live KEGG.
Database and package versions can affect results.

Script 03 requires six historical CIBERSORTx files in
`data/cibersortx_historical/` and a saved preprocessing object. Two archived results and two email-redacted
runtime logs are supplied directly. Restore the preprocessing object and the
following two mixture files from the included parts by running
`source("restore_large_files.R")` from the repository root; see the
[restore instructions](../RESTORE_LARGE_FILES.md):

* `CIBERSORT_input_PD_linear.txt` (119,780,832 bytes)
* `CIBERSORT_input_OSCC_linear.txt` (103,006,748 bytes)

The restoration script checks file sizes and SHA-256 hashes. Restoration and
hash checks were tested locally. Rights and privacy review remain unresolved,
including when files are distributed as binary parts. Do not substitute the
20,848-gene matrices from Script 01: the historical mixture files differ.
Script 05 also needs Script 03's `03_cell_zero_QC.csv`, Script 01's expression
objects and Script 02's DEG lists. The stages require the R environment,
packages, network downloads and earlier outputs described above; the draft
has not been tested end-to-end.

## Archived logs and integrity

Runtime-log copies here replace email addresses with `[REDACTED_EMAIL]`;
the original logs remain unchanged in the working project. The file
`source_manifest.csv` records the MD5 hashes and byte sizes of the **original**
six historical files, not the redacted copies. `redacted_log_manifest.csv`
separately records each original MD5 and the corresponding edited copy MD5.
Script 03 records current input hashes in its own output manifest and checks
for the expected settings strings in the logs; its code does not read the
original `source_manifest.csv` or compare logs byte-for-byte. Script 03 completed
in a local isolated test with saved inputs and selected results matched archived
comparisons; this does not establish the precise identity of the files uploaded
for the historical CIBERSORTx jobs.

## Scientific scope and unresolved checks

* Script 03 reads historical CIBERSORTx results; it does not rerun deconvolution.
  Archived logs record LM22, quantile normalization enabled and one permutation;
  this is insufficient for reliable permutation-based CIBERSORTx P values.
  The exact relationship of the historical mixture files to the submitted jobs
  remains to be documented; their identity must not be assumed.
* The default Script 02 PD analysis reproduces an unblocked biopsy-level model;
  patient-blocked PD analysis is a separately labelled sensitivity branch.
* An earlier manuscript file reported 16 PD / 788 OSCC significant GO BP terms;
  the recorded repeat found 18 / 729. The author subsequently showed revised
  text with 18 / 729. Verify the saved final manuscript and Figure 4 against
  the results; the original reason for the difference is not established.
* Script 05 distinguishes the historical multiple-testing calculation from
  correction across all tested pairs. Historical PD analyses use biopsy-level
  observations and do not account for repeated samples from patients.
* Sample maps and deconvolution results are sample-level data, and the archived
  result/log redistribution status is not yet determined. Review before sharing.
  Record environment versions and document how to obtain all inputs.

Original sources are preserved in the working project. No R workspace files
were copied into this draft.
