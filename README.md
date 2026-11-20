# PD–OSCC analysis code: reviewer draft

This draft collects R analysis code and selected supporting materials for the
periodontitis (PD) and oral squamous cell carcinoma (OSCC) study. **It is not a
complete raw-data-to-results workflow.** Three large files needed for some
checks are deliberately omitted; their names and SHA-256 hashes are in
[`large_files_manifest.csv`](large_files_manifest.csv). No download links for
those files have yet been provided. Data sharing and third-party terms still
require review. The package has not been tested on another computer.

## Where to start

| Folder | Material | Scope and prerequisites |
| --- | --- | --- |
| `early/` | Scripts 01–05 and optional Figure 5 QC; sample maps, reference tables and archived CIBERSORTx outputs. | Script 01 downloads GEO datasets GSE10334 and GSE30784 and GPL570 annotation. Subsequent steps need its saved outputs. Script 03 additionally expects the two omitted CIBERSORTx mixture files. |
| `step06/` | Historical WGCNA reconstruction, feature and sample manifests, DEG object and archived comparison records. | The 195 MB saved preprocessing RDS is omitted. Input checks and network rebuilding cannot run until the matching file is supplied in the path recorded in `large_files_manifest.csv`. |
| `step07/` | Fixed-panel cross-validation and external OSCC validation, with source scripts and selected expression and sample data. | The provided data include processed expression matrices and the external cohort ExpressionSet; inspect the sample-level data and redistribution terms. See `step07/README_RU.md`. |
| `steps08_10/` | Description of the manual Step 08, R scripts 09–10, analysis-ready subsets and references. | Step 08 is a methods record, not an executable script. See `09_README.md` and `10_README.md` in that directory. |

### Files not included in this draft

The following files are not present in this repository and no download location
is currently available. Their exact expected paths, sizes and SHA-256 hashes
are recorded in [`large_files_manifest.csv`](large_files_manifest.csv).

| Missing file | Role |
| --- | --- |
| `step06/results/01_preprocessing/01_preprocessed_esets.rds` | Saved preprocessing object required by the historical Step 06 reconstruction and its input check. |
| `early/data/cibersortx_historical/CIBERSORT_input_PD_linear.txt` | Archived PD mixture required by Script 03's historical-input checks. |
| `early/data/cibersortx_historical/CIBERSORT_input_OSCC_linear.txt` | Archived OSCC mixture required by Script 03's historical-input checks. |

Script 03 also requires the preprocessing object at
`early/results/01_preprocessing/01_preprocessed_esets.rds`; generating or
obtaining that object and the two mixture files is necessary before Script 03
can run. Step 07's included processed expression matrices are different inputs
and are not substitutes for this object. A matching saved preprocessing object
can support downstream checks, but does not independently verify Script 01's
raw-data-to-result workflow.

Each top-level folder is a separate working root, not part of one automated
pipeline. Keep the `09_assets` and `10_assets` folders next to their scripts.
Consult the script headers and step-specific documents for packages, inputs,
output locations, options and analysis modes. Installing R and the required
packages and supplying missing files are the user's responsibility. A
complete, version-locked environment is not provided.

For example, after setting the R working directory to `step07/`, run
`source("07_machine_learning_reproduction.R")` for its default input-only
check. Set `PD_OSCC_07_MODE` to `run` to run the fixed-panel analyses, and
unset it afterward; see [`step07/README_RU.md`](step07/README_RU.md). For
Steps 09 and 10, set the working directory to `steps08_10/` and follow
[`09_README.md`](steps08_10/09_README.md) and
[`10_README.md`](steps08_10/10_README.md). Step 10 needs access to KEGG for
its enrichment run. The included step-specific READMEs describe the run
modes and where outputs are written; Script 09 runs its analysis by default.

## What has been checked

The earlier local draft was tested on the author's Windows system, not from
this smaller distribution on another computer. Scripts 02, 03, 05 and 05B
completed in an isolated local test using saved inputs; selected results
matched their archived comparisons. Script 04 produced 18 PD and 729 OSCC GO
terms in that test, but its KEGG portion did not complete there. Script 07
reproduced five archived AUC summary values locally. Input checks for the full
Step 06, 09 and 10 drafts were performed locally; Step 06's full input is not
present here. These checks are not independent validation of the whole
workflow or guarantees of identical results on another machine.

The author's later connection test reached GEO and KEGG. However, downloading
the GSE10334 matrix timed out after receiving approximately 3 MB of a 44.5 MB
file. Thus end-to-end execution of Script 01 from GEO remains unverified in
that test. GEO downloads and current KEGG annotations can depend on connection
speed, service availability and database versions.

## Important provenance and interpretation limits

- Step 03 reads archived CIBERSORTx results and does not rerun deconvolution.
  Archived logs record LM22, quantile normalization and one permutation. The
  exact identity of the archived mixture files with the historical job uploads
  has not been independently established.
- Step 06 reconstructs from archived probe-row manifests; the original
  feature-selection rule has not been recovered. The historical PD
  module–trait calculation treats biopsies as independent despite repeated
  sampling from patients. Its historical hypergeometric background counts
  probe rows although the lists are gene-level. Archived STRING evidence is a
  manual result; the original settings have not been fully recovered.
- Step 07 uses a panel selected before cross-validation. The CV AUCs are
  conditional on this fixed panel and do not independently validate panel
  selection. External OSCC validation is separate; an overlapping dataset
  cannot serve as independent PD validation. Earlier SVM results are not
  reproduced by this wrapper.
- Step 08 records manual cBioPortal/GEPIA2 steps; web settings, exports and
  dates are not fully documented by the text file alone.
- Step 09 retains historical correlations and a separate patient-cluster
  bootstrap sensitivity analysis. Saved references are comparison targets,
  not independent validation, and the script does not recreate the final
  Figure 9 layout. Bootstrap intervals are not adjusted for 65 comparisons.
- Step 10 is exploratory enrichment of the five-gene panel, separate from
  DEG enrichment in Step 04. The historical five-gene GO object was not
  recovered. KEGG and GO results depend on annotation versions.

Revised manuscript text shown by the author reports 18 PD and 729 OSCC GO
terms, as in the recorded Step 04 rerun. The saved final manuscript and Figure
4 have not been verified together. The specific causes of earlier differing
GO counts have not been established. Historical and repeated analyses use
different expression matrices in some steps, so they should not be described
as interchangeable.

## Data and release status

`large_files_manifest.csv` describes three files not included here. Before
giving a reviewer access, arrange an appropriate way to supply any needed
files, verify their hashes, and check whether sample-level inputs, archived
outputs and third-party material can be shared. The included files also
contain sample-level data and derived CIBERSORTx results; inclusion here does
not establish privacy approval or redistribution rights. This draft was
assembled locally and has not been published or submitted.
