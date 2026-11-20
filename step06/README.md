# Step 06: historical WGCNA reconstruction and STRING evidence

This draft preserves the historical recipe for audit. Run from this `step06` directory. It includes saved inputs from steps 01 and 02, archived feature and sample manifests, and selected records from an earlier completed reconstruction. It does **not** establish the complete original preprocessing and feature-selection provenance. Review data-sharing rights and metadata before public release. The approximately 195 MB preprocessing RDS must first be reconstructed from the included binary parts: from the repository root, follow the [restore instructions](../RESTORE_LARGE_FILES.md). Restoration verifies file size and SHA-256 and was tested locally, not on another computer.

## Contents and prerequisites

- `scripts/06_WGCNA_PPI_hubgenes.R`: input checks and optional historical WGCNA rebuild.
- `results/01_preprocessing/`: PD/OSCC sample maps and, after restoration, saved input ExpressionSets.
- `results/02_DEG/`: saved differential-expression object and summary, with `historical_unblocked` PD model.
- `data/reference/`: archived feature order, probe-to-gene labels, sample order, and numeric fingerprints.
- `audit_reference/`: selected results from the completed 18 November 2026 run, package versions, and the recovered STRING export and image.

The scripts require R with Biobase and digest for checking the inputs, plus WGCNA for a rebuild. The prior reconstruction used R 4.6.1 and WGCNA 1.74 on Windows; see `audit_reference/06_sessionInfo.txt` for the recorded environment. Software versions and hardware can affect recomputation. The included input files are snapshots; use the prior steps and source datasets to investigate earlier processing, not this folder alone.

## Checking and rebuilding

Open R with the working directory set to this `step06` folder. The default check reads inputs, checks cohort order, groups, PD patient mapping, complete saved DEG lists, historical row manifests and numeric SHA-256 fingerprints; it builds no networks and writes no outputs. If set in the environment, `PD_OSCC_PROJECT_ROOT` and `PD_OSCC_06_REFERENCE_DIR` override these default paths and must point to the intended input locations.

```r
tryCatch({
  Sys.setenv(PD_OSCC_06_MODE = "check_inputs")
  source("scripts/06_WGCNA_PPI_hubgenes.R")
}, finally = Sys.unsetenv("PD_OSCC_06_MODE"))
```

The optional `rebuild` mode repeats lengthy network calculations and writes a fresh timestamped directory under `results/06_WGCNA/historical_reconstruction/`. It can require substantial time and memory. It does not overwrite the archived audit reference.

```r
tryCatch({
  Sys.setenv(PD_OSCC_06_MODE = "rebuild")
  source("scripts/06_WGCNA_PPI_hubgenes.R")
}, finally = Sys.unsetenv("PD_OSCC_06_MODE"))
```

## Evidence and limits

The historical input matrices comprise 28,655 PD and 28,704 OSCC **probe rows with gene labels**, including duplicate labels. The original rule selecting those rows was not recovered. The archived manifests reproduce the saved input rows but do not reconstruct that selection rule. No top-5,000 gene-selection rule is asserted here.

The saved reconstruction reported 29 and 64 non-grey modules, 32 and 662 DEG/module intersections, and the same five shared candidates as the archived panel. At fixed PD beta 6 and OSCC beta 3, recalculated scale-free fit R² values were 0.8255 and 0.8758 respectively. These values refer to the historical probe-row matrices. The separate PD analysis of 20,848 unique genes gave R² 0.7797972 at beta 6 and approximately 0.817 at beta 7; these are different inputs and must not be presented as the same diagnosis. See the two comparison CSVs in `audit_reference/`.

The PD module–trait p-values use 247 biopsies as independent observations although they come from 90 patients. The historical hypergeometric calculation uses a probe-row background for gene-level lists; its arithmetic reproduction does not justify that background. Computational agreement is not validation of those assumptions.

The STRING TSV and image are archived evidence of a manual step, not output created by this script. The TSV has four undirected edges represented in both directions; the image also shows CHST2 as an isolated node. The original STRING version, date, and exact query settings have not been independently established. The five shared candidates are not necessarily network-centrality hubs merely because they intersect. This code does not rerun STRING.

The earlier network computation preceded subsequent documentation-only changes to the copy of this script. On 20 November 2026, the standalone draft's `check_inputs` mode completed successfully on the original Windows computer with R 4.6.1 and the installed Biobase and digest packages: both archived matrix fingerprints matched, and the PD/OSCC checks reported 247/212 samples and 57/1,266 complete DEG entries. A before/after comparison of file paths, sizes and modification times found no changed, new or deleted files in this folder. This check did not rebuild networks. Distribution and execution on another machine have not yet been tested.
