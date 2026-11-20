Historical input manifests: provenance and limits

These manifests inventory the actual saved matrices archive_env$mat_pd and archive_env$mat_oscc in R_workspace_23_10_2026.RData.
They record the original GPL570 probe ID, first annotated gene symbol, and exact row order. Repeated gene symbols are preserved.
The original filtering rule that produced these saved matrix subsets has NOT been recovered. These manifests do not claim that a first-symbol-only algorithm reconstructs the historical matrices.
Applying only first Gene.symbol before ///, trim, and removal of missing/empty names to either the current or archived ExpressionSets retains 45118 rows per cohort; the archived matrices instead contain 28655 PD and 28704 OSCC rows.
Each archived matrix row was mapped to exactly one input probe using the first four expression values, and then ALL expression values, sample order, and gene labels were required to be strictly identical. No ambiguous matches were accepted.
The current and archived ExpressionSet expression matrices and Gene.symbol vectors were strictly identical.
Network_row denotes the saved input matrix row before any subsequent WGCNA sample/gene quality filtering; it is not automatically a position in the final saved colors vector.
The saved PD network colors have 28655 entries with names identical to the saved PD matrix row names.
The saved OSCC network colors have 28703 entries. Their names exactly match saved OSCC matrix row names after removal of row 97 (DEFB106B; probe 1552412_a_at).
The reason for the removal of OSCC row 97 is not preserved in net_oscc$goodGenes: that vector already has length 28703 and is entirely TRUE. The omitted row has finite and varying values, so it is not an all-missing or exactly constant row.
The manifests and QC were derived only by reading and comparing saved objects. WGCNA, limma, module-trait correlations, and biological gene-list intersections were not rerun.
File basenames, object references, and MD5 hashes are recorded in the QC CSV. Original source files were not modified.
