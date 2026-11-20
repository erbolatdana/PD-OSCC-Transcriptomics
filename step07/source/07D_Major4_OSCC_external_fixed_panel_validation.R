# ============================================================================
# PD–OSCC REANALYSIS | STEP 07D — OSCC EXTERNAL VALIDATION AUDIT
# Strictly fixed historical 5-gene panel: CXCL1, CXCL6, CXCR4, RGS1, CHST2.
# Train ONLY on GSE30784 Step 01D; evaluate GSE37991 ONLY once.
# No external-label-based probe choice, feature selection, fitting, threshold,
# normalization, hyperparameter tuning, or ROC-direction optimization.
# GPL6883 (Illumina) differs from GPL570 (Affymetrix). Direct coefficient
# transport is an exploratory cross-platform test, NOT clinical validation.
# For CXCR4's 3 GPL6883 probes, use the unweighted mean on log2 scale;
# specified BEFORE external labels are considered.
# AUC CI uses 40 matched patient-pair bootstrap (not independent samples).
# Permutation swaps OSCC/control labels WITHIN each patient pair.
# The 5 genes were selected from complete discovery datasets: neither this
# external check nor the earlier fixed-panel CV removes discovery-stage leakage.
# Run after steps 07A, 07B and 07C; NO rerun of GEO/DEG/WGCNA/CIBERSORT.
# ============================================================================

project_dir <- .script07_data_root
if (!dir.exists(project_dir)) stop("Project folder is missing: ", project_dir)
for (pkg in c("Biobase", "pROC")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install package: ", pkg)
}

panel <- c("CXCL1", "CXCL6", "CXCR4", "RGS1", "CHST2")
B <- 2000L
out <- file.path(project_dir, "results", "07_Major4_external_validation", "03_OSCC_fixed_panel_external")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

in_train <- file.path(project_dir, "data", "processed", "GSE30784_01D_log2_gene_matrix.rds")
map_train <- file.path(project_dir, "data", "processed", "GSE30784_01D_sample_map.csv")
orig_ext <- file.path(project_dir, "data", "external_validation", "original_GEO", "GSE37991_ORIGINAL_GEO_ExpressionSet.rds")
map_ext <- file.path(project_dir, "data", "external_validation", "original_GEO", "GSE37991_VERIFIED_sample_map.csv")
ann_ext <- file.path(project_dir, "results", "07_Major4_external_validation", "02_mapping_overlap_QC", "GSE37991_five_gene_probe_annotation.csv")
inputs <- c(in_train, map_train, orig_ext, map_ext, ann_ext)
if (!all(file.exists(inputs))) stop("Missing required input files:\n", paste(inputs[!file.exists(inputs)], collapse="\n"))

output_main <- file.path(out, "Major4_OSCC_GSE37991_external_AUC_summary.csv")
if (file.exists(output_main)) stop("A completed summary already exists. Do not silently overwrite:\n", output_main)

log_path <- file.path(out, "Major4_07D_validation_log.txt")
writeLines(c("STEP 07D — historical 5-gene panel external OSCC validation",
             paste("Time:", Sys.time()), paste("R:", R.version.string)), log_path)
logmsg <- function(...) {
  msg <- paste0(..., collapse="")
  cat(msg, "\n", sep="")
  cat(msg, "\n", file=log_path, append=TRUE, sep="")
}

# ---- Train data from frozen Step 01D (GSE30784, GPL570)
tr <- readRDS(in_train)
mt <- utils::read.csv(map_train, check.names=FALSE, stringsAsFactors=FALSE)
if (!is.matrix(tr) || anyDuplicated(rownames(tr)) || anyDuplicated(colnames(tr))) stop("Invalid training gene expression matrix")
if (!all(panel %in% rownames(tr))) stop("Missing training genes: ",paste(setdiff(panel,rownames(tr)),collapse=", "))
if (!all(c("GSM","Group") %in% names(mt)) || anyDuplicated(mt$GSM)) stop("Invalid training labels")
idx <- match(colnames(tr), mt$GSM)
if (anyNA(idx)) stop("GSM mismatch in discovery training labels")
mt <- mt[idx,,drop=FALSE]
if (!all(mt$Group %in% c("OSCC","Control")) || ncol(tr)!=212L || sum(mt$Group=="OSCC")!=167L || sum(mt$Group=="Control")!=45L) {
  stop("Unexpected GSE30784 sample groups/counts")
}
xtr <- t(tr[panel,,drop=FALSE]); storage.mode(xtr) <- "double"
if (any(!is.finite(xtr))) stop("Non-finite training data")
ytr <- as.integer(mt$Group=="OSCC")

# ---- External data: 40 OSCC + 40 non-OSCC tissues, matched by person
es <- readRDS(orig_ext)
if (!methods::is(es,"ExpressionSet")) stop("External RDS is not an ExpressionSet")
ext <- Biobase::exprs(es)
me <- utils::read.csv(map_ext,check.names=FALSE,stringsAsFactors=FALSE)
if (anyDuplicated(colnames(ext)) || anyDuplicated(rownames(ext))) stop("Duplicate GSM/probe IDs in external series")
if (!all(c("GSM","Group","Patient_ID") %in% names(me)) || anyDuplicated(me$GSM)) stop("External sample map invalid")
idx <- match(colnames(ext),me$GSM)
if(anyNA(idx))stop("External GSM mismatch")
me <- me[idx,,drop=FALSE]
if (ncol(ext)!=80L || sum(me$Group=="OSCC")!=40L || sum(me$Group=="Control")!=40L || length(unique(me$Patient_ID))!=40L) {
  stop("Unexpected external cohort sample counts")
}
paired <- table(me$Patient_ID,me$Group)
if(!identical(sort(colnames(paired)),sort(c("OSCC","Control"))) || any(paired!=1L)) stop("External pair design invalid")

# ---- Platform-native signal -> log2 only if clearly on large linear scale.
# Rule is based on marginal expression distribution ONLY, never labels or AUC.
# Do not quantile-normalize platforms together or center using external labels.
p99_before <- as.numeric(stats::quantile(as.numeric(ext),0.99,na.rm=TRUE))
if (!is.finite(p99_before)) stop("Bad external expression distribution")
transform_rule <- if (p99_before > 100) "log2(x + 1) applied because external marginal P99 > 100" else "as deposited (already low-scale / likely log2)"
if (p99_before>100) {
  if (any(ext < 0,na.rm=TRUE)) stop("External has negative raw values; automatic log2 transform is unsafe")
  ext <- log2(ext + 1)
}
logmsg("GSE37991 expression transformation: ", transform_rule,"; original P99=", round(p99_before,3))

# ---- Fixed external probe mapping without inspecting clinical group differences.
ann <- utils::read.csv(ann_ext,check.names=FALSE,stringsAsFactors=FALSE)
if (!all(c("Probe_ID","Gene_Symbol","Present_in_GSE37991") %in% names(ann))) stop("Annotation columns missing")
ann$Gene_Symbol <- toupper(trimws(as.character(ann$Gene_Symbol)))
ann$Probe_ID <- trimws(as.character(ann$Probe_ID))
ann <- ann[ann$Gene_Symbol %in% panel & ann$Probe_ID %in% rownames(ext),,drop=FALSE]
ann <- unique(ann[,c("Probe_ID","Gene_Symbol"),drop=FALSE])
if (anyDuplicated(ann$Probe_ID)) stop("Ambiguous probe-to-gene mapping")
ct <- table(factor(ann$Gene_Symbol,levels=panel))
if (!identical(as.integer(ct),c(1L,1L,3L,1L,1L))) stop("Probe counts changed. Expected CXCL1=1,CXCL6=1,CXCR4=3,RGS1=1,CHST2=1; got ",paste(ct,collapse=","))
xext <- vapply(panel,function(g) {
  probes <- ann$Probe_ID[ann$Gene_Symbol==g]
  colMeans(ext[probes,,drop=FALSE])
},numeric(ncol(ext)))
rownames(xext) <- colnames(ext)
if (!identical(colnames(xext),panel) || any(!is.finite(xext))) stop("Bad external mapped panel")

# ---- Training-only preprocessing, fit ONCE on full GSE30784.
center <- colMeans(xtr); sdev <- apply(xtr,2,stats::sd)
if (any(!is.finite(sdev)) || any(sdev<=0)) stop("Training expression SD invalid")
ztr <- sweep(sweep(xtr,2,center,"-"),2,sdev,"/")
zext <- sweep(sweep(xext,2,center,"-"),2,sdev,"/")
train_df <- as.data.frame(ztr); train_df$y <- ytr
ext_df <- as.data.frame(zext)
model <- suppressWarnings(stats::glm(y ~ ., data=train_df,family=stats::binomial()))
if(any(!is.finite(stats::coef(model))))stop("Full-training GLM coefficients nonfinite")
prob <- as.numeric(stats::predict(model,newdata=ext_df,type="response"))
if (length(prob)!=80L || any(!is.finite(prob))) stop("Model predictions missing/nonfinite; inspect cross-platform shift")

# ---- Outcome labels used ONLY for evaluation, not for train/test transforms.
y <- as.character(me$Group)
auc_fun <- function(labels, probabilities) {
  as.numeric(pROC::auc(pROC::roc(response=factor(labels,levels=c("Control","OSCC")),
                        predictor=probabilities,levels=c("Control","OSCC"),direction="<",quiet=TRUE)))
}
obs_auc <- auc_fun(y,prob)
pairs <- split(seq_along(y),as.character(me$Patient_ID))
if (length(pairs)!=40L || any(lengths(pairs)!=2L)) stop("Pair index invalid")

# ---- Paired patient-level bootstrap of fixed external predictions.
set.seed(20261023)
boot <- replicate(B, {
  chosen <- sample(names(pairs),length(pairs),replace=TRUE)
  ii <- unlist(pairs[chosen],use.names=FALSE)
  auc_fun(y[ii],prob[ii])
})
if (any(!is.finite(boot))) stop("Bootstrap AUC failed")
ci <- unname(stats::quantile(boot,c(.025,.975)))

# ---- Paired exchangeability permutation test: hold model/scores fixed.
set.seed(20261024)
perm <- replicate(B, {
  yy <- y
  for (ii in pairs) if (runif(1)<0.5) yy[ii] <- rev(yy[ii])
  auc_fun(yy,prob)
})
p_one_sided <- (sum(perm>=obs_auc)+1)/(B+1)

# ---- Paired discrimination: scores higher in tumor for matched person.
pair_score_difference <- vapply(pairs,function(ii) {
  prob[ii[y[ii]=="OSCC"]] - prob[ii[y[ii]=="Control"]]
},numeric(1))

# Distribution audit; large shift calls for careful interpretation.
feature_qc <- data.frame(
  Gene=panel, Training_mean=as.numeric(center),Training_SD=as.numeric(sdev),
  External_mean=as.numeric(colMeans(xext)),External_SD=as.numeric(apply(xext,2,stats::sd)),
  External_mean_in_training_SD=as.numeric(colMeans(zext)),
  External_min_Z=as.numeric(apply(zext,2,min)),External_max_Z=as.numeric(apply(zext,2,max)),
  External_probes=as.integer(ct),stringsAsFactors=FALSE)
utils::write.csv(feature_qc,file.path(out,"Major4_GSE37991_cross_platform_feature_QC.csv"),row.names=FALSE)
utils::write.csv(data.frame(GSM=me$GSM,Patient_ID=me$Patient_ID,Group=me$Group,
    Predicted_OSCC_probability=prob,stringsAsFactors=FALSE),
    file.path(out,"Major4_GSE37991_fixed_panel_predictions.csv"),row.names=FALSE)
utils::write.csv(data.frame(Patient_ID=names(pair_score_difference),
     Tumor_minus_control_probability=as.numeric(pair_score_difference)),
    file.path(out,"Major4_GSE37991_within_pair_scores.csv"),row.names=FALSE)

result <- data.frame(
  Validation="GSE37991",Training="GSE30784",Panel=paste(panel,collapse=";"),
  Method="glm binomial; train-only center/scale",N_external_samples=80L,
  N_external_patient_pairs=40L,AUC=obs_auc,Paired_bootstrap_CI95_low=ci[1],
  Paired_bootstrap_CI95_high=ci[2],Bootstrap_reps=B,
  Paired_permutation_one_sided_P=p_one_sided,Permutation_reps=B,
  Proportion_pairs_tumor_higher=mean(pair_score_difference>0),
  External_transformation=transform_rule,CXCR4_probe_handling="mean of 3 log2 probes",
  Discovery_caveat="Full-data feature selection preceded original CV; OSCC external evaluation is exploratory and cross-platform",
  stringsAsFactors=FALSE)
utils::write.csv(result,output_main,row.names=FALSE)

pROC::roc(response=factor(y,levels=c("Control","OSCC")),predictor=prob,
          levels=c("Control","OSCC"),direction="<",quiet=TRUE)->roc
png(file.path(out,"Major4_GSE37991_external_ROC.png"),width=1900,height=1800,res=240)
plot(roc,main=sprintf("OSCC external GSE37991 | AUC=%.3f",obs_auc),col="black",lwd=2)
abline(a=0,b=1,lty=2,col="gray")
dev.off()

logmsg("Training GSE30784: ",nrow(xtr)," samples, 5 fixed genes")
logmsg("External GSE37991: ",nrow(xext)," samples / ",length(pairs)," pairs")
logmsg("CXCR4 probe handling: average of 3 pre-specified exact symbol probes on log2 scale")
logmsg("AUC = ",sprintf("%.4f",obs_auc),"; patient-paired 95% bootstrap CI [",sprintf("%.4f",ci[1]),", ",sprintf("%.4f",ci[2]),"]")
logmsg("Paired label-swap permutation one-sided P = ",sprintf("%.5f",p_one_sided))
logmsg("Tumor score higher than matched control in ",sum(pair_score_difference>0),"/",length(pairs)," pairs")
logmsg("External feature means in training SD units: ",paste(paste(panel,round(colMeans(zext),2),sep="="),collapse="; "))
logmsg("WARNING: AUC does NOT validate biological progression, clinical screening, or unbiased selection. Cross-platform data may have major calibration/scale shifts.")
logmsg("WARNING: This analysis MAY NOT reproduce original AUC 0.726 because original fitting, probe collapse and processing may have differed. Do not tune steps to match that number.")
logmsg("Output summary: ",output_main)
utils::write.csv(data.frame(Input=inputs,MD5=unname(tools::md5sum(inputs))),
                 file.path(out,"Major4_07D_input_manifest.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"Major4_07D_sessionInfo.txt"))
cat("\nDONE. Please send TWO CSV files:\n",output_main,"\n",
    file.path(out,"Major4_GSE37991_cross_platform_feature_QC.csv"),"\n",sep="")
