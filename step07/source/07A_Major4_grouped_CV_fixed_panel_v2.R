# MAJOR COMMENT 4 — STEP 07A
# Group-aware 10-fold out-of-fold audit for the PRE-SPECIFIED original five-gene panel.
# This is a corrected CROSS-VALIDATION OF A FIXED PANEL, NOT nested feature selection.
# It does NOT demonstrate that full-data DEG/WGCNA feature discovery was leakage-free.
# It does NOT estimate external validation AUC (GSE16134/GSE37991 not needed here).
# It DOES NOT modify originals; creates new files only under results/07_Major4_fixed_panel_CV.
# Requires R packages: caret, pROC.
# V2: uses the existing Step 01D sample maps (same as successful Major3 audit);
# does NOT require the optional 01B APPROVED CSV copies in data/metadata.

project_dir <- .script07_data_root
if (!dir.exists(project_dir)) stop("Project folder missing: ", project_dir)
if (!requireNamespace("caret", quietly = TRUE) || !requireNamespace("pROC", quietly = TRUE)) {
  stop("Install packages first: install.packages(c('caret','pROC'))")
}

out_dir <- file.path(project_dir, "results", "07_Major4_fixed_panel_CV")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

panel <- c("CXCL1", "CXCL6", "CXCR4", "RGS1", "CHST2")

read_dataset <- function(gse, disease_name) {
  mat_path <- file.path(project_dir, "data", "processed", paste0(gse, "_01D_log2_gene_matrix.rds"))
  lab_path <- file.path(project_dir, "data", "processed", paste0(gse, "_01D_sample_map.csv"))
  if (!file.exists(mat_path)) stop("Missing processed 01D matrix: ", mat_path)
  if (!file.exists(lab_path)) stop("Missing Step 01D sample map: ", lab_path)
  mat <- readRDS(mat_path)
  lab <- read.csv(lab_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!is.matrix(mat) || is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop(gse, " 01D object is not a gene x GSM matrix with dimnames")
  }
  if (!all(panel %in% rownames(mat))) {
    stop(gse, " missing genes: ", paste(setdiff(panel, rownames(mat)), collapse=", "))
  }
  if (anyDuplicated(rownames(mat))) stop(gse, " duplicate gene rows")
  if (anyDuplicated(colnames(mat))) stop(gse, " duplicate GSM columns")
  req <- c("GSM", "Group")
  if (!all(req %in% names(lab)) || anyDuplicated(lab$GSM)) {
    stop(gse, " invalid/duplicate Step 01D sample map")
  }
  # Align BY GSM, never by row ordering.
  idx <- match(colnames(mat), lab$GSM)
  if (anyNA(idx)) stop(gse, " missing Step 01D labels for some expression samples")
  lab <- lab[idx, , drop=FALSE]
  if (!identical(as.character(lab$GSM), colnames(mat))) stop(gse, " alignment failed")
  if (anyNA(lab$Group)) stop(gse, " has missing group labels")
  observed <- as.character(lab$Group)
  if (!all(observed %in% c(disease_name, "Control"))) {
    stop(gse, " unexpected class labels; check dysplasia removal: ",
         paste(unique(observed), collapse=", "))
  }
  x <- t(mat[panel, , drop=FALSE])
  if (any(!is.finite(x))) stop(gse, " nonfinite gene expression values")
  df <- as.data.frame(x, check.names=FALSE)
  df$y <- factor(ifelse(observed == disease_name,"Disease","Control"),
                 levels=c("Disease","Control"))
  rownames(df) <- colnames(mat)
  if (gse == "GSE10334") {
    if (!"Patient_ID" %in% names(lab) || anyNA(lab$Patient_ID) ||
        any(!nzchar(as.character(lab$Patient_ID)))) {
      stop("GSE10334 missing valid Patient_ID; cannot do grouped CV safely")
    }
    patient <- as.character(lab$Patient_ID)
    if (length(unique(patient)) != 90L) stop("PD expected 90 patients, found ",length(unique(patient)))
    if (nrow(df)!=247L || sum(df$y=="Disease")!=183L || sum(df$y=="Control")!=64L) {
      stop("GSE10334 sample/group counts unexpected")
    }
  } else {
    patient <- rownames(df)  # sample-wise CV; independent people are not verified here
    if (nrow(df)!=212L || sum(df$y=="Disease")!=167L || sum(df$y=="Control")!=45L) {
      stop("GSE30784 sample/group counts unexpected")
    }
  }
  cat(gse, ": ", nrow(df), " samples (",sum(df$y=="Disease")," disease; ",
      sum(df$y=="Control")," control); ",length(unique(patient)),
      " CV grouping IDs\n",sep="")
  list(df=df,patient=patient)
}

pd <- read_dataset("GSE10334", "PD")
oscc <- read_dataset("GSE30784", "OSCC")

# Produce patient-disjoint 10-fold CV, retrying seeds until ALL held-out folds contain both classes.
make_pd_indices <- function(y, patient, k=10) {
  for (seed in 42:242) {
    set.seed(seed)
    idx <- caret::groupKFold(patient, k=k) # training indices; holds out entire patients
    outs <- lapply(idx, function(tr) setdiff(seq_along(y), tr))
    ok <- all(vapply(seq_along(outs), function(i) {
      te <- outs[[i]]; tr <- idx[[i]]
      length(unique(as.character(y[te])))==2L &&
        length(unique(as.character(y[tr])))==2L &&
        !length(intersect(unique(patient[tr]), unique(patient[te])))
    }, logical(1)))
    if (ok) {
      names(idx) <- paste0("Fold", seq_along(idx))
      names(outs) <- names(idx)
      return(list(index=idx,indexOut=outs,seed=seed))
    }
  }
  stop("Could not make 10 grouped folds with both labels; STOP rather than leak")
}

make_oscc_indices <- function(y,k=10) {
  set.seed(42)
  idx <- caret::createFolds(y, k=k, returnTrain=TRUE)
  outs <- lapply(idx,function(tr) setdiff(seq_along(y),tr))
  if (!all(vapply(outs, function(te) length(unique(as.character(y[te])))==2L, logical(1)))) {
    stop("OSCC has a held-out fold with one class")
  }
  list(index=idx,indexOut=outs,seed=42)
}

# Within-fold scaling: caret fits center/scale parameters using ONLY training rows.
# If using label-selected five genes, resulting AUC remains conditional on prior full-data selection.
fit_cv <- function(dat,label,folds) {
  ctrl <- caret::trainControl(method="cv", number=length(folds$index),
     index=folds$index, indexOut=folds$indexOut,
     classProbs=TRUE, summaryFunction=caret::twoClassSummary,
     savePredictions="final", allowParallel=FALSE)
  set.seed(20261023)
  model <- suppressWarnings(caret::train(y ~ ., data=dat$df,
     method="glm", family=binomial(), trControl=ctrl, metric="ROC",
     preProcess=c("center","scale")))
  pred <- model$pred
  if (is.null(pred) || nrow(pred)!=nrow(dat$df)) {
    stop(label, " prediction count invalid: expected exactly one held-out prediction per sample")
  }
  if (anyDuplicated(pred$rowIndex) || !setequal(pred$rowIndex,seq_len(nrow(dat$df)))) {
    stop(label, " OOF prediction indices duplicated/missing")
  }
  pred$GSM <- rownames(dat$df)[pred$rowIndex]
  pred$Patient_ID <- dat$patient[pred$rowIndex]
  pred$True_label <- as.character(dat$df$y[pred$rowIndex])
  if (!identical(as.character(pred$obs), pred$True_label)) stop(label," OOF labels misaligned")
  roc <- pROC::roc(response=pred$obs, predictor=pred$Disease,
     levels=c("Control","Disease"),direction="<",quiet=TRUE)
  pooled_auc <- as.numeric(pROC::auc(roc))
  # Bootstrap ONLY the fixed held-out predictions: uncertainty is conditional on the fitted CV.
  # PD group resampling by Patient_ID; OSCC by individual sample (identity).
  unique_ids <- unique(pred$Patient_ID)
  B <- 1000L
  set.seed(1107)
  boot <- replicate(B, {
    draw <- sample(unique_ids,length(unique_ids),replace=TRUE)
    ii <- unlist(lapply(draw, function(id) which(pred$Patient_ID==id)),use.names=FALSE)
    if (length(unique(pred$obs[ii]))<2L) return(NA_real_)
    as.numeric(pROC::auc(pROC::roc(response=pred$obs[ii],predictor=pred$Disease[ii],
       levels=c("Control","Disease"),direction="<",quiet=TRUE)))
  })
  boot <- boot[is.finite(boot)]
  if (length(boot)<950) stop(label, " excessive degenerate bootstrap samples")
  ci <- as.numeric(stats::quantile(boot,c(0.025,0.975),names=FALSE))
  # Save new audit ONLY, no overwrite of previous stage 07.
  new_pred_path <- file.path(out_dir,paste0(label,"_OOF_predictions.csv"))
  if(file.exists(new_pred_path)) stop("Output already exists; move/rename before rerun: ",new_pred_path)
  write.csv(pred[,c("GSM","Patient_ID","True_label","Disease","Control","Resample","rowIndex")],
            new_pred_path,row.names=FALSE)
  fold_map <- do.call(rbind,lapply(seq_along(folds$indexOut),function(i){
    te <- folds$indexOut[[i]]
    data.frame(GSM=rownames(dat$df)[te],Patient_ID=dat$patient[te],
               Fold=names(folds$indexOut)[i],stringsAsFactors=FALSE)
  }))
  write.csv(fold_map,file.path(out_dir,paste0(label,"_folds.csv")),row.names=FALSE)
  report <- data.frame(
    Dataset=label, Samples=nrow(dat$df), CV_grouping_units=length(unique(dat$patient)),
    Panel=paste(panel,collapse=";"), Model="glm logistic regression",
    CV_Folds=10, Fold_seed=folds$seed, Preprocessing="center/scale WITHIN folds",
    Mean_fold_AUC=model$results$ROC, SD_fold_AUC=model$results$ROCSD,
    Pooled_OOF_AUC=pooled_auc, OOF_AUC_bootstrap_95_low=ci[1],
    OOF_AUC_bootstrap_95_high=ci[2], Bootstrap_n=length(boot),
    Patient_grouped_CV=(label=="PD"),
    Caveat="Panel selected using complete discovery datasets; CV AUC remains potentially optimistic. Bootstrap CI conditional on OOF predictions.",
    stringsAsFactors=FALSE)
  cat("\n",label," grouped/stratified 10-fold CV (FIXED panel):\n",sep="")
  print(report[,c("Mean_fold_AUC","SD_fold_AUC","Pooled_OOF_AUC",
                 "OOF_AUC_bootstrap_95_low","OOF_AUC_bootstrap_95_high")],row.names=FALSE)
  report
}

fold_pd <- make_pd_indices(pd$df$y,pd$patient)
fold_oscc <- make_oscc_indices(oscc$df$y)

# Output guarded so old files cannot silently be overwritten.
summary_path <- file.path(out_dir,"Major4_fixed_panel_corrected_CV_summary.csv")
if(file.exists(summary_path)) stop("Existing output: ",summary_path,
      "; move/rename it before another run")
res_pd <- fit_cv(pd,"PD",fold_pd)
res_oscc <- fit_cv(oscc,"OSCC",fold_oscc)
write.csv(rbind(res_pd,res_oscc),summary_path,row.names=FALSE)
writeLines(c(
  paste("Time:",Sys.time()),
  paste("Project:",project_dir),
  "NOTE: This is correct cross-validation for a FIXED, preselected panel.",
  "Feature discovery on full samples before CV is not corrected by this step.",
  "GSE37991 and GSE16134 external AUC not computed in this step.",
  capture.output(sessionInfo())),
  file.path(out_dir,"Major4_session_info.txt"))
cat("\nDONE. Send this summary to reviewer analysis:\n",summary_path,"\n",sep="")
