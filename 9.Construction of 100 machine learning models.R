# 1. Loading R packages ----
library(mlr3)
library(mlr3learners)
library(mlr3tuning)
library(mlr3pipelines)
library(mlr3extralearners)
library(mlr3viz)
library(mlr3filters)
library(mlr3fselect)
library(caret)
library(pROC)
library(PRROC)
library(dplyr)
library(xgboost)
library(randomForest)
library(e1071)
library(lightgbm)
library(catboost)
library(keras)
library(ada)
library(naivebayes)
library(data.table)
library(ggplot2)
library(Boruta)
ibrary(infotheo)

# 2. Data preparation ----
# training set GSE186143 gene pair matrix
exprset_REO_rank 
dim(exprset_REO_rank)  #1573770  56
matadata
# all gene pairs associated with irAE severity
GPs_Fatal 
dim(GPs_Fatal)   #14284  16
# all irAE severity-associated cell pairs and their corresponding gene pairs
phyper_CC_Pvalue_ssc
phyper_CC_GPs_ssc 
length(phyper_CC_GPs_ssc)  #67

# 3. Initial feature selection ----
C_C_p<-phyper_CC_Pvalue_ssc$Cell_Type_1 %in% c("CD56_dim_CD16_hi_NK_sp")
table(C_C_p) #13 cell-cell pairs
sum(phyper_CC_Pvalue_ssc$Num_FatalGPs[C_C_p]) 
a<-phyper_CC_GPs_ssc[C_C_p]
enrichGP <- unlist(a)
enrichGP <- unique(enrichGP)
mps_enrich <- exprset_REO_rank[enrichGP,4:56]

# 4. Preprocessing of the training and validation sets ----
# Preprocessing of training set
mps_enrich
matadata
data_train <- t(mps_enrich)
colnames(data_train) <- rownames(mps_enrich)
data_train <- as.data.frame(data_train)
data_train[] <- lapply(data_train, as.numeric)
data_train$irAE_severity <- as.factor(matadata$irAE_grade_binary[match(rownames(data_train), matadata$Sample_names)])
data_train <- na.omit(data_train)
dim(data_train) # [53, 1086]

# Preprocessing of validation set 1 GSE186143 gene pair matrix
Test_1_bulk
Test_1_bulk_metadata
data_test <- t(Test_1_bulk[enrichGP,4:ncol(Test_1_bulk)])
data_test <- as.data.frame(data_test)
data_test[] <- lapply(data_test, as.numeric) 
data_test$irAE_severity <- as.factor(Test_1_bulk_metadata$irAE_grade_binary[match(rownames(data_test), Test_1_bulk_metadata$Sample_names)])
data_test <- na.omit(data_test);dim(data_test) #[7,1086]
data_val1<-data_test

# Preprocessing of validation set 2 GSE189125 gene pair matrix
Test_2_bulk
Test_2_bulk_metadata
data_test <- t(Test_2_bulk)
data_test <- as.data.frame(data_test)
data_test[] <- lapply(data_test, as.numeric)
data_test$irAE_severity <- as.factor(Test_2_bulk_metadata$irAE_grade_binary[match(rownames(data_test), Test_2_bulk_metadata$V1)])
data_test <- na.omit(data_test);dim(data_test) #[16 1086]
data_val2.2<-data_test

# Preprocessing of validation set 3 GSE253720 gene pair matrix
Test_3_bulk_metadata
Test_3_bulk
data_test <- t(Test_3_bulk)
data_test <- as.data.frame(data_test)
data_test[] <- lapply(data_test, as.numeric)
data_test$irAE_severity <- as.factor(Test_3_bulk_metadata$irAE_grade_binary[match(rownames(data_test), Test_3_bulk_metadata$V1)])
data_test <- na.omit(data_test);dim(data_test) #[8 1086]
data_val3.3<-data_test

# 5. Create an mlr3 task----
if(TRUE){
mi_scores <- apply(data_train[, !names(data_train) %in% "irAE_severity"], 2, 
                   function(x) mutinformation(x, data_train$irAE_severity))
top_features <- names(sort(mi_scores, decreasing=TRUE))[1:100]
data_train_100 <- data_train[, c(top_features, "irAE_severity")]
task_train <- TaskClassif$new(id="irAE_train", backend=data_train_500, target="irAE_severity", positive="1")
}
colnames(data_train)<-make.names(colnames(data_train))
colnames(data_val1)<-make.names(colnames(data_val1))
colnames(data_val2.2)<-make.names(colnames(data_val2.2))
colnames(data_val3.3)<-make.names(colnames(data_val3.3))

task_train <- TaskClassif$new(id="irAE_train", backend=data_train, target="irAE_severity", positive="1")
task_val1 <- TaskClassif$new(id="irAE_val1", backend=data_val1, target="irAE_severity", positive="1")
task_val2.2 <- TaskClassif$new(id="irAE_val2.2", backend=data_val2.2, target="irAE_severity", positive="1")
task_val3.3 <- TaskClassif$new(id="irAE_val3.3", backend=data_val3.3, target="irAE_severity", positive="1")

# 6. Define the learner ----
mlr_learners$keys(pattern="classif") 
learners <- list(
  lrn("classif.xgboost", id="xgboost", predict_type="prob", eta=0.05, max_depth=3, lambda=1, nrounds=50),
  lrn("classif.ranger", id="random_forest", predict_type="prob", num.trees=100, min.node.size=5),
  lrn("classif.svm", id="svm", predict_type="prob", kernel="radial", type="C-classification", cost=1),
  lrn("classif.log_reg", id="logistic_regression", predict_type="prob"),
  lrn("classif.kknn", id="knn", predict_type="prob", k=5),
  lrn("classif.lda", id="lda", predict_type="prob"),
  lrn("classif.lightgbm", id="lightgbm", predict_type="prob", num_leaves=7, learning_rate=0.05, min_data_in_leaf=3, num_iterations=50),
  lrn("classif.catboost", id="catboost", predict_type="prob", iterations=50, depth=3),
  lrn("classif.nnet", id="neural_network", predict_type="prob", size=3, maxit=100,MaxNWts = 1000000),
  lrn("classif.naive_bayes", id="naive_bayes", predict_type="prob")
)

# 7. Feature selection ----
feature_subsets <- list()
set.seed(123)
for (learner in learners) {
  cat("Processing learner:", learner$id, "\n")
  fselector <- fs("random_search", max_features=15, batch_size = 50) 
  instance <- fselect(
    fselector = fselector,
    task = task_train,
    learner = learner,
    resampling = rsmp("cv", folds = 5)$instantiate(task_train),
    measure = msr("classif.auc"), 
    terminator = trm("combo", terminators = list(
      trm("evals", n_evals = 1000),  
      trm("stagnation", iters = 50)  
    ))
  )
  selected_features <- instance$result_feature_set
  if (length(selected_features) > 20) {
    selected_features <- selected_features[1:20]  
  }
  feature_subsets[[learner$id]] <- selected_features
  cat("Selected features£¨", learner$id, "£©£º", selected_features, "\n")
}

write.csv(data.frame(
  Learner=names(feature_subsets),
  Features=sapply(feature_subsets, paste, collapse=", ")
), "feature_subsets2.csv", row.names=FALSE)

# 8. Train the model ----
models <- list()
set.seed(123)
for (learner in learners) {
  learner_id <- learner$id
  for (feature_learner_id in names(feature_subsets)) {
    selected_features <- feature_subsets[[feature_learner_id]]
    model_id <- paste0(learner_id, "_features_", feature_learner_id)
    task_subset <- TaskClassif$new(
      id=paste0("train_", model_id),
      backend=data_train[, c(selected_features, "irAE_severity")],
      target="irAE_severity",
      positive="1"      
    )
    learner_clone <- learner$clone(deep=TRUE) 
    learner_clone$train(task_subset)
    models[[model_id]] <- learner_clone
  }
}

# 9. Evaluation of model predictive performance in the validation set ----
validation_results <- list()
validation_tasks<-list(task_train,task_val1,task_val2.2,task_val3.3)
set.seed(123)
for (val_task in validation_tasks) {
  val_id <- val_task$id
  cat("Evaluating on validation set:", val_id, "\n")
  results <- data.frame(Model=character(), AUC=numeric(), F1=numeric(), Recall=numeric(), stringsAsFactors=FALSE)
  for (learner in learners) {
    learner_id <- learner$id
    for (feature_learner_id in names(feature_subsets)) {
      model_id <- paste0(learner_id, "_features_", feature_learner_id)
      selected_features <- feature_subsets[[feature_learner_id]]
      val_data <- as.data.frame(val_task$data())
      missing_features <- setdiff(selected_features, names(val_data))
      if (length(missing_features) > 0) {
        val_data[missing_features] <- 0
      }
      val_task_subset <- TaskClassif$new(
        id=paste0(val_id, "_", model_id),
        backend=val_data[, c(selected_features, "irAE_severity")],
        target="irAE_severity",
        positive="1"
      )
      pred <- models[[model_id]]$predict(val_task_subset)
      auc_score <- tryCatch({
        roc(val_data$irAE_severity, pred$prob[, "1"],levels = c("0", "1"),direction = "<")$auc
      }, error=function(e) NA)
      pred_class <- ifelse(pred$prob[, "1"] > 0.5, 1, 0)
      cm <- tryCatch({
        confusionMatrix(factor(pred_class), factor(val_data$irAE_severity), positive="1")
      }, error=function(e) list(byClass=c(F1=NA, Recall=NA)))
      f1_score <- cm$byClass["F1"]
      recall <- cm$byClass["Recall"]
      
      results <- rbind(results, data.frame(
        Model=model_id,
        AUC=auc_score,
        F1=f1_score,
        Recall=recall
      ))
    }
  }
  validation_results[[val_id]] <- results
  write.csv(results, paste0("validation_results_", val_id, ".csv"), row.names=FALSE)
}
saveRDS(validation_results,"5.validation_results")

# 10.Summary of predictive performance results from 100 models ----
learners_names <- c("xgboost", "random_forest", "svm", "logistic_regression", "knn", "lda", "lightgbm", "catboost", "neural_network", "naive_bayes")
feature_names <- names(feature_subsets)
model_ids <- paste0(
  rep(learners_names, each=10), "_features_", rep(feature_names, times=10)
)
validation_names <- names(validation_results) 
auc_matrix <- matrix(NA, nrow=100, ncol=length(validation_names), dimnames=list(model_ids, validation_names))
for (val_id in validation_names) {
  results <- validation_results[[val_id]]
  if (any(is.na(results$AUC))) {
    cat("Warning: NA values found in AUC for", val_id, "\n")
    results$AUC[is.na(results$AUC)] <- mean(results$AUC, na.rm=TRUE) 
  }
  for (i in 1:nrow(results)) {
    model_id <- results$Model[i]
    auc_matrix[model_id, val_id] <- results$AUC[i]
  }
}
saveRDS(auc_matrix,"6.auc_matrix")
write.csv(auc_matrix, "6.auc_matrix.csv", row.names=TRUE)


