# 1. Data preparation ----
TrainData_TPM  #training set GSE186143 bulk transcriptome data
TrainData_metadata

Test_1_bulk_TPM #validation set 1 GSE186143 bulk transcriptome data
Test_1_bulk_metadata

Test_2.2_bulk_count #validation set 2 GSE189125 single-cell¨Cderived pseudo-bulk transcriptome data
Test_2.2_bulk_metadata

Test_3.3_bulk_count #validation set 3 GSE253720 Single-cell¨Cderived pseudo-bulk transcriptome data
Test_3.3_bulk_metadata

# 2. Evaluation of the predictive performance of five biomarkers in the training set ----
## PMID£º33009409
# Formula£º0.37¡ÁLCP1+0.70¡ÁADPGK¨C9.10
c("LCP1","ADPGK") %in%  rownames(TrainData_TPM)
lcp1_expr <- as.numeric(TrainData_TPM["LCP1", ])
adpgk_expr <- as.numeric(TrainData_TPM["ADPGK", ])
composite_score <- 0.37 * lcp1_expr + 0.70 * adpgk_expr - 9.10
roc_final <- roc(TrainData_metadata$irAE_grade_binary, composite_score, quiet = TRUE)

## PMID£º31611368
# Formula£ºCD74 and GNAL expression
c("CD74","GNAL") %in%  rownames(TrainData_TPM)
CD74_expr <- as.numeric(TrainData_TPM["CD74", ])
GNAL_expr <- as.numeric(TrainData_TPM["GNAL", ])
roc_final <- roc(TrainData_metadata$irAE_grade_binary, CD74_expr, quiet = TRUE)
roc(TrainData_metadata$irAE_grade_binary, GNAL_expr, quiet = TRUE)

## PMID£º30409824
# Formula£ºGeometric mean of 11 cytokine genes
cytokine_genes <- c("CSF3", "CSF2", "CX3CL1", "FGF2", "IFNA2", "IL12A", "IL1A", "IL1B", "IL1RA", "IL2", "IL13")
cytokine_genes %in% rownames(TrainData_TPM)
cytokine_genes<-cytokine_genes[cytokine_genes %in% rownames(TrainData_TPM)]
cytokine_expr <- TrainData_TPM[cytokine_genes, ]
geometric_mean <- function(x) {
  exp(mean(log(x + 1), na.rm = TRUE)) - 1
}
cytox_scores <- apply(cytokine_expr, 2, geometric_mean)
roc(TrainData_metadata$irAE_grade_binary, cytox_scores, quiet = TRUE)

## PMID£º39115425
# 5-gene SYK-GEP expression (CD22, PAG1, CD33, HNRNPU, and FCGR2C)
SYK_GEP<-c("CD22", "PAG1", "CD33", "HNRNPU", "FCGR2C")
SYK_GEP %in% rownames(TrainData_TPM)
SYK_GEP_expr <- TrainData_TPM[SYK_GEP, ]
geometric_mean <- function(x) {
  exp(mean(log(x + 1), na.rm = TRUE)) - 1
}
SYK_GEP_scores <- apply(SYK_GEP_expr, 2, geometric_mean)
roc(TrainData_metadata$irAE_grade_binary, SYK_GEP_scores, quiet = TRUE)

