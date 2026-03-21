# 1. Risk score calculation for individual cell pairs ----
Single_CellPair_irAE_Risk<-matrix(,length(phyper_CC_GPs_ssc),4+nrow(matadata))
colnames(Single_CellPair_irAE_Risk)<-c("Cell_Cell_Type","Gene_pair_TrueNum","Gene_pair_TestNum","Gene_pair_Loss",matadata$Sample_names)
for(i in 1:length(phyper_CC_GPs_ssc)){
  Cell_Cell_pair<-phyper_CC_Pvalue_ssc$Cell_Cell_Type[i]
  Gene_pair<-as.character(phyper_CC_GPs_ssc[[i]])
  Gene_pair_TrueNum<-length(Gene_pair)
  
  Gene_pair_Test<-Gene_pair[Gene_pair %in% exprset_REO_rank[,3]]
  Gene_pair_TestNum<-sum(Gene_pair %in% exprset_REO_rank[,3])
  
  Gene_pair_Loss<- paste0(Gene_pair[!(Gene_pair %in% exprset_REO_rank[,3])],"/")
  
  a<-exprset_REO_rank[Gene_pair_Test,4:56]
  b<-matrix(as.numeric(a), nrow = nrow(a), ncol = ncol(a))
  b[is.na(b)]=0
  irAE_Risk<-as.numeric(colSums(b)/Gene_pair_TestNum)
  
  Single_CellPair_irAE_Risk[i,1]<-Cell_Cell_pair
  Single_CellPair_irAE_Risk[i,2]<-Gene_pair_TrueNum
  Single_CellPair_irAE_Risk[i,3]<-Gene_pair_TestNum
  Single_CellPair_irAE_Risk[i,4]<-Gene_pair_Loss
  Single_CellPair_irAE_Risk[i,5:57]<-irAE_Risk
}

# 2. ROC-based evaluation of risk score predictive performance ----
Single_CellPair_irAE_Risk<-as.data.frame(Single_CellPair_irAE_Risk)
Single_CellPair_irAE_Risk$AUC<-""
for(i in 1:nrow(Single_CellPair_irAE_Risk)){
  true_labels <- matadata$irAE_grade
  predicted_scores<-as.numeric(Single_CellPair_irAE_Risk[i,5:57])
  roc_obj <- roc(true_labels, predicted_scores, levels = c("TOX_low", "TOX_high"),direction = "<") ##control VS case
  Single_CellPair_irAE_Risk$AUC[i]<-round(auc(roc_obj), 3)
}
