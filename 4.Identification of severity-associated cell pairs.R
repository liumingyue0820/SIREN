# 1. Construction of possible gene-pair sets between any two cell subtypes ----
gene$CellType<-factor(gene$CellType,levels = c("Naive_T_cells_sp","Early_Activated_CD4_T_Cells_sp","Treg_sp","CD8_Te_sp","CD8_Tem_sp","Exhausted_CD8_T_sp","Mucosal_Associated_invariant_T_cells_sp",
        "Classical_Monocyte_sp","Non_classical_Monocyte_sp","Inflammatory_macrophage_sp",
        "CD56_bright_CD16_lo_NK_sp","CD56_dim_CD16_hi_NK_sp",
        "Memory_B_cell_sp","Plasma_cell_sp",
        "Myeloid_dendritic_cells_sp","Plasmacytoid_dendritic_cells_sp"))
Cell_Cell_list<-list()
names(table(gene$CellType))
s=1
for(i in 1:length(table(gene$CellType))){
  for(j in 1:length(table(gene$CellType))){
    if(i==j){next}
    if(i!=j){
      Cell1_marker<-gene[gene$CellType==names(table(gene$CellType))[i],2]
      Cell2_marker<-gene[gene$CellType==names(table(gene$CellType))[j],2]
      Cell_Cell_GPs<-matrix(,length(Cell1_marker)*length(Cell2_marker)+1,3)
      colnames(Cell_Cell_GPs)<-c("gene1","gene2","gene1_gene2")
      Cell_Cell_GPs[1,]<-c(names(table(gene$CellType))[i],names(table(gene$CellType))[j],paste(names(table(gene$CellType))[i],"_",names(table(gene$CellType))[j]))
      k=2
      for(m in 1:length(Cell1_marker)){
        for(n in 1:length(Cell2_marker)){
          Cell_Cell_GPs[k,]<-c(Cell1_marker[m],Cell2_marker[n],paste(Cell1_marker[m],"_",Cell2_marker[n]))
          k=k+1
        }
      }
      Cell_Cell_list[[s]]<-Cell_Cell_GPs
      s=s+1
    }
  }
}
length(Cell_Cell_list)  #Potential gene pairs across 240 cell subtype pairs

# 2. Identification of severity-associated cell pairs ----
phyper_CC_Pvalue<-matrix(,length(Cell_Cell_list),4)
colnames(phyper_CC_Pvalue)<-c("Cell_Cell_Type","P_value","Num_FatalGPs","Enrichment_Ratio")
phyper_CC_GPs<-list()
for(i in 1:length(Cell_Cell_list)){
  q<-length(intersect(GPs_Fatal$gene1_gene2,as.data.frame(Cell_Cell_list[[i]])$gene1_gene2))
  m<-nrow(GPs_Fatal)
  n<-nrow(exprset_REO_rank)-m
  k<-nrow(Cell_Cell_list[[i]])-1
  phyper_CC_Pvalue[i,1]<-as.data.frame(Cell_Cell_list[[i]])$gene1_gene2[1]
  phyper_CC_Pvalue[i,2]<-phyper(q-1, m, n, k, lower.tail=F)
  phyper_CC_Pvalue[i,3]<-q
  phyper_CC_Pvalue[i,4]<-(q/k)/(m/nrow(exprset_REO_rank))
  phyper_CC_GPs[[i]]<-intersect(GPs_Fatal$gene1_gene2,as.data.frame(Cell_Cell_list[[i]])$gene1_gene2)
}
phyper_CC_Pvalue<-as.data.frame(phyper_CC_Pvalue)

phyper_CC_Pvalue_ssc<-phyper_CC_Pvalue[which(as.numeric(phyper_CC_Pvalue$P_value)<0.05),]   
phyper_CC_Pvalue_ssc$Cell_Type_1<-""
phyper_CC_Pvalue_ssc$Cell_Type_2<-""
Cells<-strsplit(phyper_CC_Pvalue_ssc$Cell_Cell_Type, " _ ")
for(i in 1:nrow(phyper_CC_Pvalue_ssc)){
  phyper_CC_Pvalue_ssc$Cell_Type_1[i]<-Cells[[i]][1]
  phyper_CC_Pvalue_ssc$Cell_Type_2[i]<-Cells[[i]][2]
}
phyper_CC_GPs_ssc<-phyper_CC_GPs[which(as.numeric(phyper_CC_Pvalue$P_value)<0.05)]
sum(as.numeric(phyper_CC_Pvalue_ssc$Num_FatalGPs)) #6448 gene pairs

phyper_CC_Pvalue_ssc$P_value<-as.numeric(phyper_CC_Pvalue_ssc$P_value)
phyper_CC_Pvalue_ssc$Num_FatalGPs<-as.numeric(phyper_CC_Pvalue_ssc$Num_FatalGPs)

View(phyper_CC_Pvalue_ssc)  #67 cell subtype pairs






