# 1. Identification of severity-associated gene pairs ----
chisq_matrix<-matrix(,nrow(exprset_REO_rank),14)
colnames(chisq_matrix)<-c("gene1","gene2","gene1_gene2","0_TOX_high","0_TOX_low","1_TOX_high","1_TOX_low","EXP_0_TOX_high","EXP_0_TOX_low","EXP_1_TOX_high","EXP_1_TOX_low","0_TOX_high%","1_TOX_high%","P_value")
for(i in 1:nrow(exprset_REO_rank)){
  chisq_matrix[i,1:3]<-exprset_REO_rank[i,1:3]
  ssc<-cbind(exprset_REO_rank[i,4:56],matadata$irAE_grade)
  chisq_table<-matrix(c(sum(ssc[which(ssc[,1]==0),2]=="TOX_high"),sum(ssc[which(ssc[,1]==1),2]=="TOX_high"),sum(ssc[which(ssc[,1]==0),2]=="TOX_low"),sum(ssc[which(ssc[,1]==1),2]=="TOX_low")),2,2)
  chisq_matrix[i,c(4,6,5,7)]<-chisq_table[1:4]
  chisq_matrix[i,c(8,10,9,11)]<-chisq.test(chisq_table)$expected
  chisq_matrix[i,12]<-as.numeric(chisq_matrix[i,4])/(as.numeric(chisq_matrix[i,4])+as.numeric(chisq_matrix[i,5]))
  chisq_matrix[i,13]<-as.numeric(chisq_matrix[i,6])/(as.numeric(chisq_matrix[i,6])+as.numeric(chisq_matrix[i,7]))
  chisq_matrix[i,14]<-chisq.test(chisq_table)$p.value
}
ssc<-chisq_matrix[which(as.numeric(chisq_matrix[,14])<0.05),]  
GPs_Fatal<-as.data.frame(ssc[ssc[,12]<ssc[,13],])   
GPs_Fatal$Cell.type_gene1<-gene$CellType[match(GPs_Fatal$gene1,gene$Gene)]
GPs_Fatal$Cell.type_gene2<-gene$CellType[match(GPs_Fatal$gene2,gene$Gene)]
