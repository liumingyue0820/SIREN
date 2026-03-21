# 1. Construction of the gene pair matrix ----
exprset_REO<-exprset[gene$Gene,]
dim(exprset_REO)   #1255  53
exprset_REO_rank<-matrix(,nrow(exprset_REO)*(nrow(exprset_REO)-1),56)
dim(exprset_REO_rank) ##1573770  56
colnames(exprset_REO_rank)<-c("gene1","gene2","gene1_gene2",colnames(exprset_REO))
n=1
for(i in 1:nrow(exprset_REO)){
  for(j in 1:nrow(exprset_REO)){
    if(i==j){next}
    if(i!=j){
      for(k in 1:ncol(exprset_REO)){
        exprset_REO_rank[n,1]<-rownames(exprset_REO)[i]
        exprset_REO_rank[n,2]<-rownames(exprset_REO)[j]
        exprset_REO_rank[n,3]<-paste(rownames(exprset_REO)[i],"_",rownames(exprset_REO)[j])
        if(exprset_REO[i,k]>exprset_REO[j,k]){exprset_REO_rank[n,k+3]=1}  
        else if(exprset_REO[i,k]<=exprset_REO[j,k]){exprset_REO_rank[n,k+3]=0}  
      }
      n=n+1
    }
  }
}