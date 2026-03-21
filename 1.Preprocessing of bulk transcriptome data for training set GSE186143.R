#####Data description
##66 samples
##5 healthy controls;
##61 metastatic melanoma patient samples: 
  #bulk cohort 1, n = 26; ¡ª¡ªsequencing batch: B
  #bulk cohort 2, n = 27; ¡ª¡ª18 sequencing batch: A;9  sequencing batch: C
  #8 patient samples used for technical experiments.(4A,4C) ¡¾YUTRIL¡¿ severity grade data missing, consider excluding

# 1. Loading R packages ----
library(readxl)
library(sva)
library(tinyarray)
library(dplyr)

# 2. Load Expression Data (TPM) ----
data_TPM <- read.table("...\\GSE186143_bulk_TPM_samples.txt", header=TRUE,fill=TRUE)
data_TPM_matedata<-as.data.frame(read_excel("...\\GSE186143_series_matrix.xlsx"))
a<-c("YUFURL","YUTORY","YUHERN","YUROD","YUALOE","YUNANCY","WUMENLO","YUHEAL","YUBAIX","YUDONA","YUKEID","YURADIO","YUBRET","YUCUSK","YUGIM","YUGRUE","YUMAXI","YUMEDIC","YUQIP","YURAISE","YURICA","YUEVRY","YUPRAF","YUROSA","YUBIFO","YUBUMP","YUDIME","YUFROME","YUGIRL","YUHIMO","YUJADE","YUPARK","YUTACI","YUAMIGO","YUBROMO","YUCADA","YUFROZE","YULAMB","YUPOGER","YUSMAL","YUGRAPH","YUPRETE","YUBLIT","YUPESO","YUDEBIT","YURIDA","YUVITZ","YUARGE","YUCOT","YUSTE","YUZIRG","YUHONEY","YUDUAR","YUSONG","WULOUIS","WUPALO","YUCEVO","YUDNA","YUWHEY","YUADD","YUTRIL","H-1","H-2","H-3","H-4","H-5")
data_TPM_matedata$Sample_names<-a

# 3. Handle Duplicate Gene Symbols ----
countdata<-data_TPM   
dim(countdata)  #56809  67
countdata <- aggregate(countdata[,2:67], by=list(countdata$GENE), FUN=mean)
row.names(countdata)<-countdata$Group.1
countdata<-countdata[,-1]
dim(countdata)  #56809  66

# 4. Selecting samples from bulk cohort 1 and bulk cohort 2 as the training set ----
cohort_1_2<-c(data_TPM_matedata$Sample_names[grep("cohort: Bulk cohort 1",data_TPM_matedata$Sample_characteristics_ch1...10)],
              data_TPM_matedata$Sample_names[grep("cohort: Bulk cohort 2",data_TPM_matedata$Sample_characteristics_ch1...10)])
countdata<-countdata[,cohort_1_2]
dim(countdata)  #56809  53
matadata<-data_TPM_matedata[data_TPM_matedata$Sample_names%in%cohort_1_2,]
matadata$irAE_grade<-"TOX_low"
matadata$irAE_grade[matadata$Sample_characteristics_ch1...9 %in% c("highest irae grade: 3","highest irae grade: 4")]<-"TOX_high"

# 5. Filter Lowly Expressed Genes ----
countdata=countdata[rowSums(countdata>0) >= floor(3/4*ncol(countdata)),]
dim(countdata)   #20013     53

# 6. Batch Correction (ComBat) ----
exprset <- log2(countdata+0.1)
batch_info <- factor(matadata$Sample_characteristics_ch1...11)
condition <- factor(matadata$irAE_grade)
mod <- model.matrix(~ condition)  
exprset <- ComBat(dat = expr_matrix, batch = batch_info, mod = NULL, par.prior = TRUE)










