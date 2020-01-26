library("Matrix")
setwd("/storage/gaurav.ahuja/siddhant/thesis_v2/EXP0053/Normal_breast_TPM")
counts <- readMM("E-CURD-7.expression_tpm.mtx")
genes <- read.csv("E-CURD-7.expression_tpm.mtx_rows", header=FALSE,sep="\t")
gene_ids <- genes$V1
cell_ids <- read.csv("E-CURD-7.expression_tpm.mtx_cols",header=FALSE)$V1
rownames(counts) <- gene_ids
colnames(counts) <- cell_ids
counts_2<-counts
counts_2<-as.data.frame(counts_2)
View(counts_2[1:10,1:10])
counts_2<-setDT(counts_2,keep.rownames = "ID")[]
write.table(counts_2,file="Normal_breast_TPM.csv",sep=",",row.names=FALSE,quote = FALSE)
