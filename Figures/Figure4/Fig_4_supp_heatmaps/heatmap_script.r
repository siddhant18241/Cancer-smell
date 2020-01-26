#fig4
setwd("/home/siddhant/Desktop/Fig_4_supp_heatmaps")
phenotypes<-read.csv("phenotype.tsv",sep="\t")
exp52<-read.csv("exp0052.csv")
exp52<-merge(exp52,phenotypes,by.x="sample","samples")
write.table(exp52,file="EXP0052.csv",row.names = FALSE,quote = FALSE,sep=",")

exp53<-read.csv("exp0053.csv")
exp53<-merge(exp53,phenotypes,by.x="sample","samples")
write.table(exp53,file="EXP0053.csv",row.names = FALSE,quote = FALSE,sep=",")

exp54<-read.csv("exp0054.csv")
exp54<-merge(exp54,phenotypes,by.x="sample","samples")
write.table(exp54,file="EXP0054.csv",row.names = FALSE,quote = FALSE,sep=",")

#install.packages("pheatmap")
library(heatmap3)
library(RColorBrewer)
library(pheatmap)
my_palette <- colorRampPalette(c("tomato4", "snow1", "slateblue1"))(n = 299)
setwd("/home/siddhant/Desktop/Fig_4_supp_heatmaps")
pdf(file="Supp_heatmaps.pdf")
x<-read.csv("heatmap_0052.csv",row.names = 1)
#heatmap3(x,cexRow=1,cexCol = 1,main="PDX Labels vs Tumor Stages",col=my_palette,cex.main = 1)
pheatmap(x, color = my_palette,main="PDX Labels vs Tumor Stages",cluster_rows = FALSE)

x<-read.csv("heatmap_0053.csv",row.names = 1)
pheatmap(x, color = my_palette,main="Tissue Labels vs Tumor Stages",cluster_rows = FALSE)

x<-read.csv("heatmap_0054.csv",row.names = 1)
pheatmap(x, color = my_palette,main="CTC Labels vs Tumor Stages",cluster_rows = FALSE)

dev.off()
