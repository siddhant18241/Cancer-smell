setwd("/storage/gaurav.ahuja/siddhant/thesis_v2/Just_a_copy_of_exp0053")
f<-read.csv("Final_Receptor_file.csv")
f<-data.frame(f$OR1A1,f$OR2M3)
colnames(f)<-c("OR1A1","OR2M3")
x<-melt(as.data.table(f))
x <- na.omit(x)
colnames(x)<-c("Receptor","TPM")
ggplot(x, aes(x = TPM, fill = Receptor)) + geom_density(alpha = 0.3)
