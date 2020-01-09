##################CANCER-SMELL################
#The following code is the main pipeline for CancerSmell

#Initial step is to load the desired Receptor files
#f_or contains the list of Functional Olfactory Receptors
#p_or contains the list of Pseudo Olfactory Receptors
#t1r contains the list of Taste Receptors (T1Rs)
#t2r contains the list of Taste Receptors (T2Rs)

f_or<-read.csv("Receptor_files/functional_or.csv",header=TRUE)
p_or<-read.csv("Receptor_files/pseudo_or.csv",header=TRUE)
t1r<-read.csv("Receptor_files/t1r.csv",header=TRUE)
t2r<-read.csv("Receptor_files/t2r.csv",header=TRUE)
taar<-read.csv("Receptor_files/taar.csv",header=TRUE)
v1r<-read.csv("Receptor_files/v1r.csv",header=TRUE)

