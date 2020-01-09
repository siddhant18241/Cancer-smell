# Cancer-smell
## Workflow
<img src="Main_pipeline/Images/Workflow.png">

**Load essential Libraries** <br/>
The initial step is to load all the necessary libraries

```r
#Loading Essential Libraries
library("reticulate") #To use one Python Library "Scrublet"
library("AnnotationDbi") #To convert one Id to Another
library("org.Hs.eg.db") #To use Genome wide annotation for Human
library(dplyr) #For Datasets Manipulation
library(tidyr) #For Datasets Manipulation
library(data.table) #For Datasets Manipulation
library(zFPKM) #To calulate zFPKM scores
library(Seurat) #For single cell analysis
library(ggplot2) #To plot graphs
library(RColorBrewer) #To define color-palette
library("GSVA") #To compute Gene Set Variation Analysis
library(GSEABase) #To compute Gene Set Variation Analysis
library(GSVAdata) #To compute Gene Set Variation Analysis
library(msigdb) #To load gmt file
```

Once all the libraries are properly loaded, we need to provide input files. This pipeline takes into account 6 sub famlies of GPCRs. The format of all these 6 Input Files is, a sinlge column csv file that contains **gene names** as entries and column heading as **"Symbol"**  <br/>

**Input Receptor files** <br/>
```r
f_or<-read.csv("Main_pipeline/Receptor_files/functional_or.csv",header=TRUE)
p_or<-read.csv("Main_pipeline/Receptor_files/pseudo_or.csv",header=TRUE)
t1r<-read.csv("Main_pipeline/Receptor_files/t1r.csv",header=TRUE)
t2r<-read.csv("Main_pipeline/Receptor_files/t2r.csv",header=TRUE)
taar<-read.csv("Main_pipeline/Receptor_files/taar.csv",header=TRUE)
v1r<-read.csv("Main_pipeline/Receptor_files/v1r.csv",header=TRUE)
```

**Loading Expression file** <br/>
In this documntation, the sample file used is under **GSE75688** folder. The structure of the file is shown below.

| __CellName__ | __SRR2973279__ | __SRR2973280__ | __SRR2973281__ | __SRR2973282__ | __SRR2973283__ |
|-------------|------------|------------|------------|------------|------------|
| Sample         | BC01     | BC01     | BC01     | BC01     | BC01     |
| ENSG00000000003        | 0.76869 | 0.625649     | 0.546056     |  0.743321     | 0.163526    |
| ENSG00000000005         | 0 | 0     | 0     |  0     | 0    |
| ENSG00000000457         | 2.333028 | 0.593446     | 68.046238     |  8.740911     | 61.360909   |

```r
#Reading Raw tpm matrix
exp<-read.csv("Main_pipeline/GSE75688/GSE75688_PCG_beforeQC.txt",header=TRUE,
               row.names=1,sep="\t",stringsAsFactors=FALSE)
```
Since the first row contains non-numneric values, the first line of the file was removed.

```r
exp<-exp[-1,]
```

The Gene Id's present in this input file is in ENSEMBL format whereas our receptor file contains gene names. So we need to convert these ENSEMBL ID's to Gene Names

```r
exp$Gene <- mapIds(org.Hs.eg.db,keys=row.names(exp),column="SYMBOL", keytype="ENSEMBL", multiVals="first")
exp<-distinct(exp,Gene, .keep_all= TRUE)
exp<-exp %>% drop_na(Gene)
row.names(exp)<-exp$Gene
exp$Gene<-NULL
t_exp<-t(exp)
t_exp<-as.data.frame(t_exp)
```

