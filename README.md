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

**File Manipulation** <br/>
Since the first row contains non-numneric values, the first line of the file was removed.

```r
exp<-exp[-1,]
```

Now our file appears like this 

| __CellName__ | __SRR2973279__ | __SRR2973280__ | __SRR2973281__ | __SRR2973282__ | __SRR2973283__ |
|-------------|------------|------------|------------|------------|------------|
| ENSG00000000003        | 0.76869 | 0.625649     | 0.546056     |  0.743321     | 0.163526    |
| ENSG00000000005         | 0 | 0     | 0     |  0     | 0    |
| ENSG00000000457         | 2.333028 | 0.593446     | 68.046238     |  8.740911     | 61.360909   |

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
After ID conversion file will apper like this.

| __CellName__ | __SRR2973279__ | __SRR2973280__ | __SRR2973281__ | __SRR2973282__ | __SRR2973283__ |
|-------------|------------|------------|------------|------------|------------|
| TSPAN6       | 0.76869 | 0.625649     | 0.546056     |  0.743321     | 0.163526    |
| TNMD         | 0 | 0     | 0     |  0     | 0    |
| DPM1         | 2.333028 | 0.593446     | 68.046238     |  8.740911     | 61.360909   |

**Saving Files**
```r
exp<-setDT(exp, keep.rownames = "Gene_names")[]
write.table(exp,file="Main_pipeline/GSE75688/raw_input.csv",sep=",",quote=FALSE,row.names = FALSE,col.names = TRUE)
t_exp<-setDT(t_exp, keep.rownames = "Cell_names")[]
write.table(t_exp,file="Main_pipeline/GSE75688/raw_input_for_doublet.csv",sep=",",quote=FALSE,
            row.names = FALSE,col.names = TRUE)
```

**Doublet Identification** <br/>
For this purpose python package *"Scrublet"* (0.2.1) was used.

```python
import pandas as pd
import numpy as np
import scrublet as scr
data = pd.io.parsers.read_csv("Main_pipeline/GSE75688/raw_input_for_doublet.csv",sep=",",index_col=0)
indexNamesArr = data.index.values
listOfRowIndexLabels = list(indexNamesArr)
listOfRowIndexLabels
row_names=pd.DataFrame(listOfRowIndexLabels)
scrub = scr.Scrublet(data)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
score=pd.DataFrame(doublet_scores)
value=pd.DataFrame(predicted_doublets)
a=np.concatenate((row_names,score),axis=1)
b=np.concatenate((a,value),axis=1)
final=pd.DataFrame(b)
final.columns = ['Cell','Doublet_Score','Doublet']
print(final)
final.to_csv(r'Main_pipeline/GSE75688/doublet_analysis.csv',index = None, header=True)
exit()
```

**Alternatively** the above code has been saved in the file *doublet.py* and can be called in R only using **reticulate()** library. For this you need to pass the location where your python is installed in the function *use_python*. For eg

```r
use_python("/home/siddhant18241/miniconda3/bin/python3", required = T)
py_config()
source_python("Main_pipeline/GSE75688/doublet.py")
```
Doublet Analysis will output the result in "doublet_analysis.csv".
This is the three column file that contains Cell_Id, Doublet_Score, and Doublet_Status
| __Cell__ | __Doublet_Score__ | __Doublet__ |
|-------------|------------|------------|
| SRR2973279      | 0.02072538860103627 | False     |
| SRR2973280       | 0.017369727047146403 | False     |
| SRR2973281       | 0.017369727047146403 | False     |


Doublet Status will be in the form of TRUE/FALSE, only FALSE cells are further used for analysis.


