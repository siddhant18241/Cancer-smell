#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#Category wise Monocle
##MONOCole
library(ggplot2)
library(dplyr)
library(monocle3)
library("PerformanceAnalytics")
setwd("/storage/gaurav.ahuja/siddhant/thesis_v2/Monocle_categories")
#########################################BASAL##################################3
#File1
expression_matrix<-read.csv("Basal_like/Basal_expression.csv",row.names = 1)
RN<-row.names(expression_matrix)
expression_matrix<-as.matrix(expression_matrix)

#File2
cell_metadata<-read.csv("Basal_like/Basal_Cell_Meta_Data_2.csv",row.names = 1)

#File3
gene_annotation<-read.csv("Basal_like/gene_meta_data.csv")
gene_annotation$genes<-NULL
row.names(gene_annotation)<-RN

#Module1
# Make the CDS object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

pdf(file="Final_pdf_V2/Basal_like.pdf")
#Cluster and classify cells
cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds) + ggtitle("PC Variance")
cds <- reduce_dimension(cds)
cds = cluster_cells(cds, resolution=1e-5)

plot_cells(cds,cell_size = 1.5)
plot_cells(cds, color_cells_by="partition", group_cells_by="partition",cell_size = 1.5)
plot_cells(cds, color_cells_by="OR_count_per_cell_zfpkm",cell_size = 1.5)+ggtitle("OR_count_per_cell_zfpkm")
plot_cells(cds, color_cells_by="Fire_score",cell_size = 1.5)+ggtitle("Fire_score")
plot_cells(cds, color_cells_by="Cell_type",cell_size = 1.5)+ggtitle("Cell_type")
plot_cells(cds, color_cells_by="Molecular_subtype",cell_size = 1.5)+ggtitle("Molecular_subtype")
plot_cells(cds, color_cells_by="Pathological_stage",cell_size = 1.5)+ggtitle("Pathological_stage")
plot_cells(cds, color_cells_by="Immunohistochemistry",cell_size = 1.5)+ggtitle("Immunohistochemistry")
plot_cells(cds, color_cells_by="patient_id",cell_size = 1.5)+ggtitle("patient_id")
plot_cells(cds, color_cells_by="patient_age",cell_size = 1.5)+ggtitle("patient_age")
plot_cells(cds, color_cells_by="Mean_OR_TPM",cell_size = 1.5)+ggtitle("Mean_OR_TPM")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_2",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_2")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_4",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_4")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_5",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_5")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_6",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_6")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_10",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_10")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_15",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_15")
#plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_20",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_20")
plot_cells(cds, color_cells_by="Vector_size",cell_size = 1.5)+ggtitle("Vector Size")

# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="Cell_type",
#                     ordering_type="cluster_row_col",
#                     max.size=3,cell_size = 1.5)+ggtitle("Cell_type")

# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="Pathological_stage",
#                     ordering_type="cluster_row_col",
#                     max.size=3,cell_size = 1.5)+ggtitle("Pathological_stage")

# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="Immunohistochemistry",
#                     ordering_type="cluster_row_col",
#                     max.size=3,cell_size = 1.5)+ggtitle("Immunohistochemistry")

#Module3
cds <- learn_graph(cds)
#learn_graph(cds)
# plot_cells(cds,
#            color_cells_by = "Pathological_stage",
#            label_groups_by_cluster=FALSE,
#            label_leaves=FALSE,
#            label_branch_points=FALSE,cell_size = 1.5)
# plot_cells(cds,
#            color_cells_by = "Stemness",
#            label_groups_by_cluster=FALSE,
#            label_leaves=FALSE,
#            label_branch_points=FALSE,cell_size = 1.5)
# plot_cells(cds,
#            color_cells_by = "Differentiation",
#            label_groups_by_cluster=FALSE,
#            label_leaves=FALSE,
#            label_branch_points=FALSE,cell_size = 1.5)

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="High_Stemness"){
  cell_ids <- which(colData(cds)[, "State"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
pseudo<-cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
pseudo<-as.data.frame(pseudo)
colnames(pseudo)<-"PseudoTime"
#View(order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds)))
# 
# traj.plot <- plot_cell_trajectory(cds)
# point.data <- ggplot_build(cds)[["plot"]][["data"]]
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)

plot_cells(cds,
           color_cells_by = "OR_count_per_cell_zfpkm",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("OR_count_per_cell_zfpkm")

plot_cells(cds,
           color_cells_by = "Fire_score",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Fire_score")

plot_cells(cds,
           color_cells_by = "Cell_type",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Cell_type")

plot_cells(cds,
           color_cells_by = "Vector_size",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Vector_size")

plot_cells(cds,
           color_cells_by = "Pathological_stage",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Pathological_stage")

plot_cells(cds,
           color_cells_by = "Immunohistochemistry",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Immunohistochemistry")

plot_cells(cds,
           color_cells_by = "patient_id",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("patient_id")

plot_cells(cds,
           color_cells_by = "patient_age",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("patient_age")

plot_cells(cds,
           color_cells_by = "Mean_OR_TPM",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Mean_OR_TPM")
plot_cells(cds,
           color_cells_by = "Stemness",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Stemness")
plot_cells(cds,
           color_cells_by = "Differentiation",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Differentiation")
plot_cells(cds,
           color_cells_by = "Molecular_subtype",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Molecular_subtype")


ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))
test_results<-as.data.frame(ciliated_cds_pr_test_res)
test_results<-setDT(test_results,keep.rownames = "Genes")[]
write.table(test_results,"Basal_like/Test_results.csv",row.names = FALSE,sep=",",quote = FALSE)
AFD_genes<-read.csv("/storage/gaurav.ahuja/siddhant/thesis_v2/Monocle_categories/Genes.csv",header=FALSE)
AFD_genes<-as.matrix(AFD_genes)
a<-length(AFD_genes)
for (j in 1:a){
  tryCatch(
    {
      
      AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes[j],
                             colData(cds)$Molecular_subtype %in% c("Basal-like")]
      k<-plot_genes_in_pseudotime(AFD_lineage_cds,
                               color_cells_by="OR_count_per_cell_zfpkm",
                               min_expr=0.5,cell_size = 1.5)
      plot(k)
      
    },
    error=function(cond) {
      message(paste("."))
      
    },
    warning=function(cond) {
      message(paste("."))
      
    },
    finally={
      
    }
  )    
}

for (j in 1:a){
  tryCatch(
    {
      
      AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes[j],
                             colData(cds)$Molecular_subtype %in% c("Basal-like")]
      k<-plot_genes_in_pseudotime(AFD_lineage_cds,
                                  color_cells_by="Fire_score",
                                  min_expr=0.5,cell_size = 1.5)
      plot(k)
      
    },
    error=function(cond) {
      message(paste("."))
      
    },
    warning=function(cond) {
      message(paste("."))
      
    },
    finally={
      
    }
  )    
}

for (j in 1:a){
  tryCatch(
    {
      
      AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes[j],
                             colData(cds)$Molecular_subtype %in% c("Basal-like")]
      k<-plot_genes_in_pseudotime(AFD_lineage_cds,
                                  color_cells_by="Mean_OR_TPM",
                                  min_expr=0.5,cell_size = 1.5)
      plot(k)
      
    },
    error=function(cond) {
      message(paste("."))
      
    },
    warning=function(cond) {
      message(paste("."))
      
    },
    finally={
      
    }
  )    
}

#View(cell_metadata)
Meta_data_pseudo<-cbind(cell_metadata,pseudo)
metadata_with_pesudotime<-setDT(Meta_data_pseudo,keep.rownames = "ID")[]
write.table(metadata_with_pesudotime,file="Basal_like/Metadata_with_pseudotime.csv",sep=",",quote=FALSE,row.names = FALSE)
sp <- ggscatter(Meta_data_pseudo, x = "PseudoTime", y = "Stemness",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                xlab = "PseudoTime",
                ylab = "Stemness"
                
)
sp+stat_cor(
  aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
  label.x = 4,label.y = 0.2
)
sp <- ggscatter(Meta_data_pseudo, x = "PseudoTime", y = "Differentiation",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                xlab = "PseudoTime",
                ylab = "Differentiation"
                
)
sp+stat_cor(
  aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
  label.x = 4,label.y = 0.2
)

sp <- ggscatter(Meta_data_pseudo, x = "PseudoTime", y = "Mean_OR_TPM",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                xlab = "PseudoTime",
                ylab = "Mean_OR_TPM"
                
)
sp+stat_cor(
  aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
  label.x = 4,label.y = 5
)

sp <- ggscatter(Meta_data_pseudo, x = "PseudoTime", y = "Fire_score",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                xlab = "PseudoTime",
                ylab = "Fire_score"
                
)
sp+stat_cor(
  aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
  label.x = 4,label.y = 950
)

sp <- ggscatter(Meta_data_pseudo, x = "PseudoTime", y = "OR_count_per_cell_zfpkm",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                xlab = "PseudoTime",
                ylab = "OR_count_per_cell_zfpkm"
                
)
sp+stat_cor(
  aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
  label.x = 4,label.y = 5
)

ggplot(Meta_data_pseudo, aes(PseudoTime, Stemness)) +geom_bar(position="identity",stat="identity", width = 0.2)
ggplot(Meta_data_pseudo, aes(PseudoTime, Differentiation)) +geom_bar(position="identity",stat="identity", width = 0.2)


newdata<-data.frame(Meta_data_pseudo$OR_count_per_cell_zfpkm,Meta_data_pseudo$Fire_score,Meta_data_pseudo$Mean_OR_TPM,Meta_data_pseudo$Stemness,Meta_data_pseudo$Differentiation,Meta_data_pseudo$PseudoTime)
colnames(newdata)<-c("OR_count","Fire_score","Mean_TPM","Stemness","Differentiation","PseudoTime")

chart.Correlation(newdata, histogram=TRUE, pch=19)

# stem<-cell_metadata[which(cell_metadata$Cell_staus=="Stemness"),]
# dif<-cell_metadata[which(cell_metadata$Cell_staus=="Differentiation"),]
# 
# df<-data.frame(stem$Stemness,stem$Mean_OR_TPM)
# colnames(df)<-c("Stemness","Mean_OR_TPM")
# sp <- ggscatter(df, x = "Stemness", y = "Mean_OR_TPM",
#                 add = "reg.line",  # Add regressin line
#                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#                 conf.int = TRUE,
#                 xlab = "Stemness",
#                 ylab = "Mean_OR_TPM"
#                 
# )
# sp+stat_cor(
#   aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
#   label.x = 0.1,label.y = 1
# )
# 
# df<-data.frame(stem$Stemness,stem$OR_count_per_cell_zfpkm)
# colnames(df)<-c("Stemness","OR_count_per_cell_zfpkm")
# sp <- ggscatter(df, x = "Stemness", y = "OR_count_per_cell_zfpkm",
#                 add = "reg.line",  # Add regressin line
#                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#                 conf.int = TRUE,
#                 xlab = "Stemness",
#                 ylab = "OR_count_per_cell_zfpkm"
#                 
# )
# sp+stat_cor(
#   aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
#   label.x = 0,label.y = 25
# )
# 
# df<-data.frame(dif$Differentiation,dif$Mean_OR_TPM)
# colnames(df)<-c("Differentiation","Mean_OR_TPM")
# sp <- ggscatter(df, x = "Differentiation", y = "Mean_OR_TPM",
#                 add = "reg.line",  # Add regressin line
#                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#                 conf.int = TRUE,
#                 xlab = "Differentiation",
#                 ylab = "Mean_OR_TPM"
#                 
# )
# sp+stat_cor(
#   aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
#   label.x = -0.1,label.y = 7.5
# )
# 
# df<-data.frame(dif$Differentiation,dif$OR_count_per_cell_zfpkm)
# colnames(df)<-c("Differentiation","OR_count_per_cell_zfpkm")
# sp <- ggscatter(df, x = "Differentiation", y = "OR_count_per_cell_zfpkm",
#                 add = "reg.line",  # Add regressin line
#                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#                 conf.int = TRUE,
#                 xlab = "Differentiation",
#                 ylab = "OR_count_per_cell_zfpkm"
#                 
# )
# sp+stat_cor(
#   aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
#   label.x = 0,label.y = 20
# )


dev.off()

##################HER2
#########################################BASAL##################################3
#File1
expression_matrix<-read.csv("Her2/Her2_expression.csv",row.names = 1)
RN<-row.names(expression_matrix)
expression_matrix<-as.matrix(expression_matrix)

#File2
cell_metadata<-read.csv("Her2/Her2_Cell_Meta_Data_2.csv",row.names = 1)

#File3
gene_annotation<-read.csv("Her2/gene_meta_data.csv")
gene_annotation$genes<-NULL
row.names(gene_annotation)<-RN

#Module1
# Make the CDS object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

pdf(file="Final_pdf_V2/HER2.pdf")
#Cluster and classify cells
cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds) + ggtitle("PC Variance")
cds <- reduce_dimension(cds)
cds = cluster_cells(cds, resolution=1e-5)

plot_cells(cds,cell_size = 1.5)
plot_cells(cds, color_cells_by="partition", group_cells_by="partition",cell_size = 1.5)
plot_cells(cds, color_cells_by="OR_count_per_cell_zfpkm",cell_size = 1.5)+ggtitle("OR_count_per_cell_zfpkm")
plot_cells(cds, color_cells_by="Fire_score",cell_size = 1.5)+ggtitle("Fire_score")
plot_cells(cds, color_cells_by="Cell_type",cell_size = 1.5)+ggtitle("Cell_type")
plot_cells(cds, color_cells_by="Molecular_subtype",cell_size = 1.5)+ggtitle("Molecular_subtype")
plot_cells(cds, color_cells_by="Pathological_stage",cell_size = 1.5)+ggtitle("Pathological_stage")
plot_cells(cds, color_cells_by="Immunohistochemistry",cell_size = 1.5)+ggtitle("Immunohistochemistry")
plot_cells(cds, color_cells_by="patient_id",cell_size = 1.5)+ggtitle("patient_id")
plot_cells(cds, color_cells_by="patient_age",cell_size = 1.5)+ggtitle("patient_age")
plot_cells(cds, color_cells_by="Mean_OR_TPM",cell_size = 1.5)+ggtitle("Mean_OR_TPM")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_2",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_2")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_4",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_4")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_5",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_5")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_6",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_6")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_10",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_10")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_15",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_15")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_20",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_20")
plot_cells(cds, color_cells_by="Vector_size",cell_size = 1.5)+ggtitle("Vector Size")

# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="Cell_type",
#                     ordering_type="cluster_row_col",
#                     max.size=3,cell_size = 1.5)+ggtitle("Cell_type")
# 
# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="Pathological_stage",
#                     ordering_type="cluster_row_col",
#                     max.size=3,cell_size = 1.5)+ggtitle("Pathological_stage")
# 
# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="Immunohistochemistry",
#                     ordering_type="cluster_row_col",
#                     max.size=3,cell_size = 1.5)+ggtitle("Immunohistochemistry")
# 
#Module3
cds <- learn_graph(cds)
# plot_cells(cds,
#            color_cells_by = "Pathological_stage",
#            label_groups_by_cluster=FALSE,
#            label_leaves=FALSE,
#            label_branch_points=FALSE,cell_size = 1.5)

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="Moderate_Stemness"){
  cell_ids <- which(colData(cds)[, "State"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
pseudo<-cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
pseudo<-as.data.frame(pseudo)
colnames(pseudo)<-"PseudoTime"
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)

plot_cells(cds,
           color_cells_by = "OR_count_per_cell_zfpkm",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("OR_count_per_cell_zfpkm")

plot_cells(cds,
           color_cells_by = "Fire_score",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Fire_score")

plot_cells(cds,
           color_cells_by = "Cell_type",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Cell_type")

plot_cells(cds,
           color_cells_by = "Vector_size",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Vector_size")

plot_cells(cds,
           color_cells_by = "Pathological_stage",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Pathological_stage")

plot_cells(cds,
           color_cells_by = "Immunohistochemistry",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Immunohistochemistry")

plot_cells(cds,
           color_cells_by = "patient_id",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("patient_id")

plot_cells(cds,
           color_cells_by = "patient_age",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("patient_age")

plot_cells(cds,
           color_cells_by = "Mean_OR_TPM",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Mean_OR_TPM")
plot_cells(cds,
           color_cells_by = "Stemness",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Stemness")
plot_cells(cds,
           color_cells_by = "Differentiation",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Differentiation")
plot_cells(cds,
           color_cells_by = "Molecular_subtype",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Molecular_subtype")

ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))
test_results<-as.data.frame(ciliated_cds_pr_test_res)
test_results<-setDT(test_results,keep.rownames = "Genes")[]
write.table(test_results,"Her2/Test_results.csv",row.names = FALSE,sep=",",quote = FALSE)
AFD_genes<-read.csv("/storage/gaurav.ahuja/siddhant/thesis_v2/Monocle_categories/Genes.csv",header=FALSE)
AFD_genes<-as.matrix(AFD_genes)
a<-length(AFD_genes)
for (j in 1:a){
  tryCatch(
    {
      
      AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes[j],
                             colData(cds)$Molecular_subtype %in% c("HER2-enriched")]
      k<-plot_genes_in_pseudotime(AFD_lineage_cds,
                                  color_cells_by="OR_count_per_cell_zfpkm",
                                  min_expr=0.5,cell_size = 1.5)
      plot(k)
      
    },
    error=function(cond) {
      message(paste("."))
      
    },
    warning=function(cond) {
      message(paste("."))
      
    },
    finally={
      
    }
  )    
}

for (j in 1:a){
  tryCatch(
    {
      
      AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes[j],
                             colData(cds)$Molecular_subtype %in% c("HER2-enriched")]
      k<-plot_genes_in_pseudotime(AFD_lineage_cds,
                                  color_cells_by="Fire_score",
                                  min_expr=0.5,cell_size = 1.5)
      plot(k)
      
    },
    error=function(cond) {
      message(paste("."))
      
    },
    warning=function(cond) {
      message(paste("."))
      
    },
    finally={
      
    }
  )    
}

for (j in 1:a){
  tryCatch(
    {
      
      AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes[j],
                             colData(cds)$Molecular_subtype %in% c("HER2-enriched")]
      k<-plot_genes_in_pseudotime(AFD_lineage_cds,
                                  color_cells_by="Mean_OR_TPM",
                                  min_expr=0.5,cell_size = 1.5)
      plot(k)
      
    },
    error=function(cond) {
      message(paste("."))
      
    },
    warning=function(cond) {
      message(paste("."))
      
    },
    finally={
      
    }
  )    
}


Meta_data_pseudo<-cbind(cell_metadata,pseudo)
metadata_with_pesudotime<-setDT(Meta_data_pseudo,keep.rownames = "ID")[]
write.table(metadata_with_pesudotime,file="Her2/Metadata_with_pseudotime.csv",sep=",",quote=FALSE,row.names = FALSE)
sp <- ggscatter(Meta_data_pseudo, x = "PseudoTime", y = "Stemness",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                xlab = "PseudoTime",
                ylab = "Stemness"
                
)
sp+stat_cor(
  aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
  label.x = 2,label.y = 0.2
)
sp <- ggscatter(Meta_data_pseudo, x = "PseudoTime", y = "Differentiation",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                xlab = "PseudoTime",
                ylab = "Differentiation"
                
)
sp+stat_cor(
  aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
  label.x = 2,label.y = -0.15
)

sp <- ggscatter(Meta_data_pseudo, x = "PseudoTime", y = "Mean_OR_TPM",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                xlab = "PseudoTime",
                ylab = "Mean_OR_TPM"
                
)
sp+stat_cor(
  aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
  label.x = 2,label.y = 5
)

sp <- ggscatter(Meta_data_pseudo, x = "PseudoTime", y = "Fire_score",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                xlab = "PseudoTime",
                ylab = "Fire_score"
                
)
sp+stat_cor(
  aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
  label.x = 2,label.y = 950
)

sp <- ggscatter(Meta_data_pseudo, x = "PseudoTime", y = "OR_count_per_cell_zfpkm",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                xlab = "PseudoTime",
                ylab = "OR_count_per_cell_zfpkm"
                
)
sp+stat_cor(
  aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
  label.x = 4,label.y = 5
)

ggplot(Meta_data_pseudo, aes(PseudoTime, Stemness)) +geom_bar(position="identity",stat="identity", width = 0.2)
ggplot(Meta_data_pseudo, aes(PseudoTime, Differentiation)) +geom_bar(position="identity",stat="identity", width = 0.2)


newdata<-data.frame(Meta_data_pseudo$OR_count_per_cell_zfpkm,Meta_data_pseudo$Fire_score,Meta_data_pseudo$Mean_OR_TPM,Meta_data_pseudo$Stemness,Meta_data_pseudo$Differentiation,Meta_data_pseudo$PseudoTime)
colnames(newdata)<-c("OR_count","Fire_score","Mean_TPM","Stemness","Differentiation","PseudoTime")

chart.Correlation(newdata, histogram=TRUE, pch=19)

# stem<-cell_metadata[which(cell_metadata$Cell_staus=="Stemness"),]
# dif<-cell_metadata[which(cell_metadata$Cell_staus=="Differentiation"),]
# 
# df<-data.frame(stem$Stemness,stem$Mean_OR_TPM)
# colnames(df)<-c("Stemness","Mean_OR_TPM")
# sp <- ggscatter(df, x = "Stemness", y = "Mean_OR_TPM",
#                 add = "reg.line",  # Add regressin line
#                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#                 conf.int = TRUE,
#                 xlab = "Stemness",
#                 ylab = "Mean_OR_TPM"
#                 
# )
# sp+stat_cor(
#   aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
#   label.x = 0.1,label.y = 1
# )
# 
# df<-data.frame(stem$Stemness,stem$OR_count_per_cell_zfpkm)
# colnames(df)<-c("Stemness","OR_count_per_cell_zfpkm")
# sp <- ggscatter(df, x = "Stemness", y = "OR_count_per_cell_zfpkm",
#                 add = "reg.line",  # Add regressin line
#                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#                 conf.int = TRUE,
#                 xlab = "Stemness",
#                 ylab = "OR_count_per_cell_zfpkm"
#                 
# )
# sp+stat_cor(
#   aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
#   label.x = 0,label.y = 25
# )
# 
# df<-data.frame(dif$Differentiation,dif$Mean_OR_TPM)
# colnames(df)<-c("Differentiation","Mean_OR_TPM")
# sp <- ggscatter(df, x = "Differentiation", y = "Mean_OR_TPM",
#                 add = "reg.line",  # Add regressin line
#                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#                 conf.int = TRUE,
#                 xlab = "Differentiation",
#                 ylab = "Mean_OR_TPM"
#                 
# )
# sp+stat_cor(
#   aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
#   label.x = -0.1,label.y = 7.5
# )
# 
# df<-data.frame(dif$Differentiation,dif$OR_count_per_cell_zfpkm)
# colnames(df)<-c("Differentiation","OR_count_per_cell_zfpkm")
# sp <- ggscatter(df, x = "Differentiation", y = "OR_count_per_cell_zfpkm",
#                 add = "reg.line",  # Add regressin line
#                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#                 conf.int = TRUE,
#                 xlab = "Differentiation",
#                 ylab = "OR_count_per_cell_zfpkm"
#                 
# )
# sp+stat_cor(
#   aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
#   label.x = 0,label.y = 20
# )


dev.off()

##################LUMINAL-A
########################################Luminal A##################################3
#File1
expression_matrix<-read.csv("LuminalA/LuminalA_expression.csv",row.names = 1)
RN<-row.names(expression_matrix)
expression_matrix<-as.matrix(expression_matrix)

#File2
cell_metadata<-read.csv("LuminalA/LuminalA_Cell_Meta_Data_2.csv",row.names = 1)

#File3
gene_annotation<-read.csv("LuminalA/gene_meta_data.csv")
gene_annotation$genes<-NULL
row.names(gene_annotation)<-RN

#Module1
# Make the CDS object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

pdf(file="Final_pdf_V2/LuminalA.pdf")
#Cluster and classify cells
cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds) + ggtitle("PC Variance")
cds <- reduce_dimension(cds)
cds = cluster_cells(cds, resolution=1e-5)

plot_cells(cds,cell_size = 1.5)
plot_cells(cds, color_cells_by="partition", group_cells_by="partition",cell_size = 1.5)
plot_cells(cds, color_cells_by="OR_count_per_cell_zfpkm",cell_size = 1.5)+ggtitle("OR_count_per_cell_zfpkm")
plot_cells(cds, color_cells_by="Fire_score",cell_size = 1.5)+ggtitle("Fire_score")
plot_cells(cds, color_cells_by="Cell_type",cell_size = 1.5)+ggtitle("Cell_type")
plot_cells(cds, color_cells_by="Molecular_subtype",cell_size = 1.5)+ggtitle("Molecular_subtype")
plot_cells(cds, color_cells_by="Pathological_stage",cell_size = 1.5)+ggtitle("Pathological_stage")
plot_cells(cds, color_cells_by="Immunohistochemistry",cell_size = 1.5)+ggtitle("Immunohistochemistry")
plot_cells(cds, color_cells_by="patient_id",cell_size = 1.5)+ggtitle("patient_id")
plot_cells(cds, color_cells_by="patient_age",cell_size = 1.5)+ggtitle("patient_age")
plot_cells(cds, color_cells_by="Mean_OR_TPM",cell_size = 1.5)+ggtitle("Mean_OR_TPM")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_2",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_2")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_4",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_4")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_5",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_5")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_6",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_6")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_10",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_10")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_15",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_15")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_20",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_20")
plot_cells(cds, color_cells_by="Vector_size",cell_size = 1.5)+ggtitle("Vector Size")

# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="Cell_type",
#                     ordering_type="cluster_row_col",
#                     max.size=3,cell_size = 1.5)+ggtitle("Cell_type")
# 
# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="Pathological_stage",
#                     ordering_type="cluster_row_col",
#                     max.size=3,cell_size = 1.5)+ggtitle("Pathological_stage")
# 
# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="Immunohistochemistry",
#                     ordering_type="cluster_row_col",
#                     max.size=3,cell_size = 1.5)+ggtitle("Immunohistochemistry")

#Module3
cds <- learn_graph(cds)
# plot_cells(cds,
#            color_cells_by = "Pathological_stage",
#            label_groups_by_cluster=FALSE,
#            label_leaves=FALSE,
#            label_branch_points=FALSE,cell_size = 1.5)

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="High_Stemness"){
  cell_ids <- which(colData(cds)[, "State"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
pseudo<-cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
pseudo<-as.data.frame(pseudo)
colnames(pseudo)<-"PseudoTime"
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)

plot_cells(cds,
           color_cells_by = "OR_count_per_cell_zfpkm",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("OR_count_per_cell_zfpkm")

plot_cells(cds,
           color_cells_by = "Fire_score",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Fire_score")

plot_cells(cds,
           color_cells_by = "Cell_type",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Cell_type")

plot_cells(cds,
           color_cells_by = "Vector_size",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Vector_size")

plot_cells(cds,
           color_cells_by = "Pathological_stage",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Pathological_stage")

plot_cells(cds,
           color_cells_by = "Immunohistochemistry",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Immunohistochemistry")

plot_cells(cds,
           color_cells_by = "patient_id",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("patient_id")

plot_cells(cds,
           color_cells_by = "patient_age",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("patient_age")

plot_cells(cds,
           color_cells_by = "Mean_OR_TPM",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Mean_OR_TPM")
plot_cells(cds,
           color_cells_by = "Stemness",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Stemness")
plot_cells(cds,
           color_cells_by = "Differentiation",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Differentiation")
plot_cells(cds,
           color_cells_by = "Molecular_subtype",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Molecular_subtype")


ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))
test_results<-as.data.frame(ciliated_cds_pr_test_res)
test_results<-setDT(test_results,keep.rownames = "Genes")[]
write.table(test_results,"LuminalA/Test_results.csv",row.names = FALSE,sep=",",quote = FALSE)
AFD_genes<-read.csv("/storage/gaurav.ahuja/siddhant/thesis_v2/Monocle_categories/Genes.csv",header=FALSE)
AFD_genes<-as.matrix(AFD_genes)
a<-length(AFD_genes)
for (j in 1:a){
  tryCatch(
    {
      
      AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes[j],
                             colData(cds)$Molecular_subtype %in% c("Luminal A")]
      k<-plot_genes_in_pseudotime(AFD_lineage_cds,
                                  color_cells_by="OR_count_per_cell_zfpkm",
                                  min_expr=0.5,cell_size = 1.5)
      plot(k)
      
    },
    error=function(cond) {
      message(paste("."))
      
    },
    warning=function(cond) {
      message(paste("."))
      
    },
    finally={
      
    }
  )    
}

for (j in 1:a){
  tryCatch(
    {
      
      AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes[j],
                             colData(cds)$Molecular_subtype %in% c("Luminal A")]
      k<-plot_genes_in_pseudotime(AFD_lineage_cds,
                                  color_cells_by="Fire_score",
                                  min_expr=0.5,cell_size = 1.5)
      plot(k)
      
    },
    error=function(cond) {
      message(paste("."))
      
    },
    warning=function(cond) {
      message(paste("."))
      
    },
    finally={
      
    }
  )    
}

for (j in 1:a){
  tryCatch(
    {
      
      AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes[j],
                             colData(cds)$Molecular_subtype %in% c("Luminal A")]
      k<-plot_genes_in_pseudotime(AFD_lineage_cds,
                                  color_cells_by="Mean_OR_TPM",
                                  min_expr=0.5,cell_size = 1.5)
      plot(k)
      
    },
    error=function(cond) {
      message(paste("."))
      
    },
    warning=function(cond) {
      message(paste("."))
      
    },
    finally={
      
    }
  )    
}
Meta_data_pseudo<-cbind(cell_metadata,pseudo)
metadata_with_pesudotime<-setDT(Meta_data_pseudo,keep.rownames = "ID")[]
write.table(metadata_with_pesudotime,file="LuminalA/Metadata_with_pseudotime.csv",sep=",",quote=FALSE,row.names = FALSE)
sp <- ggscatter(Meta_data_pseudo, x = "PseudoTime", y = "Stemness",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                xlab = "PseudoTime",
                ylab = "Stemness"
                
)
sp+stat_cor(
  aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
  label.x = 4,label.y = 0.2
)
sp <- ggscatter(Meta_data_pseudo, x = "PseudoTime", y = "Differentiation",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                xlab = "PseudoTime",
                ylab = "Differentiation"
                
)
sp+stat_cor(
  aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
  label.x = 4,label.y = 0.2
)

sp <- ggscatter(Meta_data_pseudo, x = "PseudoTime", y = "Mean_OR_TPM",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                xlab = "PseudoTime",
                ylab = "Mean_OR_TPM"
                
)
sp+stat_cor(
  aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
  label.x = 4,label.y = 5
)

sp <- ggscatter(Meta_data_pseudo, x = "PseudoTime", y = "Fire_score",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                xlab = "PseudoTime",
                ylab = "Fire_score"
                
)
sp+stat_cor(
  aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
  label.x = 4,label.y = 950
)

sp <- ggscatter(Meta_data_pseudo, x = "PseudoTime", y = "OR_count_per_cell_zfpkm",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                xlab = "PseudoTime",
                ylab = "OR_count_per_cell_zfpkm"
                
)
sp+stat_cor(
  aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
  label.x = 4,label.y = 5
)

ggplot(Meta_data_pseudo, aes(PseudoTime, Stemness)) +geom_bar(position="identity",stat="identity", width = 0.2)
ggplot(Meta_data_pseudo, aes(PseudoTime, Differentiation)) +geom_bar(position="identity",stat="identity", width = 0.2)


newdata<-data.frame(Meta_data_pseudo$OR_count_per_cell_zfpkm,Meta_data_pseudo$Fire_score,Meta_data_pseudo$Mean_OR_TPM,Meta_data_pseudo$Stemness,Meta_data_pseudo$Differentiation,Meta_data_pseudo$PseudoTime)
colnames(newdata)<-c("OR_count","Fire_score","Mean_TPM","Stemness","Differentiation","PseudoTime")

chart.Correlation(newdata, histogram=TRUE, pch=19)

# stem<-cell_metadata[which(cell_metadata$Cell_staus=="Stemness"),]
# dif<-cell_metadata[which(cell_metadata$Cell_staus=="Differentiation"),]
# 
# df<-data.frame(stem$Stemness,stem$Mean_OR_TPM)
# colnames(df)<-c("Stemness","Mean_OR_TPM")
# sp <- ggscatter(df, x = "Stemness", y = "Mean_OR_TPM",
#                 add = "reg.line",  # Add regressin line
#                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#                 conf.int = TRUE,
#                 xlab = "Stemness",
#                 ylab = "Mean_OR_TPM"
#                 
# )
# sp+stat_cor(
#   aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
#   label.x = 0.1,label.y = 1
# )
# 
# df<-data.frame(stem$Stemness,stem$OR_count_per_cell_zfpkm)
# colnames(df)<-c("Stemness","OR_count_per_cell_zfpkm")
# sp <- ggscatter(df, x = "Stemness", y = "OR_count_per_cell_zfpkm",
#                 add = "reg.line",  # Add regressin line
#                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#                 conf.int = TRUE,
#                 xlab = "Stemness",
#                 ylab = "OR_count_per_cell_zfpkm"
#                 
# )
# sp+stat_cor(
#   aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
#   label.x = 0,label.y = 25
# )
# 
# df<-data.frame(dif$Differentiation,dif$Mean_OR_TPM)
# colnames(df)<-c("Differentiation","Mean_OR_TPM")
# sp <- ggscatter(df, x = "Differentiation", y = "Mean_OR_TPM",
#                 add = "reg.line",  # Add regressin line
#                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#                 conf.int = TRUE,
#                 xlab = "Differentiation",
#                 ylab = "Mean_OR_TPM"
#                 
# )
# sp+stat_cor(
#   aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
#   label.x = -0.1,label.y = 7.5
# )
# 
# df<-data.frame(dif$Differentiation,dif$OR_count_per_cell_zfpkm)
# colnames(df)<-c("Differentiation","OR_count_per_cell_zfpkm")
# sp <- ggscatter(df, x = "Differentiation", y = "OR_count_per_cell_zfpkm",
#                 add = "reg.line",  # Add regressin line
#                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#                 conf.int = TRUE,
#                 xlab = "Differentiation",
#                 ylab = "OR_count_per_cell_zfpkm"
#                 
# )
# sp+stat_cor(
#   aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
#   label.x = 0,label.y = 20
# )


dev.off()


##################LUMINAL-b
########################################Luminal b##################################3
#File1
expression_matrix<-read.csv("LuminalB/LuminalB_expression.csv",row.names = 1)
RN<-row.names(expression_matrix)
expression_matrix<-as.matrix(expression_matrix)

#File2
cell_metadata<-read.csv("LuminalB/LuminalB_Cell_Meta_Data_2.csv",row.names = 1)

#File3
gene_annotation<-read.csv("LuminalB/gene_meta_data.csv")
gene_annotation$genes<-NULL
row.names(gene_annotation)<-RN

#Module1
# Make the CDS object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

pdf(file="Final_pdf_V2/LuminalB.pdf")
#Cluster and classify cells
cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds) + ggtitle("PC Variance")
cds <- reduce_dimension(cds)
cds = cluster_cells(cds, resolution=1e-5)

plot_cells(cds,cell_size = 1.5)
plot_cells(cds, color_cells_by="partition", group_cells_by="partition",cell_size = 1.5)
plot_cells(cds, color_cells_by="OR_count_per_cell_zfpkm",cell_size = 1.5)+ggtitle("OR_count_per_cell_zfpkm")
plot_cells(cds, color_cells_by="Fire_score",cell_size = 1.5)+ggtitle("Fire_score")
plot_cells(cds, color_cells_by="Cell_type",cell_size = 1.5)+ggtitle("Cell_type")
plot_cells(cds, color_cells_by="Molecular_subtype",cell_size = 1.5)+ggtitle("Molecular_subtype")
plot_cells(cds, color_cells_by="Pathological_stage",cell_size = 1.5)+ggtitle("Pathological_stage")
plot_cells(cds, color_cells_by="Immunohistochemistry",cell_size = 1.5)+ggtitle("Immunohistochemistry")
plot_cells(cds, color_cells_by="patient_id",cell_size = 1.5)+ggtitle("patient_id")
plot_cells(cds, color_cells_by="patient_age",cell_size = 1.5)+ggtitle("patient_age")
plot_cells(cds, color_cells_by="Mean_OR_TPM",cell_size = 1.5)+ggtitle("Mean_OR_TPM")
#plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_2",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_2")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_4",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_4")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_5",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_5")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_6",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_6")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_10",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_10")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_15",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_15")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_20",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_20")
plot_cells(cds, color_cells_by="Vector_size",cell_size = 1.5)+ggtitle("Vector Size")
# 
# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="Cell_type",
#                     ordering_type="cluster_row_col",
#                     max.size=3,cell_size = 1.5)+ggtitle("Cell_type")
# 
# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="Pathological_stage",
#                     ordering_type="cluster_row_col",
#                     max.size=3,cell_size = 1.5)+ggtitle("Pathological_stage")
# 
# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="Immunohistochemistry",
#                     ordering_type="cluster_row_col",
#                     max.size=3,cell_size = 1.5)+ggtitle("Immunohistochemistry")

#Module3
cds <- learn_graph(cds)
# plot_cells(cds,
#            color_cells_by = "Pathological_stage",
#            label_groups_by_cluster=FALSE,
#            label_leaves=FALSE,
#            label_branch_points=FALSE,cell_size = 1.5)

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="High_Stemness"){
  cell_ids <- which(colData(cds)[, "State"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
pseudo<-cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
pseudo<-as.data.frame(pseudo)
colnames(pseudo)<-"PseudoTime"
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)

plot_cells(cds,
           color_cells_by = "OR_count_per_cell_zfpkm",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("OR_count_per_cell_zfpkm")

plot_cells(cds,
           color_cells_by = "Fire_score",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Fire_score")

plot_cells(cds,
           color_cells_by = "Cell_type",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Cell_type")

plot_cells(cds,
           color_cells_by = "Vector_size",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Vector_size")

plot_cells(cds,
           color_cells_by = "Pathological_stage",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Pathological_stage")

plot_cells(cds,
           color_cells_by = "Immunohistochemistry",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Immunohistochemistry")

plot_cells(cds,
           color_cells_by = "patient_id",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("patient_id")

plot_cells(cds,
           color_cells_by = "patient_age",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("patient_age")

plot_cells(cds,
           color_cells_by = "Mean_OR_TPM",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Mean_OR_TPM")
plot_cells(cds,
           color_cells_by = "Stemness",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Stemness")
plot_cells(cds,
           color_cells_by = "Differentiation",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Differentiation")
plot_cells(cds,
           color_cells_by = "Molecular_subtype",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Molecular_subtype")


ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))
test_results<-as.data.frame(ciliated_cds_pr_test_res)
test_results<-setDT(test_results,keep.rownames = "Genes")[]
write.table(test_results,"LuminalB/Test_results.csv",row.names = FALSE,sep=",",quote = FALSE)
AFD_genes<-read.csv("/storage/gaurav.ahuja/siddhant/thesis_v2/Monocle_categories/Genes.csv",header=FALSE)
AFD_genes<-as.matrix(AFD_genes)
a<-length(AFD_genes)
for (j in 1:a){
  tryCatch(
    {
      
      AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes[j],
                             colData(cds)$Molecular_subtype %in% c("Luminal B")]
      k<-plot_genes_in_pseudotime(AFD_lineage_cds,
                                  color_cells_by="OR_count_per_cell_zfpkm",
                                  min_expr=0.5,cell_size = 1.5)
      plot(k)
      
    },
    error=function(cond) {
      message(paste("."))
      
    },
    warning=function(cond) {
      message(paste("."))
      
    },
    finally={
      
    }
  )    
}

for (j in 1:a){
  tryCatch(
    {
      
      AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes[j],
                             colData(cds)$Molecular_subtype %in% c("Luminal B")]
      k<-plot_genes_in_pseudotime(AFD_lineage_cds,
                                  color_cells_by="Fire_score",
                                  min_expr=0.5,cell_size = 1.5)
      plot(k)
      
    },
    error=function(cond) {
      message(paste("."))
      
    },
    warning=function(cond) {
      message(paste("."))
      
    },
    finally={
      
    }
  )    
}

for (j in 1:a){
  tryCatch(
    {
      
      AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes[j],
                             colData(cds)$Molecular_subtype %in% c("Luminal B")]
      k<-plot_genes_in_pseudotime(AFD_lineage_cds,
                                  color_cells_by="Mean_OR_TPM",
                                  min_expr=0.5,cell_size = 1.5)
      plot(k)
      
    },
    error=function(cond) {
      message(paste("."))
      
    },
    warning=function(cond) {
      message(paste("."))
      
    },
    finally={
      
    }
  )    
}
Meta_data_pseudo<-cbind(cell_metadata,pseudo)
metadata_with_pesudotime<-setDT(Meta_data_pseudo,keep.rownames = "ID")[]
write.table(metadata_with_pesudotime,file="LuminalB/Metadata_with_pseudotime.csv",sep=",",quote=FALSE,row.names = FALSE)
sp <- ggscatter(Meta_data_pseudo, x = "PseudoTime", y = "Stemness",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                xlab = "PseudoTime",
                ylab = "Stemness"
                
)
sp+stat_cor(
  aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
  label.x = 0.2,label.y = 0.2
)
sp <- ggscatter(Meta_data_pseudo, x = "PseudoTime", y = "Differentiation",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                xlab = "PseudoTime",
                ylab = "Differentiation"
                
)
sp+stat_cor(
  aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
  label.x = 0.2,label.y = 0.2
)

sp <- ggscatter(Meta_data_pseudo, x = "PseudoTime", y = "Mean_OR_TPM",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                xlab = "PseudoTime",
                ylab = "Mean_OR_TPM"
                
)
sp+stat_cor(
  aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
  label.x = 0.4,label.y = 5
)

sp <- ggscatter(Meta_data_pseudo, x = "PseudoTime", y = "Fire_score",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                xlab = "PseudoTime",
                ylab = "Fire_score"
                
)
sp+stat_cor(
  aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
  label.x = 0.4,label.y = 950
)

sp <- ggscatter(Meta_data_pseudo, x = "PseudoTime", y = "OR_count_per_cell_zfpkm",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                xlab = "PseudoTime",
                ylab = "OR_count_per_cell_zfpkm"
                
)
sp+stat_cor(
  aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
  label.x = 0.4,label.y = 5
)

ggplot(Meta_data_pseudo, aes(PseudoTime, Stemness)) +geom_bar(position="identity",stat="identity", width = 0.2)
ggplot(Meta_data_pseudo, aes(PseudoTime, Differentiation)) +geom_bar(position="identity",stat="identity", width = 0.2)


newdata<-data.frame(Meta_data_pseudo$OR_count_per_cell_zfpkm,Meta_data_pseudo$Fire_score,Meta_data_pseudo$Mean_OR_TPM,Meta_data_pseudo$Stemness,Meta_data_pseudo$Differentiation,Meta_data_pseudo$PseudoTime)
colnames(newdata)<-c("OR_count","Fire_score","Mean_TPM","Stemness","Differentiation","PseudoTime")

chart.Correlation(newdata, histogram=TRUE, pch=19)

# stem<-cell_metadata[which(cell_metadata$Cell_staus=="Stemness"),]
# dif<-cell_metadata[which(cell_metadata$Cell_staus=="Differentiation"),]
# 
# df<-data.frame(stem$Stemness,stem$Mean_OR_TPM)
# colnames(df)<-c("Stemness","Mean_OR_TPM")
# sp <- ggscatter(df, x = "Stemness", y = "Mean_OR_TPM",
#                 add = "reg.line",  # Add regressin line
#                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#                 conf.int = TRUE,
#                 xlab = "Stemness",
#                 ylab = "Mean_OR_TPM"
#                 
# )
# sp+stat_cor(
#   aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
#   label.x = 0.1,label.y = 1
# )
# 
# df<-data.frame(stem$Stemness,stem$OR_count_per_cell_zfpkm)
# colnames(df)<-c("Stemness","OR_count_per_cell_zfpkm")
# sp <- ggscatter(df, x = "Stemness", y = "OR_count_per_cell_zfpkm",
#                 add = "reg.line",  # Add regressin line
#                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#                 conf.int = TRUE,
#                 xlab = "Stemness",
#                 ylab = "OR_count_per_cell_zfpkm"
#                 
# )
# sp+stat_cor(
#   aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
#   label.x = 0,label.y = 25
# )
# 
# df<-data.frame(dif$Differentiation,dif$Mean_OR_TPM)
# colnames(df)<-c("Differentiation","Mean_OR_TPM")
# sp <- ggscatter(df, x = "Differentiation", y = "Mean_OR_TPM",
#                 add = "reg.line",  # Add regressin line
#                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#                 conf.int = TRUE,
#                 xlab = "Differentiation",
#                 ylab = "Mean_OR_TPM"
#                 
# )
# sp+stat_cor(
#   aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
#   label.x = -0.1,label.y = 7.5
# )
# 
# df<-data.frame(dif$Differentiation,dif$OR_count_per_cell_zfpkm)
# colnames(df)<-c("Differentiation","OR_count_per_cell_zfpkm")
# sp <- ggscatter(df, x = "Differentiation", y = "OR_count_per_cell_zfpkm",
#                 add = "reg.line",  # Add regressin line
#                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#                 conf.int = TRUE,
#                 xlab = "Differentiation",
#                 ylab = "OR_count_per_cell_zfpkm"
#                 
# )
# sp+stat_cor(
#   aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
#   label.x = 0,label.y = 20
# )


dev.off()

#Combined
setwd("/storage/gaurav.ahuja/siddhant/thesis_v2/MOCOCLE/updated/updated_bc05_removed")
#File1
expression_matrix<-read.csv("expression_matrix.csv",row.names = 1)
expression_matrix<-as.matrix(expression_matrix)

#File2
cell_metadata<-read.csv("Cell_Meta_Data_2.csv",sep=",",row.names = 1)

#File3
gene_annotation<-read.csv("gene_meta_data.csv",row.names = 1)

#Module1
# Make the CDS object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
pdf(file="/storage/gaurav.ahuja/siddhant/thesis_v2/Monocle_categories/Final_pdf_V2/Combined.pdf")
#Cluster and classify cells
cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds) + ggtitle("PC Variance")
cds <- reduce_dimension(cds)
cds = cluster_cells(cds, resolution=1e-5)

plot_cells(cds,cell_size = 1.5)
plot_cells(cds, color_cells_by="partition", group_cells_by="partition",cell_size = 1.5)
plot_cells(cds, color_cells_by="OR_count_per_cell_zfpkm",cell_size = 1.5)+ggtitle("OR_count_per_cell_zfpkm")
plot_cells(cds, color_cells_by="Fire_score",cell_size = 1.5)+ggtitle("Fire_score")
plot_cells(cds, color_cells_by="Cell_type",cell_size = 1.5)+ggtitle("Cell_type")
plot_cells(cds, color_cells_by="Molecular_subtype",cell_size = 1.5)+ggtitle("Molecular_subtype")
plot_cells(cds, color_cells_by="Pathological_stage",cell_size = 1.5)+ggtitle("Pathological_stage")
plot_cells(cds, color_cells_by="Immunohistochemistry",cell_size = 1.5)+ggtitle("Immunohistochemistry")
plot_cells(cds, color_cells_by="patient_id",cell_size = 1.5)+ggtitle("patient_id")
plot_cells(cds, color_cells_by="patient_age",cell_size = 1.5)+ggtitle("patient_age")
plot_cells(cds, color_cells_by="Mean_OR_TPM",cell_size = 1.5)+ggtitle("Mean_OR_TPM")
#plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_2",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_2")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_4",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_4")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_5",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_5")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_6",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_6")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_10",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_10")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_15",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_15")
# plot_cells(cds, color_cells_by="OR_count_TPM_cutoff_20",cell_size = 1.5)+ggtitle("OR_count_TPM_cutoff_20")
#plot_cells(cds, color_cells_by="Vector_size",cell_size = 1.5)+ggtitle("Vector Size")
# 
# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="Cell_type",
#                     ordering_type="cluster_row_col",
#                     max.size=3,cell_size = 1.5)+ggtitle("Cell_type")
# 
# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="Pathological_stage",
#                     ordering_type="cluster_row_col",
#                     max.size=3,cell_size = 1.5)+ggtitle("Pathological_stage")
# 
# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="Immunohistochemistry",
#                     ordering_type="cluster_row_col",
#                     max.size=3,cell_size = 1.5)+ggtitle("Immunohistochemistry")

#Module3
cds <- learn_graph(cds)
# plot_cells(cds,
#            color_cells_by = "Pathological_stage",
#            label_groups_by_cluster=FALSE,
#            label_leaves=FALSE,
#            label_branch_points=FALSE,cell_size = 1.5)

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="High_Stemness"){
  cell_ids <- which(colData(cds)[, "State"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
pseudo<-cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
pseudo<-as.data.frame(pseudo)
colnames(pseudo)<-"PseudoTime"
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)

plot_cells(cds,
           color_cells_by = "OR_count_per_cell_zfpkm",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("OR_count_per_cell_zfpkm")

plot_cells(cds,
           color_cells_by = "Fire_score",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Fire_score")

plot_cells(cds,
           color_cells_by = "Cell_type",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Cell_type")

# plot_cells(cds,
#            color_cells_by = "Vector_size",
#            label_cell_groups=FALSE,
#            label_leaves=FALSE,
#            label_branch_points=FALSE,
#            graph_label_size=1.5,cell_size = 1.5)+ggtitle("Vector_size")

plot_cells(cds,
           color_cells_by = "Pathological_stage",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Pathological_stage")

plot_cells(cds,
           color_cells_by = "Immunohistochemistry",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Immunohistochemistry")

plot_cells(cds,
           color_cells_by = "patient_id",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("patient_id")

plot_cells(cds,
           color_cells_by = "patient_age",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("patient_age")

plot_cells(cds,
           color_cells_by = "Mean_OR_TPM",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Mean_OR_TPM")
plot_cells(cds,
           color_cells_by = "Stemness",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Stemness")
plot_cells(cds,
           color_cells_by = "Differentiation",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Differentiation")
plot_cells(cds,
           color_cells_by = "Molecular_subtype",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Molecular_subtype")


ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))
test_results<-as.data.frame(ciliated_cds_pr_test_res)
test_results<-setDT(test_results,keep.rownames = "Genes")[]
write.table(test_results,"Test_results.csv",row.names = FALSE,sep=",",quote = FALSE)
AFD_genes<-read.csv("/storage/gaurav.ahuja/siddhant/thesis_v2/Monocle_categories/Genes.csv",header=FALSE)
AFD_genes<-as.matrix(AFD_genes)
a<-length(AFD_genes)
for (j in 1:a){
  tryCatch(
    {
      
      AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes[j],
                             colData(cds)$Molecular_subtype %in% c("Basal-like","Luminal A","Luminal B","HER2-enriched")]
      k<-plot_genes_in_pseudotime(AFD_lineage_cds,
                                  color_cells_by="Molecular_subtype",
                                  min_expr=0.5,cell_size = 1.5)
      plot(k)
      
    },
    error=function(cond) {
      message(paste("."))
      
    },
    warning=function(cond) {
      message(paste("."))
      
    },
    finally={
      
    }
  )    
}

for (j in 1:a){
  tryCatch(
    {
      
      AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes[j],
                             colData(cds)$Molecular_subtype %in% c("Basal-like","Luminal A","Luminal B","HER2-enriched")]
      k<-plot_genes_in_pseudotime(AFD_lineage_cds,
                                  color_cells_by="OR_count_per_cell_zfpkm",
                                  min_expr=0.5,cell_size = 1.5)
      plot(k)
      
    },
    error=function(cond) {
      message(paste("."))
      
    },
    warning=function(cond) {
      message(paste("."))
      
    },
    finally={
      
    }
  )    
}

for (j in 1:a){
  tryCatch(
    {
      
      AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes[j],
                             colData(cds)$Molecular_subtype %in% c("Basal-like","Luminal A","Luminal B","HER2-enriched")]
      k<-plot_genes_in_pseudotime(AFD_lineage_cds,
                                  color_cells_by="Fire_score",
                                  min_expr=0.5,cell_size = 1.5)
      plot(k)
      
    },
    error=function(cond) {
      message(paste("."))
      
    },
    warning=function(cond) {
      message(paste("."))
      
    },
    finally={
      
    }
  )    
}

for (j in 1:a){
  tryCatch(
    {
      
      AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes[j],
                             colData(cds)$Molecular_subtype %in% c("Basal-like","Luminal A","Luminal B","HER2-enriched")]
      k<-plot_genes_in_pseudotime(AFD_lineage_cds,
                                  color_cells_by="Mean_OR_TPM",
                                  min_expr=0.5,cell_size = 1.5)
      plot(k)
      
    },
    error=function(cond) {
      message(paste("."))
      
    },
    warning=function(cond) {
      message(paste("."))
      
    },
    finally={
      
    }
  )    
}
Meta_data_pseudo<-cbind(cell_metadata,pseudo)
metadata_with_pesudotime<-setDT(Meta_data_pseudo,keep.rownames = "ID")[]
write.table(metadata_with_pesudotime,file="Metadata_with_pseudotime.csv",sep=",",quote=FALSE,row.names = FALSE)
sp <- ggscatter(Meta_data_pseudo, x = "PseudoTime", y = "Stemness",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                xlab = "PseudoTime",
                ylab = "Stemness"
                
)
sp+stat_cor(
  aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
  label.x = 0.2,label.y = 0.2
)
sp <- ggscatter(Meta_data_pseudo, x = "PseudoTime", y = "Differentiation",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                xlab = "PseudoTime",
                ylab = "Differentiation"
                
)
sp+stat_cor(
  aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
  label.x = 0.2,label.y = 0.2
)

sp <- ggscatter(Meta_data_pseudo, x = "PseudoTime", y = "Mean_OR_TPM",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                xlab = "PseudoTime",
                ylab = "Mean_OR_TPM"
                
)
sp+stat_cor(
  aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
  label.x = 0.4,label.y = 8.5
)

sp <- ggscatter(Meta_data_pseudo, x = "PseudoTime", y = "Fire_score",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                xlab = "PseudoTime",
                ylab = "Fire_score"
                
)
sp+stat_cor(
  aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
  label.x = 0.4,label.y = 950
)

sp <- ggscatter(Meta_data_pseudo, x = "PseudoTime", y = "OR_count_per_cell_zfpkm",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                xlab = "PseudoTime",
                ylab = "OR_count_per_cell_zfpkm"
                
)
sp+stat_cor(
  aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
  label.x = 0.4,label.y = 5
)

ggplot(Meta_data_pseudo, aes(PseudoTime, Stemness)) +geom_bar(position="identity",stat="identity", width = 0.2)
ggplot(Meta_data_pseudo, aes(PseudoTime, Differentiation)) +geom_bar(position="identity",stat="identity", width = 0.2)


newdata<-data.frame(Meta_data_pseudo$OR_count_per_cell_zfpkm,Meta_data_pseudo$Fire_score,Meta_data_pseudo$Mean_OR_TPM,Meta_data_pseudo$Stemness,Meta_data_pseudo$Differentiation,Meta_data_pseudo$PseudoTime)
colnames(newdata)<-c("OR_count","Fire_score","Mean_TPM","Stemness","Differentiation","PseudoTime")

chart.Correlation(newdata, histogram=TRUE, pch=19)

# stem<-cell_metadata[which(cell_metadata$Cell_staus=="Stemness"),]
# dif<-cell_metadata[which(cell_metadata$Cell_staus=="Differentiation"),]
# 
# df<-data.frame(stem$Stemness,stem$Mean_OR_TPM)
# colnames(df)<-c("Stemness","Mean_OR_TPM")
# sp <- ggscatter(df, x = "Stemness", y = "Mean_OR_TPM",
#                 add = "reg.line",  # Add regressin line
#                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#                 conf.int = TRUE,
#                 xlab = "Stemness",
#                 ylab = "Mean_OR_TPM"
#                 
# )
# sp+stat_cor(
#   aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
#   label.x = 0.1,label.y = 1
# )
# 
# df<-data.frame(stem$Stemness,stem$OR_count_per_cell_zfpkm)
# colnames(df)<-c("Stemness","OR_count_per_cell_zfpkm")
# sp <- ggscatter(df, x = "Stemness", y = "OR_count_per_cell_zfpkm",
#                 add = "reg.line",  # Add regressin line
#                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#                 conf.int = TRUE,
#                 xlab = "Stemness",
#                 ylab = "OR_count_per_cell_zfpkm"
#                 
# )
# sp+stat_cor(
#   aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
#   label.x = 0,label.y = 25
# )
# 
# df<-data.frame(dif$Differentiation,dif$Mean_OR_TPM)
# colnames(df)<-c("Differentiation","Mean_OR_TPM")
# sp <- ggscatter(df, x = "Differentiation", y = "Mean_OR_TPM",
#                 add = "reg.line",  # Add regressin line
#                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#                 conf.int = TRUE,
#                 xlab = "Differentiation",
#                 ylab = "Mean_OR_TPM"
#                 
# )
# sp+stat_cor(
#   aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
#   label.x = -0.1,label.y = 7.5
# )
# 
# df<-data.frame(dif$Differentiation,dif$OR_count_per_cell_zfpkm)
# colnames(df)<-c("Differentiation","OR_count_per_cell_zfpkm")
# sp <- ggscatter(df, x = "Differentiation", y = "OR_count_per_cell_zfpkm",
#                 add = "reg.line",  # Add regressin line
#                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#                 conf.int = TRUE,
#                 xlab = "Differentiation",
#                 ylab = "OR_count_per_cell_zfpkm"
#                 
# )
# sp+stat_cor(
#   aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
#   label.x = 0,label.y = 20
# )


dev.off()


##MONOCole
library(ggplot2)
library(dplyr)
library(monocle3)
library("PerformanceAnalytics")
setwd("/storage/gaurav.ahuja/siddhant/thesis_v2/Monocle_categories/Normal/Inputs")
#File1
expression_matrix<-read.csv("Expression.csv",row.names = 1)
expression_matrix<-as.matrix(expression_matrix)

#File2
cell_metadata<-read.csv("Cell_meta_data.csv",sep=",",row.names = 1)

#File3
gene_annotation<-read.csv("gene_metadata.csv",row.names = 1)

# Make the CDS object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

#pdf


cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds) + ggtitle("PC Variance")
cds <- reduce_dimension(cds)
cds = cluster_cells(cds, resolution=1e-5)
plot_cells(cds,cell_size = 1.5)
plot_cells(cds, color_cells_by="partition", group_cells_by="partition",cell_size = 1.5)
plot_cells(cds, color_cells_by="count_value",cell_size = 1.5)+ggtitle("OR_count_per_cell_zfpkm")
plot_cells(cds, color_cells_by="Mean_TPM",cell_size = 1.5)+ggtitle("Mean TPM")
plot_cells(cds, color_cells_by="sampling_site",cell_size = 1.5)+ggtitle("Sampling site")
plot_cells(cds, color_cells_by="individual",cell_size = 1.5)+ggtitle("Individual")
plot_cells(cds, color_cells_by="Sample_title",cell_size = 1.5)+ggtitle("Sample title")
cds <- learn_graph(cds)

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="High_Stemness"){
  cell_ids <- which(colData(cds)[, "State"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
pseudo<-cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
pseudo<-as.data.frame(pseudo)
colnames(pseudo)<-"PseudoTime"
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)

plot_cells(cds,
           color_cells_by = "count_value",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("OR_count_per_cell_zfpkm")

plot_cells(cds,
           color_cells_by = "Mean_TPM",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Mean TPM")

plot_cells(cds,
           color_cells_by = "sampling_site",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Sampling Site")


plot_cells(cds,
           color_cells_by = "individual",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Individual")

plot_cells(cds,
           color_cells_by = "Stemness",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Stemness")
plot_cells(cds,
           color_cells_by = "Differentiation",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,cell_size = 1.5)+ggtitle("Differentiation")

ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))
test_results<-as.data.frame(ciliated_cds_pr_test_res)
test_results<-setDT(test_results,keep.rownames = "Genes")[]
write.table(test_results,"Test_results.csv",row.names = FALSE,sep=",",quote = FALSE)
AFD_genes<-read.csv("genes.csv",header=FALSE)
AFD_genes<-as.matrix(AFD_genes)
a<-length(AFD_genes)

for (j in 1:a){
  tryCatch(
    {
      
      AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes[j],
                             colData(cds)$sampling_site %in% c("luminal")]
      k<-plot_genes_in_pseudotime(AFD_lineage_cds,
                                  color_cells_by="Mean_TPM",
                                  min_expr=0.5,cell_size = 1.5)
      plot(k)
      
    },
    error=function(cond) {
      message(paste("."))
      
    },
    warning=function(cond) {
      message(paste("."))
      
    },
    finally={
      
    }
  )    
}
Meta_data_pseudo<-cbind(cell_metadata,pseudo)
metadata_with_pesudotime<-setDT(Meta_data_pseudo,keep.rownames = "ID")[]
write.table(metadata_with_pesudotime,file="Metadata_with_pseudotime.csv",sep=",",quote=FALSE,row.names = FALSE)
sp <- ggscatter(Meta_data_pseudo, x = "PseudoTime", y = "Stemness",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                xlab = "PseudoTime",
                ylab = "Stemness"
                
)
sp+stat_cor(
  aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
  label.x = 0.2,label.y = 0.2
)
sp <- ggscatter(Meta_data_pseudo, x = "PseudoTime", y = "Differentiation",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                xlab = "PseudoTime",
                ylab = "Differentiation"
                
)
sp+stat_cor(
  aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
  label.x = 0.2,label.y = 0.2
)

sp <- ggscatter(Meta_data_pseudo, x = "PseudoTime", y = "Mean_TPM",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                xlab = "PseudoTime",
                ylab = "Mean_TPM"
                
)
sp+stat_cor(
  aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
  label.x = 0.4,label.y = 8.5
)

# sp <- ggscatter(Meta_data_pseudo, x = "PseudoTime", y = "Fire_score",
#                 add = "reg.line",  # Add regressin line
#                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#                 conf.int = TRUE,
#                 xlab = "PseudoTime",
#                 ylab = "Fire_score"
#                 
# )
# sp+stat_cor(
#   aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
#   label.x = 0.4,label.y = 950
# )

sp <- ggscatter(Meta_data_pseudo, x = "PseudoTime", y = "count_value",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE,
                xlab = "PseudoTime",
                ylab = "OR_count_per_cell_zfpkm"
                
)
sp+stat_cor(
  aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
  label.x = 0.4,label.y = 5
)

ggplot(Meta_data_pseudo, aes(PseudoTime, Stemness)) +geom_bar(position="identity",stat="identity", width = 0.2)
ggplot(Meta_data_pseudo, aes(PseudoTime, Differentiation)) +geom_bar(position="identity",stat="identity", width = 0.2)


newdata<-data.frame(Meta_data_pseudo$count_value,Meta_data_pseudo$Mean_TPM,Meta_data_pseudo$Stemness,Meta_data_pseudo$Differentiation,Meta_data_pseudo$PseudoTime)
colnames(newdata)<-c("OR_count","Mean_TPM","Stemness","Differentiation","PseudoTime")

chart.Correlation(newdata, histogram=TRUE, pch=19)

# stem<-cell_metadata[which(cell_metadata$Cell_staus=="Stemness"),]
# dif<-cell_metadata[which(cell_metadata$Cell_staus=="Differentiation"),]
# 
# df<-data.frame(stem$Stemness,stem$Mean_OR_TPM)
# colnames(df)<-c("Stemness","Mean_OR_TPM")
# sp <- ggscatter(df, x = "Stemness", y = "Mean_OR_TPM",
#                 add = "reg.line",  # Add regressin line
#                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#                 conf.int = TRUE,
#                 xlab = "Stemness",
#                 ylab = "Mean_OR_TPM"
#                 
# )
# sp+stat_cor(
#   aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
#   label.x = 0.1,label.y = 1
# )
# 
# df<-data.frame(stem$Stemness,stem$OR_count_per_cell_zfpkm)
# colnames(df)<-c("Stemness","OR_count_per_cell_zfpkm")
# sp <- ggscatter(df, x = "Stemness", y = "OR_count_per_cell_zfpkm",
#                 add = "reg.line",  # Add regressin line
#                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#                 conf.int = TRUE,
#                 xlab = "Stemness",
#                 ylab = "OR_count_per_cell_zfpkm"
#                 
# )
# sp+stat_cor(
#   aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
#   label.x = 0,label.y = 25
# )
# 
# df<-data.frame(dif$Differentiation,dif$Mean_OR_TPM)
# colnames(df)<-c("Differentiation","Mean_OR_TPM")
# sp <- ggscatter(df, x = "Differentiation", y = "Mean_OR_TPM",
#                 add = "reg.line",  # Add regressin line
#                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#                 conf.int = TRUE,
#                 xlab = "Differentiation",
#                 ylab = "Mean_OR_TPM"
#                 
# )
# sp+stat_cor(
#   aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
#   label.x = -0.1,label.y = 7.5
# )
# 
# df<-data.frame(dif$Differentiation,dif$OR_count_per_cell_zfpkm)
# colnames(df)<-c("Differentiation","OR_count_per_cell_zfpkm")
# sp <- ggscatter(df, x = "Differentiation", y = "OR_count_per_cell_zfpkm",
#                 add = "reg.line",  # Add regressin line
#                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#                 conf.int = TRUE,
#                 xlab = "Differentiation",
#                 ylab = "OR_count_per_cell_zfpkm"
#                 
# )
# sp+stat_cor(
#   aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
#   label.x = 0,label.y = 20
# )


dev.off()
