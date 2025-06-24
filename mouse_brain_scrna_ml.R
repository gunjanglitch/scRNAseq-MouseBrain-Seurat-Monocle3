install.packages("Matrix")
install.packages("randomForest")
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3')
devtools::install_github("satijalab/seurat-wrappers")
#load packages
#core analaysis package
library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(monocle3)
library(SingleR)
library(celldex)
library(cluster)
library(randomForest)
library(SeuratWrappers)

#setting up directory
getwd()
setwd("C:/Users/Asus/Desktop/mouse_brain")

#load data 
data_path = "C:/Users/Asus/Desktop/mouse_brain"
sc_data = Read10X(data.dir = data_path)

#creating seurat object
seurat_obj = CreateSeuratObject(counts = sc_data, project = "mousebrain")

#quality control
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
qc_plot = VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
qc_plot
ggsave("qc_voilin_plot.png", qc_plot ,width = 12, height = 8)

#before filtering
seurat_obj_raw = CreateSeuratObject(counts = sc_data, project = "mousebrain")
original_cell_count = ncol(seurat_obj_raw)

#filter low quality cells
seurat_obj = subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
write.csv(seurat_obj@meta.data, file = "filtered_cells_metadata.csv")

#normalization
seurat_obj = NormalizeData(seurat_obj)
norm_data = GetAssayData(seurat_obj, layer = "data")
norm_data_df = as.data.frame(as.matrix(norm_data))
write.csv(norm_data_df, file = "normalized_counts.csv")

#feature selection
seurat_obj = FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
var_genes = VariableFeatures(seurat_obj)
var_genes_data = norm_data[var_genes, ]
write.csv(as.data.frame(as.matrix(var_genes_data)), file = "variable_features.csv")

#scaling
seurat_obj = ScaleData(seurat_obj)
scaled_data = GetAssayData(seurat_obj, slot = "scale.data")
scaled_data = as.data.frame(as.matrix(scaled_data))
write.csv(scaled_data, file = "scaled_data.csv")
scaled_var_genes = scaled_data[var_genes, ]
write.csv(as.data.frame(as.matrix(scaled_var_genes)), file = "scaled_variable_features.csv")

#PCA
seurat_obj = RunPCA(seurat_obj)
elbow_pca = ElbowPlot(seurat_obj)
elbow_pca
ggsave("pca_elbow_plot.jpg", elbow_pca, width = 12, height = 8)

#dimensionality reduction
seurat_obj = RunUMAP(seurat_obj, dims = 1:10)
umap_plot = UMAPPlot(seurat_obj)
umap_plot
ggsave("umap_plot.jpg", umap_plot, width = 12, height = 8)

#clustering
seurat_obj = FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj = FindClusters(seurat_obj, resolution = 1.2)
dim_plot = DimPlot(seurat_obj, reduction = "umap", label = TRUE)
dim_plot
ggsave("dim_plot.png", dim_plot, width = 12, height = 8, dpi = 300)

table(Idents(seurat_obj))

#marker gene identification
markers = FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25)
write.csv(markers, "markers_mousebrain.csv")

#trajectory analysis with monocle3
cds = as.cell_data_set(seurat_obj) #convert seurat to monocle3 celldataset
cds = preprocess_cds(cds, num_dim = 10) #preprocess the data
cds = reduce_dimension(cds) #reduce dimensionality using UMAP
cds = cluster_cells(cds) #cluster cells (optional but improves graph)
cds = learn_graph(cds) #learn trajectory graph
cds = order_cells(cds) #order ccells in pseudotime
traj_ana_plot = plot_cells(cds, color_cells_by = "pseudotime")
traj_ana_plot
ggsave("trajectory_analysis_plot.png", traj_ana_plot, width = 12, height = 8, dpi = 300)

#differential expression between clusters
Idents(seurat_obj) = "seurat_clusters"
deg_0_vs_1 = FindMarkers(seurat_obj, ident.1 = 0, ident.2 = 1)
write.csv(deg_0_vs_1, "deg_clusters0_vs_1.csv")

#cell type annotation
ref = MouseRNAseqData()
seurat_counts = GetAssayData(seurat_obj, slot = "data")
singleR_result = SingleR(test = seurat_counts, ref = ref, labels = ref$label.main)
seurat_obj$SingleR = singleR_result$labels
cell_anno_plot = DimPlot(seurat_obj, group.by = "SingleR", label = TRUE)
cell_anno_plot
ggsave("cell_annotation_plot.png", cell_anno_plot, width = 12, height = 8)


#ML : RandomForest for cell type prediction
gene_data = t(as.matrix(GetAssayData(seurat_obj, slot = "scale.data")))
labels = seurat_obj$SingleR

set.seed(123)
train_idx = sample(1:nrow(gene_data), size = 0.7 * nrow(gene_data))
train_data = gene_data[train_idx, ]
test_data = gene_data[-train_idx, ]
train_lables = labels[train_idx]
test_lables = labels[-train_idx]

rf_model = randomForest(x = train_data, y = as.factor(train_lables), ntree = 100)
predictions = predict(rf_model, test_data)

conf_mat = table(predictions, test_lables)
accuracy = sum(diag(conf_mat)) / sum(conf_mat)
print(paste("Random Forest Accuracy:", round(accuracy, 3)))

cluster_labels = Idents(seurat_obj)
dist_matrix = dist(GetAssayData(seurat_obj, slot = "scale.data")[1:500, ])
sil = silhouette(as.numeric(cluster_labels)[1:500], dist_matrix)
ml_plot = plot(sil)
