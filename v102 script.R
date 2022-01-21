#SF7 v102 genome
setwd("G:/Analysis done on lab computer/single cell annotations paper/re analysis on realignment")

#load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(clustree)

#load in sf7 data
#read sf7
v102_sf7 <- Read10X(data.dir = "E:/Single Cell all data/work/pachon_SF7/outs/filtered_feature_bc_matrix")

#change ensembl codes in the matrix #change gene names in the expression matrix
dim(v102_sf7)
25489  8870

#create sf7 based on no.of genes and cells
v102_sf7 <- CreateSeuratObject(counts = v102_sf7, min.cells = 2, min.features = 50, project="SFuninjured_v102")

plot1 <- FeatureScatter(v102_sf7, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 

#draw scatter plot  intercepts to look at filtering
plot1 + geom_hline(yintercept=200) + geom_vline(xintercept = 20000) + geom_hline(yintercept=2500)

#check the Violin plots
VlnPlot(v102_sf7, features = c("nFeature_RNA", "nCount_RNA"))

#filter based on no. of features/genes and counts/UMIs
#use violin plots with histograms to set thresholds
v102_sf7 <- subset(x = v102_sf7, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & nCount_RNA <20000)

#check how many cells left after filtering
table(Idents(v102_sf7))
SFuninjured_v102 
8717 

dim(v102_sf7)
[1] 16408  8717

#normalise data, identify highly variable features and scale data
v102_sf7 <- SCTransform(v102_sf7, verbose = T)      

#run PCA
v102_sf7 <- RunPCA(v102_sf7, features=VariableFeatures(v102_sf7))

#visualise PCs
DimHeatmap(v102_sf7, dims = 1:15, cells = 500, balanced = T)
DimHeatmap(v102_sf7, dims = 15:30, cells = 500, balanced = T)
DimHeatmap(v102_sf7, dims = 30:45, cells = 500, balanced = T)
ElbowPlot(v102_sf7, ndims=50)

#run UMAP
#chose dimensions of 30
v102_sf7 <- RunUMAP(v102_sf7, dims = 1:30)

#cell clustering
# will perform hierarchical clustering
v102_sf7 <- FindNeighbors(v102_sf7, dims = 1:30)

#determing resolution for clustering
#determining clustering resolution
#running multiple resolutions to build clustertree
v102_sf7 <- FindClusters(v102_sf7, resolution = seq(0,2, 0.25))

#visualising clustree
clustree(v102_sf7)
#can run additional resoltuions to refine resolution choice
v102_sf7 <- FindClusters(v102_sf7, resolution = 1.25)
v102_sf7 <- FindClusters(v102_sf7, resolution = 1.75)

clustree(v102_sf7)
clustree(v102_sf7, node_colour = "sc3_stability")
clustree(v102_sf7, node_colour = "MPEG1", node_colour_aggr = "median")

#once chosen a resolution, rerun to store chosen resoltuion in seurat_clusters slot in the metadata
v102_sf7 <- FindClusters(v102_sf7, resolution = 1.25)

#visualise clusters 
DimPlot(v102_sf7, label = TRUE) + NoLegend()

#check how many cells per cluster
table(Idents(v102_sf7))
0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22   23   24   25 
1148  893  877  853  638  562  455  447  304  293  285  249  240  209  199  186  183  176  162   80   63   61   48   43   34   29 

#find markers for each cluster
v102_sf7.markers <- FindAllMarkers(v102_sf7, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
top5 <- v102_sf7.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
DoHeatmap(v102_sf7, features = top5$gene, size = 5,)

#create Excel file of markers for each cluster
write.csv(v102_sf7.markers, file = "v102_sf7_cluster_markers.csv", row.names = T, quote = F)

#rename clusters
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25)
new.cluster.ids <- c("Cardiomyocytes", "Cardiomyocytes","Cardiomyocytes", "Endothelium", "Cardiomyocytes", "Cardiomyocytes", "Cardiomyocytes", "Endothelium", "Cardiomyocytes", "Fibroblasts", "Endothelium", "Cardiomyocytes","Cardiomyocytes", "Smooth Muscle", "Macrophages", "T Cells", "Endo/CM", "Pacemaker Cardiomyocytes", "Erythrocytes", "Granulocytes", "HSCs", "Smooth Muscle", "Erythrocytes/CM", "Epicardium", "Thrombopoeitic Cells", "Neutrophils")
Idents(v102_sf7) <- plyr::mapvalues(Idents(v102_sf7), from = current.cluster.ids, to = new.cluster.ids)
DimPlot(v102_sf7, label = TRUE)

v102_sf7$celltype <- Idents(v102_sf7)

FeaturePlot(v102_sf7, features="arhgap27", cols = c("light grey", "red"))

VlnPlot(v102_sf7, features="nFeature_RNA")

#save seurat object
saveRDS(v102_sf7, file="v102_sf7_211221.rds")


######220118
v102 <- readRDS(file="E:/v102_sf7_211221.rds")

Idents(v102) <- "seurat_clusters"
bcells <- WhichCells(v102, expression = cd37>1)
v102 <- SetIdent(v102, cells=bcells, value="B Cells")

current.cluster.ids <- c("B Cells", 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25)
new.cluster.ids <- c("B Cells", "Cardiomyocytes", "Cardiomyocytes","Cardiomyocytes", "Endothelium", "Cardiomyocytes", "Cardiomyocytes", "Cardiomyocytes", "Endothelium", "Cardiomyocytes", "Fibroblasts", "Endothelium", "Cardiomyocytes","Cardiomyocytes", "Smooth Muscle", "Macrophages", "T Cells", "Endo/CM", "Pacemaker Cardiomyocytes", "Erythrocytes", "Endo/Leuko", "HSCs", "Smooth Muscle", "Erythrocytes/CM", "Epicardium", "Thrombopoeitic Cells", "Neutrophils")
Idents(v102) <- plyr::mapvalues(Idents(v102), from = current.cluster.ids, to = new.cluster.ids)
DimPlot(v102, label = TRUE)

v102$celltype <- Idents(v102)

VlnPlot(v102, features="nFeature_RNA")

markers <- c("myl7", "tnnc1a", "chrm2a", "kcnq1.1", "epas1b","cahz", "hbaa2","tcf21", "col1a1a", "thbs1b", "ITIH3", "wt1b", "h6pd", "ENSAMXG00005004684","c1qb", "c1qa","cd37", "runx3",
             "ctsl.1", "ncf1", "mmp9",  "vwa11")

Idents(v102) <- "celltype"

DotPlot(v102, features=markers)

v102$celltype <- factor(v102$celltype,levels = c("Cardiomyocytes", "Pacemaker Cardiomyocytes", "Endothelium", "Endo/CM", "Endo/Leuko", "Erythrocytes", "Erythrocytes/CM", "Fibroblasts", "Smooth Muscle", "Epicardium", "HSCs", "Macrophages", "B Cells", "T Cells", "Neutrophils", "Thrombopoeitic Cells"))
                       
p <- DotPlot(v102, features=markers)
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

p <- DotPlot(sf7, features=markers)
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

Idents(v102) <- "celltype"
plot1 <- FeatureScatter(v102, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 
