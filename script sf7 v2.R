#SF7 v2.0 genome
setwd("G:/Analysis done on lab computer/single cell annotations paper/re analysis on realignment")

#load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(clustree)

#load in sf7 data
#read sf7
sf7 <- Read10X(data.dir = "E:/Single Cell all data/realignment surface v103 210304/surface_E103e_SF7/outs/filtered_feature_bc_matrix")

#change ensembl codes in the matrix #change gene names in the expression matrix

dim(sf7)
27420  9018

#create sf7 based on no.of genes and cells
sf7 <- CreateSeuratObject(counts = sf7, min.cells = 2, min.features = 50, project="SFuninjured_v2")

#will create a scatter plot, good to visualise outlier cells to exclude
plot1 <- FeatureScatter(sf7, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 

#can also visualise as histograms
nUMI <- as.vector(sf7$nCount_RNA)
hist(x=nUMI,
     breaks = 100, main = "Number of unique molcular identifiers detected per cell",
     xlab = "nCount_RNA") + abline(v = 500, col = "red")

nUMI <- FetchData(sf7, vars = "nCount_RNA")
ggplot(nUMI, aes(x=nCount_RNA)) + geom_histogram(binwidth=100) + geom_vline(xintercept=15000)

nGenes <- as.vector(sf7$nFeature_RNA)
hist(x=nGenes,
     breaks = 100, main = "Number of genes detected per cell",
     xlab = "nFeature_RNA") + abline(v = 200, col = "red")

nGenes <- FetchData(sf7, vars = "nFeature_RNA")
ggplot(nGenes, aes(x=nFeature_RNA)) + geom_histogram(binwidth=20) + geom_vline(xintercept=130)

#draw scatter plot  intercepts to look at filtering
plot1 + geom_hline(yintercept=200) + geom_vline(xintercept = 28000) + geom_hline(yintercept=3500)

#check the Violin plots
VlnPlot(sf7, features = c("nFeature_RNA", "nCount_RNA"))

#filter based on no. of features/genes and counts/UMIs
#use violin plots with histograms to set thresholds
sf7 <- subset(x = sf7, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & nCount_RNA <28000)

#check how many cells left after filtering
table(Idents(sf7))
SFuninjured_v2 
8926 

#many genes removed during filtering
dim(sf7)
[1] 17528  8926

#normalise data, identify highly variable features and scale data
sf7 <- SCTransform(sf7, verbose = T)      

#run PCA
sf7 <- RunPCA(sf7, features=VariableFeatures(sf7))

#visualise PCs
DimHeatmap(sf7, dims = 1:15, cells = 500, balanced = T)
DimHeatmap(sf7, dims = 15:30, cells = 500, balanced = T)
DimHeatmap(sf7, dims = 30:45, cells = 500, balanced = T)
ElbowPlot(sf7, ndims=50)

#run UMAP
#chose dimensions of 30
sf7 <- RunUMAP(sf7, dims = 1:30)

#cell clustering
# will perform hierarchical clustering
sf7 <- FindNeighbors(sf7, dims = 1:30)


#determing resolution for clustering
#determining clustering resolution
#running multiple resolutions to build clustertree
sf7 <- FindClusters(sf7, resolution = seq(0,2, 0.25))

#visualising clustree
clustree(sf7)
#can run additional resoltuions to refine resolution choice
sf7 <- FindClusters(sf7, resolution = 1.25)
sf7 <- FindClusters(sf7, resolution = 1.75)

clustree(sf7)
clustree(sf7, node_colour = "sc3_stability")
clustree(sf7, node_colour = "MPEG1", node_colour_aggr = "median")

#once chosen a resolution, rerun to store chosen resoltuion in seurat_clusters slot in the metadata
sf7 <- FindClusters(sf7, resolution = 1.25)

#visualise clusters 
DimPlot(sf7, label = TRUE) + NoLegend()

#check how many cells per cluster
table(Idents(sf7))
0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22   23   24   25   26 
1249  906  759  732  709  491  450  443  443  347  295  264  252  239  217  215  197  185  184   73   49   49   43   43   39   30   23 

#find markers for each cluster
sf7.markers <- FindAllMarkers(sf7, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
top5 <- sf7.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
DoHeatmap(sf7, features = top5$gene, size = 5,)

#create Excel file of markers for each cluster
write.csv(sf7.markers, file = "sf7_v2_cluster_markers.csv", row.names = T, quote = F)

#rename clusters
#rename clusters
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26)
new.cluster.ids <- c("Cardiomyocytes","Cardiomyocytes","Cardiomyocytes","Cardiomyocytes", "Endothelium", "Cardiomyocytes", "Cardiomyocytes", "Endothelium", "Endothelium", "Cardiomyocytes", "Fibroblasts", "Cardiomyocytes", "Cardiomyocytes","Pacemaker Cardiomyocytes", "Erythrocytes", "Smooth Muscle", "Macrophages", "T Cells", "Endothelium", "HSCs", "B Cells", "Endo/Leuko", "Epicardium", "Erythrocytes", "Thrombopoeitic Cells", "Neutrophils", "Platelets/Megakaryocytes")
Idents(sf7) <- plyr::mapvalues(Idents(sf7), from = current.cluster.ids, to = new.cluster.ids)
DimPlot(sf7, label = TRUE)

sf7$celltype <- Idents(sf7)

FeaturePlot(sf7, features="arhgap27", cols = c("light grey", "red"))

VlnPlot(sf7, features="nFeature_RNA")

saveRDS(sf7, file="sf7_211221.rds")

markers <- c("cmlc1", "tnnt2b", "chrm2a", "kcnq1.1", "epas1b","cahz", "hbaa2","TCF21", "col1a1a", "thbs1b", "acta2", "tbx18", "wt1b", "h6pd", "c1qb", "c1qa","cd37", "runx3", "ctsl.1.11",
             "ncf1", "mmp9.1",  "ENSAMXG00000001109", "mpl")


sf7$celltype <- factor(sf7$celltype,levels = c("Cardiomyocytes", "Pacemaker Cardiomyocytes", "Endothelium", "Endo/Leuko",  "Erythrocytes", "Fibroblasts", "Smooth Muscle", "Epicardium", "HSCs", "Macrophages", "B Cells", "T Cells", "Granulocytes", "Neutrophils", "Thrombopoeitic Cells", "Platelets/Megakaryocytes"))

Idents(sf7) <- "celltype"

DotPlot(sf7, features=markers)

p <- DotPlot(sf7, features=markers)
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

Idents(sf7) <- "celltype"
plot1 <- FeatureScatter(sf7, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 
