sf7 <- readRDS(file="G:/Analysis done on lab computer/single cell annotations paper/re analysis on realignment/sf7_211221.rds")
v102_sf7 <- readRDS(file="G:/Analysis done on lab computer/single cell annotations paper/re analysis on realignment/v102_sf7_211221.rds")

sf7$genome <- "v2"
v102_sf7$genome <- "v102"

#select integration features
sf7.features <- SelectIntegrationFeatures(object.list = list(sf7,v102_sf7), nfeatures = 3000)
options(future.globals.maxSize = 2000 * 1024^2)
sf7.list <- PrepSCTIntegration(object.list = list(sf7,v102_sf7), anchor.features = sf7.features, 
                               verbose = FALSE)

#identify anchors and integrate the datasets
sf7.anchors <- FindIntegrationAnchors(object.list = sf7.list, normalization.method = "SCT", 
                                      anchor.features = sf7.features, verbose = FALSE)
sf7.integrated <- IntegrateData(anchorset = sf7.anchors, normalization.method = "SCT", 
                                verbose = FALSE)

memory.limit(size = 30000)

#proceed downstream with clustering analysis
DefaultAssay(sf7.integrated) <- "integrated"

#when using sct integration pipeline DO NOT RUN scale data
sf7.integrated <- RunPCA(sf7.integrated, verbose = FALSE)
ElbowPlot(sf7.integrated, ndims=50)
sf7.integrated <- RunUMAP(sf7.integrated, dims = 1:30)
plots <- DimPlot(sf7.integrated, group.by = "genome")
plots

#use resoltion of 1.5, same as before
sf7.integrated <- sf7.integrated %>% FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = seq(0.5,2, 0.25))
library(clustree)
clustree(sf7.integrated)

#use resolution=1.25
sf7.integrated <- sf7.integrated %>% FindClusters(resolution = 1.25) %>% 
  identity()

DimPlot(sf7.integrated, label=T)

#check how many cells/cluster
table(Idents(sf7.integrated))

Idents(sf7.integrated) <- "seurat_clusters"
0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22   23   24   25   26   27   28   29 
1699 1661 1408 1178 1077 1047  950  857  763  710  630  581  546  519  504  478  430  424  370  351  341  227  215  179  135   93   92   87   57   34 

table(Idents(sf7.integrated), sf7.integrated$genome)
v102  v2
0   877 822
1   814 847
2   710 698
3   620 558
4   544 533
5   531 516
6   459 491
7   390 467
8   376 387
9   339 371
10  313 317
11  289 292
12  267 279
13  249 270
14  249 255
15  238 240
16  209 221
17  208 216
18  183 187
19  147 204
20  173 168
21  113 114
22   98 117
23   82  97
24   64  71
25   47  46
26   40  52
27   49  38
28   28  29
29   11  23


#Find markers per cluster per genome 
head(sf7.integrated@meta.data)

#create column identifying sample origin
sf7.integrated$genome_seurat <- paste(sf7.integrated$genome, sf7.integrated$seurat_clusters, sep = "_")

Idents(sf7.integrated) <- "genome_seurat"

#once clustering is complete, switch back to RNA slot for rest of downstream analysis
DefaultAssay(sf7.integrated) <- "RNA"

#not sure if should normalise the RNA slot 
#need to normalise RNA data slot to adjust for differences in sequencing depth between cells in the integrated Seuart object
sf7.integrated <- NormalizeData(sf7.integrated, verbose = FALSE)

#find markers for each cluster, need to scale the RNA slot fist before running analysis
all.genes <- rownames(x = sf7.integrated)
sf7.integrated <- ScaleData(sf7.integrated, features = all.genes)

Idents(sf7.integrated) <- "seurat_clusters"
integrated.markers <- FindAllMarkers(sf7.integrated, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
top5 <- integrated.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
DoHeatmap(sf7.integrated, features = top5$gene, size = 5,)

Idents(sf7.integrated) <- "seurat_clusters"
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29)
new.cluster.ids <- c("Cardiomyocytes","Cardiomyocytes", "Cardiomyocytes", "Cardiomyocytes", "Cardiomyocytes", "Endothelium", "Cardiomyocytes","Cardiomyocytes", "Endothelium", "Cardiomyocytes", "Endothelium","Fibroblasts", "Cardiomyocytes", "Cardiomyocytes", "Endothelium", "Pacemaker Cardiomyocytes", "T cells", "Smooth Muscle", "Endo_CM", "Erythrocytes", "Macrophages", "Endothelium", "Cardiomyocytes", "Endo_Leuko", "HSCs", "Epicardium", "B cells", "Erythrocytes_CM", "Neutrophils", "Platelets_Megakaryocytes")
Idents(sf7.integrated) <- plyr::mapvalues(Idents(sf7.integrated), from = current.cluster.ids, to = new.cluster.ids)
DimPlot(sf7.integrated, label = TRUE)
sf7.integrated$mcte <- Idents(sf7.integrated)

#find markers for each cluster
Idents(sf7.integrated) <- "genome_seurat"
genome.markers <- FindAllMarkers(sf7.integrated, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

#create excel file of top markers for each cluster
write.csv(genome.markers, file = "sf7_genome_clustermarkers.csv", row.names = T, quote = F)

#to plot split feature plots with grey colour
remotes::install_version("Seurat", version = "3.1.05")
FeaturePlot(sf7.integrated, features="acta2", split.by="genome")
FeaturePlot(sf7.integrated, features="tnnc1a", split.by="genome")

##looking at differences in endothelial cells
sf7.integrated$mct_genome <- paste(sf7.integrated$mct, sf7.integrated$genome, sep = "_")
Idents(sf7.integrated) <- "mct_genome"


DGE_endo <- FindMarkers(sf7.integrated, ident.1="Endothelium_v2", ident.2="Endothelium_v102", test.use="LR")


###LR test
### Calculate ranks for GSEA
genes <- DGE_endo
genes$X <- rownames(DGE_endo)

genes$ranks <- sign(genes$avg_log2FC) * -log10(genes$p_val)
head(genes)

#create ranks file containing only gene name and ranks
endo_ranked <- subset(genes, select = c("X", "ranks"))
colnames(endo_ranked) <- c("GeneName","rank")

#sort ranks in decreasing order and save your ranked gene list for GSEA use
endo_ranked <- endo_ranked[order(as.numeric(endo_ranked[,2]),decreasing = TRUE),]
write.csv(endo_ranked, "endo_ranked.csv", row.names = FALSE, quote = FALSE)

#convert fish genes to mouse homologs
AMgenes <- (endo_ranked[,1])

#convert AMgenes to mouse orthologues
library(biomaRt)
#set AM genome as database
ensembl = useMart("ensembl",dataset="amexicanus_gene_ensembl")

#have to run twice to run external gene names, eg mmp9 vs ensembl codes 
matrix <- getBM(attributes=c('external_gene_name', 'mmusculus_homolog_associated_gene_name'), 
                filters="external_gene_name",
                values = AMgenes, 
                mart = ensembl)

matrix2 <- getBM(attributes=c("ensembl_gene_id", 'mmusculus_homolog_associated_gene_name'), 
                 filters="ensembl_gene_id",
                 values = AMgenes, 
                 mart = ensembl)

musgenes <- c(matrix[,2], matrix2[,2])

#call mouse ensemble database
musensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")

mouse_entrez_id <- getBM(attributes=c(
  "external_gene_name", 
  "entrezgene_id"), 
  filters = "external_gene_name",
  values = musgenes,
  mart = musensembl,
  uniqueRows = FALSE)

#need to merge matrix into one file 
colnames(matrix) <- c("fish", "mmusculus_homolog_associated_gene_name")
colnames(matrix2) <- c("fish", "mmusculus_homolog_associated_gene_name")
fish2mousemap <- merge(matrix, matrix2, all=T)

fish2mousemap <- fish2mousemap %>% 
  mutate_all(~ifelse(. %in% c(""), NA, .)) %>% 
  na.omit()

fish2mousemap <- fish2mousemap[!duplicated(fish2mousemap$fish), ]

#find matching mouse genes
mouse.genes <- fish2mousemap[match(AMgenes, fish2mousemap[, "fish"]),"mmusculus_homolog_associated_gene_name"]

#change gene names in ranked list
endo_ranked$GeneName <- mouse.genes

endo_ranked <- endo_ranked %>% 
  mutate_all(~ifelse(. %in% c(""), NA, .)) %>% 
  na.omit()

write.csv(endo_ranked, file="endo_gseagenes.csv", quote=F)

require(org.Mm.eg.db)
endo_ranked$ID <- mapIds(org.Mm.eg.db, endo_ranked$GeneName, 'ENTREZID', 'SYMBOL')

endo_ranked[1:92,2] <- 9999
endo_ranked[1697:1712,2] <- -9999

#class(genes)
GenesRank_ID <- endo_ranked$rank
names(GenesRank_ID) <- endo_ranked$ID

pathways.hallmark <- readRDS(file = "E:/Single Cell all data/Analysis/GSEA/Mm.h.all.v7.1.entrez.rds") 
pathways.hallmark %>% 
  head() %>% 
  lapply(head)

# Run the GSEA
library(fgsea)
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=GenesRank_ID)

fgseaRes[, leadingEdge := mapIdsList(
  x=org.Mm.eg.db, 
  keys=leadingEdge,
  keytype="ENTREZID", 
  column="SYMBOL")]

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

ggplot(fgseaResTidy[fgseaResTidy$padj<0.25,], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=NES<0)) +
  coord_flip() +
  labs(x=NULL, y="Normalized Enrichment Score",
       title="endo_markers- Hallmarks") + 
  scale_fill_discrete(name="Enriched\npathways",
                      labels=c("v2", "v102")) +
  theme_classic()

fwrite(fgseaResTidy, file="endo_GSEA_fgseaRes_LR.txt", sep="\t", sep2=c("", " ", ""))

#####UV DNA response damage is significant using LR test
library(EnhancedVolcano)
EnhancedVolcano(DGE_endo,
                lab = row.names(DGE_endo),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                title = "DGE endo v2 vs v102, LR test",
                selectLab = "id11",
                xlim = c(-6, 6))

DGE_endo$X <- row.names(DGE_endo)
sig <- DGE_endo %>% filter(p_val_adj < 0.05) %>% pull(X)

v2_endo_cells <- WhichCells(sf7.integrated, ident="Endothelium_v2") 
v102_endo_cells <- WhichCells(sf7.integrated, ident="Endothelium_v102") 

v102_endo_specific_cells <- setdiff(v102__endo_cells, v2_cells_trimmed)

library(stringr)
v102_endo_cells <- str_sub(v102_endo_cells, 1, 18)
v2_endo_cells <- str_sub(v2_endo_cells, 1, 18)

v102_spef_endocells <- setdiff(v102_endo_cells, v2_endo_cells)
v2_spef_endocells <- setdiff(v2_endo_cells, v102_endo_cells)

####checking if any genes are DE in both datasets
DGE_endo_2 <- DGE_endo[!duplicated(DGE_endo$X), ]

Idents(sf7)
head(sf7@meta.data)
head(v102_sf7@meta.data)

tail(sf7.integrated@meta.data)

v2_cells <- colnames(sf7)
v102_cells <- colnames(v102_sf7)

v102specific_cells <- setdiff(v102_cells, v2_cells)

v2specific_cells <- setdiff(v2_cells, v102_cells)

v102specific_genes <- setdiff(v102_cells, v2_cells_trimmed)

v2specific_genes <- setdiff(v2_cells_trimmed, v102_cells)

v2_cells <- paste(v2specific_cells, "_1", sep="")
v102_cells <- paste(v102specific_cells, "_2", sep="")

all_missing_cells <- c(v102_cells, v2_cells)

DimPlot(sf7.integrated, cells.highlight = v2_cells)
DimPlot(sf7.integrated, cells.highlight = v102_cells)

DimPlot(sf7.integrated, cells.highlight = all_missing_cells)

DimPlot(sf7.integrated, cells.highlight = list(v102_cells, v2_cells), cols.highlight = c("#44C0C4","#F39F9B"), cols= "grey")

####how many assembly-specific genes are there
v2_genes <- row.names(sf7)
v102_genes <- row.names(v102_sf7)

v2_spef_genes <- setdiff(v2_genes, v102_genes)
v102_spef_genes <- setdiff(v102_genes, v2_genes)
#v102 specific #4,311, v2 specific 5,638

thromb_v1 <- WhichCells(sf7.integrated, expression=vwa11>1)
thromb_v2 <- WhichCells(sf7.integrated, expression=vwa11.1>1)

thromb_all <- c(thromb_v1, thromb_v2)

Idents(sf7.integrated) <- "seurat_clusters"
sf7.integrated <- SetIdent(sf7.integrated, cells=thromb_all, value="Thrombopoeitic Cells")
current.cluster.ids <- c("Thrombopoeitic Cells", 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29)
new.cluster.ids <- c("Thrombopoeitic Cells","Cardiomyocytes","Cardiomyocytes", "Cardiomyocytes", "Cardiomyocytes", "Cardiomyocytes", "Endothelium", "Cardiomyocytes","Cardiomyocytes", "Endothelium", "Cardiomyocytes", "Endothelium","Fibroblasts", "Cardiomyocytes", "Cardiomyocytes", "Endothelium", "Pacemaker Cardiomyocytes", "T cells", "Smooth Muscle", "Doublets", "Erythrocytes", "Macrophages", "Endothelium", "Cardiomyocytes", "Doublets", "HSCs", "Epicardium", "B cells", "Doublets", "Neutrophils", "Platelets_Megakaryocytes")
Idents(sf7.integrated) <- plyr::mapvalues(Idents(sf7.integrated), from = current.cluster.ids, to = new.cluster.ids)
DimPlot(sf7.integrated, label = TRUE)

sf7.integrated$celltype <- Idents(sf7.integrated)

sf7.integrated <- SetIdent(sf7.integrated, cells=thromb_all, value="30")

Idents(sf7.integrated) <- "seurat_clusters"
current.cluster.ids <- c("30", 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29)
current.cluster.ids <- c("30", 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29)
Idents(sf7.integrated) <- plyr::mapvalues(Idents(sf7.integrated), from = current.cluster.ids, to = new.cluster.ids)
DimPlot(sf7.integrated, label = TRUE)

CM <- FindConservedMarkers(sf7.integrated, ident.1="Cardiomyocytes", grouping.var = "genome", verbose = FALSE)
endo <- FindConservedMarkers(sf7.integrated, ident.1="Endothelium", grouping.var = "genome", verbose = FALSE)
Tcells <- FindConservedMarkers(sf7.integrated, ident.1="T cells", grouping.var = "genome", verbose = FALSE)
macro <- FindConservedMarkers(sf7.integrated, ident.1="Macrophages", grouping.var = "genome", verbose = FALSE)
fibro <- FindConservedMarkers(sf7.integrated, ident.1="Fibroblasts", grouping.var = "genome", verbose = FALSE)
sm <- FindConservedMarkers(sf7.integrated, ident.1="Smooth Muscle", grouping.var = "genome", verbose = FALSE)
pCM <- FindConservedMarkers(sf7.integrated, ident.1="Pacemaker Cardiomyocytes", grouping.var = "genome", verbose = FALSE)
epi <- FindConservedMarkers(sf7.integrated, ident.1="Epicardium", grouping.var = "genome", verbose = FALSE)
rbcs <- FindConservedMarkers(sf7.integrated, ident.1="Erythrocytes", grouping.var = "genome", verbose = FALSE)
Bcells <- FindConservedMarkers(sf7.integrated, ident.1="B cells", grouping.var = "genome", verbose = FALSE)
platelets <- FindConservedMarkers(sf7.integrated, ident.1="Platelets_Megakaryocytes", grouping.var = "genome", verbose = FALSE)
neutro <- FindConservedMarkers(sf7.integrated, ident.1="Neutrophils", grouping.var = "genome", verbose = FALSE)
hsc <- FindConservedMarkers(sf7.integrated, ident.1="HSCs", grouping.var = "genome", verbose = FALSE)
thromb <- FindConservedMarkers(sf7.integrated, ident.1="Thrombopoeitic Cells", grouping.var = "genome", verbose = FALSE)


write.csv(CM, file = "CM_conserved_genome_markers.csv", row.names = T, quote = F)
write.csv(endo, file = "endo_conserved_genome_markers.csv", row.names = T, quote = F)
write.csv(Tcells, file = "Tcells_conserved_genome_markers.csv", row.names = T, quote = F)
write.csv(macro, file = "macro_conserved_genome_markers.csv", row.names = T, quote = F)
write.csv(fibro, file = "fibro_conserved_genome_markers.csv", row.names = T, quote = F)
write.csv(endo/CM, file = "endo/CM_conserved_genome_markers.csv", row.names = T, quote = F)
write.csv(sm, file = "sm_conserved_genome_markers.csv", row.names = T, quote = F)
write.csv(pCM, file = "pCM_conserved_genome_markers.csv", row.names = T, quote = F)
write.csv(epi, file = "epicardium_conserved_genome_markers.csv", row.names = T, quote = F)
write.csv(rbcs, file = "rbcs_conserved_genome_markers.csv", row.names = T, quote = F)
write.csv(Bcells, file = "Bcells_conserved_genome_markers.csv", row.names = T, quote = F)
write.csv(platelets, file = "platelets_conserved_genome_markers.csv", row.names = T, quote = F)
write.csv(neutro, file = "neutro_conserved_genome_markers.csv", row.names = T, quote = F)
write.csv(hsc, file = "hsc_conserved_genome_markers.csv", row.names = T, quote = F)
write.csv(thromb, file = "thromb_conserved_genome_markers.csv", row.names = T, quote = F)

FeaturePlot(sf7.integrated, features="nFeature_RNA", split.by="genome")
Idents(sf7.integrated) <- "genome"
DimPlot(sf7.integrated, group.by="genome")

doublet_markers <- c("hbba2", "hbaa2", "hbae1.1", "hbba1", "cd74a", "lcp1", "cxcr4b", "f8", "lfng", "epas1b",  "tnni4b.1", "myl7", "cox4i1l")

doublets <- WhichCells(sf7.integrated, idents =c("0","5", "18","20", "27","23"))

sf7.integrated$seurat_clusters2 <- Idents(sf7.integrated)

sf7.integrated$seurat_clusters2 <- factor(sf7.integrated$seurat_clusters2, levels = c("0", "1", "2", "3", "4", "5", "20", "18", "23", "27", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "19", "21", "22", "24", "25", "26", "28", "29"))

Idents(sf7.integrated) <- "seurat_clusters2"
DoHeatmap(sf7.integrated, features = doublet_markers, size = 5,cells=doublets)


head(sf7@meta.data)

DimPlot(sf7.integrated, cells.highlight = list(v102cells_2, v2cells_1), cols.highlight = c("#44C0C4","#F39F9B"), cols= "grey")
