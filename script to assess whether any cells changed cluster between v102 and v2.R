###########exploring if any cells changed cluster between v102 and v2
#####tidyverse
library(tidyverse)

sf7_t <- sapply(sf7cellident, "length<-", max(lengths(sf7cellident))) 

v102_t <- sapply(v102cellident, "length<-", max(lengths(v102cellident))) 

sf7_t <- as.tibble(sf7_t)

names(sf7_t)[names(sf7_t) == "Platelets/Megakaryocytes"] <- "Platelets_Megakaryocytes"
names(sf7_t)[names(sf7_t) == "Endo/Leuko"] <- "Endo_Leuko"
names(sf7_t)[names(sf7_t) == "Thrombopoeitic Cells"] <- "Thrombopoeitic_Cells"


sf7_t <- sf7_t %>% 
  pivot_longer(Cardiomyocytes:Platelets_Megakaryocytes,
      names_to = "v2", 
    
    values_drop_na = TRUE
  )

#gives warning message but works
sf7_t <- remove_rownames(sf7_t)
sf7_t <- column_to_rownames(sf7_t, var = "value")

#####for v102
v102_t <- as.tibble(v102_t)

names(v102_t)[names(v102_t) == "Endo/Leuko"] <- "Endo_Leuko"
names(v102_t)[names(v102_t) == "Erythrocytes/CM"] <- "Erythrocytes_CM"
names(v102_t)[names(v102_t) == "Thrombopoeitic Cells"] <- "Thrombopoeitic_Cells"
names(v102_t)[names(v102_t) == "Endo/CM"] <- "Endo_CM"

v102_t <- v102_t %>% 
  pivot_longer(Cardiomyocytes:Thrombopoeitic_Cells,
               names_to = "v102", 
               values_drop_na = TRUE)
  
v102_t <- remove_rownames(v102_t)
v102_t <- column_to_rownames(v102_t, var = "value")

####merge columns v2 and v102 together, matching row names and filling in NA if necessary

b <- merge(v102_t, sf7_t)

b$Same <- b$v102 == b$v2

b$asint <- as.integer(b$Same)

result <- filter(b, asint == 0)


#######
head(sf7.integrated@meta.data)
Idents(sf7.integrated) <- "mct"
sf7integrated_cellident <- CellsByIdentities(sf7.integrated)

sf7integrated_cellident <- sapply(sf7integrated_cellident, "length<-", max(lengths(sf7integrated_cellident))) 

sf7integrated_cellident <- as.tibble(sf7integrated_cellident)

names(sf7integrated_cellident)[names(sf7integrated_cellident) == "T cells"] <- "T Cells"
names(sf7integrated_cellident)[names(sf7integrated_cellident) == "B cells"] <- "B Cells"

sf7integrated_cellident <- sf7integrated_cellident %>% 
  pivot_longer(Cardiomyocytes:Platelets_Megakaryocytes,
               names_to = "integrated", 
               
               values_drop_na = TRUE
  )

###remove v102 or v2 specifier in cell name ID
library(stringr)

sf7integrated_cellident$value <- str_sub(sf7integrated_cellident$value, 1, 18)

c <- merge(b, sf7integrated_cellident)

c$Same_v102_int <- c$v102 == c$integrated

c$Same_v2_int <- c$v2 == c$integrated

c$asint_Same_v102_int <- as.integer(c$Same_v102_int)

c$asint_Same_v2_int <- as.integer(c$Same_v2_int)

result2 <- filter(c, asint_Same_v102_int == 0)
result3 <- filter(c, asint_Same_v2_int == 0)


####assessing doublets
Idents(sf7.integrated) <- "seurat_clusters"
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29)
new.cluster.ids <- c("Cardiomyocytes","Cardiomyocytes", "Cardiomyocytes", "Cardiomyocytes", "Cardiomyocytes", "Endothelium", "Cardiomyocytes","Cardiomyocytes", "Endothelium", "Cardiomyocytes", "Endothelium","Fibroblasts", "Cardiomyocytes", "Cardiomyocytes", "Endothelium", "Pacemaker Cardiomyocytes", "T cells", "Smooth Muscle", "Doublets", "Erythrocytes", "Macrophages", "Endothelium", "Cardiomyocytes", "Doublets", "HSCs", "Epicardium", "B cells", "Doublets", "Neutrophils", "Platelets_Megakaryocytes")
Idents(sf7.integrated) <- plyr::mapvalues(Idents(sf7.integrated), from = current.cluster.ids, to = new.cluster.ids)
DimPlot(sf7.integrated, label = TRUE)
sf7.integrated$doublets <- Idents(sf7.integrated)

Idents(sf7.integrated) <- "doublets"
sf7integrated_doub_cellident <- CellsByIdentities(sf7.integrated)

sf7integrated_doub_cellident <- sapply(sf7integrated_doub_cellident, "length<-", max(lengths(sf7integrated_doub_cellident))) 

sf7integrated_doub_cellident <- as.tibble(sf7integrated_doub_cellident)

names(sf7integrated_doub_cellident)[names(sf7integrated_doub_cellident) == "T cells"] <- "T Cells"
names(sf7integrated_doub_cellident)[names(sf7integrated_doub_cellident) == "B cells"] <- "B Cells"

sf7integrated_doub_cellident <- sf7integrated_doub_cellident %>% 
  pivot_longer(Cardiomyocytes:Platelets_Megakaryocytes,
               names_to = "integrated", 
               
               values_drop_na = TRUE
  )

sf7integrated_doub_cellident$value <- str_sub(sf7integrated_doub_cellident$value, 1, 18)

v102cells <- result2$value
v2cells <- result3$value

result4 <- sf7integrated_doub_cellident[sf7integrated_doub_cellident$value %in% v102cells, ]
result5 <- sf7integrated_doub_cellident[sf7integrated_doub_cellident$value %in% v2cells, ]

result4$doub_check <- result4$integrated == "Doublets"

result4_1 <- filter(result4, doub_check == TRUE)

result5 <- sf7integrated_doub_cellident[sf7integrated_doub_cellident$value %in% v2cells, ]

result5$doub_check <- result5$integrated == "Doublets"

result5_1 <- filter(result5, doub_check == TRUE)

cells_1 <- result4_1$value

v102cells_2 <- paste(v102cells, "_2", sep="")
v2cells_1 <- paste(v2cells, "_1", sep="")

DimPlot(sf7.integrated, cells.highlight=v102cells_2)
DimPlot(sf7.integrated, cells.highlight=v2cells_1)

DimPlot(sf7.integrated, cells.highlight = list(v102cells_2, v2cells_1), cols.highlight = c("#44C0C4","#F39F9B"), cols= "grey")

cells <- paste(cells_1, "_2", sep="")
