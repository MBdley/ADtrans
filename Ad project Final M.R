## libraries #######################
library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)
library(grid)
library(viridis)


# AD
## DATA PREP #######################
## Load the AD dataset
setwd('C:\\Users\\mhkbr\\Desktop\\adsnMtx\\mtx')
ad.lists <- c('AD1_AD2','AD3_AD4','AD5_AD6')
ad.data <- Read10X(data.dir = ad.lists)
ad <- CreateSeuratObject(counts = ad.data, project = "ad_scRNAseq", min.cells = 3, min.features = 200)
ad$stim <- "AD"

control.lists <- c('Ct1_Ct2','Ct3_Ct4','Ct5_Ct6')
control.data <- Read10X(data.dir = control.lists)
control <- CreateSeuratObject(counts = control.data, project = "control_scRNAseq", min.cells = 3, min.features = 200)
control$stim <- "CTRL"

ad.combine <- merge(ad, y = control, add.cell.ids = c("AD", "C"), project = "combine_scRNAseq")
ad.combine
head(colnames(ad.combine))
table(ad.combine$orig.ident)

## QC #######################################################
ad.combine [["percent.mt"]] <- PercentageFeatureSet(ad.combine, pattern = "^MT-")

# Visualize QC metrics as a violin plot
#VlnPlot(ad.combine, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(ad.combine, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ad.combine, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

##ad.combine <- subset(ad.combine, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 5) # see @meta.data

ad.combine
nFet <- data.matrix(ad.combine[['nFeature_RNA']], rownames.force = NA)
pctl <- quantile(nFet, c(.05, .95)) 
ad.combine <- subset(ad.combine, subset = nFeature_RNA > unname(pctl)[1] & nFeature_RNA < unname(pctl)[2] & percent.mt < 10)

## normalization ################################################
ad.combine <- NormalizeData(ad.combine, normalization.method = "LogNormalize", scale.factor = 10000)

## identification of highly variable features #####################
###ad.combine <- FindVariableFeatures(ad.combine, selection.method = "vst", nfeatures = 2000)

ad.combine <- FindVariableFeatures(object = ad.combine, selection.method = "vst", mean.cutoff = 
     c(0.0125, 3), dispersion.cutoff = c(0.5, Inf), nfeatures = 2000)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

ad.combine<- FindVariableFeatures(ad.combine, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(ad.combine), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(ad.combine)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
## scaling ######################################################
all.genes <- rownames(ad.combine)
ad.combine <- ScaleData(ad.combine, features = all.genes)
## PCA #########################################################
ad.combine <- RunPCA(ad.combine, features = VariableFeatures(object = ad.combine))

## Clustering #######################################################
ElbowPlot(ad.combine)
ad.combine <- FindNeighbors(ad.combine, dims = 1:10)
ad.combine <- FindClusters(ad.combine, resolution = 0.1)

## UMAP/tSNE #########################################################
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
ad.combine <- RunUMAP(ad.combine, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(ad.combine, reduction = "umap", label = FALSE)

 ## Differential gene expression analysis ############################

# find markers for every cluster compared to all remaining cells, report only the positive ones
ad.markers <- FindMarkers(ad.combine, ident.1 = "AD", ident.2 = "CTRL", min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox", only.pos = TRUE)
ad.markers %>% top_n(n = 100, wt = avg_logFC)

ad.markers <- FindAllMarkers(ad.combine)
# OR individual celltype Wilcox test
Oligo.markers <- FindMarkers(ad.combine, ident.1 = "Oligo_AD", indent.2 = "Oligo_CTRL", logfc.threshold = 0.25, test.use = "wilcox", only.pos = FALSE)
Endo.markers <- FindMarkers(ad.combine, ident.1 = "Endo_AD", ident.2 = "Endo_CTRL", logfc.threshold = 0.25, test.use = "wilcox", only.pos = FALSE)
Micro.markers <- FindMarkers(ad.combine, ident.1 = "Micro_AD", ident.2 = "Micro_CTRL", logfc.threshold = 0.25, test.use = "wilcox", only.pos = FALSE)
OPC.markers <- FindMarkers(ad.combine, ident.1 = "OPC_AD", ident.2 = "OPC_CTRL", logfc.threshold = 0.25, test.use = "wilcox", only.pos = FALSE)
Astro.markers <- FindMarkers(ad.combine, ident.1 = "Astro_AD", ident.2 = "Astro_CTRL", logfc.threshold = 0.25, test.use = "wilcox", only.pos = FALSE)
Neuron.markers <- FindMarkers(ad.combine, ident.1 = "Neurons_AD", ident.2 = "Neurons_CTRL", logfc.threshold = 0.25, test.use = "wilcox", only.pos = FALSE)


ad.combine$celltype.stim <- paste(Idents(ad.combine), ad$stim, sep = "_")
ad.combine$celltype <- Idents(ad.combine)
Idents(ad.combine) <- "celltype.stim"
Idents(ad.combine) <- "seurat_clusters"
Idents(ad.combine) <- "stim"

ad.new <-subset(ad.combine, idents = c("Microglia_AD", "Microglia_CTRL", "Endothelial_AD", "Endothelial_CTRL", "Neurons_AD", "Neurons_CTRL","Astrocytes_AD", "Astrocytes_CTRL", "OPC_AD", "OPC_CTRL","Oligodendrocytes_AD", "Oligodendrocytes_CTRL"))

ad.combine <- subset(ad.combine, features = rownames(counts))
                          
#Plots of genes Identified as DEGs
UP = c ( "MT-ND4", "MT-CO2", "MT-ATP6", "LINGO1", "HSPA1A", "CRYAB")
EndothelialUP = c("SPP1", "CD163","NEAT1")
Downregulated = c ()

VlnPlot(ad.new, features = c(Allup2))
VlnPlot(ad.new, features = c(MicrogliaUP))
VlnPlot(ad.new, features = c("BIN1", "CTNNA3","GFAP"))
VlnPlot(ad.new, features = c(Downregulated))
                                
## Makers for gene identification
micro = c("CD74")
Astro = c("AQP4")
Neurons = c("NDUFA4")
oligo = c("MBP")
endo = c("CEMIP")
OPC = c("PCDH15")


AD.cc <- subset(ad.combine, idents = c("Oligo", "Astro", "Oligo", "OPC",
                                        "Neurons", "Micro", "Endo"))

VlnPlot(AD.cc, features = UP, pt.size = 0, combine = TRUE, ncol = 3)  + theme(legend.position = 'none')
VlnPlot(AD.cc, features = c(EndothelialUP), pt.size = 0, combine = TRUE, ncol = 3)  + theme(legend.position = 'none')
VlnPlot(AD.cc, features = c(Downregulated), pt.size = 0, combine = TRUE, ncol = 3)  + theme(legend.position = 'none')

levels(ad.new) <- c("Microglia_AD", "Microglia_CTRL", "Endothelial_AD", "Endothelial_CTRL", "Neurons_AD", "Neurons_CTRL","Astrocytes_AD", "Astrocytes_CTRL", "OPC_AD", "OPC_CTRL","Oligodendrocytes_AD", "Oligodendrocytes_CTRL")
levels(ad.combine) <- c(0, 2, 1, 3, 4 , 5, 6, 7)

#microglia 
VlnPlot(ad.combine, features = micro)
#astrocytes
VlnPlot(ad.combine, features = Astro)
#neurons
VlnPlot(ad.combine, features = Neurons)
#oligo
VlnPlot(ad.small, features = oligo)
#OPC
VlnPlot(ad.combine, features = OPC)
#Endo
VlnPlot(ad.combine, features = endo)

Markers <- VlnPlot(ad.combine, features = c(micro, Astro,Neurons, oligo, OPC, endo), 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = Markers, ncol = 3, label = FALSE)+ theme(legend.position="none")



# on UMAP plot
FeaturePlot(ad.combine, features = c("AQP4"))
FeaturePlot(ad.combine, features = c("CD74"))
FeaturePlot(ad.combine, features = c("MBP"))
FeaturePlot(ad.combine, features = c("PCDH15"))
FeaturePlot(ad.combine, features = c("FLT1"))
FeaturePlot(ad.combine, features = c("RBFOX1"))

new.cluster.ids <- c("Oligodendrocyte", "Astrocyte", "Oligodendrocyte", "OPC",
                     "Neurons", "Micro", "Unknown", "Endothelium")

new.cluster.ids <- c(1,2,3,4,5,6,7,8)

names(new.cluster.ids) <- levels(ad.combine)
ad.combine <- RenameIdents(ad.combine, new.cluster.ids)




# AD vs Control #####
DimPlot(ad.combine, reduction = 'umap', group.by = 'stim')
plot_grid(P1, P2)
# Subset microglia #######
microglia.c0 <- subset(ad.combine, idents = "Micro")
microglia.c0 <- NormalizeData(microglia.c0, normalization.method = "LogNormalize", scale.factor = 10000)
microglia.c0 <- FindVariableFeatures(microglia.c0, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(microglia.c0)
microglia.c0 <- ScaleData(microglia.c0, features = all.genes)
microglia.c0 <- RunPCA(microglia.c0, features = VariableFeatures(object = microglia.c0))

microglia.c0 <- FindNeighbors(microglia.c0, dims = 1:10)
microglia.c0 <- FindClusters(microglia.c0, resolution = 0.2)

micro.markers <- FindMarkers(microglia.c0, ident.1= "AD", ident.2 = "CTRL", test.use = "wilcox",  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
micro.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
top10 <- micro.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
Idents(microglia.c0) <- "stim"


Idents(object = microglia.c0) <- "stim"

geneOfInterest = c ( "MT-ND4", "MT-CO2", "MT-ATP6", "SPP1", "CD163","NEAT1", "MT-CYB", "MT-ND2", "BIN1")

VlnPlot(microglia.c0, features = geneOfInterest)

########## Sub Cluster UMAP########
DimPlot(ad.combine, reduction = "umap", label = FALSE, pt.size = 0.5, cols = c("Other Cell Types"='#e8eaee', "m1"='#0CB702', "m2"='#FF61CC', "m3"='#8494FF', "m4"='#ffff66', "m5"='#25aff5', "m6"='#00ffff'))


allsub.cluster.ids <- c("Other Cell Types", "Other Cell Types","Other Cell Types","Other Cell Types","Other Cell Types", "m2", "Other Cell Types","m5", "m3", "m1", "m4")
names(allsub.cluster.ids) <- levels(ad.combine)
ad.combine <- RenameIdents(ad.combine, allsub.cluster.ids)

Idents(object = ad.combine) <- "seurat_clusters"
Idents(object = microglia.c0) <- "stim"
levels(ad.combine)
