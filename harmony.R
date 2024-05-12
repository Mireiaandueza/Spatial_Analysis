# HARMONY INTEGRATION FOR FULL DATASET

# LIBRARIES
library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(rhdf5)
library(devtools)
library(harmony)
library(cowplot)


# DIRECTORIES
samples_dir <- "./ST_files"
out_dir <-  "./HARMONY"
setwd(samples_dir)

# LOAD DATA
ST <- readRDS(file = "ST_merge.rds")

# TRANSFORM AND RUN HARMONY

# Create object with only count data
ST_harm <- CreateSeuratObject(ST@assays$Spatial$counts, assay = "Spatial")
ST_harm <- AddMetaData(
  object = ST_harm,
  metadata = ST@meta.data[c("orig.ident","nCount_Spatial" ,"nFeature_Spatial", 
                            "Genotype","Timepoint", "Genotype_Timepoint", 
                            "Sample", "Slide")],
  
)


#ST_harm <- NormalizeData(ST_harm, assay="Spatial")
#ST_harm <- FindVariableFeatures(ST_harm, assay="Spatial")
#ST_harm <- ScaleData(ST_harm, assay="Spatial")

# SCT Normalize
ST_harm <- SCTransform(ST_harm, assay = "Spatial")
DefaultAssay(ST_harm) <- "SCT"

# Run PCA
ST_harm <- RunPCA(ST_harm , verbose = TRUE)

# Run HARMONY
ST_harm <- RunHarmony(object = ST_harm, group.by.vars="Slide", verbose = T, assay.use="SCT")


# Run clustering and UMAP at resolution 0.2
ST_harm <- ST_harm %>%
    FindNeighbors(reduction = "harmony") %>%
    FindClusters(resolution = 0.2) 

ST_harm <- ST_harm %>%
    RunUMAP(reduction = "harmony",  dims = 1:20)


# SAVE
setwd(out_dir)
saveRDS(ST_harm, "ST_harmony.rds")


# PLOT RESULTS
# Clusters
pdf("UMAP_Clustering.pdf")
p2 <- DimPlot(ST_harm, reduction = "umap", label = TRUE,  pt.size = .1)
p2
dev.off()
# Clusters colored by Genotype_Timepoint
pdf("UMAP_Clustering_GT.pdf")
p1 <- DimPlot(ST_harm, reduction = "umap", group.by = "Genotype_Timepoint", pt.size = .1)
p1
dev.off()
# Clusters colored by Genotype
pdf("UMAP_Clustering_G.pdf")
p1 <- DimPlot(ST_harm, reduction = "umap", group.by = "Genotype", pt.size = .1)
p1
dev.off()
# Clusters colored by Slide
pdf("UMAP_Clustering_S.pdf")
p1 <- DimPlot(ST_harm, reduction = "umap", group.by = "Slide", pt.size = .1)
p1
dev.off()
# Clusters colored by Sample
pdf("UMAP_Clustering_Sample.pdf")
p1 <- DimPlot(ST_harm, reduction = "umap", group.by = "Sample", pt.size = .1)
p1
dev.off()


# Evaluate integration
Per226 <- subset(ST_harm, subset = Sample == "PER226")
Per189 <- subset(ST_harm, subset = Sample == "PER189")
p1 <- DimPlot(Per226, reduction = "umap", label = TRUE, group.by = "Slide")
pdf("Per226_by_Slide.pdf")
p1 
dev.off()
p1 <- DimPlot(Per189, reduction = "umap", label = TRUE, group.by = "Slide")
pdf("Per189_by_Slide.pdf")
p1 
dev.off()


# MARKERS for clusters
pbmc.markers <- FindAllMarkers(ST_harm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC) %>% print(n = 30)

# Save table
write.table(pbmc.markers, file ="markers0_2.csv" , sep = ",")

