
#CODE FOR MERGING DATASET AND DOWNSTREAM ANALYSIS: normalization, 
  #dimension reduction and clustering

# LIBRARIES
library(Matrix)
library(Seurat)
library(hdf5r)
library(ggplot2)
library(patchwork)
library(dplyr)
library(rhdf5)
library(SeuratData)
library(devtools)
library(clustree)
library(anndata)
library(stringr)
library(purrr)

# DIRECTORIES
samples_dir <-"./ST_files"
out_path <- "./MERGE"
setwd(samples_dir)

dir.create(file.path(getwd(), "MERGE"))

# FUNCTIONS
source("norm_and_clustering.R")

# LOAD SAMPLES
files <- c("ST1.rds","ST2.rds","ST3.rds","ST4.rds", "ST5.rds", "ST6.rds",
           "ST7.rds", "ST8.rds","ST09.rds","ST10.rds","ST11.rds","ST12.rds")
rds_data <- map(files, readRDS)


# Add slide metadata
rds_data <- map2(rds_data, files, function(obj, filename) {
  # Extract the numeric part from the file name
  slide_number <- as.numeric(str_extract(filename, "\\d+"))
  	
  # Transform the slide name in the object
  obj$Slide <- paste("ST", slide_number, sep = "")
 
  # Return the modified object
  return(obj)
})


# Add batch metadata
batches <- c("1","1","1","1","1","1","1","1","2","2","2","2")
rds_data <- map2(rds_data, batches, function(file, bt){file$batch <-bt 
return(file)})


# MERGE
ST_merged  <- merge(x = rds_data[[1]], y = c(rds_data[[2]], rds_data[[3]], 
                                             rds_data[[4]], rds_data[[5]], 
                                             rds_data[[6]], rds_data[[7]], 
                                             rds_data[[8]],rds_data[[9]],
                                             rds_data[[10]],rds_data[[11]]), 
                    add.cell.ids = c("ST1","ST2","ST3","ST4","ST5","ST6","ST7",
                                     "ST8","ST9","ST10","ST11"), 
                    merge.data = TRUE)

# Image names
names(ST_merged@images) <- c("ST1", "ST2", "ST3", "ST4", "ST5", "ST6", "ST7", 
                             "ST8", "ST9", "ST10", "ST11")
ST_merged$Slide <- as.factor(ST_merged$Slide)
ST_merged <- JoinLayers(ST_merged,  assay = "Spatial")


# FILTER OUT SAMPLES FROM OTHER STUDIES
ST_merged <- subset(x = ST_merged, subset = Sample == "PEO688", invert = TRUE)


# FILTER OUT REPLICATES (ST11)
ST_filtered <- subset(x = ST_merged, subset = Slide == "ST11", invert = TRUE)
ST_filtered@images <-ST_merged@images[-which(names(ST_merged@images) == "ST11")]


# FILTER OUT GANGLION 
ST_merged <- subset(x = ST_filtered, subset = Genotype == "Ganglion", 
                    invert = TRUE)

##########################################

# SOME QUALITY MARKERS
setwd(out_path)

# QUALITY nCount vs nFeature
pdf("quality.pdf")
plot2 <- FeatureScatter(ST_merged, feature1 = "nCount_Spatial", 
                        feature = "nFeature_Spatial", group.by = "Genotype")
plot2
dev.off()

##########################################

# CALL DOWNSTREAM ANALYSIS 
ST_SCT <- norm_and_clustering(ST_merged)

#Save dataset
saveRDS(ST_SCT, file = "ST_merge.rds")

##########################################

# DIMPLOT OF FEATURES
p1 <- DimPlot(ST_SCT, reduction = "umap", label = TRUE, group.by = "Genotype", 
              cols = c("darkorchid", "mediumorchid1","grey")) + 
  ggtitle("Clusters colored by genotype")
pdf("clustering_by_Genotype.pdf")
p1 
dev.off()

p1 <- DimPlot(ST_SCT, reduction = "umap", label = TRUE, group.by = "Timepoint",
              cols = c("darkorchid", "mediumorchid1","grey"))+
  ggtitle("Clusters colored by timepoint")
pdf("clustering_by_Timepoint.pdf")
p1 
dev.off()

p1 <- DimPlot(ST_SCT, reduction = "umap", label = TRUE, group.by = "Sample")+ 
  ggtitle("Clusters colored by sample")
pdf("clustering_by_Sample.pdf")
p1 
dev.off()

p1 <- DimPlot(ST_SCT, reduction = "umap", label = TRUE, group.by = "Slide")+ 
  ggtitle("Clusters colored by capture area")
pdf("clustering_by_Slide.pdf")
p1 
dev.off()

##########################################

# PLOT SOME MARKERS

#ACINAR 
pdf("Acinar_umap.pdf")
FeaturePlot(ST_SCT, features =c("Agr2", "Cel","Pnlip", "Ctrb1"), 
            reduction = "umap")
dev.off()

#ENDOCRINE
pdf("Endocrine.pdf")
FeaturePlot(ST_SCT, features = c("Sst", "Gcg","Ins1", "Ins2"), 
            reduction = "umap")
dev.off()

#EPITELIAL
pdf("Epitelial.pdf")
FeaturePlot(ST_SCT, features = c("Pecam1"), reduction = "umap")
dev.off()

#DUCTAL
pdf("Ductal_umap.pdf")
FeaturePlot(ST_SCT, features =c("Krt19", "Sox9"), reduction = "umap")
dev.off()


##########################################

# FIND CLUSTER MARKERS
Idents(ST_SCT) <- ST_SCT$SCT_snn_res.0.2
ST_SCT$seurat_clusters <- Idents(ST_SCT)

pbmc.markers <- FindAllMarkers(ST_SCT, only.pos = TRUE, min.pct = 0.25, 
                               logfc.threshold = 0.25)
pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC) %>% print(n = 30)
write.table(pbmc.markers, file ="markers0_2.csv" , sep = ",")

names(cluster_names) <- levels(ST_SCT)
ST_SCT <- RenameIdents(ST_SCT, cluster_names)
ST_SCT$seurat_clusters <- Idents(ST_SCT)
ST$Genotype <- factor(ST$Genotype, levels = c("WT","HET","DD"))

ST_SCT$seurat_clusters <- factor(ST_SCT$seurat_clusters, 
                                 levels = c("Acinar WT (I)","Acinar WT (II)",
                                            "Acinar HET (I)","Acinar HET (II)",
                                            "KRAS WT","KRAS DD", "Immune",
                                            "Endocrine","Per 106", 
                                            "Per 152 (Islet Mice)", "Per 222", 
                                            "Per 135", "Per 189",  "Cluster 8", 
                                            "Cluster 14"))

# STACKED BARPLOT OF CLUSTERS
# Convert cluster information to a data frame
cluster_data <- as.data.frame(table(ST_SCT$Genotype,ST_SCT$Timepoint, 
                                    ST_SCT$seurat_clusters))
colnames(cluster_data) <- c("Genotype","Timepoint", "Cluster", "Count")


# Create barplot
pdf("Stacked_Cluster.pdf")
ggplot(cluster_data, aes(x = Cluster, y = Count, fill = as.factor(Genotype))) +
  geom_bar(stat = "identity") +
  labs(title = "Cluster Distribution by Genotype") + 
  theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off()
