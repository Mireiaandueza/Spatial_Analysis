# INTEGRATION OF SINGLE CELL DATA FROM MANUEL SERRANO AND SCHLEISINGER
# USING SEURAT FIND-ANCHORS

#LIBRARIES
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


#PATHS
spatial_data <- "./VISION/" #DIRECTORY OF SPATIAL DATA
single_cell_data <- "./single_cell_data/"


#LOAD SPATIAL DATA (WITH SIGNATURES)
setwd(spatial_data)
ST <- readRDS("ST_signatures.rds")

#LOAD SINGLE CELL DATA
setwd(single_cell_data)

#Serrano
load("WT_CER_i4F_epithelial.rda")
load("i4F.rda")
load("CER.rda")
load("WT.rda")

#Schleisinger
load("sc_OP.rda")

# UPDATE OBJECTS
WT_CER_i4F_epithelial <- UpdateSeuratObject(WT_CER_i4F_epithelial)
i4F <- UpdateSeuratObject(i4F)
CER <- UpdateSeuratObject(CER)
WT <- UpdateSeuratObject(WT)
sc_OP <- UpdateSeuratObject(sc_OP)

# RENAME epithelial clusters (original dataset not annotated)
cluster_names <- c("Acinar_CER", "Acinar_Low", "Ductal_WT_LoW", 
"Acinar_WT","Ductal_CER_High","Acinar_High", "Epitelial_Reprogrammed", 
"Endocrine", "Macrophage_like", 
"Mesenchymal_like", "10")

names(cluster_names) <- levels(WT_CER_i4F_epithelial)
WT_CER_i4F_epithelial <- RenameIdents(WT_CER_i4F_epithelial, cluster_names)

# Create a list
datasets <- c("WT_CER_i4F_epithelial", "i4F","CER","WT","sc_OP")
single_cells <- list(WT_CER_i4F_epithelial = WT_CER_i4F_epithelial, 
                     i4F = i4F, CER = CER, WT = WT, sc_OP = sc_OP)

# Normalize query
ST <- NormalizeData(ST, assay = "Spatial")
ST <- FindVariableFeatures(ST, assay = "Spatial")
ST <- ScaleData(ST, assay = "Spatial")

# Define a function for data transfer
transfer_data <- function(ST, ref_data, dataset_name) {

  ref_data <- NormalizeData(ref_data, assay = "RNA")
  ref_data <- FindVariableFeatures(ref_data, assay = "RNA")
  ref_data <- ScaleData(ref_data, assay = "RNA")

  print("Normalized and Scaled")
  #Find anchors
  anchors <- FindTransferAnchors(reference = ref_data, query = ST, 
                                 normalization.method = "LogNormalize", 
                                 reference.assay = "RNA", 
                                 query.assay = "Spatial")
  # Transfer data
  predictions.assay <- TransferData(anchorset = anchors, 
                                    refdata = Idents(ref_data), 
                                    prediction.assay = FALSE, 
                                    weight.reduction = ST[["pca"]], 
                                    dims = 1:30)
  # Add as metaData
  ST <- AddMetaData(object = ST, metadata = predictions.assay$predicted.id, col.name =dataset_name)

  # Return
  return(ST)
}


# Create a new file named SINGLECELL_INTEGRATION
mainDir <- getwd()
dir.create(file.path(mainDir, "SINGLECELL_INTEGRATION"))
setwd(file.path(mainDir, "SINGLECELL_INTEGRATION"))

# Apply function for all datasets
for (dataset_name in datasets) {
  # Transfer data
  ST <- transfer_data(ST, single_cells[[dataset_name]], dataset_name)

  # Plot integration
  pdf(paste0(dataset_name, "_integrated_clusters.pdf"))
  Plot <- DimPlot(ST, group.by = dataset_name)
  print(Plot)
  dev.off()
}

# SAVE INTEGRATED DATASET
saveRDS(ST, file = "ST_singlecell.rds")

#PLOT ORIGINAL DATASETS
lapply(single_cells, function(single_cell){
	pdf(paste0(dataset_name, ".pdf"))
	Plot <- DimPlot(single_cell)
	print(Plot)
	dev.off()
})

