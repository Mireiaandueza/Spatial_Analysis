# INTEGRATION OF SINGLE CELL DATA FROM MANUEL SERRANO AND SCHLEISINGER
# USING SPACEXR TOOL RCTD

#LIBRARIES
library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(rhdf5)
library(devtools)
library(spacexr)
library(Matrix)


# DIRECTORIES (single cell rda file to integrate)
spatial_data <- "./VISION/" #DIRECTORY OF SPATIAL DATA (with signature scores)
single_cell_data <- "./single_cell_data/WT_CER_i4F_epithelial.rda"
# i4F.rda,CER.rda,WT.rda

# Create a new file named RCTD_INTEGRATION
mainDir <- getwd()
dir.create(file.path(mainDir, "RCTD_INTEGRATION"))
out_dir <- file.path(mainDir, "RCTD_INTEGRATION")

file_name <- tools::file_path_sans_ext(basename(single_cell_data))

# LOAD SPATIAL DATA (WITH SIGNATURES)
setwd(spatial_data)
ST <- readRDS("ST_signatures.rds") 


#LOAD SINGLE CELL DATA
X <- load(single_cell_data)
single_cell <- get(X)
rm(X)

# UPDATE OBJECTS
single_cell <- UpdateSeuratObject(single_cell)

# CREATE RCTD OBJECTS

# REFERENCE DATA
counts <- single_cell@assays$SCT@counts	
rownames(counts) <- rownames(single_cell)
colnames(counts) <- colnames(single_cell)

cell_types <- Idents(single_cell)
names(cell_types) <- colnames(single_cell)
cell_types <- as.factor(cell_types) 
nUMI <- single_cell$nCount_SCT; names(nUMI) <- colnames(single_cell) 

reference <- Reference(counts, cell_types, nUMI)
print("referenced")

# SPATIAL DATA
ST_list <- SplitObject(ST, split.by = "Slide")

# FUNCTION to run deconvolution
run_deconvolution <- function(ST, single_cell, slide){
	coords<- Map(function(x){x@coordinates[,4:5]}, ST@images)
	coords <- do.call(rbind, coords)
	rownames(coords) <- colnames(ST)
	colnames(coords) <- c("x", "y")


	counts <- ST@assays$SCT@counts
	rownames(counts) <- rownames(ST)
	colnames(counts) <- colnames(ST)

	spaceRNA <- SpatialRNA(coords, counts)



	myRCTD <- create.RCTD(spaceRNA, reference, max_cores = 1,
	                      CELL_MIN_INSTANCE = 15)
	print("create.RCTD")

	myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
	print("run.RCTD")


	dir.create(file.path(out_dir, file_name, slide))
	setwd(file.path(out_dir, file_name, slide))
	saveRDS(myRCTD,  paste0(file_name,slide,".rds"))
	resultsdir <- file.path(out_dir, file_name, slide)

	results <- myRCTD@results
	# normalize the cell type proportions to sum to 1.
	norm_weights = normalize_weights(results$weights) 
	cell_type_names <- myRCTD@cell_type_info$info[[2]] s
	spatialRNA <- myRCTD@spatialRNA


	plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights) 
	# Plots all weights for each cell type as in full_mode. 
	#'results/cell_type_weights.pdf'

	plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 
	# Plots the weights for each cell type as in doublet_mode. 
	# 'results/cell_type_weights_doublets.pdf'

	plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, 
                     results$results_df) 
	# Plots the number of confident pixels of each cell type in 'full_mode'.
	# 'results/cell_type_occur.pdf'

	plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)

}

dir.create(file.path(out_dir, file_name))
setwd(file.path(out_dir, file_name))

#Run deconvolution on all slides at once
run_deconvolution(ST, slide = "all", single_cell)

#Run deconvolution per slide 
#lapply(names(ST_list), function(slide) {
#  run_deconvolution(ST_list[[slide]], single_cell, slide)

#})


