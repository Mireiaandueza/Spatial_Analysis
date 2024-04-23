#CODE FOR ADDING METADATA FROM LOUPE BROWSER

#LIBRARIES
library(Matrix)
library(Seurat)
library(hdf5r)


# Set path variables
cloupe_files <- "./CLOUPES"
spaceranger_files <- "./spaceranger_out"
save_dir <- "./ST_files"

# Samples and files
samples <- c("ST1", "ST2", "ST3", "ST4", "ST5", "ST6", "ST7", "ST8","ST09", "ST10", "ST11", "ST12")

files <- c("ST1-PER226-189","ST2-PER182-185","ST3-PER158-157-155","ST4-PER174-152","ST5-PER104-111","ST6-PER124-106","ST7-PER099-112","ST8-PER197-222","ST09_PEO688-PER135","ST10_PER156-160","ST11_PER226-189","ST12_PIR615-636")

study <- cbind(samples = samples, files = files)

#Loop 
for (i in 1:nrow(study)) {

	# New directory
	dir <- paste0(cloupe_files, study[[i,1]])
	setwd(dir)

	# Load metadata
	Genotype <- read.csv("Genotype.csv")
	Timepoint <- read.csv("Timepoint.csv")
	Genotype_Timepoint <- read.csv("Genotype_Timepoint.csv")
	Sample <- read.csv("Sample.csv")


	# LOAD SPACERANGER DATA H5
	path = paste0(spaceranger_files, study[[i,2]],"/outs/")
	ST <- Load10X_Spatial(data.dir = path, filename ="filtered_feature_bc_matrix.h5",assay = "Spatial")
	print(dir)

	# Add metadata to object 
	ST_data <- AddMetaData(
  		object = ST,
  		metadata = Genotype$Genotype,
  		col.name = 'Genotype'
	)
	ST_data <- AddMetaData(
  		object = ST_data,
  		metadata = Timepoint$Timepoint,
  		col.name = 'Timepoint'
	)
	ST_data <- AddMetaData(
  		object = ST_data,
  		metadata = Genotype_Timepoint$Genotype_Timepoint,
  		col.name = 'Genotype_Timepoint'
	)
	ST_data <- AddMetaData(
  		object = ST_data,
  		metadata = Sample$Sample,
  		col.name = 'Sample'
	)

	# Save RDS object
	if (!dir.exists(save_dir)) {
    		dir.create(save_dir)
  	}

	filename <- file.path(save_dir, paste0(study[[i, 1]], ".rds"))
	saveRDS(ST_data, file = filename)
}

