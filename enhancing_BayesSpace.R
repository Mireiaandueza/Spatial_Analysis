# RUN BAYESSPACE ON ST DATA TO ENHANCE SPOT RESOLUTION
# Computes BayesSpace on each slide separately

#LIBRARIES
library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)
library(Matrix)
library(Seurat)
library(dplyr)
library(purrr)


# DIRECTORIES
spatial_data <- "./ST_FILES" # Files with no preprocessing
setwd(spatial_data)
ST <- readRDS("ST_signatures.rds") 
# Create a new directory for results
out_dir <- dir.create(file.path(getwd(), "BayesSpace"))


# Files and clusters
files <- c("ST1.rds","ST2.rds","ST3.rds","ST4.rds", "ST5.rds", 
           "ST6.rds","ST7.rds", "ST8.rds","ST09.rds","ST10.rds",
           "ST11.rds","ST12.rds")
clusters<-c(8,7,6,8,8,7,4,4,7,6,6,8) #Decided after running qTune() individually

data <- data.frame(file = files, clusters = clusters)
data$file <- as.character(data$file)
data$clusters <- as.character(data$clusters)

# FUNCTION
run_bayes <- function(file_dir, out_dir, file, clusters){

	# Input:
		# file_dir: input file directory
		# file: input file 
		# out_dir: output file directory
		# clusters: number of clusters
	# Output: enhanced seurat object saved in out_dir
	
	print(paste0(file_dir, file))
	ST <- readRDS(paste0(file_dir, file))
	slide <- sub("\\.rds$", "", file)

	names(ST@images) <- slide
	# Filtering on a single file due to other unwanted sample and ganglion
	if (slide == "ST9"){
	ST <- subset(x = ST, subset = Sample == "PEO688", invert = TRUE)
	ST <- subset(x = ST, subset = Genotype == "Ganglion", invert = TRUE)}

	#Convert to SCE
	diet.seurat = Seurat::DietSeurat(ST, graphs = "pca")
	sce = as.SingleCellExperiment(diet.seurat)

	colData(sce) = cbind(colData(sce),ST@images[[slide]]@coordinates)
	sce <- spatialPreprocess(sce,platform = "Visium", n.PCs = 50, log.normalize = T) 	
	
	# Uncomment if number of clusters not decided
	#sce <- qTune(sce, qs=seq(2, 15))
	#pdf(paste0(slide,"_likelihood.pdf"))
	#pls<-qPlot(sce)
	#print(pls)
	#dev.off()
	
	n_cluster <- as.numeric(clusters)
	sce <- spatialCluster(sce, nrep = 1000, burn.in = 100, q = n_cluster )

	# Plot
	pdf(paste0(slide,"_cluster.pdf"))
	plt<-clusterPlot(sce)
	print(plt)
	dev.off()

	# ENHANCED RESOLUTION OF CLUSTERS
	print(paste0("Enhancing clusters for slide ", slide))
	sce.enhanced <- spatialEnhance(sce, q= n_cluster, 
                                    model="t", gamma=2,
                                    jitter_prior=0.3, jitter_scale=3.5,
                                    nrep=1000, burn.in=100)

	print("Enhanced!")
	pdf(paste0(slide,"_Enhanced.pdf"))
	plt <-clusterPlot(sce.enhanced)
	print(plt)
	dev.off()

	# ENHANCED RESOLUTION OF MARKERS
	print(paste0("Enhancing markers for slide ", slide, " can take a while..."))
	markers <- rownames(ST)
	
	sce.enhanced <- enhanceFeatures(sce.enhanced, sce,
                                     feature_names=markers,
                                     nrounds=0)
	print("Markers enhanced!")

	#CREATE SEURAT ENHANCED OBJECT
  	seurat_enh <- as.Seurat(sce.enhanced, counts = NULL, 
				data = "logcounts",project = "Spatial" )

	#Add metadata
	A <- rbind(ST@meta.data, ST@meta.data, ST@meta.data, 
		   ST@meta.data, ST@meta.data, ST@meta.data)
	seurat_enh <- AddMetaData(seurat_enh, A)
	seurat_enh@images[[slide]] <- ST@images[[1]] 
	
	#Save
	setwd(out_dir)
	print("Saving data")
	seurat_enh <- SetIdent(seurat_enh, value = "spatial.cluster")
	saveRDS(seurat_enh, file = paste0(slide,"_enhanced.rds"))

}

#RUN
Results<- pmap(data, run_bayes, file_dir = spatial_data, out_dir = out_dir)

