norm_and_clustering <- function(data){
	
	# FUNCTION FOR SCT TRANSFORM NORMALIZATION, UMAP AND CLUSTERING

	# Input:
		# data: data to transform (Seurat object)
	# Output: Seurat object with SCT normalization, UMAP reduction and clustered
	
	# ST TRANSFORM
	#Remove subset with no UMIS and correct errors
	ST_SCT <- subset(data, nCount_Spatial>0)

	# NORMALIZE WITH SCT
	ST_SCT <- SCTransform(ST_SCT, assay = "Spatial", verbose = TRUE, 
	                      vars.to.regress = c("nCount_Spatial", "Sample"), 
	                      return.only.var.genes = FALSE)

	# Dimensionality reduction, clustering, and visualization
	ST_SCT <- RunPCA(ST_SCT, assay = "SCT",slot = "scale.data", verbose = FALSE)
	ST_SCT <- FindNeighbors(ST_SCT, assay = "Spatial",
	                        reduction = "pca", dims = 1:30)

	# CLUSTREE SCT
	for (res in seq(0,1,by = 0.1)) {
    		ST_SCT <- FindClusters(ST_SCT, resolution = res, algorithm = 3)
	}

	ST_SCT <- RunUMAP(ST_SCT, reduction = "pca", dims = 1:30)

	# Plot clustree
	pdf("clustree_st_SCT.pdf" , width = 8, height = 10)
	print(clustree(ST_SCT, prefix = "SCT_snn_res.", node_colour = "nCount_Spatial", node_colour_aggr = "mean"))
	dev.off()
	
	# CLUSTERING 0.3
	Idents(ST_SCT) <- ST_SCT$SCT_snn_res.0.3

	p1 <- DimPlot(ST_SCT, reduction = "umap", label = TRUE)
	pdf("clustering_0.3.pdf")
	print(p1 )
	dev.off()

	# CLUSTERING 0.2
	Idents(ST_SCT) <- ST_SCT$SCT_snn_res.0.2

	p1 <- DimPlot(ST_SCT, reduction = "umap", label = FALSE)
	pdf("clustering_0.2.pdf")
	print(p1 )
	dev.off()

	return(ST_SCT) #Return normalized and clustered object
}
