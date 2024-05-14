find_correlation <- function(ST , feature1, feature2, out_dir){
	# Input:
		# ST: Spatial transcritpmics seurat object
		# feature1: feature in x-axis
		# feature2: feature in y-axis
		# out_dir: output directory
	# Output: correlation plot


	# CAUTION: LOG NORMALIZED SEURAT OBJECT
	setwd(out_dir)
	
	# Subset of spots with gene expression > 0
	filtered <- eval(parse(text = paste0(
	  "subset(x = ST, subset = ", feature1,"> 0)")))
	filtered <- eval(parse(text = paste0(
	  "subset(x = filtered, subset = ", feature2,"> 0)")))
	
	# Split by genotype
	splitted <- SplitObject(filtered, split.by ="Genotype")
	
	# Function to compute correlation and pvalue
	cor_test_func <- function(data) {
		DefaultAssay(data)<- "Spatial"
		exp_feature1 <- FetchData(data, vars = feature1)[, 1]
		exp_feature2 <- FetchData(data, vars = feature2)[, 1]
	  	cor_result <- cor.test(exp_feature1, exp_feature2)
	  	return(c(Pearson_Correlation = cor_result$estimate, 
	  	         P_Value = cor_result$p.value))
	}
	
	# Apply function to each group/genotype
	results <- sapply(splitted, cor_test_func)
	
	exp1 <- FetchData(filtered, vars = feature1)[, 1]
	exp2 <- FetchData(filtered, vars = feature2)[, 1]
	total_cor_result <- cor.test(exp1, exp2)
	
	
	# SCATTERPLOT
	pdf(paste0(feature1,"_", feature2, ".pdf"))
	plot2 <- FeatureScatter(filtered, feature1 = feature1, feature = feature2, 
	                        group.by = "Genotype",  cols = c("darkorchid4", 
	                                                         "mediumorchid1",
	                                                         "grey")) +
	  xlab(feature1) + 
	  ylab(feature2) + 
	  ggtitle(paste0("WT:", round(results[1,2],2),"  DD:", round(results[1,1],2),
	                 "  HET:",round(results[1,3],2)," Total:", 
	                 round(total_cor_result$estimate,2))) 
	print(plot2)
	dev.off()
}


# Function for plotting feature spatially in ALL slides for 2 genotypes
spatial_features <- function(ST, slide, feature1,feature2, out_dir){

	# Input:
		# ST: Spatial transcritpmics seurat object
		# slide: slide to plot
		# feature1: feature in x-axis
		# feature2: feature in y-axis
		# out_dir: output directory
	# Output: both features ploted on top of the slide, side by side
	
	library(scales)
	setwd(out_dir)
	exp_feature1 <- FetchData(ST, vars = feature1, layer = "data")[, 1]
	exp_feature2 <- FetchData(ST, vars = feature2, layer = "data")[, 1]
	
	scale_range1 <- range(exp_feature1)
	scale_range2 <- range(exp_feature2)

	ST_sub <-subset(x = ST, subset = Slide == slide)

	# Plots 
	pdf(paste0(slide,"_", feature1,"_", feature2, ".pdf"))
	plot1 <- SpatialFeaturePlot(ST_sub, features = feature1, images = slide, 
	                            pt.size.factor = 1.8, image.alpha = 0) + 
	  theme(legend.position = "top") +
	  scale_fill_gradientn(colors = c( "blue","yellow","red"),
	                       limits = c(scale_range1[1],scale_range1[2]), 
	                       values = rescale(c(scale_range1[1],
	                                          (scale_range1[1]+scale_range1[2])/2,
	                                          scale_range1[2])))
	plot2 <- SpatialFeaturePlot(ST_sub, features = feature2, images = slide, 
	                            pt.size.factor = 1.8, image.alpha = 0) + 
	  theme(legend.position = "top") +
	  scale_fill_gradientn(colors = c( "blue","yellow","red"),
	                       limits = c(scale_range2[1],scale_range2[2]), 
	                       values = rescale(c(scale_range2[1],
	                                          (scale_range2[1]+scale_range2[2])/2,
	                                          scale_range2[2])))
 	# Combine two features  
	combined_plots <- plot1 + plot2
	print(combined_plots)
	dev.off()

}

