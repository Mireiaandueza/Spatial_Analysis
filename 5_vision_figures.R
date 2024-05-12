# PLOT SIGNATURE/FEATURE STATISTICS AND SIGNATURES/FEATURES IN SITU
# (code to use when KRAS samples are included is commented)

#LIBRARIES
library(devtools)
library(VISION)
library(Seurat)
library(ggplot2)
library(Matrix)
library(scales)
library(ggridges)
library(purrr)
library(ggpubr)

# DIRECTORIES
read_dir <- "./VISION/"
save_dir <- "./FIGURES/"
setwd(read_dir)
source("find_correlation.R")

#LOAD DATA
ST <- readRDS(file = paste0(read_dir,"ST_signatures.rds"))

#SET WORKING DIRECTORY
setwd(save_dir)

#FUNCTION DEFINITION

violinplot_genotype <- function(ST, firma){
	
	#comparisons <- list(c("DD", "HET"), c("HET", "WT"),c("DD", "WT"),c("KRAS WT", "WT"),c("KRAS DD", "DD"), c("KRAS DD", "KRAS WT"))
	
	#VIOLINPLOT by Genotype
	pdf(paste0(firma, "_violin_" ,".pdf"))
	plot <- VlnPlot(object = ST, features = firma, group.by = "Genotype" ,
	                pt.size = FALSE,y.max = 4) + NoLegend() + 
	  stat_compare_means(comparisons = comparisons, label = "p.signif") +
	  stat_summary(fun = median, geom='point', size = 25,
	               colour = "black", shape = 95)+ 
	  ggtitle(paste0(firma, " signature by genotype"))+ 
  	scale_fill_manual(values = cols1)
	print(plot)
	dev.off()

}

plt_feature_splitted_violin <- function(ST, firma){
	
	comparisons <- list(c("DD", "HET"), c("HET", "WT"),c("DD", "WT"))
	
	# FEATUREPLOT
	pdf(paste0(firma, ".pdf"))
	plot2 <- FeaturePlot(ST, features = firma, reduction = "umap") +
	  ggtitle(paste0("Expression of ",firma))
	print(plot2)
	dev.off()

	# VIOLIN PER GENOTYPE_TIMEPOINT
	pdf(paste0(firma, "_violin.pdf"))
	plt <- VlnPlot(object = ST, features = firma, group.by = "Genotype", 
	               split.by = "Genotype_Timepoint",pt.size = FALSE, 
	               y.max = 20)  + 
	  stat_compare_means(comparisons = comparisons, label = "p.signif") + 
	  coord_cartesian(ylim = c(0,2.5)) +
  	scale_fill_manual(values = cols2) +
	  ggtitle(paste0("Expression of ",firma)) +
	  xlab("Genotype and Timepoint") +
	  ylab("Signature score")
	print(plt)
	dev.off()

}

spatial_features <- function(ST, slide, feature){
	
	scale_range <- range(FetchData(object = ST, vars = feature))
	ST_sub <-subset(x = ST, subset = Slide == slide)

	# PLOT spatialfeatureplot
	pdf(paste0(slide,"_", feature, ".pdf"))
	plot2 <- SpatialFeaturePlot(ST_sub, features = feature, images = slide,
	                            pt.size.factor = 1.8, image.alpha = 0) +
	  theme(legend.position = "right") +
	  scale_fill_gradientn(colors = c( "blue","yellow","red"),
                   limits = c(scale_range[1],scale_range[2]), 
                   values = rescale(c(scale_range[1],
                                      (scale_range[1]+scale_range[2])/2,
                                      scale_range[2])))
	print(plot2)
	dev.off()

}

#################################


# SPATIAL FEATURE plot

#Signatures

firmas <- c("Signature_acinar","Signature_DD","DD_3months", "DD_6months", "DD_1.5years", "Metaplastic_Schlesinger","REACTOME_DIABETES",
"Prasad_et_al_Panin1_2_up","ER_Stress","PERK_UPR", "IRE1_UPR", "ATF6_UPR","ACUTE_INFLAMMATORY_RESPONSE", "Reprogrammed")


ST$seurat_clusters <- ST@active.ident

#ST$Genotype_Timepoint <- factor(ST$Genotype_Timepoint, levels = c("WT 3 months", "WT 6 months", "WT 1.5 years","HET 3 months", "HET 6 months","HET 1.5 years","DD 3 months", "DD 6 months", "DD 1.5 years", "KRAS WT 6 months","KRAS DD 6 months"))
ST$Genotype_Timepoint <- factor(ST$Genotype_Timepoint, levels = c("WT 3 months", "WT 6 months", "WT 1.5 years","HET 3 months", "HET 6 months","HET 1.5 years","DD 3 months", "DD 6 months", "DD 1.5 years"))

#ST$Genotype <- factor(ST$Genotype, levels = c("WT","HET","DD","KRAS WT","KRAS DD"))
ST$Genotype <- factor(ST$Genotype, levels = c("WT","HET","DD"))

Groups <- c("seurat_clusters", "Genotype_Timepoint", "Genotype")

# Colors 
#cols1 <- c("grey","hotpink","darkorchid4","skyblue", "skyblue4")
cols1 <- c("grey","mediumorchid1","darkorchid4")
#cols2 <- c("gray80","gray60", "gray40","pink","hotpink","deeppink2","mediumorchid1","darkorchid", "darkorchid4", "skyblue", "skyblue4")
cols2 <- c("gray80","gray60", "gray40","pink","hotpink","deeppink2","mediumorchid1","darkorchid", "darkorchid4")


# Apply violinplot to each signature
map(firmas, violinplot_genotype, ST = ST)

# Apply featureplot and violinplot divided by genotype-timepoint
map(firmas, plt_feature_splitted_violin, ST = ST)


#BACK TO SLIDES

# Create a new file for SLIDES (warning if it already exists)
mainDir <- getwd()
dir.create(file.path(mainDir, "SLIDES"))
setwd(file.path(mainDir, "SLIDES"))

# SPATIAL FEATURE plot
#names(ST@images) <- c("ST1", "ST2", "ST3", "ST4", "ST5", "ST6", "ST7", "ST8", "ST9", "ST10", "ST12")
names(ST@images) <- c("ST1", "ST2", "ST3", "ST4", "ST5", "ST6", "ST7", "ST8", "ST9", "ST10")

# Create colour palette
Idents(ST) <- ST$SCT_snn_res.0.2
ST$seurat_clusters <- Idents(ST)
default_palette <- hue_pal()(length(levels(ST$seurat_clusters)))
ST@meta.data$seurat_clusters_colors <- default_palette[ST@meta.data$seurat_clusters]

Signatures <- c("Signature_Acinar_modified","Signature_DD","ER_Stress")
#,"DD_3months", "DD_6months", "DD_1.5years", "Metaplastic_Schlesinger","REACTOME_DIABETES", "ACUTE_INFLAMMATORY_RESPONSE", #"Prasad_et_al_Panin1_2_up","PERK_UPR", "IRE1_UPR", "ATF6_UPR"
features <- c("Agr2", "Ins2","Cel" ,"Sst","Gcg") 
#"Ins1", "Sst","Gcg", "Neurog3","Ern2","Ern1", "Pdx1" )
slides <- names(ST@images)

# PLOT FEATURES/GENES IN SITU

# Create all combinations of feature and capture area
combinations <- expand.grid(feature = features, slide = slides)
combinations$slide <- as.character(combinations$slide)
combinations$feature <- as.character(combinations$feature)

# Apply function to each combination using pmap
pmap(combinations, spatial_features, ST = ST)

# PLOT SIGNATURES IN SITU

# Create all combinations of feature and capture area
combinations <- expand.grid(feature = Signatures, slide = slides)
combinations$slide <- as.character(combinations$slide)
combinations$feature <- as.character(combinations$feature)

# Apply function to each combination using pmap
pmap(combinations, spatial_features, ST = ST)


# PLOT CLUSTERS IN SITU

slides <- names(ST@images)
lapply(slides, function(slide) {
	ST_sub <- subset(x = ST, subset = Slide == slide)
	colors  <- default_palette[sort(match(unique(ST_sub@meta.data$seurat_clusters_colors), default_palette))]
	
	pdf(paste0(slide,"_cluster.pdf"))
	plot2 <- SpatialDimPlot(ST_sub, label = TRUE, label.size = 1.5, images = slide, pt.size.factor = 2, image.alpha = 0)+
 	scale_fill_manual(values= colors)+NoLegend()
	print(plot2)
	dev.off()

})



#CORRELATION PLOTS

DefaultAssay(ST) <- "Spatial"
ST <- NormalizeData(ST, assay="Spatial")
ST <- FindVariableFeatures(ST)
ST <- ScaleData(ST)

# Ins1 vs Gcg
find_correlation(ST, "Ins1", "Gcg", getwd())





