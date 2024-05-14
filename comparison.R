comparison <- function(ST , group, signature, out_dir){

# FUNCTION TO COMPUTE DIFFERENTIAL EXPRESSION FROM SPOTS, 
  # that have high vs low expresion of the signature. 

  	# Input:
		  # ST: Spatial transcritpmics seurat object
		  # group: genotype_timepoint subgroup to analyze
		  # signature: signature to compare
		  # out_dir: output directory
	# Output: It saves a cvs file in the ouput directory.

setwd(out_dir)

filtered_cells <- subset(ST, subset = Genotype_Timepoint == group)
filtered_cells <- eval(parse(text = paste0("subset(x = filtered_cells,
                                           subset = ",signature,"> 0)")))


signature_values <- FetchData(filtered_cells, 
                              vars = signature, layer = "data")[, 1]

DefaultAssay(filtered_cells) <- "Spatial"
pdf(paste0("Density", signature,gsub(" ", "", group), ".pdf"))

# Distribution plot of "SIGNATURE"
density_plot <- density(signature_values)
plot(density_plot, main = paste0("Distribution of " , signature," for ",
                                 group), xlab = signature)

# Quantiles
quantiles <- quantile(signature_values, c(0.25, 0.5, 0.75))
abline(v = quantiles, col = c("red", "green", "blue"), lty = 2)
legend("topright", legend = c("Q1", "Median", "Q3"), 
       col = c("red", "green", "blue"), 
       lty = 2)
dev.off()

# Down and up spots
Down <- eval(parse(text = paste0("subset(filtered_cells, subset =",
                                 signature,"< quantiles[[1]])")))
Up <- eval(parse(text = paste0("subset(filtered_cells, subset =",
                               signature,"> quantiles[[3]])")))
Down$Value <- "Down"
Up$Value <- "Up"

merged <- merge(Down, Up, add.cell.ids = c("Signature_Down", "Signature_Up"), 
                project = "merged_seurat")

# Find Markers
Idents(merged) <- "Value"
merged <- JoinLayers(merged)
pbmc.markers <- FindMarkers(merged, ident.1 = "Up", ident.2 = "Down")
pbmc.markers %>%
    slice_max(n = 20, order_by = avg_log2FC)

# Save
group <- gsub(" ", "", group)
write.table(pbmc.markers, 
            file =paste0(signature, group, "markers.csv"), 
            sep = ",")
}


# LIBRARIES
library(devtools)
library(Seurat)
library(ggplot2)
library(ggridges)
library(patchwork)
library(purrr)
library(scales)


# DIRECTORY
read_dir <- "./VISION"
setwd(read_dir)
out_dir < "./COMPARISONS"

# DATA
ST <- readRDS("ST_signatures.rds")

# COMPARISONS
comparison(ST , "DD 3 months", "Agr2", out_dir)
comparison(ST , "DD 6 months", "Agr2", out_dir)
comparison(ST , "DD 1.5 years", "Agr2", out_dir)

