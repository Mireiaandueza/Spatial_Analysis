# RUN VISION SIGNATURES AND CELL CYCLE 


#LIBRARIES
library(devtools)
library(VISION)
library(Seurat)
library(ggplot2)
library(Matrix)
library(scales)
library(ggridges)
library(dplyr)


# DIRECTORIES
samples_dir <- "./ST_files/"
sig_dir <- "./signatures/"
out_dir <- "./VISION/"

# LOAD DATA
setwd(samples_dir)
ST <- readRDS(file = "ST_merge.rds")


#SET WORKING DIRECTORY
setwd(out_dir)


#CREATE PROJECTION
ST_projection <- ST@reductions$umap@cell.embeddings




# Define function
Runvision <- function(ST,signature, col_name) {
	
	
	ST_vis <- Vision(ST@assays$SCT$data,signatures = signature, 
	                 meta = ST@meta.data, pool=FALSE, 	
	                 min_signature_genes=1, sig_gene_threshold = 0.0001)

	ST_vis <- addProjection(ST_vis, "UMAP", ST_projection)
	ST_vis <- analyze(ST_vis)

	ST <- AddMetaData(
  		object = ST,
  		metadata = ST_vis@SigScores,
  		col.name = col_name)
	return(ST)

}

# RUN VISION SIGNATURES

#Define signature filenames
Signatures <- c("Metaplastic_Schlesinger_et_al.gmt",  
                "Prasad_et_al_Panin1_2_up.gmt", "ERSTRESS.gmt",
                "DD_overall_Timecourse.gmt",
                "GOBP_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE.v2023.2.Mm.gmt", 
                "GOBP_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE.v2023.2.Mm.gmt", 
                "GOBP_ATF6_MEDIATED_UNFOLDED_PROTEIN_RESPONSE.v2023.2.Mm.gmt",
                "REACTOME_DIABETES_PATHWAYS.v2023.2.Hs.gmt", 
                "KEGG_TYPE_II_DIABETES_MELLITUS.v2023.2.Hs.gmt", 
                "KEGG_TYPE_I_DIABETES_MELLITUS.v2023.2.Hs.gmt", 
                "KEGG_MATURITY_ONSET_DIABETES_OF_THE_YOUNG.v2023.2.Hs.gmt",
                "HALLMARK_PANCREAS_BETA_CELLS.v2023.2.Hs.gmt", 
                "GOBP_ACUTE_INFLAMMATORY_RESPONSE.v2023.2.Hs.gmt", 
                "geneset_acinar_enzymes_modified.gmt","reprogramed.gmt")

#Defin signature names for data column
sig_names <- c('Metaplastic_Schlesinger', 'Prasad_et_al_Panin1_2_up', 
               'ER_Stress','Signature_DD','PERK_UPR', 'IRE1_UPR', 'ATF6_UPR',
               'REACTOME_DIABETES', 'KEGG_DIABETES_II', 'KEGG_DIABETES_I', 
               'KEGG_MATURITY_ONSET_DIABETES','BETA_CELLS', 
               'ACUTE_INFLAMMATORY_RESPONSE', 'Signature_Acinar','Reprogrammed')


# RUN signatures
for (i in seq_along(Signatures)) {
  ST <- Runvision(ST, paste0(sig_dir,"/",Signatures[i]), sig_names[i])

}


# SAVE
setwd(out_dir)
saveRDS(ST, file = "ST_signatures.rds")

