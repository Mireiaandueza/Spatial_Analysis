# SCENIC PIPELINE FOR REGULON ANALYSIS

library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SCENIC)
library(scales)
library(AUCell)

dbFiles <- c("https://resources.aertslab.org/cistarget/databases/mus_musculus/\
             mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather",
"https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/\
mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather")

# DIRECTORIES
out_dir <- dir.create("SCENIC_HET")
setwd(out_dir)
dbfiles <- "./SCENIC/"

# DATA
ST_1 <- readRDS("./ST_files/ST1.rds")
ST_2 <- readRDS("./ST_files/ST2.rds")

# MERGE
ST <- merge(x =ST_1, y = ST_2, add.cell.ids = c("ST1","ST2"), merge.data = TRUE)
names(ST@images) <- c("ST1", "ST2")
ST <- JoinLayers(ST,  assay = "Spatial")

#SETTINGS
org <- "mgi" 
dbDir <- dbfiles
myDatasetTitle <- "SCENIC OR ST" 
data(defaultDbNames)
dbs <- defaultDbNames[[org]]

data(list="motifAnnotations_mgi_v9", package="RcisTarget")
motifAnnotations_mgi <- motifAnnotations_mgi_v9

scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, 
                                  datasetTitle=myDatasetTitle, nCores=10) 
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
Idents(ST) <- ST$Genotype
cellInfo <- data.frame(seuratCluster=Idents(ST))


exprMat <-ST@assays$Spatial$counts
exprMat <- as.matrix(exprMat)

# FILTERING
#genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
 #                          minCountsPerGene=3*.01*ncol(exprMat),
  #                         minSamples=ncol(exprMat)*.01)

#genesKept <- loadInt(scenicOptions, "genesKept")
#check
#interestingGenes <- c("Agr2", "Fetub", "Ins2", "Ctrb1")
#interestingGenes[which(!interestingGenes %in% genesKept)]

#Filter
#exprMat_filtered <- exprMat[genesKept, ]
#dim(exprMat_filtered)
#rm(exprMat)

#Run in all genes
exprMat_filtered <- exprMat
rm(exprMat)

runCorrelation(exprMat_filtered, scenicOptions)

# Norm
exprMat_filtered <- log2(exprMat_filtered+1) 

# Run GENIE3
runGenie3(exprMat_filtered, scenicOptions)

# Build and score the GRN
exprMat_log <- log2(exprMat_filtered +1)
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) 
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)


scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") 

# Export
export2loom(scenicOptions, exprMat_filtered)
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

# REGULON ANALYSIS
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

regulonActivity_byCellType <- sapply(split(rownames(cellInfo), 
                                           cellInfo$seuratCluster),
                                     function(cells) 
                                       rowMeans(getAUC(regulonAUC)[,cells]))

regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), 
                                             center = T, scale=T))
# HEATMAP ORDERED 
pdf("heatmap_SCENIC.pdf")
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, 
                        name="Regulon activity")
dev.off()

regulonActivity_byCellType <- as.data.frame(regulonActivity_byCellType)
regulonActivity_byCellType$absolute_difference <- abs(
  regulonActivity_byCellType$DD - regulonActivity_byCellType$WT)

regulons_sorted <- regulonActivity_byCellType[
  order(-regulonActivity_byCellType$absolute_difference), ]

# Select top 10 regulons with largest differences
top_10_regulons <- head(regulons_sorted, 10)
pdf("heatmap_top10_SCENIC.pdf")
ComplexHeatmap::Heatmap(top_10_regulons, 
                        name="Regulon activity")
dev.off()

## PLOTS

# GET OASIS/Creb3l1 expression
creb <-getAUC(regulonAUC)["Creb3l1 (416g)",]
ST$Creb3l1_regulon <-  creb
scale_range <- range(creb)

# PLOT OASIS/Creb3l1 REGULON
p1 <- SpatialFeaturePlot(ST, features = "Creb3l1_regulon", images = "ST1", 
                         pt.size.factor = 1.8, image.alpha = 0) + 
  theme(legend.position = "right")+
  scale_fill_gradientn(colors = c( "blue","darkseagreen4","khaki","red"),
                       limits = c(0.2,scale_range[2]), values = rescale(c(0.2,0.4,0.6,0.7)))
pdf("ST1_Creb3l1_regulon.pdf")
p1 
dev.off()

p1 <- SpatialFeaturePlot(ST, features = "Creb3l1_regulon", images = "ST2", 
                         pt.size.factor = 1.8, image.alpha = 0) + 
  theme(legend.position = "right")+
  scale_fill_gradientn(colors = c( "blue","darkseagreen4","khaki","red"),
                       limits = c(0.2,scale_range[2]), 
                       values = rescale(c(0.2,0.4,0.6,0.7)))
pdf("ST2_Creb3l1_regulon.pdf")
p1 
dev.off()

# PLOT Creb3l1 TF
ST<- NormalizeData(ST, assay="Spatial")
ST <- FindVariableFeatures(ST, assay="Spatial")
ST <- ScaleData(ST, assay="Spatial", features = rownames(ST))

p1 <- SpatialFeaturePlot(ST, features = "Creb3l1", images = "ST1",
                         pt.size.factor = 1.8, image.alpha = 0) + 
  theme(legend.position = "right")+
  scale_fill_gradientn(colors = c( "blue","darkseagreen4","khaki","red"),
                       limits = c(0,5), values = rescale(c(0,1,2,4.83)))
pdf("ST1_Creb3l1.pdf")
p1 
dev.off()

p1 <- SpatialFeaturePlot(ST, features = "Creb3l1", images = "ST2", 
                         pt.size.factor = 1.8, image.alpha = 0) + 
  theme(legend.position = "right")+
  scale_fill_gradientn(colors = c( "blue","darkseagreen4","khaki","red"),
                       limits = c(0,5), values = rescale(c(0,1,2,4.83)))
pdf("ST2_Creb3l1.pdf")
p1 
dev.off()

