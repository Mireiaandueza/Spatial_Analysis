# ADD RCTD RESULTS TO SEURAT SPATIAL DATA WITH SIGNATURES TO MERGE ALL 
# AVAILABLE INFO

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
library(ggridges)
library(tidyr)

# DIRECTORIES
spatial_data <- "./VISION"
rctd_data <- "./RCTD_INTEGRATION/"
setwd(rctd_data)


# RCTD RESULTS
CER <- readRDS(paste0("./CER/all/","CERall.rds"))
WT <- readRDS(paste0("./WT/all/","WTall.rds"))
i4F <- readRDS(paste0("./i4F/all/","i4Fall.rds"))
sc_OP <- readRDS(paste0("./sc_OP/all/","sc_OPall.rds"))
WT_CER_i4F_epithelial <- readRDS(paste0(
  "./WT_CER_i4F_epithelial/all/","WT_CER_i4F_epithelialall.rds"))


datasets <- list(WT = WT, CER = CER, i4F = i4F, sc_OP = sc_OP,
                 WT_CER_i4F_epithelial = WT_CER_i4F_epithelial)
names <- c("WT", "CER", "i4F", "sc_OP", "WT_CER_i4F_epithelial")

# ST DATA
setwd(spatial_data)
ST <- readRDS("ST_signatures.rds")


# Create a function to add metadata to ST
add_metadata <- function(dataset, name) {
  metadata <- dataset@results$weights
  first_type <- dataset@results$results_df$first_type
  col_names <- paste0(name, "_", colnames(metadata))
  
  ST <<- AddMetaData(object = ST, metadata = first_type, 
                     col.name = paste0(name, "_first_type"))
  ST <<- AddMetaData(object = ST, metadata = metadata, col.name = col_names)
}

# Apply the function to each dataset
Map(add_metadata, datasets, names)


#PLOT
Data <- c("WT_CER_i4F_epithelial_first_type","sc_OP_first_type" ,
          "i4F_first_type","CER_first_type","WT_first_type")

for (i in Data){
df_data <- ST@meta.data[i]
df_data$Genotype <- ST$Genotype
df <- data.frame(df_data)

df_percent <- eval(parse(text = paste0(
  "df_percent <- df %>% group_by( Genotype,",i,")%>%
  summarize(count = n()) %>%
  group_by(Genotype) %>%
  mutate(Percent = (count / sum(count)) * 100)")))

cols <- c("darkorchid4","hotpink","grey")
pdf(paste0(i,"_Stacked.pdf"))

# Create a stacked barplot
Pls <- eval(parse(text = paste0(
  "ggplot(df_percent, aes(x = ",i,", y = Percent, fill = Genotype))")))
Pls <- Pls +
  geom_bar(stat = "identity") #+ scale_fill_manual(values = cols)+
  labs(title = paste0("Cell Type Distribution by Genotype for " ,i," Serrano")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(Pls)
dev.off()
}

# SAVE ALL IN A .RDS FILE
saveRDS(ST, file = "ST_RCTD_annotated.rds")
