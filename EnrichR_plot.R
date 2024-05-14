# CODE TO PLOT ENRICHR RESULTS FROM AGR2 HIGH VS LOW SPOTS

# LIBRARIES
library(ggplot2)
library(lattice)
library(dplyr)
library(patchwork)
library(ggpubr)

# DIRECTORY
dir <- "./AGR_HIG_VS_LOW/"
setwd(dir)


# Functions
genratio_count <- function(data){
  overlap_values <- strsplit(as.character(data$Overlap), "/")
  # Convert the values to numeric
  overlap_values <- sapply(overlap_values, as.numeric)
  # Compute the GeneRatio
  GeneRatio <- overlap_values[1, ] / overlap_values[2, ]
  Count <- overlap_values[1, ] 
  data$Count <- Count
  data$GeneRatio <- GeneRatio
  return(data)}

files <- c( "3m_UP_DD_Reactome_2022_table.txt",
            "6m_UP_DD_Reactome_2022_table.txt")

# files <- c( "3m_DN_DD_Reactome_2022_table.txt",
#             "6m_DN_DD_Reactome_2022_table.txt",
#             "1.5y_DN_DD_Reactome_2022_table.txt")
#files <- c( "3m_UP_DD_BioPlanet_2019_table.txt","6m_UP_DD_BioPlanet_2019_table.txt")

#files <- c( "3m_DN_DD_BioPlanet_2019_table.txt",
#            "6m_DN_DD_BioPlanet_2019_table.txt",
#           "1.5y_DN_DD_BioPlanet_2019_table.txt")



genes_3m <- genratio_count(read.table(file = files[1], 
                                      header = TRUE, sep ="\t"))
genes_6m <- genratio_count(read.table(file =  files[2],
                                      header = TRUE, sep ="\t"))
#genes_18m <- genratio_count(read.table(file =  files[3], 
                                      #header = TRUE, sep ="\t"))

genes_3m$time <- "3 months"
genes_6m$time <- "6 months"
#genes_18m$time <- "1.5 years"


genes <- bind_rows(genes_3m,genes_6m)
#genes <- bind_rows(genes_3m,genes_6m,genes_18m)

top_10 <- genes %>%
  group_by(time) %>%
  slice_min(order_by = Adjusted.P.value, n = 5) 

genes_10<-genes[genes$Term%in%top_10$Term,]
#genes_10$time <- factor(genes_10$time, 
#levels = c("3 months", "6 months", "1.5 years"))
genes_10$time <- factor(genes_10$time, levels = c("3 months", "6 months"))
genes_10 <-  arrange(genes_10,Adjusted.P.value)


data3 <- genes_10[genes_10$time == "3 months",]
missing_terms <- unique(genes_10$Term)[!(unique(genes_10$Term) %in% data3$Term)]
missing_rows <- data.frame(Term = missing_terms, Overlap  =NA,P.value=NA ,
                           Adjusted.P.value = NA, 
                           Old.P.value = NA, Old.Adjusted.P.value = NA,
                           Odds.Ratio = NA, Combined.Score = NA,
                           Genes = NA, Count = NA, GeneRatio = 0,
                           time = NA, stringsAsFactors = FALSE)

# Combine data with missing rows
data3 <- rbind(data3, missing_rows)
data3 <- arrange(data3, GeneRatio)
data3$Term <- factor(data3$Term, levels = unique(data3$Term))

data6 <- genes_10[genes_10$time == "6 months",]
missing_terms <- unique(genes_10$Term)[!(unique(genes_10$Term) %in% data6$Term)]
data6 <- arrange(data6, GeneRatio)
data6$Term <- factor(data6$Term, levels = unique(data3$Term))



#data18 <- genes_10[genes_10$time == "1.5 years",]

#missing_terms <- unique(genes_10$Term)[!(unique(genes_10$Term) %in% data18$Term)]
#missing_rows <- data.frame(Term = missing_terms, Overlap  =NA,P.value=NA ,
                          # Adjusted.P.value = NA, 
                          # Old.P.value = NA, Old.Adjusted.P.value = NA,
                          # Odds.Ratio = NA, Combined.Score = NA,
                          # Genes = NA, Count = NA, GeneRatio = 0,
                          # time = NA, stringsAsFactors = FALSE)
#data18 <- rbind(data18, missing_rows)
#data18 <- arrange(data18, GeneRatio)
#data18$Term <- factor(data18$Term, levels = unique(data3$Term))

# plot: dot plot
create_plot <- function(data, labels, title){
  plt <- ggplot(data =data, aes(x = GeneRatio, y = Term, 
                        color = Adjusted.P.value, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + scale_y_discrete(labels = labels) + xlim(c(0,0.3))
  theme(axis.text.x = element_text(size = 14), 
          axis.text.y = element_text(size = 14),
          axis.title = element_blank()) +   
    ylab("") + 
  xlab("")+ggtitle(title)
  return(plt)
}

plt3 <- create_plot(data3, labels = data3$Term, title = "3 months")
plt6 <- create_plot(data6, labels = NULL, title = "6 months")
#plt18<- create_plot(data18, labels = NULL, title = "1.5 years")
ggarrange(
  #plt3, plt6, plt18,
  plt3, plt6,
  common.legend = TRUE, legend = "right", widths = c(4.5,1),ncol = 2
)







