convert_gmt_files <- function(file, output_file){

  # Converts a gmt file containing human genes to mouse genes
  # Saves the new file in directory output_file
  
  # LIBRARIES
  library(devtools)
  library(Seurat)
  library(dplyr)
  library(qusage)
  
  # READ
  Signature <- read.gmt(file)
  
  # CONVERT
  mouse_signature <- convert_human_to_mouse(Signature[[1]])
  
  # SAVE
  
  gmt_file <- file(output_file, "w")
  cat(names(Signature)[1], "\t", sep = "", file = gmt_file)
      
  # Write gene symbols
  cat(paste(mouse_signature, collapse = "\t"), "\n", sep = "", file = gmt_file)
    
  # Close the connection
  close(gmt_file)
}

convert_human_to_mouse <- function(gene_list) {
    
    # Converts a list of human genes to mouse ortholog genes from Biomart
    output = c()
    mouse_human_genes = read.csv(
      "https://www.informatics.jax.org/downloads/reports/\
      HOM_MouseHumanSequence.rpt",sep="\t")

    for(gene in gene_list) {
          class_key = (mouse_human_genes %>% filter(
            Symbol == gene & Common.Organism.Name == "human"))[['DB.Class.Key']]
          
          if( !identical(class_key, integer(0)) ) {
            human_genes = (mouse_human_genes %>% filter(
              DB.Class.Key == class_key & 
                Common.Organism.Name=="mouse, laboratory"))[,"Symbol"]
            
            for(human_gene in human_genes) {
                output = rbind(c(human_gene), output)
            }
          }
     }
     return (output)
}
