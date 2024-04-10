library(stringr)

normalize <- function(path, type){
  paths <- paste0(path, c("/out5p_matrix.csv", "/in5p_matrix.csv"))

  for(file in paths){
    dat <- read.csv(file)
    #--- normalization factor from expected 4mer distribution 
    if(type=="MSK"){
      norm <- read.csv("/n/data1/hms/dbmi/park/ctDNA_loci_project/fragmentomics/shared_tools/MSK_hg38_4mer.csv")
    }
    if(type=="Gerberg"){
      norm <- read.csv("/n/data1/hms/dbmi/park/ctDNA_loci_project/fragmentomics/shared_tools/Gerberg_hg38_4mer.csv")
    }
    
    norm_1 <- norm[,1:3]
    norm_2 <- norm[,4:ncol(norm)]
    ind <- order(colnames(norm_2))
    norm_2 <- norm_2[, ind]
    norm_factor <- colSums(norm_2)
    #---
    
    norm_matrix <- sweep(dat[,2:ncol(dat)], 1, norm_factor, "/")
    scale_factor <- colSums(dat[,2:ncol(dat)])/colSums(norm_matrix)
    
    scaled_matrix <- sweep(norm_matrix, 2, scale_factor, '*')
    
    scaled_matrix <- cbind(dat$X, scaled_matrix)
    colnames(scaled_matrix)[1] <- X
    write.csv(scaled_matrix, paste0(str_remove(basename(file), "_matrix.csv"), "_norm_matrix.csv"), quote=FALSE, row.names=FALSE)  
  }
}

args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
type <- args[2]
normalize(path, type)
