library(data.table)
library(stringr)
library(tidyverse)

bind_pos_neg <- function(data, bin){
  data$X <- rownames(data)
  df_long <- data %>%
    pivot_longer(cols = -X, names_to = "sample", values_to = "value") %>%
    mutate(id = gsub("_(neg|pos)_.*", "", sample), # Extract common identifier
           type = ifelse(grepl("neg", sample), "neg", "pos")) %>%
    select(-sample) %>%
    arrange(id, type)
  
  df_summed <- df_long %>%
    group_by(id, X) %>%
    summarise(sum_value = sum(value, na.rm = TRUE), .groups = 'drop')
  
  df_wide <- df_summed %>%
    pivot_wider(names_from = id, values_from = sum_value)
  
  ind <- match(bin, df_wide$X)
  df_wide <- df_wide[ind,]
  return(df_wide)
}

normalize <- function(sub, dataset, strand, pattern){
  
  if(dataset == "MSK"){
    if(strand=="pos"){
      if(pattern == "in5p.csv"){ref <- read.csv("/n/data1/hms/dbmi/park/ctDNA_loci_project/fragmentomics/shared_tools/panel_4mer_counts/forward.csv")}
      if(pattern == "out5p.csv"){ref <- read.csv("/n/data1/hms/dbmi/park/ctDNA_loci_project/fragmentomics/shared_tools/panel_4mer_counts/forward_flip.csv")}
    }else if(strand=="neg"){
      if(pattern == "in5p.csv"){ref <- read.csv("/n/data1/hms/dbmi/park/ctDNA_loci_project/fragmentomics/shared_tools/panel_4mer_counts/rev_comp.csv")}
      if(pattern == "out5p.csv"){ref <- read.csv("/n/data1/hms/dbmi/park/ctDNA_loci_project/fragmentomics/shared_tools/panel_4mer_counts/forward_base.csv")}      
    }
  }
  if(dataset == "Gerberg"){
    if(strand=="pos"){
      if(pattern == "in5p.csv"){ref <- read.csv("/n/data1/hms/dbmi/park/ctDNA_loci_project/fragmentomics/shared_tools/panel_4mer_counts/Gerberg_forward.csv")}
      if(pattern == "out5p.csv"){ref <- read.csv("/n/data1/hms/dbmi/park/ctDNA_loci_project/fragmentomics/shared_tools/panel_4mer_counts/Gerberg_forward_flip.csv")}
    }else if(strand=="neg"){
      if(pattern == "in5p.csv"){ref <- read.csv("/n/data1/hms/dbmi/park/ctDNA_loci_project/fragmentomics/shared_tools/panel_4mer_counts/Gerberg_rev_comp.csv")}
      if(pattern == "out5p.csv"){ref <- read.csv("/n/data1/hms/dbmi/park/ctDNA_loci_project/fragmentomics/shared_tools/panel_4mer_counts/Gerberg_forward_base.csv")}      
    }
  }
  
  if(!identical(ref$x, rownames(sub))){
    print("warning check ordering of normalization factor")
  }

  norm_sub <- sweep(sub, 1, ref$y, "/")
  scale_factor <- colSums(sub)/colSums(norm_sub)
  
  scaled_matrix <- sweep(norm_sub, 2, scale_factor, '*')
  return(scaled_matrix)
}

make_matrix <- function(input_dir, output_dir, dataset, param1, param2){
  
  param1 <- tolower(param1) == "true"
  param2 <- tolower(param2) == "true"
  
  bin1 <- readRDS("/n/data1/hms/dbmi/park/ctDNA_loci_project/fragmentomics/shared_tools/length_bins.Rds")
  bin1 <- bin1[2:length(bin1)]
  bin2 <- readRDS("/n/data1/hms/dbmi/park/ctDNA_loci_project/fragmentomics/shared_tools/combi_4mer.Rds")
  
  for(i in 1:3){
    patterns <- c("length.csv", "in5p.csv", "out5p.csv")
    segs <- c("length.csv", "in5p.csv","out5p.csv")
    
    files <- list.files(input_dir, pattern=patterns[i], full.names=TRUE) # Corrected variable from 'path' to 'input_dir'
    if(i==1){bin <- bin1}else if(i %in% c(2,3)){bin <- bin2}
    
    comb <- data.frame(matrix(nrow=length(bin), ncol=length(files)))
    rownames(comb) <- bin
    
    for(j in 1:length(files)){
      dat <- fread(files[j], header=FALSE)
      
      if(nrow(dat) != length(bin)){print("WARNING, mismatch")}
      
      comb[,j] <- dat$V2
      colnames(comb)[j] <- str_remove(basename(files[j]), segs[i])
    }
    
    if(param1){
      if(patterns[i] != "length.csv"){
        pos <- comb[,str_detect(colnames(comb), "pos")]
        pos_norm <- normalize(pos, dataset, "pos", patterns[i])
        
        neg <- comb[,str_detect(colnames(comb), "neg")]
        neg_norm <- normalize(neg, dataset, "neg", patterns[i])
        
        norm <- cbind(X=rownames(comb), pos_norm, neg_norm)
      }else{
        norm <- cbind(X=rownames(comb), comb)
      }
      out <- paste0(output_dir, "/", str_remove(patterns[i], ".csv"), "_norm_matrix.csv")
    }else{
      print("Skipping normalize")
      norm <- cbind(X=rownames(comb), comb)
      out <- paste0(output_dir, "/", str_remove(patterns[i], ".csv"), "_matrix.csv")
    }
    
    if(param2){
      final <- bind_pos_neg(norm, bin)
    }else{
      colnames(norm) <- str_replace(colnames(norm), "(neg|pos).*", "\\1")
      final <- norm
    }
    
    write.csv(final, out, quote=FALSE, row.names=FALSE)
  }
}

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]
dataset <- args[3] #MSK or Gerberg
param1 <- args[4] #normalize
param2 <- args[5] #combine pos/neg
make_matrix(input_dir, output_dir, dataset, param1, param2)
