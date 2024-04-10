library(stringr)
library(readr)
library(tidyverse)
library(ggplot2)

get_df <- function(dat, lab, pat){
  sample_ids <- unique(sub("-cfDNA.*", "", basename(dat)))
  
  if(pat == "length.csv"){
    df <- data.frame(matrix(nrow=37, ncol=length(sample_ids)))
  }else{
    df <- data.frame(matrix(nrow=256, ncol=length(sample_ids)))
  }

  for(i in 1:length(sample_ids)){
    sub <- dat[grepl(basename(sample_ids[i]), dat)]
    sub_data <- lapply(sub, function(f) {
      read_csv(f, col_names = FALSE, show_col_types=FALSE)
    })
    
    second_columns <- lapply(sub_data, function(df) df[[2]])
    
    combined_column <- Reduce(`+`, second_columns)
    
    df[,i] <- combined_column 
    colnames(df)[i] <- paste0(basename(sample_ids[i]),"_",lab)
  }
  return(df)
}

bind_pos_neg <- function(data, bin){
  data$X <- bin #
  df_long <- data %>%
    pivot_longer(cols = -X, names_to = "sample", values_to = "value") %>%
    mutate(id = gsub("_(neg|pos)*", "", sample), #
           type = ifelse(grepl("neg", sample), "neg", "pos")) %>%
    select(-sample) %>%
    arrange(id, type)
  
  df_summed <- df_long %>%
    group_by(id, X) %>%
    summarise(sum_value = sum(value, na.rm = TRUE), .groups = 'drop')
  
  df_wide <- df_summed %>%
    pivot_wider(names_from = id, values_from = sum_value)
  
  #ind <- match(bin, df_wide$X)
  #df_wide <- df_wide[ind,]
  return(df_wide)
}

run <- function(path, outdir){
  pats <- c("out5p.csv", "in5p.csv", "length.csv")
  for(pat in pats){
    files <- list.files(path, pattern=pat, full.names=TRUE)
    neg <- files[str_detect(files, "neg")] 
    pos <- files[str_detect(files, "pos")]
    
    neg_df <- get_df(neg,"neg", pat)
    pos_df <- get_df(pos,"pos", pat)
    
    df <- cbind(neg_df, pos_df)
    
    if(pat == "length.csv"){
      bin <- readRDS("/n/data1/hms/dbmi/park/ctDNA_loci_project/fragmentomics/shared_tools/length_bins.Rds")
      bin <- bin[2:length(bin)]
    }else{
      bin <- readRDS("/n/data1/hms/dbmi/park/ctDNA_loci_project/fragmentomics/shared_tools/combi_4mer.Rds")
    }
    
    combined <- bind_pos_neg(df, bin)
    save_path <- paste0(outdir,"/",str_remove(pat, ".csv"),"_matrix.csv")
    write.csv(combined, save_path, row.names=FALSE)
    print(paste0("saved: ",save_path))
  }  
}

args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
outdir <- args[2]
#path <- "/n/data1/hms/dbmi/park/ctDNA_loci_project/fragmentomics/Wyatt/allreads/output_baseline_shifted/csv"
#outdir <- "/n/data1/hms/dbmi/park/ctDNA_loci_project/fragmentomics/Wyatt/allreads/matrix/baseline_shifted"
run(path, outdir)
