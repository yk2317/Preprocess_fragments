library(data.table)
library(stringr)
library(tidyverse)

bind_pos_neg <- function(data){
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
    return(df_wide)
}


make_matrix <- function(input_dir, output_dir){
    bin1 <- readRDS("/n/data1/hms/dbmi/park/ctDNA_loci_project/fragmentomics/shared_tools/length_bins.Rds")
    bin2 <- readRDS("/n/data1/hms/dbmi/park/ctDNA_loci_project/fragmentomics/shared_tools/combi_4mer.Rds")

    for(i in 2:3){
        patterns <- c("length.csv", "in5p.csv", "out5p.csv")
        segs <- c("length.csv", "in5p.csv","out5p.csv")

        files <- list.files(input_dir, pattern=patterns[i], full.names=TRUE) # Corrected variable from 'path' to 'input_dir'
        if(i==1){
            bin <- bin1
        }else{
            bin <- bin2
        }

        comb <- data.frame(matrix(nrow=length(bin), ncol=length(files)))
        rownames(comb) <- bin

        for(j in 1:length(files)){
            dat <- fread(files[j], header=FALSE)
         
            if(nrow(dat) != length(bin)){
                print("WARNING, mismatch")
            } 

            comb[,j] <- dat$V2
            colnames(comb)[j] <- str_remove(basename(files[j]), segs[i])
            #print(j)
        }
        
        final <- bind_pos_neg(comb)         

        out <- paste0(output_dir, "/", str_remove(patterns[i], ".csv"), "_matrix.csv")
        write.csv(final, out, quote=FALSE, row.names=FALSE)
    }
}

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]
make_matrix(input_dir, output_dir)
