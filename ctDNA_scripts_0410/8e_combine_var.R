# duplex: /n/data1/hms/dbmi/park/DATA/ctDNA_datasets/MSK-GRAIL/fgbio_clara/step6_out/no_mark_duplicates_hg38/.PreProcessing/
# simplex: /n/data1/hms/dbmi/park/DATA/ctDNA_datasets/MSK-GRAIL/fgbio_clara/step6b_out/.PreProcessing/

combine <- function(path){
  types <- c("out5p_matrix.csv", "in5p_matrix.csv", "length_matrix.csv")
  
  for(type in types){
    type <- "out5p_matrix.csv"
    dat1 <- read.csv(paste0(path,"/duplex/matrix/variants/",type))
    dat2 <- read.csv(paste0(path,"/simplex/matrix/variants/",type))
    
    ind <- match(colnames(dat1), colnames(dat2))
    print("in duplex, but not in simplex")
    print(colnames(dat1)[is.na(ind)])
    ind <- match(colnames(dat2), colnames(dat1))
    print("in simplex, but not in duplex")
    print(colnames(dat2)[is.na(ind)])
    
    shared <- intersect(colnames(dat1), colnames(dat2))
    unique_dat1 <- setdiff(colnames(dat1), colnames(dat2))
    unique_dat2 <- setdiff(colnames(dat2), colnames(dat1))
    
    shared <- shared[shared != "X"]
    combined_common <- dat1[, shared] + dat2[, shared]
    
    df_unique_dat1 <- dat1[, unique_dat1]
    df_unique_dat2 <- dat2[, unique_dat2]
    
    combined_dat <- cbind(combined_common, df_unique_dat1, df_unique_dat2)
    outs <- paste0(path, "/combined/matrix/combined/", type)
    write.csv(combined_dat, outs)
    print(paste0("file created: ", outs))
  }
}

args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
combine(path)
