library(rtracklayer)
library(GenomicRanges)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)

freq <- function(path, outdir, dir){
  granges <- import.bed(path)
  
  if(dir == "forward"){
    seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, granges)
  }
  if(dir == "reverse"){
    seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, granges)
    seqs <- reverseComplement(seqs)
  }
  
  count_4mers <- function(sequence) {
    nucleotides <- c("A", "C", "G", "T")
    all_4mers <- expand.grid(nucleotides, nucleotides, nucleotides, nucleotides)
    four_mer_names <- apply(all_4mers, 1, paste, collapse = "")
    counts <- rep(0, length(four_mer_names))
    names(counts) <- four_mer_names
    
    for(i in 1:(nchar(sequence) - 3)){
      four_mer <- substring(as.character(sequence), i, i+3)
      if(nchar(sequence) > 4){
        counts[four_mer] <- counts[four_mer] + 1
      }
    }
    return(counts)  
  }
  
  frequency_list <- lapply(as.character(seqs), count_4mers)
  rds_path <- str_replace(outdir, "\\.csv$", ".Rds")
  saveRDS(frequency_list, rds_path)
 
  frequency_df <- do.call(rbind, frequency_list)
  row.names(frequency_df) <- paste("Region", seq_along(frequency_list))
  
  coordinates <- as.data.frame(seqnames(granges))
  coordinates$start <- start(granges)
  coordinates$end <- end(granges)
  #coordinates$gene <- df_bed$V4
  frequency_df <- cbind(coordinates, frequency_df)
  
  write.csv(frequency_df, outdir, row.names = FALSE)  
}

args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
outdir <- args[2]
dir <- args[3]
freq(path, outdir, dir)
