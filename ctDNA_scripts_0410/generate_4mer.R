library(rtracklayer)
library(GenomicRanges)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)

freq <- function(path, outdir, direction){
  basename <- gsub(pattern = ".bed$", replacement = "", x = basename(path))
  out_csv <- paste0(outdir, "/", basename,"_", direction,  "_4mer_counts.csv")
  out_rds <- paste0(outdir, "/", basename,"_", direction, "_4mer_counts.Rds")

  granges <- import.bed(path)

  if(direction == "forward"){
    seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, granges)
  }
  if(direction == "reverse"){
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
      if(nchar(four_mer) == 4 && !grepl("N", four_mer)){
          counts[four_mer] <- counts[four_mer] + 1
      }
    }
    return(counts)
  }

  frequency_list <- lapply(as.character(seqs), count_4mers)
  print("up to count_4mer, lapply")
  #saveRDS(frequency_list, out_rds)

  expected_names <- names(frequency_list[[1]])
  if(!all(sapply(frequency_list, function(x) all(names(x) == expected_names)))) {
    stop("Not all elements have consistent names.")
  }

  frequency_df <- do.call(rbind, frequency_list)
  print("up to do.call(rbind)")
  row.names(frequency_df) <- paste("Region", seq_along(frequency_list))

  coordinates <- as.data.frame(seqnames(granges))
  coordinates$start <- start(granges)
  coordinates$end <- end(granges)
  frequency_df <- cbind(coordinates, frequency_df)

  write.csv(frequency_df, out_csv, row.names = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
outdir <- args[2]
direction <- args[3]
freq(path, outdir, direction)
