library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)

plt_motif <- function(data, title){
  data <- data %>% mutate(alpha_group = ((row_number() - 1) %/% 16) %% 2)
  alpha_levels <- c(0.3,1)
  data$gp <- substr(data$V1, 1,1)
  
  p <- ggplot(data, aes(x=V1, y=V2, fill=gp)) + geom_bar(stat = "identity", aes(alpha = factor(alpha_group))) +
    scale_alpha_manual(values = alpha_levels) +
    theme(
      panel.background = element_rect(fill = "white", color=NA),  
      plot.background = element_rect(fill = "white", color=NA),
      panel.border = element_blank(), 
      panel.grid.major = element_line(color="white"),              
      panel.grid.minor = element_line(color="white"),               
      axis.text.x = element_blank(),  
      axis.ticks.x = element_blank(), 
      axis.text.y = element_blank(),  
      axis.ticks.y = element_blank(), 
      legend.position = "none"
    ) +    
    labs(x = "", y = "", title = title)
  print(p)
}

get_data <- function(file){
  dat1 <- fread(file, select=c(6,9,12,7))
  colnames(dat1) <- c("string", "original", "corrected", "outward")
  dat1$inward <- str_sub(dat1$string, 3,10)
  
  #---
  pos <- read.csv("~/Desktop/count_8mer_MSK_pos.csv", header=FALSE)
  top1 <- pos$V1[1:10] # flip1
  top2 <- pos$V1[1:1000] # flip2
  top3 <- pos$V1[1:100000] # flip3 
  #---
  
  flip <- function(seq) {
    split_seq <- strsplit(seq, "")[[1]]
    reversed_seq <- paste(rev(split_seq), collapse = "")
    return(reversed_seq)
  }
  
  dat1$top1 <- dat1$original
  dat1$top1[dat1$original %in% top1] <- dat1$inward[dat1$original %in% top1]
  
  dat1$top2 <- dat1$original
  dat1$top2[dat1$original %in% top2] <- dat1$inward[dat1$original %in% top2]
  
  dat1$top3 <- dat1$original
  dat1$top3[dat1$original %in% top3] <- dat1$inward[dat1$original %in% top3]
  
  return(dat1)
}

check <- function(file, column, outdir){
  
  dat1 <- get_data(file)
  
  col <- dat1[[column]]
  out5p <- sapply(str_sub(col,1,4), flip)
  in5p <- str_sub(col,5,8)
  
  tmp <- data.frame(col, out5p, in5p)
  
  # Count for out5p
  count_out5p <- data.frame(table(tmp[,2]))
  colnames(count_out5p) <- c("V1", "V2")
  title_out5p <- paste0(column, "_", "out5p")
  
  # Count for in5p
  count_in5p <- data.frame(table(tmp[,3]))
  colnames(count_in5p) <- c("V1", "V2")
  title_in5p <- paste0(column, "_", "in5p_", basename(file))
  
  # Start PDF device
  path <- paste0(outdir,"/",column, "_", basename(file),".pdf")
  pdf(path, width = 14, height = 7)
  par(mfrow=c(1,2)) # Set layout for 2 plots side by side
  
  # Plotting with plt_motif
  plt_motif(count_out5p, title_out5p)  
  plt_motif(count_in5p, title_in5p)
  
  dev.off() # Close the PDF device
}

args <- commandArgs(trailingOnly = TRUE)
#file <- "~/Desktop/filtered_MSK11172ACHraw_pos_pass.bed_pymaster.txt"
#outdir <- "~/Desktop"

file <- args[1]
outdir <- args[2]
check(file, "original", outdir)
check(file, "corrected", outdir) # selective correction 
check(file, "outward", outdir) 
check(file, "inward", outdir) 
check(file, "top1", outdir) # from original - messing up the most frequent 10
check(file, "top2", outdir) # from original - messing up the most frequent 1000
check(file, "top3", outdir) # from original - messing up the most frequent 100000



