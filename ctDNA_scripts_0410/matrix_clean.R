library(stringr)

dir <- "/n/data1/hms/dbmi/park/ctDNA_loci_project/fragmentomics/MSK_hg38_old/matrix_var_wbc"
files <- list.files(dir, pattern=".csv", full.names = T)

for(i in 1:length(files)){
  name <- str_remove(basename(files[i]), ".csv")
  dat <- read.csv(files[i])
  x <- dat$X
  dat <- dat[,2:ncol(dat)]
  
  #sums <- colSums(dat)
  #sums[sums<500]
  #sum(sums<500)
  
  dat <- dat[,colSums(dat)>500]
  dat <- cbind(x, dat)
  path <- paste0(dir, "/", name, "_cleaned.csv")
  write.csv(dat, path, quote=FALSE, row.names=FALSE)
}

