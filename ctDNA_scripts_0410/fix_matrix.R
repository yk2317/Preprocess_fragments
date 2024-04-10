setwd("/n/data1/hms/dbmi/park/ctDNA_loci_project/fragmentomics/MSK_hg38/duplex/matrix/norm_variant_sep/")

#----
input <- "./original/out5p_norm_matrix.csv"
output <- "out5p_norm_matrix.csv"

input <- "./original/in5p_norm_matrix.csv"
output <- "in5p_norm_matrix.csv"

input <- "./original/length_norm_matrix.csv"
output <- "length_norm_matrix.csv"

#---- 
input <- "./original/out5p_matrix.csv"
output <- "out5p_matrix.csv"

input <- "./original/in5p_matrix.csv"
output <- "in5p_matrix.csv"

input <- "./original/length_matrix.csv"
output <- "length_matrix.csv"

dat <- read.csv(input)
rownames(dat) <- dat$X
dat <- dat[,-1]

print("data")
few <- colnames(dat)[which(colSums(dat) <500)]
print(few)
colSums(dat)[which(colSums(dat) <500)]

dat <- dat[,!colnames(dat) %in% few]
X <- rownames(dat)
dat <- cbind(X, dat)

write.csv(dat, output, row.names=FALSE)

