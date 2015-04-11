# this script is continue from script "dataProcess1"
# Here we filtering out rows that contain low variance ...

setwd("~/CloudChen/Projects/MOGSA/")

library(Biobase)
library(preprocessCore)
library(genefilter)

data <- readRDS(file = "Res/2014.10.21_BLCA_mRNA_CNV_repro/DataProc/original_data_matched_samples.RDS")

# =================================================================================
# =                                                                               =
# =                             confirm                                           =
# =                                                                               =
# =================================================================================
dim(data$data$cnv)
dim(data$data$rna)
identical(colnames(data$data$cnv), colnames(data$data$rna))

# =================================================================================
# =                                                                               =
# =                    filtering copy number variation                            =
# =                                                                               =
# =================================================================================


# check distributions
hist(data$data$cnv, breaks = 100)
# slow:
# boxplot(data$data$cnv, cex=0.3, pch=20, axes=FALSE)

# calculate the MAD for each gene
dcnv <- data$data$cnv
cnv_rsd <- apply(dcnv, 1, mad, na.rm=TRUE)
cnv_rsd[is.na(cnv_rsd)] <- 0
hist(cnv_rsd, breaks = 100)

# select the upper 50%
idx <- cnv_rsd > median(cnv_rsd)
sum(idx)
dcnv <- dcnv[idx, ]
dcnv_des <- data$des$cnv[idx, ]

# write.table(dcnv_des, file = "temp.txt", row.names = FALSE, quote=FALSE, sep="\t")
# boxplot(dcnv)

# check abnormal values
sum(is.na(dcnv))
sum(is.infinite(dcnv))
hist(dcnv, breaks = 100)

dcnv[is.na(dcnv)] <- 0

# check dimensions
dim(dcnv_des)
dim(dcnv)

# =================================================================================
# =                                                                               =
# =                     filtering RNA seq data                                    =
# =                                                                               =
# =================================================================================

drna <- log10(data$data$rna+1)
hist(drna, breaks = 100)
dim(drna)

rsum <- rowSums(drna)
rsd <- rowSds(drna)
rmad <- apply(drna, 1, mad, na.rm=TRUE)

hist(rsum)
hist(rsd)
hist(rmad, breaks = 100)

idx <- rsum > 300 & rmad > 0.1

sum(idx)

drna_trim <- drna[idx, ]
dim(drna_trim)
drna_trim_des <- data$des$rna[idx, ]

hist(drna_trim, breaks = 100)

# centering patients
drna_trim <- sweep(drna_trim, 2, colMeans(drna_trim), "-")
drna_trim <- normalize.quantiles(drna_trim)
# boxplot(drna_trim)
colnames(drna_trim) <- colnames(drna)


# check abnormal data
sum(is.na(drna_trim))
sum(is.infinite(drna_trim))

# 
ttd <- as.data.frame(data$des$rna)
ttd$sd <- rsd
ttd$sum <- rsum

rownames(drna_trim) <- paste(drna_trim_des[, 1], 1:nrow(drna_trim_des), sep = "_")
rownames(dcnv) <- paste(dcnv_des[, 1], 1:nrow(dcnv_des), sep = "_")



# =================================================================================
# =                                                                               =
# =                            save data                                          =
# =                                                                               =
# =================================================================================
d <- list(expr=list(cnv=dcnv, rna=drna_trim),
          des=list(cnv=dcnv_des, rna=drna_trim_des))

saveRDS(d, file = "Res/2014.10.21_BLCA_mRNA_CNV_repro/DataProc/data_matched_processed.RDS")



# =================================================================================
# =                                                                               =
# =                            venn diagram                                       =
# =                                                                               =
# =================================================================================

library(VennDiagram)
layout(matrix(1:2, 1, 2))
hist(dcnv, breaks = 100)
hist(drna_trim, breaks = 100)

vd <- venn.diagram(list(CNV=dcnv_des[, 1], mRNA=drna_trim_des[, 1]), filename = NULL)
grid.draw(vd)




