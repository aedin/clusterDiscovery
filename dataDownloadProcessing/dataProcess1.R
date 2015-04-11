# this script selects common participants in the two datasets and
# reorder the expression/CNV matrix so the columns in two datasets
# are matched. The result is saved.


# setwd
setwd("~/CloudChen/Projects/MOGSA/")

# =================================================================================
# =                                                                               =
# =               data loading and processing, find the intersection              =
# =                                                                               =
# =================================================================================

load("Dat/BLCA/Assembly_process/cnv_snp_6_nosnp_hg19.txt.rda")
cnv <- list(data = Data, des = Des)
load("Dat/BLCA/Assembly_process/rnaseqv2.txt.rda")
rnaseqv2 <- list(data = Data, des = Des)

samples <- lapply(list(cnv, rnaseqv2), function(x) colnames(x$data))
samples <- lapply(samples, function (x) substr(x, 1, 15))

participants <- lapply(samples, function (x) substr(x, 1, 12))

common_parti <- Reduce("intersect", participants)
common_samples <- Reduce("intersect", samples)

length(common_parti)
length(common_samples)

# find normal matched tissue
cs2 <- common_samples[!grepl("11", common_samples)]

# =================================================================================
# =                                                                               =
# =               select intersect patients across the two datasets               =
# =                                                                               =
# =================================================================================

# CNV
Dcnv <- cnv$data[, samples[[1]] %in% cs2]
# remove duplicates
Dcnv <- Dcnv[, ! colnames(Dcnv) %in% c("TCGA-BL-A13J-01A-11D-A273-01", "TCGA-BL-A0C8-01A-11D-A273-01" )]
dim(Dcnv)
sort(table(samples[[1]][samples[[1]] %in% common_samples]))

# RNA seq
Drna <- rnaseqv2$data[, samples[[2]] %in% cs2]
dim(Drna)


# =================================================================================
# =                                                                               =
# =                   reorder patients in the two datasets                        =
# =                                                                               =
# =================================================================================

# reorder patients
samples_sel <- cbind(sort(colnames(Dcnv)), sort(colnames(Drna)))
samples_sel <- as.data.frame(samples_sel)

colorder <- substr(samples_sel[, 1], 1, 12)
samples_sel$parti <- colorder
colnames(samples_sel) <- c("cnv", "rna", "participant")

# =================================================================================
# =                                                                               =
# =               rearrange fo the order of the columns                           =
# =                                                                               =
# =================================================================================
colnames(Dcnv) <- substr(colnames(Dcnv), 1, 12)
colnames(Drna) <- substr(colnames(Drna), 1, 12)

Dcnv <- Dcnv[, colorder]
Drna <- Drna[, colorder]

# =================================================================================
# =                                                                               =
# =                         create output data                                    =
# =                                                                               =
# =================================================================================
# 
ori_data <- list(cnv = Dcnv, rna=Drna)
des <- list(cnv=cnv$des, rna=rnaseqv2$des)
data <- list(data=ori_data, des=des, sample=samples_sel)


saveRDS(data, file = "Res/2014.10.21_BLCA_mRNA_CNV_repro/DataProc/original_data_matched_samples.RDS")