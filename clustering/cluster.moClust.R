setwd("~/CloudChen/Projects/MOGSA/")
library(corpcor)
library(mogsa)
source("~/CloudChen/Projects/IntegrativeClustering/R/pack/deflat.R") # the directory is changed
source("~/CloudChen/Projects/IntegrativeClustering/R/pack/mbpca.R")
source("~/CloudChen/Projects/IntegrativeClustering/R/pack/nipalsSoftK.R")
source("~/CloudChen/Projects/IntegrativeClustering/R/pack/softK.R")
source("~/CloudChen/Projects/IntegrativeClustering/R/pack/msvd.R")
source("~/CloudChen/Projects/IntegrativeClustering/R/pack/processOpt.R")
source("~/CloudChen/Projects/IntegrativeClustering/R/pack/toMoa.R")
source("~/CloudChen/Projects/IntegrativeClustering/R/pack/bootMbpca.R")

# loading gene expression data and supplementary data
data <- readRDS(file = "Res/2014.10.21_BLCA_mRNA_CNV_repro/DataProc/data_matched_processed.RDS")


expr <- sapply(data$expr, t)

# no sparsity, for comparison with permutation later
res <- mbpca(x = expr, ncomp = 20, k = "all", method = "globalScore", option = "inertia", 
             center=TRUE, scale=FALSE, moa = TRUE, svd.solver = "fast")
plot(res, value="eig", type=2)


# x - a list of expression matrices or data.frame, the rows are samples columns are variables
# ncomp - the number of latent variable want to calculate
# k - the parameter controls the sparsity, "all" means no sparsity, it could be number between (0, 1), indicate the 
#   proportion of non-zero coefficients or integers, indicates the absolute number of non-zero coefficients


# use scree plot, we can select 4 latent variables

# # or you can try permutation evaluate the axes for coherent structure
# perm.sample <- bootMbpca(moa = res, ncomp = 8, mc.cores = 2, B = 30, resample = "sample", svd.solver="fast", replace = FALSE)

# introduce sparsity on variable loading matrix
res.sparse <- mbpca(x = expr, ncomp = 4, k = 0.1, method = "globalScore", option = "inertia", 
                    center=TRUE, scale=FALSE, moa = TRUE, svd.solver = "fast")

layout(matrix(1:2, 1, 2))
plot(res, value="eig", type=2)
plot(res.sparse, value="eig", type=2)
latentVar <- res.sparse@fac.scr

# the latentVar could be used as in cluster.mogsa

# =====
layout(matrix(1:4, 2, 2))
for (i in 1:4) {
  plot(res@fac.scr[, i], res.sparse@fac.scr[, i])
}










