setwd("~/CloudChen/Projects/MOGSA/")
library(mogsa)
library(biomaRt)
library(pheatmap)
library(NMF)
library(cluster)
library(class)
library(densityClust)
library(ConsensusClusterPlus)

# loading gene expression data and supplementary data
data <- readRDS(file = "Res/2014.10.21_BLCA_mRNA_CNV_repro/DataProc/data_matched_processed.RDS")

# =================================================================================
# =                                                                               =
# =                       run multiple omics data analysis                        =
# =                                                                               =
# =================================================================================

# multi-omics analysis
moana_mfa <- moa(data$expr, proc.row = "center_ssq1", w.data = "lambda1", statis = FALSE)
# show eigenvalues
plot(moana_mfa, value="eig", type=2, n=30)
# show normed eigenvalues, contribution are separated according to datasets
plot(moana_mfa, value="tau", type=2, n=30)

# top 5 PC capturing the unique strutures

# # another preprocessing method, STATIS, did not used
# moana_statis <- moa(data$expr, proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
# plot(moana_statis, value="eig", type=2, n=30)
# plot(moana_statis, value="tau", type=2, n=30)

# =================================================================================
# =                                                                               =
# =       basic explroe by hierachical clustering and silhouette plot             =
# =                                                                               =
# =================================================================================
plotClusterSilhouette <- function(x, d="pearson", link.method="ward.D", return=NULL) {
  
  if (d %in% c("pearson", "spearman")) {
    dst <- as.dist(1-cor(t(x), method = d))
  } else {
    dst <- dist(x, method = d)
  }
  
  hcl <- hclust(dst, method = link.method)
  layout(matrix(1:6, 2, 3, byrow = TRUE))
  for (i in 2:7 ) {
    cls <- cutree(hcl, k = i)
    sil <- silhouette(cls, dst)
    plot(sil, col="gray25", border="gray25")
    if (i == return)
      res <- cls
  }
  return(invisible(res))
}

latentVar <- moana_mfa@fac.scr[, 1:5]

# explore
plotClusterSilhouette(latentVar, d = "euclidean", link.method = "ward.D")
plotClusterSilhouette(latentVar, d = "pearson", link.method = "ward.D")
plotClusterSilhouette(latentVar, d = "spearman", link.method = "ward.D")

# select pearson and ward.D method, different number of cluster
cls2 <- plotClusterSilhouette(latentVar, d = "pearson", link.method = "ward.D", return = 2)
cls3 <- plotClusterSilhouette(latentVar, d = "pearson", link.method = "ward.D", return = 3)
cls4 <- plotClusterSilhouette(latentVar, d = "pearson", link.method = "ward.D", return = 4)
cls5 <- plotClusterSilhouette(latentVar, d = "pearson", link.method = "ward.D", return = 5)
cls6 <- plotClusterSilhouette(latentVar, d = "pearson", link.method = "ward.D", return = 6)

# =================================================================================
# =                                                                               =
# =                       quick prediction strength check                         =
# =                                                                               =
# =================================================================================
pred.strength <- function(x, cls, n=20, nnk=3) {
  if (!is.factor(cls))
    cls <- as.factor(cls)
  require(parallel)
  r <- mclapply(1:n, mc.cores = 12, function(u) {
    idx <- lapply(levels(cls), function(vv) sample(which(cls == vv)))
    idxTrain <- unlist(lapply(idx, function(xx) xx[1:ceiling(length(xx)/2)]))
    idxTest <- (1:length(cls))[! (1:length(cls)) %in% idxTrain]
    
    cvI <- knn(train = x[idxTrain, ], test = x[idxTest, ], cl = cls[idxTrain], k = nnk)
    min(tapply(cvI==cls[idxTest], cls[idxTest], function(x) sum(x)/length(x)))
  })
  return(unlist(r))
}

ps2 <- pred.strength(x = latentVar, cls = cls2, n=100, nnk=3)
ps3 <- pred.strength(x = latentVar, cls = cls3, n=100, nnk=3)
ps4 <- pred.strength(x = latentVar, cls = cls4, n=100, nnk=3)
ps5 <- pred.strength(x = latentVar, cls = cls5, n=100, nnk=3)
ps6 <- pred.strength(x = latentVar, cls = cls6, n=100, nnk=3)

boxplot(cbind(ps2, ps3, ps4, ps5, ps6))

# =================================================================================
# =                                                                               =
# =                            consnesus clustering                               =
# =                                                                               =
# =================================================================================
distance <- as.dist(1-cor(t(latentVar), method = "pearson"))
ccl <- ConsensusClusterPlus(distance, maxK = 6)
ccl_cls <- ccl[[3]]$consensusClass


# =================================================================================
# =                                                                               =
# =                            consnesys clustering                               =
# =                                                                               =
# =================================================================================
dc <- densityClust(distance, gaussian = TRUE)
plot(dc)
dcr <- findClusters(dc, rho = 6, delta = 0.5, plot = TRUE)
dc_cls <- dcr$clusters

# =================================================================================
# =                                                                               =
# =                                   comparison                                  =
# =                                                                               =
# =================================================================================

sil <- silhouette(cls3, distance)
rownames(sil) <- names(cls3)
sil <- sil[order(-sil[, "cluster"], sil[, "sil_width"], decreasing = TRUE), ]

cls_all <- cbind(cls3, ccl_cls, dc_cls)
cls_all <- cls_all[rownames(sil), ]

image(cls_all, col=2:4)









