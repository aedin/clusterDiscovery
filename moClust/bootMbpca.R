# bootstrapping kernel
# data - the data to bootstrap
# replace - Boolean; sampling with or without replacement
# B - number of bootstrap
# mc.cores - Integer; number of cores used in bootstrap. This value is passed to function mclapply
# resample - Could be either "sample" or "total". "sample" means sample wise resampling. "total" means total resampling. See "detail"
# ncomp - passed to mbpca
# method - passed to mbpca
# k - passed to mbpca
# center - passed to mbpca
# scale - passed to mbpca
# option - passed to mbpca
# maxiter - passed to mbpca

bootMbpcaK <- function(data, replace, B=100, mc.cores=1, resample = c("sample", "total", "gene"), 
                       ncomp, method, k, 
                       center=FALSE, scale=FALSE, option="uniform", 
                       maxiter=1000, svd.solver=c("svd", "fast.svd", "propack")) {
  
  
  resampleMbpca <- function(d, ncomp, method, k, center, scale, option, 
                            maxiter, replace, resample, svd.solver) {
    rsd <- switch(resample,
                  "sample" = lapply(d, function(x) x[sample(1:nrow(x), replace = replace), ]),
                  "gene" = lapply(d, function(x) apply(x, 2, sample, replace=replace)),
                  "total" = lapply(d, function(x) t(apply(x, 1, sample, replace=replace))))
                  # "total" = lapply(d, function(x) x[sample(1:nrow(x), replace = replace), sample(1:ncol(x), replace = replace)]))
    res <- mbpca(x = rsd, verbose = FALSE, moa=FALSE, 
                 ncomp=ncomp, method=method, k=k, center=center, 
                 scale=scale, option=option, maxiter=maxiter, svd.solver)
    diag(crossprod(res$t))
  }
  
  svd.solver <- match.arg(svd.solver)
  resample <- match.arg(resample)
  r <- mclapply(1:B, mc.cores = mc.cores, function(x) 
    resampleMbpca(data, ncomp, method, k, center, scale, option, maxiter, replace, resample, svd.solver))
  do.call("rbind", r)
  
}


# bootstrapping and plot
# moa - an object of "moa" returned by mbpca
# mc.cores - Integer; number of cores used in bootstrap. This value is passed to function mclapply
# B - number of bootstrap
# replace - Boolean; sampling with or without replacement
# resample - Could be either "sample" or "total". "sample" means sample wise resampling. "total" means total resampling. See "detail"
# ncomp - passed to mbpca
# k - passed to mbpca
# method - passed to mbpca
# maxiter - passed to mbpca 

bootMbpca <- function(moa, mc.cores=1, B = 100, replace=TRUE, resample=c("sample", "gene", "total"), log="y",
                      ncomp = NULL, method = NULL, maxiter=1000, svd.solver=c("svd", "fast.svd", "propack"), plot=TRUE) {
  data <- moa@data
  call <- moa@call
  
  if (is.null(ncomp)) {
    ncomp <- call$ncomp
    if (!is.integer(ncomp))
      stop("Invalid ncomp, please set by yourself.")
    cat(paste("ncomp is set to ", ncomp, ".\n", sep = ""))
  }
  if (is.null(method)) {
    method <- call$method
    cat(paste("method is set to '", method, "'.\n", sep = ""))
  }
  if (!is.null(call$maxiter) && call$maxiter > 0 && is.integer(call$maxiter)) {
    maxiter <- call$maxiter
    cat(paste("maxiter is set to ", maxiter, ".\n", sep = ""))
  }  
  if (call$k != "all") {
    call$k <- "all"
    call$verbose <- FALSE
    moa <- eval(call)
  }
  
  ncomp <- min(c(ncomp), length(moa@eig))
  svd.solver <- match.arg(svd.solver)
  resample <- match.arg(resample)
  btk <- bootMbpcaK(data, B = B, mc.cores=mc.cores, replace = replace, resample=resample,
                    option = "uniform", center = FALSE, scale = FALSE,
                    ncomp = ncomp, k = "all", method = method, maxiter=maxiter, svd.solver=svd.solver)
  
  if (plot) {
    sc <- min(ncol(btk), length(moa@eig))
    isc <- 1:sc
    boxplot(rbind(moa@eig[isc], btk[, isc]), col=NA, border = NA, log=log)
    boxplot(btk, add = TRUE, pch=NA)
    lines(isc, moa@eig[isc], pch=20)
    points(isc, moa@eig[isc], pch=20)
  }
  
  btk
}
