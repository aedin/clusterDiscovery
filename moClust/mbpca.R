# multiple block PCA, main called function

# ncomp - number of components
# ... - currently nothing
# method - one of c("globalScore", "blockScore", "blockLoading")
# k - number por percentage of retained elements, could be"all"
# center - logical; if the data should be centered
# scale - logical; if the data should be scaled, passed to function "scale"
# option - how each of the dataset should be normalized, one of c("uniform", "inertia", "lambda1")
# maxiter - integer, maximum number of allow iteration in the NAPLS algorithm
# moa - logical, if the result should be in a class of "moa"
# verbose - logical; if the process (# of PC) should be printed

mbpca <- function(x, ..., ncomp, method, k="all", center=TRUE, 
                  scale=FALSE, option="uniform", maxiter=1000, 
                  moa=TRUE, verbose=TRUE, svd.solver=c("svd", "fast.svd", "propack")) {
  
  #   x <- list(x, ...)
  call <- match.call()
  x <- processOpt(x, center=center, scale=scale, option = option)
  prddata <- x
  nc <- sapply(x, ncol)
  keepAll <- k == 'all' | k > min(nc)
  
  ssl <- match.arg(svd.solver)
  svdf <- switch(ssl, 
                 "svd" = svd,
                 "fast.svd" = fast.svd,
                 "propack" = function(X) propack.svd(X, neig = 1))
  
  for (i in 1:ncomp) {
    if (verbose)
      cat(paste("calculating component ", i, " ...\n", sep = ""))
    if (keepAll)
      r <- msvd(x, svd.sol=svdf) else
        r <- nipalsSoftK(x, maxiter=maxiter, k=k)
    x <- deflat(x, r$t, r$tb, r$pb, method)
    ge <<- x
    if (i == 1)
      res <- r else {
        res$t <- cbind(res$t, r$t)
        res$w <- cbind(res$w, r$w)
        res$tb <- mapply(cbind, res$tb, r$tb, SIMPLIFY = FALSE)
        res$pb <- mapply(cbind, res$pb, r$pb, SIMPLIFY = FALSE)
      }
  }
  if (moa) 
    res <- toMoa(prddata, res, call=call)
  return(res)
}