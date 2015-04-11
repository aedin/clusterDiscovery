# internal function call by "mbpca"

# x - the input matrix, rows are observations, columns are variables
# maxiter - # of maximum interation the algorithm can run
# k - the number (>=1) or proportion (<1) of variable want to keep

# find the tb, pb, t, w for a matrix
# result comparable with msvd

nipalsSoftK <- function(x, maxiter, k) {
  # NIPALS
  # t <- x[[1]][, 1]
  t <- svd(do.call("cbind", x))$u[, 1]
  
  regproj <- function(xb, t, k) { 
    pb <- t(xb) %*% t / c(t(t) %*% t)
    pb <- pb/sqrt(sum(pb^2))
    pb <- softK(pb, k)  # soft-thresholding
    tb <- xb %*% pb
    list(tb=tb, pb=pb) # t-score, p-loading
  }
  
  for (i in 1:maxiter) {
    told <- t
    rp <- lapply(x, regproj, t, k)
    tm <- sapply(rp, "[[", "tb")
    w <- t(tm) %*% t / c(t(t) %*% t)
    w <- w/sqrt(sum(w^2))
#     w <- w/sum(w)
    t <- tm %*% w
    if (all.equal(c(t), c(told))[1] == TRUE)
      break
    if (i == maxiter)
      warning("Algrithm does not converge!")
  }
  
  res <- list(tb = lapply(rp, "[[", "tb"),
              pb = lapply(rp, "[[", "pb"),
              t = t,
              w = w)
  
  return(res)
}



