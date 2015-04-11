# internal function called by mbpca

# xb - block b
# t - global score
# tb - block score
# pb - block loading
# method - deflation with respect to ...

#   t, tb, pb could be returned by function "nipalsSoftK" or "msvd"

# return a deflated matrix 



deflat <- function(x, t, tb, pb, method="globalScore") {
  # globalScore, blockScore, blockLoading
  switch(method,
         "globalScore" = lapply(x, function(xb) { xb - t %*% t(t) %*% xb / c(t(t) %*% t) }),
         "blockLoading" = mapply(SIMPLIFY = FALSE, function(xb, pb) {xb - xb %*% pb %*% t(pb)}, x, pb),
         "blockScore" = mapply(SIMPLIFY = FALSE, function(xb, tb) {xb - tb %*% t(tb) %*% xb / c(t(tb) %*% tb)}, x, tb))
}