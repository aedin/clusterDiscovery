# convert to MOA
# internal function called by "mbpca"
# data - processed data in mbpca
# x - the object calculated in mbpca
# call - the call of mbpca

toMoa <- function(data, x, call) {
  
  dname <- names(data)
  x$pb <- x$pb[dname]
  x$tb <- x$tb[dname]
  
  eig <- diag(crossprod(x$t))
  tab.dim <- sapply(data, dim)
  partial.eig <- sweep(x$w^2, 2, eig, "*")
  eig.vec <- sweep(x$t, 2, sqrt(eig), "/")
  
  rn <- unlist(lapply(x$pb, rownames))
  rn <- paste(rn, rep(names(data), tab.dim[2, ]), sep = "_")
  loading <- do.call("rbind", x$pb)
  rownames(loading) <- rn
  
  rownames(x$t) <- rownames(x$tb[[1]])
  fac.scr <- x$t
  partial.fs <- x$tb
  ctr.obs <- NA # sweep(abs(x$t), 2, apply(abs(x$t), 2, sum), "/")
  ctr.var <- NA # sweep(abs(loading), 2, apply(abs(loading), 2, sum), "/")
  ctr.tab <- NA # x$w^2
  RV <- pairwise.rv(data, match = "row")
  w.row <- NA
  
  res <- new("moa", 
             eig = eig,
             tau = NA, 
             partial.eig = partial.eig,
             eig.vec = eig.vec,
             loading = loading,
             fac.scr = fac.scr,
             partial.fs = partial.fs,
             ctr.obs = NA,
             ctr.var = NA,
             ctr.tab = NA,
             RV= RV,
             w.row = NA,
             data = data,
             tab.dim =tab.dim,
             call=call)
  return(res)
}

