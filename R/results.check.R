
#' @title GenomeHetFilter()
GenomeHetFilter <- function(obs, mode.res, max.non.clonal, max.neg.genome,Q, verbose=FALSE) {
  ## calculate provisional seg_Z_tab and filter out modes that imply > 50% genome het.
  ## and filter out modes with > 2.5% het genome < 0
  mode.tab <- mode.res[["mode.tab"]]
  ## init both to all zeros
  frac.het <- frac.neg.het <- rep(0, nrow(mode.tab))
  
  for (i in seq_len(nrow(mode.tab))) {
    delta <- mode.tab[i, "delta"]
    comb <-  GetCopyRatioComb(Q, delta)
    seg.z <- mode.res[["seg.qz.tab"]][i, , Q+1]
    frac.het[i] <- sum(seg.z * obs$w)
    frac.neg.het[i] <- sum((obs$w * seg.z)[obs$ratio < comb[1]])
  }
  
  if (max.non.clonal > 0) {
    nc.ix <- (frac.het > max.non.clonal)
    
    if (verbose) {
      print(paste("removing ", sum(nc.ix), " / ", length(nc.ix),
                  " modes with > ", (max.non.clonal * 100), "% genome non-clonal.", sep=""))
    }
  } else {
    nc.ix <- rep(FALSE, length(frac.het))
  }
  
  if (max.neg.genome > 0) {
    neg.mode.ix <- (frac.neg.het > max.neg.genome) & (!nc.ix)
    
    if (verbose) {
      print(paste("removing ", sum(neg.mode.ix), " / ", length(neg.mode.ix),
                  " modes with >", (max.neg.genome * 100) ,
                  "% genome non-clonal < 0 copies.", sep=""))
    }
  } else {
    neg.mode.ix <- rep(FALSE, length(frac.het))
  }
  
  ## return the 'bad' indices
  return(nc.ix | neg.mode.ix)
}

#' @title ReorderModeRes()
ReorderModeRes = function(mode.res, ix, DROP=FALSE) {
  mode.res[["mode.tab"]] = mode.res[["mode.tab"]][ix,, drop=DROP]
  mode.res[["seg.z.tab"]] = mode.res[["seg.z.tab"]][ix,, drop=DROP]
  mode.res[["seg.qz.tab"]] = mode.res[["seg.qz.tab"]][ix,,, drop=DROP]
  mode.res[["seg.q.tab"]] = mode.res[["seg.q.tab"]][ix,,, drop=DROP]
  mode.res[["theta.q.tab"]] = mode.res[["theta.q.tab"]][ix, ,drop=DROP]
  mode.res[["theta.qz.hat"]] = mode.res[["theta.qz.hat"]][ix, ,drop=DROP]
  mode.res[["chr.arm.tab"]] = mode.res[["chr.arm.tab"]][ix ,,, , drop=DROP]
  mode.res[["mode.clust.p"]] = mode.res[["mode.clust.p"]][ix , , drop=DROP]
  mode.res[["mode.posts"]] = mode.res[["mode.posts"]][ix]
  return(mode.res)
}
