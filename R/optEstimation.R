
#' @title FitSample()
#' @export
FitSample = function(seg.dat, Q, pi_theta_qz, sigma.h, tau.dom,sigma.h.dom, verbose=FALSE) {
  require(numDeriv)
  kThetaQz = rep(1 / (Q + 1), Q + 1)
  obs = seg.dat
  res = FindLocationModes(obs, Q, kThetaQz, sigma.h,tau.dom, verbose=verbose)
  if (!is.null(res[["mode.flag"]])) { 
    return(list("mode.flag"=res[["mode.flag"]]))
  } else {
    mode.tab = res[["mode.tab"]]
  }
  if (all(is.na(mode.tab)))  {
    return(list(mode.flag="ERROR"))
  }
  n.modes = nrow(mode.tab)
  #theta.qz.hat = matrix(NA, nrow=n.modes, ncol=Q+1)
  theta.q.tab = array(NA, dim=c(n.modes, Q))
  if (verbose) {
    print(paste("Optimizing LL(data, theta.qz, sigma.h | comb) for ",
                n.modes, " modes: ", sep=""))
  }
  
  seg.z.tab = array(NA, dim=c(n.modes, length(obs$w)))
  seg.qz.tab = array(NA, dim=c(n.modes, length(obs$w), Q+1))
  seg.q.tab = array(NA, dim=c(n.modes, length(obs$w), Q))
  
  #chr.arm.tab = array(NA, dim=c(n.modes, 1, nrow(chr.arms.dat), Q))
  #dimnames(chr.arm.tab)[[3]] = rownames(chr.arms.dat)
  mode.tab = cbind(mode.tab, NA, NA, NA, NA, NA, NA, NA, NA)
  colnames(mode.tab)[c((ncol(mode.tab)-7) : ncol(mode.tab))] = c("genome mass", "sigma.h.hat", "theta.z.hat", "frac.het", "log.partition", "entropy", "clust_LL", "post_LL")
  
  
  for (i in 1:n.modes) {
    delta = mode.tab[i, "delta"]
    ## can be NA for non-array platforms
    comb =  GetCopyRatioComb(Q, delta)
    ## optimize
    res = OptThetaQzSigmaH(obs, comb, sigma.h.dom, pi_theta_qz, verbose=verbose)
    
    mode.tab[i, "sigma.h.hat"] = res[["sigma.h.hat"]]
    mode.tab[i, "log.partition"] = res[["LL"]]
    ## initial version with data-fit term only
    mode.tab[i, "post_LL"] = res[["LL"]]  
    lambda.qz.res = res[["lambda.qz.res"]]
    mode.tab[i, "theta.z.hat"] = res$theta_qz_hat[(Q + 1)]
    
    ## posterior distribution over segment CNs
    res = GetSegQzPost(obs$ratio, obs$sd,obs$w, comb, NA, mode.tab[i,"sigma.h.hat"],theta.qz=res[["theta_qz_hat"]])
    seg.z.tab[i, ] = res[["QZ"]][, (Q+1)]
    seg.qz.tab[i, , ] = res[["QZ"]]
    seg.q.tab[i, , ] = res[["Q"]]
    theta.q.tab[i, ] = colSums(res[["Q"]] * obs$w)
    
    mode.tab[i, "theta.z.hat"] = sum(seg.qz.tab[i, , (Q + 1)] * obs$w)
    
    ## compute % non-clonal genome
    mode.tab[i, "frac.het"] = sum(obs$w * seg.z.tab[i,])  
    
    mode.tab[i,"genome mass"] = 1 * sum(c((1:Q)-1) * colSums(res[["Q"]] * obs$w))
    
    
    ## weighted entropy average over segs
    mode.tab[i, "entropy"] = CalcFitEntropy(obs, res[["QZ"]])
    
    #chr.arm.tab[i, , , ] = CalcChrArmDistr(seg.obj, res[["Q"]], chr.arms.dat)
  }
  
  return(list(mode.tab=mode.tab, mode.posts=NA,theta.q.tab=theta.q.tab,
              seg.z.tab=seg.z.tab, seg.qz.tab=seg.qz.tab,
              seg.q.tab=seg.q.tab,mode.flag=NA))
}


#' @title OptThetaQzSigmaH()
#' @export

OptThetaQzSigmaH <- function(obs, comb, sigma.h.dom, pi_theta_qz, verbose=FALSE) {
  objective <- function(par, obs, comb, lambda) {
    sigma.h <- par
    max.q <- length(comb)
    theta.qz.hat <- GetThetaQzPost(obs$ratio, obs$sd,obs$w, comb, sigma.h, pi_theta_qz)
    theta.qz.map <- theta.qz.hat
    theta.qz.map <- theta.qz.map / sum(theta.qz.map)
    theta.qz.hat <- theta.qz.map
    if (verbose) {
      cat("S")
    }
    LL <- CalcNormLoglik(obs$ratio, obs$sd, obs$w, comb, NA, theta.qz.hat, sigma.h)
    if (is.nan(LL)) {
      LL <- -Inf
    }
    return(LL)
  }
  
  max_q = length(comb)
  ## find a value of Lambda to initialize with
  use.sigma.h <- mean(sigma.h.dom)
  theta.qz.hat <- GetThetaQzPost(obs$ratio, obs$sd,obs$w, comb, use.sigma.h, pi_theta_qz)
  lambda.qz.res <- GetLambda(theta.qz.hat, obs$w, max.q + 1,verbose=verbose)
  init.lambda <- lambda.qz.res[["Lambda"]]
  
  ## now optimize sigma.h
  res <- optimize(f = objective, interval = sigma.h.dom, tol = 0.001, maximum = TRUE, 
                  obs = obs, comb = comb, lambda = init.lambda)
  
  sigma.h.hat <- res[["maximum"]]
  
  ## now get LL, Lambda, and theta.qz.hat conditional on optimized sigma.h
  theta.qz.hat <- GetThetaQzPost(obs$ratio, obs$sd,obs$w, comb, sigma.h.hat, pi_theta_qz)
  
  lambda.qz.res <- GetLambda(theta.qz.hat, obs$w, max.q + 1, verbose=verbose)
  LL <- CalcNormLoglik(obs$ratio, obs$sd, obs$w,comb, lambda.qz.res, NA, sigma.h.hat)
  
  if (is.nan(LL)) {
    LL <- -Inf
    if (verbose) {
      print("Warning: NaN loglik in opt_theta_Z_sigma_H")
    }
  }
  
  if (verbose) {
    cat(paste(" sigma.h.hat=", round(sigma.h.hat, 5), sep = ""))
    cat("\n")
  }
  
  return(list(LL = LL, sigma.h.hat = sigma.h.hat, theta_qz_hat = theta.qz.hat, 
              lambda.qz.res = lambda.qz.res))
}


#' @title GetThetaQzPost()
#' @export

GetThetaQzPost <- function(d, sigma, w, comb, sigma.h, pi_theta_qz, theta.qz = NA, comb.s = 1) {
  q <- length(comb)
  unif.qz <- rep(1 / (q + 1), q + 1)
  log.ss.p.1 <- sapply(comb, dnorm, d,
                       sqrt(sigma^2 + comb.s^2 * sigma.h^2), log = TRUE)
  
  if (is.na(theta.qz)) {
    log.pqz <- matrix(log(unif.qz), nrow = length(w), ncol = q + 1, byrow = TRUE)
  } else {
    log.pqz <- matrix(log(theta.qz), nrow = length(w), ncol = q + 1, byrow = TRUE)
  }
  
  theta.q <- unif.qz[c(1:q)]
  theta.q <- theta.q / sum(theta.q)
  ##    log_pQ <- matrix(log(theta.q), nrow = length(w), ncol = q, byrow = TRUE)
  
  log.comb.seg.mat <- cbind(log.pqz[, c(1:q)] + log.ss.p.1, log(1/7) +
                              log.pqz[, q + 1])  
  seg_qz <- exp(log.comb.seg.mat - LogAdd(log.comb.seg.mat))    
  
  prior_w = pi_theta_qz$PC
  pi_w = pi_theta_qz$W
  pi_w = pi_theta_qz$W / sum(pi_w)
  
  post_theta_qz = (1 - prior_w) * colSums(w * seg_qz) + (prior_w * pi_w)
  
  return(post_theta_qz)
}


#' @title GetSegQzPost()
#' @export

GetSegQzPost <- function(d, sigma, w, comb, lambda.qz.res, sigma.h, 
                         theta.qz=NA, comb.s=1) {
  if (is.list(lambda.qz.res)) {
    log.pqz <- t((LambdaExpr(lambda.qz.res[["Lambda"]], w))) +
      log(lambda.qz.res[["norm.term"]])
    
    ## length of lambda is N_states - 1 ;  == Q
    q <- length(lambda.qz.res[["Lambda"]])
    comb <- comb[1:q]
    log.pq <- log.pqz[, c(1:q)]
    log.pq <- log.pq - LogAdd(log.pq)  ## renormalize
  } else {
    q <- length(comb)   
    theta.q <- theta.qz[c(1:q)]
    theta.q <- theta.q / sum(theta.q)
    log.pqz <- matrix(log(theta.qz), nrow = length(w), ncol = q + 1, byrow = TRUE)
    log.pq <- matrix(log(theta.q), nrow = length(w), ncol = q, byrow = TRUE)
  }
  
  log.ss.p.1 <- sapply(comb, dnorm, d,
                       sqrt(sigma^2 + comb.s^2 * sigma.h^2), log = TRUE)
  
  log.comb.seg.mat <- cbind(log.pqz[, c(1:q)] + log.ss.p.1,
                            log(1/7) + log.pqz[, q + 1])
  
  seg_qz <- exp(log.comb.seg.mat - LogAdd(log.comb.seg.mat))
  
  log.comb.seg.mat <- log.pq + log.ss.p.1  ## Nseg x Q
  seg.Q <- exp(log.comb.seg.mat - LogAdd(log.comb.seg.mat))
  
  return(list(QZ = seg_qz, Q = seg.Q))
}


#' @title CalcFitEntropy()
#' @export
CalcFitEntropy <- function(obs, seg.qz) {
  H <- function(p) {
    v <- -p * log(p)
    v[p == 0] <- 0
    return(sum(v))
  }
  
  seg.h <- apply(seg.qz, 1, H)
  res <- sum(obs$w * seg.h )
  
  return(res)
}


#' @title CalcChrArmDistr()
#' @export
CalcChrArmDistr = function(seg.obj, seg_q, chr_arms_dat) {
  n_arm = nrow(chr_arms_dat)
  chr_arm_tab = array(NA, dim=c(1, n_arm, ncol(seg_q)))
  
  for (i in seq_len(n_arm)) {
    chr_dat = get_tcr_chr_arm_segs(seg.obj, chr_arms_dat[i, ])      
    
    if (length(chr_dat$int_W) == 0 ) { 
      next 
    }
    
    chr_arm = array(NA, dim=c(1, ncol(seg_q)))
    chr_arm[1, ] = colSums(seg_q[chr_dat$ix, , drop=FALSE] * chr_dat$int_W)
    
    chr_arm_tab[, i, ] = 0
    chr_arm_tab[1, i, which.max(chr_arm[1,])] = 1
  }
  
  return(chr_arm_tab)
}


#' @title doCNV()
#' @export
doCNV <- function(seg.dat){
  Q = 8 # maximum integer copy number
  CNVres=list()
  for (Q in 3:8){
    kTauDom = c(0.95, Q) #range of ploidy
    kSigmaHDom = c(0, sd(seg.dat$ratio)) #range of variation at sample level
    #kPiSomThetaQ = c(100, 50, rep(2, (Q - 2)))
    W = c(1, 25, 100, 25, 10, 5)
    if (Q < 6){
      pi_theta_qz = list(W = c(W[1:Q], 10), PC = 0.05)
    }else{
      pi_theta_qz = list(W = c(W, rep(1, Q - 6), 10), PC = 0.05)
    }
    sigma.p = 0
    sigma.h = sigma.p ##define  by user
    mode.res = FitSample(seg.dat, Q, pi_theta_qz, sigma.h, kTauDom, kSigmaHDom, verbose=FALSE)
    
    if (inherits(mode.res, "try-error")) {
      mode.res = list(mode.flag = "FAIL")
    }
    if (is.na(mode.res[["mode.flag"]])) {
      bad.ix = GenomeHetFilter(seg.dat, mode.res, max.non.clonal=0,max.neg.genome=0, Q, verbose=FALSE)
      if (sum(bad.ix) == nrow(mode.res[["mode.tab"]])) {
        mode.res = list(mode.flag="TAU_DOM")
      } else {
        mode.res = ReorderModeRes(mode.res, !bad.ix)
      }
    }
    mode.tab=mode.res$mode.tab
    CNVindex=order(-mode.tab[,8],exp(mode.tab[,12]),decreasing = T)[1]
    #prob=exp(-mode.tab[,1]*mode.tab[,5]/sum(mode.tab[,5])*mode.tab[,7])
    #probOrder=order(prob,decreasing = T)
    #CNVindex=probOrder[which.min(mode.tab[probOrder[1:2],8])]
    tau=mode.tab[CNVindex,1]
    delta=mode.tab[CNVindex,2]
    frac.het=mode.tab[CNVindex,8]
    post_LL=mode.tab[CNVindex,12]
    thetaCNV=mode.res$theta.q.tab[CNVindex,]
    qz=mode.res$seg.qz.tab[CNVindex,,]
    integerCNV=apply(qz, 1, which.max)-1
    CNVres[[Q-2]]=list(tau=tau,delta=delta,frac.het=frac.het,post_LL=post_LL,thetaCNV=thetaCNV,qz=qz,integerCNV=integerCNV)
  }
  res=c()
  for (i in 1:length(CNVres)){
    res=rbind(res,c(CNVres[[i]]$tau,CNVres[[i]]$delta,CNVres[[i]]$frac.het,CNVres[[i]]$post_LL))
  }
  CNVindex=order(-res[,3],exp(res[,4]),decreasing = T)[1]
  CNVestimate=CNVres[[CNVindex]]
  seg.dat$CNV=CNVestimate$delta*CNVestimate$integerCNV
  seg.dat$integerCNV=CNVestimate$integerCNV
  return(list(CNVestimate=CNVestimate,seg.dat=seg.dat))
}


#' @title EstimateCNVPost()
#' @export
EstimateCNVPost <- function(seg.dat, Q,delta,sigma.h){
  d=seg.dat$ratio
  sigma=seg.dat$sd
  w=seg.dat$w
  comb.s=1
  comb=GetCopyRatioComb(Q, delta)
  log.ss.p.1 <- sapply(comb, dnorm, d,sqrt(sigma^2 + comb.s^2 * sigma.h^2), log = TRUE)
  seg.Q <- apply(log.ss.p.1,1,which.max)-1
  return(seg.Q)
}

#' @title EstimateCNVPost1()
#' @export
##return different dat
EstimateCNVPost1 <- function(seg.dat, CNVestimate){
  d=seg.dat$ratio
  sigma=seg.dat$sd
  w=seg.dat$w
  Q=CNVestimate$Q
  delta=CNVestimate$delt
  sigma.h=CNVestimate$sigma.h
  comb.s=1
  comb=GetCopyRatioComb(Q, delta)
  log.ss.p.1 <- sapply(comb, dnorm, d,sqrt(sigma^2 + comb.s^2 * sigma.h^2), log = TRUE)
  seg.Q <- apply(log.ss.p.1,1,which.max)-1
  seg.dat$integerCNV<-seg.Q
  return(seg.dat)
}



#' @title doCNV1()
#' @export

doCNV1 <- function(seg.dat){
  Q = 8 # maximum integer copy number
  CNVres=list()
  kTauDom = c(0.95, Q) #range of ploidy
  kSigmaHDom = c(0, sd(seg.dat$ratio)) #range of variation at sample level
  W = c(1, 25, 100, 25, 10, 5)
  if (Q < 6){
    pi_theta_qz = list(W = c(W[1:Q], 10), PC = 0.05)
  }else{
    pi_theta_qz = list(W = c(W, rep(1, Q - 6), 10), PC = 0.05)
  }
  sigma.p = 0
  sigma.h = sigma.p ##define  by user
  mode.res = FitSample(seg.dat, Q, pi_theta_qz, sigma.h, kTauDom, kSigmaHDom, verbose=FALSE)
  
  if (inherits(mode.res, "try-error")) {
    mode.res = list(mode.flag = "FAIL")
  }
  if (is.na(mode.res[["mode.flag"]])) {
    bad.ix = GenomeHetFilter(seg.dat, mode.res, max.non.clonal=0,max.neg.genome=0, Q, verbose=FALSE)
    if (sum(bad.ix) == nrow(mode.res[["mode.tab"]])) {
      mode.res = list(mode.flag="TAU_DOM")
    } else {
      mode.res = ReorderModeRes(mode.res, !bad.ix)
    }
  }
  mode.tab=mode.res$mode.tab
  subclonal=apply(mode.res$seg.z.tab,1,sum)
  mode.tab=cbind(mode.tab,subclonal)
  diff=sapply(1:dim(mode.tab)[1], function(i,mode.tab,seg.dat){
    tau=mode.tab[i,1]
    delta=mode.tab[i,2]
    sigma.h=mode.tab[i,6]
    integerCNV=EstimateCNVPost(seg.dat, Q,delta,sigma.h)
    expectCNV=delta*integerCNV
    diff=lossCNV(integerCNV,tau)
    return(diff)
  },mode.tab,seg.dat)
  mode.tab=cbind(mode.tab,diff)
  mode.tab=as.data.frame(mode.tab)
  
  #
  mode.tab$score=(mode.tab$subclonal+mode.tab$diff)/2
  CNVindex=which.min(mode.tab$score)
  tau=mode.tab$tau[CNVindex]
  delta=mode.tab$delta[CNVindex]
  frac.het=mode.tab$theta.z.hat[CNVindex]
  sigma.h=mode.tab$sigma.h.hat[CNVindex]
  subclonal.het=mode.tab$subclonal[CNVindex]
  thetaCNV=mode.res$theta.q.tab[CNVindex,]
  qz=mode.res$seg.qz.tab[CNVindex,,]
  integerCNV=apply(qz, 1, which.max)-1
  heterogeneity=mode.tab$diff[CNVindex]
  CNVestimate=list(tau=tau,delta=delta,frac.het=frac.het,sigma.h=sigma.h,subclonal.het=subclonal.het,thetaCNV=thetaCNV,qz=qz,integerCNV=integerCNV,Q=Q,heterogeneity=heterogeneity)
  CNVestimate$integerCNV=EstimateCNVPost(seg.dat,Q,delta,sigma.h)
  seg.dat$CNV=CNVestimate$delta*CNVestimate$integerCNV
  seg.dat$integerCNV=CNVestimate$integerCNV
  return(list(CNVestimate=CNVestimate,seg.dat=seg.dat))
}


#' @title doCNV1_v2()
#' @export
#change the method of find the best delta (distance of copy number change)
doCNV1_v2 <- function(seg.dat,Determ_method = "Post",Q = 8){
 # Q = 8 # maximum integer copy number
  CNVres=list()
  kTauDom = c(0.95, Q) #range of ploidy
  kSigmaHDom = c(0, sd(seg.dat$ratio)) #range of variation at sample level
  W = c(1, 25, 100, 25, 10, 5)
  if (Q < 6){
    pi_theta_qz = list(W = c(W[1:Q], 10), PC = 0.05)
  }else{
    pi_theta_qz = list(W = c(W, rep(1, Q - 6), 10), PC = 0.05)
  }
  sigma.p = 0
  sigma.h = sigma.p ##define  by user
  mode.res = FitSample(seg.dat, Q, pi_theta_qz, sigma.h, kTauDom, kSigmaHDom, verbose=FALSE)
  
  if (inherits(mode.res, "try-error")) {
    mode.res = list(mode.flag = "FAIL")
  }
  if (is.na(mode.res[["mode.flag"]])) {
    bad.ix = GenomeHetFilter(seg.dat, mode.res, max.non.clonal=0,max.neg.genome=0, Q, verbose=FALSE)
    if (sum(bad.ix) == nrow(mode.res[["mode.tab"]])) {
      mode.res = list(mode.flag="TAU_DOM")
    } else {
      mode.res = ReorderModeRes(mode.res, !bad.ix)
    }
  }
  mode.tab=as.data.frame(mode.res$mode.tab)  #！！！！DEBUG
  
  # ----------------------------------------------------------------------------
  # -- Determine Copy Number 
  # ----------------------------------------------------------------------------
  #(Posterior probability method)
 # if(Determ_method=="Post"){
    subclonal=apply(mode.res$seg.z.tab,1,sum)
    mode.tab=cbind(mode.tab,subclonal)
    diff=sapply(1:dim(mode.tab)[1], function(i,mode.tab,seg.dat){
      tau=mode.tab[i,1]
      delta=mode.tab[i,2]
      sigma.h=mode.tab[i,6]
      integerCNV=EstimateCNVPost(seg.dat, Q,delta,sigma.h)
      expectCNV=delta*integerCNV
      diff=lossCNV(integerCNV,tau)
      return(diff)
    },mode.tab,seg.dat)
    mode.tab=cbind(mode.tab,diff)
    mode.tab=as.data.frame(mode.tab)
    mode.tab$score=(mode.tab$subclonal+mode.tab$diff)/2
    #CNVindex=which.min(mode.tab$score)
    #CNVindex=which.min(mode.tab$diff)
 # }else{
    # #(SSE Method)
    # outerRaw         = seg.dat$ratio %o% mode.tab$delta
    # outerRound       = round(outerRaw)
    # outerDiff        = (outerRaw - outerRound) ^ 2
    # outerColsums = colSums(outerDiff, na.rm = FALSE, dims = 1)
    # CNmult      = mode.tab$delta[order(outerColsums)]
    # CNerror     = round(sort(outerColsums), digits=2)
    # 
    # CNVindex=which.min(outerColsums)
    # 
    raw <- matrix(seg.dat$ratio,nrow(seg.dat),nrow(mode.tab))
    outerRaw <-matrix(NA,nrow(seg.dat),nrow(mode.tab))

    for(k in 1:nrow(mode.tab)){
      tau=mode.tab$tau[k]
      delta=mode.tab$delta[k]
      frac.het=mode.tab$theta.z.hat[k]
      sigma.h=mode.tab$sigma.h.hat[k]
      integerCNV=EstimateCNVPost(seg.dat,Q,delta,sigma.h)
      outerRaw[,k]=delta*integerCNV
    }
    outerDiff = (outerRaw - raw) ^ 2
    outerColsums = colSums(outerDiff, na.rm = FALSE, dims = 1)
    mode.tab$Colsums <-outerColsums
    mode.tab$score2=(mode.tab$subclonal+mode.tab$diff+mode.tab$Colsums)/3
    #CNVindex=which.min(mode.tab$score2)
    CNVindex=which.min(mode.tab$Colsums)
  #}

  tau=mode.tab$tau[CNVindex]
  delta=mode.tab$delta[CNVindex]
  frac.het=mode.tab$theta.z.hat[CNVindex]
  sigma.h=mode.tab$sigma.h.hat[CNVindex]
  subclonal.het=mode.tab$subclonal[CNVindex]
  thetaCNV=mode.res$theta.q.tab[CNVindex,]
  qz=mode.res$seg.qz.tab[CNVindex,,]   #对应CN state 的后验概率
  integerCNV=apply(qz, 1, which.max)-1  ###
  heterogeneity=mode.tab$diff[CNVindex]
  CNVestimate=list(tau=tau,delta=delta,frac.het=frac.het,sigma.h=sigma.h,subclonal.het=subclonal.het,thetaCNV=thetaCNV,qz=qz,integerCNV=integerCNV,Q=Q,heterogeneity=heterogeneity)
  CNVestimate$integerCNV=EstimateCNVPost(seg.dat,Q,delta,sigma.h)
  seg.dat$CNV=CNVestimate$delta*CNVestimate$integerCNV
  seg.dat$integerCNV=CNVestimate$integerCNV
  return(list(CNVestimate=CNVestimate,seg.dat=seg.dat))
}


#' @title lossCNV()
#' @export
lossCNV <- function(observed,tau){
  lower=0
  upper=max(observed)+1
  da=density(observed,from=lower,to=upper)
  diff=sapply(1:100,function(i,observed,lower,upper,da){
    db=density(sample(c(1:upper),length(observed),replace = T),from=lower,to=upper)
    d=data.frame(x=da$x,a=da$y,b=db$y)
    diff=sum((d$b[d$b>d$a&d$x<tau][-1]-d$a[d$b>d$a&d$x<tau][-1])*diff(d$x[d$b>d$a&d$x<tau]))
    return(diff)
  },observed,lower,upper,da)
  return(mean(diff))
}


#' @title iterEst()
#' @export

iterEst <- function(seg.dat){
  output = lapply(1:5, function(i,seg.dat){
    output<-doCNV1(seg.dat)
    #output<-doCNV1_v2(seg.dat)
    return(output)
  },seg.dat)
  tau=unlist(lapply(1:length(output), function(i,output){
    return(output[[i]]$CNVestimate$tau)
  },output))
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  index=which(tau==getmode(tau))[1]
  return(output[[index]])
}