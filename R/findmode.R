
#' @title GetLambda()
#' @export

GetLambda <- function(theta, W, n.states, lambda.init=NA, w.thresh=0,verbose=FALSE) {
  if (all(theta == theta[1])) {
    lambda.res <- list()
    lambda.res[["Lambda"]] <- rep(0, (n.states-1))
    lambda.res[["norm.term"]] <- rep(1 / n.states, length(W))
    lambda.res[["resid"]] <- 0 
    return(lambda.res)
  }
  
  if (!is.na(lambda.init)) {
    lambda.res <- CalcLambda(theta, W, n.states, first.lambda.guess=lambda.init,
                             verbose=verbose)
    resid <- lambda.res[["resid"]]
    
    if (resid < 0.05) {
      return(lambda.res)
    }
  }
  wsix <- sort(W, index.return=TRUE)[["ix"]]
  orig.w <- W
  if (w.thresh > 0) {
    n.small <- max(which(cumsum(W[wsix]) < w.thresh))
    small.ix <- wsix[c(1:n.small)]
    
    W <- W[-small.ix]
    W <- W / sum(W)
  }
  
  theta[theta < .Machine$double.xmin] = .Machine$double.xmin
  
  resid <- Inf #the minimum value of loss function
  wix <- sort(W, decreasing=TRUE, index.return=TRUE)[["ix"]]
  iix <- 1
  init.ix <- wix[iix]
  
  while (resid > 0.05 & iix <= length(W)) {
    lambda.init <- -1/W[init.ix] * log(theta[c(2:length(theta))] /theta[1]) / c(1:(length(theta)-1))
    lambda.init[lambda.init == Inf] = .Machine$double.xmax / 2
    lambda.init[lambda.init == -Inf] = -.Machine$double.xmax / 2
    
    lambda.res <- CalcLambda(theta, W, n.states, first.lambda.guess=lambda.init,
                             verbose=verbose)
    resid <- lambda.res[["resid"]]
    iix <- iix+1
    init.ix <- wix[iix]
  }
  
  #if (resid > 0.05) {
  #  if (verbose) {
  #    print(paste("Warning: Lambda resid = ", resid, sep=""))
  #  }
  #  stop()
  #}
  
  if (w.thresh > 0) {
    norm.term <- rep(NA, length(orig.w))
    norm.term[small.ix] <- 1/n.states
    norm.term[-small.ix] <- lambda.res[["norm.term"]]
    lambda.res[["norm.term"]] <- norm.term
  }
  return(lambda.res)
}


#' @title CalcLambda()
#' @export
##nlm estimation to get lambda parameter and likelihood
CalcLambda <- function(theta, W, n.states, first.lambda.guess,verbose=FALSE) {
  if (all(!is.na(first.lambda.guess))) {
    cur.lambda <- first.lambda.guess
  } else {
    ## init guess
    cur.lambda <-  runif(c(n.states-1), -1e3, 1e3)    
  }
  
  while (!is.finite(L2Loss(cur.lambda, n.states, theta, W))) {
    cur.lambda <- GetInitGues(W, n.states)
    if (verbose) {
      cat("!")
    }
  }
  
  stol <- 1e-6
  gtol <- 1e-6
  
  res  <- try(nlm(f=L2Loss, p=cur.lambda, n.states=n.states, theta=theta,
                  W=W, steptol=stol, gradtol=gtol), silent=TRUE)
  if (inherits(res, "try-error")) {
    if (verbose) {
      cat("!")
    }
    return(list(resid=Inf))
  }
  
  lambda <- res[["estimate"]]
  resid=res[["minimum"]]
  ## calc normalization factor (for each seg)
  log.p.q.terms <- LambdaExpr(lambda, W)
  p.q.terms <- exp(log.p.q.terms)
  seg.totals <- apply(p.q.terms, 2, sum)
  lambda.norm.term <- 1 / seg.totals
  return(list(Lambda=lambda, norm.term=lambda.norm.term, resid=resid))
}


#' @title L2Loss()
#' @export
###define L2 loss function
L2Loss <- function(lambda, n.states, theta, W) {
  res <- LambdaEval(lambda, W, n.states) 
  loss <- sqrt(sum((theta - res)^2))
  if (!is.finite(loss)) {
    stop()
  }  
  return(loss)
}

#' @title GetInitGues()
#' @export
##generate initial lambda if likelihood is infinite
GetInitGues <- function(W, n.states) {
  l0 <- seq(1, 0, length=(n.states-1))
  N <- length(W)
  
  l0 <- l0 * runif(1, N, 3*N)
  sgn <- ifelse(rbinom(1, 1, 0.5) == 1, -1, 1)
  return(l0 * sgn)
}


#' @title LambdaEval()
#' @export
LambdaEval <- function(lambda, W, n.states) {
  log.p.q.terms <- LambdaExpr(lambda, W)  
  p.q.terms <- exp(log.p.q.terms)
  p.q.terms[p.q.terms == Inf] <- .Machine[["double.xmax"]] / 2
  seg.totals <- apply(p.q.terms, 2, sum)
  seg.densities <-  t(p.q.terms) / seg.totals
  g1 <- apply(W * seg.densities, 2, sum) 
  if (any(!is.finite(g1))) {
    stop()
  }
  
  return(g1)
}


#' @title LambdaExpr()
#' @export

LambdaExpr <- function(lambda, W) {
  n.states <- length(lambda) + 1
  return(-c(0, lambda) * ((c(1:n.states)-1) %*% matrix(W, nrow=1)))
}


#' @title MargModeFinder()
#' @export

MargModeFinder <- function(obs, dom2, Q, lambda.qz.res, pz, sigma.h, tau_dom, d.res=0.025, verbose=FALSE) {
  #d.res=0.125
  delta_dom = log(c(1 / tau_dom[2] - 0.05, 1))
  mode.tab = run_1d_opt(obs, delta_dom, d.res, lambda.qz.res, pz, sigma.h, Q, verbose=verbose)  #line 244
  if (verbose) {
    print("1d mode opt: ")
    print(res_1d)
  }
  
  ## find unique modes in table
  ix <- !is.na(mode.tab[, 1])
  mode.list <- mode.tab[ix, 1, drop = FALSE]
  umodes = unique(round(mode.list, 2))
  
  if (verbose) {
    print(paste(nrow(umodes), " unique modes found", sep = ""))
  }
  umodes <- umodes[umodes[, 1] >= dom2[1], , drop = FALSE]
  if (verbose) {
    print(paste(nrow(umodes), " modes in delta range.", sep = ""))
  }
  return(umodes)
}


#' @title run_1d_opt()
#' @export

run_1d_opt = function(obs, dom, d_res, lambda_qz_res, pz, sigma_h, Q, verbose=FALSE) {
  comb_1d_ll = function(par, Q, obs, lambda_qz_res, pz, sigma_h) {
    delta = exp(par)
    comb =  GetCopyRatioComb(Q, delta) #relative ratio corresponding to each integer copy number
    LL = CalcNormLoglik(obs$ratio, obs$sd, obs$w, comb, lambda_qz_res, pz, sigma_h) #lian 219
    return(-LL)
  }
  
  d_grid = seq(dom[1], dom[2], d_res) #ploidy distance grid
  mode_tab = array(NA, dim=c(length(d_grid), 2))
  
  for (i in seq_along(d_grid)) {
    opt = nlm(f=comb_1d_ll, p=d_grid[i], Q=Q, obs=obs, lambda_qz_res=lambda_qz_res, pz=pz, sigma_h=sigma_h)
    mode_tab[i, ] = c(opt$estimate, -opt$minimum)
    if (verbose) {
      cat("-")
    }
  }
  
  if (verbose) {
    cat("\n")
  }
  
  return(mode_tab)
}


#' @title GetCopyRatioComb()
#' @export

GetCopyRatioComb = function(Q, delta) {
  xx = (delta * (c(1:Q) - 1))
  #means = Atten(xx, error_model$fit.at)
  #tx_means = TxData(error_model, means)
  
  return(xx)
}


#' @title CalcNormLoglik()
#' @export
CalcNormLoglik <- function(d, sigma, w, comb, lambda.qz.res, pz, sigma.h, comb.s=1) {
  if (is.list(lambda.qz.res)) {
    log.pqz <- t((LambdaExpr(lambda.qz.res[["Lambda"]], w))) +
      log(lambda.qz.res[["norm.term"]])
    q <- length(lambda.qz.res[["Lambda"]])
    comb <- comb[1:q]
  } else {
    q <- length(comb)
    log.pqz <- matrix(log(pz), nrow = length(w), ncol = q + 1, byrow = TRUE)
  }
  
  log.ss.p.1 <- sapply(comb, dnorm, d, sqrt(sigma^2 + comb.s^2 * sigma.h^2),
                       log=TRUE)
  log.comb.seg.mat <- cbind(log.pqz[, c(1:q)] + log.ss.p.1,
                            log(1 / 7) + log.pqz[, q + 1])  
  LL <- sum(LogAdd(log.comb.seg.mat))
  
  return(LL)
}


#' @title LogAdd()
#' @export

LogAdd <- function(X) {
  ##  Calculates log( sum(exp(x)) )  without 'leaving' log space
  if (is.vector(X)) {
    mix <- which.max(X)
    max <- X[mix]
    res <- max + log(sum(exp(X - max)))
  }
  
  if (is.matrix(X)) {
    mv <- apply(X, 1, max)
    res <- mv + log(rowSums(exp(X - mv)))
  }
  
  return(res)
}


#' @title CombLL()
#' @export

CombLL <- function(par, Q, obs, dom2, lambda.qz.res, pz, sigma.h) {
  if (any(is.na(par))) {
    cat("$")
    return(Inf)
  }
  
  if (par[1] < log(0.07) | par[1] > dom2[2]) {
    return(Inf)
  }
  delta <- exp(par[1])
  comb <- GetCopyRatioComb(Q, delta)
  LL <- CalcNormLoglik(obs$ratio, obs$sd, obs$w, comb,lambda.qz.res, pz, sigma.h)
  return(-LL)
}


#' @title CalcModeLogCurv()
#' @export
CalcModeLogCurv <- function(mode, mode.hess, verbose=FALSE) {
  hess.mat <- mode.hess
  curvature <- abs(det(hess.mat / (2 * pi)))
  mode.curv <- log((curvature)^(-1/2) ) - log(2)
  if ((!is.finite(mode.curv)) && verbose) {
    print("WARNING: NON-FINITE log_evidence")
  }
  return(mode.curv)
}


#' @title GetTau()
#' @export
GetTau = function(delta) {
  tau = 1 / delta
  return(list(tau=tau)) 
}


#' @title FindLocationModes()
#' @export
##FindLocationModes: find the integer coly number distance
FindLocationModes <- function(obs, Q, theta.qz, sigma.h, tau.dom, verbose=FALSE) {
  kDom2 <- log(c(0.08, 1.05)) 
  ## lambda estimation
  lambda.qz.res <- GetLambda(theta.qz, obs$w, Q + 1, verbose=verbose)
  pz <- theta.qz[Q + 1]
  
  #####
  mode.tab <- MargModeFinder(obs, kDom2, Q, lambda.qz.res, pz, sigma.h,tau.dom, verbose=verbose)
  #############################################
  if (nrow(mode.tab) == 0) {
    return(list(mode.flag = "DELTA_DOM"))
  }
  
  mode.tab <- cbind(mode.tab, NA, NA)
  
  for (i in 1:nrow(mode.tab)) {
    par <- mode.tab[i, 1]
    cur.par <- par
    mode.params <- list(delta = par)
    
    LL <- CombLL(cur.par, Q = Q, obs = obs, dom2 = kDom2,lambda.qz.res, pz, sigma.h)
    
    mode.hess <- hessian(CombLL, x = cur.par, method = "Richardson", Q = Q, 
                         obs = obs, dom2 = kDom2,
                         lambda.qz.res = lambda.qz.res, pz = pz, 
                         sigma.h = sigma.h)
    
    if (!is.na(mode.hess[1])) {
      mode.curv <- CalcModeLogCurv(par, mode.hess, verbose=verbose)
    } else {
      LL <- NA
      mode.curv <- NA
    }
    mode.tab[i, ] <- c(mode.params[["delta"]],LL, mode.curv)
  }
  delta <- exp(mode.tab[, 1])
  res <- GetTau(delta)
  tau <- res[["tau"]]
  LL <- mode.tab[, 2]
  mode.curv <- mode.tab[, 3]
  mode.tab <- cbind(tau, delta, LL, mode.curv)
  colnames(mode.tab) <- c("tau","delta", "LL", "mode_curv")
  
  if (verbose) {
    print(mode.tab)
  }
  
  mode.ix <- (mode.tab[, "tau"] >= tau.dom[1] &mode.tab[, "tau"] <= tau.dom[2])
  
  if (verbose) {
    print(paste("removing ", sum(!mode.ix), " / ",
                length(mode.ix), " modes outside of tau range.", 
                sep = ""))
  }
  mode.tab <- mode.tab[mode.ix, , drop = FALSE]
  
  if (nrow(mode.tab) == 0) {
    return(list(mode.flag = "TAU_DOM"))
  }
  return(list(mode.tab = mode.tab))
}

