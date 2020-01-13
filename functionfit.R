################
DataSetzeroNB <- function (n, p, K, param, sdsignal, DE,pzero=NA) 
{
  if (n < 4 * K) 
    stop("We require n to be at least 4*K.")
  
  q0 <- rexp(p, rate = 1/25)
  isDE <-  runif(p)< DE    
  classk <- matrix(NA, nrow = K, ncol = p)
  for (k in 1:K) {
    lfc <- rnorm(p, mean=0, sd = sdsignal)
    classk[k, ] <- ifelse(isDE, q0 * exp(lfc), q0)
  }
  truesf <- runif(n)*20+20
  truesfte <- runif(n)*20+20
  conds <- sample(c(rep(1:K, 4), sample(1:K, n - 4 * K, replace = TRUE)))
  condste <- sample(c(rep(1:K, 4), sample(1:K, n - 4 * K, replace = TRUE)))
  x <- xte <- matrix(NA, nrow = n, ncol = p)
  if (is.null(pzero)) {
    prob<- runif(n,0.05,0.5)
  }else{
    prob<-rep(pzero,n)
  }
  
  for (i in 1:n) {
    for (k in 1:K) {
      if (conds[i] == k){
        x[i, ] <- rnbinom(p, mu = truesf[i] * classk[k, ], size = param)
        index<-which(x[i,]>0)
        if(sum(x[i,]>0)>(p*prob[i])){
          ind0<-sample(index,round(p*prob[i]), replace = FALSE)
          x[i,ind0]<-0
        }else{
          x[i,]<-0
          x[i,sample(1:p,1)]<-1
        }
      }
      
      if (condste[i] == k){
        xte[i, ] <- rnbinom(p, mu = truesfte[i] * classk[k, ], size = param)
        index<-which(xte[i,]>0)
        if(sum(xte[i,]>0)>(p*prob[i])){
          ind0<-sample(index,round(p*prob[i]), replace = FALSE)
          xte[i,ind0]<-0
        }else{
          xte[i,]<-0
          xte[i,sample(1:p,1)]<-1
          
        }
      }
    }
  }

  rm <- apply(x, 2, sum) == 0
  return(list(x = x[, !rm], xte = xte[, !rm], y = conds, yte = condste, 
              truesf = truesf, truesfte = truesfte, dkg=classk))
}

#####################
estimatepZIP<-function(x,y,xte=NULL,yte,beta=1,type=c("mle","deseq","quantile"), prior=NULL){
  if(is.null(xte)){
    xte <- x
    warning("Since no xte was provided, testing was performed on training data set.")
  }
  type <- match.arg(type)
  if(is.null(prior)) prior <- rep(1/length(unique(y)), length(unique(y)))
  
  null.out <- NullModel(x, type=type)
  ns <- null.out$n
  nste <- NullModelTest(null.out,x,xte,type=type)$nste
  uniq <- sort(unique(y))
  ds <-  GetDn(ns,x,y,beta)
  
  mu <- matrix(NA, nrow=nrow(x), ncol=length(x[1,]))
  mute <- matrix(NA, nrow=nrow(xte), ncol=length(xte[1,]))
  for(i in 1:nrow(x)){
    dstar = ds[y[i],]
    mu[i,] <- ns[i,]*dstar
  }
  
  for(i in 1:nrow(xte)){
    dstar = ds[yte[i],]
    mute[i,] <- nste[i,]*dstar
  }
  
  
  G <- length(x[1,])
  lib <- rowSums(x)
  x1 <- t(x)
  
  mu1 <- as.vector(t(mu))
  librep <- rep(lib,rep(G,length(lib)))/(lib[1])
  x2 <- as.vector(x1)
  
  y <- x2
  y[y!=0] <- 1
  xreg <- cbind(y,librep,mu1)
  glm.out <- glm(y~ librep+ mu1,family=binomial("logit"),data=data.frame(xreg))
  summary(glm.out)
  
  coef <- as.matrix(glm.out$coefficients)
  inter <- rep(1,G)
  muu <- t(mute)
  libte <- rowSums(xte)
  
  xte1 = t(xte)
  p <- xte1
  for(i in 1:length(xte1[1,])){
    libsize <- rep(sum(xte1[,i]),G)/libte[1]
    estx1 <- cbind(inter,libsize,muu[,i])
    dd <- estx1%*% coef
    dd[dd>50] <- 50
    dd[dd<(-50)] <- -50
    p1<-exp(dd)
    
    p[,i] <- ((1-exp(-muu[,i])*(1+p1))/((1+p1)*(1-exp(-muu[,i]))))
  }
  p[p<0] <- 0
  pp <- p
  return(pp)   
}



estimatepZINB <- function(x,y,xte=NULL,yte,phihat,beta=1,type=c("mle","deseq","quantile"), prior=NULL){
  if(is.null(xte)){
    xte <- x
    warning("Since no xte was provided, testing was performed on training data set.")
  }
  type <- match.arg(type)
  if(is.null(prior)) prior <- rep(1/length(unique(y)), length(unique(y)))
  
  null.out <- NullModel(x, type=type)
  ns <- null.out$n
  nste <- NullModelTest(null.out,x,xte,type=type)$nste
  uniq <- sort(unique(y))
  ds <-  GetDn(ns,x,y,beta)
  
  md <- matrix(NA, nrow=nrow(x), ncol=length(x[1,]))
  mdte <- matrix(NA, nrow=nrow(xte), ncol=length(xte[1,]))
  for(i in 1:nrow(x)){
    dstar = ds[y[i],]
    md[i,] <- ns[i,]*dstar
  }
  
  for(i in 1:nrow(xte)){
    dstar = ds[yte[i],]
    mdte[i,] <- nste[i,]*dstar
  }
  
  
  G <- length(x[1,])
  lib <- rowSums(x)
  x1 <- t(x)

  md1 <- as.vector(t(md))
  librep <- rep(lib,rep(G,length(lib)))/(lib[1])
  x2 <- as.vector(x1)
  
  y <- x2
  y[y!=0] <- 1
  xreg <- cbind(y,librep,md1)
  glm.out <- glm(y~ librep+ md1,family=binomial("logit"),data=data.frame(xreg))
  summary(glm.out)
  
  coef <- as.matrix(glm.out$coefficients)
  inter <- rep(1,G)
  mdd <- t(mdte)
  libte <- rowSums(xte)
  
  xte1 = t(xte)
  p <- xte1
  for(i in 1:length(xte1[1,])){
    
    libsize <- rep(sum(xte1[,i]),G)/libte[1]
    estx1 <- cbind(inter,libsize,mdd[,i])
    dd <- estx1%*% coef
    dd[dd > 50] <- 50
    dd[dd < (-50)] <- -50
    p1 <- (1/(1+mdd[,i]*phihat))^(1/phihat)
    p2 <- exp(dd)
    p[,i] <- (1-(1+p2)*p1)/((1+p2)*(1-p1))
  }
  p[p<0] <- 0
  pp <- p
  return(pp)   
}
######################
#ZIPLDA
ZIPDA.cv <-
  function(x,y,rhos=NULL,beta=1,nfolds=5,prob0=NULL,type=c("mle","deseq","quantile"),folds=NULL,transform=TRUE, alpha=NULL, prior=NULL){
    type <- match.arg(type)
    if(!transform && !is.null(alpha)) stop("You have asked for NO transformation but have entered alpha.")
    if(transform && is.null(alpha)) alpha <- FindBestTransform(x)
    if(transform){
      if(alpha<=0 || alpha>1) stop("alpha must be between 0 and 1")
      x <- x^alpha
    }
    if(is.null(rhos)){
      ns <- NullModel(x,type=type)$n
      uniq <- sort(unique(y))
      maxrho <- rep(NA, length(uniq))
      for(k in 1:length(uniq)){
        a <- colSums(x[y==uniq[k],])+beta
        b <- colSums(ns[y==uniq[k],])+beta
        maxrho[k] <- max(abs(a/b-1)*sqrt(b),na.rm=TRUE)
      }
      rhos <- seq(0, max(maxrho,na.rm=TRUE)*(2/3), len=30)
    }
    if(is.null(folds)) folds <- balanced.folds(y,nfolds=nfolds)
    nfolds <- length(folds)
    errs <- nnonzero <- matrix(NA, nrow=nfolds, ncol=length(rhos))
    for(i in 1:nfolds){
      cat(i,fill=FALSE)
      tr <- -folds[[i]]
      te <- folds[[i]]
      out <- ZIPLDA(x[tr,],y[tr],x[te,],rhos=rhos,beta=beta,prob0=prob0,type="mle", prior=prior, transform=FALSE) # Have already power-transformed x, so don't need to do it again!!!
      for(j in 1:length(rhos)){      
        errs[i,j] <- sum(out[[j]]$ytehat!=y[te])
        nnonzero[i,j] <- sum(colSums(out[[j]]$ds!=1)!=0)
      }
    }
    cat(fill=TRUE)
    save <- list(errs=errs, bestrho=rhos[max(which(colMeans(errs)==min(colMeans(errs))))], rhos=rhos, nnonzero=nnonzero,folds=folds, alpha=alpha,type=type)
    return(save)
  }

######################
ZIPLDA<-
  function(x,y,xte=NULL,rho=0,beta=1,rhos=NULL,prob0=NULL,type=c("mle","deseq","quantile"), prior=NULL, transform=TRUE, alpha=NULL){
    if(is.null(xte)){
      xte <- x
      warning("Since no xte was provided, testing was performed on training data set.")
    }
    if(!is.null(rho) && length(rho)>1) stop("Can only enter 1 value of rho. If you would like to enter multiple values, use rhos argument.")
    type <- match.arg(type)
    if(!transform && !is.null(alpha)) stop("You have asked for NO transformation but have entered alpha.")
    if(transform && is.null(alpha)) alpha <- FindBestTransform(x)
    if(transform){
      if(alpha<=0 || alpha>1) stop("alpha must be between 0 and 1")
      x <- x^alpha
      xte <- xte^alpha
    }  
    if(is.null(prior)) prior <- rep(1/length(unique(y)), length(unique(y)))
    if(is.null(rho)&&is.null(rhos)) stop("Must enter rho or rhos.")
    null.out <- NullModel(x, type=type)
    ns <- null.out$n
    nste <- NullModelTest(null.out,x,xte,type=type)$nste
    uniq <- sort(unique(y))
    signx1<-sign(xte==0)
    if(is.null(rhos)){
      ds <- GetD(ns,x,y,rho,beta)
      discriminant <- matrix(NA, nrow=nrow(xte), ncol=length(uniq))
      for(k in 1:length(uniq)){
        for(i in 1:nrow(xte)){
          dstar = ds[k,]
          part2=nste[i,]*dstar 
          part1=prob0[i,]+(1-prob0[i,])*exp(-part2) 
          part1[part1==0]=1
          discriminant[i,k] <-sum(signx1[i,]*log(part1))+sum(xte[i,]*(1-signx1[i,])*log(dstar))-sum((1-signx1[i,])*part2)+log(prior[k])
        }
      }
      save <- list(ns=ns,nste=nste,ds=ds,discriminant=discriminant,ytehat=apply(discriminant,1,which.max), alpha=alpha, rho=rho,x=x,y=y,xte=xte,alpha=alpha,type=type)
      return(save)
    } else {
      save <- list()
      ds.list <- GetD(ns,x,y,rho=NULL, rhos=rhos,beta)
      for(rho in rhos){
        ds <- ds.list[[which(rhos==rho)]]
        discriminant <- matrix(NA, nrow=nrow(xte), ncol=length(uniq))
        for(k in 1:length(uniq)){
          for(i in 1:nrow(xte))   {
           
            dstar = ds[k,]
            part2=nste[i,]*dstar 
            part1=prob0[i]+(1-prob0[i])*exp(-part2) 
            part1[part1==0]=1
            discriminant[i,k] <-sum(signx1[i,]*log(part1))+sum(xte[i,]*(1-signx1[i,])*log(dstar))-sum((1-signx1[i,])*part2)+log(prior[k])
            
          }
          
        }
        save[[which(rhos==rho)]] <- (list(ns=ns,nste=nste,ds=ds,discriminant=discriminant,ytehat=apply(discriminant,1,which.max), alpha=alpha, rho=rho,x=x,y=y,xte=xte,alpha=alpha,type=type))
      }
      return(save)
    }
  }


####################
#PLDA
PLDA.cv <-
  function(x,y,rhos=NULL,beta=1,nfolds=5,type=c("mle","deseq","quantile"),folds=NULL,transform=TRUE, alpha=NULL, prior=NULL){
    type <- match.arg(type)
    if(!transform && !is.null(alpha)) stop("You have asked for NO transformation but have entered alpha.")
    if(transform && is.null(alpha)) alpha <- FindBestTransform(x)
    if(transform){
      if(alpha<=0 || alpha>1) stop("alpha must be between 0 and 1")
      x <- x^alpha
    }
    if(is.null(rhos)){
      ns <- NullModel(x,type=type)$n
      uniq <- sort(unique(y))
      maxrho <- rep(NA, length(uniq))
      for(k in 1:length(uniq)){
        a <- colSums(x[y==uniq[k],])+beta
        b <- colSums(ns[y==uniq[k],])+beta
        maxrho[k] <- max(abs(a/b-1)*sqrt(b),na.rm=TRUE)
      }
      rhos <- seq(0, max(maxrho,na.rm=TRUE)*(2/3), len=30)
    }
    if(is.null(folds)) folds <- balanced.folds(y,nfolds=nfolds)
    nfolds <- length(folds)
    errs <- nnonzero <- matrix(NA, nrow=nfolds, ncol=length(rhos))
    for(i in 1:nfolds){
      cat(i,fill=FALSE)
      tr <- -folds[[i]]
      te <- folds[[i]]
      out <- Classify(x[tr,],y[tr],x[te,],rhos=rhos,beta=beta,type="quantile", prior=prior, transform=FALSE) # Have already power-transformed x, so don't need to do it again!!!
      for(j in 1:length(rhos)){      
        errs[i,j] <- sum(out[[j]]$ytehat!=y[te])
        nnonzero[i,j] <- sum(colSums(out[[j]]$ds!=1)!=0)
      }
    }
    cat(fill=TRUE)
    save <- list(errs=errs, bestrho=rhos[max(which(colMeans(errs)==min(colMeans(errs))))], rhos=rhos, nnonzero=nnonzero,folds=folds, alpha=alpha,type=type)
    return(save)
  }

########################
PLDA <-
  function(x,y,xte=NULL,rho=0,beta=1,rhos=NULL,type=c("mle","deseq","quantile"), prior=NULL, transform=TRUE, alpha=NULL){
    if(is.null(xte)){
      xte <- x
      warning("Since no xte was provided, testing was performed on training data set.")
    }
    if(!is.null(rho) && length(rho)>1) stop("Can only enter 1 value of rho. If you would like to enter multiple values, use rhos argument.")
    type <- match.arg(type)
    if(!transform && !is.null(alpha)) stop("You have asked for NO transformation but have entered alpha.")
    if(transform && is.null(alpha)) alpha <- FindBestTransform(x)
    if(transform){
      if(alpha<=0 || alpha>1) stop("alpha must be between 0 and 1")
      x <- x^alpha
      xte <- xte^alpha
    }  
    if(is.null(prior)) prior <- rep(1/length(unique(y)), length(unique(y)))
    if(is.null(rho)&&is.null(rhos)) stop("Must enter rho or rhos.")
    null.out <- NullModel(x, type=type)
    ns <- null.out$n
    nste <- NullModelTest(null.out,x,xte,type=type)$nste
    uniq <- sort(unique(y))
    if(is.null(rhos)){
      ds <- GetD(ns,x,y,rho,beta)
      discriminant <- matrix(NA, nrow=nrow(xte), ncol=length(uniq))
      for(k in 1:length(uniq)){
        discriminant[,k] <- rowSums(scale(xte,center=FALSE,scale=(1/log(ds[k,])))) - rowSums(scale(nste,center=FALSE,scale=(1/ds[k,]))) + log(prior[k])
      }
      save <- list(ns=ns,nste=nste,ds=ds,discriminant=discriminant,ytehat=apply(discriminant,1,which.max), alpha=alpha, rho=rho,x=x,y=y,xte=xte,alpha=alpha,type=type)
      return(save)
    } else {
      save <- list()
      ds.list <- GetD(ns,x,y,rho=NULL, rhos=rhos,beta)
      for(rho in rhos){
        ds <- ds.list[[which(rhos==rho)]]
        discriminant <- matrix(NA, nrow=nrow(xte), ncol=length(uniq))
        for(k in 1:length(uniq)){
          discriminant[,k] <- rowSums(scale(xte,center=FALSE,scale=(1/log(ds[k,])))) - rowSums(scale(nste,center=FALSE,scale=(1/ds[k,]))) + log(prior[k])
        }
        save[[which(rhos==rho)]] <- (list(ns=ns,nste=nste,ds=ds,discriminant=discriminant,ytehat=apply(discriminant,1,which.max), alpha=alpha, rho=rho,x=x,y=y,xte=xte,alpha=alpha,type=type))
      }
      return(save)
    }
  }
###################
#NBLDA
NBLDA.cv <-
  function(x,y,rhos=NULL,beta=1,nfolds=5,phihat=0,type=c("mle","deseq","quantile"),folds=NULL, prior=NULL){
    type <- match.arg(type)
    
    if(is.null(rhos)){
      ns <- NullModel(x,type=type)$n
      uniq <- sort(unique(y))
      maxrho <- rep(NA, length(uniq))
      for(k in 1:length(uniq)){
        a <- colSums(x[y==uniq[k],])+beta
        b <- colSums(ns[y==uniq[k],])+beta
        maxrho[k] <- max(abs(a/b-1)*sqrt(b),na.rm=TRUE)
      }
      rhos <- seq(0, max(maxrho,na.rm=TRUE)*(2/3), len=30)
    }
    if(is.null(folds)) folds <- balanced.folds(y,nfolds=nfolds)
    nfolds <- length(folds)
    errs <- nnonzero <- matrix(NA, nrow=nfolds, ncol=length(rhos))
    for(i in 1:nfolds){
      cat(i,fill=FALSE)
      tr <- -folds[[i]]
      te <- folds[[i]]
      out <- NBLDA(x[tr,],y[tr],x[te,],rhos=rhos,phihat=phihat,beta=beta,type="mle", prior=prior) # Have already power-transformed x, so don't need to do it again!!!
      for(j in 1:length(rhos)){      
        errs[i,j] <- sum(out[[j]]$ytehat!=y[te])
        nnonzero[i,j] <- sum(colSums(out[[j]]$ds!=1)!=0)
      }
    }
    cat(fill=TRUE)
    save <- list(errs=errs, bestrho=rhos[max(which(colMeans(errs)==min(colMeans(errs))))], rhos=rhos, nnonzero=nnonzero,folds=folds,type=type)
    return(save)
  }
###################
NBLDA <-
  function(x,y,xte=NULL,rho=0,beta=1,rhos=NULL,phihat=0,type=c("mle","deseq","quantile"), prior=NULL){
    if(is.null(xte)){
      xte <- x
      warning("Since no xte was provided, testing was performed on training data set.")
    }
    
    if(is.null(prior)) prior <- rep(1/length(unique(y)), length(unique(y)))
    if(is.null(rho)&&is.null(rhos)) stop("Must enter rho or rhos.")
    null.out <- NullModel(x, type=type)
    ns <- null.out$n
    nste <- NullModelTest(null.out,x,xte,type=type)$nste
    uniq <- sort(unique(y))
    if(is.null(rhos)){
      ds <- GetD(ns,x,y,rho,beta)
      discriminant <- matrix(NA, nrow=nrow(xte), ncol=length(uniq))
      for (k in 1:length(uniq)) {
        
        for(l in 1:nrow(xte))   {
          
          dstar = ds[k,]
          part2=1+nste[l,]*dstar*phihat 
          part1=dstar/part2 
          
          
          discriminant[l, k]<- sum(xte[l,]*log(part1))-sum((1/phihat)*log(part2))+log(prior[k])
          
        }
      }
      save <- list(ns=ns,nste=nste,ds=ds,discriminant=discriminant,ytehat=apply(discriminant,1,which.max), rho=rho,x=x,y=y,xte=xte,type=type)
      return(save)
    } else {
      save <- list()
      ds.list <- GetD(ns,x,y,rho=NULL, rhos=rhos,beta)
      for(rho in rhos){
        ds <- ds.list[[which(rhos==rho)]]
        discriminant <- matrix(NA, nrow=nrow(xte), ncol=length(uniq))
        for (k in 1:length(uniq)) {
          
          for(l in 1:nrow(xte))   {
            
            dstar = ds[k,]
            part2=1+nste[l,]*dstar*phihat 
            part1=dstar/part2 
            
            
            discriminant[l, k]<- sum(xte[l,]*log(part1))-sum((1/phihat)*log(part2))+log(prior[k])
            
          }
        }
        save[[which(rhos==rho)]] <- (list(ns=ns,nste=nste,ds=ds,discriminant=discriminant,ytehat=apply(discriminant,1,which.max), rho=rho,x=x,y=y,xte=xte,type=type))
      }
      return(save)
    }
  }
########################
GetDn <- function(ns, x, y, beta){
  uniq <- sort(unique(y))
  ds <- matrix(1, nrow=length(uniq), ncol=ncol(x))
  for(k in 1:length(uniq)){
    a <- colSums(x[y==uniq[k],])+beta
    b <- colSums(ns[y==uniq[k],])+beta
    ds[k,] <- a/b
  }
  return(ds)
}


###############
#estimated parameter based on the ZINBLDA model
fun <- function(x,mu,disperhatmean){
  
  
  signx <- sign(x==0)
  zinb <- function(p) {
    
    res <- sum(sign(x==0)*log(p[1]+(1-p[1])*(1/(1+p[3]*p[2]))^(1/p[2]))+(1-sign(x==0))*(log(1-p[1])+lgamma(x+1/p[2])-lgamma(x+1)-lgamma(1/p[2])+x*log(p[3]*p[2])-x*log(1+p[3]*p[2])-(1/p[2])*log(1+p[3]*p[2])))
    
    return(-res)
    
  }
  nlminb(c(0,1,mu),zinb,lower=c(0.01,0.001,0.01),upper=c(0.999,disperhatmean,9999),control = list(step.max=0.2))$par
}

################
#ZINBLDA
ZINB.cv <-
  function(x,y,rhos=NULL,beta=1,nfolds=5,phihat=0,prob0=NULL,type=c("mle","deseq","quantile"),folds=NULL, prior=NULL){
    type <- match.arg(type)
    if(is.null(rhos)){
      ns <- NullModel(x,type=type)$n
      uniq <- sort(unique(y))
      maxrho <- rep(NA, length(uniq))
      for(k in 1:length(uniq)){
        a <- colSums(x[y==uniq[k],])+beta
        b <- colSums(ns[y==uniq[k],])+beta
        maxrho[k] <- max(abs(a/b-1)*sqrt(b),na.rm=TRUE)
      }
      rhos <- seq(0, max(maxrho,na.rm=TRUE)*(2/3), len=30)
    }
    if(is.null(folds)) folds <- balanced.folds(y,nfolds=nfolds)
    nfolds <- length(folds)
    errs <- nnonzero <- matrix(NA, nrow=nfolds, ncol=length(rhos))
    for(i in 1:nfolds){
      cat(i,fill=FALSE)
      tr <- -folds[[i]]
      te <- folds[[i]]
      out <-ZINBLDA(x[tr,],y[tr],x[te,],rhos=rhos,phihat=phihat,beta=beta,prob0=prob0,type="mle", prior=prior) # Have already power-transformed x, so don't need to do it again!!!
      for(j in 1:length(rhos)){      
        errs[i,j] <- sum(out[[j]]$ytehat!=y[te])
        nnonzero[i,j] <- sum(colSums(out[[j]]$ds!=1)!=0)
      }
    }
    cat(fill=TRUE)
    save <- list(errs=errs, bestrho=rhos[max(which(colMeans(errs)==min(colMeans(errs))))], rhos=rhos, nnonzero=nnonzero,folds=folds,type=type)
    return(save)
  }

######################
ZINBLDA<-
  function(x,y,xte=NULL,rho=0,beta=1,rhos=NULL,phihat=0,prob0=NULL,type=c("mle","deseq","quantile"), prior=NULL){
    if(is.null(xte)){
      xte <- x
      warning("Since no xte was provided, testing was performed on training data set.")
    }
    if(is.null(prior)) prior <- rep(1/length(unique(y)), length(unique(y)))
    if(is.null(rho)&&is.null(rhos)) stop("Must enter rho or rhos.")
    null.out <- NullModel(x, type=type)
    ns <- null.out$n
    nste <- NullModelTest(null.out,x,xte,type=type)$nste
    uniq <- sort(unique(y))
    signx3<-sign(xte==0)
    if(is.null(rhos)){
      ds <- GetD(ns,x,y,rho,beta)
      discriminant <- matrix(NA, nrow=nrow(xte), ncol=length(uniq))
      for(k in 1:length(uniq)){
        for(i in 1:nrow(xte)){
          dstar = ds[k,]
          part2=nste[i,]*dstar 
          part3<-(1/(1+part2*phihat))^(1/phihat)
          discriminant[i,k] <-sum(signx3[i,]*log(prob0[i,]+(1-prob0[i,])*part3))+sum(xte[i,]*(1-signx3[i,])*(log(dstar)-log(1+part2*phihat)))-sum((1-signx3[i,])*(1/phihat)*log(1+part2*phihat))+log(prior[k])
          
        }
      }
      save <- list(ns=ns,nste=nste,ds=ds,discriminant=discriminant,ytehat=apply(discriminant,1,which.max),rho=rho,x=x,y=y,xte=xte,type=type)
      return(save)
    } else {
      save <- list()
      ds.list <- GetD(ns,x,y,rho=NULL, rhos=rhos,beta)
      for(rho in rhos){
        ds <- ds.list[[which(rhos==rho)]]
        discriminant <- matrix(NA, nrow=nrow(xte), ncol=length(uniq))
        for(k in 1:length(uniq)){
          for(i in 1:nrow(xte))   {
            
            dstar = ds[k,]
            part2=nste[i,]*dstar 
            part3<-(1/(1+part2*phihat))^(1/phihat)
            discriminant[i,k] <-sum(signx3[i,]*log(prob0[i,]+(1-prob0[i,])*part3))+sum(xte[i,]*(1-signx3[i,])*(log(dstar)-log(1+part2*phihat)))-sum((1-signx3[i,])*(1/phihat)*log(1+part2*phihat))+log(prior[k])
            
          }
        }
        save[[which(rhos==rho)]] <- (list(ns=ns,nste=nste,ds=ds,discriminant=discriminant,ytehat=apply(discriminant,1,which.max),rho=rho,x=x,y=y,xte=xte,type=type))
      }
      return(save)
    }
  }
###################################
calcNormFactors <- function(dataMatrix, refColumn=1, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10) {
  if( !is.matrix(dataMatrix) )
    stop("'dataMatrix' needs to be a matrix")
  if( refColumn > ncol(dataMatrix) )
    stop("Invalid 'refColumn' argument")
  apply(dataMatrix,2,.calcFactorWeighted,ref=dataMatrix[,refColumn], logratioTrim=logratioTrim, 
        sumTrim=sumTrim, doWeighting=doWeighting, Acutoff=Acutoff)
}

.calcFactorWeighted <- function(obs, ref, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10) {
  
  if( all(obs==ref) )
    return(1)
  
  nO <- sum(obs)
  nR <- sum(ref)
  logR <- log2((obs/nO)/(ref/nR))          # log ratio of expression, accounting for library size
  absE <- (log2(obs/nO) + log2(ref/nR))/2  # absolute expression
  v <- (nO-obs)/nO/obs + (nR-ref)/nR/ref   # estimated asymptotic variance
  
  # remove infinite values, cutoff based on A
  fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)
  
  logR <- logR[fin]
  absE <- absE[fin]
  v <- v[fin]
  
  # taken from the original mean() function
  n <- sum(fin)
  loL <- floor(n * logratioTrim) + 1
  hiL <- n + 1 - loL
  loS <- floor(n * sumTrim) + 1
  hiS <- n + 1 - loS
  
  keep <- (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
  if (doWeighting) 
    2^( sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE) )
  else
    2^( mean(logR[keep], na.rm=TRUE) )
}
##########
balanced.folds <- function(y, nfolds = min(min(table(y)), 10)){
  totals <- table(y)
  fmax <- max(totals)
  nfolds <- min(nfolds, fmax)
  # makes no sense to have more folds than the max class size
  folds <- as.list(seq(nfolds))
  yids <- split(seq(y), y)
  # nice way to get the ids in a list, split by class
  ###Make a big matrix, with enough rows to get in all the folds per class
  bigmat <- matrix(NA, ceiling(fmax/nfolds) * nfolds, length(totals))
  for(i in seq(totals)) {
    bigmat[seq(totals[i]), i] <- sample(yids[[i]])
  }
  smallmat <- matrix(bigmat, nrow = nfolds) # reshape the matrix
  ### Now do a clever sort to mix up the NAs
  smallmat <- permute.rows(t(smallmat)) ### Now a clever unlisting
  x <- apply(smallmat, 2, function(x) x[!is.na(x)])
  if(is.matrix(x)){
    xlist <- list()
    for(i in 1:ncol(x)){
      xlist[[i]] <- x[,i]
    }
    return(xlist)
  }
  return(x)
}
############################