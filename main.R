library(foreach)
library(doParallel)
library("e1071")
library("class")
library("sSeq")
library("limma")
library("edgeR")
library("PoiClaClu")

source("E:/RNA-seq/yb/论文/functionfit.r")


#parameters of simulated data
dispersion <- 0.001   # dispersion parameter φ
DE <- 0.2             # differential expression ratio
n = seq(8,64,8)       # the training and testing sample size
pzero = 0.2           # the percentage of excess of zeros in the data
sdsignal = 0.5        # The variances of the log(d)
p <- 100              # Number of genes
k = 2                 # Number of classes

simtime = 50          # The simulation number of times
variable = n
countP <- countNB <- countZIP <- countZINB <- rep(0,length(variable))




myfun <- function(i){
  for(j in 1:length(variable)){
    
    #Generate simulated data
    dat <- DataSetzeroNB(n=n[j],p=p,sdsignal=sdsignal,K=k,param=1/dispersion,DE=DE,pzero=pzero)
    
    #ZIP model estimates the probability of zero in data
    prob0_ZIP_train <- estimatepZIP(dat$x,dat$y,xte=dat$x,dat$y,beta=1,type=c("mle","deseq","quantile"), prior=NULL)
    prob0_ZIP_test <- estimatepZIP(dat$x,dat$y,dat$xte,dat$yte,beta=1,type=c("mle","deseq","quantile"), prior=NULL)
    prob0_ZIP_train <- t(prob0_ZIP_train); prob0_ZIP_test <- t(prob0_ZIP_test)
   
    
    #ZIPLDA
    zipldacv.out2 <- ZIPDA.cv(dat$x,dat$y,prob0=prob0_ZIP_train,type="mle")
    ZIPLDA.out2 <- ZIPLDA(dat$x,dat$y,dat$xte,transform=FALSE,prob0=prob0_ZIP_test,rho=zipldacv.out2$bestrho,type="mle") 
    ZIPLDAaa2 = length(which(ZIPLDA.out2$ytehat!=dat$yte))
    countZIP[j] = ZIPLDAaa2 / length(dat$yte)
  
  
    #PLDA
    cv.out <- PLDA.cv(dat$x,dat$y)
    out <- PLDA(dat$x,dat$y,dat$xte,transform=FALSE,rho=cv.out$bestrho) 
    aa = length(which(out$ytehat!=dat$yte))
    countP[j] = aa / length(dat$yte)
    
    
    #calculate the method-of-moment estimation on dispersions;
    X = t(dat$x)
    tt = getT(X,sizeFactors=rep(1,ncol(X)))$target  
    rM = rowMeans(X);
    rV = rowVars(X);
    disp = (rV-rM) / rM^2
    disp[which(disp<0)] <- 0
    disperhat = getAdjustDisp(disp,shrinkTarget=tt)$adj
    
    #NBLDA
    NBLDAcv.out2 <- NBLDA.cv(dat$x,dat$y,phihat=disperhat)
    NBLDA.out2 <- NBLDA(dat$x,dat$y,dat$xte,phihat=disperhat,rho=NBLDAcv.out2$bestrho) 
    bb = length(which(NBLDA.out2$ytehat!=dat$yte))
    countNB[j] = bb / length(dat$yte)
    
    #estimated mean mu
    ns <- NullModel(dat$x, type="mle")$n
    uniq <- sort(unique(dat$y))
    ds <-  GetDn(ns,dat$x,dat$y,beta=1)
    mu <- mean(ns * ds[dat$y,])
    #estimated  dispersion parameter based on the ZINBLDA model
    disperhatmean <- mean(disperhat)
    par <- t(apply(dat$x,2,fun,mu=mu,disperhatmean=disperhatmean))
    disperhatZINB <- par[,2]

    #Estimate the probability of zero in the data based on ZINBLDA model
    prob0_ZINB_train <- estimatepZINB(dat$x,dat$y,dat$x,dat$y,phihat=disperhatZINB,beta=1,type="mle", prior=NULL)
    prob0_ZINB_test <- estimatepZINB(dat$x,dat$y,dat$xte,dat$yte,phihat=disperhatZINB,beta=1,type="mle", prior=NULL)
    prob0_ZINB_train <- t(prob0_ZINB_train); prob0_ZINB_test <- t(prob0_ZINB_test)
    
    #ZINBLDA
    zinbldacv.out2 <- ZINB.cv(dat$x,dat$y,prob0=prob0_ZINB_train,phihat=disperhatZINB)
    ZInbLDA.out2 <- ZINBLDA(dat$x,dat$y,dat$xte,phihat=disperhatZINB,prob0=prob0_ZINB_test,rho=zinbldacv.out2$bestrho) 
    ZInbLDAaa2 = length(which(ZInbLDA.out2$ytehat!=dat$yte))
    countZINB[j] = ZInbLDAaa2/length(dat$yte)
    
  }
  rbind(countP,countNB,countZIP,countZINB)
}


#Parallel computing step
com.myfun <- cmpfun(myfun)
cl <- makeCluster(4)
registerDoParallel(cl)
ptm <- proc.time()
x <- foreach(j=1:simtime,.combine='rbind',.export = c("countP","countNB","countZIP","countZINB"),
             .packages = c("e1071", "class","sSeq","limma","edgeR","PoiClaClu")) %dopar% myfun(i)

stopCluster(cl)


#All simulation results of the four models are averaged respectively
PLDA = colMeans(x[seq(1,nrow(x),4),]) 
NBLDA = colMeans(x[seq(2,nrow(x),4),]) 
ZIPLDA = colMeans(x[seq(3,nrow(x),4),]) 
ZINBLDA = colMeans(x[seq(4,nrow(x),4),]) 

PLDA
NBLDA
ZIPLDA
ZINBLDA



lwd <- 2
plot( variable,PLDA,type="b",lty=1,lwd=lwd,col = "black",ylim=c(0,0.3),xlab="n",ylab="Misclassification rate",main=paste('pzero=',pzero,' dispersion=',dispersion,' DE=',DE))
lines( variable,NBLDA,type="b",lty=2,lwd=lwd,col = "green") 
lines(variable, ZIPLDA,type="b",lty=3,lwd=lwd,col = "blue") 
lines(variable, ZINBLDA,type="b",lty=4,lwd=lwd,col = "red") 
legend("topright", c("PLDA","NBLDA","ZIPLDA","ZINBLDA"),
       col=c("black","green","blue","red"),
       lty=c(1,2,3,4), lwd=1, cex=0.7)


proc.time() - ptm  