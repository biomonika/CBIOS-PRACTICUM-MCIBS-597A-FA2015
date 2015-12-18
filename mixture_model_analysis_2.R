library(ggplot2)


ipd_p1=read.table(file="ipds_P8.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)

###Obtain IPD with 8X coverage
ipd8=ipd_p1[which(ipd_p1$coverage==8),]
ipd8_A=ipd8[which(ipd8$base=="A"),]
ipd8_T=ipd8[which(ipd8$base=="T"),]
ipd8_G=ipd8[which(ipd8$base=="G"),]
ipd8_C=ipd8[which(ipd8$base=="C"),]

###FUNCTIONS FOR PLOTS

library(mixtools)


plot.normal.components <- function(mixture,component.number,...) {
  curve(mixture$lambda[component.number] *
          dnorm(x,mean=mixture$mu[component.number],
                sd=mixture$sigma[component.number]), add=TRUE, ...)
}

plot.gamma.components <- function(mixture,component.number,...) {
  curve(mixture$lambda[component.number] *
          dgamma(x,shape=mixture$gamma.pars[1,][component.number],
                 scale=mixture$gamma.pars[2,][component.number]), add=TRUE, ...)
}


pnormmix <- function(x,mixture) {
  lambda <- mixture$lambda
  k <- length(lambda)
  pnorm.from.mix <- function(x,component) {
    lambda[component]*pnorm(x,mean=mixture$mu[component],
                            sd=mixture$sigma[component])
  }
  pnorms <- sapply(1:k,pnorm.from.mix,x=x)
  return(rowSums(pnorms))
}

pgammamix <- function(x,mixture) {
  lambda <- mixture$lambda
  k <- length(lambda)
  pgamma.from.mix <- function(x,component) {
    lambda[component]*pgamma(x,shape=mixture$gamma.pars[1,][component],
                             scale=mixture$gamma.pars[2,][component])
  }
  pgammas <- sapply(1:k,pgamma.from.mix,x=x)
  return(rowSums(pgammas))
}



dnormalmix <- function(x,mixture,log=FALSE) {
  lambda <- mixture$lambda
  k <- length(lambda)
  # Calculate share of likelihood for all data for one component
  like.component <- function(x,component) {
    lambda[component]*dnorm(x,mean=mixture$mu[component],
                            sd=mixture$sigma[component])
  }
  # Create array with likelihood shares from all components over all data
  likes <- sapply(1:k,like.component,x=x)
  # Add up contributions from components
  d <- rowSums(likes)
  if (log) {
    d <- log(d)
  }
  return(d)
}

loglike.normalmix <- function(x,mixture) {
  loglike <- dnormalmix(x,mixture,log=TRUE)
  return(sum(loglike))
}



########ANALYSIS

###NORMAL MIXTURE MODEL


###2 Mixtures
### more than ones and take the max 
ipd8_tm_A.k2_100<-vector('list',100)
for(i in 1:100)
{
  ipd8_tm_A.k2_100[[i]]<- normalmixEM(ipd8_A$tMean,k=2,maxit=100,epsilon=0.01) 
}
max <- which.max(sapply(ipd8_tm_A.k2_100, function(x){x$loglik}))
min <- which.min(sapply(ipd8_tm_A.k2_100, function(x){x$loglik}))
ipd8_tm_A.k2_100[[max]]$loglik
ipd8_tm_A.k2_100[[min]]$loglik



windows()
plot(hist(ipd8_A$tMean,breaks=100),col="white",border="grey",freq=FALSE,xlab="tMean",main='MAXloglik(blue) and MINloglik(red)')
#lines(density(ipd8_A$tMean),lty=2)
sapply(1:2,plot.normal.components,mixture=ipd8_tm_A.k2_100[[max]], col='blue', lwd='2')
sapply(1:2,plot.normal.components,mixture=ipd8_tm_A.k2_100[[min]], col='red', lwd='2')

windows()
loglik2<-sapply(ipd8_tm_A.k2_100,function(x){x$loglik})
loglik2
plot(loglik2, col='blue', type='p', ylim = c(-155,-140))


##QQ plot 2 component max
windows()
distinct.ipd8_tm_A <- sort(unique(ipd8_A$tMean))
tcdfs <- pnormmix(distinct.ipd8_tm_A,mixture=ipd8_tm_A.k2_100[[max]])
ecdfs <- ecdf(ipd8_A$tMean)(distinct.ipd8_tm_A)
plot(tcdfs,ecdfs,col='blue',xlab="Theoretical CDF",ylab="Empirical CDF",xlim=c(0,1),
     ylim=c(0,1))
abline(0,1)

##QQ plot 2 component min
windows()
distinct.ipd8_tm_A <- sort(unique(ipd8_A$tMean))
tcdfs <- pnormmix(distinct.ipd8_tm_A,mixture=ipd8_tm_A.k2_100[[min]])
ecdfs <- ecdf(ipd8_A$tMean)(distinct.ipd8_tm_A)
plot(tcdfs,ecdfs,col='red',xlab="Theoretical CDF",ylab="Empirical CDF",xlim=c(0,1),
     ylim=c(0,1))
abline(0,1)





###3 mixtures
### more than ones and take the max 
ipd8_tm_A.k3_100<-vector('list',100)
for(i in 1:100)
{
  ipd8_tm_A.k3_100[[i]]<- normalmixEM(ipd8_A$tMean,k=3,maxit=100,epsilon=0.01) 
}
max <- which.max(sapply(ipd8_tm_A.k3_100, function(x){x$loglik}))
min <- which.min(sapply(ipd8_tm_A.k3_100, function(x){x$loglik}))
ipd8_tm_A.k3_100[[max]]$loglik
ipd8_tm_A.k3_100[[min]]$loglik

ipd8_tm_A.k3_100[[1]]

windows()
plot(hist(ipd8_A$tMean,breaks=100),col="white",border="grey",freq=FALSE,xlab="tMean",main="IPD 8x A")
#lines(density(ipd8_A$tMean),lty=2)
sapply(1:3,plot.normal.components,mixture=ipd8_tm_A.k3_100[[1]], col='red')
sapply(1:3,plot.normal.components,mixture=ipd8_tm_A.k3_100[[2]], col='blue')
sapply(1:3,plot.normal.components,mixture=ipd8_tm_A.k3_100[[3]], col='green')
sapply(1:3,plot.normal.components,mixture=ipd8_tm_A.k3_100[[4]], col='yellow')


windows()
plot(hist(ipd8_A$tMean,breaks=100),col="white",border="grey",freq=FALSE,xlab="tMean",main='MAXloglik(blue) and MINloglik(red)')
#lines(density(ipd8_A$tMean),lty=2)
sapply(1:3,plot.normal.components,mixture=ipd8_tm_A.k3_100[[max]], col='blue', lwd='2')
sapply(1:3,plot.normal.components,mixture=ipd8_tm_A.k3_100[[min]], col='red', lwd='2')

windows()
loglik3<-sapply(ipd8_tm_A.k3_100,function(x){x$loglik})
loglik3
plot(loglik3, col='blue', type='p', ylim = c(-155,-140))


##QQ plot 3 component max
windows()
distinct.ipd8_tm_A <- sort(unique(ipd8_A$tMean))
tcdfs <- pnormmix(distinct.ipd8_tm_A,mixture=ipd8_tm_A.k3_100[[max]])
ecdfs <- ecdf(ipd8_A$tMean)(distinct.ipd8_tm_A)
plot(tcdfs,ecdfs,col='blue',xlab="Theoretical CDF",ylab="Empirical CDF",xlim=c(0,1),
     ylim=c(0,1))
abline(0,1)

##QQ plot 3 component min
windows()
distinct.ipd8_tm_A <- sort(unique(ipd8_A$tMean))
tcdfs <- pnormmix(distinct.ipd8_tm_A,mixture=ipd8_tm_A.k3_100[[min]])
ecdfs <- ecdf(ipd8_A$tMean)(distinct.ipd8_tm_A)
plot(tcdfs,ecdfs,col='red',xlab="Theoretical CDF",ylab="Empirical CDF",xlim=c(0,1),
     ylim=c(0,1))
abline(0,1)



##### number of mixtures plot


ipd8_tm_A.k2_100<-vector('list',100)
ipd8_tm_A.k3_100<-vector('list',100)
ipd8_tm_A.k4_100<-vector('list',100)
ipd8_tm_A.k5_100<-vector('list',100)
ipd8_tm_A.k6_100<-vector('list',100)
ipd8_tm_A.k7_100<-vector('list',100)
ipd8_tm_A.k8_100<-vector('list',100)
ipd8_tm_A.k9_100<-vector('list',100)
ipd8_tm_A.k10_100<-vector('list',100)

                 
for(i in 1:100)
{
  ipd8_tm_A.k2_100[[i]]<- normalmixEM(ipd8_A$tMean,k=2,maxit=100,epsilon=0.01) 
  ipd8_tm_A.k3_100[[i]]<- normalmixEM(ipd8_A$tMean,k=3,maxit=100,epsilon=0.01) 
  ipd8_tm_A.k4_100[[i]]<- normalmixEM(ipd8_A$tMean,k=4,maxit=100,epsilon=0.01) 
  ipd8_tm_A.k5_100[[i]]<- normalmixEM(ipd8_A$tMean,k=5,maxit=100,epsilon=0.01) 
  ipd8_tm_A.k6_100[[i]]<- normalmixEM(ipd8_A$tMean,k=6,maxit=100,epsilon=0.01) 
  ipd8_tm_A.k7_100[[i]]<- normalmixEM(ipd8_A$tMean,k=7,maxit=100,epsilon=0.01) 
  ipd8_tm_A.k8_100[[i]]<- normalmixEM(ipd8_A$tMean,k=8,maxit=100,epsilon=0.01) 
  ipd8_tm_A.k9_100[[i]]<- normalmixEM(ipd8_A$tMean,k=9,maxit=100,epsilon=0.01) 
  ipd8_tm_A.k10_100[[i]]<- normalmixEM(ipd8_A$tMean,k=10,maxit=100,epsilon=0.01) 
}



lg2 <- rep(0,100)
lg3 <- rep(0,100)
lg4 <- rep(0,100)
lg5 <- rep(0,100)
lg6 <- rep(0,100)
lg7 <- rep(0,100)
lg8 <- rep(0,100)
lg9 <- rep(0,100)
lg10 <- rep(0,100)



for(i in 1:100)
{
  lg2[i]<-ipd8_tm_A.k2_100[[i]]$loglik  
  lg3[i]<-ipd8_tm_A.k3_100[[i]]$loglik
  lg4[i]<-ipd8_tm_A.k4_100[[i]]$loglik  
  lg5[i]<-ipd8_tm_A.k5_100[[i]]$loglik  
  lg6[i]<-ipd8_tm_A.k6_100[[i]]$loglik  
  lg7[i]<-ipd8_tm_A.k7_100[[i]]$loglik  
  lg8[i]<-ipd8_tm_A.k8_100[[i]]$loglik  
  lg9[i]<-ipd8_tm_A.k9_100[[i]]$loglik  
  lg10[i]<-ipd8_tm_A.k10_100[[i]]$loglik
} 

bic2 <- -2*lg2 + 5*log(100) 
bic3 <- -2*lg3 + 8*log(100)
bic4 <- -2*lg4 + 11*log(100)
bic5 <- -2*lg5 + 14*log(100)
bic6 <- -2*lg6 + 17*log(100)
bic7 <- -2*lg7 + 20*log(100)
bic8 <- -2*lg8 + 23*log(100)
bic9 <- -2*lg9 + 26*log(100)
bic10 <- -2*lg10 + 29*log(100)

windows()
boxplot(lg2,lg3,lg4,lg5,lg6,lg7,lg8,lg9,lg10)

windows()
boxplot(bic2,bic3,bic4,bic5,bic6,bic7,bic8,bic9,bic10)


# for each data i compute the maximum of the posterior prob of each trial

max_2 <- matrix(0,177,100)
max_3 <- matrix(0,177,100)
max_4 <- matrix(0,177,100)
max_5 <- matrix(0,177,100)
max_6 <- matrix(0,177,100)
max_7 <- matrix(0,177,100)
max_8 <- matrix(0,177,100)
max_9 <- matrix(0,177,100)
max_10 <- matrix(0,177,100)

for(i in 1:177){
  for(j in 1:100){
    max_2[i,j] <- max(ipd8_tm_A.k2_100[[j]]$posterior[i,])
    max_3[i,j] <- max(ipd8_tm_A.k3_100[[j]]$posterior[i,])
    max_4[i,j] <- max(ipd8_tm_A.k4_100[[j]]$posterior[i,])
    max_5[i,j] <- max(ipd8_tm_A.k5_100[[j]]$posterior[i,])
    max_6[i,j] <- max(ipd8_tm_A.k6_100[[j]]$posterior[i,])
    max_7[i,j] <- max(ipd8_tm_A.k7_100[[j]]$posterior[i,])
    max_8[i,j] <- max(ipd8_tm_A.k8_100[[j]]$posterior[i,])
    max_9[i,j] <- max(ipd8_tm_A.k9_100[[j]]$posterior[i,])
    max_10[i,j] <- max(ipd8_tm_A.k10_100[[j]]$posterior[i,])
  }
}



# for each data, I compute the mean of the max posterior of each trial
mean_max_2 <- rep(0,177) 
mean_max_3 <- rep(0,177) 
mean_max_4 <- rep(0,177) 
mean_max_5 <- rep(0,177) 
mean_max_6 <- rep(0,177) 
mean_max_7 <- rep(0,177) 
mean_max_8 <- rep(0,177) 
mean_max_9 <- rep(0,177) 
mean_max_10 <- rep(0,177) 

for(i in 1:177){
    mean_max_2[i] <- mean(max_2[i,])   
    mean_max_3[i] <- mean(max_3[i,])
    mean_max_4[i] <- mean(max_4[i,])
    mean_max_5[i] <- mean(max_5[i,])
    mean_max_6[i] <- mean(max_6[i,])
    mean_max_7[i] <- mean(max_7[i,])
    mean_max_8[i] <- mean(max_8[i,])
    mean_max_9[i] <- mean(max_9[i,])
    mean_max_10[i] <- mean(max_10[i,])
}  

windows()
boxplot(mean_max_2, mean_max_3, mean_max_4, mean_max_5, mean_max_6,
        mean_max_7, mean_max_8, mean_max_9, mean_max_10 )

windows()
plot(mean_max_2)

windows()
plot(mean_max_3)



mean_loglik_2 <- mean(sapply(ipd8_tm_A.k2_100, function(x){x$loglik}))
mean_loglik_3 <- mean(sapply(ipd8_tm_A.k3_100, function(x){x$loglik}))
mean_loglik_4 <- mean(sapply(ipd8_tm_A.k4_100, function(x){x$loglik}))
mean_loglik_5 <- mean(sapply(ipd8_tm_A.k5_100, function(x){x$loglik}))
mean_loglik_6 <- mean(sapply(ipd8_tm_A.k6_100, function(x){x$loglik}))
mean_loglik_7 <- mean(sapply(ipd8_tm_A.k7_100, function(x){x$loglik}))
mean_loglik_8 <- mean(sapply(ipd8_tm_A.k8_100, function(x){x$loglik}))
mean_loglik_9 <- mean(sapply(ipd8_tm_A.k9_100, function(x){x$loglik}))
mean_loglik_10 <- mean(sapply(ipd8_tm_A.k10_100, function(x){x$loglik}))

mean_loglik_k <- c(mean_loglik_2,mean_loglik_3, mean_loglik_4, mean_loglik_5,
                   mean_loglik_6, mean_loglik_7, mean_loglik_8,
                   mean_loglik_9, mean_loglik_10)


par <- c(5,8,11,14,17,20,23,26,29)
bic <- -2*mean_loglik_k + par*log(100)
bic

windows()
plot(c(2,3,4,5,6,7,8,9,10),mean_loglik_k, pch = 19, col='red' )

windows()
plot(c(2,3,4,5,6,7,8,9,10),bic, pch = 19, col='blue')









################## GAMMA MIXTURE MODEL

###3 mixtures
### more than ones and take the max 
ipd8_tm_A.k3_100_g<-vector('list',100)
for(i in 1:100)
{
  ipd8_tm_A.k3_100_g[[i]]<- gammamixEM(ipd8_A$tMean,k=3,maxit=100,epsilon=0.01) 
}
max <- which.max(sapply(ipd8_tm_A.k3_100_g, function(x){x$loglik}))
min <- which.min(sapply(ipd8_tm_A.k3_100_g, function(x){x$loglik}))
ipd8_tm_A.k3_100_g[[max]]$loglik



windows()
plot(hist(ipd8_A$tMean,breaks=100),col="white",border="grey",freq=FALSE,xlab="tMean",main='MAXloglik(blue) and MINloglik(red)')
#lines(density(ipd8_A$tMean),lty=2)
sapply(1:3,plot.gamma.components,mixture=ipd8_tm_A.k3_100_g[[max]], col='blue', lwd='2')
sapply(1:3,plot.gamma.components,mixture=ipd8_tm_A.k3_100_g[[min]], col='red', lwd='2')

windows()
loglik3_g<-sapply(ipd8_tm_A.k3_100,function(x){x$loglik})
loglik3_g
plot(loglik3_g, col='blue', type='p')








###################### T,C,G #############à

###T

ipd8_tm_T.k3_100<-vector('list',100)
for(i in 1:100)
{
  ipd8_tm_T.k3_100[[i]]<- normalmixEM(ipd8_T$tMean,k=3,maxit=100,epsilon=0.01) 
}
max <- which.max(sapply(ipd8_tm_T.k3_100, function(x){x$loglik}))
ipd8_tm_T.k3_100[[max]]$loglik


windows()
plot(hist(ipd8_T$tMean,breaks=100),col="white",border="grey",freq=FALSE,xlab="tMean",main='MAXloglik mixture T')
#lines(density(ipd8_A$tMean),lty=2)
sapply(1:3,plot.normal.components,mixture=ipd8_tm_T.k3_100[[max]], col='blue', lwd='2')


###G
ipd8_tm_G.k3_100<-vector('list',100)
for(i in 1:100)
{
  ipd8_tm_G.k3_100[[i]]<- normalmixEM(ipd8_G$tMean,k=3,maxit=100,epsilon=0.01) 
}
max <- which.max(sapply(ipd8_tm_G.k3_100, function(x){x$loglik}))
ipd8_tm_G.k3_100[[max]]$loglik


windows()
plot(hist(ipd8_G$tMean,breaks=100),col="white",border="grey",freq=FALSE,xlab="tMean",main='MAXloglik mixture G')
#lines(density(ipd8_A$tMean),lty=2)
sapply(1:3,plot.normal.components,mixture=ipd8_tm_G.k3_100[[max]], col='blue', lwd='2')


##### C

ipd8_tm_C.k3_100<-vector('list',100)
for(i in 1:100)
{
  ipd8_tm_C.k3_100[[i]]<- normalmixEM(ipd8_C$tMean,k=3,maxit=100,epsilon=0.01) 
}
max <- which.max(sapply(ipd8_tm_C.k3_100, function(x){x$loglik}))
ipd8_tm_C.k3_100[[max]]$loglik


windows()
plot(hist(ipd8_C$tMean,breaks=100),col="white",border="grey",freq=FALSE,xlab="tMean",main='MAXloglik mixture C')
#lines(density(ipd8_A$tMean),lty=2)
sapply(1:3,plot.normal.components,mixture=ipd8_tm_C.k3_100[[max]], col='blue', lwd='2')

