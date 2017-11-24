library(truncnorm)
library(mvtnorm)

#functions

ll <- function(theta,r,d,x)
{
  E0 <- theta[1]
  Emax <- theta[2]
  ED50 <- theta[3]
  lambda <- theta[4]
  sigmasq <- exp(theta[5])
  beta1 <- theta[6]
  n <- length(r)
  
  if((ED50 + beta1) > 0)
  {
    mu <- E0 + ((d^lambda)*Emax)/((d^lambda)+(ED50+beta1*x)^lambda)
    out <-  (-n/2)*(log(2*pi)+log(sigmasq))-(1/(2*sigmasq))*sum((r-mu)^2)
  } else{
    out <- -Inf
  }
  return(out)
}

pri <- function(theta)
{
  E0 <- theta[1]
  Emax <- theta[2]
  ED50 <- theta[3]
  lambda <- theta[4]
  sigmasq <- exp(theta[5])
  beta1 <- theta[6]  
  
  E0pri <- dnorm(E0,0,10)
  Emaxpri <- dnorm(Emax,100,10)  
  pED50 <- ED50/200
  pED50pri <- dbeta(pED50,2.5,5)
  lambdapri <- dunif(lambda,0.5,3)
  sigmasqpri <- dtruncnorm(sigmasq,a=0,mean=3,sd=5)
  betapri <- dnorm(beta1,10,4)
  
  ind <- ((ED50 + beta1)>0)
  
  return(log(E0pri*Emaxpri*pED50pri*lambdapri*sigmasqpri*betapri*ind))
}

#read in data
test<-read.csv("CW18.csv")
#test<-test[test$bio1==1,]


#Metropolis Hastings
set.seed(123)
theta=c(4.63,114,75,3,4,20)
plotting <- 3 ## Which theta[i] to plot
n.rep <- 35000
prop.stdDev<-c(4.63,114,75,3,0.0125,0.005)
#n.burnin <- 5000
close.screen(all = TRUE)
split.screen(c(3,1))

y <- 1:60
ll0 <- ll(theta,test["r"],test["d"],test["bio1"])+pri(theta)
all.alpha <- rep(0,n.rep)
all.accept <- rep(0,n.rep)
th <- matrix(0,6,n.rep)
th[,1] <- theta
for (i in 2:n.rep) {
  
  #th[,i]<-th[,i-1]+rmvnorm(1,mean=c(1,0,0,0,0,0))*prop.stdDev

  ## proposal...
  th[1,i] <- th[1,i-1]+rt(1,df=10)*5
  ## comment out next line for no time trend...
  th[2,i] <- th[2,i-1]+rt(1,df=2)*5.2
  th[3,i] <- th[3,i-1]+rt(1,df=5)*2
  th[4,i] <- th[4,i-1]+rt(1,df=2)*0.0115
  th[5,i] <- th[5,i-1]+rt(1,df=5)*0.0125
  th[6,i] <- th[6,i-1]+rt(1,df=3)*0.0105
  
  ## get ll proposal
  ll1 <- ll(th[,i],test["r"],test["d"],test["bio1"])+pri(th[,i])
  acc <- min(1,exp(ll1-ll0))
  all.alpha[i] <- acc
  if (runif(1) <= acc) { ## accept
    all.accept[i] <- 1
    ll0 <- ll1 ## Keep ll0 in sync with th
  } else { ## reject
    th[,i] <- th[,i-1]
  }
    
  if (i %% 35000 == 0) { ## Do some nice plotting...
#    op=par(mfrow=c(2,2))
#    tmp = acf(th[plotting,1:i],plot=FALSE,lag.max=1)
#    acf(th[plotting,1:i],
#        main=paste("Lag-one correlation = ",signif(tmp$acf[2],4),sep=""))
#    hist(all.alpha[1:i],
#         main=paste("Mean alpha = ",signif(mean(all.alpha[1:i]),4),sep=""),
#         xlab="alpha")
#    plot(th[plotting,1:i],type="l",
#         main=paste("Mean theta_",plotting," = ",
#                    signif(mean(th[plotting,1:i]),4)),
#         ylab=(paste("theta_",plotting,sep="")))
#    inter = quantile(th[plotting,1:i],c(0.025,0.975))
#    hist(th[plotting,1:i],
#         main=(paste("Credible Int. = (",signif(inter[1],4),", ",
#                     signif(inter[2],4),")",sep="")),
#         xlab=paste("theta_",plotting,sep=""))
#    par(op)
	    screen(1)
     plot(th[4,1:i],type="l",
         ylab=(paste("theta_",4,sep="")))

	    screen(2)     
     plot(th[5,1:i],type="l",
          ylab=(paste("theta_",5,sep="")))

	screen(3)     
     plot(th[6,1:i],type="l",
          ylab=(paste("theta_",6,sep="")))
  }
}

mean(all.accept)
mean(all.alpha)

th.mod<-th[,-(1:5000)]
V<-cov(t(th.mod))
mu<-rowMeans(th.mod)
cor.test(th.mod[5,],th.mod[6,],method="spearman")

##Adjusted MH including correlation between the parameters
# V is the variance-covariance matrix of the parameters
# mu is the vector containing the mean of each parameter
theta=c(4.63,114,75,3,4,20)
ll0 <- ll(theta,test["r"],test["d"],test["bio1"])+pri(theta)
n.rep <- 35000
all.alpha <- rep(0,n.rep)
all.accept <- rep(0,n.rep)
th.cor <- matrix(0,6,n.rep)
th.cor[,1] <- theta
for (i in 2:n.rep){
	th.cor[,i]<-th.cor[,i-1]+rmvnorm(1,sigma=V)*(2.4/sqrt(6))^2
	ll1 <- ll(th.cor[,i],test["r"],test["d"],test["bio1"])+pri(th.cor[,i])
  	acc <- min(1,exp(ll1-ll0))
  	all.alpha[i] <- acc
	if (runif(1) <= acc) { ## accept
    		all.accept[i] <- 1
    		ll0 <- ll1 ## Keep ll0 in sync with th
  	} else { ## reject
    		th.cor[,i] <- th.cor[,i-1]
  }

}

par(mfrow=c(3,2))
for (i in 1:6){
plot(th.cor[i,],type='l',ylab=(paste("theta_",i,sep="")))
}
