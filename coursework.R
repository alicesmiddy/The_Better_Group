library("truncnorm", lib.loc="~/R/win-library/3.3")
library("mvtnorm", lib.loc="~/R/win-library/3.3")

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

trial <- read.csv("C:/Users/Alice Smiddy/Documents/Bath/Applied Statistical inference/coursework/CW18.csv")

#explatory analysis

close.screen(all = TRUE)
split.screen(c(2,1))
split.screen(c(1,2),1)
screen(3)
with(trial, plot(d,r))
screen(4)
with(trial, plot(age,r))
screen(2)
with(trial, hist(r, breaks=15))

#metropolis hastings

theta=c(0,100,50,0.5,3,10)
plotting <- 1 ## Which theta[i] to plot
n.rep <- 100000
#n.burnin <- 5000
close.screen(all = TRUE)
split.screen(c(3,1))

y <- 1:60
ll0 <- ll(theta,trial["r"],trial["d"],trial["bio1"])+pri(theta)
all.alpha <- rep(0,n.rep)
all.accept <- rep(0,n.rep)
th <- matrix(0,6,n.rep)
th[,1] <- theta
for (i in 2:n.rep) {
  ## proposal...
  th[1,i] <- th[1,i-1]+rnorm(1,mean=1)
  ## comment out next line for no time trend...
  th[2,i] <- th[2,i-1]+rnorm(1)*0.001
  th[3,i] <- th[3,i-1]+rnorm(1)*0.001
  th[4,i] <- th[4,i-1]+rnorm(1)*0.1
  th[5,i] <- th[5,i-1]+rnorm(1)*0.1
  th[6,i] <- th[6,i-1]+rnorm(1)*0.001
  
  ## get ll proposal
  ll1 <- ll(th[,i],trial["r"],trial["d"],trial["bio1"])+pri(th[,i])
  all.alpha[i] <- exp(ll1-ll0)
  if (runif(1) < exp(ll1-ll0)) { ## accept
    all.accept[i] <- 1
    ll0 <- ll1 ## Keep ll0 in sync with th
  } else { ## reject
    th[,i] <- th[,i-1]
    ll1 <- ll0
  }

  
  # 
  # #this bit needs to be looked into... May or may not have finished the metropolis hastings algorithm.
  # #plot the graph to see if you get as expected
  # #So far you've only looked at E0 (I think?). Double check and code up for the rest of theta if correct
  # 
  if (i %% 2000 == 0) { ## Do some nice plotting...
  #   op=par(mfrow=c(2,2))
  #   tmp = acf(th[plotting,1:i],plot=FALSE,lag.max=1)
  #   acf(th[plotting,1:i],
  #       main=paste("Lag-one correlation = ",signif(tmp$acf[2],4),sep=""))
  #   hist(all.alpha[1:i],
  #        main=paste("Mean alpha = ",signif(mean(all.alpha[1:i]),4),sep=""),
  #        xlab="alpha")

    screen(1)
     plot(th[plotting,1:i],type="l",
         ylab=(paste("theta_",plotting,sep="")))
     
     # plot(th[2,1:i],type="l",
     #      ylab=(paste("theta_",2,sep="")))
     
     # plot(th[3,1:i],type="l",
     #      main=paste("Mean theta_",3," = ",
     #                 signif(mean(th[3,1:i]),4)),
     #      ylab=(paste("theta_",3,sep="")))
     screen(2)     
     plot(th[4,1:i],type="l",
          ylab=(paste("theta_",4,sep="")))
     
     # plot(th[5,1:i],type="l",
     #      main=paste("Mean theta_",5," = ",
     #                 signif(mean(th[5,1:i]),4)),
     #      ylab=(paste("theta_",5,sep="")))
     screen(3)     
     plot(th[6,1:i],type="l",
          ylab=(paste("theta_",6,sep="")))
  #   inter = quantile(th[plotting,1:i],c(0.025,0.975))
  #   hist(th[plotting,1:i],
  #        main=(paste("Credible Int. = (",signif(inter[1],4),", ",
  #                    signif(inter[2],4),")",sep="")),
  #        xlab=paste("theta_",plotting,sep=""))
  #   par(op)
  }
}



