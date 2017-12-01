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
# trial <-trial[trial$bio1==1,]

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

theta=c(4.63,115,75,3,3.66,11.8)
n.rep <- 1000000
#The burn in might change for each variable in theta. Look at this later
#n.burnin <- 5000
close.screen(all = TRUE)
split.screen(c(3,2))


ll0 <- ll(theta,trial["r"],trial["d"],trial["bio1"])+pri(theta)
all.alpha <- rep(0,n.rep)
all.accept <- rep(0,n.rep)
avetheta1 <- rep(0,n.rep)
avetheta2 <- rep(0,n.rep)
avetheta3 <- rep(0,n.rep)
avetheta4 <- rep(0,n.rep)
avetheta5 <- rep(0,n.rep)
avetheta6 <- rep(0,n.rep)
th <- matrix(0,6,n.rep)
th[,1] <- theta
for (i in 2:n.rep) {
  ## proposal...
  th[1,i] <- th[1,i-1]+rnorm(1,0,1.9)
  ## comment out next line for no time trend...
  th[2,i] <- th[2,i-1]+rnorm(1,0,2.4)
  th[3,i] <- th[3,i-1]+rnorm(1,0,1.5)
  th[4,i] <- th[4,i-1]+rnorm(1,0,0.04)
  th[5,i] <- th[5,i-1]+rnorm(1,0,0.04)
  th[6,i] <- th[6,i-1]+rnorm(1,0,1.4)
  
  if(th[4,i]<0){
    th[,i] <- th[,i-1]
  } else {
    ## get ll proposal
    ll1 <- ll(th[,i],trial["r"],trial["d"],trial["bio1"])+pri(th[,i])
    all.alpha[i] <- exp(ll1-ll0)
    avetheta1[i]=(avetheta1[i-1]*(i-1)+th[1,i])/i
    avetheta2[i]=(avetheta2[i-1]*(i-1)+th[2,i])/i
    avetheta3[i]=(avetheta3[i-1]*(i-1)+th[3,i])/i
    avetheta4[i]=(avetheta4[i-1]*(i-1)+th[4,i])/i
    avetheta5[i]=(avetheta5[i-1]*(i-1)+th[5,i])/i
    avetheta6[i]=(avetheta6[i-1]*(i-1)+th[6,i])/i
    if (runif(1) < exp(ll1-ll0)) { ## accept
      all.accept[i] <- 1
      ll0 <- ll1 ## Keep ll0 in sync with th
    } else { ## reject
      th[,i] <- th[,i-1]
      ll1 <- ll0
    }
  }
  
  
  if (i %% 5000 == 0) { ## Do some nice plotting...
    
    screen(1)
     plot(th[1,1:i],type="l",
         ylab=(paste("theta_",1,sep="")))
    #  
    #  screen(2)     
    #  plot(avetheta1[1:i],type="l",
    #       ylab=(paste("average_theta_",1,sep="")))
    #  
    #  screen(3)  
    #  inter = quantile(th[1,1:i],c(0.025,0.975))
    #  hist(th[1,1:i],
    #       main=(paste("Credible Int. = (",signif(inter[1],4),", ",
    #                   signif(inter[2],4),")",sep="")),
    #       xlab=paste("theta_",1,sep=""))
    
    screen(2)
    plot(th[2,1:i],type="l",
         ylab=(paste("theta_",2,sep="")))
    
    # screen(3)     
    # plot(avetheta2[1:i],type="l",
    #      ylab=(paste("average_theta_",2,sep="")))
    # 
    # screen(5)  
    # inter = quantile(th[2,1:i],c(0.025,0.975))
    # hist(th[2,1:i],
    #      main=(paste("Credible Int. = (",signif(inter[1],4),", ",
    #                  signif(inter[2],4),")",sep="")),
    #      xlab=paste("theta_",2,sep=""))
    
    screen(3)
    plot(th[3,1:i],type="l",
         ylab=(paste("theta_",3,sep="")))
    # 
    # screen(4)
    # plot(avetheta3[1:i],type="l",
    #      ylab=(paste("average_theta_",3,sep="")))
    # 
    # screen(6)
    # inter = quantile(th[3,1:i],c(0.025,0.975))
    # hist(th[3,1:i],
    #      main=(paste("Credible Int. = (",signif(inter[1],4),", ",
    #                  signif(inter[2],4),")",sep="")),
    #      xlab=paste("theta_",3,sep=""))


    screen(4)
    plot(th[4,1:i],type="l",
         ylab=(paste("theta_",4,sep="")))
    # 
    # screen(4)
    # plot(avetheta4[1:i],type="l",
    #      ylab=(paste("average_theta_",4,sep="")))
    # 
    # screen(6)
    # inter = quantile(th[4,1:i],c(0.025,0.975))
    # hist(th[4,1:i],
    #      main=(paste("Credible Int. = (",signif(inter[1],4),", ",
    #                  signif(inter[2],4),")",sep="")),
    #      xlab=paste("theta_",4,sep=""))
    
    screen(5)
    plot(th[5,1:i],type="l",
         ylab=(paste("theta_",5,sep="")))
    # 
    # screen(4)
    # plot(avetheta5[1:i],type="l",
    #      ylab=(paste("average_theta_",5,sep="")))
    # 
    # screen(6)
    # inter = quantile(th[5,1:i],c(0.025,0.975))
    # hist(th[5,1:i],
    #      main=(paste("Credible Int. = (",signif(inter[1],4),", ",
    #                  signif(inter[2],4),")",sep="")),
    #      xlab=paste("theta_",5,sep=""))
    
    screen(6)
    plot(th[6,1:i],type="l",
         ylab=(paste("theta_",6,sep="")))
    
    # screen(4)
    # plot(avetheta6[1:i],type="l",
    #      ylab=(paste("average_theta_",6,sep="")))
    # 
    # screen(6)
    # inter = quantile(th[6,1:i],c(0.025,0.975))
    # hist(th[6,1:i],
    #      main=(paste("Credible Int. = (",signif(inter[1],4),", ",
    #                  signif(inter[2],4),")",sep="")),
    #      xlab=paste("theta_",6,sep=""))
  }
}



