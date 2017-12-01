#Library imports

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

#The initial theta

theta=c(4.63,115,75,3,3.66,11.8)

#n.rep and n.burnin are arbitrary and could be lowered if needed

n.rep <- 1000000
n.burnin <- 5000

#Remove exploratory plots and prepare screen for trace plots

close.screen(all = TRUE)
split.screen(c(3,2))

#Calculate the loglikelihood of initial theta

ll0 <- ll(theta,trial["r"],trial["d"],trial["bio1"])+pri(theta)

#Vectors to store the alpha and number of accepts

all.alpha <- rep(0,n.rep)
all.accept <- rep(0,n.rep)

#matrix to store theta proposals
th <- matrix(0,6,n.rep)
th[,1] <- theta

#the main loop

for (i in 2:n.rep) {
  #set proposal
  th[1,i] <- th[1,i-1]+rnorm(1,0,1.9)
  th[2,i] <- th[2,i-1]+rnorm(1,0,2.4)
  th[3,i] <- th[3,i-1]+rnorm(1,0,1.5)
  th[4,i] <- th[4,i-1]+rnorm(1,0,0.04)
  th[5,i] <- th[5,i-1]+rnorm(1,0,0.04)
  th[6,i] <- th[6,i-1]+rnorm(1,0,1.4)
  

  ## get ll proposal and calculate alpha
  
  ll1 <- ll(th[,i],trial["r"],trial["d"],trial["bio1"])+pri(th[,i])
  all.alpha[i] <- exp(ll1-ll0)
  
  if (runif(1) < exp(ll1-ll0)) { ## accept
    all.accept[i] <- 1
    ll0 <- ll1 ## Keep ll0 in sync with th
  } else { ## reject
    th[,i] <- th[,i-1]
    ll1 <- ll0
  }

  
  
  if (i %% 5000 == 0) { ## Do some nice plotting...
    
    #Trace plots
    
    screen(1)
     plot(th[1,1:i],type="l",
         ylab=(paste("theta_",1,sep="")))
     
   screen(2)
   plot(th[2,1:i],type="l",
        ylab=(paste("theta_",2,sep="")))
   
   screen(3)
   plot(th[3,1:i],type="l",
        ylab=(paste("theta_",3,sep="")))
   
   screen(4)
   plot(th[4,1:i],type="l",
        ylab=(paste("theta_",4,sep="")))
   
   screen(5)
   plot(th[5,1:i],type="l",
        ylab=(paste("theta_",5,sep="")))
   
   screen(6)
   plot(th[6,1:i],type="l",
        ylab=(paste("theta_",6,sep="")))
  }
}



