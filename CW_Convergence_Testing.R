#The corrplot library is used to produce these nice correlation plots
# https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
library(corrplot)
library("truncnorm")
library("mvtnorm")
library(nortest)

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

trial <- read.csv("CW18.csv")
trial.biomarker<-subset(trial,bio1==1)
trial.no.biomarker<-subset(trial,bio1==0)
plot(trial$d,trial$r,xlab="Dose Amount",ylab="Response",main="Dose response for population with/without biomarker")
lines(trial.biomarker$d,trial.biomarker$r,col="red")
lines(trial.no.biomarker$d,trial.no.biomarker$r,col="blue")
##I haven't worked out how to make the legend work ¯\(°_o)/¯
legend(x=0,y=100,legend=c("Response","Biomarker","No Biomarker"),col=c("black","red",blue"))
hist(trial.biomarker$r,freq=FALSE,main="Histogram of response to treatment",xlab="Response")
hist(trial.no.biomarker$r,freq=FALSE,main="Histogram of response to treatment",xlab="Response")


theta=c(4.63,115,75,3,3.66,11.8)
n.rep <- 100000
#The burn in might change for each variable in theta. Look at this later
n.burnin <- 5000
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
mean(all.accept)
th<-th[,-(1:n.burnin)]

##After running the first version of Metropolis-Hastings
#Say our parameter matrix AFTER removing the burn-in period is called th
V<-cor(t(th)) #Correlation matrix of the parameters
VCOV<-cov(t(th))
corrplot(V,method="circle") #Plot the correlation matrix with a color bar
#We will see if any of the parameters is correlated

#To decide wether the correlation is significant or not, we perform the Pearson
# correlation test on each pair of parameters
#If the p-value is <5%, then the correlation is significant
cor.test(th[1,],th[2,],method="pearson")
cor.test(th[1,],th[3,],method="pearson")
cor.test(th[1,],th[4,],method="pearson")
cor.test(th[1,],th[5,],method="pearson")
cor.test(th[1,],th[6,],method="pearson")
cor.test(th[2,],th[3,],method="pearson")
cor.test(th[2,],th[4,],method="pearson")
cor.test(th[2,],th[5,],method="pearson")
cor.test(th[2,],th[6,],method="pearson")
cor.test(th[3,],th[4,],method="pearson")
cor.test(th[3,],th[5,],method="pearson")
cor.test(th[3,],th[6,],method="pearson")
cor.test(th[4,],th[5,],method="pearson")
cor.test(th[4,],th[6,],method="pearson")
cor.test(th[5,],th[6,],method="pearson")

# IF the correlation between the pairs is significant, we should try to modify
#the Metropolis-Hastings sampler by generating multivariate normal RV's
#centered on the previous VECTOR of parameters and with correlation matrix V
# i.e. in each iterate generate th.cor[,i]<-th.cor[,i-1]+rmvnorm(1,sigma=V)*std_dev
#Of course, the parameter std_dev should be tuned in order to get a ~25% acceptancece
#
theta=c(4.63,115,75,3,3.66,11.8)
n.rep <- 100000
#The burn in might change for each variable in theta. Look at this later
n.burnin <- 5000
close.screen(all = TRUE)
split.screen(c(3,2))
ll0 <- ll(theta,trial["r"],trial["d"],trial["bio1"])+pri(theta)
all.alpha <- rep(0,n.rep)
all.accept <- rep(0,n.rep)
th.cor <- matrix(0,6,n.rep)
th.cor[,1] <- theta
std_dev<-0.9 
for(i in 2:n.rep) {
	th.cor[,i]<-th.cor[,i-1]+rmvnorm(1,sigma=VCOV)*std_dev
	ll1 <- ll(th.cor[,i],trial["r"],trial["d"],trial["bio1"])+pri(th.cor[,i])
      all.alpha[i] <- exp(ll1-ll0)
	if(th.cor[4,i]<0){
    th.cor[,i] <- th.cor[,i-1]
  } else {
	if (runif(1) < exp(ll1-ll0)) { ## accept
      all.accept[i] <- 1
      ll0 <- ll1 ## Keep ll0 in sync with th
    } else { ## reject
      th.cor[,i] <- th.cor[,i-1]
      ll1 <- ll0
    }
}
	
	if (i %% 5000 == 0) { ## Do some nice plotting...
    close.screen(all = TRUE)
	split.screen(c(3,2))
    screen(1)
     plot(th.cor[1,1:i],type="l",
         ylab=(paste("theta_",1,sep="")))

    
    screen(2)
    plot(th.cor[2,1:i],type="l",
         ylab=(paste("theta_",2,sep="")))
    

    
    screen(3)
    plot(th.cor[3,1:i],type="l",
         ylab=(paste("theta_",3,sep="")))



    screen(4)
    plot(th.cor[4,1:i],type="l",
         ylab=(paste("theta_",4,sep="")))

    
    screen(5)
    plot(th.cor[5,1:i],type="l",
         ylab=(paste("theta_",5,sep="")))
       
    screen(6)
    plot(th.cor[6,1:i],type="l",
         ylab=(paste("theta_",6,sep="")))
    
     }

}
mean(all.accept)
th.cor<-th.cor[,-(1:n.burnin)]
#ACF lengths
#The acf length has to be computed for each paramter in order to draw independent
#subsamples

th1.ac<-acf(th[1,])[[1]][,,1]
th1.acl<-2*sum(th1.ac)-1

th2.ac<-acf(th[2,])[[1]][,,1]
th2.acl<-2*sum(th2.ac)-1

th3.ac<-acf(th[3,])[[1]][,,1]
th3.acl<-2*sum(th3.ac)-1

th4.ac<-acf(th[4,])[[1]][,,1]
th4.acl<-2*sum(th4.ac)-1

th5.ac<-acf(th[5,])[[1]][,,1]
th5.acl<-2*sum(th5.ac)-1

th6.ac<-acf(th[6,])[[1]][,,1]
th6.acl<-2*sum(th6.ac)-1

#Drawing independent subsamples for each parameter
th1.ind<-th[1,][seq(1,length(th[1,]),th1.acl)]
th2.ind<-th[2,][seq(1,length(th[2,]),th2.acl)]
th3.ind<-th[3,][seq(1,length(th[3,]),th3.acl)]
th4.ind<-th[4,][seq(1,length(th[4,]),th4.acl)]
th5.ind<-th[5,][seq(1,length(th[5,]),th5.acl)]
th6.ind<-th[6,][seq(1,length(th[6,]),th6.acl)]

#To decide wether each parameter converged or not, perform the two-sample
#Kolmogorov-Smirnov test by splitting into two each independent subsample of
#parameters
#We expect to observe a p-value GREATER than 0.05 to conclude that the parameter
#converged, i.e., each subsample has the same distribution}
#For simplicity, I'll just split each vector into two halves

#Theta1
H<-length(th1.ind)/2
th1.ind.A<-th1.ind[1:H]
th1.ind.B<-th1.ind[H+1:length(th1.ind)]
ks.test(th1.ind.A,th1.ind.B)

#Theta2
H<-length(th2.ind)/2
th2.ind.A<-th2.ind[1:H]
th2.ind.B<-th2.ind[H+1:length(th2.ind)]
ks.test(th2.ind.A,th2.ind.B)

#Theta3
H<-length(th3.ind)/2
th3.ind.A<-th3.ind[1:H]
th3.ind.B<-th3.ind[H+1:length(th3.ind)]
ks.test(th3.ind.A,th3.ind.B)

#Theta4
H<-length(th4.ind)/2
th4.ind.A<-th4.ind[1:H]
th4.ind.B<-th4.ind[H+1:length(th4.ind)]
ks.test(th4.ind.A,th4.ind.B)

#Theta5
H<-length(th5.ind)/2
th5.ind.A<-th5.ind[1:H]
th5.ind.B<-th5.ind[H+1:length(th5.ind)]
ks.test(th5.ind.A,th5.ind.B)

#Theta6
H<-length(th6.ind)/2
th6.ind.A<-th6.ind[1:H]
th6.ind.B<-th6.ind[H+1:length(th6.ind)]
ks.test(th6.ind.A,th6.ind.B)

#We perform the same procedure on the samples from the second Metropolis-Hastings
#loop, i.e., using the multivariate normal distribution AFTER discarding a certain
#burn-in period
#Recall that I named the output of this M-H sampling: th.cor
#My notation might become a bit cumbersome ¯\_(?)_/¯ ¯\(°_o)/¯

#ACF Lengths
th1.cor.ac<-acf(th.cor[1,])[[1]][,,1]
th1.cor.acl<-2*sum(th1.cor.ac)-1

th2.cor.ac<-acf(th.cor[2,])[[1]][,,1]
th2.cor.acl<-2*sum(th2.cor.ac)-1

th3.cor.ac<-acf(th.cor[3,])[[1]][,,1]
th3.cor.acl<-2*sum(th3.cor.ac)-1

th4.cor.ac<-acf(th.cor[4,])[[1]][,,1]
th4.cor.acl<-2*sum(th4.cor.ac)-1

th5.cor.ac<-acf(th.cor[5,])[[1]][,,1]
th5.cor.acl<-2*sum(th5.cor.ac)-1

th6.cor.ac<-acf(th.cor[6,])[[1]][,,1]
th6.cor.acl<-2*sum(th6.cor.ac)-1

#Drawing independent subsamples for each parameter
th1.cor.ind<-th.cor[1,][seq(1,length(th.cor[1,]),th1.cor.acl)]
th2.cor.ind<-th.cor[2,][seq(1,length(th.cor[2,]),th2.cor.acl)]
th3.cor.ind<-th.cor[3,][seq(1,length(th.cor[3,]),th3.cor.acl)]
th4.cor.ind<-th.cor[4,][seq(1,length(th.cor[4,]),th4.cor.acl)]
th5.cor.ind<-th.cor[5,][seq(1,length(th.cor[5,]),th5.cor.acl)]
th6.cor.ind<-th.cor[6,][seq(1,length(th.cor[6,]),th6.cor.acl)]

#Testing for convergence
#Theta1 w/correlation
H<-length(th1.cor.ind)/2
th1.cor.ind.A<-th1.cor.ind[1:H]
th1.cor.ind.B<-th1.cor.ind[H+1:length(th1.cor.ind)]
ks.test(th1.cor.ind.A,th1.cor.ind.B)

#Theta2 w/correlation
H<-length(th2.cor.ind)/2
th2.cor.ind.A<-th2.cor.ind[1:H]
th2.cor.ind.B<-th2.cor.ind[H+1:length(th2.cor.ind)]
ks.test(th2.cor.ind.A,th2.cor.ind.B)

#Theta3 w/correlation
H<-length(th3.cor.ind)/2
th3.cor.ind.A<-th3.cor.ind[1:H]
th3.cor.ind.B<-th3.cor.ind[H+1:length(th3.cor.ind)]
ks.test(th3.cor.ind.A,th3.cor.ind.B)

#Theta4 w/correlation
H<-length(th4.cor.ind)/2
th4.cor.ind.A<-th4.cor.ind[1:H]
th4.cor.ind.B<-th4.cor.ind[H+1:length(th4.cor.ind)]
ks.test(th4.ind.A,th4.ind.B)

#Theta5 w/correlation
H<-length(th5.cor.ind)/2
th5.cor.ind.A<-th5.cor.ind[1:H]
th5.cor.ind.B<-th5.cor.ind[H+1:length(th5.cor.ind)]
ks.test(th5.cor.ind.A,th5.cor.ind.B)

#Theta6 w/correlation
H<-length(th6.cor.ind)/2
th6.cor.ind.A<-th6.cor.ind[1:H]
th6.cor.ind.B<-th6.cor.ind[H+1:length(th6.cor.ind)]
ks.test(th6.cor.ind.A,th6.cor.ind.B)

##ASSUMING both th and th.cor converged nicely, i.e., every parameter converged
#How do we decide which one is a better model?
#Calculate the DIC for both estimates. Keep the one with lowest DIC
#Check Core Statistic 6.3, but it does not explain well how the DIC is obtained
#Also check slides 54-56 from https://www.ukdataservice.ac.uk/media/307220/presentation4.pdf

##DIC for the First M-H Approach
th.bar<-rowMeans(th)
D.th.bar<- (-2)*ll(th.bar,trial["r"],trial["d"],trial["bio1"])
D.bar<-0
K<-length(th[1,])
for (i in 1:K){
D.bar=D.bar+ (-2)*ll(th[,i],trial["r"],trial["d"],trial["bio1"])
}
D.bar=D.bar/K
DIC.th<-2*D.bar-D.th.bar

##DIC for the Second M-H Approach
#Again, sorry for the cumbersome notation ¯\_(?)_/¯ ¯\(°_o)/¯
th.cor.bar<-rowMeans(th.cor)
D.th.cor.bar<- (-2)*ll(th.cor.bar,trial["r"],trial["d"],trial["bio1"])
D.cor.bar<-0
K.cor<-length(th.cor[1,])
for (i in 1:K.cor){
D.cor.bar=D.cor.bar+ (-2)*ll(th.cor[,i],trial["r"],trial["d"],trial["bio1"])
}
D.cor.bar=D.cor.bar/K.cor
DIC.th.cor<-2*D.cor.bar-D.th.cor.bar

#Comparing the DIC for each model
DIC.th
DIC.th.cor

#We keep the model that includes the correlation between parameters
#Now, we calculate Credible Intervals for each parameter

#CI for E0
quantile(th.cor[1,],c(0.025,0.975))

#CI for Emax
quantile(th.cor[2,],c(0.025,0.975))

#CI for ED50
quantile(th.cor[3,],c(0.025,0.975))

#CI for lambda
quantile(th.cor[4,],c(0.025,0.975))

#CI for sigma
quantile(th.cor[5,],c(0.025,0.975))

#CI for Beta
quantile(th.cor[6,],c(0.025,0.975))

##Testing for normality
##qqnorm plots may provide visual evidence that the posterior distribution
#for each parameter is normal
par(mfrow=c(3,2))
qqnorm(th.cor[1,])
qqnorm(th.cor[2,])
qqnorm(th.cor[3,]) 
qqnorm(th.cor[4,]) #based on this, lambda might not have a normal posterior distribution
qqnorm(th.cor[5,])
qqnorm(th.cor[6,])
#However, the rest of the parameters MIGHT be normally distributed
#To confirm my suspicion, perform the Cramer von Mises test. (Install the nortest package)
#We expect to observe p-values greater than 5%, as the Null Hypothesis states that
#the data is normally distributed
cvm.test(th.cor[1,])
cvm.test(th.cor[2,])
cvm.test(th.cor[3,])
cvm.test(th.cor[4,])
cvm.test(th.cor[5,])
cvm.test(th.cor[6,])
#After performing the tests, neither of the parameters turned out to be normally distributed
#Perhaps this has to do with some behaviour around the tails of the posterior distribution

#Density plots
plot(density(th.cor[1,]),main="Kernel density for E0")
plot(density(th.cor[2,]),main="Kernel density for Emax")
plot(density(th.cor[3,]),main="Kernel density for ED50")
plot(density(th.cor[4,]),main="Kernel density for lambda")
plot(density(th.cor[5,]),main="Kernel density for sigma^2")
plot(density(th.cor[6,]),main="Kernel density for Beta")

##Final parameter estimates
E0.hat<-mean(th.cor[1,])
Emax.hat<-mean(th.cor[2,])
ED50.hat<-mean(th.cor[3,])
lambda.hat<-median(th.cor[4,])
sigmasq.hat<-mean(th.cor[5,])
beta.hat<-mean(th.cor[6,])

#Predicted response
r.hat<-E0.hat+((trial$d^(lambda.hat))*Emax.hat)/((trial$d^(lambda.hat))+(ED50.hat+beta.hat*trial$bio1)^(lambda.hat))
r.bio.hat<-E0.hat+((trial.biomarker$d^(lambda.hat))*Emax.hat)/((trial.biomarker$d^(lambda.hat))+(ED50.hat+beta.hat*trial.biomarker$bio1)^(lambda.hat))
r.no.bio.hat<-E0.hat+((trial.no.biomarker$d^(lambda.hat))*Emax.hat)/((trial.no.biomarker$d^(lambda.hat))+(ED50.hat+beta.hat*trial.no.biomarker$bio1)^(lambda.hat))

#Comparing observed vs Predicted
plot(trial$d,trial$r,xlab="Dose Amount",ylab="Response",main="Predicted Dose Response for population with/without biomarker")
lines(trial.biomarker$d,r.bio.hat,col="red")
lines(trial.no.biomarker$d,r.no.bio.hat,col="blue")
#Again, we should add legends showing which line is which
#The red line is the predicted response for the population WITH the biomarker
#The blue line is the predicted response for the population WITHOUT the biomarker

#Residuals

res<-trial$r-r.hat
par(mfrow=c(1,3))
plot(res)
qqnorm(res)
acf(res)
#Testing for normal distribution in residuals
#We expect to observe a p-value greater than 5%
cvm.test(res)
res.bio<-trial.biomarker$r-r.bio.hat
res.no.bio<-trial.no.biomarker$r-r.no.bio.hat
#Residual scatter plots
par(mfrow=c(3,1))
plot(res)
plot(res.bio)
plot(res.no.bio)


#Calculation of ED30
ED30<-ED50.hat*(0.3/0.7)^(1/lambda.hat)
