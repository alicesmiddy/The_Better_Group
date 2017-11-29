#The corrplot library is used to produce these nice correlation plots
# https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
library(corrplot)

##After running the first version of Metropolis-Hastings
#Say our parameter matrix AFTER removing the burn-in period is called th
V<-cor(t(th)) #Correlation matrix of the parameters
corrplot(M,method="circle") #Plot the correlation matrix with a color bar
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
th4.ind.B<-th5.ind[H+1:length(th5.ind)]
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
th1.cor.ind<-th[1,][seq(1,length(th.cor[1,]),th1.cor.acl)]
th2.cor.ind<-th[2,][seq(1,length(th.cor[2,]),th2.cor.acl)]
th3.cor.ind<-th[3,][seq(1,length(th.cor[3,]),th3.cor.acl)]
th4.cor.ind<-th[4,][seq(1,length(th.cor[4,]),th4.cor.acl)]
th5.cor.ind<-th[5,][seq(1,length(th.cor[5,]),th5.cor.acl)]
th6.cor.ind<-th[6,][seq(1,length(th.cor[6,]),th6.cor.acl)]

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
th4.cor.ind.B<-th5.cor.ind[H+1:length(th5.cor.ind)]
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
D.th.bar<- (-2)*ll(th.bar,r,d,x)
D.bar<-0
K<-length(th[1,])
for (i in 1:K){
D.bar=D.bar+ (-2)*ll(th[,i],r,d,x)
}
D.bar=D.bar/K
DIC.th<-2*D.bar-D.th.bar

##DIC for the Second M-H Approach
#Again, sorry for the cumbersome notation ¯\_(?)_/¯ ¯\(°_o)/¯
th.cor.bar<-rowMeans(th.cor)
D.th.cor.bar<- (-2)*ll(th.cor.bar,r,d,x)
D.cor.bar<-0
K.cor<-length(th.cor[1,])
for (i in 1:K.cor){
D.cor.bar=D.cor.bar+ (-2)*ll(th.cor[,i],r,d,x)
}
D.cor.bar=D.cor.bar/K.cor
DIC.th.cor<-2*D.cor.bar-D.th.cor.bar
