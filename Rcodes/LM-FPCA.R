###############################################################
## R codes for Landmark analysis with FPCA scores (LM-FPCA) 
## using simlation data  
## created by Bin Shi
## Reference: Shi, Wei and Huang (2020) Statistics in Medicine
##############################################################

rm(list = ls())
library(survival)
library(fdapace)
library(mvtnorm)
library(MASS)
library(fpca)

setwd("C:\\Users\\binsh\\Desktop\\majorrevising\\submissioncodes")
source("dataset_sim.R")

###perform FPCA for total cholestrol###
##reshape wide to long format
exam = dat_sim
w<-exam[order(exam$shareid),]
l<-reshape(w,varying = list(c("tc1","tc2","tc3","tc4","tc5","tc6","tc7"),c("ageAtexam1","ageAtexam2","ageAtexam3","ageAtexam4","ageAtexam5","ageAtexam6","ageAtexam7")),v.names =c("tc_val","examage_val"),direction = "long") 
l.sort <- l[order(l$shareid),]
tc<-l.sort[,1:7]
fda=tc[tc$examage_val <= 75 & tc$examage_val >= 25,]
#dim(fda)
#[1] 11319     7
file<-fda[,c(1,7,6)]

##to standardize the measurement and to make time interval(0---1)
dat=na.omit(file)
dat$tc_val=(dat$tc_val-mean(dat$tc_val))/sd(dat$tc_val)
dat$examage_val=dat$examage_val/75.5
dat=cbind(dat$shareid,dat$tc_val,dat$examage_val)
colnames(dat)=c("ID","measurement","time")

##parameters for fpca
M.set<-c(4,5,6,7) # candidate model number
r.set<- 3 # eigen value number
ini.method="EM"
basis.method="bs"
sl.v=rep(0.5,10)
max.step=50
grid.l=seq(0,1,0.01)
grids=c(seq(0,1,0.002))

result<-fpca.mle(dat, M.set,r.set,ini.method, basis.method,sl.v,max.step,grid.l,grids)
#summary(result)  
#result$eigenvalues
#dim(result$eigenfunctions)

##rescaled grid
M<-result$selected_model[1] 
r<-result$selected_model[2]
evalest<-result$eigenvalues ## estimated  
sig2est<-result$error_var ## estimated 
eigenfest<-result$eigenfunctions 
muest<-result$fitted_mean 
##look at the CV scores and convergence
#result$cv_scores #CV 
#result$converge #convergence

###Landmark analysis for FPCA scores###
range_up=matrix(rep(NA,7*3),ncol=3)
range_low=matrix(rep(NA,7*3),ncol=3)
sigma=matrix(rep(NA,7*3),ncol=3)

coef<-matrix(rep(NA,7*3),ncol=3)
upper<-matrix(rep(NA,7*3),ncol=3)
lower<-matrix(rep(NA,7*3),ncol=3)
time<-c(40,45,50,55,60,65,70)
k=1 # count number of landmark datasets

for(j in c(40,45,50,55,60,65,70)){
  #exam age before j and no CHD event or no loss of followup before j(= large or equal to j)
  dat2<-subset(fda,examage_val<=j & !(chd==1 & chddate <j) & (chddate >= j))
  dat2$examage_val=dat2$examage_val/75.5
  ##rescaled grid
  grids.new<-result$grid[result$grid >= min(dat2$examage_val,na.rm=T) & result$grid <= max(dat2$examage_val,na.rm=T)] 
  gridname = paste0("grid", as.character(j))
  assign(gridname, result$grid)
  ##derive fpc scores and look at the predicted curve 
  fpcs<-fpca.score(dat,grids.new,muest,evalest,eigenfest,sig2est,r) 
  #get predicted trajectories on the observed measurement points
  fpca3<-cbind(unique(na.omit(fda$shareid)),fpcs)
  colnames(fpca3)<-c("shareid","gamma1","gamma2","gamma3")
  chddata<-merge(unique(dat2[,1:4]),fpca3, by="shareid")
  dataname = paste0("dat", as.character(j))
  
  cox<-coxph(Surv(lage,chddate, chd) ~ gamma1+gamma2+gamma3, chddata)
  sigma[k,]=apply(chddata[,5:7],2,sd)
  
  ##store beta1-3 for separate models
  coef[k,]<-cox$coef[1:3]
  upper[k,]<-confint(cox)[,1][1:3]
  lower[k,]<-confint(cox)[,2][1:3]
  #get the pt value
  chddata[,"predtime"]<- j
  chddata[,"pt"]<- (chddata$predtime-25)/100
  chddata[,"T"]<- chddata$chddate - j
  
  dataname = paste0("data", as.character(j))
  assign(dataname, chddata)
  coxname = paste0("cox", as.character(j))
  assign(coxname, cox)
  
  k=k+1
}


coef_s=coef*sigma
upper_s=upper*sigma
lower_s=lower*sigma

##plot three eigenfunctions
plot(result$grid*75.5,result$eigenfunctions[1,], type = "l", col = "red", lwd = 3,lty=1,xlab = "Exam age (years)", ylab="First  three eigen components",ylim=c(-5,5))#,cex.lab=1.5)
points(result$grid*75.5,result$eigenfunctions[2,], type = "l", col = "blue", lwd = 3,lty=2)
points(result$grid*75.5,result$eigenfunctions[3,], type = "l", col = "green", lwd = 3,lty=3)
abline(h=0,lty=4)
legend("topleft", c("eigen1", "eigen2", "eigen3"), lty = c(1, 2, 3), pch = 1, col = c("red","blue","green"),title = "legend")


##super models
data<-rbind(data40,data45,data50,data55,data60,data65,data70)
cox1<-coxph(Surv(T, chd) ~ gamma1+ I(gamma1*log(pt))+gamma2+ I(gamma2*log(pt))+gamma3+ I(gamma3*log(pt)) + strata(predtime), data)
sigma=apply(sigma,2,mean)

##plot separate models and super models together
layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE))
par(mar=c(4,4,2,1.5),cex=0.8)

##eigen I
##plot separate model for eigen I
sapply(LETTERS[1], function(x) { 
  plot(time,coef_s[,1], type = "l", col = "red", lwd = 3,lty=1,xlab = "Exam age (years)", ylab=expression(beta[1](L)), ylim=c(-3.1, 2))#,cex.lab=1.5)
  points(time,upper_s[,1], type = "l", col = "blue", lwd = 1,lty=2)
  points(time,lower_s[,1], type = "l", col = "blue", lwd = 1,lty=2)
  abline(h=0,lty=4)
  fig_label(x, cex=2) 
})

##plot smooth curves eigen I
coef1<-rep(NA,100)
coef2<-rep(NA,100)
se<-rep(NA,100)
ages<-seq(40,70,length.out=100)
i=1
for(k in ages){
  coef1[i]<-(cox1$coef[1]+cox1$coef[2]*log((k-25)/100))*sigma[1]
  se[i]<-sqrt(cox1$var[1,1]+cox1$var[2,2]*(log((k-15)/100))^2 + 2*log((k-15)/100)*cox1$var[1,2])
  i=i+1
}
sapply(LETTERS[2], function(x) { 
  plot(ages,coef1, type = "l", col = "red", lwd = 1,xlab = "Exam age (years)", ylab = NA, ylim=c(-3.1, 2))#,cex.lab=1.5)
  points(ages, coef1 - 1.96*se, type = "l", col = "blue", lwd = 1,lty=2)
  points(ages, coef1 + 1.96*se, type = "l", col = "blue", lwd = 1,lty=2)
  abline(h=0,lty=4)
  fig_label(x, cex=2) 
})


#eigen II
##plot separate model for eigen II
sapply(LETTERS[3], function(x) { 
  plot(time,coef_s[,2], type = "l", col = "red", lwd = 3,xlab = "Exam age (years)", ylab=expression(beta[2](L)), ylim=c(-3.1, 2))#,cex.lab=1.5)
  points(time,upper_s[,2], type = "l", col = "blue", lwd = 1,lty=2)
  points(time,lower_s[,2], type = "l", col = "blue", lwd = 1,lty=2)
  abline(h=0,lty=4)
  fig_label(x, cex=2) 
})


##plot smooth curves for eigen II
coef3<-rep(NA,100)
coef4<-rep(NA,100)
se<-rep(NA,100)
ages<-seq(40,70,length.out=100)
i=1
for(k in ages){
  coef3[i]<-(cox1$coef[3]+cox1$coef[4]*log((k-25)/100))*sigma[2]
  se[i]<-sqrt(cox1$var[1,1]+cox1$var[2,2]*(log((k-15)/100))^2 + 2*log((k-15)/100)*cox1$var[1,2])
  i=i+1
}
sapply(LETTERS[4], function(x) { 
  plot(ages,coef3, type = "l", col = "red", lwd = 1,xlab = "Exam age (years)", ylab= NA, ylim=c(-3.1, 2))#,cex.lab=1.5)
  points(ages, coef3 - 1.96*se, type = "l", col = "blue", lwd = 1,lty=2)
  points(ages, coef3 + 1.96*se, type = "l", col = "blue", lwd = 1,lty=2)
  abline(h=0,lty=4)
  fig_label(x, cex=2) 
})


#eigen III
##plot separate model for eigen III
sapply(LETTERS[5], function(x) { 
  plot(time,coef_s[,3], type = "l", col = "red", lwd = 3,lty=1,xlab = "Exam age (years)", ylab=expression(beta[3](L)), ylim=c(-3.1, 2))#,cex.lab=1.5)
  points(time,upper_s[,3], type = "l", col = "blue", lwd = 1,lty=2)
  points(time,lower_s[,3], type = "l", col = "blue", lwd = 1,lty=2)
  abline(h=0,lty=4)
  fig_label(x, cex=2) 
})


##plot smooth curves for eigen III
coef5<-rep(NA,100)
coef6<-rep(NA,100)
se<-rep(NA,100)
ages<-seq(40,70,length.out=100)
i=1
for(k in ages){
  coef5[i]<-(cox1$coef[5]+cox1$coef[6]*log((k-25)/100))*sigma[3]
  se[i]<-sqrt(cox1$var[1,1]+cox1$var[2,2]*(log((k-15)/100))^2 + 2*log((k-15)/100)*cox1$var[1,2])
  i=i+1
}
sapply(LETTERS[6], function(x) { 
  plot(ages,coef5, type = "l", col = "red", lwd = 1,xlab = "Exam age (years)", ylab= NA, ylim=c(-3.1, 2))
  points(ages, coef5 - 1.96*se, type = "l", col = "blue", lwd = 1,lty=2)
  points(ages, coef5 + 1.96*se, type = "l", col = "blue", lwd = 1,lty=2)
  abline(h=0,lty=4)
  fig_label(x, cex=2) 
})

layout(matrix(c(1,1,1,1,1,1), 3, 2, byrow = TRUE))
#dev.off()