###############################################################
## R codes for Landmark analysis (LM-only) using simlation data  
## created by Bin Shi
## Reference: Shi, Wei and Huang (2020) Statistics in Medicine
##############################################################

rm(list = ls())
library(survival)
setwd("C:\\Users\\binsh\\Desktop\\majorrevising\\submissioncodes")
source("dataset_sim.R")

######landmark analysis######
library(survival)
coef<-rep(NA,7)
upper<-rep(NA,7)
lower<-rep(NA,7)
time<-c(35,40,45,50,55,60,65)
k=1 #count number of landmark datasets
for(j in c(35,40,45,50,55,60,65)){
  chd=dat_sim
  ages<-chd[,12:18]
  tcs<-chd[,5:11]
  tcs[,8]=NA
  
  #to get the most recent TC before chd
  for (i in 1:dim(chd)[1]) {
    index <- which(ages[i,]<=j)
    if (length(as.vector(index)) ==0) {
      tcs[i,8]=NA
    } else {
      last <- max(index)
      tcs[i,8]=tcs[i,last]
    }
  }
  
  ###separate model
  chdbeta<-cbind(chd[,c(1:4)],tcs[,8])
  colnames(chdbeta)[5]<-"tc"
  #head(chdbeta)
  dataname = paste0("beta", as.character(j))
  coxname = paste0("cox", as.character(j))
  
  #exam age before j and no CHD event no loss of followup  before j, which makes with landmark time increasing, dataset size is decreasing
  dat<-subset(chdbeta,lage<=j & !(chd==1 & chddate <j) & (chddate >= j) )
  dat[,"predtime"]<- j
  dat[,"pt"]<- (dat$predtime-25)/100
  dat[,"T"] = dat$chddate-j
  assign(dataname, dat)
  
  cox<-coxph(Surv(lage,chddate, chd) ~ tc, dat)
  assign(coxname, cox)
  coef[k]<-cox$coef[1]
  upper[k]<-confint(cox)[1]
  lower[k]<-confint(cox)[2]
  k=k+1
}

layout(matrix(c(1,2), 2, 1, byrow = TRUE))
##plot the separate model
sapply(LETTERS[1], function(x) { 
  plot(time,coef, type = "l", col = "red", lwd = 3,xlab = "Exam age (years)", ylab="Landmark cox regression coefficients",ylim=c(-0.015,0.03))#,cex.lab=1.5)
  points(time,upper, type = "l", col = "blue", lwd = 1,lty=2)
  points(time,lower, type = "l", col = "blue", lwd = 1,lty=2)
  abline(h=0,lty=4)
  fig_label(x, cex=2) 
})

###super model
ages<-chd[,12:18]
tcs<-chd[,5:11]
wchd<-cbind(chd[,1:4],tcs,ages)
exam<- wchd[order(wchd$chddate),]

#construct super datasets
data<-rbind(beta35,beta40,beta45,beta50,beta55,beta60,beta65)
data[data$lage==data$chddate, "chddate"]<- data[data$lage==data$chddate,"chddate"]+1#if lage==chddate then chddate +1

cox1<-coxph(Surv(T, chd) ~ tc + I(tc*log(pt)) + strata(predtime), data)
#summary(cox1)

##smoothing coeffiencts of super models by 100 points
coef_sm<-rep(NA,100)
se<-rep(NA,100)
ages<-seq(35,65,length.out =100)
i=1
for(k in ages){
  coef_sm[i]<-cox1$coef[1]+cox1$coef[2]*log((k-25)/100)
  se[i]<-sqrt(cox1$var[1,1]+cox1$var[2,2]*(log((k-15)/100))^2 + 2*log((k-15)/100)*cox1$var[1,2])
  i=i+1
}

##plot the super model
sapply(LETTERS[2], function(x) { 
  plot(ages,coef_sm, type = "l", col = "red", lwd = 1,xlab = "Exam age (years)", ylab="Cox regression TC coefficients",ylim=c(-0.015,0.03))#,cex.lab=1.5)
  points(ages, coef_sm - 1.96*se, type = "l", col = "blue", lwd = 1,lty=2)
  points(ages, coef_sm + 1.96*se, type = "l", col = "blue", lwd = 1,lty=2)
  abline(h=0,lty=4)
  fig_label(x, cex=2) 
})

layout(matrix(c(1,1), 2, 1, byrow = TRUE))
#dev.off()