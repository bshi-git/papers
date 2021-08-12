###############################################################
## R codes using simlation data for calculating Brier scores 
## and AUC values
## created by Bin Shi
## Reference: Shi, Wei and Huang (2020) Statistics in Medicine
###############################################################

rm(list = ls())
library(survival)
library(rockchalk)
library(fpca)
library(msm)
library(pec)
library(tdROC)

Q = 500
brier.lm_only = rep(NA,Q)
brier.lm_FPCA = rep(NA,Q)
auc.lm_only = rep(NA,Q)
auc.lm_FPCA = rep(NA,Q)

setwd("C:\\Users\\binsh\\Desktop\\majorrevising\\submissioncodes\\dataset_gen")

for(q in 1:Q) {  
  set.seed(q)
  
  ####survival time simulation####
  ##parameters for the survival model
  n=1669
  b=c(-0.27,0.30,-0.48)
  S = c(1043.3, 702.3,  234.1)
  genSurvData <- function(beta = b,
                          r = 0,
                          id.study = NA){
    #Scale parameter
    lambda <- exp(-7)
    #Shape parameter
    nue <-1.5   
    n <- n
    p <- length(beta)
    beta <- matrix(beta, ncol = 1)
    mu <- rep(0, p)
    X <- mvrnorm(n, mu, diag(S*40))
    
    #Calculate survival times
    T <- (-log(runif(n)) / (lambda * exp(X %*% beta)))^(1/nue)
    T = round(T + 35, 0)# add 35 to move min age to 35
    #observations are censored by outside of unif(n,45,75)
    Ctimes <- round(runif(n, 45, 75), 0)
    Time <- pmin(T, Ctimes)
    event <- as.numeric(T <= Ctimes)
    
    #Reorder data frame: T, event, covariates
    dat <- data.frame(T = Time, X, event = event)
    dat$id.study <- id.study <- seq(1,n)    
    tmp.names <- names(dat)
    colnames(dat) = c("T", "gamma1", "gamma2", "gamma3","event","id.study" )
    #Returning a matrix speeds-up things
    dat <- as.matrix(dat)
    return(dat)
  }
  
  dat.sim <- genSurvData()
  dat.sim <- as.data.frame(dat.sim)
  
  ####generate measure-time grids####
  mu_age=c(34.5, 42.7, 47.1, 50.3, 54.0, 58.1, 60.8)
  sigma_age =c(81.9, 80.9, 82.1, 80.9, 80.7, 79.2, 78.4)
  exam_ages = matrix(NA, n, 7)
  lage = round(rtnorm(n, mu_age[1], sigma_age[1],15,65), 0)
  exam_ages[,1] = lage
  preage = lage
  for (j in c(2:7)){           
    for (i in c(1:n)) {
      postage = round(rtnorm(1, mu_age[j], sigma_age[j],15,85), 0)
      while (preage[i] > postage) {
        postage = round(rtnorm(1, mu_age[j], sigma_age[j],15,85), 0)
      }
      exam_ages[i,j] = postage
      
    }
    
    
    preage = exam_ages[,j]
  }
  
  ####generate longitudinal biomarker values####
  sigma_wn = 0.1
  sigma = 0.1
  load("tc_mu.RData")
  load("tc_FPC.RData")
  load("phi.RData")
  tc_eigenFUN = t(phi)
  X = tc_FPC %*% tc_eigenFUN
  Y <- tc_mu + X + matrix(rnorm(length(X), 0, sigma_wn), nrow=nrow(X)) # white noise
  
  ##produce simulated dataset
  dat_sim <- matrix(,length(dat.sim$T),18)
  dat_sim[,1] = dat.sim$id.study
  dat_sim[,2] = dat.sim$event
  dat_sim[,3] = lage
  dat_sim[,4] = dat.sim$T
  dat_sim[,5:11] = round(Y[,1:7],1)
  dat_sim[,12:18] = round(exam_ages,0)
  dat_sim = as.data.frame(dat_sim)
  colnames(dat_sim) = c("shareid","chd","lage","chddate","tc1","tc2","tc3","tc4","tc5","tc6","tc7","ageAtexam1","ageAtexam2","ageAtexam3","ageAtexam4","ageAtexam5","ageAtexam6","ageAtexam7")
  ##remove observations which chddate <= lage
  index<-which(dat_sim$chddate>dat_sim$lage)
  dat_sim = dat_sim[index,]
  
  
  #####################LM only models######################       
  ##reshape wide to long format
  w<-dat_sim
  l1<-reshape(w,varying = list(c("tc1","tc2","tc3","tc4","tc5","tc6","tc7"),c("ageAtexam1","ageAtexam2","ageAtexam3","ageAtexam4","ageAtexam5","ageAtexam6","ageAtexam7")),v.names =c("tc_val","examage_val"),direction = "long") 
  l1.sort <- l1[order(l1$shareid),]
  
  
  ##landmark analysis   
  coef<-rep(NA,7)
  upper<-rep(NA,7)
  lower<-rep(NA,7)
  time<-c(35,40,45,50,55,60,65)
  k=1
  for(j in time){           
    chd1 = dat_sim
    ages1<-chd1[,12:18]
    tcs1<-chd1[,5:11]
    tcs1[,8]=NA
    
    #to get the most recent TC before chd (LVCF)
    for (i in 1:dim(w)[1]) {
      
      index <- which(ages1[i,]<=j)
      #if(length(as.vector(index)) !=0){ last <- max(index)}
      if (length(as.vector(index)) ==0) {
        tcs1[i,8]=NA
      } else {
        last <- max(index)
        tcs1[i,8]=tcs1[i,last]
      }
    }        
    
    chdbeta1<-cbind(chd1[,1:4],tcs1[,8])
    colnames(chdbeta1)[5]<-"tc"    
    dataname = paste0("beta", as.character(j))
    coxname = paste0("cox", as.character(j))
    # exam age before j and no CHD event no loss of followup  before j
    dat1<-subset(chdbeta1,lage<=j & !(chd==1 & chddate <j) & (chddate >= j) )        
    dat1[,"predtime"]<- j
    dat1[,"pt"]<- (dat1$predtime-25)/100
    dat1[,"T"] = dat1$chddate-j        
    assign(dataname, dat1)
    
    try(cox<-coxph(Surv(chddate, chd) ~ tc, dat1), silent = TRUE)      
    assign(coxname, cox) 
    coef[k]<-cox$coef
    upper[k]<-confint(cox)[1]
    lower[k]<-confint(cox)[2]
    k=k+1
  } 
  
  ##construct super model from super datasets
  data1<-rbind(beta35,beta40,beta45,beta50,beta55,beta60,beta65)
  data1[data1$lage==data1$chddate, "chddate"]<- data1[data1$lage==data1$chddate,"chddate"]+1 #if lage==chddate then chddate +1     
  
  try(cox1<-coxph(Surv(T, chd) ~ tc + I(tc*log(pt)) + strata(predtime), data1, x= TRUE, ties = "breslow", singular.ok=TRUE, robust=TRUE), silent = TRUE) 
  
  ##### brier score of Survival data     
  dat1 = na.omit(data1)
  try(Models <- list("Cox.X1"=coxph(Surv(T, chd) ~ tc + I(tc*log(pt)) + strata(predtime),data=dat1,x=TRUE,y=TRUE)), silent = TRUE)
  # compute the apparent prediction error
  try(PredError <- pec(object=Models,
                       formula=Surv(T, chd) ~ tc + I(tc*log(pt)) + strata(predtime),		     
                       data=dat1,
                       exact=TRUE,
                       cens.model="marginal",
                       splitMethod="none",
                       B=0,
                       verbose=TRUE), silent = TRUE)       
  
  try(fm <- tdROC( X = predict(cox1,type="lp"), Y = dat1$T, delta = dat1$chd, tau = 5, span = 0.1, nboot = 0, alpha = 0.05, n.grid = 1000, cut.off = 5:9 ), silent = TRUE)
  
  auc.lm_only[q]=fm$AUC$value
  brier.lm_only[q] = summary(PredError, times =5)$AppErr[4]
  
  ##################FPCA-LM################################    
  # to perform FPCA for total cholestrol
  #reshape wide to long format
  w<-dat_sim
  l<-reshape(w,varying = list(c("tc1","tc2","tc3","tc4","tc5","tc6","tc7"),c("ageAtexam1","ageAtexam2","ageAtexam3","ageAtexam4","ageAtexam5","ageAtexam6","ageAtexam7")),v.names =c("tc_val","examage_val"),direction = "long") 
  l.sort <- l[order(l$shareid),]
  
  tc<-l.sort[,1:7]
  fda=tc[tc$examage_val <= 75 & tc$examage_val >= 25,]
  file<-fda[,c(1,7,6)]
  #to standardize the measurement and to make time interval(0---1)
  datf=na.omit(file)
  datf$tc_val=(datf$tc_val-mean(datf$tc_val))/sd(datf$tc_val)
  datf$examage_val=datf$examage_val/75.5
  dt=cbind(datf$shareid,datf$tc_val,datf$examage_val)
  colnames(dt)=c("ID","measurement","time")
  ##parameters for fpca
  M.set<-c(4,5,6,7)
  r.set<-3   
  ini.method="EM"
  basis.method="bs"
  sl.v=rep(0.5,10)
  max.step=50
  grid.l=seq(0,1,0.01)
  grids=c(seq(0,1,0.002))    
  
  result<-fpca.mle(dt, M.set,r.set,ini.method, basis.method,sl.v,max.step,grid.l,grids)  
  ##rescaled grid
  M<-result$selected_model[1] 
  r<-result$selected_model[2]
  evalest<-result$eigenvalues ## estimated  
  sig2est<-result$error_var ## estimated 
  eigenfest<-result$eigenfunctions 
  muest<-result$fitted_mean
  
  ##Landmark analysis for FPCA acores
  range_up=matrix(rep(NA,7*3),ncol=3)
  range_low=matrix(rep(NA,7*3),ncol=3)
  sigma=matrix(rep(NA,7*3),ncol=3)
  
  coef<-matrix(rep(NA,7*3),ncol=3)
  upper<-matrix(rep(NA,7*3),ncol=3)
  lower<-matrix(rep(NA,7*3),ncol=3)   
  k=1  
  for(j in c(40,45,50,55,60,65,70)){      
    dat2<-subset(fda,examage_val<=j & !(chd==1 & chddate <j) & (chddate >= j))#exam age before j and no CHD event or no loss of followup before j(= large or equal to j)
    dat2$examage_val=dat2$examage_val/75.5
    ##rescaled grid
    grids.new<-result$grid[result$grid >= min(dat2$examage_val,na.rm=T) & result$grid <= max(dat2$examage_val,na.rm=T)] 
    gridname = paste0("grid", as.character(j))
    assign(gridname, result$grid)
    ## derive fpc scores and look at the predicted curve---integral fpc scores with different time interval 
    fpcs<-fpca.score(dt,grids.new,muest,evalest,eigenfest,sig2est,r) 
    #get predicted trajectories on the observed measurement points
    fpca3<-cbind(unique(na.omit(fda$shareid)),fpcs)
    colnames(fpca3)<-c("shareid","gamma1","gamma2","gamma3")
    chddata<-merge(unique(dat2[,1:4]),fpca3, by="shareid")               
    dataname = paste0("dat", as.character(j))        
    
    #perform cox analysis
    try(cox<-coxph(Surv(chddate, chd) ~ gamma1+gamma2+gamma3, chddata), silent = TRUE) 
    
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
  
  data.bind<-rbind(data40,data45,data50,data55,data60,data65,data70)
  try(cox1<-coxph(Surv(T, chd) ~ gamma1+ I(gamma1*log(pt))+gamma2+ I(gamma2*log(pt))+gamma3+ I(gamma3*log(pt))+ strata(predtime), data = data.bind, ties = "breslow", singular.ok=TRUE, robust=TRUE), silent = TRUE)   
  
  ##### brier score of Survival data    
  data = na.omit(data.bind)
  try(Models <- list("Cox.X2"=coxph(Surv(T, chd) ~ gamma1+ I(gamma1*log(pt))+gamma2+ I(gamma2*log(pt))+gamma3+ I(gamma3*log(pt))+ strata(predtime),data=data,x=TRUE,y=TRUE)), silent = TRUE)
  # compute the apparent prediction error
  try(PredError <- pec(object=Models,
                       formula=Surv(T, chd) ~ gamma1+ I(gamma1*log(pt))+gamma2+ I(gamma2*log(pt))+gamma3+ I(gamma3*log(pt))+ strata(predtime),
                       data=data,
                       exact=TRUE,
                       cens.model="marginal",
                       splitMethod="none",
                       B=0,
                       verbose=TRUE), silent = TRUE)   
  
  try(fm <- tdROC( X = predict(cox1,type="lp"), Y = data$T, delta = data$chd, tau = 5, span = 0.1, nboot = 0, alpha = 0.05, n.grid = 1000, cut.off = 5:9 ), silent = TRUE)
  auc.lm_FPCA[q]=fm$AUC$value
  brier.lm_FPCA[q] = summary(PredError, times =5)$AppErr[4]
  #return(cbind(brier.lm_only, brier.lm_FPCA, auc.lm_only, auc.lm_FPCA))
  
}#end of forloop

res = cbind(brier.lm_only, brier.lm_FPCA, auc.lm_only, auc.lm_FPCA)
res
#save(res, file= "res.RData")
apply(res,2, mean,na.rm = TRUE)
apply(res,2, sd,na.rm = TRUE)
apply(res,2,quantile,c(0.025,0.975),na.rm = TRUE)


