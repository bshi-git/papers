###############################################################
## R codes for dataset generation based on the three FPC score   
## estimation and betas from survival analysis of FHS data
## created by Bin Shi
## Reference: Shi, Wei and Huang (2020) Statistics in Medicine
##############################################################

rm(list = ls())
library(survival)
library(fdapace)
library(mvtnorm)
library(MASS)
library(msm)

setwd("C:\\Users\\binsh\\Desktop\\majorrevising\\submissioncodes\\dataset_gen")

####survival time simulation####
##parameters for the survival model
set.seed(8828)
n=1669
b=c(-0.27,0.30,-0.48)
S = diag(c(1043.3, 702.3,  234.1))
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
  X <- mvrnorm(n, mu, S)
  
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
#1- sum(dat.sim$event)/n
#hist(dat.sim$T)
#head(dat.sim) 

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
#head(exam_ages)

####generate longitudinal biomarker values####
sigma_wn = 0.1
sigma = 0.1

#tc_mu = fpcaObjFlies$mu
#tc_mu = c(180.6225, 181.1963, 181.5945, 182.0185, 182.6073, 183.4407, 184.5404, 185.8607, 187.2865, 188.6668, 189.8807, 190.9002, 191.7989, 192.6953, 193.6790, 194.7723, 195.9366, 197.1140, 198.2738, 199.4212, 200.5699, 201.7229,202.8793, 204.0453, 205.2255, 206.4027, 207.5301, 208.5406, 209.3709, 209.9872, 210.3961, 210.6290, 210.7201,210.7016, 210.6141, 210.5076, 210.4264, 210.3924, 210.4018, 210.4349, 210.4609, 210.4369, 210.3073, 210.0117,209.5057, 208.7885, 207.9264, 207.0581, 206.3758, 206.0776, 206.3068)
#save(tc_mu, file = "tc_mu.RData")
load("tc_mu.RData")
#tc_FPC = fpcaObjFlies$xiEst[,1:3]
#save(tc_FPC, file = "tc_FPC.RData")
load("tc_FPC.RData")
#phi = fpcaObjFlies$phi[,1:3]
#save(phi, file = "phi.RData")
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
#head(dat_sim)
#  shareid chd lage chddate   tc1   tc2   tc3   tc4   tc5   tc6   tc7 ageAtexam1 ageAtexam2
#1       1   0   45      68 243.5 254.9 232.6 205.3 217.4 216.5 202.4         45         46
#3       3   0   49      63 122.9 165.8 190.8 200.2 226.8 228.3 208.9         49         56
#4       4   0   38      70 207.4 223.1 215.9 197.1 206.7 205.0 189.1         38         82
#5       5   1   17      57 221.5 234.9 220.4 194.9 201.3 201.3 191.6         17         61
#6       6   0   59      62 114.2 169.6 190.8 192.8 193.5 232.9 233.5         59         79
#7       7   0   34      51 208.8 227.8 224.3 207.0 192.1 219.8 213.1         34         52
#dim(dat_sim)
#[1] 1394   18

write.table(dat_sim, "dat_sim.csv")

####fig_lable function####
fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 300/100
  sh <- strheight(text, cex=cex) * 300/100
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}