#packages
library(NADA)
library(glmnet)
library(hdnom)
#library(bujar)
library(parcor)
library(mboost)
library(gbm)
library(elasticnet)
library(mpath)

#Download useless function
source("Bujar/bstfit_PS.R")
source("Bujar/bujar_PS.R")
source("Bujar/predval_PS.R")
source("Bujar/bjboost.R")
source("Bujar/cvmboost.R")
source("Bujar/fpartial.R")
source("Bujar/gcv_enet.R")
source("Bujar/important_inter.R")
source("Bujar/interactions.R")
source("Bujar/mypartialPlot.R")
source("Bujar/rif.R")
source("Bujar/vimpint.R")
source("GaussBJ/cvGaussBJ.R")
source("GaussBJ/GaussBJ.R")

#Data
data(Recon)
head(Recon)

#Description
#We add 100 noise variables
n <- dim(Recon)[1]
p <- 100
Xnoise <- matrix(rnorm(n*p, 0, 1), ncol = p)
Xtrue <- Recon[,1:8]
X <- as.matrix(cbind(Xtrue, Xnoise))
Y <- Recon$AtraConc
#the detection limit is 0.05 ug/L to have 11% of left-censoring data
#Y[which(Recon$AtraConc <= 0.05)] <- 0.05
Cens <- Y
Cens[Y==0.05] <- 0
Cens[Y>0.05] <- 1
s <- rep(0, n) 
s[Cens == 0] <- 0.05 #Threshold variable

#Return the variable to use in survival analysis
M = max(Y) + 1 
Z = M - Y

#Fold CV
foldid <- rep(0, n)
foldid[which(Cens == 0)] <- sample(length(which(Cens == 0))) %% 5 + 1
foldid[which(Cens == 1)] <- sample(length(which(Cens == 1))) %% 5 + 1

#Fit method

#Lasso LOD
fit_LOD = glmnet::cv.glmnet(X, Y, nfolds = 5, foldid = foldid, alpha=1)

#Lasso NonParaBJ <- modified version from bujar package
lambda <- seq(0, 5, by = 0.2)
fit_NonParBJ <- bujar(x = X, y = Z, cens = Cens, learner = "lasso", valdata = cbind(Z, Cens, X), lambda = lambda,
                      nfold = 5, foldid = foldid, cv = TRUE, tuning = TRUE, trace = TRUE)

#Quantile regression
#fit_QuantReg <-  rqPen::cv.rq.pen(X, Y, tau=.5, lambda=NULL, weights=NULL, penalty="LASSO",
#                                  nfolds = 5,foldid=foldid)

#Gauss BJ
lambda <- seq(0, 5, by = 0.2)
fit_GaussBJ <- cvGaussBJ(Y = Y, X = X, s = s, Cens = Cens,
                         lambda = lambda, groups = foldid)
