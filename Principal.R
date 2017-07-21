library(NADA)
library(glmnet)
library(hdnom)
library(bujar)
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
#We choose the detection limit at 0.05 ug/L to have 10% of left-censoring data
Y[which(Recon$AtraConc <= 0.05)] <- 0.05
Cens <- Y
Cens[Y==0.05] <- 0
Cens[Y>0.05] <- 1
s <- rep(0, n) 
s[Cens == 0] <- 0.05 #Threshold variable

#Fold CV
foldid <- rep(0, n)
foldid[which(Cens == 0)] <- sample(length(which(Cens == 0))) %% 5 + 1
foldid[which(Cens == 1)] <- sample(length(which(Cens == 1))) %% 5 + 1

#Lasso LOD
fit_LOD = glmnet::cv.glmnet(X, Y, nfolds = 5, foldid = foldid, alpha=1)

#Lasso NonParaBJ <- version modifié du package bujar
lambda <- seq(0, 5, by = 0.2)
fit_NonParBJ <- bujar::bujar(x = X, y = Z, cens = Cens, learner = "lasso", lambda = lambda,
                      nfold = 5, foldid = foldid, cv = TRUE, tuning = TRUE, trace = TRUE)

#Quantile regression
fit_QuantReg <-  rqPen::cv.rq.pen(X, Y, tau=.5, lambda=NULL, weights=NULL, penalty="LASSO",
                                  nfolds = 5,foldid=foldid)

#Gauss BJ
lambda <- seq(0, 5, by = 0.2)
fit_GaussBJ <- cvGaussBJ(Y = Y, X = X, s = s, Cens = Cens,
                         lambda = lambda, groups = foldid)
