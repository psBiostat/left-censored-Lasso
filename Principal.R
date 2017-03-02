library(NADA)
#Data
data(TCEReg)
head(TCEReg)

#Description
with(TCEReg, cenboxplot(TCEConc, TCECen))

#Add 100 noise variables
n <- dim(TCEReg)[1]
p <- 100
Xnoise <- matrix(rnorm(n*p, 0, 1), ncol = p)

X <- cbind(TCEReg$LandUse, TCEReg$PopDensity, TCEReg$PctIndLU, TCEReg$Depth, Xnoise)
Y <- TCEReg$TCEConc
Cens <- TCEReg$TCECen
Cens[TCEReg$TCECen == TRUE] <- 0
Cens[TCEReg$TCECen == FALSE] <- 1
s <- rep(0, n) 
s[Cens == 0] <- 5 #Threshold variable

#Fold CV
foldid <- rep(0, dim(TCEReg)[1])
foldid[which(TCEReg$TCECen == TRUE)] <- sample(length(which(TCEReg$TCECen == TRUE))) %% 5 + 1
foldid[which(TCEReg$TCECen == FALSE)] <- sample(length(which(TCEReg$TCECen == FALSE))) %% 5 + 1

#Lasso LOD
fit_LOD = glmnet::cv.glmnet(X, Y, nfolds = 5, foldid = foldid, alpha=1)

#Lasso Logistic
fit_Logistic = glmnet::cv.glmnet(X, Cens, family = "binomial", alpha = 1, nfolds = 5, foldid = foldid)

#Lasso ReverseCox
M = max(Y) + 1 
Z = M - Y
Surv <- cbind(Z, Cens)
colnames(Surv) <- c("time", "status")
fit_ReverseCox = glmnet::cv.glmnet(X, Surv, family = "cox", alpha = 1, nfolds = 5, foldid = foldid)

#Lasso NonParaBJ <- version modifié du package bujar
lambda <- seq(0, 5, by = 0.2)
fit_NonParBJ <- bujar(x = X, y = Z, cens = Cens, learner = "lasso", lambda = lambda,
                      nfold = 5, foldid = foldid, cv = TRUE, tuning = TRUE)

#Quantile regression
fit_QuantReg <-  rqPen::cv.rq.pen(X, Y, tau=.5, lambda=NULL, weights=NULL, penalty="LASSO",
                                     nfolds = 5,foldid=foldid)

#Gauss BJ
fit_GaussBJ <- cvGaussBJ(Y = Y, X = X, s = s, Cens = Cens, lambda = lambda, groups = foldid)