cvGaussBJ = function(Y, X, s, Cens, zero = 10^(-10), maxIter = 1500,
                   lambda,  penalty.factor = rep(1, dim(X)[2]), groups){
  
  if(sum(search()=="package:glmnet")==0){library(glmnet)}
  if(sum(search()=="package:foreach")==0){library(foreach)}
  if(sum(search() == "package:parcor") == 0){library(parcor)}
  
  MSE <- c()
  
  for(ii in 1:length(lambda)){
    mse <- c()
    for(l in 1:max(groups)){	
      
      #Initialisation
      iter = 1              
      test = c(1, 1, 1)
      beta.hat <- matrix(ncol = dim(X)[2] + 2)
      sigma2.hat <- c()
      
      Y.star <- Y[groups != l]
      glmnet.est = glmnet(as.matrix(X[groups != l,]), Y.star, lambda = lambda[ii], alpha = 1, 
                          intercept = TRUE, standardize = TRUE,
                          penalty.factor =  penalty.factor)
      beta.hat[iter,] = c(glmnet.est$a0, as.vector(glmnet.est$beta))
      sigma2.hat[iter] = sum((Y.star - cbind(1,as.matrix(X[groups != l,])) %*% beta.hat[iter,])^2) / length(Y.star)
      
      #Iteration
      while(sum(test < zero) < 3 && iter < maxIter && cycleperiod < maxcycle){ 
        #     cat(iter, "\n")
        iter = iter + 1
        beta.old = beta.hat[iter - 1, ]
        sigma2.old = sigma2.hat[iter - 1]
        Y.star.old = Y.star
        # Imputation of Y:
        C <- which(Cens[groups != l] == 0)
        if (length(C) > 0) {
          Y.star[C] = cbind(1,X[groups != l,])[C,] %*% beta.old - (sigma2.old * dnorm(s[C] - cbind(1,X[groups != l,])[C,] %*% beta.old, sd = sigma2.old))/pnorm(s[C] - cbind(1,X[groups != l,])[C,]%*%beta.old, sd = sigma2.old)  
        } else break
        glmnet.est = glmnet(X[groups != l,], Y.star, lambda = lambda[ii], alpha = 1, intercept = TRUE, standardize = TRUE,
                            penalty.factor = penalty.factor)
        
        beta.hat <- rbind(beta.hat, c(glmnet.est$a0, as.vector(glmnet.est$beta)))
        sigma2.hat[iter] = sum((Y.star - cbind(1, X[groups != l,]) %*% beta.hat[iter,])^2) / length(Y.star)
        test = c(supbeta = max(abs(beta.hat[iter,] - beta.old)), sigma2 = max(abs(sigma2.hat[iter] - sigma2.old)),
                 Yhat = max(abs(Y.star - Y.star.old)))
        cycleperiod  = length(which((round(sigma2.hat[iter],4) == round(sigma2.hat[-iter],4)) == "TRUE"))
        
      }
      C <- which(Cens[groups == l] == 0)
      if (length(C) > 0) {
        Y[which(groups == l)][C] = cbind(1,X[which(groups == l),])[C,] %*% beta.hat[iter,] - (sigma2.hat[iter] * dnorm(s[which(groups == l)][C] - cbind(1,X[which(groups == l),])[C,] %*% beta.hat[iter,],  sigma2.hat[iter]))/pnorm(s[which(groups == l)][C] - cbind(1,X[which(groups == l),])[ C,]%*%beta.hat[iter,],  sigma2.hat[iter])
        mse[l] = sum((Y[groups == l] - cbind(1,X[which(groups == l),])%*%beta.hat[iter,])^2)   
      } else {
        mse[l] = sum((Y[groups == l] - cbind(1,X[which(groups == l),])%*%beta.hat[iter,])^2)
      }
      
    }
    MSE[ii] <- mean(mse)
  }
  
  # browser()
  inds = which.min(MSE)
  lambda.opt <- lambda[inds]
  
  Res <- list(CV = MSE, lambda = lambda, lambda.opt = lambda.opt)
  
  
  return(Res)
}
