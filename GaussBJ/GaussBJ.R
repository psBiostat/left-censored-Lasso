GaussBJ = function(Y, X, s, Cens, zero = 10^(-10), maxIter = 1500,
                      lambda,  penalty.factor = rep(1, dim(X)[2]), groups){
  
  if(sum(search()=="package:glmnet")==0){library(glmnet)}
  if(sum(search()=="package:foreach")==0){library(foreach)}
  if(sum(search() == "package:parcor") == 0){library(parcor)}
      
      #Initialisation
      iter = 1            	
      test = c(1, 1, 1)	
      beta.hat <- matrix(ncol = dim(X)[2] + 2)
      sigma2.hat <- c()
      Y.star <- Y 
      glmnet.est = glmnet(as.matrix(X), Y.star, lambda = lambda, alpha = 1, 
                          intercept = TRUE, standardize = TRUE,
                          penalty.factor =  penalty.factor)
      beta.hat[iter,] = c(glmnet.est$a0, as.vector(glmnet.est$beta))
      sigma2.hat[iter] = sum((Y.star - cbind(1,as.matrix(X)) %*% beta.hat[iter,])^2) / length(Y.star)
      
      #Iteration
      while(sum(test < zero) < 3 && iter < maxIter && cycleperiod < maxcycle){ 
        #     cat(iter, "\n")
        iter = iter + 1
        beta.old = beta.hat[iter - 1, ]
        sigma2.old = sigma2.hat[iter - 1]
        Y.star.old = Y.star
        # Imputation of Y:
        C <- which(Cens == 0)
        if (length(C > 0)) {
          Y.star[C] = cbind(1,X)[C,] %*% beta.old - (sigma2.old * dnorm(s[C] - cbind(1,X)[C,] %*% beta.old, sd = sigma2.old))/pnorm(s[C] - cbind(1,X)[C,]%*%beta.old, sd = sigma2.old)  
        } else break
        glmnet.est = glmnet(X, Y.star, lambda = lambda, alpha = 1, intercept = TRUE, standardize = TRUE,
                            penalty.factor = penalty.factor)
        
        beta.hat <- rbind(beta.hat, c(glmnet.est$a0, as.vector(glmnet.est$beta)))
        sigma2.hat[iter] = sum((Y.star - cbind(1, X ) %*% beta.hat[iter,])^2) / length(Y.star)
        test = c(supbeta = max(abs(beta.hat[iter,] - beta.old)), sigma2 = max(abs(sigma2.hat[iter] - sigma2.old)),
                 Yhat = max(abs(Y.star - Y.star.old)))
        cycleperiod  = length(which((round(sigma2.hat[iter],4) == round(sigma2.hat[-iter],4)) == "TRUE"))
        
      }
  
  Res <- list(BJ = list(beta = beta.hat, sigma2 = sigma2.hat, iter = iter),
              Beta.Est = as.vector(beta.hat[iter,]),
              Sigma2.Est = sigma2.hat[iter],
              Y.imput = Y.star)
  
  
  return(Res)
}
