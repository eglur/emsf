source("dp.R")

pisf <- function(Da, K, rb, gamma, pi = NULL)
{
   num.actions <- dim(Da)[3]
   num.states  <- ncol(K)
   m <- nrow(K)
 
   if (is.null(pi)) pi <- sample(1:num.actions, num.states, TRUE)
   pi2 <- pi
   pi2[1] <- pi2[1] + 1; if (pi2[1] > num.actions) pi2[1] <- 1 # to make sure pi and pi2 are different
   
   PIS <- pi

   Qt <- NULL
   while (sum(pi != pi2) > 0)
   {
      Dpi <- matrix(0, num.states, m)
      for (i in 1:num.states) Dpi[i,] <- Da[i,,pi[i]]
      Pb <- K %*% Dpi
      vb <- solve((diag(nrow(Pb)) - gamma * Pb), rb)
      
      Qt <- matrix(0, num.states, num.actions)
      
      for (i in 1:num.actions) Qt[,i] <- Da[,,i] %*% vb

      pi2 <- pi
      pi <- apply(Qt, 1, which.max)
      PIS <- rbind(PIS, pi)
   }
   list(pi = pi, Q = Qt, PIS = PIS)
}

test.sf.normalized <- function(n, m, num.actions, gamma = 0.95)
{
   Da <- array(runif(n*m * num.actions), c(n, m, num.actions))
   for (a in 1:num.actions) 
   {
      Da[,,a] <- Da[,,a] / apply(Da[,,a], 1, sum)
#       pd(Da[,,a])
#       print(apply(Da[,,a],1,sum))
#       print(""); print(""); print("")
   }
   
   K <- matrix(runif(n*m), m, n)
   K <- K / apply(K, 1, sum)
   
   rb <- runif(m)
#    rb[1] <- -1000
   
   Pa <- array(0, c(n, n, num.actions))
   ra <- matrix(0, n, num.actions)
   
   Pab <- array(0, c(m, m, num.actions))
   rab <- matrix(0, m, num.actions)

   for (a in 1:num.actions)
   {

      Pa[,,a] <- Da[,,a] %*% K
      ra[,a] <- Da[,,a] %*% rb

      Pab[,,a] <- K %*% Da[,,a] 
      rab[,a] <- rb 
      
#       print(ra[,a])
#       print(rab[,a])
   }
   

   Y <- policy.iteration(ra, Pa, gamma)
    ## PISF ##
   X <- pisf(Da, K, rb, gamma)
   pit2 <- X$pi
   
    ## POLICY ITERATION ##
   X <- policy.iteration(rab, Pab, gamma)
   Qt <- matrix(0, n, num.actions)
   
#    X$Q <- X$Q / apply(X$Q, 1, sum)  ## nao funciona 
   for (a in 1:num.actions) 
   {
      
       X$Q[,a] <- X$Q[,a] - mean(X$Q[,a]) ## funciona razoavelmente (mas nao perfeitamente)
#       X$Q[,a] <- X$Q[,a] / sum(X$Q[,a]) ## nao funciona 
      Qt[,a] <- Da[,,a] %*% X$Q[,a]
#       Qt[,a] <- Qt[,a] - mean(Qt[,a])
   }
   
#    V <- apply(Y$Q, 1, max)
#    Vt <- apply(Qt, 1, max)
#    
#    matplot(cbind(V, Vt), t="l")
   
   pit <- apply(Qt, 1, which.max)
   
#    print(Y$pi)
#    print(pit)
   
   print(c(paste("Erro PI normalizado:" , sum(pit != Y$pi)), (paste("Erro PISF:" , sum(pit2 != Y$pi)))))
}