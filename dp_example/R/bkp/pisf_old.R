source("data.manipulation.R")
source("dp.R")

solve.mp.exactly <- function(R, P, df, precision = 1e-6) {
   solve((diag(nrow(P)) - df*P), R) 
   }


pisf <- function(SF, df, pi = NULL, max.iter = Inf,
                 solve.mp.function = solve.mp.exactly) {
# df is the discount factor
# m is the number of archetypes
# pi is the initial policy
# max.iter is the maximum number of iterations 
# solve.mp.function is the function used to compute the value of a policy
#    (this can be used to implement the "modified policy iteration algorihtm",
#     which is the method recommended by Puterman (p.186))
# factor.function is the function used to approximate the SFDP stochastically
   
   num.actions <- dim(SF$D)[3]
   num.states <-  ncol(SF$K)
   m <- nrow(SF$K)

   if (is.null(pi)) pi <- sample(1:num.actions, num.states, TRUE)

   pi.old <- array(0, num.states)
   Q <- matrix(0, num.states, num.actions)
      

   Dpi <- matrix(0, num.states, m)
   Ppi <- matrix(0,m,m)
   it <- 0
   if (max.iter == Inf) max.iter <- length(pi)
   
   while (it < max.iter && !identical(pi, pi.old)) {
      
      # policy evaluation
      for (i in 1:nrow(Dpi)) Dpi[i,] <- SF$D[i,,pi[i]]
      Ppi <- SF$K %*% Dpi
      V <- solve.mp.function(SF$r, Ppi, df)
    
      for (a in 1:num.actions) {
         Q[,a] <- SF$D[,,a] %*%  V  
         }

      pi.old <- pi
      pi <- apply(Q,1,which.max)
      
      it <- it + 1
      }
   
   list(pi = pi, Q = Q)
   }


## Regular factorization
##------------------------

nmf.kmeans.smooth <- function(W, m, tau, max.iter, ...) {
# K-means factorization
   A <- kmeans(normalize(W), m, iter.max = max.iter, ...)
   D <- matrix(0, nrow(W), m)
   for (i in 1:nrow(D)) {
      for (j in 1:ncol(D)) {
         D[i,j] <- exp(-sum((W[i,] - A$centers[j,])^2) / tau)
         }
      }
   D <- D / apply(D,1,sum)
   
   list(D = D, K = A$centers[,2:ncol(A$centers)], r = A$centers[,1])
   }


factor.mdp <- function(R, P, m, tau, 
                        max.iter = 5, 
                        factor.function = nmf.kmeans.smooth, ...) {
# stochastic approximation of an MPD
# R is |S| x |A|
# P is |S| x |S| x |A|
     
     num.states <- nrow(R)
     num.actions <- ncol(R)
     
     P2 <- P[,,1]
     R2 <- as.matrix(R[,1])
     if (num.actions > 1) {
      for (i in 2:num.actions) {
            R2 <- rbind(R2, as.matrix(R[,i]))
            P2 <- rbind(P2, P[,,i])
            }
       }
       
     W <- cbind(R2,P2)
     
     Q <- factor.function(W, m, tau, max.iter, ...)
     D <- array(0, c(num.states, m, num.actions))
     for (i in 1:num.actions) {
          b <- (i-1)* num.states + 1
          e <- b + num.states - 1
          D[,,i] <- Q$D[b:e, ]
          }
     list(D=D, K = Q$K, r = Q$r)
     }



## Factorization in the space of features
##------------------------

kb.factor.mdp <- function(R, P, S, mda, tau, 
                          max.iter = 5, 
                          factor.function = kmeans, ...) {
   num.states  <- nrow(R)
   num.actions <- ncol(R)
   
   A <- factor.function(normalize(S), mda, iter.max = max.iter, ...)
   T <- matrix(0, num.states, mda)
   for (i in 1:nrow(T)) {
      for (j in 1:ncol(T)) {
         T[i,j] <- exp(-sum((S[i,] - A$centers[j,])^2) / tau)
         }
      }
   T <- T / apply(T,1,sum)
   
   T2 <- T
   ind <- array(0, mda)
   for (i in 1:mda) {
      ind[i] <- which.min(T2[,i])
      T2[ind[i],] <- Inf # so it won't be picked twice
      }

   D <- array(0, c(num.states, mda * num.actions, num.actions))
   K <- NULL
   r <- NULL
   for (a in 1:num.actions) {
      K <- rbind(K, P[ind,,a])
      r <- c(r, R[ind,a])
      b <- (a-1) * num.actions + 1
      e <- b + mda - 1
      D[,b:e,a] <- T
      }
  
   list(D = D, K = K, r = r)
   }


pisf.new <- function(Da, K, rb, gamma, pi = NULL)
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


test.bound.pisf <- function(n, m, num.actions, gamma = 0.99, num.runs = 100)
{
   for (run in 1:num.runs)
   {
     Da <- array(runif(n * m * num.actions), c(n, m, num.actions)) 
     for (a in 1:num.actions) for (i in 1:n) Da[i,sample(1,1:m),a] <- 10000 ## remove?
     for (i in 1:num.actions) Da[,,i] <- Da[,,i] / apply(Da[,,i], 1, sum)
     K <- matrix(runif(n * m), m, n)
     for (i in 1:m) K[i,sample(1,1:n)] <- 10000 ## remove?
     K <- K / apply(K,1,sum)
      
     rb <- runif(m)
#      rb[1] <- 100

     Qt <- pisf.new(Da, K, rb, gamma)$Q
     
     P <- array(0, c(n,n,num.actions))
     r <- matrix(0, n, num.actions)
     for (i in 1:num.actions) 
     {
        P[,,i] <- Da[,,i] %*% K
        r[,i] <- Da[,,i] %*% rb
     }
     
     Q <- policy.iteration(r, P, gamma)$Q
    
    diff <- sum(Q - Qt)^2 / length(Q)
    print(paste("Diff PISF and MT", diff))
    if (diff > 1e-10) 
    {
       print("Something wrong in the relation of PISF and MT")
       break
    }
    
    pi <- sample(1:num.actions, n, TRUE)
    PIS <- pisf.new(Da, K, rb, gamma, pi)$PIS
    PIS2 <- policy.iteration(r, P, gamma, pi)$PIS
    diff <- sum(PIS != PIS2)
    print(paste("POLICY HISTORIES",diff))
    
    if (diff > 0) 
    {
       print("Something wrong in the relation of PISF and MT; policy histories do not match")
       break
    }

     P3 <- array(0, c(m,m,num.actions))
     r3 <- matrix(0, m, num.actions)
     for (i in 1:num.actions) 
     {
        P3[,,i] <- K %*% Da[,,i]
        r3[,i] <- rb
     }
     
     Qtmp <- policy.iteration(r3, P3, gamma)$Q
     vtmp <- apply(Qtmp, 1, max)
     Q3 <- matrix(0, n, num.actions) 
    for (i in 1:num.actions) Q3[,i] <- Da[,,i] %*% vtmp
    Q3 <- Q3 - mean(Q3)
    Q4 <- Q - mean(Q)
    diff <- sum(Q4 - Q3)^2 / length(Q)
    print(paste("Diff PISF and separated factorizations", diff))
    if (diff > 1e-10) 
    {
       print("Something wrong in the relation of PISF and separated factorizations")
       break
    }



    P2 <- array(runif(n*n*num.actions), c(n,n,num.actions))
    for (i in 1:num.actions) P2[,,i] <- P2[,,i] / apply(P2[,,i], 1, sum)
    r2 <- matrix(runif(n*num.actions), n, num.actions)
    r2[1,1] <- 200
    
    Q2 <- policy.iteration(r2, P2, gamma)$Q
    
    Cp <- -Inf
    Cr <- -Inf
    for (i in 1:num.actions) 
    {
       C <- max(apply(abs(P[,,i] - P2[,,i]), 1, sum))
       if (C > Cp) Cp <- C
       
       C <- max(abs(r[,i] - r2[,i]))
       if (C > Cr) Cr <- C
          
    }
    
    deltar <- max(max(r), max(r2)) - min(min(r), min(r2))
    
    bound <- 1 / (1 - gamma * (1 - 0.5 * Cp)) * (Cr + gamma / (2 * (1 - gamma)) * Cp *deltar)
    
    v <- apply(Q,1,max)
    v2 <- apply(Q2, 1, max)
    diff <- max(abs(v - v2))
   
   print(paste("Diff MT and M", round(diff, 2), "   Bound", round(bound,2)))
      
   if (diff > bound) 
   {
      print("Something wrong with the bound")
      break
   }
   }
}

print("pisf.R loaded")