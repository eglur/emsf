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
   rb[1] <- -1000
   
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
#    X <- pisf(Da, K, rb, gamma)
#    pit <- X$pi
   
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
   
   V <- apply(Y$Q, 1, max)
   Vt <- apply(Qt, 1, max)
   
   matplot(cbind(V, Vt), t="l")
   
   pit <- apply(Qt, 1, which.max)
   
#    print(Y$pi)
#    print(pit)
   
   print(sum(pit != Y$pi))
}

generate.stochastic.matrix <- function(nrows, ncols)
{
   A <- matrix(runif(nrows * ncols), nrows, ncols)
   A <- A / apply(A, 1, sum)
   A
}

sample.from.dist <- function(D)
{
   v <- runif(1)
   vv <- 0
   ind <- 0
   while (vv < v) 
   {
      ind <- ind + 1
      vv <- vv + D[ind]
   }
   ind
}


soft.counting.outer <- function(n, T, num.batches)
{
   P <- generate.stochastic.matrix(n,n)
   Q <- generate.stochastic.matrix(n,n)
   
   Z <- Q %*% P %*% Q
   
   counts <- matrix(0, n, n)
   soft.counts <- matrix(0, n, n)
   error <- NULL
   serror <- NULL
   
   for (b in 1:num.batches)
   {

      s <- sample(1:n, 1)
      for (i in 1:T)
      {
         
         sp <- sample.from.dist(Q[s,])
         sp <- sample.from.dist(P[sp,])
         s2 <- sample.from.dist(Q[sp,])         
         
         # conventional
         counts[s,s2] <- counts[sp,s2] + 1

         #proposed 
         soft.counts <- soft.counts + Q[s,] %o% Q[s2,]
         
         s <- s2
      }
      
      Zt <- counts / apply(counts, 1, sum)
      sZt <- soft.counts / apply(soft.counts, 1, sum)
      
  
      error <- c(error, sum((Zt - Z)^2))
      serror <- c(serror, sum((sZt - Z)^2))
      
     
      matplot(cbind(error, serror), t="l")
      
   }

   error
}


em.sf.imp <- function(M, max.it = 10)
{
   n <- M$n
   m <- M$m
   na <- M$na
   P <- M$P
   OBS <- M$y
   ACT <- M$a 
   
   
# Assumes that each observation and each action has appeared at least once
#     n <- length(unique(unlist(OBS)))
#     na <- length(unique(unlist(ACT)))
   
    D <- generate.stochastic.matrices(n,m,na)
    K <- generate.stochastic.matrices(m,n,na)
    
    G  <- array(0, c(n,n,na))
    rr <- matrix(rep(1:n, na), n, na) 
    
    D2 <- D
    K2 <- K
    
    C <- array(0, c(n,n,na))
    for (batch in 1:length(OBS))
    {
    
        y <- OBS[[batch]]
        a <- ACT[[batch]]

        T <- length(y)
        
        for (i in 1:(T-1)) C[y[i],y[i+1],a[i]] <- C[y[i],y[i+1],a[i]] + 1 ## counting
        
    }
    
#     cnt <- order(apply(C,1,sum))
    
   
#     baseline.ll <- frobenius(P, Pt)
    
    old.score <- 0 # eps + 1 
    score <- 0
    it <- 0
    
    error <- NULL
    error2 <- NULL
       
    while (it < max.it) # abs(score - old.score) > eps && 
    {
        
        old.score <- score 
        score <- 1
        
        for (a in 1:na)
        {
            # original version # 
            Dh <- matrix(0, n, m)
            Kh <- matrix(0, m, n)
            
            for (i in 1:n)
            {
                for (j in 1:n)
                {
                    if (C[i,j,a] != 0) 
                    {
                        tmp <- D[i,,a] * K[,j,a]
                        g <- sum(tmp)

                        w <- C[i,j,a] * tmp / g
                        Dh[i,] <- Dh[i,] + w
                        Kh[,j] <- Kh[,j] + w

                    }
                }
            }
            
            
            
            Dh <- Dh / apply(Dh, 1, sum)
            Kh <- Kh / apply(Kh, 1, sum)
            
            D[,,a] <- Dh
            K[,,a] <- Kh

            # improved version # 
            Dh <- matrix(0, n, m)
            Kh <- matrix(0, m, n)
            
            for (i in 1:n)
            {
                for (j in 1:n)
                {
                    if (C[i,j,a] != 0) 
                    {
                        tmp <- D2[rr[i,a],,a] * K2[,j,a]
                        G[i,j,a] <- sum(tmp)

                        w <- C[i,j,a] * tmp / G[i,j,a]
                        Dh[i,] <- Dh[i,] + w
                        Kh[,j] <- Kh[,j] + w
                     
                    }
                    
                }

            }
            
            Kh <- Kh / apply(Kh, 1, sum)
            Dh <- Dh / apply(Dh, 1, sum)
    
            D2[,,a] <- Dh 
            K2[,,a] <- Kh
            
            
            ## swap order of the rows of Da to improve likelhood 
            rr[,a] <- 1:n # the assignment above "resets" the order
            h <- array(0, n)
            for (i in 1:n) h[i]  <- sum(C[i,,a] * G[rr[i,a],,a])
#             score <- sum(H)
            for (i in 1:(n-1))
            {
               for (j in (i+1):n)
               {
                  hij <- sum(C[i,,a] * G[rr[j,a],,a])
                  hji <- sum(C[j,,a] * G[rr[i,a],,a])
                  
                  if (hij  + hji > h[rr[i,a]] + h[rr[j,a]]) 
                  { # swap
                     tmp <- rr[i,a]
                     rr[i,a] <- rr[j,a]
                     rr[j,a] <- tmp
                  }
               }
            
#             print(rr[,a])   
            }     
        }
      

        Pt <- array(0, c(n,n,na))
        for (a in 1:na) Pt[,,a] <- D[,,a] %*% K[,,a]
        e <- log.likelihood(Pt, OBS, ACT)
#         e <- frobenius(P, Pt)
        error <- c(error, e)
        
        Pt2 <- array(0, c(n, n, na))
        for (a in 1:na) 
        {
           for (i in 1:n) Pt2[i,,a] <- D2[rr[i,a],,a] %*% K2[,,a] #  right order at this point
        }
        
        e <- log.likelihood(Pt2, OBS, ACT)
#         e <- frobenius(P, Pt2) 
        error2 <- c(error2, e)
        
        matplot(cbind(error, error2), t="l")
        
        print(error2[length(error2)])
        
        it <- it + 1

   }
        
    list(D = D, K = K)
}

