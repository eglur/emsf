source("dp.R")
source("pisf.R")
source("markov.chains.R")

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


kld <- function(d1, d2)
{
  err <- 0
  for (i in 1:length(d1)) if (d1[i] > 0 && d2[i] > 0) err <- err + d1[i] * log(d1[i] / d2[i])
  err
}

KL <- function(P, Pt, dist = NULL)
{
  if (is.null(dist)) dist <- 1/nrow(P)
  err <- 0
  for (i in 1:nrow(P))
  {
    err <- err + dist[i] * kld(P[i,], Pt[i,])
  }
  err
}

emsf <- function(P, m, ta, tc, 
                 alpha = 1,
                 max.it = 10^10
                )
{
   
   n <- dim(P)[1]
   na <- dim(P)[3]
   
   dist <- matrix(0, n, na)
   for (a in 1:na) dist[,a] <- stationary.distribution(P[,,a]) # this is wrong if na > 1
  
   D <- array(runif(n * m * na), c(n, m, na))
   for (u in 1:na) D[,,u] <- D[,,u] / apply(D[,,u], 1, sum)
      
   K  <- matrix(runif(n * m), m, n)
   K <- K / apply(K, 1, sum)
   
   Dh <- array(0, c(n, m, na))
   Kh  <- matrix(0, m, n)
   
   x <- matrix(0, n, na)
   y  <- array(0, m)
   
   s <- sample(1:n, 1)
#    pi <- sample(1:na, n, TRUE)
   C <- array(0, c(n,n,na)) # sparse in practice
   
   ll <- 0
   it <- 1
   error <- NULL
   while (it <= max.it)
   {
#       a <- pi[s]
#       s2 <- trasition.function(s, a)
      
      a <- sample(1:na, 1)
      s2 <- sample.from.dist(P[s,,a])
         
      C[s,s2,a] <- C[s,s2,a] + 1
      
      if (it %% ta == 0)
      {
         for (u in 1:na)
         {
            for (i in 1:n)
            {
               for (j in 1:n)
               {
                  if (C[i,j,u] != 0)
                  {
                     g <- D[i,,u] %*% K[,j]
                     
                      if (g > 0) ## THINK
                     {
                        w <- (C[i,j,u] / g) * D[i,,u] * K[,j]
                        
                        Dh[i,,u] <- Dh[i,,u] + w
                        x[i,u] <- x[i,u] + sum(w)
                        
                        Kh[,j] <- Kh[,j] + w
                        y <- y + w
                     
                        ll <- ll + C[i,j,u] * log(g)
                     }
                     
                  }
               }
            }
             C[,,u] <- 0
         }
      }
      
      if (it %% tc == 0)
      {
        
         
         for (i in 1:m) if (y[i] != 0) Kh[i,] <- Kh[i,] / y[i]
         
         K <- (1 - alpha) * K + alpha * Kh

         Kh  <- matrix(0, m, n)
         y  <- array(0, m)      


         for (u in 1:na)
         {
            for (i in 1:n) 
            {
                if (x[i,u] != 0) D[i,,u] <- (1 - alpha) * D[i,,u] + alpha * Dh[i,,u] / x[i,u]
            }
            
         }
         
         Dh <- array(0, c(n, m, na))
         x <- matrix(0, n, na)

         err <- 0
         for (a in 1:na) err <- err + KL(P[,,a], D[,,a] %*% K, dist[,a])
         error <- c(error, err)
#          plot(error)
      }
   
   
   s <- s2
   it <- it + 1
   }
   error   
}






online.emsf <- function(P, m, tc, 
                 max.it = 10^10,
                 alpha = 1e-4, 
                 beta = 1e-4 
                )
{
   
   n <- dim(P)[1]
   na <- dim(P)[3]
   
   D <- array(runif(n * m * na), c(n, m, na))
   for (u in 1:na) D[,,u] <- D[,,u] / apply(D[,,u], 1, sum)
      
   K  <- matrix(runif(n * m), m, n)
   K <- K / apply(K, 1, sum)
   
   cd <- matrix(0, n, na)
   ck <- array(0, n)

   s <- sample(1:n, 1)
   
#    pi <- sample(1:na, n, TRUE)
   ll <- 0
   it <- 1
   while (it <= max.it)
   {
#       a <- pi[s]
#       s2 <- trasition.function(s, a)
      
      a <- sample(1:na, 1)
      s2 <- sample.from.dist(P[s,,a])
         
      w <- (D[s,,a] * K[,s2]) / (D[s,,a] %*% K[,s2])
      
      cd[s,a] <- cd[s,a] + 1
      ck[s2]  <- ck[s2] + 1
      
#       alpha <-  1 / cd[s,a]
#       beta <- 1 / ck[s2]
      
      D[s,,a] <- (1 - alpha) * D[s,,a] + alpha * w 
      
      
      K[,s2] <- (1 - beta) * K[,s2]  + beta * w 
     
     
      if (it %% tc == 0)
      {
         
         K <- K / apply(K, 1, sum)
         
         err <- 0
         for (i in 1:na) err <- err + sum((P[,,i] - D[,,i] %*% K)^2)
         print(err)
      
      }
   
   
   s <- s2
   it <- it + 1
   }
   
}



emsf.pisf <- function(P, rp, df, m, tc, 
                 epsilon = 0.1,
                 max.it = 10^10,
                 alpha = 1e-4, 
                 beta = 1e-4
                 
                )
{
   
   n <- dim(P)[1]
   na <- dim(P)[3]
   
   r <- matrix(0, n, na)
   for (u in 1:na) r[ , u] <- P[,,u] %*% rp
      
   F <- policy.iteration(r, P, df)
   Q <- F$Q
   op <- F$pi
   
   D <- array(runif(n * m * na), c(n, m, na))
   for (u in 1:na) D[,,u] <- D[,,u] / apply(D[,,u], 1, sum)
      
   K  <- matrix(runif(n * m), m, n)
   K <- K / apply(K, 1, sum)
   
#    rt <- matrix(0, n, na)
#    ct <- matrix(0, n, na)
      
   s <- sample(1:n, 1)
   
   pi <- sample(1:na, n, TRUE)
   
   se <- NULL
   it <- 1
   while (it <= max.it)
   {
      if (runif(1) < epsilon) a <- sample(1:na, 1)
      else a <- pi[s]
          
      s2 <- sample.from.dist(P[s,,a])
      
      g <- D[s,,a] * K[,s2]
      w <- g / sum(g)
      
      D[s,,a] <- (1 - alpha) * D[s,,a] + alpha * w 

      K[,s2] <- (1 - beta) * K[,s2]  + beta * w 
     
#       rt[s,a] <- rt[s,a] + r[s2,a]
#       ct[s,a] <- ct[s,a] + 1
      
      if (it %% tc == 0)
      {
         
         K <- K / apply(K, 1, sum)
         rb <- K %*% rp
         
         F <- pisf(D, K, rb, df)
         Qt <- F$Q
         pi <- apply(Qt, 1, which.max)

#          err <-  sum((Q - Qt)^2) # sum(op != pi) #sum((Q - Qt)^2)
         
         err <- 0
         for (u in 1:na) err <- err + sum((P[,,u] - D[,,u] %*% K)^2)
         
         print(sum(op != pi))
         
         se <- c(se, err)
         plot(se)
      }
   
   ## BACK!!!
   s <-  sample(1:n, 1) # s2  # sample(1:n, 1)
   it <- it + 1
   }
   
}