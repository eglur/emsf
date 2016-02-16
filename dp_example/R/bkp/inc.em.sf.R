
source("em.sf.R")

update.model <- function(D, K, Dh, Kh, x, y, llh, C)
{
   for (i in 1:nrow(C))
   {
      for (j in 1:ncol(C))
      {
         if (C[i,j] != 0)
         {
            g <- D[i,] %*% K[,j]
            w <- (C[i,j] / g) * (D[i,] * K[,j]) 

            Dh[i,] <- Dh[i,] +  w
            x[i] <- x[i] + sum(w) 
           
            Kh[,j] <- Kh[,j] + w 
            y <- y + w 

            llh <- llh + C[i,j] * log(g)
         }
      }
   }
   
   
   list(Dh = Dh, Kh = Kh, x = x, y = y, llh = llh)
}



inc.emsf <- function(n, m, na, P, OBS, ACT, eps = 1e-20, max.it = 10)
{

    D <- generate.stochastic.matrices(n,m,na)
    K <- generate.stochastic.matrices(m,n,na)
    
   
    C <- array(0, c(n,n,na))
    
    llh <- 0
    llhp <- Inf
    
    error <- NULL 
    
    while (abs(llh - llhp) > eps)
    {
        llhp <- llh
        llh <- 0
        
        Dh <- array(0, c(n, m, na))
        Kh <- array(0, c(m, n, na))
        x <- matrix(0, n, na)
        y <- matrix(0, m, na)        
        
        # counts transitions per batch (this is different in the paper)
        batch <- 1        
        while (batch <= length(OBS))
        {
           yy <- OBS[[batch]]
           aa <- ACT[[batch]]

           T <- length(yy)
        
           for (i in 1:(T-1)) C[yy[i],yy[i+1],aa[i]] <- C[yy[i],yy[i+1],aa[i]] + 1 ## counting
           
           for (a in 1:na)
           {
              R <- update.model(D[,,a], K[,,a], Dh[,,a], Kh[,,a], x[,a], y[,a], llh, C[,,a])
              
              Dh[,,a] <- R$Dh
              Kh[,,a] <- R$Kh
              x[,a]   <- R$x
              y[,a]   <- R$y
              llh <- R$llh
              
              C[,,a] <- 0
           }
           
           batch <- batch + 1
        }
        
       for (a in 1:na)
       {
          D[,,a] <- Dh[,,a] / x[,a]
          K[,,a] <- Kh[,,a] / y[,a]
       }
       
   e <- 0
   for (i in 1:na) e <- e + sum((D[,,i] %*% K[,,i] - P[,,i])^2)
   
   error <- c(error, e)
   plot(error)
       
   }
     
}



on.line.emsf <- function(P, pi, mu, na, m, tm, alpha = 0.01, max.it = 1000000)
{
   
   n <- nrow(P[,,1])
   D <- generate.stochastic.matrices(n,m, na)
   K <- generate.stochastic.matrices(m,n, na)
   C <- array(0, c(n, n, na))
   
   s <- sample(1:n, 1)
   
   llh <- 0
   error <- NULL
   for (t in 1:max.it)
   {
      a <- sample.from.dist(pi[s,])
      s2 <- sample.from.dist(P[s,,a])
      C[s,s2,a] <- C[s,s2,a] + 1
         
      if (!(t %%  tm)) # is it time to update model?
      {
         e <- 0
      
         for (a in 1:na)
         {
            Dh <- matrix(0, n, m)
            Kh <- matrix(0, m, n)
            x <- array(0, n)
            y <- array(0, m)
            R <- update.model(D[,,a], K[,,a], Dh, Kh, x, y, llh, C[,,a])
            Dh <- R$Dh
            Kh <- R$Kh
            x  <- R$x
            y  <- R$y
            llh <- R$llh
            
            for (i in 1:length(x)) 
            {
               if (x[i] != 0) # only updates the rows of D that have changed
               {
                  D[i,,a] <- (1 - alpha) * D[i,,a] + alpha * Dh[i,] / x[i]
               }
            }

            dim(y) <- NULL
            K[,,a] <- (1 - alpha) * K[,,a] + alpha * Kh / y 
            
            C[,,a] <- 0
            
#             print(D[,,a])
            
            e <- e + sum((P[,,a] - D[,,a] %*% K[,,a])^2)
         }
         error <- c(error,e)
         plot(error, t="l")
      }
      
      
      s <- s2
   }
      
}


run.experiment <- function(n = 10, m = 3, na = 2, T = 1000, num.batches = 100, max.it = 30)
{
    M <- generate.model.and.data(n, m, na, T, num.batches)

#     print("Running batch version of EMSF...")
#     em.sf2(M$n, M$m, M$na, M$P, M$y, M$a, max.it = max.it)
#     
#     x11()
    print("Running incremental version of EMSF...")
    inc.emsf(M$n, M$m, M$na, M$P, M$y, M$a, max.it = max.it)
    
    print("Done.")
}


run.experiment.online <- function(n = 10, na = 2, m = 3, tm = 1000, max.it = 10000000)
{
  
   M <- generate.model(n, m, na)

    print("Running online version of EMSF...")
    on.line.emsf(M$P, M$pi, M$mu, na, m, tm, max.it = max.it)
    
    print("Done.")
}

    
    
print("inc.em.sf.R loaded") 
    