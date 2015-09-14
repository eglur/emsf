on.line.emsf <- function(P, m, max.it = 1000)
{
   
   n <- nrow(P)
   D <- generate.stochastic.matrix(n,m)
  
   K <- generate.stochastic.matrix(m,n)
   cc <- array(0, n)
   cr <- array(0, n)
   sr <- apply(K, 1, sum)
   
   s <- sample(1:n, 1)
   
   error <- NULL
   for (i in 1:max.it)
   {
      s2 <- sample.from.dist(P[s,])
      w <- D[s,] * K[,s2]
      pij <- sum(w)
      w <- w / pij
      
      # update s2-ith column of K
      oc <- K[,s2]
      K[,s2] <-  (K[,s2] * cc[s2] + w) 
      cc[s2] <- cc[s2] + 1
      K[,s2] <- K[,s2] / cc[s2]
      sr <- sr - oc + K[,s2]
      
     
      # update s-th row of D
      D[s,] <- (D[s,] * cr[s] + w)
      cr[s] <- cr[s] + 1
      D[s,] <- D[s,] / cr[s]
      D[s,] <- D[s,] / sum(D[s,])
      
      K2 <- K / sr # apply(K, 1, sum)
      e <- sum((P - D %*% K2)^2)
      error <- c(error,e)
      plot(error, t="l")
      
      s <- s2
      
   }
      
}
