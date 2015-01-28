

generate.stochastic.matrix <- function(nrows, ncols)
{
   A <- matrix(runif(nrows * ncols), nrows, ncols)
   A <- A / apply(A, 1, sum)
   A
}

generate.stochastic.matrices <- function(nrows, ncols, na)
{
   A <- array(0, c(nrows, ncols, na))
   for (a in 1:na) A[,,a] <- generate.stochastic.matrix(nrows, ncols)
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


generate.model <- function(n, m, na)
{
    D <- generate.stochastic.matrices(n, m, na)
    K <- generate.stochastic.matrices(m, n, na)
    
    P <- array(0, c(n,n,na))
    for (a in 1:na) P[,,a] <- D[,,a] %*% K[,,a]
        
    mu <- runif(n)
    mu <- mu / sum(mu)
    pi <- generate.stochastic.matrix(n, na)
    
    list(P = P, pi = pi, mu = mu)
}

generate.data <- function(P, mu, pi, T)
{
    y <- array(0, T)
    a <- array(0, T-1)
    
    y[1] <- sample.from.dist(mu)
    for (i in 1:(T-1))
    {
        a[i] <- sample.from.dist(pi[y[i],])
        y[i + 1] <- sample.from.dist(P[y[i],, a[i]])
    }
    
    list(y = y, a = a)
}

generate.batch.data <- function(P, mu, pi, T, num.batches)
{
   y <- NULL
   a <- NULL
   
   for (i in 1:num.batches) 
   {
      dt <- generate.data(P, mu, pi, T)
      y <- c(y, list(dt$y))
      a <- c(a, list(dt$a))
   }
   
   list(y = y, a = a)
}


generate.model.and.data <- function(n, m, na, T, num.batches)
{
   md <- generate.model(n,m,na)
   dt <- generate.batch.data(md$P, md$mu, md$pi, T, num.batches)
   list(P = md$P, pi = md$pi, mu = md$mu, y = dt$y, a = dt$a)
}
    
em.sf <- function(P, OBS, ACT, m, mu, pi, eps = 1e-20, max.it = 10)
{
# Assumes that each observation and each action has appeared at least once
    n <- length(unique(unlist(OBS)))
    na <- length(unique(unlist(ACT)))
   
    D <- generate.stochastic.matrices(n,m,na)
    K <- generate.stochastic.matrices(m,n,na)
    

    old.score <- -Inf 
    score <- 0
    it <- 0
    
    error <- NULL
       
    while (abs(score - old.score) > eps || it < max.it)
    {
        
       old.score <- score 
       score <- 0
        
        
        D2 <- array(0, c(n, m, na))
        K2 <- array(0, c(m, n, na))

       
        
       for (batch in 1:length(OBS))
       {
        
        y <- OBS[[batch]]
        a <- ACT[[batch]]

        T <- length(y)
        
        A <- matrix(0, T-1, m)
        B <- matrix(1, T-1, m)
        NF <- array(0, T-1)
        
        
        A[1,] <- mu[y[1]] * pi[y[1], a[1]] * D[y[1],,a[1]]
        NF[1] <- sum(A[1,])
        A[1,] <- A[1,] / NF[1]
        
       
        for (t in 2:(T-1))
        {
            aa <- sum(A[t-1,] * K[, y[t], a[t-1]])
            A[t,] <- aa * pi[y[t],a[t]] * D[y[t],,a[t]]
            NF[t] <- sum(A[t,])
            A[t,] <- A[t,] / NF[t] 
        }
        
                
#         B[T-1, ] <- 1 / NF[T-1]
        for (t in (T-2):1)
        {
            bb <- sum(B[t+1,] * D[y[t+1], , a[t+1]])
            B[t,] <- bb * pi[y[t+1], a[t+1]] * K[, y[t+1], a[t]]
            B[t,] <- B[t,] / NF[t+1] 
        }
        
         
        C <- A * B 
#         C <- C / apply(C, 1, sum)
        
        for (t in 1:(T-1))
        {
            D2[y[t], , a[t]] <- D2[y[t], , a[t]] + C[t,]
            K2[, y[t+1], a[t]] <- K2[, y[t+1], a[t]] + C[t,]
        }
        
        score <- score - sum(log(1/NF)) ## is this correct?
      }

        
      for (i in 1:na) 
      {
         D2[,,i] <- D2[,,i] / apply(D2[,,i], 1, sum)
         K2[,,i] <- K2[,,i] / apply(K2[,,i], 1, sum)
      }
      
      D <- D2
      K <- K2
      
      e <- 0
      for (i in 1:na) e <- e + sum((D[,,i] %*% K[,,i] - P[,,i])^2)
      
      error <- c(error, e)
      plot(error)

      print(score)
       
   }
        
    list(D = D, K = K)
}
    

print("em.sf.R loaded") 
    
