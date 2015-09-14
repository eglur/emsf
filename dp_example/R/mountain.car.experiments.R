source("mountain.car.R")
source("dp.R")
source("emsf.R")
source("util.R")
source("data.manipulation.R")


mc.ml <- function(cells, num.episodes = 5e3, epsilon = 0.15, tc = 10, df = 0.999, max.iter = 300,
  test.set = rbind(gp(seq(-1,0.15,l=5), seq(-0.07,0.02,l=5)), c(0,0)),
  max.steps = 300)  
{
  
  n <- cells[1] * cells[2] + 1
  na <- 3
  
  P <- array(0, c(n, n, na))
  r <- matrix(0, n, na)
  P[n,n,] <- 1
  
  C <- array(0, c(n, n, na)) # counting
  rs <- matrix(0, n, na) # sums

  pi <- sample(1:na, n, TRUE)
     
  R <- array(0, num.episodes %/% tc)
  T <- array(0, num.episodes %/% tc)
  
  start.time <- proc.time()[1]
  time.eval <- 0
  err <- NULL
  for (ep in 1:num.episodes)
  {
    
    s <- test.set[sample(1:nrow(test.set), 1), ]
    j <- mc.state.index(s, cells)
    
    t <- NULL
    st <- 1
    g <- FALSE
    while (!g && st <= max.steps)
    {
     
      i <- j
      if (runif(1) < epsilon) a <- sample(1:na, 1)
      else a <- pi[i]
      
      t <- mc.transition(s, a)
     
      s <- t$s
      j <- mc.state.index(s, cells)
      
      C[i,j,a] <- C[i,j,a] + 1
      rs[i,a] <- rs[i,a] + t$r
      
      g <- t$g
      st <- st + 1
    }
    
   if (ep %% tc == 0)
   {
      for (u in 1:na)
      {
         cnt <- apply(C[,,u], 1, sum)
         for (k in 1:length(cnt))
         {
            if (cnt[k] > 0)
            {
               P[k,,u] <- C[k,,u] / cnt[k]
               r[k,u] <- rs[k,u] / cnt[k]
            }
            else P[k,,u] <- 1 / ncol(P[,,u])
         }
      }
      

      Y <- policy.iteration(r, P, df = df, max.iter = max.iter)
      pi <- Y$pi
            
      V <- apply(Y$Q, 1, which.max)
      V <- V[-length(V)]
      V <- matrix(V, sqrt(length(V)), sqrt(length(V)))
      image(V)
      
      
      
      tt <- proc.time()[1]
      R[ep / tc] <- mc.evaluate.policy(pi, cells, test.set, max.steps)
      print(R[ep / tc])
      err <- c(err, R[ep / tc])
#       plot(err)
      time.eval <- time.eval + (proc.time()[1] - tt)
      
      T[ep / tc] <- (proc.time()[1] - start.time) - time.eval
   }

    
  }
  
  list(R = R, T = T)
}	



mc.emsf <- function(cells, m, alpha = 1e-1, beta = NULL, 
                    num.episodes = 5e3, epsilon = 0.15, tc = 10, df = 0.999, max.iter = 300, 
                    test.set = rbind(gp(seq(-1,0.15,l=5), seq(-0.07,0.02,l=5)), c(0,0)),
                    max.steps = 300
                   ) 
{
 
  n <- cells[1] * cells[2] + 1
  na <- 3

  if (is.null(beta)) beta <- alpha * 1/na
  
  D  <- array(runif(n * m * na), c(n, m, na))
  for (u in 1:na) D[,,u] <- D[,,u] / apply(D[,,u], 1, sum)
  
  K  <- matrix(runif(m * n), m, n)
  K <- K / apply(K, 1, sum)
  
  D[n,,] <- 0
  D[n,m,] <- 1
  K[m,] <- 0
  K[m, n] <- 1
  
  
  rp <- array(0, n) 
  rc <- array(0, n) # counting

  pi <- sample(1:na, n, TRUE)
     
  R <- array(0, num.episodes %/% tc)
  T <- array(0, num.episodes %/% tc)
  
  start.time <- proc.time()[1]
  time.eval <- 0
  err <- NULL
  for (ep in 1:num.episodes)
  {
    
#     s <- c(runif(1, -1.2, 0.4), runif(1,-0.07, 0.07)) #t est.set[sample(1:nrow(test.set), 1), ]
    s <- test.set[sample(1:nrow(test.set), 1), ]
    j <- mc.state.index(s, cells)
    
    t <- NULL
    st <- 1
    g <- FALSE
    while (!g && st <= max.steps)
    {
     
      i <- j
      if (runif(1) < epsilon) a <- sample(1:na, 1)
      else a <- pi[i]
      
      t <- mc.transition(s, a)
     
      s <- t$s
      j <- mc.state.index(s, cells)
      
      w <- D[i,,a] * K[,j]
      w <- w / sum(w)
      
      D[i,,a] <- (1 - alpha) * D[i,,a] + alpha * w
      K[,j] <- (1 - beta) * K[,j] + beta * w
      
      rp[j] <- rp[j] + t$r
      rc[j] <- rc[j] + 1
 
      g <- t$g
      st <- st + 1
      
#       s <- c(runif(1, -1.2, 0.4), runif(1,-0.07, 0.07)) #t est.set[sample(1:nrow(test.set), 1), ]
#       j <- mc.state.index(s, cells)
      
    }
 
 
   if (ep %% tc == 0)
   {
      
      K <- K / apply(K, 1, sum)
      
#       print(K)
      
      rt <- array(0, n)
      for (i in 1:length(rc)) if (rc[i] > 0) rt[i] <- rp[i] / rc[i]
      rb <- K %*% rt

     
#       plot(rt)
#       matplot(cbind(D[,,1] %*% rb, D[,,2] %*% rb, D[,,3] %*% rb))
      
      Y <- pisf(D, K, rb, df = df, max.iter = max.iter)
      pi <- Y$pi
#       
#       V <- apply(Y$Q, 1, which.max)
#       V <- V[-length(V)]
#       V <- matrix(V, sqrt(length(V)), sqrt(length(V)))
#       image(V)
#       
      tt <- proc.time()[1]
      R[ep / tc] <- mc.evaluate.policy(pi, cells, test.set, max.steps)
      print(R[ep / tc] )
      err <- c(err, R[ep / tc])
      plot(err)
      time.eval <- time.eval + (proc.time()[1] - tt)
      
      T[ep / tc] <- (proc.time()[1] - start.time) - time.eval
   }
   
  }
  
  list(R = R, T = T)
}



mc.emsf.full <- function(cells, m, alpha = 1e-1, beta = NULL, 
                    num.episodes = 5e3, epsilon = 0.15, tc = 10, df = 0.999, max.iter = 300, ne = 1e6,
                    test.set = rbind(gp(seq(-1,0.15,l=5), seq(-0.07,0.02,l=5)), c(0,0)),
                    max.steps = 300
                   ) 
{
  if (is.null(beta)) beta <- alpha * 0.5
 
  n <- cells[1] * cells[2] + 1
  na <- 3
     
  D  <- array(runif(n * m * na), c(n, m, na))
  for (u in 1:na) D[,,u] <- D[,,u] / apply(D[,,u], 1, sum)
  
  K  <- matrix(runif(m * n), m, n)
  K <- K / apply(K, 1, sum)
  
  D[n,,] <- 0
  D[n,m,] <- 1
  K[m,] <- 0
  K[m, n] <- 1
  
  Dh <- array(0, c(n, m, na))
  Kh  <- matrix(0, m, n)
  Dh[n,m,] <- 1
  Kh[m, n] <- 1
  
  rp <- array(0, n) 
  rc <- array(0, n) # counting

  pi <- sample(1:na, n, TRUE)
     
  R <- array(0, num.episodes %/% tc)
  T <- array(0, num.episodes %/% tc)
  
  start.time <- proc.time()[1]
  err <- NULL
  time.eval <- 0
  for (ep in 1:num.episodes)
  {
    
    s <- c(runif(1, -1.2, 0.4), runif(1,-0.07, 0.07)) #t est.set[sample(1:nrow(test.set), 1), ]
#     s <- test.set[sample(1:nrow(test.set), 1), ]
    j <- mc.state.index(s, cells)
    
    t <- NULL
    st <- 1
    g <- FALSE
    while (!g && st <= max.steps)
    {
     
      i <- j
      if (runif(1) < epsilon) a <- sample(1:na, 1)
      else a <- pi[i]
      
      t <- mc.transition(s, a)
     
      s <- t$s
      j <- mc.state.index(s, cells)
      
      w <- D[i,,a] * K[,j]
      w <- w / sum(w)
      
      Dh[i,,a] <- Dh[i,,a] + w
      Kh[,j] <- Kh[,j] + w
      
      rp[j] <- rp[j] + t$r
      rc[j] <- rc[j] + 1
 
      g <- t$g
      st <- st + 1
      
#       s <- c(runif(1, -1.2, 0.4), runif(1,-0.07, 0.07)) #t est.set[sample(1:nrow(test.set), 1), ]
#       j <- mc.state.index(s, cells)
      
    }
 
 
   if (ep %% tc == 0)
   {
      
      for (u in 1:na) 
      {
        for (i in 1:nrow(Dh[,,u])) 
        {
          ss <- sum(Dh[i,,u])
          if (ss > 0) D[i,,u] <- (1 - alpha) *  D[i,,u] + alpha * Dh[i,,u] / ss
        }
      }
      
      K <- (1 - alpha) * K + alpha * Kh / apply(Kh, 1, sum)
      
      rt <- array(0, n)
      for (i in 1:length(rc)) if (rc[i] > 0) rt[i] <- rp[i] / rc[i]
      rb <- K %*% rt
      
#       image(D[1:(nrow(D)-1),,3] %*% K)
#       print(apply(D[,,1], 1, sum))
#     
      
#       plot(rt)
#       matplot(cbind(D[,,1] %*% rb, D[,,2] %*% rb, D[,,3] %*% rb))
      
      Y <- pisf(D, K, rb, df = df, max.iter = max.iter)
      pi <- Y$pi
      
#       V <- apply(Y$Q, 1, which.max)
#       V <- matrix(V, sqrt(length(V)), sqrt(length(V)))
#       image(V)
      
      tt <- proc.time()[1]
      R[ep / tc] <- mc.evaluate.policy(pi, cells, test.set, max.steps)
      print(R[ep / tc] )
      err <- c(err, R[ep / tc])
      plot(err)
      time.eval <- time.eval + (proc.time()[1] - tt)
      
      T[ep / tc] <- (proc.time()[1] - start.time) - time.eval

      Dh[,,] <- 0
      Kh[,]  <- 0
      Dh[n,m,] <- 1
      Kh[m, n] <- 1
  }
   
  }
  
  list(R = R, T = T)
}


emsf <- function(C, D, K, num.it)
{
  
  n <- dim(D)[1]
  m <- dim(D)[2]
  na <- dim(D)[3]

  
#   P <- C
#   for (u in 1:na) P[,,u] <- P[,,u] / apply(P[,,u], 1, sum)
#   
  
  for (it in 1:num.it)
  {
    
    Dh <- array(0, c(n, m, na))
    Kh <- matrix(0, m, n)
    
    for (u in 1:na) 
    {   
      for (i in 1:(n-1)) ## attention to '-1'
      {
        for (j in 1:m)
        {
          
          if (C[i,j,u] > 0)
          {
            w <- D[i,,u] * K[,j]
#             print(D[i,,u])
            w <- C[i,j,u] * w / sum(w)
            Dh[i,,u] <- Dh[i,,u] + w
            Kh[, j] <- Kh[,j] + w
          }
        }
      }
    }
       
   
#     for (u in 1:na) Dh[,,u] <- Dh[,,u] / apply(Dh[,,u], 1, sum)
#     K <- K / apply(K, 1, sum)

    for (u in 1:na) 
    {
      for (i in 1:(n-1)) ## attention to '-1'
      {
        ss <- sum(Dh[i,,u])
        if (ss > 0) Dh[i,,u] <- Dh[i,,u] / ss
      }
    }

    for (i in 1:m)
    {
      ss <- sum(Kh[i,])
      if (ss > 0) Kh[i,] <- Kh[i,] / ss
    }
    
    D <- Dh
    D[n,m,] <- 1

    K <- Kh
    K[m,n] <- 1
    
#     sse <- 0
#     for (u in 1:na) 
#     {
# #       print(paste("u", u))
# #       print("P")
# #       print(max(P[-n,n,u]))
#       sse <- sse + sum((D[,,u] %*% K - P[,,u])^2)
# #       print("DK")
# #       print(round(D[,,u] %*% K[,n], d = 2))
# #       print("")
#     }
#     print(paste("sse", sse))

  }
  
  list (D = D, K = K)
}


mc.emsf.comp <- function(cells, m, num.it, num.episodes = 5e3, epsilon = 0.15, tc = 10, df = 0.999, max.iter = 300,
  test.set = rbind(gp(seq(-1,0.15,l=5), seq(-0.07,0.02,l=5)), c(0,0)),
  max.steps = 300)  
{
  
  n <- cells[1] * cells[2] + 1
  na <- 3
  
  D  <- array(runif(n * m * na), c(n, m, na))
  for (u in 1:na) D[,,u] <- D[,,u] / apply(D[,,u], 1, sum)
  
  K  <- matrix(runif(m * n), m, n)
  K <- K / apply(K, 1, sum)

  D[n,,] <- 0
  D[n,m,] <- 1
  K[m,] <- 0
  K[m, n] <- 1

  C <- array(1, c(n, n, na)) # counting ## THINK
  C[n,,]  <- 0
  C[n,n,] <- 1

  rp <- array(0, n) 
  rc <- array(0, n) # counting

  pi <- sample(1:na, n, TRUE)
     
  R <- array(0, num.episodes %/% tc)
  T <- array(0, num.episodes %/% tc)
  
  start.time <- proc.time()[1]
  time.eval <- 0
  err <- NULL
  for (ep in 1:num.episodes)
  {
    
    s <- test.set[sample(1:nrow(test.set), 1), ]
    j <- mc.state.index(s, cells)
    
    t <- NULL
    st <- 1
    g <- FALSE
    while (!g && st <= max.steps)
    {
     
      i <- j
      if (runif(1) < epsilon) a <- sample(1:na, 1)
      else a <- pi[i]
      
      t <- mc.transition(s, a)
     
      s <- t$s
      j <- mc.state.index(s, cells)
      
      C[i,j,a] <- C[i,j,a] + 1
      rp[j] <- rp[j] + t$r
      rc[j] <- rc[j] + 1
      
      g <- t$g
      st <- st + 1
    }
    
   if (ep %% tc == 0)
   {
      
      SF <- emsf(C, D, K, num.it)
      
      D <- SF$D
      K <- SF$K
      
      rt <- array(0, n)
      for (i in 1:length(rc))
      {
        if (rc[i] > 0) rt[i] <- rp[i] / rc[i]
      }
      rb <- K %*% rt 
      
#       print(D)
#       pd(K)
      
      Y <- pisf(D, K, rb, df = df, max.iter = max.iter)
      pi <- Y$pi
#             
#       V <- apply(Y$Q, 1, which.max)
#       V <- V[-length(V)]
#       V <- matrix(V, sqrt(length(V)), sqrt(length(V)))
#       image(V)
      
      tt <- proc.time()[1]
      R[ep / tc] <- mc.evaluate.policy(pi, cells, test.set, max.steps)
      print(R[ep / tc])
#       err <- c(err, R[ep / tc])
#       plot(err)
      time.eval <- time.eval + (proc.time()[1] - tt)
      
#       P <- C
#       sse <- 0
#       for (u in 1:dim(P)[3]) 
#       {
#         P[,,u] <- P[,,u] / apply(P[,,u], 1, sum)
# #       print(paste("u", u))
# #       print(paste(max(P[-n,n,u]), which.max(P[-n,n,u])))
#         sse <- sse+ sum((D[,,u] %*% K - P[,,u])^2)
#       }
#       
#       err <- c(err, sse)
#       plot(err)
#       image(C[,,1])
      
      T[ep / tc] <- (proc.time()[1] - start.time) - time.eval
   }

    
  }
  
  list(R = R, T = T)
} 



print("mountain.car.experiments.R loaded")