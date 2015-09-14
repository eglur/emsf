## new version in the file emsf.R

generate.stochastic.matrix <- function(nrows, ncols)
{
   A <- matrix(runif(nrows * ncols), nrows, ncols)
 
# #    This block is to make things more structured
#    for (i in 1:nrows) 
#    {
#       ind <- sample(1:ncols, 2)
#       for (j in ind) A[i,j] <- 1000
#    }
      
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
   list(P = md$P, pi = md$pi, mu = md$mu, y = dt$y, a = dt$a, n = n, m = m, na = na)
}
    

em.sf <- function(n, m, na, P, OBS, ACT, mu, pi, eps = 1e-20, max.it = 10)
{
# Assumes that each observation and each action has appeared at least once
#     n <- length(unique(unlist(OBS)))
#     na <- length(unique(unlist(ACT)))
   
    D <- generate.stochastic.matrices(n,m,na)
    K <- generate.stochastic.matrices(m,n,na)
    

    old.score <- -Inf 
    score <- 0
    it <- 0
    
    error <- NULL
       
    while (abs(score - old.score) > eps && it < max.it)
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
        
        score <- score + sum(log(NF)) + log(A[T-1,] %*% K[, y[T], a[T-1]])
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
      it <- it + 1       
   }
        
    list(D = D, K = K)
}
    

em.sf2 <- function(n, m, na, P, OBS, ACT, eps = 1e-20, max.it = 10)
{
# Assumes that each observation and each action has appeared at least once
#     n <- length(unique(unlist(OBS)))
#     na <- length(unique(unlist(ACT)))
   
    D <- generate.stochastic.matrices(n,m,na)
    K <- generate.stochastic.matrices(m,n,na)
    
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
    
#     for (a in 1:na) C[,,a] <- C[,,a] / apply(C[,,a], 1, sum) ## REMOVE
    
    Pt <- array(0, c(n,n,na))

    for (a in 1:na) Pt[,,a] <- C[,,a] / apply(C[,,a], 1, sum)
#     baseline.ll <- log.likelihood(Pt, OBS, ACT)
    baseline.ll <- frobenius(P, Pt)
    
    old.score <- eps + 1 
    score <- 0
    it <- 0
    
    error <- NULL
    error2 <- NULL
       
    while (abs(score - old.score) > eps && it < max.it) 
    {
        
        old.score <- score 
        score <- 1
        
        for (a in 1:na)
        {
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

                        score <- score +  C[i,j,a] * log(g)
                        
                    }
                }
            }
            
            Dh <- Dh / apply(Dh, 1, sum)
            Kh <- Kh / apply(Kh, 1, sum)
            
            D[,,a] <- Dh
            K[,,a] <- Kh

            # element-by-element
            Dh <- matrix(0, n, m)
            Kh <- matrix(0, m, n)
#             PC <- C[,,a] / apply(C[,,a], 1, sum)
            PP <- D2[,,a] %*% K2[,,a]
            
            for (i in 1:n)
            {
                for (j in 1:m)
                {
                     Dh[i,j] <- D2[i,j,a] * sum(K2[j,,a] * C[i,,a] / PP[i,]) 
                }
            }
            
            for (i in 1:m)
            {
                for (j in 1:n)
                {
                     Kh[i,j] <- K2[i,j,a] * sum((D2[,i,a] * C[,j,a]) / PP[,j]) 
                }
            }

            Dh <- Dh / apply(Dh, 1, sum)
            Kh <- Kh / apply(Kh, 1, sum)
            
            D2[,,a] <- Dh
            K2[,,a] <- Kh
            
        }
      

        for (a in 1:na) Pt[,,a] <- D[,,a] %*% K[,,a]
#         e <- log.likelihood(P, OBS, ACT) # sum((D[,,i] %*% K[,,i] - P[,,i])^2)
        
        e <- frobenius(P, Pt)
        error <- c(error, e)
        
        Pt2 <- array(0, c(n, n, na))
        for (a in 1:na) Pt2[,,a] <- D2[,,a] %*% K2[,,a]

        e <- frobenius(P, Pt2) 
        error2 <- c(error2, e)
        
        matplot(cbind(error, error2), t="l")
        
        print(error2[length(error2)])
        
        it <- it + 1

   }
        
    list(D = D, K = K)
}


is.emsf.nmf <- function(n, m, na, P, OBS, ACT, eps = 1e-20, max.it = 10)
{
# Assumes that each observation and each action has appeared at least once
#     n <- length(unique(unlist(OBS)))
#     na <- length(unique(unlist(ACT)))
   
    D <- generate.stochastic.matrices(n,m,na)
    K <- generate.stochastic.matrices(m,n,na)

    D2 <- D
    K2 <- K
    
    D3 <- D2
    K3 <- K2
    
    tau <- 0
    C <- array(0, c(n,n,na))
    for (batch in 1:length(OBS))
    {
    
        y <- OBS[[batch]]
        a <- ACT[[batch]]

        T <- length(y)
        tau <- tau + T
        
        for (i in 1:(T-1)) C[y[i],y[i+1],a[i]] <- C[y[i],y[i+1],a[i]] + 1 ## counting
        
    }
    
#     for (a in 1:na) C[,,a] <- C[,,a] / apply(C[,,a], 1, sum) ## REMOVE
    
    Pt <- array(0, c(n,n,na))
    for (a in 1:na) Pt[,,a] <- C[,,a] / apply(C[,,a], 1, sum)

    Jt <- array(0, c(n,n,na))
    for (a in 1:na) Jt[,,a] <- C[,,a] / (tau - 1)

      
    #     baseline.ll <- log.likelihood(Pt, OBS, ACT)
    baseline.ll <- frobenius(P, Pt)
    
    old.score <- eps + 1 
    score <- 0
    it <- 0
    
    error <- NULL
    error2 <- NULL
    error3 <- NULL
       
    while (it < max.it)  ## abs(score - old.score) > eps && 
    {
        
        old.score <- score 
        score <- 1
        
        for (a in 1:na)
        {
            # EMSF
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

                        score <- score +  C[i,j,a] * log(g)
                        
                    }
                }
            }
            
            Dh <- Dh / apply(Dh, 1, sum)
            Kh <- Kh / apply(Kh, 1, sum)
            
            D[,,a] <- Dh
            K[,,a] <- Kh
        
#             print("SUMS")
#             print(apply(D[,,a], 1, sum))
#             print(apply(K[,,a], 1, sum))
            
#          NMF (Lee and Seung)
            ## KL divergence (Lee and Seung)
            Dh <- D2[,,a]
            Kh <- K2[,,a]
            

            P2 <- Dh %*% Kh
            for (i in 1:m)
            {
                  for (j in 1:n)
                  {
                        Kh[i,j] <- K2[i,j,a] * sum(D2[,i,a] * Jt[,j,a] / P2[,j]) / sum(D2[,i,a]) 
                  }
            }
                  
            for (i in 1:n)
            {
                  for (j in 1:m)
                  {
                        Dh[i,j] <- D2[i,j,a] * sum(K2[j,,a] * Jt[i,,a] / P2[i,]) / sum(K2[j,,a]) 
                  }
            }
            D2[,,a] <- Dh
            K2[,,a] <- Kh
            
#             print("SUMS")
#             print(apply(D2[,,a], 1, sum))
            print(apply(K2[,,a], 1, sum))
            
            ## Frobenius (Lee and Seung)
            P3 <- D3[,,a] %*% K3[,,a]
            
            K3[,,a] <- K3[,,a] * (t(D3[,,a]) %*% Pt[,,a]) / (t(D3[,,a]) %*% P3)
            D3[,,a] <- D3[,,a] * (Pt[,,a] %*% t(K3[,,a])) / (P3 %*% t(K3[,,a]))
            
#             print("SUMS")
#             print(apply(D3[,,a], 1, sum))
#             print(apply(K3[,,a], 1, sum))
#             
        }
        
        Pe <- array(0, c(n,n,na))
        for (a in 1:na) Pe[,,a] <- D[,,a] %*% K[,,a]
           
        Pn <- array(0, c(n,n,na))
        for (a in 1:na) 
        {
           Pn[,,a] <- D2[,,a] %*% K2[,,a]
           Pn[,,a] <- Pn[,,a] / apply(Pn[,,a], 1, sum)
        }
        
        Pf <- array(0, c(n,n,na))
        for (a in 1:na) 
        {
           Pf[,,a] <- D3[,,a] %*% K3[,,a]
           Pf[,,a] <- Pf[,,a] / apply(Pf[,,a], 1, sum)
        }

        #         e <- log.likelihood(P, OBS, ACT) # sum((D[,,i] %*% K[,,i] - P[,,i])^2)
#       e <- frobenius(P, Pe)
        e <- kld(P, Pe)
        error <- c(error, e)

#       e <- frobenius(P, Pn)
        e <- kld(P, Pn)
        error2 <- c(error2, e)

#       e <- frobenius(P, Pf)
        e <- kld(P, Pf)
        error3 <- c(error3, e)
        
        bb <- 1
        if (it > 100) bb <- 100
        
       #         plot(error,t="l")
        matplot(cbind(error[bb:(it+1)], error2[bb:(it+1)], error3[bb:(it+1)]), t="l")
        
#         print(score)
        
        it <- it + 1

   }
        
    list(D = D, K = K)
}



full.on.line.emsf <- function(P, m, max.it = 1000)
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



log.likelihood <- function(P, OBS, ACT)
{
   
   ll <- 0
   for (batch in 1:length(OBS))
   {
    
      y <- OBS[[batch]]
      a <- ACT[[batch]]
      T <- length(y)
      
      for (t in 1:(T-1))
      {
   #       print(log(P[y[t], y[t+1], a[t]]))
         ll <- ll + log(P[y[t], y[t+1], a[t]])
      }
   }
   ll
}



frobenius <- function(P, Pt)
{
#    print("P")
#    print("----------")
#    print("")
#    print(round(P, d =1))
#    print("Pt")
#    print("----------")
#    print("")
#    print(round(Pt, d =1))
   sse <- 0
   for (a in 1:dim(P)[3]) sse <- sse + sum((P[,,a] - Pt[,,a])^2)
   sse
}



kl <- function(p, pt) 
{
  err <- 0
  for (i in 1:length(p)) err <- err + p[i] * log(p[i] / pt[i])
  err
}


kld <- function(P, Pt)
{
   err <- 0
   for (a in 1:dim(P)[3]) 
   {
      for (i in 1:nrow(P[,,a])) err <- err + kl(P[i,,a], Pt[i,,a])
   }
   err
}


run.experiment <- function(n = 10, m = 3, na = 2, T = 1000, num.batches = 10, max.it = 30)
{
    M <- generate.model.and.data(n, m, na, T, num.batches)

    print("Running original version of EMSF...")
    em.sf(n, na, M$P, M$y, M$a, m, M$mu, M$pi, max.it = max.it)

    x11()
    print("Running new version of EMSF...")
    em.sf2(n, na, M$P, M$y, M$a, m, max.it = max.it)
    
    print("Done.")
}



print("em.sf.R loaded") 
    