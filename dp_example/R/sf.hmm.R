require(Matrix)

## CONCLUSIONS SO FAR
# - Given Pb and K compute D as
# 1. D <- t(K) 
# 2. for (i in 1:ncol(D)) D[,i] <- D[,i] * pib[i]
# 3. D <- D / apply(D,1,sum)
# - Then, DK will give the transition probabilities between "observations" s_i
# - However, in general DK is different from Pb
# - Here's the good part: let Pb2 = KD; then Pb2 has the same stationary distribution as Pb
# - Then, Pb2 is the "true" reduced form of P
# - One way of looking at this is the following: there are (infinitely?) many matrices Pb
# whose stationary distribution pib satisfies pib K = pi, where pi is the stationary distribution
# of P. However, there is a single one given by KD, where D is computed as above.
# It's like a fixed point!! (see test.fixed.point.idea() below)
# For the "fixed point" Pb2, it seems like P(state | observation, state) = 
# P(state | observation), that is, the current observation retains all the information for the dynamics
# From a sampling perspective, the chain defined by Pb2 is such that if we sample from D, then from K, 
# then from D, and so on, we would observe the same distribution of states and observations (if 
# we "censor" one or another we don't see any difference)


euc.dist <- function(x,y) sum((x - y)^ 2)

local.sf <- function(n = 20, taus = c(100, 10, 1))
{
   A <- matrix(runif(2 * n), n, 2)
   y <- matrix(runif(2), 1, 2)
   
   for (i in 1:length(taus))
   {
      d <- matrix(0, 1, n)
      
      for (j in 1:ncol(d)) d[1,j] <- exp(-euc.dist(y, A[j,]) / taus[i])
      d <- d / sum(d)
      
      B <- d %*% A
      
      x11()
      plot(A)
      points(B, col = "RED")
      points(y, col = "BLUE")
   }

}

exp.nz <- function(n, m) (1 - 0.75^(m-1)/2) * n^2 + 0.75^(m-1)/2 * n


g <- function(x, c = 0, sigma = 1) exp(-(x - c)^2 / 2*sigma)

ng <- function(x1, x2, c = 0, sigma = 1)
g(x1, c, sigma) / (g(x1, c, sigma) + g(x2, c, sigma))

dif <- function(x1, x2, n1, n2, c = 0, sigma = 1)
abs(ng(x1,x2, c, sigma) - ng(n1,n2, c, sigma))

Delta <- function(sigma = 1) exp(-1/(2*sigma^2))

dif2 <- function(x1, n1, c = 0, sigma = 1) abs(g(x1,c,sigma) / 2*Delta(sigma) - 
g(n1, c, sigma) / 2 * Delta(sigma))
 
bound <- function(x1,x2,n1,n2,c=0, sigma=1)
print(paste(dif(x1,x2,n1,n2,c, sigma), dif2(x1,n1,c,sigma)))

stat.dist <- function(P, eps = 1e-6)
{
   pi <- matrix(runif(nrow(P)), 1, nrow(P))
   pi <- pi / sum(pi)
   error <- Inf
   while (error > eps)
   {
      pi2 <- pi      
      pi <- pi2 %*% P 
      error <- sum(abs(pi - pi2))
   }
   pi
}

test.hmm.idea <- function(n, m)
{
   Pb <- matrix(runif(m*m), m, m)
   for (i in 1:nrow(P)) Pb[i, sample(1:ncol(P), 1)] <- 100
   Pb <- Pb / apply(Pb, 1, sum)
   
   K <- matrix(runif(m*n), m, n)
   K <- K / apply(K, 1, sum)
   
   pib <- stat.dist(Pb)
   D <- t(K)
   for (i in 1:ncol(D)) D[,i] <- D[,i] * pib[i]
   
#   D <- matrix(runif(n*m), n, m)
  D <- D / apply(D, 1, sum)
   
   P <- D %*% K
   pi <- stat.dist(P)
   
   pit <- pib %*% K
   
   print(sum(abs(pi - pit)))
   matplot(cbind(t(pi), t(pit)), t="l")
   
   Pb2 <- K %*% D
   pib2 <- stat.dist(Pb2)
   
   print(sum(abs(pib - pib2)))
   matplot(cbind(t(pib), t(pib2)), t="l")
   
   pit2 <- pib2 %*% K
   
   print(sum(abs(pi - pit2)))
   matplot(cbind(t(pi), t(pit2)), t="l")
   
   
}


compute.D <- function(Pc, K, pic)
{
   
   n <- ncol(K)
   m <- nrow(K)
   
   Z <- K
   for (i in 1:ncol(Z)) 
   {
      Z[,i] <- Z[,i] * pic
   }
   
   
   D <- matrix(0, n, m)
   for (i in 1:n)
   {
      for (j in 1:m)
      {
         sum.zi <- sum(Z[,i])
         if (sum.zi > 0)  D[i,j] <- sum(Pc[,j] * Z[,i]) / sum.zi
      }
   }
   
   D
}

generate.stochastic.matrix <- function(nrows, ncols)
{
   A <- matrix(runif(nrows * ncols), nrows, ncols)
   A <- A / apply(A, 1, sum)
   A
}

test.fixed.point.idea <- function(n, m, length = 100)
{
   Pb <- matrix(runif(m*m), m, m)
   Pb <- Pb / apply(Pb, 1, sum)
   
   K <- matrix(runif(m*n), m, n)
   K <- K / apply(K, 1, sum)

   pib <- stat.dist(Pb)
   D <- compute.D(Pb, K, pib)
   
   
   err <- array(0, length)
   i <- 1
   dif <- Inf
   while (dif > 0 && i <= length)
   {
      Pb <- K %*% D
#       pib2 <- stat.dist(Pb)
      D2 <- compute.D(Pb, K, pib)
      dif <- sum(abs(D -  D2))
#       dif2 <- sum(abs(pib -  pib2)) 
      D <- D2
#       pib <- pib2
      err[i] <- dif
      i <- i + 1
      print(paste("D - D2:", dif))
#       print(paste("pib - pib2:", dif2))
   }
  err
}

start.from.sf <- function(n, m, length = 100)
{
   D <- generate.stochastic.matrix(n,m)
   K <- generate.stochastic.matrix(m,n)
   P  <- D %*% K
   
   e <- array(100)
   for (i in 1:100)
   {
      
      Pb <- K %*% D
   
      pib <- stat.dist(Pb)
      D <- compute.D(Pb, K, pib)
      
      e[i] <- sum(abs(D %*% K - P))
      print(e)
      
   }
   e
}


test.relationship.D.and.K <- function(n, m)
{
   K <- matrix(runif(m*n), m, n)
   K <- K / apply(K, 1, sum)

   D <- matrix(runif(n*m), n, m)
   D <- D / apply(D, 1, sum)

   Pb <- K %*% D
   pib <- stat.dist(Pb)
   
   D2 <- t(K)
   for (i in 1:ncol(D2)) D2[,i] <- D2[,i] * pib[i]
   D2 <- D2 / apply(D2, 1, sum)
   
   print(sum(abs(D2 - D)))
   
   print(D)
   print(D2)
}


start.from.hmm <- function(n, m)
{
   Pb <- matrix(m*m, m, m)
   Pb <- Pb / apply(K, 1, sum)
   K  <- matrix(n*m, m, n)
   K <- K / apply(K, 1, sum)
   
   pib <- stat.dist(Pb)
   D <- t(K) 
   for (i in 1:ncol(D)) D[,i] <- D[,i] * pib[i]
   D <- D / apply(D,1,sum)
   
   Pb2 <- K %*% D
   
}


test.proposition <- function(n, m)
{
   stop <- FALSE
   while (!stop)
   {
      K <- matrix(runif(n*m), m, n)
      K <- K / apply(K, 1, sum)
      
      stop <- rankMatrix(K) == m
   }
   
   D <- matrix(runif(n*m), n, m)
   D <- D / apply(D, 1, sum)
   
#    mu <- runif(m)
#    mu <- mu / sum(mu)
   mu <- stat.dist(K%*%D)
   
   A <- t(K)
   for (i in 1:m) A[,i] <- A[,i] * mu[i]
   A <- A / apply(A, 1, sum)
#    print(apply(A, 1, sum))
   
   Pc <- solve(t(A) %*% A) %*% t(A) %*% D
   
   #    print(Pc)
#     print(sum(abs(t(A) %*% D - t(A) %*% A %*% Pc)))
    print(sum(abs(D - A %*% Pc)))

   
}