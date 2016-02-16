require("hmm.discnp")
source("data.plot.R")
source("util.R")

compute.value.function <- function(P, r, gamma = 0.99)
{
  v <- solve(diag(nrow(P)) - gamma * P, r)
}

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

sf.hmm <- function(n, m, k, p, generate.method, nrun = 1)
{
   # generate.method: h = start from HMM; s = start from SF; p = generate P 
   
   if (generate.method != "h" && generate.method != "s" && generate.method != "p") 
   {
     print("Invalid generate method")
     break
   }
     
   if (generate.method == "p") P <- generate.stochastic.matrix(n,n)
   else
   {
     # Generate D and K
    K <- generate.stochastic.matrix(m,n)

    D <- NULL
  
    if (generate.method == "s") D <- generate.stochastic.matrix(n,m)
    else
    {
	Pc <- generate.stochastic.matrix(m,m)
	pic <- stat.dist(Pc)
	D <- compute.D(Pc, K, pic)
    }
    
    P <- D %*% K
   }
     
   pi <- stat.dist(P)
   
   # Generate data 
   o <- NULL
   for (i in 1:k)
   {
  
      zo <- array(0,p)
      s <- sample.from.dist(pi)
      zo[1] <- s
      for (j in 2:p) 
      {
         s <- sample.from.dist(P[s,])
         zo[j] <- s
      }
      
      o <- c(o, list(zo))
   }
   
  # Compute HMM
  
   H <- NULL
   log.like <- -Inf
   
   for (i in 1:nrun)
   {
      H2 <- hmm(o, yval = 1:n, K = m, stationary = TRUE, cis = TRUE, rand.start = list(tpm = TRUE, Rho = FALSE))
      if (H2$log.like > log.like) H <- H2
   }
   
   pih <- H$ispd ## compute pih from T?
   T   <- H$tpm
   O   <- t(H$Rho)
   
   D2 <- compute.D(T, O, pih)
   
   Pt <- D2 %*% O
   
   # Counting
   Ptt <- matrix(0, n, n)
   for (i in 1:length(o))
   {
      s <- o[[i]][1]
      
      for (j in 2:length(o[[i]]))
      {
         sp <- o[[i]][j]
         Ptt[s,sp] <- Ptt[s,sp] + 1
         s <- sp
      }
   }
   Ptt <- Ptt / apply(Ptt, 1, sum)
   
   list(ml = sum((P - Ptt)^2), hmm = sum((P - Pt)^2))
   
   # compare value functions
#    r.tmp <- runif(n)
#    r  <- P %*%   r.tmp
#    rb <- Pt %*%  r.tmp
#    rc <- Ptt %*% r.tmp
#    
#    v <- compute.value.function(P, r)
#    vt <- compute.value.function(Pt, rb)
#    vtt <- compute.value.function(Ptt, rc)
#    
#    c(sum((v - vtt)^2), sum((v - vt)^2))
   
}

run.experiment <- function(num.obs, num.avg, n, m, generate.method = "s", nrun = 1, save = TRUE)
{
    nexp <- length(num.obs)
    res <- array(0, c(nexp, num.avg, 2))
    
    for (i in 1:nexp)
    {
        for (j in 1:num.avg) 
        {
	    res[i,j,1] <- sf.hmm(n, m, 1, num.obs[i], generate.method, nrun)$ml
	    res[i,j,2] <- sf.hmm(n, m, 1, num.obs[i], generate.method, nrun)$hmm
	}
    }
    if (save) save.res(res, num.obs, n, m, generate.method)
    res
}

plot.res <- function(res, x = NULL)
{
  
  R <- matrix(0, dim(res)[1], 2)
  Yui <- R
  Yli <- R
  
  for (i in 1:dim(res)[1])
  {
    R[i,1] <- mean(res[i,,1])
    R[i,2] <- mean(res[i,,2])
    
   
    ci <- sd(res[i,,1]) / sqrt(length(res[i,,1]))
    Yui[i,1] <- R[i,1] + ci
    Yli[i,1] <- R[i,1] - ci
    
    ci <- sd(res[i,,2]) / sqrt(length(res[i,,2]))
    Yui[i,2] <- R[i,2] + ci
    Yli[i,2] <- R[i,2] - ci
  }

  if (is.null(x)) x <- 1:nrow(res)
  mp(x, R, Yui, Yli, ylab = "||P - estimate(P)||", xlab = "Number of observations")
  leg("topright", leg = c("Maximum Likelihood", "HMM"))
}
        
save.res <- function(res, num.obs, n, m, generate.method, filename = "./sf_hmm/mse_P")
{
    if (dim(res)[1] != length(num.obs)) print("Sizes of res and num.obs are inconsistent")
    else
    {
	for (i in 1:length(num.obs)) wt(res[i,,], paste(paste(filename, n, m, num.obs[i], generate.method, sep = "_"),".txt"))
    }      
    print("res saved")  
}

load.res <- function(num.obs, num.avg, n, m, generate.method, filename = "./sf_hmm/mse_P")
{
    
    res <- array(0, c(length(num.obs), num.avg, 2))
    for (i in 1:length(num.obs)) 
    {
        res[i,,] <- as.matrix(read.table(paste(paste(filename, n, m, num.obs[i], generate.method, sep = "_"),".txt")))
    }
    res
}

print("sf_hmm.R loaded")