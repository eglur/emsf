# Analysis of data from the "multi component replacement" task 

source("util.R")
source("ml.R")
source("cr.data.analysis.R")

fit.pi.times <- function(red = FALSE, nc = 2:5, lt = 10,
   dir = "~/cr/experiments/", ...)
{
  suffix <- ""
  if (red) suffix = "red_"
   
   p <- length(nc)
   
   bf <- 2^nc # branching factor
   
   X <- matrix(0, p, 3)
   for (i in 1:p)
   {
     X[i,1] <- bf[i]
     X[i,2] <- bf[i]^2
     X[i,3] <- bf[i]^3
   }
    
   D <- matrix(0, 1, 3)
   m <- 2^(nc[length(nc)] + 1)
   D[1,1] <- m
   D[1,2] <- m^2
   D[1,3] <- m^3
   
  
   K <- matrix(0, 2, 1)
   for (t in c(1,2))
   {
      y<- array(0, p)
      for (i in 1:p)
      {
         T <- read.table(ps(dir, "pi_", suffix, "time_", nc[i], "_", lt, ".txt"))
         y[i] <- T[t,1]
      }
      
      w <- least.squares.params(X, y)
           
      K[t,1] <- least.squares.batch(D, w)
     
      plot(c(nc,m), c(y, K[t,1]), ylab="Seconds", xlab="Number of components", ...)
    }
    
    wt(K, ps(dir, "pi_", suffix, "time_", nc[length(nc)]+1, "_", lt, ".txt"))
}



print("cr.R loaded")