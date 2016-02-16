source("data.plot.R")

euc.dist <- function(x,y) sum((x - y)^ 2)

local.sf <- function(n = 20, taus = c(100, 10, 1), thr = 0.1)
{

   set.par()
   
   A <- matrix(runif(2 * n), n, 2)
   d <- matrix(runif(n), 1, n)
   d <- d / sum(d)
   y <- matrix(d %*% A, 1, 2)
   
   for (i in 1:length(taus))
   {
      d <- matrix(0, 1, n)
      
      for (j in 1:ncol(d)) d[1,j] <- exp(-euc.dist(y, A[j,]) / taus[i])
      d <- d / sum(d)
     
      B <- d %*% A
      
      x11()
      
      plot(A, pch = get.pch(1)[1], col = "BLACK", bg = "BLACK", xlab = NA, ylab = NA, cex = 1.5)
      dgt <- d > thr
      C <- A[dgt, , drop = FALSE]
      ind <- chull(C)
      polygon(C[ind,1], C[ind,2], col = "GREY", border = NA, cex = 1.5)
      points(A[-dgt,], pch = get.pch(1)[1], col = "BLACK", bg = "BLACK", cex = 1.5)
      points(C[ind,], pch = get.pch(1)[1], col = "GREEN", bg = "GREEN", cex = 1.5)
      
      points(B, pch = get.pch(3)[3], col = "RED", bg = "RED", cex = 1.5)
      points(y, pch = get.pch(4)[4], col = "BLUE", bg = "BLUE", cex = 1.5)
   }

}

print("pisf.jair.R loaded")