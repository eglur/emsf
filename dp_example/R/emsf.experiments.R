library(gtools)
source("emsf.R")
source("util.R")
source("data.plot.R")


emsf.experiment.generate.P <- function(n, alpha.dir, na = 1)
{
  P <- array(0, c(n, n, na))
  for (a in 1:na) P[,,a] <- rdirichlet(n,rep(alpha.dir, n))
  P 
}

emsf.experiment.fixed.tc <- function(n, alpha.dir, m, tc, alphas, max.it, na = 1, num.avg = 50, dir = "./files/", idx = "")
{
  
  R <- NULL
  D <- NULL
  
  for (i in 1:length(alphas))
  {
    L <- NULL
    
    for (j in 1:num.avg)
    {
      print(paste("Running alpha = ", alphas[i], "  run", j))
      P <- emsf.experiment.generate.P(n, alpha.dir, na)
      L <- cbind(L, emsf(P, m, tc, tc, alphas[i], max.it))
    }
    
    R <- cbind(R, apply(L, 1, mean))
    colnames(R)[ncol(R)] <- paste(alphas[i])
    
    D <- cbind(D, apply(L, 1, sd) / num.avg)
    colnames(D)[ncol(D)] <- paste(alphas[i])
    
    prefix <- ps(dir, paste("emsf", n, m, tc, max.it, alpha.dir, sep="_"))
    write.table(R, ps(prefix, idx, "_kld.txt"), row.names=FALSE, col.names=TRUE, quote=FALSE)
    write.table(D, ps(prefix, idx, "_ser.txt"), row.names=FALSE, col.names=TRUE, quote=FALSE)
  }
  
  list(R = R, D = D)
}


make.leg.alpha <- function(tc, alpha) {
# makes a legend with labels epsilon = dfs[1, 2, ...]
   l <- expression()
   for (a in tc)
   {
    for(d in alpha) {
      l <- c(l, substitute(expression(t[c] == a2 & alpha == d2), list(a2=a, d2=d))[[2]]) 
      }
   }
   print(l)
 l
}


emsf.experiment.plot.results <- function(alphas = c(0.1, 0.4, 1), num.points = 50)
{
  
  Z  <- as.matrix(read.table("./files/emsf_100_30_1000_1e+06_1_kld.txt", header = TRUE))
  Z2 <- as.matrix(read.table("./files/emsf_100_30_10000_1e+06_1_kld.txt", header = TRUE))
  Z <- Z[seq(1, nrow(Z), length = num.points), ]
  Z2 <- Z2[seq(1, nrow(Z2), length = num.points), ]
  
  D  <- as.matrix(read.table("./files/emsf_100_30_1000_1e+06_1_ser.txt", header = TRUE))
  D2 <- as.matrix(read.table("./files/emsf_100_30_10000_1e+06_1_ser.txt", header = TRUE))
  D <- D[seq(1, nrow(D), length = num.points), ]
  D2 <- D2[seq(1, nrow(D2), length = num.points), ]

  R <- cbind(Z,Z2)
  D <- cbind(D,D2)
  
  mp(seq(1, 2e6, l = nrow(R)), R, R + D, R - D, xlab = expression(tau), ylab = expression("KL"[rho]*"(P, DK)"))
  leg("topright", make.leg.alpha(c(10^3, 10^4), alphas))
}
  


