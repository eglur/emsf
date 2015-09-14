source("dist.mat.R")

fuzzy.c.means <- function(X, c, memb.exp=2, max.it = nrow(D), max.dif = 1e-5) {
	D <- matrix(runif(nrow(X)*c), nrow(X), c)
	D <- D / apply(D,1,sum)

   it <- 0
	dif <- Inf
	
	while(it < max.it && dif > max.dif) {
      U <- t(D^memb.exp)
      U <- U / apply(U,1,sum)
      W <- (U %*% X)
      
      B <- dist.mat(X,W)^(2/(memb.exp-1)) 
       
      for (i in 1:nrow(D)) {
      	for (j in 1:ncol(D)) {
      		D[i,j] <- 1/sum(B[i,j]/B[i,])
      		}
      	}

      it <- it + 1
      }
   list(D=D, W=W)
   }
      
print("fuzzy.c.means.R loaded")      
      
   
   

