# This is a (very) simplified version of the file 'detour.nmf.R'
# it has only the functions used in the thesis, and also has the same name
# used for the variables (which makes it easier to understand it)

# Non-negative matrix factorization
#-------------------------------------------------------------------------------
# All these functions perform a non-negative matrix factorization (NPF), that
# is, given a (n x m) matrix A and an integer k, they return two matrices, 
# K (n x k) and D (k x m) so that ||A - LK|| is minimized in some sense (the
# norm ||.|| being  minimized varies over the functions).

source("proj.sum.one.R")
source("rbfn.gaussian.R")
source("dist.mat.R")
source("fuzzy.c.means.R")

nmf.mult <- function(P, m, max.iter = 50, epsilon = -Inf, precision = 1e-10,
                        hist = FALSE, plot = FALSE, P.only=TRUE) {
# Lee and Seung's multiplicative rule (projected)
# This function tries to minimize the Euclidean distance
# In order to avoid a division by zero, P should not have a zero column 
# 'hist' indicates whether the history should be returned
	
      # Every element of K and D should be > 0
      #K <- matrix(runif(m * ncol(P), 0.1, 1), m, ncol(P))
      #K <- K/ apply(K, 1, sum)
      #Another possibility:
      K <- P[sample(1:nrow(P), m), ]
	
      D <- matrix(runif(nrow(P) * m, 0.1, 1), nrow(P), m)
	 D <- D / apply(D, 1, sum)

      h <- NULL
      dif <- Inf
	 it <- 0
	 while (dif > epsilon && it < max.iter) {
		K.old <- K		
		
		# update D
		N <- P %*% t(K)
		Dn <- D %*% K %*% t(K)
          # I added this to avoid division-by-zero:
          #min.d <- max(N) * precision
          ind <- Dn == 0
          N[ind] <- 1 
          Dn[ind] <- 1 
		
          D <- D * N / Dn
		D <- proj.sum.one(D)
		
		# update K
		N <- t(D) %*% P
		Dn <- t(D) %*% D %*% K
          # I added this to avoid division-by-zero:
          ind <- Dn == 0
          N[ind] <- 1 
          Dn[ind] <- 1 
		
          K <- K * N / Dn
		if (P.only) K <- proj.sum.one(K)
          else K <- cbind(K[,1],proj.sum.one(K[,2:ncol(K)]))        
          
		
		# compute changes
		dif <- max(abs(K - K.old))
		it <- it + 1
          
          if (hist) {
            h <- c(h, sse(D%*% K, P))
		  if (plot) plot(h, t="l")
            }        
		}
     
     if (P.only) list(D = D, K = K, hist = h)
     else list(D = D, K = K[,2:ncol(K)], r=K[,1], hist = h)
	}


kl <- function(P,Q) {
# Kullback_leibler divergence
   sum(P * log(P/Q) - P + Q)
   }
   
   
nmf.mult.kl <- function(R, P, k, max.iter = 50, K = NULL, D = NULL,   
                        epsilon = -Inf, precision = 1e-6, hist = FALSE,
                        plot = FALSE) {
# Dee and Seung's multiplicative rule 
# This function tries to minimize the Kullback-Leibler divergence
     
     
     if (is.null(K)) {
          # Every element of K and D should be > 0
          #K <- matrix(runif(k * ncol(P), 0.1, 1), k, ncol(P))
          #K <- K / apply(K, 1, sum)
          K <- P[sample(1:nrow(P), k), ]
          }
          
     if (is.null(D)) {
          D <- matrix(runif(nrow(P) * k, 0.1, 1), nrow(P), k)
          D <- D / apply(D, 1, sum)
          }
   
     dif <- Inf
     it <- 0
     h <- NULL
     while (dif > epsilon && it < max.iter) {
          K.old <- K          
          
           # update K
          Dn <- D %*% K 
          Dn[Dn==0] <- 1 # I added this to avoid a division-by-zero
          A <- P / Dn
          N <- t(D) %*% A
          Dn <- apply(D, 2, sum)
          Dn <- matrix(rep(Dn, ncol(K)), nrow(K), ncol(K), byrow = FALSE)
          Dn[Dn==0] <- 1 # I added this to avoid a division-by-zero
          K <- K * N / Dn
         
          # update D
         Dn <- D %*% K 
          Dn[Dn==0] <- 1 # I added this to avoid a division-by-zero
          A <- P / Dn
	  N <- A %*% t(K)
          Dn <- apply(K, 1, sum)
          Dn <- matrix(rep(Dn, nrow(D)), nrow(D), ncol(D), byrow = TRUE)
          Dn[Dn==0] <- 1 # I added this to avoid a division-by-zero
          D <- D * N / Dn
             
         # compute changes
          dif <- max(abs(K - K.old))
          it <- it + 1
          
          if (hist) {
            h <- c(h, kl(D%*% K, P))
            if (plot) plot(h, t="l")
            }        
          }
     
     
     K <- proj.sum.one(K)
     D <- proj.sum.one(D)
     
     DS <- svd(D)
     # remove the small singular values
     d.min <- max(DS$d) * precision
     DS$d[DS$d < d.min] <- 0
     non.zero.ind <- DS$d != 0
     DS$d[non.zero.ind] <- 1 / DS$d[non.zero.ind]
     DS$d <- diag(DS$d)
     # solve the system
     R2 <- DS$v %*% DS$d %*% t(DS$u) %*% R
     
     list(D = D, K = K, r = R2, hist = h)
     }


nmf.kmeans <- function(P, m, max.iter = 10, P.only = TRUE, ...) {
# K-means factorization
## o certo aqui seria normalizar no caso em que P.only == FALSE
   A <- kmeans(P, m, iter.max = max.iter, ...)
   D <- matrix(0, nrow(P), m)
   for (i in 1:nrow(D)) D[i,A$cluster[i]] <- 1
   if (P.only) list(D= D, K = A$centers)
   else list(D = D, K = A$centers[,2:ncol(A$centers)], r = A$centers[,1])
   }
   


nmf.fast.heuristic <- function(P, m, max.iter = 10, P.only = TRUE, n = nrow(P), 
                              dif.min = 0, ...) {
# Fast
   K <- NULL
   D <- NULL

   max <- apply(P,2,which.max)
   min <- apply(P,2,which.min)
   keep <- array(FALSE, length(max))
   for (i in 1:length(keep)) {
     if ((P[max[i],i] - P[min[i],i]) > dif.min) keep[i] <- TRUE
     }

   sel <- unique(c(max[keep],min[keep]))

   if (length(sel) >= n) { # Nothing happens
     sel <- 1:nrow(P)
     D <-  diag(nrow(P))
     }
   else {
     K <- P[sel,]
     D <- matrix(0, nrow(P), nrow(K))
     D[sel,1:length(sel)] <- diag(length(sel))
     D[-sel,] <- proj.sum.one(P[-sel,] %*% t(K) %*% solve(K %*% t(K)))
     }
   if (P.only) list(D= D, K = P[sel,])
   else list(D = D, K = P[sel,2:ncol(P)], r = P[sel,1])
   }


nmf.fast.heuristic2 <- function(P, m, max.iter = Inf, P.only = TRUE, n =
                              nrow(P), dif.min = 0, ...) {
# Deal with rewards!!!
   K <- matrix(0, ncol(P), ncol(P))

   max <- apply(P,2,max)
   remain <- array(1, nrow(K))
   used <- array(FALSE, nrow(K))

   for (col in 1:nrow(K)) {
      ddif <- max(max[col] - dif.min, 0)
      if (ddif > 0) {
          remain2 <- remain - ddif
          remain2[remain2 < 0] <- Inf
          row <- which.min(remain2)
          K[row,col] <- ddif 
          remain[row] <- remain2[row]
          used[row] <- TRUE
          }
      }

    K <- K[used,]
    remain <- remain[used]
    for (i in 1:nrow(K)) {
      cols <- K[i,] > 0
      K[i, cols] <- K[i, cols] + remain[i] / sum(cols)
      }
 
    D <- proj.sum.one(P %*% t(K) %*% solve(K %*% t(K)))
    
   if (P.only) list(D= D, K = K)
   else list(D = D, K = K[,2:ncol(K)], r = K[,1])
   }



nmf.fcm <- function(P, m, max.iter = 10, P.only = TRUE, memb.exp = 1.2,...) {
# Fuzzy C-means factorization
   A <- fuzzy.c.means(P, m, memb.exp, max.it = max.iter)
   if (P.only) list(D= A$D, K = A$W)
   else list(D = A$D, K = A$W[,2:ncol(A$W)], r = A$W[,1])
   }
   

nmf.fanny <- function(P, m, max.iter = 10, P.only = TRUE, memb.exp = 1.01,...) {
# Factorization using the FANNY algorithm
   A <- fanny(P, m, maxit = max.iter, memb.exp=memb.exp, stand=TRUE,
   			 cluster.only = TRUE, keep.diss=FALSE, keep.data=FALSE, ...)
   D <- A$membership
   U <- t(D^memb.exp)
   P <- (U %*% P) / apply(U,1,sum)
   
   if (P.only) list(D= D, K = P)
   else list(D = D, K = P[,2:ncol(P)], r = P[,1])
   }
   

   
nmf.kmeans.prop <- function(P, m, max.iter = 10, P.only = TRUE, tau = 0.2,...){
# K-means factorization with proportional assignment 
      A <- kmeans(P, m, iter.max = max.iter, ...)
      D <- matrix(0, nrow(P), m)
      #uses a RBF network to assign states to archetypes
      rbfn <- make.null.rbfn()
      rbfn$c <- A$centers
      rbfn$s <- rep(tau, nrow(rbfn$c))
      D <- rbfn.norm.design.matrix(rbfn, P)
      
      if (P.only) list(D= D, K = A$centers)
      else list(D = D, K = A$centers[,2:ncol(A$centers)], r = A$centers[,1])
      }
   


nmf.pam <- function(P, m, max.iter = 10, P.only = TRUE, ...) {
# Factorization using PAM algorithm
## tá normalizando...
   A <- pam(P, m, metric= "euclidean", keep.data = FALSE, cluster.only = FALSE,
               stand = TRUE, ...)
   D <- matrix(0, nrow(P), m)
   for (i in 1:nrow(D)) D[i,A$clustering[i]] <- 1
     
   K <- P[A$id.med,]
   if (P.only) list(D= D, K = K)
   else list(D = D, K = K[,2:ncol(K)], r = K[,1])
   }


nmf.clara <- function(P, m, max.iter = 10, P.only = TRUE, ...) {
# Factorization using CLARA algorithm
## tá normalizando...
   A <- clara(P, m, metric= "euclidean", keep.data = FALSE, medoids.x=FALSE,
               stand = TRUE, ...)
   D <- matrix(0, nrow(P), m)
   for (i in 1:nrow(D)) D[i,A$clustering[i]] <- 1
     
   K <- P[A$i.med,]
   if (P.only) list(D= D, K = K)
   else list(D = D, K = K[,2:ncol(K)], r = K[,1])
   }

      
nmf.S.D <- function(P, m, max.iter = 10, P.only = TRUE, tau=0.1, S, ...){
# proportional assignment performed in S and K determined by clara
## tá normalizando...
   A <- clara(P, m, metric= "euclidean", keep.data = FALSE, medoids.x=FALSE,
   stand = TRUE, ...)
   K <- P[A$i.med,]
   
   #uses a RBF network to assign states to archetypes
   rbfn <- make.null.rbfn()
   s <- A$i.med %% nrow(S)
   s[s==0] <- nrow(S)
   rbfn$c <- S[s,]
   rbfn$s <- rep(tau, nrow(rbfn$c))
   
   S <- rbind(S,S,S)
   D <- rbfn.norm.design.matrix(rbfn, S)
   
   if (P.only) list(D= D, K = K)
   else list(D = D, K = K[,2:ncol(K)], r = K[,1])
   }

    
nmf.S.DK <- function(P, m, max.iter = 10, P.only = TRUE, tau=0.1, S, ...){
# D and K determined in S 
   num.actions <- 3 # mountain
   num.states <- nrow(P) %/% num.actions
         
   k <- m %/% num.actions 
   m <- k * num.actions
   ind <- select.grid(S, k, "euclidean")
   K <- matrix(0, m, ncol(P))
   for (a in 1:num.actions) {
         b <- (a-1) * k + 1
         e <- b + k - 1
         delta <- (a-1) * num.states
         K[b:e,] <- P[ind + delta,]   
         }

   d <- as.matrix(dist(normalize(S)))
   Do <- alloc.prop(d[,ind],tau)
   D <- matrix(0, nrow(P), m)
   for (a in 1:num.actions) {
         b <- (a-1) * k + 1
         e <- b + k - 1
         b2 <- (a-1) * num.states + 1
         e2 <- b2 + num.states - 1
         D[b2:e2,b:e] <- Do 
         }
         
   if (P.only) list(D= D, K = K)
   else list(D = D, K = K[,2:ncol(K)], r = K[,1])
   }

   
nmf.ls <- function(P, k, max.iter = 50, P.only = TRUE, epsilon = 1e-6, precision
= 1e-6) {
# Non-negative projected least-squares
# No convergence guarantees
   
   K <- P[sample(1:nrow(P), k), ]
   D <- matrix(0, nrow(P), k)
   
   # h <- NULL
   dif <- Inf
   it <- 0
   while (it <= max.iter && dif > epsilon) {
      ##  FIX K AND COMPUTE D ##
      # singular value decomposition
      KS <- svd(t(K))
      # remove the small singular values
      d.min <- max(KS$d) * precision
      KS$d[KS$d < d.min] <- 0
      non.zero.ind <- KS$d != 0
      KS$d[non.zero.ind] <- 1 / KS$d[non.zero.ind]
      KS$d <- diag(KS$d)
      
      # solve the system 
      D <- proj.sum.one(t(KS$v %*% KS$d %*% t(KS$u) %*% t(P)))
   
      
      K.old <- K
      
      ##  FIX D AND COMPUTE K ##
      # singular value decomposition
      LS <- svd(D)
      # remove the small singular values
      d.min <- max(LS$d) * precision
      LS$d[LS$d < d.min] <- 0
      non.zero.ind <- LS$d != 0
      LS$d[non.zero.ind] <- 1 / LS$d[non.zero.ind]
      LS$d <- diag(LS$d)
      # solve the system
      K <- LS$v %*% LS$d %*% t(LS$u) %*% P 
      K <- cbind(K[,1], proj.sum.one(K[,2:ncol(K)])) # considering P includes R
      
      # h <- c(h, sse(D%*% K, P))
      # plot(h)
      dif <- max(abs(K - K.old))
      it <- it +1 
      }
   #print(h[length(h)])
   list(D = D, K = K[,2:ncol(K)], r = K[,1])
   }
   

nmf.subset <- function(P, m, max.iter = 10, P.only = TRUE, precision = 1e-6,...) {
	t <- system.time({
 	ind <- sample(1:nrow(P), m)
 	K <- P[ind,]
   
   D <- matrix(runif(m*nrow(P)), nrow(P), m)
   D <- D / apply(D,1,sum)}
   , TRUE)[1]
   print(paste("PS:", t))
   
   if (P.only) list(D= D, K = K)
   else list(D = D, K = K[,2:ncol(K)], r = K[,1])
   }
   

nmf.bucket <- function(P, m, max.iter = 10, P.only = TRUE, 
							precision = 1e-6, sigma = 0.1, ...) {

	
	if (!P.only) {
		R <- P[,1]
		P <- P[,2:ncol(P)]
		}
	
	S <- apply(P, 2, max)
	o <- order(S, decreasing = TRUE)
		
	# compute K
	K <- matrix(0, m, ncol(P))
	for (i in 1:ncol(P)) {
		ind <- i %% m
		if (ind == 0) ind <- m
		K[ind,o[i]] <- S[o[i]]
		}
	K <- K / apply(K,1,sum)
	
	# compute D 
	rbfn <- make.null.rbfn()
	rbfn$c <- K
	rbfn$s <- array(sigma, nrow(K))
	D <- rbfn.norm.design.matrix(rbfn,P)
		
	Rb <- NULL
	if (!P.only) {
      # compute Rb
      Rb <- solve(t(D)%*%D)%*%t(D) %*% R
      }
 
   if (P.only) list(D = D, K = K)
   else list(D = D, K = K, r = Rb[,1])
   }
   

nmf.dist <-function(P, m, max.iter = 10, P.only = TRUE, precision = 1e-10,...){
# Factorization based on the distance between vectors...
   
   D <- dist.mat(P,P)
   D <- apply(D, 1, sum)
   o <- order(D, decreasing = TRUE)
   K <- P[o[1:m],]
   # singular value decomposition
   KS <- svd(t(K))
   # remove the small singular values
   d.min <- max(KS$d) * precision
   KS$d[KS$d < d.min] <- 0
   non.zero.ind <- KS$d != 0
   KS$d[non.zero.ind] <- 1 / KS$d[non.zero.ind]
   KS$d <- diag(KS$d)
   # solve the system 
   D <- proj.sum.one(t(KS$v %*% KS$d %*% t(KS$u) %*% t(P)))
   
	if (P.only) list(D = D, K = K)
   else list(D = D, K = K[,2:ncol(K)], r = K[,1])
   }



nmf.simple <- function(P, m, max.iter = 10, P.only = TRUE, 
							precision = 1e-6, ...) {

	
# 	if (!P.only) {
# 		R <- P[,1]
# 		P <- P[,2:ncol(P)]
# 		}
	
	S <- apply(P, 2, max)
	o <- order(S, decreasing = TRUE)
	I <- apply(P,2,which.max)
	
	# compute K
	K <- P[unique(I)[(o)[1:m]],]
	
   # singular value decomposition
   KS <- svd(t(K))
   # remove the small singular values
   d.min <- max(KS$d) * precision
   KS$d[KS$d < d.min] <- 0
   non.zero.ind <- KS$d != 0
   KS$d[non.zero.ind] <- 1 / KS$d[non.zero.ind]
   KS$d <- diag(KS$d)
   # solve the system 
   D <- proj.sum.one(t(KS$v %*% KS$d %*% t(KS$u) %*% t(P)))
   
# 	Rb <- NULL
# 	if (!P.only) {
#       # compute Rb
#       LS <- svd(D)
#       # remove the small singular values
#       d.min <- max(LS$d) * precision
#       LS$d[LS$d < d.min] <- 0
#       non.zero.ind <- LS$d != 0
#       LS$d[non.zero.ind] <- 1 / LS$d[non.zero.ind]
#       LS$d <- diag(LS$d)
#       # solve the system
#       Rb <- LS$v %*% LS$d %*% t(LS$u) %*% R
#       }
 
   if (P.only) list(D = D, K = K)
   else list(D = D, K = K[,2:ncol(K)], r = K[,1])
   }


nmf.heuristic <- function(P, m, max.iter = 10, P.only = TRUE, 
							precision = 1e-6, ...) {

	
 	if (!P.only) {
 		R <- P[,1]
 		P <- P[,2:ncol(P)]
 		}
	
	S <- apply(P, 2, max)
	o <- order(S, decreasing = TRUE)
	I <- apply(P,2,which.max)
	
	# compute K
	K <- P[unique(I)[(o)[1:m]],]
	
   # singular value decomposition
   KS <- svd(t(K))
   # remove the small singular values
   d.min <- max(KS$d) * precision
   KS$d[KS$d < d.min] <- 0
   non.zero.ind <- KS$d != 0
   KS$d[non.zero.ind] <- 1 / KS$d[non.zero.ind]
   KS$d <- diag(KS$d)
   # solve the system 
   D <- proj.sum.one(t(KS$v %*% KS$d %*% t(KS$u) %*% t(P)))
   
 	Rb <- NULL
 	if (!P.only) {
       # compute Rb
       LS <- svd(D)
       # remove the small singular values
       d.min <- max(LS$d) * precision
       LS$d[LS$d < d.min] <- 0
       non.zero.ind <- LS$d != 0
       LS$d[non.zero.ind] <- 1 / LS$d[non.zero.ind]
       LS$d <- diag(LS$d)
#       # solve the system
       Rb <- LS$v %*% LS$d %*% t(LS$u) %*% R
       }
 
   if (P.only) list(D = D, K = K)
   else list(D = D, K = K, r = Rb[,1])
   }


nmf.rbfn <- function(P, m, max.iter = 10, P.only = TRUE, 
							precision = 1e-6, tau = 0.1, ...) {

	
	rbfn <- make.null.rbfn()
	S <- apply(P, 2, max)
	o <- order(S, decreasing = TRUE)
	rbfn$c <- P[o[1:m],]
	
	D <- rbfn.norm.design.matrix(rbfn,P)
	P <- rbfn$c
	
	if (P.only) list(D = D, K = P)
   else list(D = D, K = P[,2:ncol(P)], r = P[,1])
   }


print("nmf.R loaded")