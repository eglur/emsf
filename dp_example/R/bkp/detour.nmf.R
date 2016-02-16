# O único KL que tá certo é o 'nmf.mult.kl.no.proj'
# TODOS OS OUTROS ESTÃO ERRADOS!!!

# Non-negative matrix factorization
#-------------------------------------------------------------------------------
# All these functions perform a non-negative matrix factorization (NMF), that
# is, given a (n x m) matrix A and an integer k, they return two matrices, 
# K (n x k) and L (k x m) so that ||A - LK|| is minimized in some sense (the
# norm ||.|| being  minimized varies over the functions).


source("util.R") # for 'euclidean.norm()'
source("proj.sum.one.R")

nmf.mult <- function(W, k, iter.max = 50, K = NULL, L = NULL,  R =
                     NULL, epsilon = -Inf, precision = 1e-10, hist = FALSE, plot
                     = FALSE) {
# Lee and Seung's multiplicative rule (projected)
# This function tries to minimize the Euclidean distance
# In order to avoid a division by zero, W should not have a zero column 
# 'hist' indicates whether the history should be returned
# IMPORTANT: W is cbind(R, P) - the first column of W is not projected!
	
	if (is.null(K)) {
		# Every element of K and L should be > 0
		K <- matrix(runif(k * ncol(W), 0.1, 1), k, ncol(W))
		K <- cbind(K[,1], proj.sum.one(K[,2:ncol(K)]))
		#Another possibility:
		#K <- W[sample(1:nrow(W), k), ]
		}
	else K <- cbind(R,K)
	
	if (is.null(L)) { #why not find the least-squares solution here?
		L <- matrix(runif(nrow(W) * k, 0.1, 1), nrow(W), k)
		L <- L / apply(L, 1, sum)
		#KS <- svd(t(K))
		# remove the small singular values
		#d.min <- max(KS$d) * precision
		#KS$d[KS$d < d.min] <- 0
		#non.zero.ind <- KS$d != 0
		#KS$d[non.zero.ind] <- 1 / KS$d[non.zero.ind]
		#KS$d <- diag(KS$d)
		# solve the system	
		#L <- proj.sum.one(t(KS$v %*% KS$d %*% t(KS$u) %*% t(W)))
		}

     h <- NULL

	dif <- Inf
	it <- 0
	while (dif > epsilon && it < iter.max) {
		K.old <- K		
		
		# update L
		N <- W %*% t(K)
		D <- L %*% K %*% t(K)
         	          
          
          # I added this to avoid division-by-zero:
          min.d <- max(N) * precision
          ind <- D == 0
          N[ind] <- 1 
          D[ind] <- 1 
		
          L <- L * N / D
		L <- proj.sum.one(L)
		
		# update K
		N <- t(L) %*% W
		D <- t(L) %*% L %*% K
          
          # I added this to avoid division-by-zero:
          min.d <- max(N) * precision
          ind <- D == 0
          N[ind] <- 1 
          D[ind] <- 1 
		
          K <- K * N / D
		K <- cbind(K[,1], proj.sum.one(K[,2:ncol(K)]))
		
		# compute changes
		dif <- max(abs(K - K.old))
		it <- it + 1
          
          if (hist) {
            h <- c(h, sse(L%*% K, W))
		  if (plot) plot(h, t="l")
            }        
		}
     list(D = L, K = K[,2:ncol(K)], r = K[,1], hist = h)
	}




nmf.kmeans <- function(W, k, iter.max = 10, L = NULL, K = NULL,  
                  R = NULL,...) {
# K-means factorization
# L is here for compatibility only
   if (is.null(K)) seed <- k # only gives the desired number of centers 
   else seed <- cbind(R,K) # gives the initial centers
   D <- kmeans(W, seed, iter.max = iter.max, ...)
   L <- matrix(0, nrow(W), k)
   for (i in 1:nrow(L)) L[i,D$cluster[i]] <- 1
   list(L= L, K = D$centers[,2:ncol(D$centers)], R = D$centers[,1])
   }
   






nmf.mult.Ralone <- function(R, P, k, iter.max = 50, K = NULL, L = NULL,
                              epsilon = 1e-6, precision = 1e-6) {
# Lee and Seung's multiplicative rule (projected)
# This function tries to minimize the Euclidean distance
# In order to avoid a division by zero, W should not have a zero column 
# IMPORTANT: W is cbind(R, P) - the first column of W is not projected!
     
     if (is.null(K)) {
          # Every element of K and L should be > 0
          K <- matrix(runif(k * ncol(P), 0.1, 1), k, ncol(P))
          K <- K / apply(K, 1, sum)
          #Another possibility:
          #K <- P[sample(1:nrow(P), k), ]
          }
     
     if (is.null(L)) { #why not find the least-squares solution here?
          L <- matrix(runif(nrow(P) * k, 0.1, 1), nrow(P), k)
          L <- L / apply(L, 1, sum)
          KS <- svd(t(K))
          # remove the small singular values
          d.min <- max(KS$d) * precision
          KS$d[KS$d < d.min] <- 0
          non.zero.ind <- KS$d != 0
          KS$d[non.zero.ind] <- 1 / KS$d[non.zero.ind]
          KS$d <- diag(KS$d)
          # solve the system  
          L <- proj.sum.one(t(KS$v %*% KS$d %*% t(KS$u) %*% t(P)))
          }

     #h <- NULL

     dif <- Inf
     it <- 0
     while (dif > epsilon && it < iter.max) {
          K.old <- K          
          
          # update L
          N <- P %*% t(K)
          D <- L %*% K %*% t(K)
          D[D==0] <- 1 # I added this to avoid division-by-zero
          L <- L * N / D
          #L <- proj.sum.one(L)
          
          # update K
          N <- t(L) %*% P
          D <- t(L) %*% L %*% K
          D[D==0] <- 1 # I added this to avoid division-by-zero
          K <- K * N / D
          #K <- proj.sum.one(K)
          
          # compute changes
          dif <- max(abs(K - K.old))
          it <- it + 1
          #h <- c(h, sse(L%*% K, P))
          #plot(h)
          #plot(h, t="l")
          }
      
      LS <- svd(L)
      # remove the small singular values
      d.min <- max(LS$d) * precision
      LS$d[LS$d < d.min] <- 0
      non.zero.ind <- LS$d != 0
      LS$d[non.zero.ind] <- 1 / LS$d[non.zero.ind]
      LS$d <- diag(LS$d)
      # solve the system
      R2 <- LS$v %*% LS$d %*% t(LS$u) %*% R
     
     print(sum(apply(L,1,sum)))
     print(sum(apply(K,1,sum)))
          
     list(L = L, K = K, R = R2)
     }



nmf.mult.kl <- function(W, k, iter.max = 50, K = NULL, L = NULL,  R =
                     NULL, epsilon = 1e-6, precision = 1e-6) {
# Lee and Seung's multiplicative rule (projected)
# This function tries to minimize the Kullback-Leibler divergence
# IMPORTANT: W is cbind(R, P) - the first column of W is not projected!
# In order to avoid a division by zero, W should not have a zero column 
	
     if (is.null(K)) {
          # Every element of K and L should be > 0
          K <- matrix(runif(k * ncol(W), 0.1, 1), k, ncol(W))
          K <- K / apply(K, 1, sum)
          #Another possibility:
          #K <- W[sample(1:nrow(W), k), ]
          }
     else K <- cbind(R,K)
     
     if (is.null(L)) { #why not find the least-squares solution here?
          L <- matrix(runif(nrow(W) * k, 0.1, 1), nrow(W), k)
          L <- L / apply(L, 1, sum)
          #KS <- svd(t(K))
          # remove the small singular values
          #d.min <- max(KS$d) * precision
          #KS$d[KS$d < d.min] <- 0
          #non.zero.ind <- KS$d != 0
          #KS$d[non.zero.ind] <- 1 / KS$d[non.zero.ind]
          #KS$d <- diag(KS$d)
          # solve the system  
          #L <- proj.sum.one(t(KS$v %*% KS$d %*% t(KS$u) %*% t(W)))
          }

 	
 	dif <- Inf
	it <- 0
	while (dif > epsilon && it < iter.max) {
		K.old <- K		
		
          D <- L %*% K        
          D[D==0] <- 1 # I added his to avoid division-by-zero
          
          # update K
		A <- W / D
		N <- t(L) %*% A
		D <- apply(L, 2, sum)
          D <- matrix(rep(D, ncol(K)), nrow(K), ncol(K), byrow = FALSE)
          D[D==0] <- 1 # I added his to avoid division-by-zero
		K <- K * N / D
		K <- cbind(K[,1], proj.sum.one(K[,2:ncol(K)]))
		
		# update L
          N <- t(K %*% t(A))
          #N <- A %*% t(K)
		D <- apply(K, 1, sum)
          D <- matrix(rep(D, nrow(L)), nrow(L), ncol(L), byrow = TRUE)
	     D[D==0] <- 1 # I added his to avoid division-by-zero
       	L <- L * N / D 
		L <- proj.sum.one(L)
          
		# compute changes
		dif <- max(abs(K - K.old))
		it <- it + 1
		}
	
		
     list(L = L, K = K[,2:ncol(K)], R = K[,1])
	}
	


nmf.mult.kl.no.proj <- function(R, P, k, iter.max = 50, K = NULL, L = NULL,  
                  epsilon = 1e-6, precision = 1e-6) {
# Lee and Seung's multiplicative rule (NOT projected)
# This function tries to minimize the Kullback-Leibler divergence
# In order to avoid a division by zero, W should not have a zero column 
     
       if (is.null(K)) {
          # Every element of K and L should be > 0
          K <- matrix(runif(k * ncol(P), 0.1, 1), k, ncol(P))
          K <- K / apply(K, 1, sum)
          #K <- P[sample(1:nrow(P), k), ]
          }
          
     if (is.null(L)) {
          L <- matrix(runif(nrow(P) * k, 0.1, 1), nrow(P), k)
          L <- L / apply(L, 1, sum)
          }
   
     dif <- Inf
     it <- 0
     while (dif > epsilon && it < iter.max) {
          K.old <- K          
          
          # update K
          D <- L %*% K 
          D[D==0] <- 1
          A <- P / D
          N <- t(L) %*% A
          D <- apply(L, 2, sum)
          D <- matrix(rep(D, ncol(K)), nrow(K), ncol(K), byrow = FALSE)
          K <- K * N / D
         
          # update L
          #N <- t(K %*% t(A))
          N <- A %*% t(K)
          D <- apply(K, 1, sum)
          D <- matrix(rep(D, nrow(L)), nrow(L), ncol(L), byrow = TRUE)
          L <- L * N / D
          
          # compute changes
          dif <- max(abs(K - K.old))
          it <- it + 1
          
#           NZL <- apply(L, 2, sum) != 0
#           NZK <- apply(K, 1, sum) != 0
#         
#           L <- L[,NZL]
#           K <- K[NZL,]
#           L <- L[,NZK]
#           K <- K[NZK,]
      }
     
     L <- proj.sum.one(L)
     K <- proj.sum.one(K)
     
     LS <- svd(L)
     # remove the small singular values
     d.min <- max(LS$d) * precision
     LS$d[LS$d < d.min] <- 0
     non.zero.ind <- LS$d != 0
     LS$d[non.zero.ind] <- 1 / LS$d[non.zero.ind]
     LS$d <- diag(LS$d)
     # solve the system
     R2 <- LS$v %*% LS$d %*% t(LS$u) %*% R
      
     list(L = L, K = K, R = R2)
     }



nmf.mult.kl.Ralone <- function(R, P, k, iter.max = 50, K = NULL, L = NULL,   
                              epsilon = 1e-6, precision = 1e-6) {
# Lee and Seung's multiplicative rule (projected)
# This function tries to minimize the Kullback-Leibler divergence
# IMPORTANT: W is cbind(R, P) - the first column of W is not projected!
# In order to avoid a division by zero, W should not have a zero column 
     
       if (is.null(K)) {
          # Every element of K and L should be > 0
          K <- matrix(runif(k * ncol(P), 0.1, 1), k, ncol(P))
          K <- K / apply(K, 1, sum)
          #K <- P[sample(1:nrow(P), k), ]
          }
          
     if (is.null(L)) {
          L <- matrix(runif(nrow(P) * k, 0.1, 1), nrow(P), k)
          L <- L / apply(L, 1, sum)
          }
   
     dif <- Inf
     it <- 0
     while (dif > epsilon && it < iter.max) {
          K.old <- K          
          
          # update K
          D <- L %*% K
          D[D==0] <- 1 # I added his to avoid division-by-zero
          A <- P / D
          N <- t(L) %*% A
          D <- apply(L, 2, sum)
          D[D==0] <- 1 # I added his to avoid division-by-zero
             
          K <- K * N / D
          K <- proj.sum.one(K)
          
          # update L
          N <- t(K %*% t(A))
          D <- apply(K, 1, sum)
          D[D==0] <- 1 # I added his to avoid division-by-zero
          L <- L * N / D 
          L <- proj.sum.one(L)
          
          # compute changes
          dif <- max(abs(K - K.old))
          it <- it + 1
          
          }
     
     LS <- svd(L)
     # remove the small singular values
     d.min <- max(LS$d) * precision
     LS$d[LS$d < d.min] <- 0
     non.zero.ind <- LS$d != 0
     LS$d[non.zero.ind] <- 1 / LS$d[non.zero.ind]
     LS$d <- diag(LS$d)
     # solve the system
     R2 <- LS$v %*% LS$d %*% t(LS$u) %*% R
      
     list(L = L, K = K, R = R2)
     }



nmf.pg <- function(W, k, iter.max = 50, K = NULL, L = NULL, R = NULL, 
            epsilon = 1e-6, precision = 0.001, iter.max.inner = 1000,
            epsilon.dec = 0.1) {
# NMF by alternative non-negative least squares using projected gradients
# Author: Chih-Jen Lin, National Taiwan University
# Converted from Matlab to R by Andre' M. S. Barreto, COPPE / UFRJ

# L,K: initial solution
# epsilon: tolerance for a relative stopping condition
# iter.max: limit of iterations
   
   if (is.null(K)) {
         # Every element of K and L should be > 0
         #K <- matrix(runif(k * ncol(W), 0.1, 1), k, ncol(W))
         #K <- K / apply(K, 1, sum)
         K <- W[sample(1:nrow(W), k), ]
         }
         
   if (is.null(L)) { #why not find the least-squares solution here?
         #L <- matrix(runif(nrow(W) * k, 0.1, 1), nrow(W), k)
         #L <- L / apply(L, 1, sum)
         KS <- svd(t(K))
         # remove the small singular values
         d.min <- max(KS$d) * precision
         KS$d[KS$d < d.min] <- 0
         non.zero.ind <- KS$d != 0
         KS$d[non.zero.ind] <- 1 / KS$d[non.zero.ind]
         KS$d <- diag(KS$d)
         # solve the system	
         L <- proj.sum.one(t(KS$v %*% KS$d %*% t(KS$u) %*% t(W)))
         }
   
   
   grad.L <- L %*% (K %*% t(K)) - W %*% t(K)
   grad.K <- (t(L) %*% L) %*% K - t(L) %*% W
      
   grad <- c(grad.L, grad.K)
   norm.initial.grad <- euclidean.norm(grad)
   epsilon.L <- max(precision,epsilon) * norm.initial.grad
   epsilon.K <- epsilon.L
   
   #norm.proj.grad <- Inf
   dif <- Inf
   it <- 0
   
   while (dif > epsilon && it < iter.max){
      # proj.grad <- c(grad.L[grad.L < 0 | L > 0], grad.K[grad.K < 0 | K >0])
      # norm.proj.grad <- euclidean.norm(proj.grad)
      
      K.old <- K
      
      R <- nnls(W, L, K, epsilon.K, iter.max.inner)
      K <- proj.sum.one(R$mat)
      grad.K <- R$grad
      if (R$it == 1) epsilon.K <- epsilon.dec * epsilon.K
      
      R <- nnls(t(W),t(K),t(L),epsilon.L,iter.max.inner) 
      L <- proj.sum.one(t(R$mat))
      grad.L <- t(R$grad)
      if (R$it == 1) epsilon.L <- epsilon.dec * epsilon.L
   
      dif <- max(abs(K - K.old))
	 it <- it + 1
      }
      
  list(L = L, K = K)
  }
   
   

nnls <- function(W, L, K, epsilon, iter.max, alpha = 1,
               beta = 0.1, max.iter.search = 20) {

# Auxiliar function of the previous nmf function
# K, grad: output solution and gradient
# iter: #iterations used
# W, L: constant matrices
# Hinit: initial solution
# epsilon: stopping epsilonerance
# iter.max: limit of iterations

   LtW <- t(L) %*% W
   LtL <- t(L) %*% L
   
   it <- 0
   norm.proj.grad <- Inf
   grad <- NULL
   
   while (norm.proj.grad > epsilon && it < iter.max) {
      grad <- LtL %*% K - LtW
      #proj.grad <- grad[grad < 0 | K > 0]
      proj.grad <- grad
      dim(proj.grad) <- length(proj.grad)
      norm.proj.grad <- euclidean.norm(proj.grad)
      
      # search step size 
      Kp <- NULL
      for (inner_iter in 1:max.iter.search) {
         Kn <- K - alpha * grad
         #Kn[Kn < 0] <- 0
         
         d <- Kn - K
         gradd <- sum(grad * d)
         
         dQd <- sum((LtL %*% d) * d)

         suff_decr <- (0.99 * gradd + 0.5 * dQd) < 0
         
         if (inner_iter==1) {
            decr_alpha <- !suff_decr
            Kp <- K
            }
         
         if (decr_alpha) {
            if (suff_decr) {
               K <- Kn
               break
               }
            else alpha <- alpha * beta
            }
         else
            if (!suff_decr || identical(Kp,Kn)) {
               K <- Kp
               break
               }
            else {
               alpha <- alpha / beta
               Kp <- Kn
               }
         }
     it <- it + 1
     }
   list(mat = K, grad = grad, it = it)  
   }
   
   
   

nmf.ls <- function(W, k, iter.max = 50, L = NULL, K = NULL,  
                  R = NULL, epsilon = 1e-6, precision = 1e-6) {
# Non-negative projected least-squares
# No convergence guarantees
   if (is.null(K)) K <- W[sample(1:nrow(W), k), ]
   if (is.null(L)) L <- matrix(0, nrow(W), k)
   
   # h <- NULL
   dif <- Inf
   it <- 0
   while (it <= iter.max && dif > epsilon) {
      ##  FIX K AND COMPUTE L ##
      # singular value decomposition
      KS <- svd(t(K))
      # remove the small singular values
      d.min <- max(KS$d) * precision
      KS$d[KS$d < d.min] <- 0
      non.zero.ind <- KS$d != 0
      KS$d[non.zero.ind] <- 1 / KS$d[non.zero.ind]
      KS$d <- diag(KS$d)
      
      # solve the system	
      L <- proj.sum.one(t(KS$v %*% KS$d %*% t(KS$u) %*% t(W)))
   
      
      K.old <- K
      
      ##  FIX L AND COMPUTE K ##
      # singular value decomposition
      LS <- svd(L)
      # remove the small singular values
      d.min <- max(LS$d) * precision
      LS$d[LS$d < d.min] <- 0
      non.zero.ind <- LS$d != 0
      LS$d[non.zero.ind] <- 1 / LS$d[non.zero.ind]
      LS$d <- diag(LS$d)
      # solve the system
      K <- LS$v %*% LS$d %*% t(LS$u) %*% W 
      K <- cbind(K[,1], proj.sum.one(K[,2:ncol(K)])) # considering P includes R
      
      # h <- c(h, sse(L%*% K, W))
      # plot(h)
      dif <- max(abs(K - K.old))
      it <- it +1 
      }
   #print(h[length(h)])
   list(L = L, K = K[,2:ncol(K)], R = K[,1])
   }
   
   
nmf.ls.Ralone <- function(R, P, k, iter.max = 50, L = NULL, K = NULL,  
                  epsilon = 1e-6, precision = 1e-6) {
# Non-negative projected least-squares
# No convergence guarantees
   if (is.null(K)) K <- P[sample(1:nrow(P), k), ]
   if (is.null(L)) L <- matrix(0, nrow(P), k)
   
   # h <- NULL
   dif <- Inf
   it <- 0
   LS <- NULL
   while (it <= iter.max && dif > epsilon) {
      ##  FIX K AND COMPUTE L ##
      # singular value decomposition
      KS <- svd(t(K))
      # remove the small singular values
      d.min <- max(KS$d) * precision
      KS$d[KS$d < d.min] <- 0
      non.zero.ind <- KS$d != 0
      KS$d[non.zero.ind] <- 1 / KS$d[non.zero.ind]
      KS$d <- diag(KS$d)
      
      # solve the system 
      L <- proj.sum.one(t(KS$v %*% KS$d %*% t(KS$u) %*% t(P)))
   
      
      K.old <- K
      
      ##  FIX L AND COMPUTE K ##
      # singular value decomposition
      LS <- svd(L)
      # remove the small singular values
      d.min <- max(LS$d) * precision
      LS$d[LS$d < d.min] <- 0
      non.zero.ind <- LS$d != 0
      LS$d[non.zero.ind] <- 1 / LS$d[non.zero.ind]
      LS$d <- diag(LS$d)
      # solve the system
      K <- proj.sum.one(LS$v %*% LS$d %*% t(LS$u) %*% P)
      
      # h <- c(h, sse(L%*% K, P))
      # plot(h)
      dif <- max(abs(K - K.old))
      it <- it +1 
      }
   
   # solve the system using the previously computed 'LS'
   R2 <- LS$v %*% LS$d %*% t(LS$u) %*% R
   
   list(L = L, K = K, R = R2)
   }
   
   

# EXPERIMENTS

compare.ls.mult <- function(num.avg = 20, size = 100, sd = 5, perc = 0.3,
                     iter.max = 50) {
   res <- matrix(0, num.avg, 4)
   for (i in 1:num.avg) {
      P <- normal.transition.matrix(size,sd = sd)
      k <- round(perc * size)
      
      t <- system.time(R <- nmf.mult(P, k, iter.max = iter.max))
      res[i,1] <- t[1]
      R$K <- cbind(R$R, R$K)
      res[i,3] <- sse(R$L %*% R$K, P)
      
      t <- system.time(R <- nmf.ls(P, k, iter.max = iter.max))
      res[i,2] <- t[1]
      R$K <- cbind(R$R, R$K)
      res[i,4] <- sse(R$L %*% R$K, P)
      
      print(res[i,])
      }
   as.data.frame(res)
   }
      

compare.functions.dp <- function(functions, num.avg = 20, size = 100, sd =   
                           5, perc = 0.3, iter.max = 50, df = 0.9) {
   res <- matrix(0, num.avg, 2)
   for (i in 1:num.avg) {
      P <- normal.transition.matrix(size,sd = sd)
      R <- runif(size)
      k <- round(perc * size)
      
      V <- value.iteration(R, P, df = df)$Q
      
      
      W <- cbind(R,P)
      S <- functions[[1]](W, k, iter.max = iter.max)
      Pb <- S$K %*% S$L
      Vb <- value.iteration(S$R, Pb, df = df)$Q
      
      S2 <- functions[[2]](W, k, iter.max = iter.max)
      Pb2 <- S2$K %*% S2$L
      Vb2 <- value.iteration(S2$R, Pb2, df = df)$Q
      
      V <- V - mean(V)
      Vb <- S$L %*% Vb
      Vb <- Vb - mean(Vb)
      Vb2 <- S2$L %*% Vb2
      Vb2 <- Vb2 - mean(Vb2)
      
      res[i,1] <- sse(V, Vb)
      res[i,2] <- sse(V, Vb2)
      
      }
  res
  }
      


compare.Ralone <- function(functions, num.avg = 20, size = 100, sd = 5, 
                    perc = 0.3, iter.max = 50, df = 0.9) {
   res <- matrix(0, num.avg, 2)
   for (i in 1:num.avg) {
      P <- normal.transition.matrix(size,sd = sd)
      R <- runif(size)
      k <- round(perc * size)
      
      V <- value.iteration(R, P, df = df)$Q
      
      
      W <- cbind(R,P)
      S <- functions[[1]](W, k, iter.max = iter.max)
      Pb <- S$K %*% S$L
      Vb <- value.iteration(S$R, Pb, df = df)$Q
      
      S2 <- functions[[2]](R, P, k, iter.max = iter.max)
      Pb2 <- S2$K %*% S2$L
      Vb2 <- value.iteration(S2$R, Pb2, df = df)$Q
      
      V <- V - mean(V)
      Vb <- S$L %*% Vb
      Vb <- Vb - mean(Vb)
      Vb2 <- S2$L %*% Vb2
      Vb2 <- Vb2 - mean(Vb2)
      
      res[i,1] <- sse(V, Vb)
      res[i,2] <- sse(V, Vb2)
      
      }
  res
  }
      
  
      
      
 	
print("detour.nmf.R loaded")