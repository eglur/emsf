## Sem bias ##

source("data.manipulation.R")
source("rbfn.gaussian.R")
source("matrix.outer.product.fast.R")
#source("matrix.outer.product.bellman.R") 
source("rbfn.rl.R")
source("util.R")


rbfn.features <- function(X, a, rbfn, rbfn.dm= rbfn.gaussian.design.matrix) {
# rbfn.dm defines the way the rbfn will compute the design matrix;
	H <- rbfn.dm(rbfn, X) #design matrix
    
	B <- matrix(0, nrow(H), ncol(rbfn$w) * ncol(H))
	for (i in 1:nrow(H)) {
        B[i, ((a[i] -1) * ncol(H) + 1) :  (a[i] * ncol(H))] <- H[i,]
        }
	B
   }


lstd.q.rbfn <- function(ft.s, ft.s2, a, a2, r, rbfn, df = 0.9, precision = 1e-6,
                  verbose = FALSE) {
	# "ft.s" is the matrix whose columns are the features of 's'
	# "r" is a vector with the rewards
	# "df" is the discount factor
	# "rbfn" is a rbf network with |A| output units

	if (verbose) print("Computing A and b...")
	ft.dif <- ft.s - df * ft.s2

	k <- (nrow(rbfn$c)) * ncol(rbfn$w) 
	# if you want a bias:
	# k <- (nrow(rbfn$c) + 1) * ncol(rbfn$w) 
	
	b <- matrix(0, k, 1)
	
	# Least-square fixed-point approximation (Lagoudakis & Parr, 2003)
	A <- matrix.outer.product.fast(ft.s, ft.dif, nrow(rbfn$c), a, a2) 
	# if you want a bias:
	#A <- matrix.outer.product.fast(ft.s, ft.dif, nrow(rbfn$c) + 1, a, a2) 
     for (i in 1:length(r)) b <- b + ft.s[i,] * r[i]
	
	# Bellman residual minimizing approximation for deterministic environments 
     # (see Lagoudakis & Parr, 2003)
	#A <- matrix.outer.product.bellman(ft.dif, ft.dif, nrow(rbfn$c), a, a2) 
     #for (i in 1:length(r)) b <- b + ft.dif[i,] * r[i]
  	
	if (verbose) print("Solving the system...")
	# singular value decomposition
	A <- svd(A)
	# remove the small singular values
	d.min <- max(A$d) * precision
	A$d[A$d < d.min] <- 0
	non.zero.ind <- A$d != 0
	A$d[non.zero.ind] <- 1 / A$d[non.zero.ind]
	A$d <- diag(A$d)
	# solve the system	
	rbfn$w <- matrix(A$v %*% A$d %*% t(A$u) %*% b, nrow(rbfn$w), ncol(rbfn$w))
	
	rbfn
	}


lspi.rbfn <- function(samples, rbfn, df = 0.9, num.iterations = 10, precision = 
               1e-6, verbose = FALSE, rbfn.dm = rbfn.gaussian.design.matrix) {
	# "samples" is a matrix containing the data in the format 
	# l = list(s,a,r,s2)
	
	it <- 0	
	if (verbose) print("Computing features of S...")
	ft.s <- rbfn.features(samples$s, samples$a, rbfn, rbfn.dm = rbfn.dm)
     
  	while (it < num.iterations) {
 		# these are the only features that MUST be recomputed at  each
          # iteration
		if (verbose) print("Computing features of S2...")
		o2 <- rbfn.dm(rbfn, samples$s2) %*% rbfn$w
		a2 <- apply(o2, 1, which.max) 
		ft.s2 <- rbfn.features(samples$s2, a2, rbfn, rbfn.dm= rbfn.dm)

		rbfn <- lstd.q.rbfn(ft.s, ft.s2, samples$a, a2, samples$r, rbfn, df, 
                     precision, verbose)
		it <- it + 1
		}
	rbfn
	}



print("lspi.rbfn.R loaded")	