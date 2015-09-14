## Sem bias ##
source("lspi.rbfn.R")

lstd.q.rbfn <- function(ft.s, ft.dif, a, a2, r, rbfn, df = 0.9, precision =
								1e-6,verbose = FALSE) {
	# "ft.s" is the matrix whose columns are the features of 's'
	# "r" is a vector with the rewards
	# "df" is the discount factor
	# "rbfn" is a rbf network with |A| output units

	if (verbose) print("Computing A and b...")

	k <- (nrow(rbfn$c)) * ncol(rbfn$w) 
	# if you want a bias:
	# k <- (nrow(rbfn$c) + 1) * ncol(rbfn$w) 
	

	# Least-square fixed-point approximation (Lagoudakis & Parr, 2003)
  	n.col <- ncol(ft.s)
	A <- matrix(0, n.col, n.col)
	A <- .C("matrix_outer_product_fast", as.double(t(ft.s)),
				as.double(t(ft.dif)),  as.integer(n.col), 
				as.integer(nrow(ft.s)), as.integer(nrow(rbfn$c)), 
		  		as.integer(a), as.integer(a2), as.double(A))
	
	A <- matrix(A[[length(A)]], n.col, n.col, byrow = TRUE)
	      								
   
	# if you want a bias:
	#A <- matrix.outer.product.fast(ft.s, ft.dif, nrow(rbfn$c) + 1, a, a2) 
	b <- matrix(0, k, 1)
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
               1e-6, verbose = FALSE, rbfn.dm = rbfn.norm.design.matrix,
               tolerance = 1e-6, perc.imp = 1) {
	# "samples" is a matrix containing the data in the format 
	# l = list(s,a,r,s2)
	
	
	it <- 0	
	if (verbose) print("Computing features of S...")
	
   H <- rbfn.dm(rbfn, samples$s) #design matrix
    
	ft.s <- matrix(0, nrow(H), ncol(rbfn$w) * ncol(H))
	for (i in 1:nrow(H)) {
        ft.s[i,((samples$a[i] -1)*ncol(H)+1):(samples$a[i] * ncol(H))]<-H[i,]
		  }
  	
  	H <- rbfn.dm(rbfn, samples$s2)
  	
  	## Potentially remove (depending on whether perc.imp works)
    o2 <- H %*% rbfn$w
   a2 <- apply(o2, 1, which.max) 
   ft.s2 <- matrix(0, nrow(H), ncol(rbfn$w) * ncol(H))
   for (i in 1:nrow(H)) {
      ft.s2[i,((a2[i] -1)*ncol(H)+1):(a2[i] * ncol(H))]<-H[i,]
      }
   ft.dif <- ft.s - df * ft.s2   
  	## until here (dont cut things before nor after)
  	
  	old.w <- matrix(Inf, nrow(rbfn$w), ncol(rbfn$w))
   
   
   while (it < num.iterations && sum(abs(old.w-rbfn$w)) > tolerance) {
 		# these are the only features that MUST be recomputed at  each
          # iteration
      old.w <- rbfn$w
      
		if (verbose) print("Computing features of S2...")
		
		ind <- sample(1:nrow(H), round(perc.imp*nrow(H)))
		o2[ind,] <- H[ind,] %*% rbfn$w
		a2[ind] <- apply(o2[ind,], 1, which.max) 
      
  		ft.s2[ind,] <- 0
      for (i in ind) {
         ft.s2[i,((a2[i] -1)*ncol(H)+1):(a2[i] * ncol(H))] <- H[i,]
         }

      ft.dif[ind,] <- ft.s[ind,] - df * ft.s2[ind,]
      
		rbfn <- lstd.q.rbfn(ft.s, ft.dif, samples$a, a2, samples$r,
								rbfn, df, precision, verbose)
								
      ## REMOVE
#        v <- matrix(compute.v.rbfn(rbfn, S), 100,100)
#        pp(v)
# 
      it <- it + 1
		}
	rbfn
	}



print("lspi.rbfn.optimized.R loaded")	