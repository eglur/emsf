check.eigenvalues <- function(num.tests, n, m, sr, precision = 1e-10) {
# "sr" is the stochastic rank
# "m" must be in between n and the stochastic rank
	for (t in 1:num.tests) {
		O <- matrix(runif(sr*n), sr, n) # "kernel"
		O <- O / apply(O,1,sum)
		C <- matrix(runif(m*sr), m, sr)
		C <- C / apply(C,1,sum)
		K <- C %*% O

		D <- matrix(runif(m*n), n, m)
		D <- D / apply(D,1,sum)

		P  <- D %*% K
		Pb <- K %*% D

		v <- eigen(P)$values
  		v[abs(Re(v)) < precision] <- 0
 		v <- unique(v)
 		n <- length(v)
		r <- qr(P)$rank
		
		vb <- eigen(Pb)$values
  		vb[abs(Re(vb)) < precision] <- 0
 		vb <- unique(vb)
 		nb <- length(vb)
		rb <- qr(Pb)$rank

# 		print(paste("P:",n,"Pb:",nb, "Min", min(abs(v)), "Min_b",
# min(abs(vb)), "Rank:", r, "Rank_b", rb))
			
		if (rb < (nb-1)) {
			print(v)
			print(vb)
			S <- readline()
			}
		}
	}


find.last.row <- function(P) {
	A <- P
	A <- cbind(1,A)
	n <- ncol(A)
	v = array(0,n)
	v[1] <- 1
	v[n] <- 0.5
	s <- solve.QP(diag(n-2), array(0,n-2), A, bvec=v,meq=1,
factorized=FALSE)$solution
	rbind(P, s%*%P)
	}

# find.matrix <- function() {
# 	v <- seq(0,0.5,step=0.1)
# 	P <- 
# 	for (i in 1:length(v)) {
# 		for (i in 1:length(v)) {
# 
	
	
	

print("sf.mc.R loaded")