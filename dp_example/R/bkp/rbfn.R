dyn.load(paste("./c/compute_rbf_widths", .Platform$dynlib.ext, sep = ""))
#melhorar

##COLOCAR BIAS##
# OBS.: the function "rbfn.spread.centers" has had its name changed to
# "uniform.grid", which is defined in the file "util.R"

source("util.R")
source("dist.mat.R")

int.prod <- function(v1, v2) {
	drop(drop(v1) %*% drop(v2))
	}

euclidean.norm <- function(v) {
	sqrt(int.prod(v,v))
	}

sqr.euclidean.norm <- function(v) {
	int.prod(v,v)
	}

gaussian <- function(dist,  width = 1) {
	exp(-dist / width)
	}

pseudo.gaussian <- function(dist, width = 1) {
	r <- array(0, length(dist))
	ind <- dist < width
	r[ind] <- (1 - dist[ind] / width)**2
	r
	}



make.rbfn <- function(X, Y, num.rbfs, random.centers = FALSE, 
		      num.neighbors=max(1,num.rbfs %/% 5), 
		      scale = 2, norm= sqr.euclidean.norm, 
	              act.function = gaussian) {
# Here, the criteria to set the RBFs' widths are 'scale' and 'num.neighbors'.
# The width of the i-th RBF will be scale *  avg.dist, where 'avg_dist' is the
# average distance to the 'num.neighbors' nearest neighbours

# Y is either the output data or the dimension of the output space
# num.neighbors is the number of neighbors to be considered in the computation of the widths
# scale determines the proportion between the widths and the average distance to the neighbors

	# set the centers
	max.X <- apply(X, 2, max)
	min.X <- apply(X, 2, min)

	c <- matrix(0, num.rbfs, ncol(X))
	if (random.centers) for (j in 1:ncol(X)) c[,j] <- runif(num.rbfs, min.X[j], max.X[j])
	else c <- uniform.grid(c, max.X, min.X)

	# set the widths
	dists <- matrix(0, num.rbfs, num.rbfs)
	s <- array(0, num.rbfs)

	if (num.rbfs > 1) {
		for (i in 1:(num.rbfs-1)) {
			for (j in (i+1):num.rbfs) {
				dists[i,j] <- norm(c[i,] - c[j,])
				dists[j,i] <- dists[i,j]
				}
			}
		for (i in 1:num.rbfs) s[i] <- scale *
mean(sort(dists[i,])[2:(num.neighbors+1)])
		}
	else {
		dists <- array(0, nrow(X))
		for (i in 1:(num.rbfs-1)) dists[i] <- norm(c[1,] - X[i,])
		s[i] <- scale * max(dists)
		}
		


	# set the output weights
	if (is.null(dim(Y))) { # if Y is a scalar
            w <- matrix(rnorm(num.rbfs * Y, 0, 1), num.rbfs, Y) 
            }
     else {	 # is Y is a matrix with the data
		max.Y <- apply(Y, 2, max)
		min.Y <- apply(Y, 2, min)
		w <- matrix(0, num.rbfs, ncol(Y))
		for (i in 1:ncol(Y)) w[,i] <- runif(num.rbfs, min.Y[i], max.Y[i])
		}

	list(c = c, s = s, w = w, act.function = act.function, norm = norm)
	}

	

make.rbfn.tau <- function(X, dim.output, num.rbfs, random.centers = FALSE, 
                     tau = 0.5, norm = sqr.euclidean.norm, act.function =
                     gaussian) {
# Here, the criterion to set the RBFs' widths is tau, that is, the width
# of the i-th RBF will be 'tau' at its nearest neighbor

# dim.output is the dimension of the output weight
# tau is the overlapping activation

	# set the centers
	max.X <- apply(X, 2, max)
	min.X <- apply(X, 2, min)

	c <- matrix(0, num.rbfs, ncol(X))
	if (random.centers) {
        for (j in 1:ncol(X)) c[,j] <- runif(num.rbfs,min.X[j],max.X[j])
        }
	else c <- uniform.grid(c, max.X, min.X)

	# set the widths
	dists <- matrix(0, num.rbfs, num.rbfs)
	for (i in 1:(num.rbfs-1)) {
		for (j in (i+1):num.rbfs) {
			dists[i,j] <- norm(c[i,] - c[j,])
			dists[j,i] <- dists[i,j]
			}
		}
	s <- array(0, num.rbfs)
	for (i in 1:num.rbfs) {
		dists[i,i] <- Inf
		s[i] <- -min(dists[i,]) / log(tau)
		}

	# set the output weights
	w <- matrix(rnorm(num.rbfs * dim.output, 0, 1), num.rbfs, dim.output) 
	
	list(c = c, s = s, w = w, act.function = act.function, norm = norm)
	}
	
make.null.rbfn <- function(act.function = gaussian, norm = sqr.euclidean.norm) {
	# returns the 'blueprint' of a rbf network
	list(c = NULL, s = NULL, w = NULL, act.function = act.function, norm =
norm)
	}



make.rbfn.bias <- function(X, dim.output, num.rbfs, random.centers = FALSE, 
								tau = 0.5, norm = sqr.euclidean.norm, 
								act.function = gaussian) {
# identical to "make.rbfn.rgd", but includes a bias term
# dim.output is the dimension of the output weight
# tau is the overlapping activation

	# set the centers
	max.X <- apply(X, 2, max)
	min.X <- apply(X, 2, min)

	c <- matrix(0, num.rbfs, ncol(X))
	if (random.centers) for (j in 1:ncol(X)) c[,j] <- runif(num.rbfs, min.X[j],
max.X[j])
	else c <- uniform.grid(c, max.X, min.X)

	# set the widths
	dists <- matrix(0, num.rbfs, num.rbfs)
	for (i in 1:(num.rbfs-1)) {
		for (j in (i+1):num.rbfs) {
			dists[i,j] <- norm(c[i,] - c[j,])
			dists[j,i] <- dists[i,j]
			}
		}
	s <- array(0, num.rbfs)
	for (i in 1:num.rbfs) {
		dists[i,i] <- Inf
		s[i] <- -min(dists[i,]) / log(tau)
		}

	# set the output weights
	w <- matrix(rnorm((num.rbfs + 1) * dim.output, 0, 1), num.rbfs + 1, dim.output) 
	
	list(c = c, s = s, w = w, act.function = act.function, norm = norm)
	}
	


rbfn.output <- function(rbfn, X, return.matrices = FALSE, normalize = FALSE) {
	V <- array(0, c(nrow(X), nrow(rbfn$c), ncol(X))) # keeps the vectors x - c
	
	for (i in 1:nrow(rbfn$c)) V[,i,] <- t(apply(X, 1, "-", rbfn$c[i,]))
	
	D <- matrix(0, nrow(X), nrow(rbfn$c)) # D keeps the norm of the difference between centers and inputs
	for (i in 1:nrow(rbfn$c)) D[,i] <- apply(V[,i, , drop = FALSE] , 1, rbfn$norm) 

	H <- matrix(0, nrow(X), nrow(rbfn$c)) #design matrix
	for (i in 1:nrow(rbfn$c)) H[,i] <- apply(D[,i, drop = FALSE], 2, rbfn$act.function, width = rbfn$s[i])  
   
	if (normalize) H <- H / apply(H, 1, sum)
	
	o <- H %*% rbfn$w

	if (!return.matrices)	result <- o
	else result <- list(o = o, V = V, D = D, H = H)
	result
	}


step.bp.gaussian <- function(rbfn, x, y, lr.centers = 1e-3, lr.widths = 1e-3, lr.weights = 1e-3,
										update.hidden = TRUE) {
	
	ro <- rbfn.output(rbfn, x, return.matrices = TRUE)

     e <- matrix(y - ro$o)
     # update linear weights
	rbfn$w <- rbfn$w + lr.weights * t(ro$H) %*% e
	
	if (update.hidden) {
		S <- matrix(rep(rbfn$s, nrow(ro$D)), nrow(ro$D), ncol(ro$D), byrow = TRUE)
		E <- e %*% t(rbfn$w)
		
		# update widths
		rbfn$s <- rbfn$s + lr.widths * apply(matrix(E * ro$H * ro$D / S**3, nrow(x), nrow(rbfn$c)), 2, sum)
	
		F <- E * ro$H / S**2
		F <- array(rep(F, dim(ro$V)[3]), dim(ro$V))
		F <- F * ro$V
		# update centers
		for (i in 1:nrow(rbfn$c)) rbfn$c[i,] <- rbfn$c[i,] +  lr.centers * apply(F[,i, , drop = FALSE], 2, sum)
		}
	rbfn
	}


test.bp.gaussian <- function(X, Y, rbfn, num.epochs = 1000, ...) {
	sse <- array(0, num.epochs)
	o <- NULL
	ra <- matrix(0, nrow(x), nrow(rbfn$c))
	for (i in 1:num.epochs) {
		rbfn <- step.bp.gaussian(rbfn, X, Y, ...)
		H <- rbfn.gaussian.design.matrix(rbfn, X)
		o <- H %*% rbfn$w
		sse[i] <- sum((Y-o)**2)
		#H <- cbind(cbind(Y,o),H)
		#lty <- c(c(1,1),rep(2,ncol(H)-2))
		#W <- matrix(rep(rbfn$w, length(x)),  length(x), nrow(rbfn$w), byrow = TRUE)
		#W <- cbind(o,H*W)
 		matplot(cbind(Y,o), t="l")
		}
	sse
	}

rbfn.combine <- function(rbfn1, rbfn2) {
	rbfn <- rbfn1
	rbfn$c <- rbind(rbfn1$c, rbfn2$c)
	rbfn$s <- c(rbfn1$s, rbfn2$s)
	rbfn$w <- rbind(rbfn1$w, rbfn2$w)
	rbfn
	}



 make.rbfn.simple <- function(centers, tau, num.actions = 1, 
                           num.neighbors = min(ncol(centers)+1,nrow(centers)-1),
                           same.width = TRUE, p.centers = 1,
                           use.c.code = FALSE) {
# "p.centers" defines the percentage of the centers to be used to
# compute the widths
   rbfn <- make.null.rbfn(act.function = gaussian, norm= sqr.euclidean.norm)
   # set the centers
   num.rbfs <- nrow(centers)
   rbfn$c <- centers
   rbfn$w <- matrix(rnorm(num.rbfs * num.actions), num.rbfs, num.actions)
   
   # set the widths
      d.mean <-  array(0, num.rbfs)
      rbfn$s <- array(0, num.rbfs)
      base.inds <- 1:num.rbfs
      nd <- round(p.centers * num.rbfs)
      blocks <- array(nd, num.rbfs %/% nd)
      extra  <- sample(1:length(blocks), num.rbfs %% nd, replace=TRUE)
      for (i in 1:length(extra)) blocks[extra[i]] <- blocks[extra[i]] + 1
      for (i in 1:length(blocks)) {
         i <- sample(1:length(base.inds), blocks[i])
         
         inds <- base.inds[i]
         base.inds <- base.inds[-i]
         
         e <- min(num.neighbors, length(inds)-1)
         if (!use.c.code) {
            dists <- as.matrix(dist(rbfn$c[inds,]))^2 # must be squared 
            for (i in 1:length(inds)) {
               dists[i,i] <- Inf
               d.mean[inds[i]] <- mean(sort(dists[i,])[1:e])
               }
           if (same.width) rbfn$s[inds] <- -mean(d.mean[inds])/log(tau)
           else rbfn$s[inds] <- -d.mean[inds]/log(tau) 
           }
         else rbfn$s[inds] <-  .C("compute_rbf_widths",
            as.double(t(rbfn$c[inds,])), as.integer(length(inds)),
            as.integer(ncol(rbfn$c)), as.integer(e),
            as.double(tau), array(0, length(inds)))[[6]]
         }
   rbfn
   }


make.rbfn.really.simple <- function(centers, tau, num.actions = 1) {
   rbfn <- make.null.rbfn(act.function = gaussian, norm= sqr.euclidean.norm)
   num.rbfs <- nrow(centers)
   rbfn$c <- centers
   rbfn$w <- matrix(rnorm(num.rbfs * num.actions), num.rbfs, num.actions)
   rbfn$s <- array(tau, num.rbfs)
   rbfn
   }


print("rbfn.R loaded")	
	
