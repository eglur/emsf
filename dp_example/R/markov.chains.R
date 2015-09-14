source("mdp.R")
source("util.R")
# source("rbfn.R")


kl <- function(P, Q, epsilon = 1e-10) {
# Kullback_leibler divergence
   D <- (P + epsilon)/(Q + epsilon) # to circumvent log(0)
   sum(P * log(D) - P + Q)
   }
   

group.inverse <- function(A) {
# A = I- P, where P is a transition-probability matrix
   n <- nrow(A)
   U <- A[1:(n-1), 1:(n-1)]
   d <- A[n, 1:(n-1)]

   iU <- solve(U)
   h <- d%*% iU

   beta <- drop(1 - sum(h))

   j <- array(1, n-1)
   delta <- drop(-h %*% iU %*% j)
   F <- iU - delta/beta * diag(n-1)

   A1 <- iU + (iU %*% j %*% h %*% iU) / delta - (F %*% j %*% h %*% F) / delta
   
   rbind(cbind(A1,-(F %*% j) / beta), cbind((h%*%F) / beta, delta / beta^2))
   }




sf.kl <- function(P, m, max.iter = 50, K = NULL, D = NULL,   
                  epsilon = -Inf, hist = FALSE, plot = FALSE,
		  norm = norm.1) {
# Lee and Seung's multiplicative rule 
# This function tries to minimize the Kullback-Leibler divergence
     
     if (is.null(K)) {
          # Every element of K and D should be > 0
          K <- P[sample(1:nrow(P), m), ]
 	  K <- 0.9 * K + 0.1 * (1/ncol(P))
          }
          
     if (is.null(D)) {
          D <- matrix(runif(nrow(P) * m, 0.1, 1), nrow(P), m)
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
          dif <- norm(K - K.old)
          it <- it + 1

	 if (hist) {
            h <- c(h, kl(P, D%*% K))
            if (plot) plot(h, t="l")
            }        
          }
     
    # See corollary 2 of Ho and Dooren (2007) to understand this
    # It will make both K and D stochastic
    DS <- apply(D,2,sum)
    KS <- apply(K,1,sum)
    D <- (D / matrix(rep(DS, nrow(D)), nrow(D), ncol(D), byrow = TRUE)) %*%
diag(DS * KS)
    K <- K / matrix(rep(KS, ncol(K)), nrow(K), ncol(K), byrow = FALSE)
     
     list(D = D, K = K)
     }



sf.kl.fast <- function(P, m, max.iter = 50, 
                  epsilon = -Inf, hist = FALSE, plot = FALSE,
		  norm = norm.1) {
# Lee and Seung's multiplicative rule 
# This function tries to minimize the Kullback-Leibler divergence
# This is a faster version, in which K is "heuristically" defined and 
#  D is iteratively refined
     
	# define K, which will be kept fixed
	max.rows     <- apply(P,2,max)
	max.rows.ind <- apply(P,2,which.max)
	o <- order(max.rows, decreasing = TRUE)
	sel.rows <- unique(max.rows.ind[o])[1:m]
	K <- P[sel.rows,]
	
     # define D (the rows selected in K don't need to be optimized)
	D <- matrix(runif(nrow(P) * m - m^2, 0.1, 1), nrow(P) - m, m)
     D <- D / apply(D, 1, sum)
	P <- P[-sel.rows,]
	
	dif <- Inf
     it <- 0
     h <- NULL
     while (dif > epsilon && it < max.iter) {
          D.old <- D
# 		K <- 0.9 * K + 0.1 * 1/ncol(K)
		# update D 
		# O(nrow(D) ncol(D) ncol(K)) + O(nrow(P) ncol(P) nrow(K))
          Dn <- D %*% K 
          Dn[Dn==0] <- 1 # I added this to avoid a division-by-zero
          A <- P / Dn
		N <- A %*% t(K)
          Dn <- apply(K, 1, sum)
          Dn <- matrix(rep(Dn, nrow(D)), nrow(D), ncol(D), byrow = TRUE)
          Dn[Dn==0] <- 1 # I added this to avoid a division-by-zero
          D <- D * N / Dn
		
	     # compute changes
          dif <- norm(D - D.old)
          it <- it + 1
		

		if (hist) {
            h <- c(h, kl(P, D%*% K))
            if (plot) plot(h, t="l")
            }        
          }

	D.new <- matrix(0, ncol(K), m)
     D.new[-sel.rows,] <- D
	for (i in 1:length(sel.rows)) {
		D.new[sel.rows[i],]  <- 0
		D.new[sel.rows[i],i] <- 1
		}
	D <- D.new

# See corollary 2 of Ho and Dooren (2007) to understand this
    # It will make both K and D stochastic
    DS <- apply(D,2,sum)
    KS <- apply(K,1,sum)
    D <- (D / matrix(rep(DS, nrow(D)), nrow(D), ncol(D), byrow = TRUE)) %*%
diag(DS * KS)
    K <- K / matrix(rep(KS, ncol(K)), nrow(K), ncol(K), byrow = FALSE)
     
     list(D = D, K = K)
     }


sf.exact <- function(P) {
## This is not exact!
	# define K, which will be kept fixed
	max.rows.ind <- apply(P,2,which.max)
	min.rows.ind <- apply(P,2,which.min)
	sel.rows <- unique(c(max.rows.ind, min.rows.ind))
	K <- P[sel.rows,]
	m <- length(sel.rows)
	D <- matrix(0, nrow(P), m)
# 	D[-sel.rows,] <- t(solve(K%*%t(K))%*%K%*%t(P[-sel.rows,]))
 	D[-sel.rows,] <- P[-sel.rows,] %*% t(K) %*% solve(K %*% t(K))
 	D[sel.rows,1:m] <- diag(m)
#  	D <- D / apply(D,1,sum) # to guarantee
	list(D=D, K=K)
	}


sf.exact.fast <- function(P) {
##wrong; left here for historical reasons

	max.rows.ind <- apply(P,2,which.max)
	max.rows.col <- array(0,ncol(P))
	max.rows.col[max.rows.ind] <- 1:length(max.rows.ind)


	dup.ind <- duplicated(max.rows.ind)
	
	D <- NULL
	K <- NULL

	if (length(dup.ind)==0)  { # Nothing to remove
		K <- P
		D <- diag(nrow(P))
		}
	else {
			
		dup <- unique(max.rows.ind[dup.ind])
		uni <- setdiff(max.rows.ind, dup)
		uni.cols <- max.rows.col[uni]

		A <- P

# 		A[,max.rows.col[dup[1]]] <- A[,max.rows.col[dup[1]]] + 
# 					apply(as.matrix(A[,uni.cols], ), 1, sum)
							
		A[uni,] <- 0
		A[,uni.cols] <- 0
		A[uni, uni.cols] <- diag(length(uni)) 
	
		min.rows.ind  <- apply(A,2,which.min)
		min <- setdiff(unique(min.rows.ind), uni)

		essential.rows <- unique(c(min,dup))
		print(essential.rows)

		K <- matrix(A[essential.rows,-uni.cols], length(essential.rows),
				  ncol(P)-length(uni.cols), byrow = TRUE)
		
		print(K)
		sel.rows <- sort(c(uni,essential.rows)) ##think
		print(sel.rows)

		D <- matrix(0, nrow(P), length(sel.rows))
		print(uni)
		print(uni.cols)

		D[,uni] <- P[,uni.cols]
		D[,-uni] <- P[,-uni.cols] %*% t(K) %*% solve(K%*%t(K))
		K <- A[sel.rows,]
		print(K)
		}

	#  	D <- D / apply(D,1,sum) # to guarantee
	list(D=D, K=K)
	}
	

sf.rbf <- function(P) {
	# define K, which will be kept fixed
	max.rows.ind <- apply(P,2,which.max)
	min.rows.ind <- apply(P,2,which.min)
	sel.rows <- unique(c(max.rows.ind, min.rows.ind))
	K <- P[sel.rows,]
	rbfn <- make.null.rbfn()
	rbfn$c <- K
	rbfn$s <- array(0.01,nrow(K))
	list(D=rbfn.norm.design.matrix(rbfn, P), K=K)
	}



stationary.distribution <- function(P, epsilon = 1e-5, norm = norm.1,    
				    pi = array(1/ncol(P), ncol(P)), max.it = NULL) {
## uses the power method
   if (is.null(max.it)) max.it <- nrow(P)
   dif <- Inf
   it <- 0
   while(dif > epsilon && it <= max.it) {
      pi.old <- pi
      pi <- array(pi %*% P, ncol(P))
      dif <- norm(pi-pi.old)
# 	 print(dif)
	    it <- it + 1
      }
#    print(it)
   pi
   }


direct.stationary.distribution <- function(P) {
	n <- nrow(P)
	Q <- t(P-diag(n))
	Q[n,] <- 1
	b <- array(0,n)
	b[n] <- 1
	solve(Q,b)
	}
   
 
test <- function(P, sf.function = sf.exact,
			  sd.function = stationary.distribution, ...) {
	gc()
	t1 <- system.time(v <- sd.function(P, ...), FALSE)[1]
	t2 <- system.time({
		J <- sf.function(P)
		va <- sd.function(J$K %*% J$D, ...) %*% J$K
		}, FALSE)[1]
    va <- drop(va)
    print(t2/t1)
#     print(paste("t1:", t1, "t2:", t2))
    print(paste("Error: ", norm.1(v-va)))
    }

		

      

 factor.P <- function(factor.function = sf.kl,
                        max.iter = 100,
			epsilon.ff = -Inf,
			norm = norm.1,
                        name = "mult",
                        sizes = c(100), 
                        p.centers = c(0.1,0.3,0.5,0.7,0.9), 
                        sds = c(1,3,5), 
                        noises = 1e-4,
                        perc= 0.2,
                        num.avg = 50, 
			epsilon.sd = 1e-6,
                        load.p = TRUE, 
                        save.p  = TRUE,
                        dir ="./res_paper/",
                        dir.res = "./res_paper_mc/", ... ) {
# 'factor.functions' is a list with the functions to do the fatorization
# 'max.iter' is an array with the maximum number of iterations per function
# 'sizes' is the number of states |S| of each markov chain (MC)
# 'sds' is the standard deviation to generate the transition matrix
# 'perc' is the percentual of reduction of the MP, i.e., m = perc * |S|
# 'num.avg' is the number of times each MP will be reduced
# 'names' is the strings with the names of the method (to write to file)
# 'load.p' defines whether the MC should be created or loaded from file
# 'dir' is the directory where results should be saved

   for (s in 1:length(sizes)) {
     ncenters <- round(p.centers * sizes[s]) 
     for (nc in 1:length(ncenters)) {
      print(paste("nc", nc, "/", length(ncenters)))
      for (i in 1:length(sds)) {
            print(paste("sd", i, "/", length(sds)))
            for (n in 1:length(noises)) {
            print(paste("noise", n, "/", length(noises)))
            res <- array(0, c(length(perc), num.avg, 8))
                        
            for (j in 1:num.avg) {
               cat(paste(j, " ", sep=""))
               # create (or load) MC
               P <- NULL
               filename <- paste(dir, "P_norm_s", sizes[s], "_nc",ncenters[nc],
                                 "_sd", sds[i], "_n", noises[n], "_", j, ".txt",
                                 sep = "")
               if (!load.p) {
                  P <- normal.transition.matrix(sizes[s], ncenter =
                         ncenters[nc], sd = sds[i], noise = noises[n])
                  }
               else {
                  P <- as.matrix(read.table(filename))
                  }
               

               if (save.p) wt(P,filename)
              
	       v <- stationary.distribution(P, epsilon = epsilon.sd, norm =norm)
            
	       for (p in 1:length(perc)) {
		  m <- round(perc[p] * sizes[s])
		  # print(paste("Perc.: ", perc[p], " (m =", m,")", sep=""))
		  
		  D <- factor.function(P, m, max.iter = max.iter,
			   epsilon=epsilon.ff, norm = norm)
		  
		  PA <- D$D %*% D$K
		  DIF <- P - PA
		  
		  #1-norm
		  res[p,j,1]<- norm.1(DIF)
		  res[p,j,2]<- norm.inf(DIF)
		  res[p,j,3]<- norm.frobenius(DIF)
		  res[p,j,4]<- kl(P,PA)
		  
		  va <- drop(stationary.distribution(D$K %*% D$D, 
		     epsilon=epsilon.sd, norm =norm) %*% D$K)

		  dif = v - va
		  #1-norm
		  res[p,j,5]<- norm.1(dif)
		  res[p,j,6]<- norm.inf(dif)
		  res[p,j,7]<- norm.frobenius(dif)
  		  res[p,j,8]<- kl(v,va)
		  }
	       }

            res.mat <- matrix(0, num.avg, 8)
	    for (p in 1:length(perc)) {
	       filename <- paste(dir.res, name,"_P_s", sizes[s], 
				 "_nc", ncenters[nc], "_sd", sds[i],
				 "_n", noises[n], "_p", perc[p], ".txt",
				 sep="")
	       res.mat <- res[p,,]
	       colnames(res.mat) <- c("P_norm_1",  "P_norm_inf",  "P_norm_frob",
				 "P_kullback", "pi_norm_1", "pi_norm_inf",
				 "pi_norm_frob", "pi_kullback")  
	       write.table(res.mat, filename, quote = FALSE, row.names=FALSE,
		     col.names=TRUE)  
	       }
	    }
         print(" ")
         }
        }
      }
   print("All done!!")
   }



plot.P.script.paper.mc <- function(
                        size = 100,
                        p.centers = c(0.1,0.3,0.5,0.7,0.9),
                        sds=c(1,3,5),
                        noises = 1e-4, 
                        type="P",
			cols = 1,
                        cube = res.cube(cols = cols, perc=0.2, type=type,
                                       names="mult", noises=noises,
                                       sds=sds, p.centers = p.centers,
                                       stat.function = mean,
                                       dir="./res_paper_mc/"), 
                        dir= "./fig_paper_mc/",
                        leg.pos = "topleft"
                        ){
# plots the results of factor.experiments() stored in "cube"
   
   #ylim <- c(min(cube), max(cube))
   n.centers = size * p.centers
   l <- make.leg.widths(sds)
   plot.cube.2d(cube, c(1,1,1,0,-1,1,1), x=n.centers, t="o", 
               ylab=expression("||P-DK||"[1]),
	       xlab=expression(tilde(rho)), 
               leg=l, lwd=1.5, cex=2, cex.lab=1.2, leg.pos=leg.pos,
               text.width = NULL, y.intersp = 1.2)
   filename <- paste(dir, "res_sf_P.eps", sep="")
   dev.copy2eps(file = filename)
   }



google.matrix <- function(n, alpha = 0.85) {
	num.links <- sample(0:n, n ,TRUE)
	P <- matrix(0, n, n)
	for (i in 1:n) {
		if (num.links[i] > 0) {
			links <- sample(1:n, num.links[i], FALSE)
			P[i,links] <- 1 / num.links[i]
			}
		else P[i,] <- 1/n
		}
#	alpha * P + (1-alpha)* 1/ncol(P)
	P
	}


sf.exact.google <- function(P) {

	max.rows.ind <- apply(P,2,which.max)
	sel <- unique(max.rows.ind)
	
	# 	K <- matrix(0, length(sel), ncol(P))
	# 	for (i in 1:length(sel)) {
	# 		cols <- max.rows.ind==sel[i]
	# 		K[i,cos] <- 1 / length(cols
	# 	
	P[sel,]
	}
	

 print("markov.chains.R loaded")
