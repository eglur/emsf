random.transition.matrix <- function(nrow, ncol = nrow) {
## O certo seria calcular as probabilidades como a area do intervalo. 
# Generates a stochastic matrix without any column suming to zero
# (because if this were the case, the state could be removed from the MC)
	P <- matrix(runif(nrow * ncol), nrow, ncol)
	P <- P / apply(P, 1, sum)
	# detect columns whose sum is zero
	zeros <- apply(P, 2, sum) == 0
	if (sum(zeros) > 0) {
		rows <- sample(1:nrow(P), length(zeros))
		ind <- cbind(zeros, rows)
		P[ind] <- runif(nrow(ind),0.1, 1)
		P <- P / apply(P, 1, sum)
		}
	P
	}
	

normal <- function(x,mu,sigma, pi = 3.1416) {
## A funcao dnorm() é identica...
   1/(sqrt(2* pi)* sigma) * exp(-((x - mu)^2/(2 * sigma^2))) 
   }



normal.transition.matrix <- function(nrow, ncol = nrow, ncenters = min(nrow,ncol), 
                                     sd = 1, noise = 0, lim=NULL) {
# Generates a stochastic matrix in which the transitions are approximately
# normal
# 'ncenters' is the number of distinct centers of the distributions
     if (is.null(lim)) lim <- 1:ncol
     P <- matrix(0, nrow, ncol)
     ind <- sample(lim, ncenters)
     if (ncenters != 1)
     	ind <- c(ind, sample(ind, nrow - ncenters, replace = TRUE))
     else ind <- rep(ind, nrow)
		
# 	A <- matrix(rnorm(nrow * ncol, ind, sd), nrow, ncol)
# 	A <- round(A)
#      O <- A > ncol(A) || A < 1    
# 	A[O] <- sample(1:ncol, sum(O), replace = TRUE)
# 	for (i in 1:ncol(P)) P[,i] <- apply(A==i, 1, sum)
     for (i in 1:nrow(P)) {
         q <- runif(ncol,0,noise)
	 P[i,] <- normal(1:ncol, ind[i], sd) + q
         }
	P <- P / apply(P, 1, sum)
	P
	}
	
mdp.transition.matrix <- function(num.states, num.actions, 				
	 transition.matrix.function = normal.transition.matrix, ...) {
	P <- array(0, c(num.states, num.states, num.actions))
	for (a in 1:num.actions) {
		P[,,a] <- transition.matrix.function(num.states, ...)
		}
	P
	}
		
mdp.melekopoglou <- function(num.states) {
# Melekopoglou & Condon, 1990
	if (num.states %% 2 != 0) print("mdp.melekopoglou: The number of states in
the MDP must be even")
	else {
		P <- array(0, c(num.states, num.states, 2))
		R <- matrix(0, num.states, 2)
		
		odd <- seq(1, num.states-3, by=2)
		for (i in odd) {
			P[i,i+2,1] <- 1
			P[i,i+1,2] <- 1
			}
		
		even <- seq(2, num.states-4, by=2)
		for (i in even) {
			P[i,i+2,] <- 0.5
			P[i,i+1,] <- 0.5
			}
		P[num.states-1, num.states,] <- 1
		P[num.states-2, num.states,] <- 1
		P[num.states, num.states,]   <- 1
		 					
		R[num.states-1,] <- 1		 					
		list(P=P, R=R)
		}
	}


mdp.line <- function(num.states, num.actions) {
	P <- array(0, c(num.states, num.states, num.actions))
	R <- matrix(0, num.states, num.actions)
	
	for (i in 1:(num.states-1)) {
		P[i,i+1,num.actions] <- 1
		P[i,i,1:(num.actions-1)] <- 1
		}
		
	R[num.states-1, num.actions] <- 1
	P[num.states,num.states,] <- 1
	list(P=P, R=R)
	}
	

mdp.double.line <- function(num.states, num.actions) {
	P <- array(0, c(num.states, num.states, num.actions))
	R <- matrix(0, num.states, num.actions)
	
	P[,1,1:(num.actions-1)] <- 1
	
	for (i in seq(1,num.states-2,by=2)) {
		P[i,i+2,num.actions] <- 1
		P[i+1,i+2,num.actions] <- 1
		}
	P[num.states,num.states,] <- 1
	R[num.states-1, num.actions] <- 1
	list(P=P, R=R)
	}
	

mdp.triangle <- function(num.states, num.actions, prob.pass = 0.5) {
	P <- array(0, c(num.states, num.states, num.actions))
	R <- matrix(0, num.states, num.actions)
		
	for (i in (num.states-1):1) {
		P[i,i+1,num.actions] <- prob.pass
		v <- runif(i)
		P[i,1:i,num.actions] <- v / sum(v) * (1-prob.pass)
		
		for (a in 1:(num.actions-1)) {
			v <- runif(i)
			P[i,1:i,a] <- v / sum(v)
			}
		}
			
	P[num.states,num.states,] <- 1
	R[num.states-1, num.actions] <- 1
	list(P=P, R=R)
	}
		

mdp.straits <- function(num.states, num.actions, max.dif = 1, 
								prob.pass = 0.5) {
	P <- array(0, c(num.states, num.states, num.actions))
	for (i in (num.states-1):1) {
		# First guarantee the minimum probability of getting to the next level
		P[i,i+1,num.actions] <- prob.pass
		v <- runif(i)
		P[i,1:i,num.actions] <- v / sum(v) * (1-prob.pass)
		# Then, fill in the rest of the transition matrix, except the first action
		if (num.actions > 2) {
			for (a in 2:(num.actions-1)) {
			      v <- runif(i+1)
			      P[i,1:(i+1),a] <- v / sum(v)
			      }
			 }
		# Set the probabilities associated with action 1, which never goes up
		v <- runif(i)
		P[i,1:i,1] <- v / sum(v)
		}
	P[num.states,num.states,] <- 1
	
	R <- matrix(0, num.states, num.actions)
	R[num.states-1,num.actions] <- 1
	
	list(P=P, R=R)
	}

mdp.straits.no.goal <- function(num.states, num.actions, max.dif = 1, 
				prob.pass = 0.5){
	P <- array(0, c(num.states, num.states, num.actions))
	for (i in (num.states-1):1) {
	 # First guarantee the minimum probability of getting to the next level
		P[i,i+1,num.actions] <- prob.pass
		v <- runif(i)
		P[i,1:i,num.actions] <- v / sum(v) * (1-prob.pass)
		# Then, fill in the rest of the transition matrix, except the first action
		if (num.actions > 2) {
		  for (a in 2:(num.actions-1)) {
			      v <- runif(i+1)
			      P[i,1:(i+1),a] <- v / sum(v)
			      }
			 }
		# Set the probabilities associated with action 1, which never goes up
		v <- runif(i)
		P[i,1:i,1] <- v / sum(v)
		}
	
	for (a in 1:num.actions) {
	    v <- runif(num.states)
	    P[num.states,,a] <- v / sum(v)
	    }
	
	R <- matrix(0, num.states, num.actions)
	R[num.states,num.actions] <- 1
	
	list(P=P, R=R)
	}

bin2int <- function(x) sum(x * 2^(rev(seq_along(x)) - 1))

flip.bit.li <- function(x, pos) {
     if (pos < length(x) && sum(x[(pos+1):length(x)]) > 0) x[pos:length(x)] <- 1
     else if (x[pos] == 1) x[pos] <- 0
          else x[pos:length(x)] <- 1
     x
     }


flip.bit.expon <- function(x, pos) {
     if (pos < length(x)) {
          inter <- (pos+1):length(x)
          if (sum(x[inter]) == 0) x[pos] <- 0
          else x[pos] <- 1
          x[inter] <- 1
          }
     else x[pos] <- 0
     x     
     }


flip.bit.linear <- function(x, pos) {
     x[pos] <- 0
     if (pos < length(x)) x[(pos+1):length(x)] <- 1
     x     
     }


flip.highest <- function(x, flip.bit.function) {
     flipped <- FALSE
     i <- length(x)
     while (!x[i] && i > 0) i <- i - 1
     if (i > 0) flip.bit.function(x,i)
     else x
     }

flip.lowest <- function(x, flip.bit.function) {
     flipped <- FALSE
     i <- 1
     while (!x[i] && i <= length(x)) i <- i + 1
     if (i <= length(x)) flip.bit.function(x,i)
     else x
     }
     

mdp.bitflip <- function(num.bits = 10, flip.bit.function = flip.bit.li,
                        noise = 0) {
     n <- 2^num.bits
     P <- array(0, c(n,n,3))
     R <- matrix(0, n, 3)

     
     s <- array(0, num.bits)
     s[length(s)] <- 1
     P[1,1,] <- 1
     cont <- 2
      while (cont <= n) {
           # action "flip_highest"
           s2 <- flip.highest(s,flip.bit.function )
           j <- bin2int(s2) + 1
           P[cont,j, 1] <- 1 - noise
           P[cont, -j, 1] <- noise / (n-1)
           
           if (j == 1 && cont != 1) R[cont, 1] <- 1

           # action "flip_lowest"
           s2 <- flip.lowest(s, flip.bit.function )
           j <- bin2int(s2) + 1
           P[cont, j , 2] <- 1 - noise
           P[cont, -j, 2] <- noise / (n-1)
          if (j == 1 && cont != 1) R[cont, 2] <- 1

           # action "flip_random" 
           nz <- sum(s)
           for (i in 1:num.bits) {
               if (s[i]) {
                    s2 <- flip.bit.function(s,i)
                    j <- bin2int(s2) + 1
                    P[cont, j , 3] <- P[cont, j , 3] + ((1 - noise)/nz)
                    P[cont, -j ,3] <- P[cont, -j ,3] + (noise/nz)/(n-1)
                    
                    if (j == 1 && cont != 1) R[cont, 3] <- R[cont, 3] + 1/nz
                    }
               }
          cont <- cont + 1
          s <- inc.odometer(s)
          }
     list(P=P, R=R)
     }
          

mdp.bitflip.actions <- function(num.bits = 10) {
     n <- 2^num.bits
     P <- array(0, c(n,n,num.bits))
     R <- matrix(-1, n, num.bits)
     R[1,] <- 0
     
     s <- array(0, num.bits)
     cont <- 1
      while (cont <= n) {
           for (i in 1:num.bits) {
               s2 <- flip.bit(s,i)
               P[cont, j <- bin2int(s2) + 1 , i] <- 1
               }
          cont <- cont + 1
          s <- inc.odometer(s)
          }
     list(P=P, R=R)
     }

stationary.distribution <- function(P, precision = 1e-6) {
# Finds the stationary distribution of a transition matrix P
	dif <- Inf
	while (dif > precision) {
		P2 <- P %*% P
		dif <- max(abs(P2-P))
		P <- P2
		}
	apply(P, 2, mean)
	}
	
save.mdp.c <- function(P, r, filename) {
# To be loaded using C++, thus the indices
   num.states  <- dim(P)[1]
   num.actions <- dim(P)[3]
   wt(matrix(c(num.states, num.actions),1,2), paste(filename,"_info",sep=""))
   for (a in 1:num.actions) {
      wt(P[,,a], paste(filename,"_Pa",a-1, sep=""))
      wt(r[,a], paste(filename,"_ra",a-1, sep=""))
      }
   }
   
create.mdp.sf <- function(n, m, a, filename, noise = 0) {
   D <- array(0, c(n,m,a))
   K <- random.transition.matrix(m, n)
   rb <- array(runif(m), m)
   
   P <- array(0, c(n, n, a))
   r <- matrix(0, n, a)
   Pb <- array(0, c(m, m, a))
   
   for (i in 1:a) {
      D[,,i] <- random.transition.matrix(n, m)
      P[,,i] <- D[,,i] %*% K + runif(n^2, 0, noise)
      P[,,i] <- P[,,i] / apply(P[,,i], 1, sum)
      r[,i]  <- D[,,i] %*% rb  + runif(n, -noise, noise)
      Pb[,,i] <- K %*% D[,,i]
      wt(D[,,i], paste(filename, "_D", i, sep=""))
      }
  wt(K, paste(filename, "_K", sep=""))
  wt(rb, paste(filename, "_rb", sep=""))
  save.mdp.c(P, r, paste(filename,"_mdp",sep=""))
#   save.mdp.c(Pb, rb, paste(filename,"_mdpb",sep=""))
  }
  
   
      
      
   
   
   



print("mdp.R loaded")