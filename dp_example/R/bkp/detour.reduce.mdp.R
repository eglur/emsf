# Decomposition of an MC and an MDP

source("detour.nmf.R")

reduce.mdp.block <- function(R, P, perc, factor.function, ...) {
     num.states <- nrow(R)
     num.actions <- ncol(R)
     reduced.num.states <- round(perc * num.states)
     
     R2 <- as.matrix(R[,1])
     P2 <- P[,,1]
     for (i in 2:num.actions) {
          R2 <- rbind(R2, as.matrix(R[,i]))
          P2 <- rbind(P2, P[,,i])
          }
     W <- cbind(R2,P2)
     Q <- factor.function(W, reduced.num.states,  ...)
     
     # Now, unstack the corresponding D's
     D <- array(0, c(num.states, reduced.num.states, num.actions))
     P <- array(0, c(reduced.num.states, reduced.num.states, num.actions))
     for (i in 1:num.actions) {
          b <- (i-1)* num.states + 1
          e <- b + num.states - 1
          D[,,i] <- Q$D[b:e, ]
          P[,,i] <-  Q$K %*% D[,,i]
          }
     R <- matrix(rep(Q$R, num.actions), reduced.num.states, num.actions)
     list(D =D, K = Q$K, R = R, P = P)
     }


reduce.mdp <- function(R, P, perc, factor.function, ...) {
	num.states <- nrow(R)
	num.actions <- ncol(R)
	reduced.num.states <- round(perc * num.states)
	
	R2 <- matrix(0, reduced.num.states, num.actions)
	P2 <- array(0, c(reduced.num.states, reduced.num.states, num.actions))
	L <- array(0, c(num.states, reduced.num.states, num.actions))
	K <- array(0, c(reduced.num.states, num.states, num.actions))
	
	for (i in 1:num.actions) {
		W <- cbind(R[,i], P[,,i])
		D <- factor.function(W, reduced.num.states, ...)
		R2[,i] <- D$R
		P2[,,i] <- D$K %*% D$L
		L[,,i] <- D$L
		K[,,i] <- D$K
		}
	
	list(L =L, K = K, R = R2, P = P2)
	}
	




expand.Q <- function(D, Q) 	{
	Q2 <- matrix(0, nrow(D$L), ncol(D$R))
	for (i in 1:ncol(Q2)) {
		Q2[,i] <- D$L[,,i] %*% Q[,i]
		}
	Q2
	}	




reduce.mc.parallel <- function(R, P, p, b=1, factor.function=
                        factor.nmf.mult, ...) {
# R is a |S|-dimensional vector with the rewards
# P is a |S| x |S| matrix with the transition probabilities
# 'p' is the proportion that defines the "internal dimension" m = p * b
# (m * b defines the number of states of the reduced MC)
# 'b' is the number of blocks W = cbind(R, P) should be divided into 
# (b = 1 means sequential processing; b > 1 means "parallel" processing)
	W <- cbind(R,P)
	K <- NULL
	L <- NULL
	R2 <- NULL
	
	bs <- array(0, b) # block sizes
	bs[] <- length(R) %/% b 			
     mod <- length(R) %% b
     if (mod) bs[1:mod] <- bs[1:mod] + 1     # redistribute the remaining
     
	for (i in 1:b) {
		# sel <- sample(1:nrow(W), bs[i])
		B <- W[1:bs[i],]  ## think if it's a good idea to shuffle (or cluster?) this
		W <- W[-(1:bs[i]),, drop = FALSE]
		
		m <- round(p * nrow(B))
		D <- factor.function(B, m, ...)
		
		K <- rbind(K, D$K)
		R2 <- c(R2, D$R)
		if (!is.null(L)) {
			T <- matrix(0, nrow(L), m)
			L <- cbind(L, T)
			T <- matrix(0, nrow(D$L), ncol(L) - ncol(T))
			L <- rbind(L, cbind(T, D$L))
			}
		else L <- D$L
		}
	list(L = L, K = K, R = R2)
	}
	
		

reduce.mdp.parallel <- function(R, P, p, b = 1, factor.function=
                           factor.nmf.mult, ...) {
# R is a |S| x |A| matrix with the rewards associated with the (s,a) transitions
# P is a |S| x |S| x |A| matrix with the transition probabilities P(s_i| s_j, a)
# 'p' is the proportion that defines the "internal dimension" m = p * bs (bs =block size)
# (m * b defines the number of states of the reduced MDP)
# 'b' is the number of blocks W  should be divided into 
# (here, W = cbind(rbind(R_1, R_2, ... R_|A|), rbind(P_1, P_2, ...P_|A|)))
# (b = 1 means sequential processing; b > 1 means "parallel" processing)
	
	# First, stack R's and P's 
	num.actions <- ncol(R)
	num.states <- nrow(R)
	
	p <- p / num.actions        ## think
	
	R2 <- as.matrix(R[,1])
	P2 <- P[,,1]	
	for (i in 2:num.actions) {
		R2 <- rbind(R2, as.matrix(R[,i]))
		P2 <- rbind(P2, P[,,i])
		}
	
	D <- decompose.mc(R2, P2, p, b, factor.function, ...)
	reduced.num.states <- length(D$R)
	# Now, unstack the corresponding L's
	L <- array(0, c(num.states, reduced.num.states, num.actions))
	P <- array(0, c(reduced.num.states, reduced.num.states, num.actions))
	for (i in 1:num.actions) {
		b <- (i-1)* num.states + 1
		e <- b + num.states - 1
		L[,,i] <- D$L[b:e, ]
		P[,,i] <-  D$K %*% L[,,i]
		}
	
	R <- matrix(rep(D$R, num.actions), reduced.num.states, num.actions)
	list(L =L, K = D$K, R = R, P = P)
	}
	


print("detour.reduce.mdp.R loaded")	