# Functions related to dynamic programming

solve.mp <- function(R, P, df, precision = 1e-6) {
# solve a Markov process exactly 
# 	P <- df * P
# 	diag(P) <- diag(P) - 1
# 	solve(P,-R)

   solve((diag(nrow(P)) - df*P), R) 
   # singular value decomposition
#    A <- svd(diag(nrow(P)) - df*P)
#    # remove the small singular values
#    d.min <- max(A$d) * precision
#    A$d[A$d < d.min] <- 0
#    non.zero.ind <- A$d != 0
#    A$d[non.zero.ind] <- 1 / A$d[non.zero.ind]
#    A$d <- diag(A$d)
#    # solve the system   
#    A$v %*% A$d %*% t(A$u) %*% R
	}


value.iteration <- function(R, P, df, iter.max = Inf, 
			epsilon = 1e-6, span =TRUE, 
			G = matrix(1, nrow(R), ncol(R)), Q = NULL) {
# Value iteration algorithm as suggested by Puterman, 1994, p. 205
# R is a |S| x |A| matrix with the rewards associated with the (s,a) transitions
# P is a |S| x |S| x |A| matrix with the transition probabilities P(s_i| s_j, a)
# df is the discount factor
# iter.max is the maximum number of iterations 
# epsilon defines the quality of the final policy, which will be epsilon-optimal
# span: should the 'span stop criterion' be used? (Puterman, 6.6.1)
# 	(generally, if span == TRUE the algorithm converges faster)
# G is a |S| x |A| a matrix with transitions that end up in a goal position
# Q is the initial action-value function

	#set up the vectors if it's a MC (not a MDC)
	if (is.na(dim(P)[3])) { 
		dim(R) <- c(length(R), 1)
		dim(P) <- c(nrow(P), ncol(P), 1)
		if (!is.null(Q)) dim(Q) <- c(length(Q), 1)
		}
	if (is.null(Q)) Q <- R

	V <- apply(Q, 1, max)
	it <- 0
     max.dif <- Inf
	tolerance <- epsilon * (1 - df) / (2 * df)
	while (max.dif > tolerance && it <= iter.max) {
		V.old <- V
		for (a in 1:ncol(Q)) {
			V.next <- df * P[, , a] %*% V
			V.next <- V.next * G[,a]
			Q[,a] <- R[,a] + V.next
			}
		V <- apply(Q, 1, max)

		V.old <- V - V.old
		if (span) max.dif <- max(V.old) - min(V.old)
		else max.dif <- max(abs(V.old)) 
		
		it <- it + 1
      	}
	list(Q = Q, it = it)
	}
	
	
policy.iteration <- function(R, P, df, pi = NULL, max.iter = Inf,
                     solve.mp.function = solve.mp, 
                     iterative = FALSE, ...) {
# Policy iteration algorithm as suggested by Puterman, 1994, p. 174
# R is a |S| x |A| matrix with the rewards associated with the (s,a) transitions
# P is a |S| x |S| x |A| matrix with the transition probabilities P(s_i| s_j, a)
# df is the discount factor
# pi is the initial policy
# max.iter is the maximum number of iterations 
# solve.mp.function is the function used to compute the value of a policy
#    (this can be used to implement the "modified policy iteration algorihtm",
#     which is the method recommend by Puterman (p.186))
   
   num.actions <- ncol(R)
   num.states <- nrow(R)

   if (is.null(pi)) pi <- sample(1:num.actions, num.states, TRUE)
   PIS <- pi
   
   pi.old <- array(0, num.states)
   P.pi <- matrix(0, num.states, num.states)
   R.pi <- array(0, num.states)
   Q <- matrix(0, num.states, num.actions)
   
   it <- 0
   while (it < max.iter && !identical(pi, pi.old)) {
	   # policy evaluation
      for (i in 1:num.states) {
         R.pi[i] <- R[i,pi[i]]
         P.pi[i, ] <- P[i,, pi[i]]
         }
      V <- solve.mp.function(R.pi, P.pi, df, ...)

   # This version doesn't use tensors
      # policy improvement
	for (a in 1:num.actions) {
	      Q[,a] <- R[,a] + df * P[,,a] %*% V
	       }
		
		# This version uses tensors
#       Q <- R + df * matrix(tensor(P, V, 2, 1), nrow(Q), ncol(Q))
      pi.old <- pi
      pi <- apply(Q,1,which.max)
      it <- it + 1

      PIS <- rbind(PIS, pi)
   	}
   list(pi = pi, Q = Q, PIS = PIS)
   }
	

load.mdp <- function(na = 11, transition.matrix.name = "P", reward.matrix.name = "R", dir = "~/tmp/dp/")
{
   # first check number of states
   P <- read.table(paste(dir,transition.matrix.name, "0.txt", sep = ""))
   ns <- nrow (P)
   
   print(paste("Detected", ns, "states"))
   
   # load matrices
   P <- array(0, c(ns, ns, na))
  
   for (a in 1:na)
   {
      P[,,a] <- as.matrix(read.table(paste(dir,transition.matrix.name, a-1, ".txt", sep = "")))
   }
   
   # load rewards
   r <- matrix(0, ns, na)
   T <- as.matrix(read.table(paste(dir,reward.matrix.name, ".txt", sep = "")))
   for (a in 1:na)
   {
      r[,a] <- T[a,]
   }

   list(P = P, r = r)   
}
	

value.iteration.finite.hor <- function(R, P, hor) {
# Value iteration algorithm as suggested by Puterman, 1994, p. 205
# R is a |S| x |A| matrix with the rewards associated with the (s,a) transitions
# P is a |S| x |S| x |A| matrix with the transition probabilities P(s_i| s_j, a)
# Q is the initial action-value function

     ns <- nrow(R)
     na <- ncol(R)
     
     pi <- matrix(0, ns, hor) 
                  
     V <- 0:(ns - 1)
    
     for (i in 1:hor) 
          {
           print(paste("Running iteration", i))
          
          Q <- matrix(0, ns, na)
          
          for (a in 1:na) 
          {
             V.next <- P[, , a] %*% V
             Q[,a] <- R[,a] + V.next
          }
          
          V <- apply(Q, 1, max)
          pi[,i] <- apply(Q, 1, which.max)
          
          }
          
     list(pi = pi, V = V)
     }


print("dp.R loaded")