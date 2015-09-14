source("dp.R")
source("data.manipulation.R")
source("rbfn.gaussian.R")
source("util.R")

make.rbfn.kbrl <- function(centers, num.actions, tau, num.neighbors =
                     ncol(centers) + 1) {
	rbfn <- make.null.rbfn(act.function = gaussian, norm = sqr.euclidean.norm)
	
	# set the centers
	num.rbfs <- nrow(centers)
	rbfn$c <- centers
	rbfn$w <- matrix(rnorm(num.rbfs * num.actions), num.rbfs, num.actions)
	
	# set the widths
	dists <- matrix(0, num.rbfs, num.rbfs)
	for (i in 1:(num.rbfs-1)) {
		for (j in (i+1):num.rbfs) {
			dists[i,j] <- rbfn$norm(rbfn$c[i,] - rbfn$c[j,])
			dists[j,i] <- dists[i,j]
			}
		}
	for (i in 1:num.rbfs) {
		dists[i,i] <- Inf
		rbfn$s[i] <- -mean(sort(dists[i,])[1:(num.neighbors)]) / log(tau) 
		}
	rbfn
	}
	


detour.rbfn <- function(lsars, tau, df, prop.centers = 1, center.function =
                  center.grid, num.neighbors = ncol(lsars[[1]]$s) + 1, 
                  goal.reward = 1, epsilon.vi = 1e-6, iter.pi = 10, 
                  run.value.iteration = TRUE, ...) {
# Notice that here 'prop.centers' refers to the proportion of all the data
    
     num.actions <- lsars$num.actions
     
     S <- NULL
     for (a in 1:num.actions) S <- rbind(S, lsars[[a]]$s)
     
     num.centers <- round(prop.centers * nrow(S))
     C <- center.function(lsars, num.centers,...)
     rbfn <- make.rbfn.kbrl(C, num.actions, tau, num.neighbors)
     
     K <- NULL
     for (a in 1:num.actions) {
         H <- t(rbfn.norm.design.matrix(rbfn, lsars[[a]]$s))
         H <- H / apply(H, 1, sum)
         K <- c(K, list(H)) 
         }
     
     R <- matrix(0, num.centers, num.actions)
     P <- array(0, c(num.centers, num.centers, num.actions))
     G <- matrix(1, num.centers, num.actions)
     
     for (a in 1:num.actions) {
          R[,a] <- K[[a]] %*% lsars[[a]]$r
          P[,,a] <- K[[a]] %*% rbfn.norm.design.matrix(rbfn, lsars[[a]]$s2)
          G[,a] <- goal.reward - R[,a]
          }
               
      if (run.value.iteration) {   
         rbfn$w <- value.iteration(R, P, G = G, epsilon = epsilon.vi, df = df)$Q
         }
      else rbfn$w <- policy.iteration(R, P, df = df, iter.max = iter.pi)$Q  
      
     rbfn
     }



detour.rbfn.same <- function(lsars, tau, df, prop.centers = 1, center.function =
                  center.grid, num.neighbors = ncol(lsars[[1]]$s) + 1, 
                  goal.reward = 1, epsilon.vi = 1e-6, iter.pi = 10, 
                  run.value.iteration = TRUE, ...) {
# IMPORTANT: lsars MUST HAVE the same initial states for each a!!
# this version is faster if the initial states are the same, and can also be 
# used to run the standard KBRL when prop.centers == 1.
# Notice that here 'prop.centers' refers to the proportion of
# nrow(lsars[[1]]$s), while in the function 'detour.rbfn' this parameter refers
# to the concatenation of all lsars[[a]]$s
    
     S <- lsars[[1]]$s
## does it make sense to have more centers than points (prop.centers >1) ? 
	num.centers <- round(prop.centers * nrow(S))
     C <- NULL
     
     if (prop.centers == 1) C <- S
     else C <- center.function(lsars, num.centers,...)
    
	num.actions <- lsars$num.actions
	
	rbfn <- make.rbfn.kbrl(C, num.actions, tau, num.neighbors)
	
     K <- NULL
	if (prop.centers == 1) K <- diag(nrow(S)) # just to guarantee 
	else {
        K <- t(rbfn.norm.design.matrix(rbfn, S)) ## normalized?
        K <- K / apply(K, 1, sum)
        }
	
	R <- matrix(0, num.centers, num.actions)
	P <- array(0, c(num.centers, num.centers, num.actions))
	G <- matrix(1, num.centers, num.actions)
	
	for (a in 1:num.actions) {
  		R[,a] <- K %*% lsars[[a]]$r
  		P[,,a] <- K %*% rbfn.norm.design.matrix(rbfn, lsars[[a]]$s2)
		G[,a] <- goal.reward - R[,a]
		}
			
      if (run.value.iteration) {	
         rbfn$w <- value.iteration(R, P, G = G, epsilon = epsilon.vi, df = df)$Q
         }
      else rbfn$w <- policy.iteration(R, P, df = df, iter.max = iter.pi)$Q  
      
 	rbfn
	}




center.random <- function(lsars, num.centers, perc.goals= 0.05) {
# lsars MUST HAVE the same initial states for each a!!
	# first get transitions that end up in a goal ##  necessario?
	#goals <- lsars[[1]]$g
	#for (a in 2:lsars$num.actions) goals <- goals | lsars[[a]]$g
	
	#G <- lsars[[1]]$s[goals, ]
	C <- NULL
	#ng <- max(round(perc.goals * num.centers),1)
	#if (nrow(G) <= ng) C <- G
	#else C <- G[sample(1:nrow(G), ng), , drop = FALSE]
	#num.centers <- num.centers - nrow(C)
   
     rbind(C, lsars[[1]]$s[sample(1:nrow(lsars[[1]]$s), num.centers),])
     }
	

center.kmeans <- function(lsars, num.centers, iter.max = 10, perc.goals= 0.05){
# num.centers is the number of centers
# perc.goals is the proportion of points saved to be goals
# lsars MUST HAVE the same initial states for each a!!
	S <- lsars[[1]]$s
	
	# first get transitions that end up in a goal ##  necessario?
	#goals <- lsars[[1]]$g
	#for (a in 2:lsars$num.actions) goals <- goals | lsars[[a]]$g
	
   	#G <- S[goals, ]
	C <- NULL
	#ng <- max(round(perc.goals * num.centers),1)
	#if (nrow(G) <= ng) C <- G
	#else C <- G[sample(1:nrow(G), ng), , drop = FALSE]
	#num.centers <- num.centers - nrow(C)
   
  	# make the set S2 with terminal states
	S2 <- lsars[[1]]$s2
	for (a in 2:lsars$num.actions) S2 <- rbind(S2, lsars[[a]]$s2)
	
	C <- rbind(C,kmeans(S2, num.centers, iter.max = iter.max)$centers)
	}
	
	
center.grid <- function(lsars, num.centers, perc.goals= 0.05){
# num.centers is the number of centers
# perc.goals is the proportion of points saved to be goals
# lsars MUST HAVE the same initial states for each a!!
	S <- lsars[[1]]$s
	
	# first get transitions that end up in a goal ##  necessario?
	#goals <- lsars[[1]]$g
	#for (a in 2:lsars$num.actions) goals <- goals | lsars[[a]]$g
	
	#G <- S[goals, ]
	C <- NULL
	#ng <- max(round(perc.goals * num.centers),1)
	#if (nrow(G) <= ng) C <- G
	#else C <- G[sample(1:nrow(G), ng), , drop = FALSE]
	#num.centers <- num.centers - nrow(C)

	# make the set S2 with terminal states
	C2 <- matrix(0, num.centers, ncol(S))
	max.X <- apply(S, 2, max)
	min.X <- apply(S, 2, min)
	rbind(C,uniform.grid(C2, max.X, min.X))
 	}
	
	
## REVIEW FROM THIS POINT DOWN

cs.fuzzy <- function(k, lsars, df, m = 2, method = "euclidean", norm =
               euclidean.norm, goal.reward = 1, p = 0.05, epsilon = 1e-5, ...) {
# 'k' is the desired number of clusters
# 'm' is the 'fuzzification' parameter
# 'norm' should agree with 'method'---"euclidean" here is not squared
# 'p' is the percentage of goal states 
# lsars MUST HAVE the same initial states for each a!!
	S <- lsars[[1]]$s
	
	# first get transitions that end up in a goal ## necessario?
	#goals <- lsars[[1]]$g
	#for (a in 2:lsars$num.actions) goals <- goals | lsars[[a]]$g
	
	#G <- S[goals, , drop = FALSE]
	#ng <- max(round(p * k),1)
 	#if (nrow(G) <= ng) C <- G
	#else C <- G[sample(1:nrow(G), ng), , drop = FALSE]
	#k <- k - nrow(C)

	# make the set S2 with terminal states
	S2 <- lsars[[1]]$s2
	for (a in 2:lsars$num.actions) S2 <- rbind(S2, lsars[[a]]$s2)
	
	# cluster states in S2 using fuzzy C-MEANS
	#D <- dist(S2, method = method,...)
	F <- cmeans(S2, k, m = m)
	CL <- F$centers
	
	# compute the clusters' centers CL 
	# This is for FANNY	
	# H2 <- t(F$membership)^m ## nao tenho certeza de que m=2 no algoritmo fanny
	# H2 <- H2 / apply(H2, 1, sum)  
	# CL <- H2 %*% S2
	#CL <- rbind(C, CL) # add centers that correspond to a goal
	
	# define constants to make code more readable
	num.actions <- lsars$num.actions
	num.clusters <- nrow(CL)
	
	# compute the memberships M of the elements in S
	M <- t(fuzzy.membership(CL, S, m, norm))
	M <- M / apply(M, 1, sum)
	
	# compute the memberships H of the elements in S2
	H <-  array(0, c(nrow(S), num.clusters, num.actions))
	#M2 <- fuzzy.membership(CL, S2, m, norm) # this is only necessary because of the goal states; otherwise I could use F$membership
	#M2 <- fuzzy.membership(CL, S2, m, norm)
	
	#for (a in 1:num.actions) {
	#	b <- (a-1)*nrow(S)+1
	#	e <- b + nrow(S) - 1
	#	H[,,a] <- M2[b:e,]
	#	}
	
	# now build P, R and G	
	R <- matrix(0, num.clusters, num.actions)
	P <- array(0, c(num.clusters, num.clusters, num.actions))
	G <- matrix(1, num.clusters, num.actions)
	
	for (a in 1:num.actions) {
		R[,a] <- M %*% lsars[[a]]$r
		G[,a] <- goal.reward - R[,a]
		P[,,a] <- M %*% fuzzy.membership(CL, lsars[[a]]$s2,m , norm)
		}
		
	Q <- value.iteration(R, P, df, epsilon = epsilon, G = G)$Q
	list(c = CL, w = Q, m = m)	
	}
	
	
fuzzy.membership <- function(CL, S, m = 2, norm = euclidean.norm) {
# compute the memberships M of the elements in S
# CL contains the clusters' centers
# S are the data points
# 'm' is the 'fuzzification' coefficient
	D <- matrix(0, nrow(S), nrow(CL))
	for (i in 1:nrow(CL)) {
		T <- apply(S, 1, "-", CL[i,])
		D[,i] <- apply(T, 2, norm) 
		 }
	D <- D^(2/(m-1)) # d_ij is the distance from s_j to cl_i 
	Z <- D == 0 # to handle zeros
	ZR <- apply(Z, 1, sum) # are there any zeros at row 'i'?
	M <- matrix(0, nrow(S), nrow(CL))
	for (i in 1:nrow(M)) {
		if (ZR[i]==0) {
			for (j in 1:ncol(M)) {
				M[i,j] <- 1/sum(D[i,j] / D[i,])
				}
			}
		else   {
			M[i,Z[i,]] <- 1 / ZR[i]
			M[i,!Z[i,]] <- 0
			}
		}
	M
	}
		

print("detour.samples.R loaded")
			