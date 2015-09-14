source("rbfn.gaussian.R")
source("kbrl.R")

# Experiments to reduce an MDP using archetypes and kernels

 select.pam <- function(S, k, metric) {
# needs library("cluster")
    pam(S,k, stand = TRUE, metric = metric)
    }
 
 
 select.clara <- function(S, k, metric, stand = TRUE) { 
# needs library("cluster")
    clara(S, k, metric = metric, stand = stand, medoids.x = FALSE)$i.med
    }
	 

select.grid.mountain <- function(S, k, metric) {
#    l <- floor(sqrt(nrow(S)))
#    n <- floor(k / l)
#    ind <- NULL
#    for (i in 1:l) {
#       b <- (i-1) * l + 1
#       e <- b + l - 1
#       ind <- c(ind, seq(b,e,l=n))
#       }
   n <- sqrt(nrow(S))
   l <- floor(sqrt(k))
   cp <- seq(-1.2, 0.5, l=l)
   cv <- seq(-0.07, 0.07, l=l)
   grid <- gp(cp,cv)
   ind  <- calc.ind.2d(grid, S[n+1,1]-S[1,1], S[2,2]-S[1,2], 
                        min.pos = -1.2, min.vel = -0.07)
   
   ind <- conv.2d.1d(ind[,1], ind[,2],n) 
   c(ind, sample(1:nrow(S), k - l^2, FALSE))
   }
      
# select.pam <- function(d, k) {
#    pam(d, k, TRUE)$medoids
#    }


alloc.knn <- function(dr, tau) {
   imin <- apply(dr, 1, which.min)
   D <- matrix(0, nrow(dr), ncol(dr))
   for (i in 1:nrow(D)) D[i,imin[i]] <- 1
   D
   }
           

alloc.prop <- function(dr, tau) {
   dr <- 1/exp(dr / tau)
   dr / apply(dr, 1, sum)
   }


## Se a seleção depender do cálculo de distâncias, é melhor normalizar S antes..
arch.kernel.policy.iteration <- function(R, P, df, m, S, 
                                        pi = NULL, 
                                        max.iter = Inf,
                                        solve.mp.function = solve.mp,
                                        selection.function= select.uniform.grid,
                                        allocation.function = alloc.prop,
                                        metric = "euclidean", 
                                        verbose = FALSE, 
                                        tau = 1, 
                                        ...) {
# Policy iteration using archetype algorithm with kernels 
# R is a |S| x |A| matrix with the rewards associated with the (s,a) transitions
# P is a |S| x |S| x |A| matrix with the transition probabilities P(s_i| s_j, a)
# df is the discount factor
# m is the number of archetypes
# S is the "observational" representation of states
# pi is the initial policy
# max.iter is the maximum number of iterations 
# solve.mp.function is the function used to compute the value of a policy
#    (this can be used to implement the "modified policy iteration algorihtm",
#     which is the method recommended by Puterman (p.186))
# metric is the dissimilarity function used (see function "dist")
   
   num.actions <- ncol(R)
   num.states <- nrow(R)
   
   # determine matrices K and vectors rb
#    if (verbose) print("Computing distances...")
#    d <- dist(S, method = method)

   if (verbose) print("Selecting archetypes...")
   k <- m %/% num.actions
   m <- k * num.actions
   if (verbose) print(paste("Using", m, "archetypes"))  
   ind <- selection.function(S, k, metric,...)
			
#    ind <- C$i.med
   K <- matrix(0, m, num.states)
   rb <- array(0, m)
   for (a in 1:num.actions) {
      b <- (a-1) * k + 1
      e <- b + k - 1
      K[b:e,] <- P[ind,,a]   
      rb[b:e] <- R[ind,a]
      }
  
   # determine matrix D
   if (verbose) print("Allocating states...")
   #Do <- allocation.function(as.matrix(d)[,ind])
   
# 	## THIS SHOULD BE TEMPORARY, ALLOCATION FUNCTION IS NOT BEING USED
# 	SN <- normalize(S)
# 	rbfn <-  make.rbfn.kbrl(SN[ind,], tau) ## could be the centers
# 	Do <- rbfn.norm.design.matrix(rbfn,SN)
	## THIS IS HOW IT WAS MADE WITH THE MOUNTAIN CAR
 	d <- as.matrix(dist(normalize(S)))
   Do <- allocation.function(d[,ind],tau)
#    Do <- matrix(0, num.states, m)
#    for (i in 1:nrow(Do)) Do[i, C$clustering[i]] <- 1

   if (is.null(pi)) pi <- sample(1:num.actions, num.states, TRUE)

   pi.old <- array(0, num.states)
   Q <- matrix(0, num.states, num.actions)
      
   if (verbose) print("Running policy iteration. This may take a while...")
   if (verbose) print("Iteration: ")
   it <- 0
   while (it < max.iter && !identical(pi, pi.old)) {
      if (verbose) cat(paste(it+1," "))
      # policy evaluation
      D <- matrix(0, num.states, m)
      for (i in 1:nrow(D)) {
         b <- (pi[i] - 1) * k + 1
         e <- b + k - 1 
         D[i, b:e] <- Do[i,]
         }
      Ppi <- K %*% D
   
      V <- solve.mp.function(rb, Ppi, df)

      # policy improvement
      V <- D %*% V

      for (a in 1:num.actions) {
         Q[,a] <- R[,a] + df * P[,,a] %*% V
         }
      pi.old <- pi
      pi <- apply(Q,1,which.max)
      it <- it + 1
      }
   if (verbose) {
      print(" ")
      print("Done.")
      }
   list(pi = pi, Q = Q, mse = mse)
   }
   
     
print("arch.kernel.R loaded")       
        
   
   
   
   
   
   
   
