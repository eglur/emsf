source("nmf.R")

# Experiments to reduce an MDP used in the thesis
sa.mdp <- function(R, P, m, factor.function = nmf.kmeans, ...) {
# stochastic approximation of an MPD
# R is |S| x |A|
# P is |S| x |S| x |A|
     
     num.states <- nrow(R)
     num.actions <- ncol(R)
     
     P2 <- P[,,1]
     R2 <- as.matrix(R[,1])
     if (num.actions > 1) {
      for (i in 2:num.actions) {
            R2 <- rbind(R2, as.matrix(R[,i]))
            P2 <- rbind(P2, P[,,i])
            }
       }
       
     W <- cbind(R2,P2)
     
     Q <- factor.function(W, m,  P.only = FALSE, ...)
     m <- nrow(Q$K) # In case "factor.function" changes m
     # Now, unstack the corresponding D's
     D <- array(0, c(num.states, m, num.actions))
     for (i in 1:num.actions) {
          b <- (i-1)* num.states + 1
          e <- b + num.states - 1
          D[,,i] <- Q$D[b:e, ]
          }
     list(D=D, K = Q$K, r = Q$r)
     }
     

encode.policy <- function(pi,M) {
# encode policy 'pi' using the archetypes in 'M'
# pi is |S| x 1
# M$D is |S| x m x |A|
# M$K is  m x |S| (the archetypes)
# M$r is  m x 1
   Dpi <- matrix(0, length(pi), nrow(M$K))
   for (i in 1:nrow(Dpi)) Dpi[i,] <- M$D[i,,pi[i]]
   Ppi <- M$K %*% Dpi
   list(Dpi = Dpi, Ppi = Ppi)
   }
                 
           
arch.policy.iteration <- function(R, P, df, m, pi = NULL, max.iter = Inf,
                     solve.mp.function = solve.mp, 
                     factor.function = nmf.kmeans, 
                     max.iter.ff = 10, M = NULL, 
                     compute.MSE = FALSE, ...) {
# This is the PISF algorithm                     
# Policy iteration using archetype algorithm (Puterman, 1994, p.174)
# R is a |S| x |A| matrix with the rewards associated with the (s,a) transitions
# P is a |S| x |S| x |A| matrix with the transition probabilities P(s_i| s_j, a)
# df is the discount factor
# m is the number of archetypes
# pi is the initial policy
# max.iter is the maximum number of iterations 
# solve.mp.function is the function used to compute the value of a policy
#    (this can be used to implement the "modified policy iteration algorihtm",
#     which is the method recommended by Puterman (p.186))
# factor.function is the function used to approximate the MDP stochastically
   
   num.actions <- ncol(R)
   num.states <- nrow(R)

   if (is.null(pi)) pi <- sample(1:num.actions, num.states, TRUE)

   pi.old <- array(0, num.states)
   Q <- matrix(0, num.states, num.actions)
      
   # approximate the MDP stochastically
   if (is.null(M)) {
      M <- sa.mdp(R,P,m, factor.function, max.iter = max.iter.ff, ...)
      }
   
   m <- nrow(M$K) # In case "sa.mdp" changes m

   mse <- 0
   if (compute.MSE) {
   	sse <- 0
      for (a in 1:num.actions) {
         sse <- sse + sum((P[,,a] - M$D[,,a] %*% M$K)^2)
         sse <- sse + sum((R[,a] - M$D[,,a] %*% M$r)^2)
         }
      mse <- sse / num.states
      }
   
   it <- 0
   

   while (it < max.iter && !identical(pi, pi.old)) {
      
      # policy evaluation
      A <- encode.policy(pi,M)
      V <- solve.mp.function(M$r, A$Ppi, df, ...)
   	
      # policy improvement
      V <- A$Dpi %*% V

      for (a in 1:num.actions) {
         Q[,a] <- R[,a] + df * P[,,a] %*% V
         }
         
      pi.old <- pi
      pi <- apply(Q,1,which.max)
      it <- it + 1
      }

   list(pi = pi, Q = Q, mse = mse)
   }
 

     
print("arch.R loaded")       
        
