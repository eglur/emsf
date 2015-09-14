source("data.manipulation.R")
source("rbfn.R")
source("rbfn.gaussian.R")


kbrl.mdp <- function(lsars, tau, rbfn.function = make.rbfn.simple, 
               goal.reward = 1, same.width = TRUE, verbose = FALSE, ...) {
# Generate an MDP using the standard KBRL. Does not solve it.
# The transitions in lsars associated with each action do NOT need to begin at
# the same set of states  
# lsars should be normalized

   if (verbose) print("Building the MDP. This may take a while...")
   	   
   num.actions <- lsars$num.actions
   num.states <- 0
   S <- NULL
   for (a in 1:num.actions) {
      num.states <- num.states  + nrow(lsars[[a]]$s2)
      S <- rbind(S, lsars[[a]]$s2)
      }
   
   P <- array(0, c(num.states, num.states, num.actions))   
   R <- array(0, c(num.states, num.actions))
   G <- array(0, c(num.states, num.actions))

   rbfns <- NULL

   b <- 1
   for (a in 1:num.actions) {
      rbfn <- rbfn.function(lsars[[a]]$s, tau, num.actions = 1, same.width =
                              same.width, ...)
      D <- rbfn.norm.design.matrix(rbfn, S)
      #       D[!is.finite(D)] <- 0 
#        D[is.na(D)] <- 0 ## this is a workaround for the case when states
#        D <- D/apply(D,1,sum) ##   are too close to each other                
#       
      for (i in 1:num.states) R[i,a] <- D[i,] %*% lsars[[a]]$r
      G[,a] <- goal.reward - R[,a]
      e <- b + nrow(lsars[[a]]$s) - 1
      P[,b:e,a] <- D
      rbfns <- c(rbfns, list(rbfn))
      b <- e + 1
      }
   list(S = S, P = P, R = R, G = G, rbfns = rbfns)
   }


kbrl.rbfns <- function(lsars, tau, df, rbfn.function = make.rbfn.simple, 
               goal.reward = 1, same.width = TRUE, verbose = FALSE,
               run.value.iteration = FALSE, max.iter.pi = Inf, ...){
     M <- kbrl.mdp(lsars, tau, rbfn.function, goal.reward, same.width, 
                    verbose, ...)
     Q <- NULL
     if (run.value.iteration)  Q <- value.iteration(M$R, M$P, df, G=M$G)$Q
     else Q <- policy.iteration(M$R, M$P, df, max.iter = max.iter.pi)$Q
     b <- 1
     for (a in 1:lsars$num.actions) {
          e <- b + nrow(lsars[[a]]$s) - 1
          M$rbfns[[a]]$w <- matrix(Q[b:e,a], e - b + 1, 1)
          b <- e + 1
          }
     M$rbfns
     }
 
 
     
print("kbrl.R loaded")
