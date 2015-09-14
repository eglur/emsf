value.iteration.4D.Q <- function(min.s, max.s, n, A, transition.function,  df
                        = 0.9, epsilon = 1e-6, max.iterations = 10000, 
                        border.penalty = 0,  ...) {
# This only works in deterministic 4D environments
# returns Q instead of V
# S is a n x 4 array, where n is the number of partitions
# A is the action space
# epsilon is the maximum difference on the value function between
# two consecutive iterations
# df is the discount factor
   delta <- (max.s - min.s) / (n-1)
   
   # build the matrix with next states
   S.next <- array(0, c(n,n,n,n,length(A),4)) # keep the INDICES of
                                   # the next states (not the states themselves)
   r.next <- array(0, c(n,n,n,n,length(A)))
   
   # keep track of goal and illegal states
   goal <- array(FALSE, c(n,n,n,n))
    
   for (i in 1:n) {
      for (j in 1:n) {
         for (k in 1:n) {
           for (l in 1:n) {
               for (a in 1:length(A)) {
                  s <- min.s + c(i,j,k,l) * delta
                  t <- transition.function(s,A[a], ...)
                  S.next[i,j,k,l,a,] <- pmin(pmax((s-min.s) %/% delta + 1,1), n)
                  r.next[i,j,k,l,a] <- t$r
                  if (t$g) goal[S.next[i,j,k,l,a,]] <- TRUE
                  }
               }
           }
         }
      }
      
   Q <- array(0, c(n,n,n,n,length(A)))
   tolerance <- epsilon * (1 - df) / (2 * df)
   max.dif <- Inf
   iteration <- 0
   while (max.dif > tolerance && iteration < max.iterations) {
      print(paste("Iteration ", iteration + 1))
      max.dif <- -Inf
      for (i in 1:n) {
         for (j in 1:n) {
            for (k in 1:n) {
               for (l in 1:n) {
                     for (a in 1:length(A)) {
                        Q.old <- Q[i,j,k,l,a]
                        V.next <- 0
                        if (!goal[i,j,k,l]) {
                           ind <- S.next[i,j,k,l,a,]
                           V.next <- max(Q[ind[1], ind[2], ind[3],ind[4],])
                           }
                        Q[i,j,k,l,a] <- r.next[i,j,k,l,a] + df * V.next
                        dif <- abs(Q[i,j,k,l,a] - Q.old)
                        if (dif > max.dif) max.dif <- dif
                        }
                    }
                 }
              }
           }
        iteration <- iteration + 1           
        }
   Q
   }

print("dp.4D.R loaded")