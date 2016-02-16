
solve.mp.exactly <- function(R, P, df, precision = 1e-6) {
   solve((diag(nrow(P)) - df*P), R) 
   }


solve.mp.iteratively <- function(R, P, df, precision = 1e-6) 
{
   dif <- Inf
   v <- R
   
   tolerance <- precision * (1 - df) / (2 * df)
   
   while (dif > tolerance)
   {
      v2 <- R + df * P %*% v
      
      v <- v2 - v
      dif <- max(v) - min(v)
      
      v <- v2
   }
   v
}


pisf <- function(Da, K, r, df, 
                 pi = NULL, 
                 max.iter = Inf,
                 solve.mp.function = solve.mp.exactly) {
   
# Da = array(..., c(num.states, m, num.actions))
# K = matrix(..., m, num.states)
# r = array(..., m)
# df is the discount factor
# pi is the initial policy
# max.iter is the maximum number of iterations 
# solve.mp.function is the function used to compute the value of a policy
#    (this can be used to implement the "modified policy iteration algorithm",
#     which is the method recommended by Puterman (p.186))

   
   num.actions <- dim(Da)[3]
   num.states <-  ncol(K)
   m <- nrow(K)

   if (is.null(pi)) pi <- sample(1:num.actions, num.states, TRUE)

   pi.old <- array(0, num.states)
   Q <- matrix(0, num.states, num.actions)
      

   Dpi <- matrix(0, num.states, m)
   Ppi <- matrix(0,m,m)
   it <- 0

   while (it < max.iter && !identical(pi, pi.old)) {
      
      # policy evaluation
      for (i in 1:nrow(Dpi)) Dpi[i,] <- Da[i,,pi[i]]
      Ppi <- K %*% Dpi
      V <- solve.mp.function(r, Ppi, df)
    
      for (a in 1:num.actions) {
         Q[,a] <- Da[,,a] %*%  V  
         }

      pi.old <- pi
      pi <- apply(Q,1,which.max)
      
      it <- it + 1
      }
   
   list(pi = pi, Q = Q)
   }



print("pisf.R loaded")