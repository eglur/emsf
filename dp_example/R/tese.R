testa.hipotese <- function(num.states, num.actions, df = 0.9, epsilon = 1e-10) {
   
   P <- mdp.transition.matrix(num.states, num.actions)
   R <- matrix(runif(num.states * num.actions), num.states, num.actions)
   coef <- matrix(0, num.states-1, num.actions)
   
   #t <- runif(num.states-1)
   #t <- t / sum(t)
   for (a in 1:num.actions) {
      #coef[,a] <- t
      coef[,a] <- runif(num.states-1)
      coef[,a] <- coef[,a] / sum(coef[,a])
      
      R[1,a] <- coef[,a] %*% R[2:num.states, a]
      for (i in 1:num.states) {
         P[1,i,a] <- coef[,a] %*% P[2:num.states, i, a]
         }
      }
   # okay, MDP ready, let's go!
   Q <- value.iteration(R,P, df, epsilon = epsilon)$Q
   
          
   # reduced MDP
   P2 <- array(0, c(num.states-1, num.states-1, num.actions))
   R2 <- matrix(0, num.states - 1, num.actions)
   L <- diag(num.states)
   L[1,1] <- 0
   for (a in 1:num.actions) {
      L[1, 2:num.states] <- coef[,a]
      T <- P[,,a] %*% L
      P2[,,a] <- T[2:num.states, 2:num.states]
      R2[,a] <- R[-1,a]
      
      }
   #print(P)
   #print(P2)
   #print(R)
   #print(R2)
   
   Q2 <- value.iteration(R2, P2, df, epsilon = epsilon)$Q
   Q3 <- matrix(0, num.states, num.actions)
   for (a in 1:num.actions) {
      Q3[,a] <- c(coef[,a] %*% Q2[,a], Q2[,a])
      }
   Q3 <- Q3 - mean(Q3)
   Q <- Q - mean(Q)
   print(sum((Q3 - Q)^2))
   p <- apply(Q, 1, which.max)
   p3 <- apply(Q3, 1, which.max)
   print(sum(p != p3))
   #matplot(cbind(p,p3), t="l")
   #print(Q)
   #print(Q3)
   }
       

teste.teorema <- function(n, k, df = 0.9) {
   P <- matrix(0, n, n)
   R <- array(0, n)
   
   P[1:k,] <- runif(n*k)
   P[1:k,] <- P[1:k,] / apply(P[1:k,], 1, sum)
   R[1:k] <- runif(k)
   
   L <- matrix(0, n, n)
   for (i in 1:k) L[i,i] <- 1
   
   for (i in (k+1):n) {
      L[i,] <- c(runif(k), rep(0, n-k))
      L[i,] <- L[i,] / sum(L[i,])
      P[i,] <- L[i,] %*% P
      R[i] <- L[i,] %*% R
      }
   
   Q <- value.iteration(R, P, df)$Q
   
   P2 <- P %*% L
   Q2 <- value.iteration(R, P2, df)$Q
   
   matplot(cbind(Q,Q2), t="l")
   sum((Q - Q2)^2)
   }
   
   
teste.teorema2 <- function(n, df = 0.9) {
   P <- normal.transition.matrix(n)
   R <- runif(n)
   
   L <- matrix(runif(n^2), n, n)
   #L <- L / apply(L, 1, sum)
   
   LI <- solve(L)
   K <- LI %*% P
   R2 <- LI %*% R
   
   Q <- solve.mc(R, P, df)
   
   P2 <- K %*% L
   
   Q2 <- solve.mc(R2, P2, df)
   Q2 <- L %*% Q2
   
   #Q <- Q - mean(Q)
   #Q2 <- Q2 - mean(Q2)
   
   matplot(cbind(Q,Q2), t="l")
   sum((Q - Q2)^2)
   }
   
   
teste.teorema.mdp <- function(num.states, k, num.actions, df = 0.9) {
# testa se a política ótima para um MDP reduzido com bases distintas para cada
# ação é a mesma do MDP original
# CONCLUSÃO: NÃO É!!!
   P <- array(0, c(num.states, num.states, num.actions))
   R <- matrix(0, num.states, num.actions)
   for (a in 1:num.actions) {
      #P[1:k,,a] <- runif(k * num.states)
      #P[1:k,,a] <- P[1:k,,a] / apply(P[1:k,,a],1,sum)
      for (i in 1:k) P[i,sample(1:num.states,1),a] <- 1
      R[1:k,a] <- runif(k)
      }
  
   K <- P[1:k,,]
   R2 <- R[1:k,]
   
   L <- array(0, c(num.states, k,num.actions))
   for (a in 1:num.actions) {
      for (i in 1:k) L[i,i,a] <- 1
      for (i in (k+1):num.states) {
         L[i,,a] <- runif(k)
         L[i,,a] <- L[i,,a] / sum(L[i,,a])
         P[i,,a] <- L[i,,a] %*% P[1:k,,a]
         R[i,a]  <- L[i,,a] %*% R[1:k,a]
         }
     }
     
   
   P2 <- array(0, c(k,k,num.actions))
   for (a in 1:num.actions) {
      P2[,,a] <- K[,,a] %*% L[,,a]
      }
   
   Q  <- policy.iteration(R,   P, df)$Q
   Qt <- policy.iteration(R2, P2, df)$Q
   
   Q2 <- matrix(0, num.states, num.actions)
   
   for (a in 1:num.actions) Q2[,a] <- L[,,a] %*% Qt[,a]
    
   print(sum((Q-Q2)^2))
   p <- apply(Q,1,which.max)
   p2 <- apply(Q2,1,which.max)
   print(sum(p != p2))
   
   plot.both(p,p2)
   }
   
   
teste.teorema.mdp.base.unica <- function(num.states, k, num.actions, df = 0.9) {
# testa se a política ótima para um MDP reduzido com um base comum para todas
# as ações é a mesma do MDP original
# CONCLUSÃO: NOPE!...
   P <- array(0, c(num.states, num.states, num.actions))
   R <- matrix(0, num.states, num.actions)
   K <- matrix(0, k, num.states)
   RB <- array(0, k)
   for (i in 1:k) K[i,sample(1:num.states,1)] <- 1
   RB <- runif(k)
  
   L <- array(0, c(num.states, k,num.actions))
   for (a in 1:num.actions) {
      for (i in 1:num.states) {
         L[i,,a] <- runif(k)
         L[i,,a] <- L[i,,a] / sum(L[i,,a])
         P[i,,a] <- L[i,,a] %*% K
         R[i,a]  <- L[i,,a] %*% RB
         }
     }
     
   
   P2 <- array(0, c(k,k,num.actions))
   R2 <- matrix(0, k, num.actions)
   for (a in 1:num.actions) {
      P2[,,a] <- K %*% L[,,a]
      R2[,a] <- RB
      }
   
   Q  <- policy.iteration(R,   P, df)$Q
   Qt <- policy.iteration(R2, P2, df)$Q
   
   Q2 <- matrix(0, num.states, num.actions)
   
   for (a in 1:num.actions) Q2[,a] <- L[,,a] %*% Qt[,a]
    
   print(sum((Q-Q2)^2))
   p <- apply(Q,1,which.max)
   p2 <- apply(Q2,1,which.max)
   print(sum(p != p2))
   
   plot.both(p,p2)
   }
   