
compare.sf.strategies <- function(n,m=n,a=5, df=0.9, ncenters = floor(0.2*n), sd
= 1) {
     D <- normal.transition.matrix(n, m, ncenters = ncenters, sd =sd) 
     K <- array(0, c(m,n,a))
     Rb <- matrix(runif(m*a), m, a)
     Rb[1,1] <- 1

     P <- array(0, c(n,n,a))
     R <- matrix(0, n, a)
     Pb <- array(0, c(m,m,a))

#      Pg <- array(0, c(m,m,a))

     for (i in 1:a) {
         K[,,i]  <- normal.transition.matrix(m, n, ncenters = ncenters, sd =sd)
    
          P[,,i]  <- D %*% K[,,i] 
          R[,i] <- D %*% Rb[,i]
          Pb[,,i] <- K[,,i] %*% D
#           Pg[,,i] <- P[,,i] %*% D
          }

     pi  <- policy.iteration(R, P, df)$pi

     Q <- D %*% policy.iteration(Rb, Pb, df)$Q
     pib1 <- apply(Q,1,which.max)

     V <- D %*% apply(policy.iteration(Rb, Pb, df)$Q, 1, max)
     Q <- matrix(0,n,a)
     for (i in 1:a) {
          Q[,i] <- R[,i] + df * P[,,i] %*% V
          }
     pib2 <- apply(Q,1,which.max)
     
     H <- matrix(0, m, n)
     H[1:m, 1:m] <- diag(m)

     Q <- t(H) %*% policy.iteration(Rb, Pb, df)$Q
     pib3 <- apply(Q,1,which.max)

     V <- t(H) %*% apply(policy.iteration(Rb, Pb, df)$Q, 1, max)
     Q <- matrix(0,n,a)
     for (i in 1:a) {
          Q[,i] <- R[,i] + df * P[,,i] %*% V
          }
     pib4 <- apply(Q,1,which.max)

     list(sum(pi!=pib1), sum(pi!=pib2), sum(pi!=pib3), sum(pi!=pib4))
     }