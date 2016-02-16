

add.column.pseudo.inv <- function(A, PA, x, epsilon = 1e-9) {
# %           A: n x m matrix
# %           PA = m x n (pseudoinverse(A))
# %           x column vector, n x 1 to append to A 
     n <- nrow(A)
     m <- ncol(A)

     PAx <- PA %*% x  # O(nm)
     P <- diag(n) - A %*% PA  # O(n^2 m)
     Px <-  P %*% x   # O(n^2)
     alpha <- drop(t(x) %*% Px) # O(n)
     b <- NULL
     if (alpha < epsilon) {
          eta = t(PAx) %*% PAx # O(m)
          b = PA %*% (PAx/(1+eta)) # O(nm) 
          }
     else b = Px/alpha; # O(n)

     rbind(PA - PAx %*% t(b), t(b)) # O(m*n)
     }

print("add.pseudo.inverse.R loaded")
