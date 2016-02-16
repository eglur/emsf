source("mdp.R")
source("dp.R")


test.individual.factorization <- function(n = 100, A =2:10, num.avg = 50,
                                          m=20, df=0.9){
# Here the MDP is factored all at once, but instead of using the PISF algorithm
# the matrices Pba and the  vectors rba are generated isolated.
# Even with a exact factorization, the policy found is not optimal

	res <- matrix(0, num.avg, length(A))
	
	for (i in 1:length(A)) {
		na <- A[i]
		
      Pb <- array(0, c(m,m,na))
      
      P <- array(0, c(n,n,na))
      r <- matrix(0, n, na)
      
      D <- array(0, c(n,m,na))
      
      for (j in 1:num.avg) {
         K <-  normal.transition.matrix(m, n, m, sd = 1, noise = 0)
         rb <- runif(m)
         
         for (a in 1:na) {
            D[,,a] <- matrix(runif(n * m), n, m)
            D[,,a] <- D[,,a] / apply(D[,,a], 1, sum)
            
            P[,,a] <- D[,,a] %*% K 
            r[,a] <- D[,,a] %*% rb
            
            Pb[,,a] <- K %*% D[,,a]
            
            }
  			       
  			
  			rb <- matrix(rep(rb,na),m,na)
        
         pi  <- policy.iteration(r, P, df)$pi
         Qb <- policy.iteration(rb, Pb, df)$Q
         
         Q <- matrix(0, n, na)
         for (a in 1:na) {
         	Q[,a] <- D[,,a] %*% Qb[,a]
         	}
         pib <- apply(Q, 1, which.max)
         
         res[j,i] <- sum(pi != pib) / n
         }
     }
  res
  }
         
plot.individual.factorization <- function(res, A = 2:10,
									eps.filename = "./fig_tese/ind_fat.eps") {
   res <- apply(res, 2, mean)
   set.par.tese()
   plot(res, t="o", ylab=expression(chi[m]),x = A,
      xlab="|A|", lwd=2, cex=2, pch = 20, ps=20)

   dev.copy2eps(file = eps.filename)
   }
	

test.individual.factorization2 <- function(n = 100, A =10,
												M=seq(10,50,by=10),num.avg = 50, df=0.9){
# Like "test.individual.factorization", but here the number os archetypes varies
# instead of the number of actions 
	
	res <- matrix(0, num.avg, length(M))
   na <- A[1]	
	
	for (i in 1:length(M)) {

		m <- M[i]
			
      Pb <- array(0, c(m,m,na))
      rb <- matrix(0, m, na)
      
      P <- array(0, c(n,n,na))
      r <- matrix(0, n, na)
      
      D <- array(0, c(n,m,na))
      rb <- matrix(0, m, na)
      
      for (j in 1:num.avg) {
         K <-  normal.transition.matrix(m, n, m, sd = 1, noise = 0)
         
         for (a in 1:na) {
            D[,,a] <- matrix(runif(n * m), n, m)
            D[,,a] <- D[,,a] / apply(D[,,a], 1, sum)
            rb[,a] <- runif(m)
            
            P[,,a] <- D[,,a] %*% K 
            r[,a] <- D[,,a] %*% rb[,a]
            
            Pb[,,a] <- K %*% D[,,a]
            
            }
  			       
  				
         pi  <- policy.iteration(r, P, df)$pi
         Qb <- policy.iteration(rb, Pb, df)$Q
         
         Q <- matrix(0, n, na)
         for (a in 1:na) {
         	Q[,a] <- D[,,a] %*% Qb[,a]
         	}
         pib <- apply(Q, 1, which.max)
         
         res[j,i] <- sum(pi != pib) / n
         }
     }
  res
  }
         

plot.individual.factorization2 <- function(res, M=seq(10,50,by=10),
									eps.filename = "./fig_tese/ind_fat.eps") {
   res <- apply(res, 2, mean)
   set.par.tese()
   plot(res, t="o", ylab=expression(chi[m]),x = M,
      xlab="m", lwd=2, cex=2, pch = 20, ps=20)

   dev.copy2eps(file = eps.filename)
   }
		
		
test.subset.factorization <- function(n = 100, A =2:10, num.avg = 50,
                                          m=20, df=0.9){
# In this function each Markov process is factored  indivudually, and the
# archetypes correspond to a subset of the original states

	res <- matrix(0, num.avg, length(A))
	
	for (i in 1:length(A)) {
		na <- A[i]
		
      Pb <- array(0, c(m,m,na))
      
      P <- array(0, c(n,n,na))
      r <- matrix(0, n, na)
      
      D <- array(0, c(n,m,na))
      
      Ka <- array(0, c(m,n,na))
      rb <- matrix(0, m, na)
      
      for (j in 1:num.avg) {
         for (a in 1:na) {
            Ka[,,a] <-  normal.transition.matrix(m, n, m, sd = 1, noise = 0)
            rb[,a] <- runif(m)
            
            D[,,a] <- matrix(runif(n * m), n, m)
            D[,,a] <- D[,,a] / apply(D[,,a], 1, sum)
            
             for (k in 1:m) {
             	D[k,,a] <- 0
             	D[k,k,a] <- 1
             	}
             	
            
            P[,,a] <- D[,,a] %*% Ka[,,a] 
            r[,a] <- D[,,a] %*% rb[,a]
            
            Pb[,,a] <- Ka[,,a] %*% D[,,a]
            
            }
  			       
  			
   
         pi  <- policy.iteration(r, P, df)$pi
         Qb <- policy.iteration(rb, Pb, df)$Q
         
         Q <- matrix(0, n, na)
         for (a in 1:na) {
         	Q[,a] <- D[,,a] %*% Qb[,a]
         	}
         pib <- apply(Q, 1, which.max)
         
         res[j,i] <- sum(pi != pib) / n
         }
     }
  res
  }


test.subset.factorization2 <- function(n = 100, A =2:10, num.avg = 50,
                                          m=20, df=0.9){
# Like "test.subset.factorization", but here the MDP is factored all at once

	res <- matrix(0, num.avg, length(A))
	
	for (i in 1:length(A)) {
		na <- A[i]
		
      Pb <- array(0, c(m,m,na))
      
      P <- array(0, c(n,n,na))
      r <- matrix(0, n, na)
      
      D <- array(0, c(n,m,na))
      
      for (j in 1:num.avg) {

         K <-  normal.transition.matrix(m, n, m, sd = 1, noise = 0)
         rb <- runif(m)
         
         for (a in 1:na) {
      	print(paste(i,j,a))
      	   D[,,a] <- matrix(runif(n * m), n, m)
            D[,,a] <- D[,,a] / apply(D[,,a], 1, sum)
            
             for (k in 1:m) {
             	D[k,,a] <- 0
             	D[k,k,a] <- 1
             	}
             	
            
            P[,,a] <- D[,,a] %*% K 
            r[,a] <- D[,,a] %*% rb
            
            Pb[,,a] <- K %*% D[,,a]
            
            }
  			       
  			
  			rb <- matrix(rep(rb,na),m,na)
        
         pi  <- policy.iteration(r, P, df)$pi
         Qb <- policy.iteration(rb, Pb, df)$Q
         
         Q <- matrix(0, n, na)
         for (a in 1:na) {
         	Q[,a] <- D[,,a] %*% Qb[,a]
         	}
         pib <- apply(Q, 1, which.max)
         
         res[j,i] <- sum(pi != pib) / n
         }
     }
  res
  }
