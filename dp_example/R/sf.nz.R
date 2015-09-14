source("count.nonzeros.R")
source("util.R")

nz <- function(A) sum(A!=0)

respect.restrictions <- function(A) {
     respect <- TRUE
     sum.cols <- apply(A,2,sum)
     if (sum(sum.cols == 0) > 0) respect <- FALSE
     if (respect) {
          sum.rows <- apply(A,1,sum)
          if (sum(sum.rows == 0) > 0) respect <- FALSE
          }
     respect
     }


prop.nonzeros <- function(n, filename = "./sf/pop_non_zeros_", verbose=FALSE) {
     filename <- paste(filename, n, ".txt", sep="")

     stats <- matrix(0,6,n-1)

     for (m in 1:(n-1)) {
          print(m)

          valids <- 0
          nz.decreased <- 0
          nz.increased <- 0
          pp <- NULL
          pm <- NULL

          nm <- n*m
          D <- array(0, nm)
          K <- array(0, nm)
          
          pm  <- 0
          pbm <- 0
          D.repeated <- FALSE
          while (!D.repeated) {

               D <- inc.odometer(D)
               if (sum(D)==0) D.repeated <- TRUE
               else {
                    dim(D) <- c(n,m)
                    while(!respect.restrictions(D) && !D.repeated) {
                         dim(D) <- nm
                         D <- inc.odometer(D)
                         if (sum(D)==0) D.repeated <- TRUE
                         dim(D) <- c(n,m)
                         }
                    }

               if (!D.repeated) {
                    K.repeated <- FALSE
                    while (!K.repeated) {
                         
                         K <- inc.odometer(K)
                         if (sum(K)==0) K.repeated <- TRUE
                         else {
                              dim(K) <- c(m,n)
                              while(!respect.restrictions(K) && !K.repeated) {
                                   dim(K) <- nm
                                   K <- inc.odometer(K)
                                   if (sum(K)==0) K.repeated <- TRUE
                                   dim(K) <- c(m,n)
                                   }
                              }

                         if (!K.repeated) {
                              valids <- valids + 1
                             if (verbose && !(valids%%10^3)) {
                                   print(paste(m,valids))
                                   }

                              p <- nz(K%*%D) / nz(D%*%K)

                              if (p < 1) {
                                   nz.decreased <- nz.decreased + 1
                                   pp <- c(pp, p)
                                   }
                              else if (p > 1) {
                                   nz.increased <- nz.increased + 1
                                   pm <- c(pm, p)

                                   P <- D %*% K
                                   normal <- TRUE
                                   for (j in 1:nrow(D)) {
                                        for (k in j:nrow(D)) {
                                             if (identical(D[j,], D[k,])) {
                                                normal <- FALSE
                                                }
                                             }
                                        }
                                   if (normal) {
                                        print(D)
                                        print(K)
                                        print(P)
                                        print("Normal matrix!!!")
                                        print(paste("nz(P):", nz(P),"nz(Pb):",
                                           nz(Pb), "sum(P):",sum(P)))
                                        readline()
                                        }
                               
                                   }
#                              if (!valids %% 1000) print(nz.decreased / valids)
                              }
                         }
                    }
               }
          stats[1,m] <- valids
          stats[2,m] <- nz.decreased
          stats[3,m] <- nz.increased
          stats[4,m] <- nz.decreased / valids
          stats[5,m] <- mean(pp)
          stats[6,m] <- mean(pm)
          }
     wt(stats, filename)
     stats
     }



estimate.prop.nz <- function(n, m, gamma, p, verbose = FALSE) {
# This function estimates the proportion of stochastic factorizations with
# dimensions n x m that generates a Pb with more nonzero elements than P
# The function uses the Chernoff bound to guarantee that the estimate is within
# gamma of the true value with probability p
     increased.nz <- 0

     t <- ceiling(number.samples(gamma,p))
     if (verbose) print(paste("Generating",t,"samples"))
     it <- 0
     not.valid <- 0
     m.squared <- m^2
     while (it < t) {
          D <- matrix(0,n,m)
          K <- matrix(0,m,n)
          
          respect <- FALSE
          while (!respect) {
               D[,] <- sample(0:1, n*m, TRUE)
               respect <- respect.restrictions(D)
               }
          
          respect <- FALSE
          while (!respect) {
               K[,] <- sample(0:1, n*m, TRUE)
               respect <- respect.restrictions(K)
               }

          nzp <- nz(D%*%K)
  
          if (nzp <= m.squared) {
               if (nz(K%*%D) > nzp) increased.nz <- increased.nz + 1
               it <- it + 1
               }
          else not.valid <- not.valid + 1
          }
               
     c(t/(t + not.valid), increased.nz / t)
     }
          

script.estimate.prop.nz <- function(ns = 4:10, gamma = 0.005, p = 0.99) {
     
     res <- matrix(0, ns[length(ns)], ns[length(ns)]-2)
     for (i in ns) {
          for (j in 2:(i-1)) {
               res[i,j-1] <- estimate.prop.nz(i,j,gamma,p)
               print(paste(i,j,res[i,j-1]))
               }
         }
     res
     }

     
number.samples <- function(gamma, p) -1/(2 * gamma^2)*(log(0.5*(1-p)))




test.theory <- function(n, m, t = 1, cols.d=n, cols.k=m) {
     ce <- 0

     for (i in 1:t) {
          D <- matrix(0,n,m)
          K <- matrix(0,m,n)
     
          for (j in 1:nrow(D)) {
                D[j,sample(1:m, cols.d)] <- 1
               }
               
          for (j in 1:ncol(D)) D[sample(1:n,1),j] <- 1

          for (j in 1:nrow(K)) {
                K[j,sample(1:n, cols.k)] <- 1
               }

          for (j in 1:ncol(K)) K[sample(1:m,1),j] <- 1

          P <- D %*% K
          Pb <- K %*% D

          if (nz(Pb) > nz(P)) {
               ce <- ce + 1

#                 if (nz(Pb) > sum(P)) {
#                      print(D)
#                      print(K)
#                      print(P)
#                      print(Pb)
#                      print("Sum is bigger!!!")
#                     print(paste("nz(P):", nz(P), "nz(Pb):",nz(Pb), "sum(P):",
#                          sum(P)))
#                      readline()
#                      }
               
               normal <- TRUE
               for (j in 1:nrow(D)) {
                    for (k in j:nrow(D)) {
                         if (identical(D[j,], D[k,])) normal <- FALSE
                         }
                    }
                    
               if (normal) {
                    print(D)
                    print(K)
                    print(P)
                    print(Pb)
                    print("Normal matrix!!!")
                    print(paste("nz(P):", nz(P), "nz(Pb):",nz(Pb), "sum(P):",
                         sum(P)))
                    readline()
                    }
               }
          }
     ce / t
     }
          

estimate.frac.valids <- function(n,m)  choose(2*m*n, m^2 + m) / 2^(2*n*m)

generate.estimate.frac.valids <- function(n) {
     R <- matrix(0,n,n)
     for (i in 2:n) {
          for (j in 1:(i-1)) {
               R[i,j] <- estimate.frac.valids(i,j)
               }
          }
     R
     }


     
z <- function(n) {
     delta <- sqrt(1 + 8*n)
     d <- 2*n
     x1 <- (-1 - delta) / d
     x2 <- (-1 + delta) / d
     x2

     }


l <- function(alpha) (2 - alpha) / alpha^2




print("sf.nz.R")
