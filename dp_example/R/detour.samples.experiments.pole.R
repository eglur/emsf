source("pole.R")
source("rbfn.rl.R")
source("detour.samples.R")
source("lspi.rbfn.R")
source("util.R")
source("detour.samples.experiments.R")

pole.data <- function (size, normalize = TRUE, lsars = FALSE, max.steps =
               Inf, sd.theta = pi/4, sd.theta.vel = 1, ...){
   t <- rnorm(size, 0, sd.theta)
   v <- rnorm(size, 0, sd.theta.vel)
   s.initial <- matrix(cbind(t,v), size, 2)
   sars <- collect.episodes(size, s.initial, pole.transition, c(-50,0,50),
               max.steps, normalize, ...) 
   if (lsars) make.lsars(sars)
   else sars
   }
   


pole.double.data <- function (size, normalize = TRUE, lsars = FALSE, max.steps =
               Inf, sd.pos = 1, sd.vel = 1.3, sd.theta1 = pi/13, 
               sd.theta1.vel = 0.5, sd.theta2 = pi / 5.5, 
               sd.theta2.vel = 0.4, A = c(-1,0,1), ...){
# the standard deviations were obtained experimentally; I ran the simulator
# with a random policy for 5 times and took the maximum value of each variable
# (except the car position) 
   x   <- rnorm(size, 0, sd.pos)
   xv  <- rnorm(size, 0, sd.vel)
   t1  <- rnorm(size, 0, sd.theta1)
   t1v <- rnorm(size, 0, sd.theta1.vel)
   t2  <- rnorm(size, 0, sd.theta2)
   t2v <- rnorm(size, 0, sd.theta2.vel)
   
   s.initial <- matrix(cbind(x, xv, t1, t1v, t2, t2v), size, 6)
   sars <- collect.episodes(size, s.initial, pole.double.transition, A,
               max.steps, normalize, ...) 
   if (lsars) make.lsars(sars)
   else sars
   }
   
     
     
pole.test.policy.rbfn <- function(S, rbfn, sample.mean, sample.sd, max.it = 
                           3000, rbfn.out = rbfn.norm.output,
                           pole.function = pole.transition, 
                           A = c(-50,0,50)) {
#    cp <- seq(-1.2, 0.5, length = 100)
#    cv <- seq(-0.07, 0.07, length = 100)
#    S2 <- gp(cp, cv)
#    S2 <- normalize(S2, means = sample.mean, stdevs = sample.sd)
#    o <- rbfn.out(rbfn, S2)
#    p <- apply(o, 1, max)
#    image(-t(matrix(p,100,100)))
   ns <- array(0, nrow(S))
   for (s in 1:nrow(S)) {
      ns[s] <- control.rbfn(rbfn, S[s,], pole.function, A, 
               sample.mean, sample.sd, max.it, rbfn.out)$ns
      }
   ns
   }
   
   
               
## generalizar pra poder usar o double             
pole.lspi <- function(num.avg, sars, S.test, rbfn.size = 3^2, df = 0.95, tau =
                   1e-1, it = 10, max.it.control = 3000, verbose = FALSE,
                   rbfn.out = rbfn.norm.output, cluster = FALSE, ...) {
     res <- matrix(0, num.avg, nrow(S.test))
     for (i in 1:num.avg) {
          if (verbose) print(paste("Run #", i))
          
          C <- NULL
          if (cluster) {
            C <- kmeans(sars$s, rbfn.size)$centers
            }
          else {
            C <- matrix(0, rbfn.size, ncol(S.test))
            max.X <- apply(sars$s, 2, max)
            min.X <- apply(sars$s, 2, min)
            C <- uniform.grid(C, max.X, min.X)
            }
          
#          rbfn <- make.rbfn.kbrl(C, 3, tau)       
          
          rbfn <- lspi.rbfn(sars, rbfn, df, it, verbose = FALSE, ...)
          
          if (verbose) print("Testing the policy...")
          res[i,] <- pole.test.policy.rbfn(S.test, rbfn, sars$means,
                        sars$stdevs, max.it = max.it.control, rbfn.out =
                        rbfn.out)
          if (verbose) print(paste("Average number of steps:", mean(res[i,]))) 
          }
     res
     }
     


pole.detour <- function(num.avg, sars, S.test, rbfn.size = 3^2, df = 0.95, tau
                   = 1e-3, it = 10, max.it.control = 3000, verbose = FALSE,
                   pole.function = pole.transition, A = c(-50,0,50), ...){
   lsars <- make.lsars(sars)
   res <- matrix(0, num.avg, nrow(S.test))
   prop.centers <- rbfn.size / nrow(sars$s) 
   for (i in 1:num.avg) {
      if (verbose) print(paste("Run #", i))
      rbfn <- detour.rbfn(lsars, tau, df, prop.centers, iter.pi = it,
               run.value.iteration = FALSE, ...)
      
      if (verbose) print("Testing the policy...")
      res[i,] <- pole.test.policy.rbfn(S.test, rbfn, lsars$means,
                        lsars$stdevs, max.it = max.it.control, rbfn.out =
                        rbfn.norm.output, pole.function = pole.function, A = A)
      
      if (verbose) print(paste("Average number of steps:", mean(res[i,]))) 
     }
   res
   }
    