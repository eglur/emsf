source("dp.R")
source("pisf.R")
source("mountain.car.R")

lim.prob.normal <- function(a, sd, tau) {
# Returns the bounds beyond which the normal probability defined around 'a' with
# standard deviation 'sd' is less than 'tau'
      c <- a^2 + 2 * sd^2 * log(sqrt(2*pi) * sd * tau) 
      d <- sqrt(a^2 - c)
      c(a - d, a + d)
      }

mountain.transition.normal <- function(s, a, sd, na, tau = 1e-2, ...) {
# computes 'na' transitions from 's' using a normal distribution with mean 'a'
# and standard deviation 'sd'
# 'na' should be odd

     lim <- c(a,a)
     if (sd > 0) lim <- lim.prob.normal(a, sd,  tau)

     aa <- seq(lim[1], lim[2],l=na)

     s2  <- matrix(0, na, 2)
     r   <- array(0, na)
     for (i in 1:length(r)) {
          t <- mountain.car.transition(s, aa[i], ...)
          s2[i,] <- t$s
          r[i] <- t$r
          }
     list(s2 = s2, r = r)
     }


calc.ind.2d <- function(s, delta.x, delta.y, min.pos = -1.2, min.vel = -0.07) {
# computes the indices as if the state space were 2D
     ind1 <- round((s[,1] - min.pos) / delta.x) + 1
     ind2 <- round((s[,2] - min.vel) / delta.y) + 1
     cbind(ind1, ind2)
     }


conv.2d.1d <- function(row, col, n) {
# converts the 2D indices to 1D indices
     (row - 1) * n + col
     }

prob.int.normal <- function(inf, sup, mean = 0, sd = 1) {
## CONFERIR
# gives the probability p(inf <= x <= sup)
     if (inf >= sup) 0
     else pnorm(sup, low=TRUE) - pnorm(inf, low=TRUE)     
     }


discretize.mountain <- function(n = 30, sd = 1, na = 101, tau = 1e-2, 
                                min.vel = -0.07, max.vel = 0.07, min.pos = -1.2,
                                max.pos =  0.5, regular.reward = 0, 
                                goal.reward = 1, A = c(-1,0,1),
                                verbose = FALSE) {
# creates a discrete MDP with the mountain-car dynamics over an nxn grid
# 'na' should be odd

     if (verbose) print("Allocating P...")
     P <- array(0, c(n^2, n^2, 3))
     if (verbose) print("Allocating R...")
     R <- matrix(0, n^2, 3)

     # computes the probability of the interval associated with each action
     probs <- array(1/na, na)
     
     if (sd > 0) {
          lim <- lim.prob.normal(0, 1, tau)
          x <- seq(lim[1], lim[2], l=na)
     
          b <- -Inf
          for (i in 1:length(probs)) { 
               if (i < length(probs)) e <- x[i] + 0.5 * (x[i+1] - x[i])
               else e <- Inf
               probs[i] <- prob.int.normal(b,e,0,1)
               b <- e
               }
          }

    
     if (verbose) print("Generating S...")
     cp <- seq(min.pos, max.pos, l=n)
     cv <- seq(min.vel, max.vel, l=n)
     delta.x <- cp[2] - cp[1]
     delta.y <- cv[2] - cv[1]
     if (verbose) print("Computing P and R. This may take a while...")
     for (i in 1:length(cp)) {
          for (j in 1:length(cv)) {
               s <- c(cp[i], cv[j])
               ind <- conv.2d.1d(i,j,n)
               for (a in 1:3) {
                    t <- mountain.transition.normal(s, A[a], sd, na, tau,
                           min.pos = min.pos, max.pos = max.pos, 
                           min.vel = min.vel, max.vel = max.vel,
                           regular.reward = regular.reward, 
                           goal.reward = goal.reward)
                    
                    R[ind,a] <- sum(probs * t$r)
                    ind.2d <- calc.ind.2d(t$s2, delta.x, delta.y, min.pos,
                                          min.vel)
                    indices <- conv.2d.1d(ind.2d[,1], ind.2d[,2],n)
                    for (k in 1:length(indices)) {
                         P[ind, indices[k],a] <- P[ind, indices[k],a] + probs[k]
                         }
                    }
               
               }
          }
     if (verbose) print("All done!")
     list(P = P, R = R, S = gp(cp, cv))
     }                   

policy.1d.2d <- function(pi) {
# converts a policy in 1d to 2d
     l <- sqrt(length(pi))
     matrix(pi, l, l, byrow = TRUE)
     }

mountain.test.policy.n <- function(policy, S, n, sd = 0, max.steps) {
   total <- 0
   for (i in 1:n) total <- total + mean(mountain.car.test.policy(policy, S, sd=sd, max.steps=max.steps))
   total / n
   }


mountain.experiment <- function( sqrt.n = 30, 
                                 ms = seq(30,270,by=30), 
                                 taus = c(100,10,1,0.1,0.01),
                                 sd=10,
                                 num.avg = 50,
                                 num.avg.test = 1,
                                 max.steps = 300,
                                 df = 0.99,
                                 max.iter = 15,
                                 max.iter.ff = 5,
                                 ST =rbind(gp(seq(-1,0.15,l=5),
                                            seq(-0.07,0.02,l=5)), c(0,0)),
                                 names = c("standard", "archetype", "kernel"),
                                 run = c(TRUE, TRUE, TRUE),
                                 dir = "./pisf/mountain/",
                                 factor.mdp.function = factor.mdp, 
                                 kb.factor.mdp.function = kb.factor.mdp, 
                                 na = 101,
                                 verbose = FALSE) {
# 'run' is a flag to (de)activate each method
   
   # solve the problem discretizing the domain
   if (verbose) print("Discretizing the domain...")
   M <- discretize.mountain(sqrt.n, sd, na=na, verbose = verbose)
   n <- sqrt.n^2
   num.actions <- ncol(M$R)
   
   if (run[1]) {
      if (verbose) print("Solving the problem using the entire grid...")
      res   <- matrix(0, num.avg, 1)
      times <- matrix(0, num.avg,1)
      for (i in 1:num.avg) {
         times[i,1] <- system.time(
            D <- policy.iteration(M$R, M$P, df = df, max.iter = max.iter)
            , TRUE)[1]
         pi <- policy.1d.2d(D$pi)
         res[i,1] <- mountain.test.policy.n(pi, ST, num.avg.test, sd, max.steps)  
         if (verbose) print(paste("Steps:", res[i,1], " Time", times[i,1]))      
         }
      wt(res, paste(dir, names[1],"_n",n,".txt", sep=""))
      wt(times, paste(dir, names[1],"_n",n,"_times.txt", sep=""))
      }

   #now, use the archetypes in its standard form
   
   if (run[2]) {
      if (verbose)  print("Solving the problem using standard archetype algorithm...")
      for (m in ms) {
         if (verbose) print(paste("m= ",m))
         res   <- matrix(0, num.avg, length(taus))
         times <- matrix(0, num.avg, length(taus))
         for (i in 1:num.avg) {
            for (j in 1:length(taus)) {
               if (verbose) cat(paste("Run",i," Tau", taus[j]))
               times[i,j] <- system.time(
                  {
                  SF <- factor.mdp.function(M$R, M$P, m, taus[j], max.iter = max.iter.ff)
                  D <- pisf(SF, df, max.iter = max.iter)
                  }
                , TRUE)[1]
               pi <- policy.1d.2d(D$pi)
               res[i,j] <- mountain.test.policy.n(pi, ST, num.avg.test, sd, max.steps)
               if (verbose) cat(paste("  Steps:", res[i,j], " Time:", times[i,j], "\n"))
            }
         }
         wt(res, paste(dir,names[2],"_n", n, "_m",m,".txt", sep=""))
         wt(times, paste(dir,names[2],"_n", n, "_m",m,"_times.txt", sep=""))
         }
      }
   

   #Finally, use the kernalized archetype algorithm 
   if (run[3]) {
      if (verbose)  print("Solving the problem using kernalized archetype algorithm...")
      for (m in ms) {
         if (verbose) print(paste("m= ",m))
         res   <- matrix(0, num.avg, length(taus))
         times <- matrix(0, num.avg, length(taus))
         for (i in 1:num.avg) {
            for (j in 1:length(taus)) {
               if (verbose) cat(paste("Run",i," Tau", taus[j]))
               times[i,j] <- system.time(
                  {
                  SF <- kb.factor.mdp.function(M$R, M$P, M$S, m %/% num.actions, taus[j], max.iter = max.iter.ff)
                  D <- pisf(SF, df, max.iter = max.iter)                  }
                , TRUE)[1]
               pi <- policy.1d.2d(D$pi)
               res[i,j] <- mountain.test.policy.n(pi, ST, num.avg.test, sd, max.steps)
               if (verbose) cat(paste("  Steps:", res[i,j], " Time:", times[i,j], "\n"))         }
         wt(res, paste(dir,names[3],"_n", n, "_m",m,".txt", sep=""))
         wt(times, paste(dir,names[3],"_n", n, "_m",m,"_times.txt", sep=""))
         }
      }

   if (verbose) print("All done.")
   }
}   

print("pisf.experiments.R loaded")