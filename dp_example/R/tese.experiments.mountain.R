# Experiments with archetypes in the mountain-car task

source("mountain.car.R")
source("tese.plot.R")
source("mdp.R")
source("dp.R")
source("data.manipulation.R")
source("arch.R")
source("arch.kernel.R")
source("dp.2D.R")
source("kbrl.R")

generate.persp.map <- function(vm = NULL, r=25, n = 250, row = 8, nc=5,
                  file.name="./fig_tese/mountain"){
# 'r' is the resolution   
   cp <- seq(-1.2, 0.5, l=n)
   cv <- seq(-0.07, 0.07, l=n)

   if (is.null(vm)) {
      vm <- value.iteration.2D(cp, cv, c(-1,0,1), mountain.car.transition,
                              0.995)
      vm <- fix.mountain(vm, 0)
      }

  # generate persp
   vmm <- matrix(0, r, r)
   ind <- floor(seq(1,ncol(vm),l=r))
   
   for (i in 1:length(ind)) vmm[i,] <- vm[ind[i], ind]
   
   par(mar=c(0,0,0,0) + 0.1)
   persp(-vmm, xlab="posição", ylab="velocidade", zlab="-V*",
		shade=0.75, theta=120, phi=50, border="BLACK", col="WHITE")
   dev.copy2eps(file=paste(file.name,"_persp.eps",sep=""))

   set.par.tese()
   image(cp, cv, vm, xlab=expression(x), ylab=expression(dot(x)),
         col=grey.colors(50))
   cp <- seq(-1.1, 0.4, l=nc)
   cv <- seq(-0.06, 0.06, l=nc)
   c <- gp(cp,cv)   

   points(matrix(c[row,],1,2), pch=21, col="black", bg="blue", cex=2)
   points(matrix(c(-0.9, -0.01), 1, 2), pch=22, col="black", bg="red", cex=2)
   
   text(locator(1),expression(s[i]))
   text(locator(1),expression(s[j]))
   
   dev.copy2eps(file=paste(file.name,"_map.eps",sep=""))  
   }


generate.map.clusters <- function(vm = NULL, n = 250, nc = 5, 
                                  rows = c(8,2,3,7), 
                                  file.name="./fig_tese/mountain"){

   cp <- seq(-1.2, 0.5, l=n)
   cv <- seq(-0.07, 0.07, l=n)

   if (is.null(vm)) {
      vm <- value.iteration.2D(cp, cv, c(-1,0,1), mountain.car.transition,
                              0.995)
      vm <- fix.mountain(vm, 0)
      }

   set.par.tese()
   image(cp, cv, vm, xlab=expression(x), ylab=expression(dot(x)),
         col=grey.colors(50))
   
   cpp <- seq(-1.1, 0.4, l=nc)
   cvp <- seq(-0.06, 0.06, l=nc)
   c <- gp(cpp,cvp)   
   p <- matrix(c(-0.9, -0.01),1,2)
   
   lines(rbind(p,c[rows[1],]), lty=1,lwd=4)
   points(c, pch=21, col="black", bg="blue",cex=2)
   points(p, pch=22, col="black", bg="red", cex=2)

   dev.copy2eps(file=paste(file.name,"_map_clusters_1.eps",sep=""))

   image(cp, cv, vm, xlab=expression(x), ylab=expression(dot(x)),
         col=grey.colors(50))
   lines(rbind(p,c[rows[1],]), lty=2,lwd=4)
   lines(rbind(p,c[rows[2],]), lty=2,lwd=4)
   lines(rbind(p,c[rows[3],]), lty=2,lwd=4)
   lines(rbind(p,c[rows[4],]), lty=2,lwd=4)
   points(c, pch=21, col="black", bg="blue",cex=2)
   points(p, pch=22, col="black", bg="red", cex=2)

   dev.copy2eps(file=paste(file.name,"_map_clusters_4.eps",sep=""))
   }
   

generate.kernel.exs <- function(x = seq(-3,3,l=100), c = 0, sigma = 1,
                               filename="./fig_tese/kernel_exs.eps") {
   gaussian <- function(d) exp(-1/2 * d^2)/(sqrt(2) *pi)
   linear <- function(d) pmax(1 - abs(d),0)
   epanechnikov <- function(d) pmax(3/4*(1-d^2), 0)
   cosseno <- function(d) pmax(0, pi/4 * cos(pi/2 * d))
   uniforme <- function(d) {r <- array(0, length(d)); r[d <= 1] <- 1/2; r}

# Generate unidimensional RBFs
   d <- sqrt((x - c)^2)/sigma
   r <- cbind(uniforme(d), linear(d), gaussian(d), epanechnikov(d))
   matplot(x, r, lwd = 3, col = col.tese(ncol(r)), t="l", cex=2.5,
           ylab=expression(kappa))
    leg(c("Uniforme", "Linear", "Gaussiana", "Epanechnikov"),
        lwd = 2, col = col.tese(ncol(r)), x=locator(1), inset = 0.02, cex=1.2)
   dev.copy2eps(file = filename)
   }   

generate.kernel.pert <- function(n = 250, 
                               cp = seq(-1, 1,l=n),
                               cv = seq(-1, 1,l=n),
                               x = gp(cp,cv),
                               nc = 5, 
                               filename="./fig_tese/kernel_pert.eps") {
   
   gaussian <- function(x, C=c(0,0), sigma=0.25) exp(-apply(t(t(x)-C)^2,1,sum) /
                                             sigma^2)
   
   cpp <- seq(-0.9, 0.9, l=nc)
   cvp <- seq(-0.9, 0.9, l=nc)
   c <- gp(cpp,cvp)   

   y <- array(0, nrow(x)) 
   for (i in 1:nrow(c)) y <- y + gaussian(x, c[i,])
   y <- matrix(y, n, n)
 
   set.par.tese()

   image(seq(-1.2, 0.5, l=n), seq(-0.07, 0.07, l=n), y,
         xlab=expression(x), ylab=expression(dot(x)),col=grey.colors(50))
   
   cp <- seq(-1.1, 0.4, l=nc)
   cv <- seq(-0.06, 0.06, l=nc)
   c <- gp(cp,cv)   
   points(c, pch=21, col="black", bg="blue", cex=2)
   dev.copy2eps(file = filename)
   }

      

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
                    ind.2d <- calc.ind.2d(t$s, delta.x, delta.y, min.pos,
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

Q.1d.2d <- function(Q) {
# converts a Q-function in 1d to 2d
     l <- sqrt(length(Q[,1]))
     Q2d <- array(0, c(l,l,ncol(Q)))
     for (a in 1:ncol(Q)) {
          Q2d[,,a] <- matrix(Q[,a], l, l, byrow = TRUE)
          }
     Q2d
     }
          
Q.1d.V.2d <- function(Q) {
# converts a Q-function in 1d to a V-function in 2d
     l <- sqrt(length(Q[,1]))
     V <- apply(Q,1,max)
     matrix(V, l, l, byrow = TRUE)
     }

mountain.test.policy.n <- function(policy, S, n, sd = 0) {
   res <- matrix(0, n, nrow(S))
   for (i in 1:n) res[i,] <- mountain.car.test.policy(policy,S, sd=sd,
                                                            max.steps=200)
   res
   }

		
mountain.experiment <- function( n = 30, 
                                 ma = c(400, 500, 600, 700), 
                                 mk = seq(30,270,by=30), 
                                 tau = c(10,1,0.1,0.01),
                                 sd=1,
                                 num.avg = 20,
                                 df = 0.99,
                                 max.iter = 15,
                                 max.iter.ff = 20,
                                 ST =rbind(gp(seq(-1,0.15,l=5),
                                            seq(-0.07,0.02,l=5)), c(0,0)),
                                 names = c("standard", "archetype", "kernel"),
                                 run = c(TRUE, TRUE, TRUE),
                                 dir = "./res_tese/mountain/",
                                 factor.function = nmf.kmeans, 
                                 select.function = select.grid.mountain,
                                 alloc.function = alloc.prop,
                                 na = 101,
                                 verbose = FALSE) {
# 'run' is a flag to (de)activate each method
   
   # solve the problem discretizing the domain
   if (verbose) print("Discretizing the domain...")
   M <- discretize.mountain(n, sd, na=na)
   
   if (run[1]) {
      if (verbose) print("Solving the problem using the entire grid...")
      time <- system.time(
         D <- policy.iteration(M$R, M$P, df = df, max.iter = max.iter)
         , TRUE)[1]
      pi <- policy.1d.2d(D$pi)
      res <- mountain.test.policy.n(pi, ST, num.avg, sd=sd)
      wt(res, paste(dir, names[1],"_n",n,".txt", sep=""))
      wt(time, paste(dir, names[1],"_n",n,"_time.txt", sep=""))
      }

   #now, use the archetypes in its standard form
   
   if (run[2]) {
      ## REMEMBER TO ADD TIME
      if (verbose) 
         print("Solving the problem using standard archetype algorithm...")
      for (m in ma) {
         if (verbose) cat(paste("m= ",m, "( "))
         res <- matrix(0, num.avg^2, nrow(ST))
         for (i in 1:num.avg) {
            if (verbose) cat(paste(i,""))
            time <- system.time(
               D <- arch.policy.iteration(M$R, M$P, df = df, m = m, 
                                       max.iter=max.iter,
                                       max.iter.ff=max.iter.ff,
                                       factor.function = factor.function)
                               , TRUE)[1]
            pi <- policy.1d.2d(D$pi)
         
            b <- (i-1) * num.avg + 1
            e <- b + num.avg -1
            res[b:e,] <- mountain.test.policy.n(pi, ST, num.avg, sd=sd)
            }
         if (verbose) cat(")\n")
         wt(res, paste(dir,names[2],"_m",m,".txt", sep=""))
         }
      }
   

   #Finally, use the kernalized archetype algorithm 
   if (run[3]) {
      if (verbose) 
        print("Solving the problem using the kernalized archetype algorithm...")
      for (m in mk) {
         if (verbose) cat(paste("m= ",m, "( "))
         for (t in tau) {
            cat(paste("tau =", t, " "))
            time <- system.time(            
            D <- arch.kernel.policy.iteration(M$R, M$P, df = df, m = m, S = M$S,
                                    max.iter=max.iter, select = select.function,
                                    alloc = alloc.function, tau=t), TRUE)[1]
            pi <- policy.1d.2d(D$pi)
            res <- mountain.test.policy.n(pi, ST, num.avg, sd=sd)
            wt(res, paste(dir,names[3],"_m",m,"_t", t, ".txt", sep=""))
            wt(time, paste(dir,names[3],"_m",m,"_t", t, "_time.txt", sep=""))
            }
         if (verbose) cat(")\n")
         }
      }

   if (verbose) print("All done.")
   }


plot.results.standard <- function(vm,
                          ST = gp(seq(-1,0.15,l=5),seq(-0.07,0.02,l=5)),
                          filename ="./res_tese/mountain/standard_n30.txt",
                          eps.filename ="./fig_tese/mountain_map_results.eps"){
   R <- mean(read.table(filename))
   R <- R[-length(R)]
   n <- nrow(vm)
   cp <- seq(-1.2, 0.5, l=n)
   cv <- seq(-0.07, 0.07, l=n)
   set.par.tese()
   image(cp, cv, vm,xlab=expression(x), ylab=expression(dot(x)),
         col=grey.colors(50))
   points(ST, pch=23, col="black", bg="green3", cex=2)
   c <- ST + matrix(rep(c(0,0.007), nrow(ST)),nrow(ST),ncol(ST), byrow=TRUE)
   text(c, format(R,nsmall=2))
   dev.copy2eps(file=eps.filename)
   }
   
 


generate.mp.sa <- function(n = 30, 
                           ma = c(400, 500, 600, 700), 
                           names = c("standard", "archetype"),
                           dir = "./res_tese/mountain/",
                           eps.filename = "./fig_tese/mountain_res_sa.eps"
                           ) {
# generates the matplots for standard and archetypes stochastic approximation
   R <- NULL
   for (m in ma) {
      R <-cbind(R, mean(read.table(paste(dir,names[2],"_m",m,".txt",sep="")))) 
      }
   R <- cbind(R, mean(read.table(paste(dir, names[1],"_n",n,".txt", sep=""))))
   R <- R[-nrow(R),] 
   leg <- c(paste("PISF-M(900,",ma,")",sep=""), "PI(900)")
   # make "black" the last color and "circle" the last pch
   col <- col.tese(ncol(R))
   col <- c(col[2:(length(col))],col[1])
   pch <- 21:(20+ncol(R))
   pch <- c(pch[2:(length(pch))],pch[1])
 
   mp(R, xlab="Estado do conjunto de teste", ylab="Passos", t="o", bty="o",
      leg.pos="bottomright", leg=leg, ps = 13, col=col, pch=pch)
   dev.copy2eps(file = eps.filename)
   }
   


generate.mp.sk <- function(n = 30, 
                           mk = c(60, 90,120), 
                           names = c("standard", "kernel"),
                           dir = "./res_tese/mountain/",
                           eps.filename = "./fig_tese/mountain_res_sk.eps"
                           ) {
# generates the matplots for standard and kernalized stochastic approximation
   R <- NULL
   for (m in mk) {
      R <-cbind(R, 
         mean(read.table(paste(dir,names[2],"_m",m,"_t0.1.txt",sep="")))) 
      }
   R <- cbind(R, mean(read.table(paste(dir, names[1],"_n",n,".txt", sep=""))))
   R <- R[-nrow(R),] 
   # make "black" the last color and "circle" the last pch
   col <- col.tese(ncol(R))
   col <- c(col[2:(length(col))],col[1])
   pch <- 21:(20+ncol(R))
   pch <- c(pch[2:(length(pch))],pch[1])
  
   leg <- c(paste("PISF-S(900,",mk,")", sep=""), "PI(900)")
   mp(R, xlab="Estado do conjunto de teste", ylab="Passos", t="o", bty="o",
      leg.pos="topleft", leg=leg, ps = 15, inset=0.01, col=col, pch=pch)
   dev.copy2eps(file = eps.filename)
   }


generate.mp.sk2 <- function(n = 30, 
                           mk = seq(30,270,by=30), 
                           names = c("standard", "kernel"),
                           dir = "./res_tese/mountain/",
                           eps.filename = "./fig_tese/mountain_res_sk2.eps"
                           ) {
# plots results of PI and PISA together
   R <- matrix(0,length(mk), 2)
   for (m in 1:length(mk)) {
      R[m,1] <- mean(
         mean(read.table(paste(dir,names[2],"_m",mk[m],"_t0.1.txt",sep="")))
         ) 
      }
   R[,2]<-mean(mean(read.table(paste(dir,names[1],"_n30.txt",sep="")))) 
   # make "black" the last color and "circle" the last pch
   col <- col.tese(ncol(R))
   col <- c(col[2:(length(col))],col[1])
   pch <- 21:(20+ncol(R))
   pch <- c(pch[2:(length(pch))],pch[1])
   
   leg <- c("PISF-S(900,.)","PI(900)")
   mp(x = mk, R, xlab="m", ylab="Passos", t="o", bty="o",
      leg.pos="topright", leg=leg, ps = 15, inset=0.05, col=col, pch=pch)
   dev.copy2eps(file = eps.filename)
   }


print.info <- function(D) {
   D <- as.matrix(D)
   D1 <- matrix(D, nrow(D) * ncol(D), 1)
   Z <- D[,ncol(D)]
   D <- D[, -ncol(D)]
   print(paste(round(mean(D),d=2), round(max(D),d=2), 
               round(min(D),d=2), round(sd(D1),d=2),
               round(mean(Z),d=2), round(max(Z),d=2), 
               round(min(Z),d=2), round(sd(Z),d=2),
               sep = " & "))
   }

generate.table.res <- function( n = 30, 
                                 ma = c(400, 500, 600, 700), 
                                 mk = seq(30,270,by=30), 
                                 tau = c(0.1),
                                 names = c("standard", "archetype", "kernel"),
                                 dir = "./res_tese/mountain/",
                                 run = c(TRUE, TRUE, TRUE)
                                 ) {
   if (run[1]) {
      D <- read.table(paste(dir, names[1],"_n",n,".txt", sep=""))   
      print.info(D)
      }

   if (run[2]) {
      for (m in ma) {
            D <- read.table(paste(dir,names[2],"_m",m,".txt", sep=""))
            print.info(D)
            }
      }

   if (run[3]) {
      for (m in mk) {
         for (t in tau) {
            D <- read.table(paste(dir,names[3],"_m",m,"_t", t, ".txt",sep=""))
            print.info(D)
            }
         }
      }
   }



mountain.data <- function (size, normalize = TRUE, p = 0, random = TRUE,
                  lsars = TRUE, sd = 0, ...){
   # make sure there is at least one goal in the dataset
   num.goals <- max(round(p*size), 1)
   # The state s = (0.5,0.01) will reach the goal even if a = -1
   s <- gp(0.5, seq(0.07, 0.01, length = num.goals))  
   remaining <- size - num.goals
		
   if (random) {
      # generate the dataset with size-num.goals points
      cp <- runif(remaining, -1.2, 0.5)
      cv <- runif(remaining,-0.07, 0.07)
      s <- rbind(s, cbind(cp,cv))
      }
   else {
      # uniform grid
      grid.side <- floor(sqrt(remaining))
      cp <- seq(-1.2, 0.5, length = grid.side)
      cv <- seq(-0.07, 0.07,length = grid.side)
      s <- rbind(s, gp(cp,cv))
      # complete the dataset with random points
      remaining <- size - nrow(s)
      if (remaining > 0) {
         cp <- runif(remaining, -1.2, 0.5)
         cv <- runif(remaining,-0.07, 0.07)
         s <- rbind(s, cbind(cp,cv))
         }
       }
	
      sa <- make.sa(s, 3)
      sars <- collect.transitions(sa, mountain.car.transition, c(-1,0,1),
               normalize = normalize, sd = sd, ...)
      if (lsars) make.lsars(sars)
      else sars
      }



mountain.kbrl.experiment <- function(
                              num.transitions = seq(100,500,by=100),
                              num.archs = c(0.5,0.3,0.1),
                              perc.archs = TRUE,
                              tau.k = 0.3, tau.q = 0.01, tau.a = 0.3,
                              random = TRUE,
                              ST=rbind(gp(seq(-1,0.15,l=5),
                                       seq(-0.07,0.02,l=5)), c(0,0)),
                              num.avg.test = 20,   
                              num.avg.runs = 20,
                              sd = 1,
                              gamma = 0.995, epsilon = 0.01, iter.max = Inf,
                              names = c("kbrl","kbrl_arch"),
                              dir = "./res_tese/mountain/",
                              run = c(TRUE, TRUE),
                              verbose = FALSE) {
# the variable `perc.archs` indicates whether `num.archs` should be interpreted
# as proportions of `num.transitions`
num.archs.original <- num.archs
for (n in num.transitions) {
         if (verbose) print(paste("Number of transitions", n))
         res.kbrl <- matrix(0, num.avg.runs * num.avg.test, nrow(ST))
         res.arch <- array(0, c(length(num.archs), num.avg.runs *
                           num.avg.test, nrow(ST)))
         if (perc.archs) num.archs <- round(num.archs.original * n)
         for (i in 1:num.avg.runs) {
            if (verbose) print(paste("Run",i, "of", num.avg.runs))
            lsars <- mountain.data(n, normalize = TRUE, p = 0, random =
                                 random, lsars = TRUE, sd = sd)
             b <- (i-1)*num.avg.test + 1
             e <- b + num.avg.test - 1
     
             if (run[1]) {
               # run KBRL
               rbfn <- kbrl(lsars, tau.k, gamma, epsilon, iter.max)
               res.kbrl[b:e,] <- mountain.car.test.policy.rbfn(rbfn, S = ST,
                                    num.avg =num.avg.test, sd =sd, means =
                                    lsars$means,stdevs=lsars$stdevs)       
               
               if (verbose) print(paste("KBRL n =",n, "t =", tau.k ,
                             "num. steps:",round(mean(res.kbrl[b:e,]),d=2)))
                }

            if (run[2]) {
               # run KBRL-A
               for (a in 1:length(num.archs)) {
                  rbfn <- kbrl.arch(lsars, tau.q= tau.q, tau.a = tau.a,
                     num.archs = num.archs[a], gamma =gamma, epsilon.vi =
                     epsilon,iter.vi =iter.max)
                  res.arch[a,b:e,] <- mountain.car.test.policy.rbfn(rbfn,S=ST,
                                       num.avg =num.avg.test, sd =sd, means =
                                       lsars$means,stdevs=lsars$stdevs)
                  if (verbose) print(paste("KBRL-A n =",n, "tq =", tau.q ,
                                       "ta =",tau.a, "a =",num.archs[a], 
                              "num.steps:", round(mean(res.arch[a,b:e,]),d=2)))
                  }
               }      
            }
         
         if (run[1]) wt(res.kbrl,paste(dir,names[1],"_n",n,"_t",tau.k,".txt",
                                       sep=""))
         if (run[2]) {
            for (a in 1:length(num.archs)) wt(res.arch[a,,], paste(dir,
                                          names[2],"_n",n,"_a", num.archs[a],
                                      "_tq",tau.q,"_ta", tau.a,".txt", sep=""))
            }
         }
   }


generate.table.kbrl<- function( num.transitions = seq(100,500,by=100), 
                                perc.transitions = c(0.5, 0.3, 0.1), 
                                tau.k = 0.3, tau.q = 0.01, tau.a = 0.3,
                                dir = "./res_tese/mountain/",
                                names = c("kbrl", "kbrl_arch")
                                 ) {
   for (n in num.transitions) {
      D <- read.table(paste(dir, names[1],"_n", n, "_t", tau.k, ".txt", sep=""))
      print.info(D)
      for (p in perc.transitions) {
         a <- p * n
         D <- read.table(paste(dir,names[2],"_n",n,"_a",a,"_tq",tau.q,"_ta",
                               tau.a,".txt",sep=""))
         print.info(D)
         }
      }
   }


generate.mp.kbrl <- function(
                           num.transitions = c(100,300,500),
                           t = 0.3, 
                           dir = "./res_tese/mountain/",
                           eps.filename = "./fig_tese/mountain_res_kbrl.eps"
                           ) {
# generates the matplots for standard KBRL
   R <- NULL
   for (n in num.transitions) {
      R <-cbind(R,mean(read.table(paste(dir,"kbrl_n",n,"_t",t,".txt",
                       sep="")))) 
      }
   R <- R[-nrow(R),] 
   leg <- paste("KBRL(",num.transitions,")",sep="")
   # make "black" the last color and "circle" the last pch
   
   mp(R, xlab="Estado do conjunto de teste", ylab="Passos", t="o", bty="o",
      leg.pos="topleft", leg=leg, ps = 13, inset=0.01)
   dev.copy2eps(file = eps.filename)
   }
   

generate.mp.kbrl.sa <- function(
                           num.transitions = 300,
                           num.archs = c(0.5, 0.1),
                           arch.perc = TRUE,
                           t = 0.3,
                           tau.q = 0.01,
                           tau.a = 0.3, 
                           dir = "./res_tese/mountain/",
                           eps.filename = "./fig_tese/mountain_res_kbrl_sa.eps"
                           ) {
# generates the matplots for standard KBRL
   if (arch.perc) num.archs <- num.archs * num.transitions

   R <- mean(read.table(paste(dir,"kbrl_n",
             num.transitions,"_t",t,".txt",sep="")))

   for (a in num.archs) {
      R <-cbind(R,mean(read.table(paste(dir,"kbrl_arch_n",num.transitions,
          "_a", a, "_tq",tau.q, "_ta", tau.a,".txt",sep="")))) 
      }
   R <- R[-nrow(R),] 
   leg <- c(paste("KBRL(",num.transitions,")", sep=""),  
            paste("KBRL-SA(",num.transitions,",",num.archs,")",sep=""))
   # make "black" the last color and "circle" the last pch
   col <- col.tese(ncol(R))
   col <- c(col[2:(length(col))],col[1])
   pch <- 21:(20+ncol(R))
   pch <- c(pch[2:(length(pch))],pch[1])
 
   mp(R, xlab="Estado do conjunto de teste", ylab="Passos", t="o", bty="o",
      leg.pos="topleft", leg=leg, ps = 13, inset=0.01)
   dev.copy2eps(file = eps.filename)
   }


plot.kbrl.inc.sample <- function(num.transitions = seq(100,500,by=50),
                                 num.archs = 100, tau.q = 0.01, tau.a = 0.3,
                                dir = "./res_tese/mountain/",                  
                                eps.filename = "./fig_tese/kbrl_inc_samp.eps"){
   d <- array(0, length(num.transitions))
   for (n in 1:length(d)) {
      d[n] <- mean(mean(read.table(paste(dir, "kbrl_arch_n", num.transitions[n],
                                 "_a", num.archs, "_tq", tau.q, "_ta",
                               tau.a,".txt",sep=""))))
      }
   
   mp(x = num.transitions, matrix(d), xlab="Número de transições", 
         ylab="Passos", t="o", bty="o", ps = 15)
   dev.copy2eps(file = eps.filename)
  }


plot.kbrl.inc.sample.english <- function(num.transitions = seq(100,500,by=50),
                        num.archs = 100, tau.q = 0.01, tau.a = 0.3,
                        dir = "./res_tese/mountain/", 
                        eps.filename ="~/tex/papers/sf_rl/kbsf_mountain.eps"){
   d <- array(0, length(num.transitions))
   for (n in 1:length(d)) {
      d[n] <- mean(mean(read.table(paste(dir, "kbrl_arch_n", num.transitions[n],
                                 "_a", num.archs, "_tq", tau.q, "_ta",
                               tau.a,".txt",sep=""))))
      }
   
   mp(x = num.transitions, matrix(d), xlab="Number of transitions n", 
         ylab="Number of steps to escape from the valley", t="o", bty="o", ps =
         15)
   lines(num.transitions, rep(107.67, length(num.transitions)), 
         col="RED", lty=2)
   text(num.transitions[8], 106.3, "KBRL(100)")
   
   lines(num.transitions, rep(76.46, length(num.transitions)), 
         col="BLUE", lty=2)
   text(num.transitions[8], 77.5, "KBRL(500)")

   text(locator(1), "KBSF(n, 100)")

dev.copy2eps(file = eps.filename)
  }


mountain.kbrl.pisa <- function(
                              num.transitions = 300,
                              num.archs = c(0.5,0.3,0.1),
                              perc.archs = TRUE,
                              taus = c(1,0.5,0.1), 
                              tau.q = 0.01, 
                              tau.a = 0.3, 
                              ST=rbind(gp(seq(-1,0.15,l=5),
                                       seq(-0.07,0.02,l=5)), c(0,0)),
                              num.avg.test = 20,   
                              num.avg.runs = 20,
                              max.iter.pi = 15,
                              select.function = select.clara,
                              alloc.function = alloc.prop,
                              arch.function = arch.clara, 
                              sd = 1,
                              gamma = 0.995, 
                              dir = "./res_tese/mountain/",
                              verbose = FALSE) {
# the variable `perc.archs` indicates whether `num.archs` should be interpreted
# as proportions of `num.transitions`
   if (perc.archs) num.archs <- num.archs * num.transitions

   for (m in num.archs) {
      res <- array(0, c(num.avg.runs * num.avg.test, nrow(ST),length(taus)))
      res.pisa <- matrix(0, num.avg.runs * num.avg.test, nrow(ST))

      for (i in 1:num.avg.runs) {
         if (verbose) print(paste("Creating MDP",i))
         lsars <- mountain.data(num.transitions, normalize = TRUE, p = 0,
                                     random=TRUE, lsars = TRUE, sd = sd)
         M <- kbrl.mdp(lsars, tau.a, goal.reward = 1)

         b <- (i-1)*num.avg.test + 1
         e <- b + num.avg.test - 1
         
         for (t in 1:length(taus)) {
            A <- arch.kernel.policy.iteration(M$R, M$P, gamma, m, M$S, 
                                          max.iter=max.iter.pi, select =
                                          select.function, alloc =
                                          alloc.function, tau = taus[t])
            rbfn <- make.rbfn.kbrl(M$S, tau.a, lsars$num.actions)
            rbfn$w <- A$Q
         
            res[b:e,,t] <- mountain.car.test.policy.rbfn(rbfn, S = ST,
                                    num.avg =num.avg.test, sd =sd, means =
                                    lsars$means, stdevs=lsars$stdevs)       
           if (verbose) print(paste("PISA-S m =", m, "t =",taus[t],":",
                           mean(mean(res[b:e,,t]))))
               
            }

         rbfn <- kbrl.arch(lsars, tau.q = tau.q, tau.a = tau.a,
                     num.archs = m, gamma =gamma, iter.pi = max.iter.pi, 
                     run.value.iteration = FALSE, arch.function = arch.function)
         res.pisa[b:e,] <- mountain.car.test.policy.rbfn(rbfn, S = ST,
                     num.avg =num.avg.test, sd =sd, means =
                     lsars$means, stdevs=lsars$stdevs)
         if (verbose) print(paste("KBRL-SA m =", m,":",
                                    mean(mean(res.pisa[b:e,]))))          
         }

      for (t in 1:length(taus)) wt(res[,,t], paste(dir,"kernel_kbrl_n",
                                 num.transitions, "_a",m,"_t",taus[t],".txt",
                                 sep=""))
      wt(res.pisa,paste(dir,"kbrl_arch_pi_n",num.transitions,"_a", m,
              "_tq",tau.q,"_ta", tau.a,".txt", sep=""))
      }

   }

print("tese.experiments.mountain.R loaded")                 