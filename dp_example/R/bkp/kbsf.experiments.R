library(cluster)
source("dp.R")
source("data.manipulation.R")
source("rbfn.rl.R")
source("kbrl.R")
source("lspi.rbfn.optimized.R")
source("kbsf.R")
source("puddle.world.R")
source("mountain.car.R")
source("pole.R")
source("pole.double.R")
source("pole.double.center.R")

##----------------- GENERAL ----------------------------------------------

generate.data.set <- function(data.function, sizes, filename,             
                                      num.avg = 30, start = 1, ...) {
   for (s in sizes) {
      for (i in start:num.avg) {
         dt <- data.function(s, ...)
         save.sars(dt, paste(filename,"_data_", s,"_",i, sep=""))
         }
      }
   }


##----------------- KBRL ----------------------------------------------

define.archetypes.kbrl <-function(lsars, num.archs, 
                                  pick.reward = NULL,
                                  include.reward=TRUE, 
                                  perc = 0.5){
# If include==TRUE, the transitions with reward==pick.reward will be picked
# otherwise, the transitions with reward != pick.reward will be picked
# perc is the maximum percentage of transitions picked
     nl <- NULL
     ms <- array(0, lsars$num.actions)
     ms[] <- num.archs %/% lsars$num.actions 
     ind <- sample(1:length(ms),num.archs %% lsars$num.actions)
     ms[ind] <- ms[ind] + 1
     for (a in 1:lsars$num.actions) {
          ind <- NULL
          if (!is.null(pick.reward)) {
               if (include.reward) ind <- lsars[[a]]$r == pick.reward
               else ind <- lsars[[a]]$r != pick.reward
               ind <- (1:length(lsars[[a]]$r))[ind]
               if (length(ind) / ms[a] > perc) {
                  sz <- floor(perc * ms[a])
                  ind <- ind[sample(1:length(ind), sz)]
                  }
               ms[a] <- ms[a] - length(ind)
               }
          ind <- c(ind, clara(lsars[[a]]$s2, ms[a])$i.med)
          s <- lsars[[a]]$s[ind,]
          r <- lsars[[a]]$r[ind]
          s2 <- lsars[[a]]$s2[ind,]
          g <- lsars[[a]]$g[ind]
          nl <- c(nl, list(list(s=s, r = r, s2 = s2, g = g)))
          }
     c(nl, list(means = lsars$means, stdevs = lsars$stdevs,
               num.actions = lsars$num.actions))
     }



exp.kbrl <- function(filename, 
                     sizes, num.archs, taus, 
                     transition.function, A, 
                     S.test, max.steps,
                     df, max.iter.pi, pick,
                     use.all.data,
                     num.avg = 30, 
                     archetypes.function = define.archetypes.kbrl, 
                     rbfn.function = make.rbfn.simple, 
                     same.width = FALSE,
                     num.neighbors = ncol(S.test) + 1,
                     run.value.iteration = FALSE,
                     include.reward =TRUE, 
                     perc = 0.3,
                     epsilon = 0.05,
                     verbose = TRUE, 
                     ...){
# if "use.all.data == T", num.archs will be ignored                     
   for (i in 1:length(sizes)) {
      lim.sup <- 1
      if (!use.all.data) lim.sup <- length(num.archs)
      for (k in 1:lim.sup) {
         na <- sizes[i]
         if (!use.all.data) na <- num.archs[k]
         if (verbose) {
            print(paste("Solving KBRL with",sizes[i],"transitions and",
                        na, "archetypes"))
            }
         res    <- matrix(0, num.avg, length(taus))          
         res.st <- matrix(0, num.avg, length(taus))          
         res.ep <- matrix(0, num.avg, length(taus))
         times  <- matrix(0, num.avg, length(taus))          
         for (j in 1:num.avg) {
            if (verbose) print(paste("Run #",j,"of",num.avg))

            lsars <- make.lsars(load.sars(paste(filename, "_data_",
                                 sizes[i],"_", j,sep="")))
            tm <- 0
            lsars2 <- NULL
            if (!use.all.data) {
               tm <- system.time(
               lsars2 <- archetypes.function(lsars, na, pick=pick,
                              include.reward=include.reward, perc=perc), 
                              TRUE)[1]
               }
           else {
               sars <- make.sars(lsars)
               # remove duplicates
               inds <- !duplicated(sars$s2) & !duplicated(sars$s)
               sars$s <- sars$s[inds,]
               sars$a <- sars$a[inds]
               sars$r <- sars$r[inds]
               sars$s2 <- sars$s2[inds,]
               lsars2 <- make.lsars(sars)
               }
      
            for (t in 1:length(taus)) {
               times[j,t] <- tm + system.time(
                  rbfns <- kbrl.rbfns(lsars2, taus[t], df, 
                               same.width = same.width,
                               run.value.iteration = run.value.iteration,
                               max.iter.pi = max.iter.pi,
                               rbfn.function=rbfn.function, 
                               num.neighbors = num.neighbors),
                               TRUE)[1]
               res.test <- array(0, nrow(S.test))
               res.test.st <- array(0, nrow(S.test))
               for (l in 1:nrow(S.test)) {
                  R <- control.rbfns(rbfns, S.test[l,],
                           transition.function, A, lsars2$means, lsars2$stdevs,
                           max.steps, df=df, rbfn.out = rbfn.norm.output,
                           epsilon = epsilon, ...)
                  res.test[l] <- R$tr
                  res.test.st[l] <- R$ns
                  }
               res[j,t] <- mean(res.test)
               res.st[j,t] <- mean(res.test.st)
               res.ep[j,t] <- sum(res.test.st == max.steps)/ length(res.test.st)
               if (verbose) print(paste("Tau =", taus[t], res[j,t],
                                        res.st[j,t], res.ep[j,t]))
               }
            if (verbose) {
               b <- which.max(res[j,])
               print(paste("Best: tau =", taus[b], res[j,b],
                              res.st[j,b],res.ep[j,b]))
               }
            }   
          if (verbose) {
            res.m <- apply(res,2,mean)
            b <- which.max(res.m)
            print(paste("Best overall: tau =", taus[b], res.m[b],
               apply(res.st,2,mean)[b], apply(res.ep,2,mean)[b],"(average)"))
           }
          wt(res, paste(filename, "_kbrl_", sizes[i],"_", na,
                           ".txt", sep=""))
          wt(res.st, paste(filename, "_kbrl_steps_", sizes[i],"_", na,
                           ".txt", sep=""))
          wt(res.ep, paste(filename, "_kbrl_episodes_", sizes[i],"_", na,
                           ".txt", sep=""))
          wt(times, paste(filename,"_kbrl_times_", sizes[i],"_",
                           na, ".txt", sep=""))
          }
       }
     if (verbose) print("All done!")
     }


##----------------- LSPI + KBSF ----------------------------------------------

define.archetypes <-function(sars, num.archs,
                              pick.reward = NULL) {
     C <- NULL                            
     if (!is.null(pick.reward)) {
        ind <- sars$r == pick.reward
        si <- sum(ind)
        if (si > 0) {
            C <- apply(matrix(sars$s[ind,], si,ncol(sars$s)),2,mean)
            num.archs <- num.archs - 1
            }
        else {
            print(paste("define.archetypes: No transitions with reward =",
                     pick.reward))
            C <- NULL
            }
        }
     rbind(C, sars$s2[clara(sars$s2, num.archs)$i.med,])
#       kmeans(sars$s2, num.archs)$centers
     }

define.archetypes.kmeans <-function(sars, num.archs,
                              pick.reward = NULL) {
     C <- NULL                            
     if (!is.null(pick.reward)) {
        ind <- sars$r == pick.reward
        si <- sum(ind)
        if (si > 0) {
            C <- apply(matrix(sars$s[ind,], si,ncol(sars$s)),2,mean)
            num.archs <- num.archs - 1
            }
        else {
            print(paste("define.archetypes: No transitions with reward =",
                     pick.reward))
            C <- NULL
            }
        }
     rbind(C, kmeans(sars$s2, num.archs)$centers)
     }


exp.lspi.kbsf <- function(filename, 
                          sizes, num.archs, taus, 
                          transition.function, A,
                          S.test, max.steps,
                          df, max.iter.pi, 
                          pick.reward,
                          tau.a, 
                          lspi.perc.imp, 
                          num.avg = 30, 
                          archetypes.function = define.archetypes,
                          rbfn.function = make.rbfn.simple,
                          same.width = FALSE, 
                          num.neighbors = ncol(S.test) + 1,
                          p.centers = 1,
                          epsilon = 0.05,
                          run.lspi = TRUE,
                          run.kbsf = TRUE,
                          taus.kbsf = NULL, # only if they're different
                          df.kbsf = NULL,   # only if they're different
                          suffix = NULL,
                          verbose = TRUE, 
                          use.c.code = TRUE,
                          load.results = FALSE,
                          start = 1,
                          ...){

   same.taus <- TRUE
   if (is.null(taus.kbsf)) taus.kbsf <- taus
   else same.taus <- FALSE
      
   if (is.null(df.kbsf)) df.kbsf <- df
   
   if (is.null(suffix)) suffix <- ""
   else suffix <- paste(suffix,"_", sep="")
   


   for (i in 1:length(sizes)) {
     for (k in 1:length(num.archs)) {
      if (verbose) {
         print(paste("Solving problem with",sizes[i],"transitions and",
                     num.archs[k],"archetypes"))
         }
       
      res.lspi    <- matrix(0, num.avg, length(taus))          
      res.kbsf    <- matrix(0, num.avg, length(taus.kbsf))          
      res.lspi.st <- matrix(0, num.avg, length(taus))          
      res.kbsf.st <- matrix(0, num.avg, length(taus.kbsf))          
      res.lspi.ep <- matrix(0, num.avg, length(taus))          
      res.kbsf.ep <- matrix(0, num.avg, length(taus.kbsf))          
      times.lspi  <- matrix(0, num.avg, length(taus))          
      times.kbsf  <- matrix(0, num.avg, length(taus.kbsf))          
      
      if (load.results) {
        if (run.lspi) {
            res.lspi <- read.table(paste(filename, "_lspi_", suffix,  sizes[i],"_",
                              num.archs[k], ".txt", sep=""))
            res.lspi.st <- read.table(paste(filename, "_lspi_steps_", suffix, sizes[i],"_",
                              num.archs[k], ".txt", sep=""))
            res.lspi.ep <- read.table(paste(filename, "_lspi_episodes_", suffix, 
                              sizes[i],"_", num.archs[k], ".txt", sep=""))
            times.lspi <- read.table(paste(filename, "_lspi_times_", suffix, sizes[i],"_",
                              num.archs[k], ".txt", sep=""))
            }
         if (run.kbsf) {
            res.kbsf <- read.table(paste(filename, "_kbsf_", suffix, sizes[i],"_",
                              num.archs[k], ".txt", sep=""))
            res.kbsf.st <- read.table(paste(filename, "_kbsf_steps_", suffix,  sizes[i],"_",
                              num.archs[k], ".txt", sep=""))
            res.kbsf.ep <- read.table(paste(filename, "_kbsf_episodes_", suffix, 
                              sizes[i],"_", num.archs[k], ".txt", sep=""))
            times.kbsf <- read.table(paste(filename, "_kbsf_times_", suffix,  sizes[i],"_",
                              num.archs[k], ".txt", sep=""))
            }  
         while (times.lspi[start,1] != 0) start <- start + 1
         }

      
      for (j in start:num.avg) {
         if (verbose) print(paste("Run #",j,"of",num.avg))
            sars <- load.sars(paste(filename, "_data_",sizes[i],"_",j,sep=""))

           tz <- system.time(
            C <- archetypes.function(sars, num.archs[k], pick.reward)
            ,TRUE)[1]
         for (t in 1:max(length(taus), length(taus.kbsf))) {
            # Define RBFNs
            tm <- tz + system.time(
            rbfn.base <- rbfn.function(C, taus[t], length(A),
                                    same.width = same.width,
                                    num.neighbors = num.neighbors)
                              ,TRUE)[1]
            tm.kbsf <- tm
            rbfn.base.kbsf <- rbfn.base
            if (!same.taus) {
               tm.kbsf <- tz + system.time(
               rbfn.base.kbsf <- rbfn.function(C, taus.kbsf[t], length(A),
                                       same.width = same.width,
                                       num.neighbors = num.neighbors)
                                 ,TRUE)[1]
               }
            
            # LSPI
            if (run.lspi && t <= length(taus)) {
               times.lspi[j,t] <- tm + system.time(
                  rbfn <- lspi.rbfn(sars, rbfn.base, df = df, 
                                 num.iterations = max.iter.pi,
                                 precision = 1e-6, verbose = FALSE, 
                                 rbfn.dm = rbfn.norm.design.matrix,
                                 perc.imp = lspi.perc.imp)
                  , TRUE)[1]

               res.lspi.test    <- array(0, nrow(S.test))
               res.lspi.test.st <- array(0, nrow(S.test))
               for (l in 1:nrow(S.test)) {
                     R <- control.rbfn(rbfn, S.test[l,],
                           transition.function, A, sars$means, sars$stdevs,
                           max.steps, df = df, rbfn.out = rbfn.norm.output,
                           epsilon = epsilon, ...)
                  res.lspi.test[l] <- R$tr
                  res.lspi.test.st[l] <- R$ns
                  }
               res.lspi[j,t] <- mean(res.lspi.test)
               res.lspi.st[j,t] <- mean(res.lspi.test.st)
               res.lspi.ep[j,t] <- sum(res.lspi.test.st==
                                    max.steps)/length(res.lspi.test.st)
               if (verbose) print(paste("LSPI: tau =", taus[t],
                           res.lspi[j,t],res.lspi.st[j,t], res.lspi.ep[j,t]))
               }
               
            # KBSF
            if (run.kbsf && t <= length(taus.kbsf)) {
               times.kbsf[j, t] <- tm.kbsf + system.time(
               rbfn <- kbsf(make.lsars(sars), rbfn.base.kbsf, tau.a = tau.a, 
                                 df = df.kbsf, 
                                 run.value.iteration = FALSE, 
                                 iter.pi = max.iter.pi,
                                 rbfn.function = rbfn.function, 
                                 num.neighbors = num.neighbors,
                                 p.centers = p.centers,
                                 same.width = FALSE,
                                 use.c.code = use.c.code)  
                                 , TRUE)[1]
                                 
               res.kbsf.test    <- array(0, nrow(S.test))
               res.kbsf.test.st <- array(0, nrow(S.test))
               for (l in 1:nrow(S.test)) {
                     R <- control.rbfn(rbfn, S.test[l,],
                           transition.function, A, sars$means, sars$stdevs,
                           max.steps, df = df, rbfn.out = rbfn.norm.output,
                           epsilon = epsilon, ...)
                  res.kbsf.test[l]    <- R$tr
                  res.kbsf.test.st[l] <- R$ns
                  }
               res.kbsf[j,t] <- mean(res.kbsf.test)
               res.kbsf.st[j,t] <- mean(res.kbsf.test.st)
               res.kbsf.ep[j,t] <- sum(res.kbsf.test.st==
                                    max.steps)/length(res.kbsf.test.st)
               if (verbose) print(paste("KBSF: tau =", taus.kbsf[t], 
                           res.kbsf[j,t], res.kbsf.st[j,t], res.kbsf.ep[j,t]))  
               }
 
           }
   
         if (run.lspi) {
            wt(res.lspi, paste(filename, "_lspi_", suffix, sizes[i],"_",
                              num.archs[k], ".txt", sep=""))
            wt(res.lspi.st, paste(filename, "_lspi_steps_", suffix,  sizes[i],"_",
                              num.archs[k], ".txt", sep=""))
            wt(res.lspi.ep, paste(filename, "_lspi_episodes_", suffix, 
                              sizes[i],"_", num.archs[k], ".txt", sep=""))
            wt(times.lspi, paste(filename, "_lspi_times_", suffix, sizes[i],"_",
                              num.archs[k], ".txt", sep=""))
            }
         if (run.kbsf) {
            wt(res.kbsf, paste(filename, "_kbsf_", suffix, sizes[i],"_",
                              num.archs[k], ".txt", sep=""))
            wt(res.kbsf.st, paste(filename, "_kbsf_steps_", suffix,  sizes[i],"_",
                              num.archs[k], ".txt", sep=""))
            wt(res.kbsf.ep, paste(filename, "_kbsf_episodes_", suffix, 
                              sizes[i],"_", num.archs[k], ".txt", sep=""))
            wt(times.kbsf, paste(filename, "_kbsf_times_", suffix, sizes[i],"_",
                              num.archs[k], ".txt", sep=""))
            }  
         }
         
      if (verbose) {
         if (run.lspi) {
            res.m    <- apply(res.lspi,2,mean)
            b <- which.max(res.m)
            print(paste("LSPI:: Best overall: tau =", taus[b],res.m[b],
                     apply(res.lspi.st,2,mean)[b],apply(res.lspi.ep,2,mean)[b]))
            }
         if (run.kbsf) {   
            res.m    <- apply(res.kbsf,2,mean)
            b <- which.max(res.m)
            print(paste("KBSF:: Best overall: tau =", taus.kbsf[b], res.m[b],
                  apply(res.kbsf.st,2,mean)[b], apply(res.kbsf.ep,2,mean)[b]))
            }      
         }

      
      }
    }
   print("All done!")
   }


##----------------- PUDDLE WORLD ----------------------------------------------

puddle.generate.data.equi <- function(size, save.data = FALSE) {
   x <- seq(0, 1, l = floor(sqrt(size %/% 4)))
   peq <- collect.transitions(make.sa(gp(x,x),4), puddle.world.transition,
            c("N", "S", "E", "W"))
   if (save.data) save.sars(peq,"./kbsf/puddle/equi/puddle_data")
   peq
   }


puddle.generate.data.uni <- function(size, save.data = FALSE,
                              p = 0.9) {
   n <- floor(p * (size/4))
   m <- ceiling((1-p) * (size/4))
   S <- cbind(runif(n), runif(n))
   S <- rbind(S, cbind(runif(m, 0.85, 1), runif(m, 0.85, 1))) # rewards
   pun <- collect.transitions(make.sa(S,4), puddle.world.transition,
            c("N", "S", "E", "W"))
   if (save.data) save.sars(pun,"./kbsf/puddle/uni/puddle_data")
   pun
   }


puddle.generate.data.exp <- function(size, save.data = FALSE, 
                              steps=300) {
   mex <- collect.episodes(puddle.world.transition, size, steps, 
         c("N", "S", "E", "W") , c(0,0), c(1,1))
   if (save.data) save.sars(mex,"./kbsf/puddle/exp/puddle_data")
   mex
   }

  

puddle.exp.kbrl <- function(filename, 
                  sizes = seq(200, 2000,by=200), 
                  num.archs = NULL,
                  taus = c(1e-4,1e-3, 1e-2, 1e-1, 0.3, 0.5, 0.7,0.9,0.95,0.99),
                  use.all.data = TRUE,
                  num.avg = 30,
                  verbose = TRUE) {
   exp.kbrl(filename, 
            sizes = sizes, 
            num.archs = num.archs,
            taus = taus,
            transition.function = puddle.world.transition, 
            A = c("N", "S", "E", "W"), 
            S.test = rbind(gp(seq(0.1,0.3,l=3), seq(0.3,0.5,l=3)), 
                           gp(seq(0.1,0.3,l=2), seq(0.9,1,l=2))), 
            max.steps = 300,
            df = 0.95, 
            max.iter.pi = 30, 
            pick = NULL, 
            use.all.data = use.all.data,
            num.avg = num.avg, 
            verbose = verbose)
            }


puddle.exp.lspi.kbsf <- function(filename, 
                  sizes = 2000, 
                  num.archs = seq(10, 100, by=10),
                  taus = c(1e-4,1e-3, 1e-2, 1e-1, 0.3, 0.5, 0.7,0.9,0.95,0.99),
                  tau.a = 0.01, 
                  num.avg = 30,
                  verbose = TRUE, ...) {
      exp.lspi.kbsf(filename,                          
            sizes = sizes, 
            num.archs = num.archs,
            taus = taus,
            transition.function = puddle.world.transition, 
            A = c("N", "S", "E", "W"), 
            S.test = rbind(gp(seq(0.1,0.3,l=3), seq(0.3,0.5,l=3)), 
                           gp(seq(0.1,0.3,l=2), seq(0.9,1,l=2))), 
            max.steps = 300, 
            df = 0.95,   
            max.iter.pi = 30, 
            pick.reward = 5,
            tau.a = tau.a, 
            num.avg = num.avg,
            lspi.perc.imp = 0.3,
            verbose = verbose, ...)
     }

      


puddle.exp.grid <- function(filename, Q, 
                        df = 0.95,
                        num.avg = 30,
                        S.test = rbind(gp(seq(0.1,0.3,l=3),seq(0.3,0.5,l=3)), 
                                        gp(seq(0.1,0.3,l=2), seq(0.9,1,l=2))),
                        max.steps = 300
                        ) {
# Q should be computed by "value.iterarion.2D.Q"
   
    res <- array(0, num.avg)
    res.st <- array(0, num.avg)
    for (i in 1:num.avg) {
      R <- puddle.test.policy(Q, S.test, df)
      res[i]    <- mean(R$tr)
      res.st[i] <- mean(R$st)
      }
    wt(res, paste(filename,"_grid_",nrow(Q),".txt",sep=""))  
    wt(res.st, paste(filename,"_grid_steps_",nrow(Q),".txt",sep=""))  
    }
   
   
   
##----------------- MOUNTAIN CAR ----------------------------------------------

mountain.generate.data.uni <- function(size, 
                                 save.data = FALSE,
                                 p = 0.9,  
                                 min.vel = -0.07, 
                                 max.vel = 0.07,
                                 min.pos = -1.2 , 
                                 max.pos = 0.5,
                                 sd = 0.1) {
   n <- floor(p * (size/3))
   m <- floor((1-p) * (size/3))
   
   cp <- runif(n, min.pos, max.pos)
   cv <- runif(n, min.vel, max.vel)
   S <- cbind(cp,cv)

   cp <- runif(m, 0.9  * max.pos, max.pos)  # (s[1] >= 0.45 && s[2] >= 0.01)
   cv <- runif(m, 0.15 * max.vel, max.vel)  # guarantees rewards
   S <- rbind(S, cbind(cp,cv)) # rewards
  
   sa <- make.sa(S,3)
   if (nrow(sa$s) < size) {
      m <- size - nrow(sa$s)
      cp <- runif(m, min.pos, max.pos)
      cv <- runif(m, min.vel, max.vel)
      sa$s <- rbind(sa$s,cbind(cp,cv))
      sa$a <- c(sa$a, sample(1:3, m, TRUE))
      }
 
   mun <- collect.transitions(sa, mountain.car.transition, c(-1,0,1), sd=sd)
            
   if (save.data) save.sars(mun,"./kbsf/mountain/uni/mountain_data")
   mun
   }



mountain.generate.data.exp <- function(size, save.data = FALSE, 
                              steps=300,
                              min.vel = -0.07, 
                              max.vel = 0.07,
                              min.pos = -1.2 , 
                              max.pos = 0.5,
                              sd = 0.1) {
   mex <- collect.episodes(mountain.car.transition, size, steps, 
         c(-1,0,1), c(min.pos,min.vel), c(max.pos,max.vel), sd=sd)
   if (save.data) save.sars(mex,"./kbsf/mountain/exp/mountain_data")
   mex
   }


mountain.exp.kbrl <- function(filename, 
                  sizes = seq(200, 2000,by=200), 
                  num.archs = NULL,
                  taus = c(1e-4,1e-3, 1e-2, 1e-1, 0.3, 0.5, 0.7,0.9,0.95,0.99),
                  use.all.data = TRUE,
                  num.avg = 30,
                  verbose = TRUE) {
   exp.kbrl(filename, 
            sizes = sizes, 
            num.archs = num.archs,
            taus = taus,
            transition.function = mountain.car.transition, 
            sd = 0.1,
            A = c(-1,0,1), 
            S.test = rbind(gp(seq(-1,0.15,l=5),
                           seq(-0.07,0.02,l=5)), c(0,0)),
            max.steps = 300,
            df = 0.99, 
            max.iter.pi = 30, 
            pick = NULL, 
            use.all.data = use.all.data,
            num.avg = num.avg, 
            verbose = verbose)
            }


mountain.exp.lspi.kbsf <- function(filename, 
                  sizes = 2000, 
                  num.archs = seq(10, 150, by=10),
                  taus = c(1e-4,1e-3, 1e-2, 1e-1, 0.3, 0.5, 0.7,0.9,0.95,0.99),
                  tau.a = 0.01, 
                  num.avg = 30,
                  verbose = TRUE, ...) {
      exp.lspi.kbsf(filename,                          
            sizes = sizes, 
            num.archs = num.archs,
            taus = taus,
            transition.function = mountain.car.transition, 
            sd = 0.1,
            A = c(-1,0,1), 
            S.test = rbind(gp(seq(-1,0.15,l=5),
                           seq(-0.07,0.02,l=5)), c(0,0)),
            max.steps = 300, 
            df = 0.99,   
            max.iter.pi = 30, 
            pick.reward = 1,
            tau.a = tau.a, 
            num.avg = num.avg,
            lspi.perc.imp = 0.3,
            verbose = verbose, ...)
     }


mountain.exp.grid <- function(filename, Q, 
                        df = 0.99,
                        num.avg = 30,
                        S.test = rbind(gp(seq(-1,0.15,l=5),
                           seq(-0.07,0.02,l=5)), c(0,0)),
                        max.steps = 300,
                        sd = 0.1
                        ) {
# Q should be computed by "value.iterarion.2D.Q"
   
    res <- array(0, num.avg)
    res.st <- array(0, num.avg)
    for (i in 1:num.avg) {
      R <- mountain.test.policy(Q, S.test, df, sd=sd)
      res[i]    <- mean(R$tr)
      res.st[i] <- mean(R$st)
      }
    wt(res, paste(filename,"_grid_",nrow(Q),".txt",sep=""))  
    wt(res.st, paste(filename,"_grid_steps_",nrow(Q),".txt",sep=""))  
    }


##----------------- POLE BALANCING----------------------------------------------

pole.generate.data.uni <- function(size, 
                                 save.data = FALSE,
                                 perc = 1,
                                 perc2 = 0,
                                 max.pos        = 2.4, 
                                 max.vel        = 2.4,
                                 max.angle      = (pi/5), 
                                 max.angle.vel  = (pi/5)) {
                                 
   n <- floor((size/3))
#    m <- floor((1-p) * (size/3))

   pos <- runif(n, perc * (-max.pos), perc * max.pos)
   vel <- runif(n, perc * (-max.vel), perc * max.vel)
   angle1     <- runif(n, perc * (-max.angle), perc * max.angle)
   angle1.vel <- runif(n, perc * (-max.angle.vel), perc * max.angle.vel)
   angle2     <- NULL
   angle2.vel <- NULL
#    if (double) {
#       angle2     <- runif(n, perc2 * (-max.angle2), perc2 * max.angle2)
#       angle2.vel <- runif(n, perc2 * (-max.angle.vel2), perc2 *
# max.angle.vel2)
#       }
      
   S <- cbind(pos, vel, angle1, angle1.vel)
#    if(double) S <- cbind(S, angle2, angle2.vel)

   sa <- make.sa(S, 2)
   
   if (nrow(sa$s) < size) {
      m <- size - nrow(sa$s)
      pos <- runif(m, perc * (-max.pos), perc * max.pos)
      vel <- runif(m, perc * (-max.vel), perc * max.vel)
      angle1     <- runif(m, perc * (-max.angle), perc * max.angle)
      angle1.vel <- runif(m, perc * (-max.angle.vel), perc * max.angle.vel)
      angle2     <- NULL
      angle2.vel <- NULL
# #       if (double) {
# #          angle2     <- runif(m, perc2 * (-max.angle2), perc2 * max.angle2)
# #          angle2.vel <- runif(m, perc2 * (-max.angle.vel2),
# perc2*max.angle.vel2)
# #          }
      S <- cbind(pos, vel, angle1, angle1.vel)
#       if(double) S <- cbind(S, angle2, angle2.vel)
      sa$s <- rbind(sa$s, S)
      sa$a <- c(sa$a, sample(1:2, m, TRUE))
      }

   pun <- collect.transitions(sa, pole.transition, c(-10, 10)) ##FIX
            
   if (save.data) save.sars(mun,"./kbsf/pole/uni/pole_data")
   pun
   }



pole.generate.data.exp <- function(size, 
                              save.data = FALSE, 
                              steps = 15000,
                              perc = 0.75,
                              noise.sd.pole = 0,
                              max.pos        = 2.4, 
                              max.vel        = 2.4,
                              max.angle      = (pi/5), 
                              max.angle.vel  = (pi/5)) {

  slim <- perc * c(max.pos,max.vel, max.angle, max.angle.vel)
#   if (double) slim <- c(slim, max.angle2, max.anlge.vel2)
  pex <- collect.episodes(pole.transition, size, steps, 
         c(-10,10), -slim, slim, noise.sd = noise.sd.pole)
         
   if (save.data) save.sars(pex,"./kbsf/pole/exp/pole_data")
   pex
   }
   

pole.make.St <- function(n.grid,
                       p = 0.75,
                       max.pos = c(2.4, 2.4, (pi/5), (pi/5)),
                       min.pos = -c(2.4, 2.4, (pi/5), (pi/5)),
                       d = 4) {
   max.pos <- p*max.pos
   min.pos <- p*min.pos
   grid <- matrix(0, n.grid^d, d)
   uniform.grid(grid, max.pos, min.pos)    ##THIS IS WRONG (THERE'S A DELTA)
   }


pole.exp.kbrl <- function(filename, 
                  sizes = seq(200, 2000,by=200), 
                  num.archs = NULL,
                  taus = c(1e-4,1e-3, 1e-2, 1e-1, 0.3, 0.5, 0.7,0.9,0.95,0.99),
                  use.all.data = TRUE,
                  num.avg = 30,
                  max.steps = 6000,
                  verbose = TRUE) {
   exp.kbrl(filename, 
            sizes = sizes, 
            num.archs = num.archs,
            taus = taus,
            transition.function = pole.transition, 
            A = c(-10,10), 
            S.test = pole.make.St(3), ##Attention DOUBLE
            max.steps = max.steps, 
            df = 0.99, 
            max.iter.pi = 30, 
            pick = NULL, 
            use.all.data = use.all.data,
            num.avg = num.avg, 
            verbose = verbose,
            epsilon = 0)
            }
   


pole.exp.lspi.kbsf <- function(filename, 
                  sizes = 150000, 
                  num.archs = seq(10, 150, by=10),
                  taus = c(0.9, 0.95, 0.99, 0.995, 0.999),
                  taus.kbsf = c(1e-4, 1e-3, 1e-2, 1e-1, 0.3),
                  tau.a = 0.9, 
                  num.avg = 30,
                  max.steps = 3000,
                  lspi.perc.imp = 0.3,
                  max.iter.pi = 30,
                  p.centers = 0.1,
                  df = 0.9,  
                  df.kbsf = 0.9,
                  noise.sd.pole = 1, 
                  verbose = TRUE,...) {
      exp.lspi.kbsf(filename,                          
            sizes = sizes, 
            num.archs = num.archs,
            taus = taus,
            taus.kbsf = taus.kbsf,
            transition.function = pole.transition, 
            A = c(-10,10), 
            S.test = pole.make.St(3), 
            max.steps = max.steps, 
            df = df,   
            df.kbsf = df.kbsf,
            max.iter.pi = max.iter.pi, 
            pick.reward = NULL,
            tau.a = tau.a, 
            num.avg = num.avg,
            lspi.perc.imp = lspi.perc.imp,
            p.centers = p.centers,
            verbose = verbose, 
            epsilon = 0, 
            noise.sd = noise.sd.pole, ...)
            }
            


pole.script <- function(filename,  filename.pre,
                  sizes = 100000, 
                  num.archs = seq(20, 100, by=20),
                  taus = c(1e-4,1e-3, 1e-2, 1e-1, 0.3, 0.9,0.95,0.99, 0.999),
                  verbose = TRUE) {
# This function uses the preliminary experiments to select the best tau 
#  for each configuration and then calls the main function
   for (s in sizes) {
      best.taus.kbsf <- array(0, length(num.archs))
      best.taus.lspi <- array(0, length(num.archs))
      for (i in 1:length(num.archs)) {
         T <- read.table(paste(filename.pre, "_lspi_", s, "_", num.archs[i],
                          ".txt", sep=""))
         best.taus.lspi[i] <- which.max(apply(T, 2, mean))
         
         T <- read.table(paste(filename.pre, "_kbsf_", s, "_", num.archs[i],
                          ".txt", sep=""))
         best.taus.kbsf[i] <- which.max(apply(T, 2, mean))
         }
      print(best.taus.lspi)
      print(best.taus.kbsf)
      for (na in 1:length(num.archs)) {
         if (verbose) {
            print(paste("Size:",s,"m:",num.archs[na]," LSPI's tau:",
             taus[best.taus.lspi[na]], " KBSF's tau:",taus[best.taus.kbsf[na]]))
             }
         pole.exp.lspi.kbsf(filename,  sizes = s, num.archs = num.archs[na],
         taus = taus[best.taus.lspi[na]], taus.kbsf=taus[best.taus.kbsf[na]])
         }
      }
    
  }
         
         
##------------DOUBLE POLE BALANCING--------------------------------------------


pole.double.generate.data.exp <- function(size, 
                              save.data = FALSE, 
                              steps = 3000,
                              perc = 0.3,
                              A = c(-10,10),
                              transition.function = pole.double.transition,
                              noise.sd.pole = 0,
                              max.pos        = 2.4, 
                              max.vel        = 2.4,
                              max.angle      = (pi/5), 
                              max.angle.vel  = (pi/5),
                              max.angle2      = 0, 
                              max.angle.vel2  = 0) {

  slim <- perc * c(max.pos,max.vel, max.angle, max.angle.vel, max.angle2,
                    max.angle.vel2)
#   if (double) slim <- c(slim, max.angle2, max.anlge.vel2)
   pex <- collect.episodes(transition.function, size, steps, 
         A, -slim, slim, noise.sd = noise.sd.pole)
         
   if (save.data) save.sars(pex,"./kbsf/double_pole/exp/pole_double_data")
   pex
   }
   

pole.double.make.St <- function(n.grid,
                       p = 0.3,
                       max.pos = c(2.4, 2.4, (pi/5), (pi/5)),
                       min.pos = -max.pos,
                       d=6) {
   max.pos <- p*max.pos
   min.pos <- p*min.pos
   grid <- matrix(0, n.grid^(d-2), d-2)
   S <- uniform.grid(grid, max.pos, min.pos) ## THIS IS WRONG (THERE'S A DELTA)
   cbind(S,0,0)
   }


pole.double.make.St.harder <- function(n.grid,
                       p = 0.3,
                       max.pos = c(2.4, 2.4, (pi/5), (pi/5), (pi/5), (pi/5)),
                       min.pos = -max.pos,
                       d=6) {
   max.pos <- p*max.pos
   min.pos <- p*min.pos
   grid <- matrix(0, n.grid^d, d)
   S <- uniform.grid(grid, max.pos, min.pos)
   S
   }


pole.double.exp.lspi.kbsf <- function(filename, 
                  sizes = 5e5, 
                  num.archs = seq(20, 250, by=20),
                  taus = c(0.9, 0.95, 0.99, 0.995, 0.999),
                  taus.kbsf = c(1e-4, 1e-3, 1e-2, 1e-1, 0.3),
                  tau.a = 0.9, 
                  num.avg = 30,
                  max.steps = 3000,
                  A = c(-10,10),
                  lspi.perc.imp = 0.3,
                  max.iter.pi = 30,
                  p.centers = 0.05,
                  df = 0.7,  
                  df.kbsf = 0.7,
                  noise.sd.pole = 0, 
                  archetypes.function = define.archetypes, 
                  S.test = pole.double.make.St(3),
                  transition.function = pole.double.transition,
                  verbose = TRUE,...) {
      exp.lspi.kbsf(filename,                          
            sizes = sizes, 
            num.archs = num.archs,
            taus = taus,
            taus.kbsf = taus.kbsf,
            transition.function = transition.function, 
            A = A, 
            S.test = S.test, 
            max.steps = max.steps, 
            df = df,   
            df.kbsf = df.kbsf,
            max.iter.pi = max.iter.pi, 
            pick.reward = NULL,
            tau.a = tau.a, 
            num.avg = num.avg,
            lspi.perc.imp = lspi.perc.imp,
            p.centers = p.centers,
            verbose = verbose, 
            epsilon = 0, 
            noise.sd = noise.sd.pole,
            archetypes.function = archetypes.function,
            ...)
            }
                          
   


   