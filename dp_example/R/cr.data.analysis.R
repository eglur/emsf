source("util.R")
source("data.plot.R")


generate.all <- function(plot = c(TRUE, TRUE, TRUE), save = c(FALSE, FALSE, FALSE))
{
   # use x11(w=8, h = 8)
   
   PM <- matrix(TRUE, 6, 6)
   PM[5:6, 1] <- FALSE
   PM[6,2] <- FALSE
#     PM[6,4] <- FALSE
   
   print(PM)
   
   if (plot[1]) # sizes
   {
      show.perf(crit = "sizes", crit.test = "sizes", thr = c(200, 400,600), load.alg = c(TRUE, TRUE, FALSE, TRUE), load.test = TRUE, nc = 2:5 , nc.test = c(6,7), load.alg.test = c(FALSE, TRUE, FALSE, TRUE), ylab = "Size", PM = PM)
      
      if (save[1]) dev.copy2eps(file = "~/pisf_jair/fig/cr_sizes.eps")
   }

   if (plot[2])
   {
      x11(w = 8, h = 8)
      show.perf(crit = "gain", crit.test = "gain_mc_5000", thr = c(200, 400,600), load.alg = c(TRUE, TRUE, TRUE, TRUE), load.test = TRUE, nc = 2:5 , nc.test = c(6,7), load.alg.test = c(FALSE, TRUE, TRUE, TRUE), ylab = expression(varphi(v^pi, v^pi[N])), PM = PM)

      if (save[2]) dev.copy2eps(file = "~/pisf_jair/fig/cr_gain.eps")
   }

   if (plot[3])
   {
      x11()
      show.perf(crit = "time", crit.test = "time", thr = c(200, 400,600), load.alg = c(TRUE, TRUE, FALSE, TRUE), load.test = TRUE, nc = 2:5 , nc.test = c(6,7), load.alg.test = c(FALSE, TRUE, FALSE, TRUE), ylab = "Hours",  PM = PM)

      if (save[3]) dev.copy2eps(file = "~/pisf_jair/fig/cr_times.eps")
   }
}

ci <- function(v, confidence.level)
{
#    qnorm(1-(1-confidence.level)/2)*sd(v)/sqrt(length(v))
      sd(v)/sqrt(length(v))
}


show.perf <- function(
               
               load.exact = TRUE,
               crit = "gain",
               nc = 2:5, 
               lt = 10,
			ma = 1:10, # max age for std
			thr = c(200,400,600), # threshold for PISF
			dir = "~/cr/default/experiments/",
			load.alg = c(TRUE, TRUE, TRUE, TRUE), # PI, PI-RED, STD, PISF
               ylim = NULL,
               confidence.level = 0.99, ## not currently used
               leg.position = "topleft",
               log  = "", 
               
               load.test = TRUE,
               crit.test = "gain_mc_5000",
               nc.test  = c(6,7),
               load.alg.test = c(FALSE, FALSE, TRUE, TRUE),

               ylab = NULL,
              
               PM = NULL, # plot matrix; logical nc x num.alg matrix with the plot pattern 
               show.hours = TRUE,
               ...
			)
{
  set.par() 
   

  if (crit == "sizes") load.alg[3] <- FALSE

  num.alg <- 4 #sum(load.alg) 
  if (load.alg[4]) num.alg <- num.alg + length(thr) - 1
  
  delta <- 0
  delta2 <- 0
  if (load.test) delta <- length(nc.test)
  if (load.exact) delta2 <- length(nc)
  E <- matrix(0, delta + delta2, num.alg)
  C <- E
  
  c <- 1
  if (crit == "time") c <- 2
  
     
   if (load.exact)
   {

      for (i in 1:length(nc))
      {

         if (load.alg[1])
         {
            
            ## Policy iteration
            fn <- ps(dir,paste("pinm",crit, nc[i], lt, sep="_"),".txt")
            
            T <- as.matrix(read.table(fn))[,c] ; pd(T); if (crit == "time" && show.hours) T <- T / 60^2
            E[i,1] <- mean(T)
            C[i,1] <- ci(T, confidence.level)
         }
         
         
         if (load.alg[2])
         {
            ## Policy iteration on reduced MDP
            fn <- ps(dir,paste("pinm_red",crit, nc[i], lt, sep="_"),".txt")
            
            T <- as.matrix(read.table(fn))[,c] ; pd(T); if (crit == "time" && show.hours) T <- T / 60^2
            E[i,2] <- mean(T)
            C[i,2] <- ci(T, confidence.level)
         }

         if (load.alg[3])
         {
            ## STD policies
            ind <- NULL
            if (crit == "time") ind <- 1:length(ma)
            else
            {
               v <- array(0, length(ma))
               for (j in 1:length(ma))
               {
                  fn <- ps(dir, paste("std", ma[j], crit, nc[i], lt, sep="_"),".txt")
                                    
                  T <- as.matrix(read.table(fn))[,c] ; pd(T); if (crit == "time" && show.hours) T <- T / 60^2
                  v[j] <- mean(T)
               }
            
               ind <- which.max(v)
            }
               
            T <- 0
            for (j in 1:length(ind))
            {
               fn <- ps(dir, paste("std", ma[ind[j]], crit, nc[i], lt, sep="_"),".txt")
              
               T <- T + as.matrix(read.table(fn))[,c]
               print(paste("Best STD for ", nc[i], " components: ", ma[ind[j]], "    Value: ", mean(T)))
            }

            E[i,3] <- mean(T)
            C[i,3] <- ci(T, confidence.level)
         }

         if (load.alg[4])
         {
            ## PISF policies
            for (j in 1:length(thr))
            {
               fn <- ps(dir,paste("pisf", thr[j], crit, nc[i], lt, sep="_"),".txt")
              
               T <- as.matrix(read.table(fn))[,c] ; pd(T); if (crit == "time" && show.hours) T <- T / 60^2
               E[i,3 + j] <- mean(T)
               C[i,3 + j] <- ci(T, confidence.level)
            }
         }
         
      }
   
      
   }
  
   if (load.test)
   {
      
      if (is.null(load.alg.test)) load.alg.test <- load.alg
         
      for (i in 1:length(nc.test)) 
      {
      
         
      if (load.alg.test[1])
      {
            
         ## Policy iteration
         fn <- ps(dir,paste("pinm",crit.test, nc.test[i], lt, sep="_"),".txt")
         if (!file.exists(fn)) print(paste("Warning: file", fn, "does not exist"))
         else
         {
            T <- as.matrix(read.table(fn))[,c] ; pd(T); if (crit == "time" && show.hours) T <- T / 60^2
            
            E[delta2 + i,1] <- mean(T)
            C[delta2 + i,1] <- ci(T, confidence.level)
         }
      }

          
      if (load.alg.test[2]) 
      {
            
         ## Reduced Policy iteration
         fn <- ps(dir,paste("pinm_red",crit.test, nc.test[i], lt, sep="_"),".txt")
         if (!file.exists(fn)) print(paste("Warning: file", fn, "does not exist"))
         else
         {
            T <- as.matrix(read.table(fn))[,c] ; pd(T); if (crit == "time" && show.hours) T <- T / 60^2
            E[delta2 + i,2] <- mean(T)
            C[delta2 + i,2] <- ci(T, confidence.level)
         }
      }

          
      if (load.alg.test[3])
      {
         ## STD policies
         ind <- NULL
         if (crit == "time_test") ind <- 1:length(ma)
         else
         {
            v <- array(0, length(ma))
            for (j in 1:length(ma))
            {
               fn <- ps(dir, paste("std", ma[j], crit.test, nc.test[i], lt, sep="_"),".txt")
               T <- as.matrix(read.table(fn))[,c] ; pd(T); if (crit == "time" && show.hours) T <- T / 60^2
               v[j] <- mean(T)
            }
         
            ind <- which.max(v)
         }
            
         T <- 0
         for (j in 1:length(ind))
         {
            fn <- ps(dir, paste("std", ma[ind[j]], crit.test, nc.test[i], lt, sep="_"),".txt")
            T <- T + as.matrix(read.table(fn))[,c]
         }
         
         E[delta2 + i,3] <- mean(T)
         print(paste(">>>", E[delta2 + i,3] ))
         C[delta2 + i,3] <- ci(T, confidence.level)

      }
      
      if (load.alg.test[4])
      {
         ## PISF policies
         for (j in 1:length(thr))
         {
#             if (thr.test[j]) ## hack
            {
               
               fn <- ps(dir,paste("pisf", thr[j], crit.test, nc.test[i], lt, sep="_"),".txt")
               if (!file.exists(fn)) print(paste("Warning: file", fn, "does not exist"))
               else
               {
                  T <- as.matrix(read.table(fn))[,c] ; pd(T); if (crit == "time" && show.hours) T <- T / 60^2
                  
                  ## hack
                  T <- T[T > 0] # for some reasons some of PISF's runs are failing, with results that are worse than the naive policy (std 0)
                  ## end hack
                  
                  E[delta2 + i, 3 + j] <- mean(T)
                  C[delta2 + i, 3 + j] <- ci(T, confidence.level)
               }
            }
            
         }
      }


      }
   }
  
  print(E)
  print("")
  print(C)
  if (log != "") E <- E + 1e-3
  C[is.na(C)] <- 0
  
  nce <- NULL
  if (load.exact) nce <- nc
  if (load.test) nce <- c(nce, nc.test)
 
  if (is.null(PM)) PM <- matrix(TRUE, length(nce), num.alg)
  for (i in 1:length(load.alg)) if (!load.alg[i]) PM[1:length(nc),i] <- FALSE
  for (i in 1:length(load.alg.test)) if (!load.alg.test[i]) PM[(length(nc)+1):length(nce),i] <- FALSE

   ylim <- c(Inf, -Inf)
   for (i in length(nce):1) # first pass is to set y limits 
   {
      ia <- (1:ncol(PM))[PM[i,]]
      if (length(ia) > 0)
      {
         E2 <- matrix(E[1:i, ia], i, length(ia))
         C2 <- matrix(C[1:i, ia], i, length(ia))
         ylim[1] <- min(ylim[1], min(E2 - C2))
         ylim[2] <- max(ylim[2], max(E2 + C2))
      }
   }
     
  
  if (is.null(ylab)) ylab <- crit
  add <- FALSE
  for (i in length(nce):1)
  {
     ia <- (1:ncol(PM))[PM[i,]]
     if (length(ia) > 0)
     {
         x <- (nce[1]):(nce[i])
         E2 <- matrix(E[1:i, ], i, ncol(E))
         C2 <- matrix(C[1:i, ], i, ncol(E))
         mp(x, E2, E2 + C2, E2 - C2, t="o", ylim = ylim, ylab = ylab, xlab = expression("Number of components (" * n[c] * ")"), log=log, inds = ia, add = add, show.shadow = FALSE, ...)
     }
     
     add <- TRUE
  }
  
  txt <- c("PI", "PI-RED", "BEST THR", "PISF")[load.alg]
  if (load.alg[4]) 
  {
     l <- paste("PISF", thr, sep="-")
     txt<- c(txt[-length(txt)], l)
  }
 
  pa <- load.alg
  if (load.alg[4]) 
  {
     pa <-c(pa, rep(TRUE, length(thr) - 1))
  }
  col = get.col(ncol(E))[pa]
  pch = get.pch(ncol(E))[pa]

  leg(leg.position, txt, col = col, pt.bg = col, pch = pch)
  E
}
    




show.perf.test <- function(
               crit = "return_test",
               nc = 7, 
               lt = 10,
               ma = 0:10, # max age for std
               thr = 300, # threshold for PISF
               hor = 5000,  # horizon used to evaluate the policies
               dir = "~/cr/default/experiments/",
               leg.position = "topleft",
               num.avg  = 91, # some runs have not finished yet
               log  = "",
               ...
               )
{
   X <- matrix(0, num.avg, length(ma) + length(thr))
   
   for (i in 1:length(ma))
   {
      fn <- ps(dir, paste("std", ma[i], crit, hor, nc, lt, sep="_"),".txt")
      X[,i] <- read.table(fn)[1:num.avg,1]
   }
   
  
   for (i in 1:length(thr))
   {
      fn <- ps(dir, paste("pisf", thr[i], crit, hor, nc, lt, sep="_"),".txt")
      T <- read.table(fn)[,1]
      X[,length(ma) + i] <- read.table(fn)[1:num.avg,1]
   }

   boxplot(X, range = 0)

   
}




show.perf.test.policy <- function(
               crit = "return_test",
               nc = 2:4, 
               lt = 10,
               ma = 4, # max age for std
               max.depth = c(5, 50, 500, 5000), # depth of UCT's trees
               C = c(1, 10, 100, 1000, 10000),       # UCT's exploration parameter
               thr = 500, # threshold for PISF
               hor = 5000,  # horizon used to evaluate the policies
               dir = "~/cr/default/experiments/",
               load.alg = c(FALSE, TRUE, TRUE),
               ylim = NULL,
               confidence.level = 0.99,
               leg.position = "topleft",
               log  = "",
               use.std.ref = TRUE,  #should the results be with respect to the naive STD-0 policy?
               ...
               )
{
  set.par() 

  num.alg <- sum(load.alg)
  if (load.alg[3]) num.alg <- num.alg + length(thr) - 1
  
  E <- matrix(0, length(nc), num.alg)
  F <- E
  
  c <- 1
  
  for (i in 1:length(nc))
  {
    
     ref <- NULL
     if (use.std.ref)
     {
         fn <- ps(dir, paste("std", 0, crit, hor, nc[i], lt, sep="_"),".txt")
         ref <- as.matrix(read.table(fn))[,c]
     }
     
    
     col <- 1 

    if (load.alg[1])
    {
      ## STD policies
      ind <- NULL

#       if (crit == "time_test") ind = 1:length(ma)
#       else
      {
         v <- array(0, length(ma))
         for (j in 1:length(ma))
         {
            fn <- ps(dir, paste("std", ma[j], crit, hor, nc[i], lt, sep="_"),".txt")
            T <- as.matrix(read.table(fn))[,c] ; pd(T); if (crit == "time" && show.hours) T <- T / 60^2
            if (use.std.ref) T <- (T - ref) / abs(ref)
            v[j] <- mean(T)
         }
      
         ind <- which.max(v)
      }
         
      T <- 0
      for (j in 1:length(ind))
      {
         fn <- ps(dir, paste("std", ma[ind[j]], crit, hor, nc[i], lt, sep="_"),".txt")
         v <- as.matrix(read.table(fn))[,c]
         if (use.std.ref) v <- (v - ref) / abs(ref)
         T <- T + v
      }

      E[i,col] <- mean(T)
      F[i,col] <- ci(T, confidence.level)

      col <- col + 1
    }

    
    if (load.alg[2])
    {
      ## UCT policies
      max.return <- -Inf # the selection is always based on the return
      best.j <- 0
      best.k <- 0
      
      max.depth2 <- max.depth
      C2 <- C
      
      if (nc[i] == 5) # this is because for nv = 2, 3, 4 the best results were found with max.depth = 5 and C = 1000, so I only ran those for nc = 5
      {
         max.depth2 <- 5
         C2 <- 1000
      }
      
      
      for (j in 1:length(max.depth2))
      {
         for (k in 1:length(C2))
         {
            fn <- ps(dir,paste("uct", max.depth2[j], C2[k], "return_test", hor, nc[i], lt, sep="_"),".txt") # the selection is always based on the return
            T <- as.matrix(read.table(fn))[,c] ; pd(T); if (crit == "time" && show.hours) T <- T / 60^2
            value <- mean(T)
            if (value > max.return)
            {
               max.return <- value
               best.j <- j
               best.k <- k
            }
         }

      }
      
      print(paste("Best depth: ", max.depth2[best.j], "   Best C: ", C2[best.k]))
      
      fn <- ps(dir,paste("uct", max.depth2[best.j], C2[best.k], crit, hor, nc[i], lt, sep="_"),".txt") 
      T <- as.matrix(read.table(fn))[,c] ; pd(T); if (crit == "time" && show.hours) T <- T / 60^2
      if (use.std.ref) T <- (T - ref) / abs(ref)
         
      E[i,col] <- mean(T)
      F[i,col] <- ci(T, confidence.level)
      col <- col + 1

   }

    
    if (load.alg[3])
    {
      ## PISF policies
      for (j in 1:length(thr))
      {
         fn <- ps(dir,paste("pisf", thr[j], crit, hor, nc[i], lt, sep="_"),".txt")
         T <- as.matrix(read.table(fn))[,c] ; pd(T); if (crit == "time" && show.hours) T <- T / 60^2
         if (use.std.ref) T <- (T - ref) / abs(ref)         
         E[i,col] <- mean(T)
         F[i,col] <- ci(T, confidence.level)
         col <- col + 1
      }
    }
    
  }
  
  if (log != "") E <- E + 1e-3
     
  mp(nc, E, E + F, E - F, t="o", ylim = ylim, ylab = ylab, xlab = "Number of components", log=log, ...)

  txt <- c("STD", "UCT", "PISF")[load.alg]
  if (load.alg[3]) 
  {
     l <- paste("PISF", thr)
     txt<- c(txt[-length(txt)], l)
  }
 
  leg(leg.position, txt)
}





show.perf.fac <- function(
               crit = "gain",
               nc = 5, 
               lt = 10,
               thr = c(200, 400, 600), # threshold for PISF
               level.fac = 1:4,
               dir = "~/cr/default/experiments/",
               load.alg = c(TRUE, TRUE, FALSE), # PI-FAC, PISF, PISF-FAC
               ylim = NULL,
               confidence.level = 0.99,
               leg.position = "bottomleft",
               log  = "",
               ylab = expression(varphi(v^pi, v*"*")), 
               ...
               )
{
  set.par() 

  num.alg <- sum(load.alg) 
#   if (load.alg[1]) num.alg <- num.alg + 1 # for pinm fac INCLUDE LATER
  if (load.alg[2]) num.alg <- num.alg + length(thr) - 1
  if (load.alg[3]) num.alg <- num.alg + length(thr) - 1
     
  E <- matrix(0, length(level.fac), num.alg)
  C <- E
  
  c <- 1
  if (crit == "time") c <- 2
  
     
  for (i in 1:length(level.fac))
  {
    
     col <- 1
    
    if (load.alg[1])
    {
       
      ## Policy iteration on factored MDP
      fn <- ps(dir,paste("pinm_fac",level.fac[i], crit, nc, lt, sep="_"),".txt")
      T <- as.matrix(read.table(fn))[,c] ; pd(T); if (crit == "time" && show.hours) T <- T / 60^2
      E[i,col] <- mean(T)
      C[i,col] <- ci(T, confidence.level)
      col <- col + 1
    }

    if (load.alg[2])
    {
      ## PISF policies
      for (j in 1:length(thr))
      {
         fn <- ps(dir,paste("pisf", thr[j], "fac_0", crit, nc, lt, sep="_"),".txt") # fac_0 is the regular PISF; it was used to compare with PINM instead of STD 0
         T <- as.matrix(read.table(fn))[,c] ; pd(T); if (crit == "time" && show.hours) T <- T / 60^2
         E[i,col] <- mean(T)
         print(E[i,col])
         C[i,col] <- ci(T, confidence.level)
         col <- col + 1
      }
    }

    if (load.alg[3])
    {
      ## PISF policies on factored MDP
      for (j in 1:length(thr))
      {
         fn <- ps(dir,paste("pisf", thr[j], "fac", level.fac[i], crit, nc, lt, sep="_"),".txt")
         T <- as.matrix(read.table(fn))[,c] ; pd(T); if (crit == "time" && show.hours) T <- T / 60^2
         E[i,col] <- mean(T)
         print(E[i,col])
         C[i,col] <- ci(T, confidence.level)
         col <- col + 1
      }
    }
    
  }
  
  if (log != "") E <- E + 1e-3
     
  ylim <- c(min(E-C), max(E+C) + 0.002)
  
  F <- matrix(E[,2:ncol(E)], nrow(E), ncol(E)-1)
  G <- matrix(C[,2:ncol(E)], nrow(E), ncol(E)-1)
  
  col <- get.col(4 + length(thr) - 1)[4:(4 + length(thr) - 1)]
  mp(level.fac, F, F + G, F - G, t="o", ylim = ylim, ylab = ylab, xlab = "Level of sparsity (z)", 
     log=log, show.shadow = TRUE, col = col, pch = NA, lty = 1, 
     show.error = FALSE, transparency = FALSE, xaxt = "n", ...)
  
  I <- matrix(E[,1], nrow(E), 1)
  J <- matrix(C[,1], nrow(E), 1)
  mp(level.fac, I, I + J, I - J, t="o", ylim = ylim, ylab = ylab, xlab = "Level of sparsity (z)", 
     log=log, show.shadow = FALSE, add = TRUE, ...)
  axis(1, at = level.fac)



  
  txt <- "PI-FAC"

  if (load.alg[2]) 
  {
     l <- paste("PISF-", thr, sep = "")
     txt<- c(txt, l)
  }

  if (load.alg[3]) 
  {
     l <- paste("PISF-FAC", thr)
     txt<- c(txt, l)
  }

#   leg(leg.position, txt)
   text(3.5, E[3,1] - 0.0015, txt[1])
   for (i in 1:ncol(F)) text(3.5, F[1,i] + 0.0015, txt[1 + i])
      
}
    

    

show.perf.truly.fac <- function(
               crit = "gain",
               nc = 5, 
               lt = 10,
               thr = c(200, 400, 600), # threshold for PISF
               ma = 1:10,
               level.fac = 1:4,
               dir = "~/cr/default/experiments/",
               load.alg = c(TRUE, TRUE), # STD, PISF
               ylim = NULL,
               confidence.level = 0.99,
               leg.position = "bottomleft",
               log  = "",
               ylab = expression(varphi(v^pi, v[z]*"*")), 
               ...
               )
{
  set.par() 

  num.alg <- sum(load.alg) 

  if (load.alg[1]) num.alg <- num.alg + length(thr) - 1
    
  E <- matrix(0, length(level.fac), num.alg)
  C <- E
  
  c <- 1
  if (crit == "time") c <- 2
  
     
  for (i in 1:length(level.fac))
  {
    
     col <- 1
    
    if (load.alg[1])
    {
      ## STD 
      
      v <- array(0, length(ma))
      for (j in 1:length(ma))
      {
         fn <- ps(dir,paste("std", ma[j], "truly_fac", level.fac[i], crit, nc, lt, sep="_"),".txt")
         T <- as.matrix(read.table(fn))[,c]
         v[i] <- mean(T)
      }

      t <- which.max(v)
      
      for (j in t)
      {
         fn <- ps(dir,paste("std", ma[j], "truly_fac", level.fac[i], crit, nc, lt, sep="_"),".txt")
         T <- as.matrix(read.table(fn))[,c]
         pd(T)
         if (crit == "time" && show.hours) T <- T / 60^2
         E[i,col] <- mean(T)
         print(E[i,col])
         C[i,col] <- ci(T, confidence.level)
         col <- col + 1
      }
      

    }

    if (load.alg[2])
    {
     
      ## PISF policies
      for (j in 1:length(thr))
      {
         fn <- ps(dir,paste("pisf", thr[j], "truly_fac", level.fac[i], crit, nc, lt, sep="_"),".txt") 
         T <- as.matrix(read.table(fn))[,c]
         pd(T)
         if (crit == "time" && show.hours) T <- T / 60^2
         E[i,col] <- mean(T)
         C[i,col] <- ci(T, confidence.level)
         col <- col + 1
      }
    }
    
  }
  
  if (log != "") E <- E + 1e-3

  col <- get.col(3 + length(thr))[3:(3 + length(thr))]
  pch <- get.pch(3 + length(thr))[3:(3 + length(thr))]
  
  mp(level.fac, E, E + C, E - C, t="o", ylim = ylim, ylab = ylab, xlab = "Level of sparsity (z)", log=log, 
     show.shadow = FALSE, xaxt = "n", col = col, pch = pch, ...)
  axis(1, at = level.fac)

  txt <- NULL

#   if (load.alg[1]) 
#   {
#      txt<- c(txt, "BEST THR")
#   }

  if (load.alg[2]) 
  {
     l <- paste("PISF", thr)
     txt<- c(txt, l)
  }

  
  col <- col[2:length(col)]
  pch <- pch[2:length(pch)]
  leg(leg.position, txt, col = col, pt.bg = col, pch = pch)
  
  text(level.fac[3] + 0.7, E[3,1] + 0.015, "BEST THR")
  
}


show.perf.uct <- function(
               crit = "gain_mc_5000",
               nc = 5, 
               lt = 10,
               thr = 600, # threshold for PISF
               max.secs = seq(1, 15, by = 2),
               C = c(1, 10, 100, 1000),
               max.depth = 5,
               dir = "~/cr/default/experiments/",
               ylim = NULL,
               confidence.level = 0.99,
               leg.position = "bottomleft",
               log  = "",
               ylab = expression(varphi(v^pi, v^pi[N])), 
               show.hours = FALSE,
               plot.vertical.time = TRUE,
               ...
               )
{
  set.par() 
     
  E <- matrix(0, length(max.secs), 1 + length(thr))
  D <- E
  
  # for times
  g <- array(0, length(thr))
  h <- g
  
  c <- 1
  if (crit == "time") c <- 2
       
  ## UCT
  for (i in 1:length(max.secs))
  {
      v <- array(0, length(C))
      for (j in 1:length(C))
      {
         fn <- ps(dir,paste("uct", max.depth, C[j], max.secs[i],crit,  nc, lt, sep="_"),".txt") 
         T <- as.matrix(read.table(fn))[,c] 
         pd(T)
         v[j] <- mean(T)
      }
      
      print(v)
      bc <- which.max(v)
      print(paste("Max secs", max.secs[i], "  Best C ", C[bc]))
      
      fn <- ps(dir,paste("uct", max.depth, C[bc], max.secs[i], crit, nc, lt, sep="_"),".txt") 
      T <- as.matrix(read.table(fn))[,c] 
      E[i,1] <- mean(T)
      D[i,1] <- ci(T, confidence.level)
   }

   ## PISF policies
   for (j in 1:length(thr))
   {
      fn <- ps(dir,paste("pisf", thr[j], crit, nc, lt, sep="_"),".txt")      
      T <- as.matrix(read.table(fn))[,c]
      if (crit == "time" && show.hours) T <- T / 60^2
      E[,1 + j] <- mean(T)
      D[,1 + j] <- ci(T, confidence.level)

      
      if (plot.vertical.time)
      {
         fn <- ps(dir,paste("pisf", thr[j], "time", nc, lt, sep="_"),".txt")      
         T <- as.matrix(read.table(fn))[,2]
         print(T)
         g[j] <- mean(T)
         h[j] <- ci(T, confidence.level)
      }
      
   }
      
#  mp(max.secs, E, E + D, E - D, ylab = ylab, xlab = "Seconds per step", log = log, ...) 
  
  col <- get.col(4 + length(thr))[4:(4 + length(thr))]
  pch <- get.pch(4 + length(thr))[4:(4 + length(thr))]
  ylim <- c(min(E-D), max(E+D) + 0.002)

  F <- matrix(E[,2:ncol(E)], nrow(E), ncol(E)-1)
  G <- matrix(D[,2:ncol(E)], nrow(E), ncol(E)-1)
  
#   col <- get.col(4 + length(thr) - 1)[4:(4 + length(thr) - 1)]

  mp(max.secs, F, F + G, F - G, t="o", ylim = ylim, ylab = ylab, xlab = expression("Seconds per step ("*t[max]*")"), 
     log=log, show.shadow = TRUE, col = col, pch = NA, lty = 1, 
     show.error = FALSE, transparency = FALSE, xaxt = "n", ...)
  
  I <- matrix(E[,1], nrow(E), 1)
  J <- matrix(D[,1], nrow(E), 1)
  mp(max.secs, I, I + J, I - J, t="o", ylim = ylim, ylab = ylab, xlab = expression("Seconds per step ("*t[max]*")"),
     log=log, show.shadow = FALSE, add = TRUE, show.error = TRUE, ...)
  axis(1, at = max.secs)



  
      txt <- "UCT"

     l <- paste("PISF-", thr, sep = "")
     txt<- c(txt, l)

     l <- paste("PISF-", thr)
     txt<- c(txt, l)

#   leg(leg.position, txt)
   text(3.5, E[3,1] - 0.0015, txt[1])
   for (i in 1:ncol(F)) text(3.5, F[1,i] - 0.04, txt[1 + i])
   
   if (plot.vertical.time)
   {
      for (j in 1:length(thr))
      {
         abline(v = g[j])
         print(paste("time", g[j]))
      }
   }

   
   E

   

}
    
    
std.values <- function(nc = 2:5, ma = 1:10, dir = "~/cr/default/experiments/")
{
   D <- matrix(0, length(ma), length(nc))
   for (i in 1:length(nc))
   {
      for (j in 1:length(ma))
      {
         fn <- ps(dir, paste("std", ma[j], "gain", nc[i], "10", sep="_"),".txt")
                           
         T <- as.matrix(read.table(fn))[,1] 
         D[j,i] <- mean(T)
      }
   }
   D
     
}

print("cr.data.analysis.R")
