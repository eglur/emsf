library(plotrix)
source("kbsf.experiments.R")

set.par <- function() {
# set parameters used  in the plots
   par(
      lwd = 1.5,
      mar = c(5, 6, 4, 2) + 0.1,
      ps  = 20)
   }

num.algs <- 5

alg.names <- function(n = num.algs) {
   names <- c("KBRL", "LSPI", "KBSF")
   names[1:n]
   }

alg.lty <- function(n = num.algs) {
   types <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
   types[1:n]
   }

alg.cols <- function(n = num.algs) {
   col <- c("BLACK", "RED", "BLUE", "GREEN3", "VIOLET")
   col[1:n]
   }

alg.pch <- function(n = num.algs) {
   seq(21, 21 + n -1, by = 1)
   }


plot.alg <- function(x, Y, alg.ind,
               Y.lb = NULL,
               Y.ub = NULL,
               plot.leg = TRUE,
               leg = NULL,
               leg.pos = "topright",
               lty = NULL,
               col = NULL,
               pch = NULL,
               y.intersp = 1.5,
               ylim = NULL,
               inset = 0,
               ...) {
# Y.lb = lower bound of Y, such as Y - sd(Y)
# Y.ub = upepr bound of Y, such as Y + sd(Y)
    set.par()
    
    plot.int <- !is.null(Y.lb) && !is.null(Y.ub)
    if (is.null(ylim)) {
      if (plot.int) ylim <-  c(min(Y.lb), max(Y.ub))
      else ylim <- c(min(Y), max(Y))
      }
   else {
      miny <- min(Y)
      maxy <- max(Y)
      if (ylim[1] > miny) ylim[1] <- miny
      if (ylim[2] < maxy) ylim[2] <- maxy
      }
      

    if (is.null(lty)) lty <- alg.lty()[alg.ind]
    if (is.null(col)) col <- alg.cols()[alg.ind]
    if (is.null(pch)) pch <- alg.pch()[alg.ind]
    
    matplot(x, Y, t="o", bg = col, col=col, ylim = ylim, lty= lty, pch=pch, ...)
  
 
   # Legend
    if (plot.leg) {
      if (is.null(leg)) leg <- alg.names()[alg.ind]
      legend(leg.pos, leg, col=col, lty=lty, pch=pch, bg="WHITE", 
             y.intersp = y.intersp, pt.bg=col, inset = inset)
      }
      
    # Intervals
    if (plot.int) {
     for (i in 1:ncol(Y)) {
               plotCI(x, Y[,i], ui=Y.ub[,i], li = Y.lb[,i], add=TRUE,
                      col=col[i], lty=lty[i], pch=pch[i], pt.bg=par("bg"))
          }
      }

   }


plot.by.size <- function(filenames, sizes, 
      num.archs = NULL,
      steps.as.criterion = FALSE,
      plot.times = FALSE,
      alg.ind = 1, 
      mult.taus = TRUE,
      plot.int = FALSE,
      leg.pos = "topleft",
      ...){
# If num.archs = NULL, the routine will consider num.archs=size[i] for each i  
 
      stat.func <- which.min
      ylab <- "Steps"
      if (!steps.as.criterion) {
         stat.func <- which.max
         ylab <- "Return"
         }
      if (plot.times) ylab <- "Seconds"
   
     res  <- matrix(0, length(sizes), length(filenames))
     rmin <- matrix(0, length(sizes), length(filenames))
     rmax <- matrix(0, length(sizes), length(filenames))
     
     sx <- NULL # suffix
     if (!is.null(num.archs)) sx <- num.archs
     
     best.taus <- array(0,length(filenames))
     if(!mult.taus && !plot.times) {
          for (j in 1:length(filenames)) {
               T <- NULL
               for (i in 1:length(sizes)) {
#                     if (is.null(num.archs)) sx <- sizes[i]
                    T <- rbind(T,apply(read.table(paste(filenames[j], "_",
                             sizes[i], "_", sx, ".txt",sep="")),2,mean))
                    }
               best.taus[j] <- stat.func(apply(T,2,mean))
               }
          }
                    
     
     for (i in 1:length(sizes)) {
          for (j in 1:length(filenames)) {
              T <- NULL
           
              if (!plot.times) {
                  T <- read.table(paste(filenames[j] , "_",
                              sizes[i],"_", sx,".txt",sep=""))
                  }
              else T <- read.table(paste(filenames[j] , "_times_",
                              sizes[i],"_", sx,".txt",sep=""))
                              
              Tm <- apply(T,2,mean)

              if (mult.taus) k <- stat.func(Tm)
              else k <- best.taus[j]

              if (k==1) {
               print(paste("Warning: the value for algorithm", filenames[j],
                         "and num_archs", sizes[i],"is minimal."))
               }
              if(k==ncol(T)) {
               print(paste("Warning: the value for algorithm", filenames[j],
                         "and num_archs", sizes[i],"is maximal."))
               }
              if (!plot.times) res[i,j] <- Tm[k]
              else res[i,j] <- mean(Tm)
#                Confidence interval   
#                -------------------
               if (shapiro.test(T[,k])$p.value < 0.01) {
                    print(paste("Warning: the data for", filenames[j],
                  "column", k, "doesn't seem to be normally distributed."))
                    }
               alpha <- 0.99
               ci <- qnorm(1-(1-alpha)/2)*sd(T[,k])/sqrt(nrow(T))
               rmax[i,j] <- res[i,j] + ci
               rmin[i,j] <- res[i,j] - ci 
#                Maximum and minimum
#                -------------------
#               rmax[i,j] <- max(T[,k])
#               rmin[i,j] <- min(T[,k])

#              Std deviation
#              -------------------
#               rmax[i,j] <- res[i,j] + sd(T[,k])
#               rmin[i,j] <- res[i,j] - sd(T[,k])
              }
          }
    if (!plot.int) {
      rmax <- NULL
      rmin <- NULL
      }
    plot.alg(sizes, res, alg.ind, rmin, rmax, ylab=ylab, xlab="n",
               leg.pos=leg.pos, ...)
    }



plot.by.archs <- function(filenames, size, num.archs, 
      steps.as.criterion = FALSE,
      plot.times = FALSE,
      alg.ind = c(2,3), 
      mult.taus = TRUE,
      plot.int = FALSE,
      leg.pos = "topleft",
      inset = 0,
      stat.func = NULL,
      ylab = NULL,
      exclude.best.column = FALSE,
      ...){
   
      if (is.null(stat.func)) {
         stat.func <- which.max
         if (steps.as.criterion) stat.func <- which.min
        }

      if (is.null(ylab)) {
         ylab <- "Return"
         if (steps.as.criterion) ylab <- "Steps"
         }

      if (plot.times) ylab <- "Seconds"
      
     res  <- matrix(0, length(num.archs), length(filenames))
     rmin <- matrix(0, length(num.archs), length(filenames))
     rmax <- matrix(0, length(num.archs), length(filenames))
     
     best.taus <- array(0,length(filenames))
     if(!mult.taus && !plot.times) {
          for (j in 1:length(filenames)) {
               T <- NULL
               for (i in 1:length(num.archs)) {
                    T <- rbind(T,apply(read.table(paste(filenames[j], "_",
                             size[j], "_",
num.archs[i],".txt",sep="")),2,mean))
                    }
               Tm <- apply(T,2,mean)
               if (exclude.best.column) Tm <- Tm[-stat.func(Tm)]
               best.taus[j] <- stat.func()
               }
          }

     for (i in 1:length(num.archs)) {
          for (j in 1:length(filenames)) {
              T <- NULL
              if (!plot.times) {
                  T <- read.table(paste(filenames[j] , "_",
                              size[j],"_",
num.archs[i],".txt",sep=""))
                  }
               else T <- read.table(paste(filenames[j] , "_times_",
                              size[j],"_",
num.archs[i],".txt",sep=""))

              Tm <- apply(T,2,mean)
               if (exclude.best.column) Tm <- Tm[-stat.func(Tm)]

               if (mult.taus) k <- stat.func(Tm)
              else k <- best.taus[j]

              if (k==1) {
               print(paste("Warning: the value for algorithm", filenames[j],
                         "and num_archs", num.archs[i],"is minimal."))
               }
              if(k==ncol(T)) {
               print(paste("Warning: the value for algorithm", filenames[j],
                         "and num_archs", num.archs[i],"is maximal."))
               }
              if (!plot.times) res[i,j] <- Tm[k]
              else res[i,j] <- mean(Tm)

#                Alternative way
#                -------------------
#               Tm <- apply(T,1,max)
#               res[i,j] <- mean(Tm)

              
#                Confidence interval   
#                -------------------
# #                if (shapiro.test(T[,k])$p.value < 0.01) {
# #                     print(paste("Warning: the data for",
# filenames[j],
# #                   "column", k, "doesn't seem to be normally
# distributed."))
# #                     }
               alpha <- 0.99
               ci <- qnorm(1-(1-alpha)/2)*sd(T[,k])/sqrt(nrow(T))
               rmax[i,j] <- res[i,j] + ci
               rmin[i,j] <- res[i,j] - ci 
               
#                Maximum and minimum
#                -------------------
#               rmax[i,j] <- max(T[,k])
#               rmin[i,j] <- min(T[,k])

#              Std deviation
#              -------------------
#               rmax[i,j] <- res[i,j] + 0.5 * sd(T[,)k])
#               rmin[i,j] <- res[i,j] - 0.5 * sd(T[,k])

#              75% percentile
#              -------------------
#               Ts <- sort(T[,k], decreasing=TRUE)[1:round(0.5 * nrow(T))]
#               rmax[i,j] <- max(Ts)
#               rmin[i,j] <- min(Ts)
              }
          }
    if (!plot.int) {
      rmax <- NULL
      rmin <- NULL
      }
    plot.alg(num.archs, res, alg.ind, rmin, rmax, leg.pos=leg.pos,
             ylab=ylab, xlab = "m", inset=inset, ...)
    }


plot.by.df <- function(filename.lspi, filename.kbsf,
      dfs = c(0.7,0.75,0.8,0.85,0.9,0.95,0.99,0.999),
      size = 150000, 
      steps.as.criterion = FALSE,
      alg.ind = c(2,3), 
      plot.int = FALSE,
      leg.pos = "topleft",
      inset = 0,
      stat.func = NULL,
      ylab = NULL,
      ...) {
   
      if (is.null(stat.func)) {
         stat.func <- which.max
         if (steps.as.criterion) stat.func <- which.min
        }

      if (is.null(ylab)) {
         ylab <- "Return"
         if (steps.as.criterion) ylab <- "Steps"
         }
     
     res  <- matrix(0, length(dfs), 2)
     rmin <- matrix(0, length(dfs), 2)
     rmax <- matrix(0, length(dfs), 2)
     
     
     for (i in 1:length(dfs)) {
         # LSPI
         L <- read.table(paste(filename.lspi , "_", dfs[i], ".txt",sep=""))
         Lm <- apply(L,2,mean)
         l <- stat.func(Lm)
         res[i,1] <- Lm[l]
         rmax[i,1] <- res[i,1] + sd(L[l])
         rmin[i,1] <- res[i,1] - sd(L[l])

         # KBSF
         L <- read.table(paste(filename.kbsf , "_", dfs[i], ".txt",sep=""))
         Lm <- apply(L,2,mean)
         l <- stat.func(Lm)
         res[i,2] <- Lm[l]
         rmax[i,2] <- res[i,2] + sd(L[l])
         rmin[i,2] <- res[i,2] - sd(L[l])
         }
         
     print(rmax)
     print(rmin)
     if (!plot.int) {
      rmax <- NULL
      rmin <- NULL
      }
      
    plot.alg(dfs, res, alg.ind, rmin, rmax, leg.pos=leg.pos,
             ylab=ylab, xlab = expression(gamma), inset=inset, ...)
    }


concat.by.archs <- function(filename, size, num.archs) {
      T <- NULL
      for (i in 1:length(num.archs)) {
            T <- rbind(T, read.table(paste(filename, "_",
                     size, "_", num.archs[i], ".txt",sep="")))
            }
      T
      }


concat.by.sizes <- function(filename, sizes) {
      T <- NULL
      for (i in 1:length(sizes)) {
            T <- rbind(T, read.table(paste(filename, "_",
                     sizes[i], "_", sizes[i], ".txt",sep="")))
            }
      T
      }


concat.by.archs.cols <- function(filename, size, num.archs) {
      T <- read.table(paste(filename, "_",
                     size, "_", num.archs[1], ".txt",sep=""))
      for (i in 2:length(num.archs)) {
            T <- cbind(T, read.table(paste(filename, "_",
                     size, "_", num.archs[i], ".txt",sep="")))
            }
      T
      }


##----------------- PUDDLE WORLD ----------------------------------------------

puddle.plot.kbrl <- function(filenames.kbrl = c("./kbsf/puddle/uni/puddle_kbrl",
                                    "./kbsf/puddle/exp/puddle_kbrl"),
                       filename.grid = "./kbsf/puddle/puddle_grid_100.txt",
                       filename.eps = "~/tex/canada/puddle_kbrl.eps",
                       sizes = seq(200,2000,by=200),
                       ylim = c(-0.3, 1.2) # based on data
                     ) {
 
 plot.by.size(filenames.kbrl, sizes, NULL, FALSE, alg.ind = c(1,1), 
      mult.taus = TRUE, plot.int = FALSE, leg.pos = "bottomright", 
      lty= alg.lty(5)[c(1,5)],
      col=alg.cols(5)[c(1,5)], pch=alg.pch(5)[c(1,5)], 
      leg=c("KBRL(n) (generative model)", "KBRL(n) (exploration policy)"), 
      ylim = ylim)
  T <- read.table(filename.grid)
  mt <- mean(T)
  lines(sizes, rep(mt, length(sizes)), col= alg.cols(4)[4],
        lty=1, lwd=2)

  text(sizes[3], mt + 0.02 * (ylim[2]-ylim[1]), 
       "\"Optimal\" policy", col=alg.cols(4)[4])  
  dev.copy2eps(file= filename.eps)
  }

puddle.plot <- function(
                       filenames.kbrl=c("./kbsf/puddle/uni/puddle_kbrl" ,
                                    "./kbsf/puddle/exp/puddle_kbrl"),
                       filenames.lspi=c("./kbsf/puddle/uni/puddle_lspi" ,
                                    "./kbsf/puddle/exp/puddle_lspi"),
                       filenames.kbsf=c("./kbsf/puddle/uni/puddle_kbsf" ,
                                    "./kbsf/puddle/exp/puddle_kbsf"),
                       filename.grid = "./kbsf/puddle/puddle_grid_100.txt",
                       size = 2000,
                       num.archs = seq(10,150,by=10)
                       ) {
 

    ## RESULTS by m#
    #------------------------------#
    # Uniformly-sampled transitions 
     T <- read.table(paste(filenames.kbrl[1],"_",size,"_",size,".txt",sep=""))
     vmax <- max(apply(T, 2, mean))
     ylim <- c(-0.3, max(apply(T,2,mean)) + 0.1) # lower lim based on data
     
     # Without KBRL
     plot.by.archs(c(filenames.lspi[1], filenames.kbsf[1]), size, num.archs,
                   FALSE,leg.pos = "bottomright", ylim=ylim,
                   leg=c("LSPI(2000,m)", "KBSF(2000,m)"))
     lines(num.archs, rep(vmax, length(num.archs)), col= alg.cols(4)[4],
        lty=alg.lty(1)[1], lwd=2)
     text(num.archs[length(num.archs)-2], vmax+0.05, "KBRL(2000)")   
     dev.copy2eps(file="~/tex/canada/fig/puddle_lspi_kbsf_uni.eps")

     # With KBRL
     ylim <- c(-0.75, max(apply(T,2,mean)) + 0.1) # lower lim based on data
    plot.by.archs(c(filenames.kbrl[1], filenames.lspi[1], filenames.kbsf[1]),
                   size, num.archs, FALSE,leg.pos =
                   "bottomleft", ylim=ylim, alg.ind = 1:3,
                   leg=c("KBRL(m)", "LSPI(2000,m)", "KBSF(2000,m)"))
     lines(num.archs, rep(vmax, length(num.archs)), col= alg.cols(4)[4],
        lty=alg.lty(1)[1], lwd=2)
     text(num.archs[length(num.archs)-2], vmax+0.05, "KBRL(2000)")   
     dev.copy2eps(file="~/tex/canada/fig/puddle_kbrl_lspi_kbsf_uni.eps")


    # Transitions sampled by an exploration policy
     T <- read.table(paste(filenames.kbrl[2],"_",size,"_",size,".txt",sep=""))
     vmax <- max(apply(T, 2, mean))

     # Without KBRL
     ylim <- c(-0.3, max(apply(T,2,mean)) + 0.1) # lower lim based on data
     plot.by.archs(c(filenames.lspi[2], filenames.kbsf[2]), size, num.archs,
                   FALSE,leg.pos = "bottomright", ylim=ylim,
                   leg=c("LSPI(2000,m)", "KBSF(2000,m)"))
     lines(num.archs, rep(vmax, length(num.archs)), col= alg.cols(4)[4],
        lty=alg.lty(1)[1])
     text(num.archs[length(num.archs)-2], vmax+0.05, "KBRL(2000)")   
     dev.copy2eps(file="~/tex/canada/fig/puddle_lspi_kbsf_exp.eps")
     
     # With KBRL
     ylim <- c(-0.65, max(apply(T,2,mean)) + 0.1) # lower lim based on data
     plot.by.archs(c(filenames.kbrl[2], filenames.lspi[2], filenames.kbsf[2]),
                   size, num.archs, FALSE,leg.pos =
                   "bottomright", ylim=ylim, alg.ind = 1:3,
                   leg=c("KBRL(m)", "LSPI(2000,m)", "KBSF(2000,m)"))
     lines(num.archs, rep(vmax, length(num.archs)), col= alg.cols(4)[4],
        lty=alg.lty(1)[1], lwd=2)
     text(num.archs[length(num.archs)-2], vmax+0.05, "KBRL(2000)")   
     dev.copy2eps(file="~/tex/canada/fig/puddle_kbrl_lspi_kbsf_exp.eps")
   
   
    ## TIMES #
    #------------------------------#
    # Uniformly-sampled transitions 
    T <- as.matrix(read.table(paste(filenames.kbrl[1], "_times_", size, "_",
                     size,".txt",sep="")))
    mv <- rep(mean(T),length(num.archs))
    # Without KBRL
    ylim <- c(0, mv[1]+5)
    plot.by.archs(c(filenames.lspi[1],filenames.kbsf[1]), size, num.archs,
               plot.times=TRUE, ylim=ylim, leg.pos="left",
                    leg=c("LSPI(2000,m)", "KBSF(2000,m)"))
    lines(num.archs, mv, col= alg.cols(4)[4],lty=1,lwd=2)
    text(num.archs[length(num.archs)-2], mv-5, "KBRL(2000)", lwd=1)   
    dev.copy2eps(file="~/tex/canada/fig/puddle_lspi_kbsf_uni_times.eps")

   # With KBRL
   plot.by.archs(c(filenames.kbrl[1],filenames.lspi[1],filenames.kbsf[1]),
               size, num.archs, plot.times=TRUE, ylim=ylim,
               leg.pos="left", alg.ind=1:3,
                leg=c("KBRL(m)", "LSPI(2000,m)","KBSF(2000,m)"))
    lines(num.archs, mv, col= alg.cols(4)[4],lty=alg.lty(1)[1],lwd=2)   
    text(num.archs[length(num.archs)-2], mv-5, "KBRL(2000)", lwd=1)   
    dev.copy2eps(file="~/tex/canada/fig/puddle_kbrl_lspi_kbsf_uni_times.eps")
    
   # Transitions sampled by an exploration policy
    T <- as.matrix(read.table(paste(filenames.kbrl[2], "_times_", size, "_",
                     size,".txt",sep="")),)
    mv <- rep(mean(T),length(num.archs))
    
    # Without KBRL
    ylim <- c(0, mv[1]+5)
    plot.by.archs(c(filenames.lspi[2],filenames.kbsf[2]), size, num.archs,
               plot.times=TRUE, ylim=ylim, leg.pos="left",
               leg=c("LSPI(2000,m)", "KBSF(2000,m)"))
    lines(num.archs, mv, col= alg.cols(4)[4],lty=alg.lty(1)[1], lwd=2)
    text(num.archs[length(num.archs)-2], mv-5, "KBRL(2000)", lwd=1)   
    dev.copy2eps(file="~/tex/canada/fig/puddle_lspi_kbsf_exp_times.eps")

    # With KBRL
    plot.by.archs(c(filenames.kbrl[2],filenames.lspi[2],filenames.kbsf[2]),
               size, num.archs, plot.times=TRUE, ylim=ylim,
               leg.pos="left", alg.ind=1:3,
                leg=c("KBRL(m)", "LSPI(2000,m)","KBSF(2000,m)"))
    lines(num.archs, mv, col= alg.cols(4)[4],lty=alg.lty(1)[1],lwd=2)   
    text(num.archs[length(num.archs)-2], mv-5, "KBRL(2000)", lwd=1)   
    dev.copy2eps(file="~/tex/canada/fig/puddle_kbrl_lspi_kbsf_exp_times.eps")
    }
    


plot.puddle.lspi.kbsf.by.n70 <- function(
                       filenames.kbrl=c("./kbsf/puddle/uni/puddle_kbrl" ,
                                    "./kbsf/puddle/exp/puddle_kbrl"),
                       filenames.lspi=c("./kbsf/puddle/uni/puddle_lspi" ,
                                    "./kbsf/puddle/exp/puddle_lspi",
                                    "./kbsf/puddle/exp/puddle_lspi_70"),
                       filenames.kbsf=c("./kbsf/puddle/uni/puddle_kbsf" ,
                                    "./kbsf/puddle/exp/puddle_kbsf", 
                                    "./kbsf/puddle/exp/puddle_kbsf_70"),
                       filename.grid = "./kbsf/puddle/puddle_grid_100.txt",
                       size = 2000,
                       num.archs = seq(10,150,by=10)
                       ) {
     # N=70 Exploration policy
     T <- read.table(paste(filenames.kbrl[2],"_",size,"_",size,".txt",sep=""))
     ylim <- c(-0.3, max(apply(T,2,mean)) + 0.1) # lower lim based on data
     plot.by.size(filenames.kbrl[2],seq(400,2000,by=200),plot.leg=FALSE,
                  ylim=ylim)
     plot.by.size(c(filenames.lspi[3],filenames.kbsf[3]), seq(400,2000,by=200),
               num.archs=70, plot.leg=FALSE, alg.ind=c(2,3), add=TRUE)
      legend("bottomright", c("KBRL(n)", "LSPI(n,70)", "KBSF(n,70)"),
             col=alg.cols(3), lty=alg.lty(3), pch=alg.pch(3), bg="WHITE", 
             y.intersp = 1.5, pt.bg=alg.cols(3))
               
      dev.copy2eps(file="~/tex/canada/fig/puddle_kbrl_lspi_kbsf_exp_n70.eps")
      }
      

plot.puddle.lspi.kbsf.by.n <- function(
                       filenames.kbrl=c("./kbsf/puddle/uni/puddle_kbrl" ,
                                    "./kbsf/puddle/exp/puddle_kbrl"),
                       filenames.lspi=c("./kbsf/puddle/uni/puddle_lspi" ,
                                    "./kbsf/puddle/exp/puddle_lspi",
                                    "./kbsf/puddle/exp/puddle_lspi_70"),
                       filenames.kbsf=c("./kbsf/puddle/uni/puddle_kbsf" ,
                                    "./kbsf/puddle/exp/puddle_kbsf", 
                                    "./kbsf/puddle/exp/puddle_kbsf_70"),
                       filename.grid = "./kbsf/puddle/puddle_grid_100.txt",
                       size = 2000,
                       num.archs = seq(10,150,by=10)
                       ) {
 
     ## RESULTS by n#
    #------------------------------#
     T <- read.table(paste(filenames.kbrl[1],"_",size,"_",size,".txt",sep=""))
     vmax <- max(apply(T, 2, mean))

      # N=40 Uniformly-sampled transitions 
     ylim <- c(-0.1, max(apply(T,2,mean)) + 0.1) # lower lim based on data
     plot.by.size(filenames.kbrl[1],seq(200,2000,by=200),plot.leg=FALSE,
                  ylim=ylim)
     plot.by.size(c(filenames.lspi[1],filenames.kbsf[1]), seq(200,2000,by=200),
               num.archs=40, plot.leg=FALSE, alg.ind=c(2,3), add=TRUE)
      legend("bottomright", c("KBRL(n)", "LSPI(n,40)", "KBSF(n,40)"),
             col=alg.cols(3), lty=alg.lty(3), pch=alg.pch(3), bg="WHITE", 
             y.intersp = 1.5, pt.bg=alg.cols(3))
               
      dev.copy2eps(file="~/tex/canada/fig/puddle_kbrl_lspi_kbsf_uni_n40.eps")

     # N=40 Exploration policy
     ylim <- c(-0.1, max(apply(T,2,mean)) + 0.1) # lower lim based on data
     plot.by.size(filenames.kbrl[2],seq(200,2000,by=200),plot.leg=FALSE,
                  ylim=ylim)
     plot.by.size(c(filenames.lspi[2],filenames.kbsf[2]), seq(200,2000,by=200),
               num.archs=40, plot.leg=FALSE, alg.ind=c(2,3), add=TRUE)
      legend("bottomright", c("KBRL(n)", "LSPI(n,40)", "KBSF(n,40)"),
             col=alg.cols(3), lty=alg.lty(3), pch=alg.pch(3), bg="WHITE", 
             y.intersp = 1.5, pt.bg=alg.cols(3))
               
      dev.copy2eps(file="~/tex/canada/fig/puddle_kbrl_lspi_kbsf_exp_n40.eps")

     # N=100 Uniformly-sampled transitions 
     ylim <- c(-0.1, max(apply(T,2,mean)) + 0.1) # lower lim based on data
     plot.by.size(filenames.kbrl[1],seq(200,2000,by=200),plot.leg=FALSE,
                  ylim=ylim)
     plot.by.size(c(filenames.lspi[1],filenames.kbsf[1]), seq(200,2000,by=200),
               num.archs=100, plot.leg=FALSE, alg.ind=c(2,3), add=TRUE)
      legend("bottomright", c("KBRL(n)", "LSPI(n,100)", "KBSF(n,100)"),
             col=alg.cols(3), lty=alg.lty(3), pch=alg.pch(3), bg="WHITE", 
             y.intersp = 1.5, pt.bg=alg.cols(3))
               
      dev.copy2eps(file="~/tex/canada/fig/puddle_kbrl_lspi_kbsf_uni_n100.eps")

     # N=100 Exploration policy
     ylim <- c(-0.3, 0.8) # lower lim based on data
     plot.by.size(filenames.kbrl[2],seq(200,2000,by=200),plot.leg=FALSE,
                  ylim=ylim, col = alg.cols(4)[4])
     plot.by.size(c(filenames.lspi[2],filenames.kbsf[2]), seq(200,2000,by=200),
               num.archs=100, plot.leg=FALSE, alg.ind=c(2,3), add=TRUE)
      legend("bottomright", c("KBRL(n)", "LSPI(n,100)", "KBSF(n,100)"),
             col=alg.cols(4)[c(4,2,3)], lty=alg.lty(3), pch=alg.pch(3), bg="WHITE", 
             y.intersp = 1.5, pt.bg=alg.cols(4)[c(4,2,3)])
               
      dev.copy2eps(file="~/tex/canada/fig/puddle_kbrl_lspi_kbsf_exp_n100.eps")



## TIMES by n#
    #------------------------------#
     # N=40 Uniformly-sampled transitions 
     ylim <- c(-0.1, max(apply(T,2,mean)) + 0.1) # lower lim based on data
     plot.by.size(filenames.kbrl[1],seq(200,2000,by=200),plot.leg=FALSE,
                  ylim=ylim, plot.times=TRUE)
     plot.by.size(c(filenames.lspi[1],filenames.kbsf[1]), seq(200,2000,by=200),
               num.archs=40, plot.leg=FALSE, alg.ind=c(2,3), add=TRUE,
               plot.times=TRUE)
      legend("topleft", c("KBRL(n)", "LSPI(n,40)", "KBSF(n,40)"),
             col=alg.cols(3), lty=alg.lty(3), pch=alg.pch(3), bg="WHITE", 
             y.intersp = 1.5, pt.bg=alg.cols(3))
               
      dev.copy2eps(file="~/tex/canada/fig/puddle_kbrl_lspi_kbsf_uni_n40_times.eps")

     # N=40 Exploration Policy
     ylim <- c(-0.1, max(apply(T,2,mean)) + 0.1) # lower lim based on data
     plot.by.size(filenames.kbrl[1],seq(200,2000,by=200),plot.leg=FALSE,
                  ylim=ylim, plot.times=TRUE)
     plot.by.size(c(filenames.lspi[1],filenames.kbsf[1]), seq(200,2000,by=200),
               num.archs=40, plot.leg=FALSE, alg.ind=c(2,3), add=TRUE,
               plot.times=TRUE)
      legend("topleft", c("KBRL(n)", "LSPI(n,40)", "KBSF(n,40)"),
             col=alg.cols(3), lty=alg.lty(3), pch=alg.pch(3), bg="WHITE", 
             y.intersp = 1.5, pt.bg=alg.cols(3))
               
      dev.copy2eps(file="~/tex/canada/fig/puddle_kbrl_lspi_kbsf_exp_n40_times.eps")

        # N=100 Uniformly-sampled transitions 
     ylim <- c(-0.1, max(apply(T,2,mean)) + 0.1) # lower lim based on data
     plot.by.size(filenames.kbrl[1],seq(200,2000,by=200),plot.leg=FALSE,
                  ylim=ylim, plot.times=TRUE)
     plot.by.size(c(filenames.lspi[1],filenames.kbsf[1]), seq(200,2000,by=200),
               num.archs=100, plot.leg=FALSE, alg.ind=c(2,3), add=TRUE,
               plot.times=TRUE)
      legend("topleft", c("KBRL(n)", "LSPI(n,100)", "KBSF(n,100)"),
             col=alg.cols(3), lty=alg.lty(3), pch=alg.pch(3), bg="WHITE", 
             y.intersp = 1.5, pt.bg=alg.cols(3))
               
      dev.copy2eps(file="~/tex/canada/fig/puddle_kbrl_lspi_kbsf_uni_n100_times.eps")

     # N=100 Exploration Policy
     ylim <- c(-0.1, max(apply(T,2,mean)) + 0.1) # lower lim based on data
     plot.by.size(filenames.kbrl[1],seq(200,2000,by=200),plot.leg=FALSE,
                  ylim=ylim, plot.times=TRUE, col=alg.cols(4)[4])
     plot.by.size(c(filenames.lspi[1],filenames.kbsf[1]), seq(200,2000,by=200),
               num.archs=100, plot.leg=FALSE, alg.ind=c(2,3), add=TRUE,
               plot.times=TRUE)
      legend("topleft", c("KBRL(n)", "LSPI(n,100)", "KBSF(n,100)"),
             col=alg.cols(4)[c(4,2,3)], lty=alg.lty(3), pch=alg.pch(3), bg="WHITE", 
             y.intersp = 1.5, pt.bg=alg.cols(4)[c(4,2,3)])
               
      dev.copy2eps(file="~/tex/canada/fig/puddle_kbrl_lspi_kbsf_exp_n100_times.eps")
      }
 

##----------------- MOUNTAIN CAR---------------------------------------------

mountain.plot.kbrl <- function(
                     filenames.kbrl = c("./kbsf/mountain/uni/mountain_kbrl",
                                 "./kbsf/mountain/exp/mountain_kbrl"),
                     filenames.grid =c("./kbsf/mountain/mountain_grid_250.txt",
                                 "./kbsf/mountain/mountain_grid_steps_250.txt"),
                       sizes = seq(200,2000,by=200),
                       ylim = c(0.07, 0.55) # based on data
                     ) {
 
 plot.by.size(filenames.kbrl, sizes, NULL, FALSE, alg.ind = c(1,1), 
      mult.taus = TRUE, plot.int = FALSE, leg.pos = "bottomright", 
      lty= alg.lty(5)[c(1,5)],
      col=alg.cols(5)[c(1,5)], pch=alg.pch(5)[c(1,5)], 
      leg=c("KBRL(n) (generative model)", "KBRL(n) (exploration policy)"), 
      ylim = ylim)
  T <- read.table(filenames.grid[1])
  mt <- mean(T)
  lines(sizes, rep(mt, length(sizes)), col= alg.cols(4)[4],
        lty=1, lwd=2)

  text(sizes[3], mt + 0.02 * (ylim[2]-ylim[1]), "\"Optimal\" policy")  
  dev.copy2eps(file= "~/tex/canada/mountain_kbrl.eps")
 

  T <- read.table(filenames.grid[2])
  mt <- mean(T)
  ylim <- c(mt-0.1*(mt), mt)
 plot.by.size(paste(filenames.kbrl,"_steps",sep=""), sizes, NULL, FALSE, 
      alg.ind = c(1,1), 
      mult.taus = TRUE, plot.int = FALSE, leg.pos = "topright", 
      lty= alg.lty(5)[c(1,5)],
      col=alg.cols(5)[c(1,5)], pch=alg.pch(5)[c(1,5)], 
      leg=c("KBRL(n) (generative model)", "KBRL(n) (exploration policy)"), 
      steps.as.criterion = TRUE, ylim=ylim)
  lines(sizes, rep(mt, length(sizes)), col= alg.cols(4)[4],
        lty=1, lwd=2)

  text(sizes[3], mt - 8,  "\"Optimal\" policy")  
  dev.copy2eps(file= "~/tex/canada/mountain_kbrl_steps.eps")
  }

plot.mountain  <- function(
                       filenames.kbrl=c("./kbsf/mountain/uni/mountain_kbrl" ,
                                    "./kbsf/mountain/exp/mountain_kbrl"),
                       filenames.lspi=c("./kbsf/mountain/uni/mountain_lspi" ,
                                    "./kbsf/mountain/exp/mountain_lspi"),
                       filenames.kbsf=c("./kbsf/mountain/uni/mountain_kbsf" ,
                                    "./kbsf/mountain/exp/mountain_kbsf"),
                       filename.grid = "./kbsf/mountain/mountain_grid_250.txt",
                       size = 2000,
                       num.archs = seq(10,150,by=10)
                       ) {
 
#     ## RETURN #
#     #------------------------------#
    # Uniformly-sampled transitions 
     T <- read.table(paste(filenames.kbrl[1],"_",size,"_",size,".txt",sep=""))
     vmax <- max(apply(T, 2, mean))
     ylim <- c(vmax, vmax + 0.02) # plot.by.archs will correct
     
     plot.by.archs(c(filenames.kbrl[1], filenames.lspi[1], filenames.kbsf[1]),
                   size, num.archs,FALSE, leg.pos ="bottomright", ylim=ylim,
                   leg=c("KBRL(m)", "LSPI(2000,m)", "KBSF(2000,m)"), 
                   alg.ind = 1:3)
     lines(num.archs, rep(vmax, length(num.archs)), col= alg.cols(4)[4],
        lty=alg.lty(1)[1], lwd=2)
      text(num.archs[length(num.archs)-2], vmax+0.01, "KBRL(2000)")   
      dev.copy2eps(file="~/tex/canada/fig/mountain_uni.eps")
      
    # Transitions sampled by an exploration policy
     T <- read.table(paste(filenames.kbrl[2],"_",size,"_",size,".txt",sep=""))
     vmax <- max(apply(T, 2, mean))
     ylim <- c(vmax, vmax + 0.02) 
     plot.by.archs(c(filenames.kbrl[2], filenames.lspi[2], filenames.kbsf[2]),
                   size, num.archs,
                   FALSE,leg.pos = "bottomright", ylim=ylim,
                   leg=c("KBRL(m)", "LSPI(2000,m)", "KBSF(2000,m)"),
                   alg.ind = 1:3)
     lines(num.archs, rep(vmax, length(num.archs)), col= alg.cols(4)[4],
        lty=alg.lty(1)[1])
     text(num.archs[length(num.archs)-2], vmax+0.01, "KBRL(2000)")   
     dev.copy2eps(file="~/tex/canada/fig/mountain_exp.eps")

#     ## STEPS #
#     #------------------------------#
    # Uniformly-sampled transitions 
     T <-read.table(paste(filenames.kbrl[1],"_steps_",
                     size,"_",size,".txt",sep=""))
     vmin <- min(apply(T, 2, mean))
     ylim <- c(vmin-5, vmin) # plot.by.archs will correct
     
     plot.by.archs(paste(c(filenames.kbrl[1], filenames.lspi[1],
                          filenames.kbsf[1]), "_steps", sep=""),
                   size, num.archs,FALSE, leg.pos ="topright", ylim=ylim,
                   leg=c("KBRL(m)", "LSPI(2000,m)", "KBSF(2000,m)"), 
                   alg.ind = 1:3, steps.as.criterion=TRUE)
     lines(num.archs, rep(vmin, length(num.archs)), col= alg.cols(4)[4],
        lty=alg.lty(1)[1], lwd=2)
      text(num.archs[length(num.archs)-2], vmin + 5, "KBRL(2000)")   
      dev.copy2eps(file="~/tex/canada/fig/mountain_steps_uni.eps")

    # Exploration policy
     T <-read.table(paste(filenames.kbrl[2],"_steps_",
                     size,"_",size,".txt",sep=""))
     vmin <- min(apply(T, 2, mean))
     ylim <- c(vmin-5, vmin) # plot.by.archs will correct
     
     plot.by.archs(paste(c(filenames.kbrl[2], filenames.lspi[2],
                          filenames.kbsf[2]), "_steps", sep=""),
                   size, num.archs,FALSE, leg.pos ="topright", ylim=ylim,
                   leg=c("KBRL(m)", "LSPI(2000,m)", "KBSF(2000,m)"), 
                   alg.ind = 1:3, steps.as.criterion=TRUE)
     lines(num.archs, rep(vmin, length(num.archs)), col= alg.cols(4)[4],
        lty=alg.lty(1)[1], lwd=2)
      text(num.archs[length(num.archs)-2], vmin - 7, "KBRL(2000)")   
      dev.copy2eps(file="~/tex/canada/fig/mountain_steps_exp.eps")

#     ## TIMES #
#     #------------------------------#
#     # Uniformly-sampled transitions 
    T <- as.matrix(read.table(paste(filenames.kbrl[1], "_times_", size, "_",
                     size,".txt",sep="")))
    mv <- rep(mean(T),length(num.archs))
    # Without KBRL
    ylim <- c(0, mv[1]+1)
    plot.by.archs(c(filenames.kbrl[1], filenames.lspi[1],filenames.kbsf[1]),
               size, num.archs,plot.times=TRUE, ylim=ylim,
               leg.pos="left",leg=c("KBRL(m)", "LSPI(2000,m)","KBSF(2000,m)"),
               alg.ind=1:3)
    lines(num.archs, mv, col= alg.cols(4)[4],lty=1,lwd=2)
    text(num.archs[length(num.archs)-2], mv-3, "KBRL(2000)", lwd=1)   
    dev.copy2eps(file="~/tex/canada/fig/mountain_uni_times.eps")

   # Transitions sampled by an exploration policy
    T <- as.matrix(read.table(paste(filenames.kbrl[2], "_times_", size, "_",
                     size,".txt",sep="")),)
    mv <- rep(mean(T),length(num.archs))
    
    # Without KBRL
    ylim <- c(0, mv[1]+1)
    plot.by.archs(c(filenames.kbrl[2],
               filenames.lspi[2],filenames.kbsf[2]), size, num.archs,
               plot.times=TRUE, ylim=ylim, leg.pos="left",
               leg=c("KBRL(m)", "LSPI(2000,m)", "KBSF(2000,m)"), alg.ind=1:3)
    lines(num.archs, mv, col= alg.cols(4)[4],lty=alg.lty(1)[1], lwd=2)
    text(num.archs[3], mv-2, "KBRL(2000)", lwd=1)   
    dev.copy2eps(file="~/tex/canada/fig/mountain_exp_times.eps")

    ## SUCCESS by m#
    #------------------------------#
    # Uniformly-sampled transitions 
     T <- read.table(paste(filenames.kbrl[1],"_",size,"_",size,".txt",sep=""))
     vmin <- rep(min(apply(T != 0, 2, sum) / nrow(T)), length(num.archs))
     ylim <- c(0.5, 1.1) # lower lim based on data
     
     res.kbrl <- array(0, length(num.archs))
     res.lspi <- array(0, length(num.archs))
     res.kbsf <- array(0, length(num.archs))
     for (na in 1:length(num.archs)) {
         T <- read.table(paste(filenames.kbrl[1], "_", size, "_", num.archs[na],
                         ".txt", sep=""))
         res.kbrl[na] <- min(apply(T != 0, 2, sum) / nrow(T))

         T <- read.table(paste(filenames.lspi[1], "_", size, "_", num.archs[na],
                         ".txt", sep=""))
         res.lspi[na] <- min(apply(T != 0, 2, sum) / nrow(T))

         T <- read.table(paste(filenames.kbsf[1], "_", size, "_", num.archs[na],
                         ".txt", sep=""))
         res.kbsf[na] <- min(apply(T != 0, 2, sum) / nrow(T))
         }
      Y <- cbind(res.kbrl, res.lspi, res.kbsf)

      plot.alg(num.archs, Y, 1:3, plot.leg = TRUE, 
               leg = c("KBRL(m)", "LSPI(2000,m)", "KBSF(2000,m)") ,
               leg.pos = "right", ylim = ylim, xlab="m", ylab="Success rate")
         
     lines(num.archs, vmin, col= alg.cols(4)[4],lty=alg.lty(1)[1], lwd=2)
     text(num.archs[3], vmin[1]+0.03, "KBRL(2000)")   
      dev.copy2eps(file="~/tex/canada/fig/mountain_uni_success.eps")

    #------------------------------#
    #Exploration policy
     T <- read.table(paste(filenames.kbrl[2],"_",size,"_",size,".txt",sep=""))
     vmin <- rep(min(apply(T != 0, 2, sum) / nrow(T)), length(num.archs))
     ylim <- c(0, vmin[1] + 0.2) # lower lim based on data
     
     res.kbrl <- array(0, length(num.archs))
     res.lspi <- array(0, length(num.archs))
     res.kbsf <- array(0, length(num.archs))
     for (na in 1:length(num.archs)) {
         T <- read.table(paste(filenames.kbrl[2], "_", size, "_", num.archs[na],
                         ".txt", sep=""))
         res.kbrl[na] <- min(apply(T != 0, 2, sum) / nrow(T))

         T <- read.table(paste(filenames.lspi[2], "_", size, "_", num.archs[na],
                         ".txt", sep=""))
         res.lspi[na] <- min(apply(T != 0, 2, sum) / nrow(T))

         T <- read.table(paste(filenames.kbsf[2], "_", size, "_", num.archs[na],
                         ".txt", sep=""))
         res.kbsf[na] <- min(apply(T != 0, 2, sum) / nrow(T))
         }
      Y <- cbind(res.kbrl, res.lspi, res.kbsf)

      plot.alg(num.archs, Y, 1:3, plot.leg = TRUE, 
               leg = c("KBRL(m)", "LSPI(2000,m)", "KBSF(2000,m)"),
               leg.pos= "topleft",
               ylim = ylim, xlab="m",ylab="Success rate")
         
     lines(num.archs, vmin, col= alg.cols(4)[4],lty=alg.lty(1)[1], lwd=2)
     text(num.archs[length(num.archs)-2], vmin[1]+0.03, "KBRL(2000)")   
     
      dev.copy2eps(file="~/tex/canada/fig/mountain_exp_success.eps")

   } 
    
                     
##----------------- POLE BALANCING---------------------------------------------


plot.pole <- function(
                       filenames.lspi=c("./kbsf/pole/uni/pole_lspi" ,
                                    "./kbsf/pole/exp/pole_lspi"),
                       filenames.kbsf=c("./kbsf/pole/uni/pole_kbsf" ,
                                    "./kbsf/pole/exp/pole_kbsf"),
                       filename.grid = "./kbsf/pole/pole_grid_250.txt",
                       sizes = c(100000,150000),
                       num.archs = seq(10,100,by=10)
                       ) {
   
   # LSPI
   #    size[1] = 100000
   plot.by.archs(paste(c(filenames.lspi[1], filenames.lspi[2]), "_steps",
                  sep=""), sizes[1], seq(10,100,by=10),
                  stat.func=which.max, leg.pos="topleft", ylab="Steps",
                  alg.ind=c(2,4), leg=c("LSPI (generative model)",
                   "LSPI (exploration policy)"))
   dev.copy2eps(file="~/tex/canada/fig/pole_lspi_steps.eps")   


   #    size[2] = 150000
   plot.by.archs(paste(c(filenames.lspi[1], filenames.lspi[2]), "_steps",
                  sep=""), sizes[2], seq(10,100,by=10),
                  stat.func=which.max, leg.pos="topleft", ylab="Steps",
                  alg.ind=c(2,4), leg=c("LSPI (generative model)",
                   "LSPI (exploration policy)"), inset=c(0,0.6))
   dev.copy2eps(file="~/tex/canada/fig/pole_lspi_steps_150.eps")   

   # KBSF
   #    size[1] = 100000
   plot.by.archs(paste(c(filenames.kbsf[1], filenames.kbsf[2]), "_steps",
                  sep=""), sizes[1], seq(10,100,by=10),
                  stat.func=which.max, leg.pos="bottomright", ylab="Steps",
                  alg.ind=c(2,4), leg=c("KBSF (generative model)",
                   "KBSF (exploration policy)"))
   dev.copy2eps(file="~/tex/canada/fig/pole_kbsf_steps.eps")   



   #    size[2] = 150000
   plot.by.archs(paste(c(filenames.kbsf[1], filenames.kbsf[2]), "_steps",
                  sep=""), sizes[2], seq(10,100,by=10),
                  stat.func=which.max, leg.pos="topright", ylab="Steps",
                  alg.ind=c(2,4), leg=c("KBSF (generative model)",
                   "KBSF (exploration policy)"), inset=c(0,0.3))
   dev.copy2eps(file="~/tex/canada/fig/pole_kbsf_steps_150.eps")   
   

# Uniformly-sampled transitions 
#   Steps
#    size[1] = 100000
     plot.by.archs(paste(c(filenames.lspi[1], filenames.kbsf[1]), "_steps",
                  sep=""), sizes[1], num.archs, leg.pos ="topleft",
                  leg=c(expression("LSPI("*10^5*",m)"),
                  expression("KBSF("*10^5*",m)")), steps.as.criterion=TRUE,
                  stat.func = which.max, ylim=c(1500,1500))
      dev.copy2eps(file="~/tex/canada/fig/pole_uni_steps.eps")
      

#    size[2] = 150000
       plot.by.archs(paste(c(filenames.lspi[1], filenames.kbsf[1]), "_steps",
                  sep=""), sizes[2], num.archs, leg.pos ="topleft",
                  leg=c(expression("LSPI(1.5 x"*10^5*",m)"),
                  expression("KBSF(1.5 x"*10^5*",m)")), steps.as.criterion=TRUE,
                  stat.func = which.max)
      dev.copy2eps(file="~/tex/canada/fig/pole_uni_steps_150.eps")
      
      #   Episodes
#    size[1] = 100000
     plot.by.archs(paste(c(filenames.lspi[1], filenames.kbsf[1]), "_episodes",
                  sep=""), sizes[1], num.archs, leg.pos ="topleft",
                  leg=c(expression("LSPI("*10^5*",m)"),
                  expression("KBSF("*10^5*",m)")), steps.as.criterion=TRUE,
                  stat.func = which.max, ylab="Successful episodes",
                  ylim=c(0.45,0.45))
      dev.copy2eps(file="~/tex/canada/fig/pole_uni_ep.eps")


      #    size[2] = 150000
       plot.by.archs(paste(c(filenames.lspi[1], filenames.kbsf[1]), "_episodes",
                  sep=""), sizes[2], num.archs, leg.pos ="topleft",
                  leg=c(expression("LSPI(1.5 x"*10^5*",m)"),
                  expression("KBSF(1.5 x"*10^5*",m)")), steps.as.criterion=TRUE,
                  stat.func = which.max, ylab="Successful episodes")
      dev.copy2eps(file="~/tex/canada/fig/pole_uni_ep_150.eps")
      
 #    Exploration policy
#    Steps
#    size[1] = 100000
     plot.by.archs(paste(c(filenames.lspi[2], filenames.kbsf[2]), "_steps",
                  sep=""), sizes[1], num.archs, leg.pos ="bottomleft",
                  leg=c(expression("LSPI("*10^5*",m)"),
                  expression("KBSF("*10^5*",m)")), steps.as.criterion=TRUE,
                  stat.func = which.max)
      dev.copy2eps(file="~/tex/canada/fig/pole_exp_steps.eps")

#    size[2] = 150000
       plot.by.archs(paste(c(filenames.lspi[2], filenames.kbsf[2]), "_steps",
                  sep=""), sizes[2], num.archs, leg.pos ="bottomleft",
                  leg=c(expression("LSPI(1.5 x"*10^5*",m)"),
                  expression("KBSF(1.5 x"*10^5*",m)")), steps.as.criterion=TRUE,
                  stat.func = which.max)
      dev.copy2eps(file="~/tex/canada/fig/pole_exp_steps_150.eps")

#   Episodes
#    size[1] = 100000
     plot.by.archs(paste(c(filenames.lspi[2], filenames.kbsf[2]), "_episodes",
                  sep=""), sizes[1], num.archs, leg.pos ="bottomleft",
                  leg=c(expression("LSPI("*10^5*",m)"),
                  expression("KBSF("*10^5*",m)")), steps.as.criterion=TRUE,
                  stat.func = which.max, ylab="Successful episodes")
      dev.copy2eps(file="~/tex/canada/fig/pole_exp_ep.eps")

#    size[2] = 150000
       plot.by.archs(paste(c(filenames.lspi[2], filenames.kbsf[2]), "_episodes",
                  sep=""), sizes[2], num.archs, leg.pos ="bottomleft",
                  leg=c(expression("LSPI(1.5 x"*10^5*",m)"),
                  expression("KBSF(1.5 x"*10^5*",m)")), steps.as.criterion=TRUE,
                  stat.func = which.max, ylab="Successful episodes")
      dev.copy2eps(file="~/tex/canada/fig/pole_exp_ep_150.eps")

     }    


lf <- function(x) x <= 1000
 
plot.pole.test <- function() {
   L <- concat.by.archs.cols("./kbsf/pole/exp/pole_lspi_steps", size=100000,
         seq(10,100,by=10))
   K <- concat.by.archs.cols("./kbsf/pole/exp/pole_kbsf_steps", size=100000,
         seq(10,100,by=10))

   Lz <- apply(apply(L,2,lf), 2, mean)
   Kz <- apply(apply(K,2,lf), 2, mean)
   
   plot.alg(seq(10,100,by=10), cbind(Lz,Kz), leg=c("LSPI", "KBSF"),
               leg.pos="topleft", xlab="m", ylab="Fraction <= 1000",
               lty=alg.lty()[c(2,3)], col=alg.cols()[c(2,3)] )

   dev.copy2eps(file = "~/tex/canada/fig/pole_failures_1000.eps")
   
   for (i in 1:length(Lz)) {
      Lz[i] <- mean(L[!lf(L[,i]),i])
      Kz[i] <- mean(K[!lf(K[,i]),i])
      }
      
   plot.alg(seq(10,100,by=10), cbind(Lz,Kz), leg=c("LSPI", "KBSF"),
               leg.pos="left", xlab="m", ylab="Mean (> 1000)",
               lty=alg.lty()[c(2,3)], col=alg.cols()[c(2,3)] )

   dev.copy2eps(file = "~/tex/canada/fig/pole_mean_suc_1000.eps")

   L <- concat.by.archs.cols("./kbsf/pole/exp/pole_lspi_steps", size=150000,
         seq(10,100,by=10))
   K <- concat.by.archs.cols("./kbsf/pole/exp/pole_kbsf_steps", size=150000,
         seq(10,100,by=10))

   Lz <- apply(apply(L,2,lf), 2, mean)
   Kz <- apply(apply(K,2,lf), 2, mean)
   
   plot.alg(seq(10,100,by=10), cbind(Lz,Kz), leg=c("LSPI", "KBSF"),
               leg.pos="topleft", xlab="m", ylab="Fraction <= 1000",
               lty=alg.lty()[c(2,3)], col=alg.cols()[c(2,3)] )

   dev.copy2eps(file = "~/tex/canada/fig/pole_failures_150_1000.eps")

   for (i in 1:length(Lz)) {
      Lz[i] <- mean(L[!lf(L[,i]),i])
      Kz[i] <- mean(K[!lf(K[,i]),i])
      }
      
   plot.alg(seq(10,100,by=10), cbind(Lz,Kz), leg=c("LSPI", "KBSF"),
               leg.pos="bottomright", xlab="m", ylab="Mean (> 1000)",
               lty=alg.lty()[c(2,3)], col=alg.cols()[c(2,3)] )

   dev.copy2eps(file = "~/tex/canada/fig/pole_mean_suc_150_1000.eps")

   } 
    


print("kbsf.data.analysis.R loaded")  