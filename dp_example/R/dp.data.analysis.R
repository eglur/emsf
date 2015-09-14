# library("psych")

# Pole n = 50.000, m = seq(10,100,by=10)
# Double Pole n = 1000000, m = seq(20, 200, by=20)
# Puddle   n = 8000, m = seq(10, 150, by = 20)
#          n = seq(1000, 10000, by = 1000), m = 100
# Mountain n = 5000, m = seq(10, 150, by = 20)
#          n = seq(1000, 10000, by = 1000), m = 100


library(plotrix)
source("util.R") 

set.par <- function() {
# set parameters used  in the plots
   par(
      lwd = 1.5,
      mar = c(5, 6, 2, 2) -1 + 0.1,
      ps  = 20)
   }

num.algs <- 5

alg.names <- function(n = num.algs) {
   names <- c("KBRL", "LSPI", "KBSF", "TREE")
   names[1:n]
   }

alg.lty <- function(n = num.algs) {
   types <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
   types[1:n]
   }


alg.cols <- function(n = num.algs, alpha = 1) {
   col <- c(rgb(0,0,0, alpha), rgb(1,0,0, alpha), rgb(0,0,1, alpha), rgb(0,1,0,alpha), rgb(0,1,1,alpha), 
      rgb(1,0,1,alpha), rgb(1,1,0,alpha))
   col[1:n]
   }

alg.pch <- function(n = num.algs) {
   seq(21, 21 + n -1, by = 1)
   }

plot.algs  <- function(x, Y, 
               alg.ind = NULL,
               Y.lb = NULL,
               Y.ub = NULL,
               plot.leg = FALSE,
               leg = NULL,
               leg.pos = "topright",
               lty = NULL,
               col = NULL,
               col.ci = NULL,
               col.sh = NULL,
               pch = NULL,
               t = "o",
               y.intersp = 1.5,
               ylim = NULL,
               inset = 0,
               show.int = NULL,
               show.shadow = NULL,
               new = TRUE,
               setpar = TRUE,
               ...) {
# Y.lb = lower bound of Y, such as Y - sd(Y)
# Y.ub = upper bound of Y, such as Y + sd(Y)

    
    Y    <- as.matrix(Y)
    plot.int <- !is.null(Y.lb) && !is.null(Y.ub)
    if (plot.int) {
      Y.lb <- as.matrix(Y.lb)
      Y.ub <- as.matrix(Y.ub)
      }

   if (is.null(show.int)) show.int <- rep(TRUE, ncol(Y))
   if (is.null(show.shadow)) show.shadow <- rep(TRUE, ncol(Y))
   
    if (is.null(alg.ind)) alg.ind <- 1:ncol(Y)
    
    if (setpar) set.par()
    
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
    if (is.null(col.ci)) col.ci <- alg.cols(alpha = 0.2)[alg.ind]
    if (is.null(col.sh)) col.sh <- rep(alg.cols(1, alpha = 0.1), length(alg.ind))
    if (is.null(pch)) pch <- alg.pch()[alg.ind]

    if (new) matplot(x, Y, t=t, bg = col, col=col, ylim = ylim, lty= lty, pch=pch, ...)
    else matlines(x, Y, t=t, bg = col, col=col, ylim = ylim, lty= lty, pch=pch, ...)
  
   # Legend
    if (plot.leg) {
      if (is.null(leg)) leg <- alg.names()[alg.ind]
      legend(leg.pos, leg, col=col, lty=lty, pch=pch, bg="WHITE", 
             y.intersp = y.intersp, pt.bg=col, inset = inset)
      }
      
    # Intervals
    if (plot.int) {
     for (i in 1:ncol(Y)) {
                if (show.int[i]) {
                  plotCI(x, Y[,i], ui=Y.ub[,i], li = Y.lb[,i], add=TRUE,
                        col=col[i], scol= col.ci[i], lty=lty[i], pch=pch[i], 
                        pt.bg=col[i], sfrac=0)
                  }
                  
               if (show.shadow[i]) {
                  px <- c(x,rev(x))
                  py <- c(Y.lb[,i], rev(Y.ub[,i]))
                  polygon(px,py, col = col.sh[1], border = NA)
                  }
          }
      }
      
    matlines(x, Y, t= t, bg = col, col=col, ylim = ylim, lty= lty, pch=pch,  ...)
   }


load.data <- function(filename, size, taus = NULL, m = NULL, taus2 = NULL, col=1, sf=".res", nm = NULL) {

   m.sx <- ""
   if (!is.null(m)) m.sx <- ps("_",m)
   
   ind.max1 <- 1
   if (!is.null(taus)) ind.max1 <- length(taus)

   ind.max2 <- 1
   if (!is.null(taus2)) ind.max2 <- length(taus2)
   
   nm.sx <- ""
   if (!is.null(nm)) nm.sx <- ps("_",nm)
   
   T <- NULL

   for (t in 1:ind.max1) {

      t.sx <- ""
      if (!is.null(taus)) { 
         t.sx <- ps("_",taus[t])
         if (taus[t] > 1e-5 && taus[t] <= 1e-4) t.sx <- ps("_", sprintf("%5.4f", taus[t])) # this is because C prints 1e-4 as "0.0004"
         }
      
      for (t2 in 1:ind.max2) {
         t.sx2 <- ""
         if (!is.null(taus2)) {
            t.sx2 <- ps("_",taus2[t2])
            if (taus2[t2] > 1e-5 && taus2[t2] <= 1e-4) t.sx <- ps("_", sprintf("%5.4f", taus2[t2])) # this is because C prints 1e-4  as          
            }
         s.sx <- size
         if (size >= 1e5 && size < 1e6) s.sx <- sprintf("%5.0f", size)
         else if (size >= 1e6 && size < 1e7) s.sx <- sprintf("%6.0f", size)
         fn <- ps(filename, "_", s.sx, nm.sx, m.sx, t.sx, t.sx2, sf)

         if (is.null(T)) {
            T <- read.matrix(fn)[,col]
            T <- matrix(T, length(T),1)
            }
         else T <- cbind(T,read.table(fn)[,col])
         colnames(T)[ncol(T)] <- paste(taus[t], t.sx2)
      }
    }
  T
  }
    

id <- function(x) x

gen.data <- function(filename, sizes, taus = NULL, ms = NULL, taus2 = NULL, col=1, pick.func = which.max, 
                     sf=".res", check.normallity = FALSE, confidence.level = 0.99,
                     transf.func = id, nm = NULL) {
   
   ind.max <- 1
   if (!is.null(ms)) ind.max <- length(ms)

   res  <- array(0, length(sizes) * ind.max)
   rmax <- array(0, length(sizes) * ind.max)
   rmin <- array(0, length(sizes) * ind.max)
   
   
   for (s in 1:length(sizes)) {
      for (m in 1:ind.max) {
                
         T <- transf.func(load.data(filename, sizes[s], taus, ms[m], taus2, col, sf, nm = nm))
         print(paste(sizes[s], ms[m], nrow(T), ncol(T)))
         
         Tm <- apply(T,2,mean)
         k <- pick.func(Tm)
         ind <- (s - 1) * ind.max + m
         res[ind] <- Tm[k]
         if (check.normallity && sum(T[,k]) > 0 && shapiro.test(T[,k])$p.value < 0.01) {
                   print(paste("Warning: the data for", filename,
                  "column", k, "doesn't seem to be normally distributed."))
                    }
         ci <- qnorm(1-(1-confidence.level)/2)*sd(T[,k])/sqrt(nrow(T))
         rmax[ind] <- res[ind] + ci
         rmin[ind] <- res[ind] - ci 
         }
      } 
   list(y = res, ub = rmax, lb = rmin)
   }


## PISF ##

plot.mountain.grid <- function( res = seq(2, 10, by = 1),  res2 = c(3, 5, 6, 10), 
    rda = seq(2,5,by=1), pisf.res = 10, taus = 0.01, col=2, cols = 2) {


#     x11(w=5, h=5)
#     y <- gen.data("~/dp/experiments/mountain/mountain_grid_policy_iteration", res^2, col = col, pick.func=which.min)
#     plot.algs(res^2, y$y, Y.ub = y$ub, Y.lb = y$lb)
  
    for (col in cols) {
#       x11(w=5, h=5)
      y <- matrix(0, length(rda), length(res2))
      for (i in 1:ncol(y)) {
         L <- read.table(ps("~/dp/experiments/mountain/mountain_grid_policy_iteration_", res2[i]^2,".res"))
         y[,i] <- rep(apply(L,2,mean)[col], length(rda))
         }
      y <- cbind(y, gen.data("~/dp/experiments/mountain/mountain_grid_pisf_s", pisf.res^2, taus =taus, m = rda^2*3 , col = col, pick.func = which.min)$y)
      leg <- c(paste("PI", res2^2), "PISF")
      plot.algs(rda^2*3, y, leg=leg, plot.leg=TRUE)
     }
}


plot.pole.grid.pi <- function( res = seq(2, 14, by = 1),  col=3, ylab="Successful episodes") {

    y <- gen.data("~/dp/experiments/pole/pole_grid_policy_iteration", res^4, col = col, pick.func=which.max)
    plot.algs(res^4, y$y, Y.ub = y$ub, Y.lb = y$lb, xlab="|S|", ylab = ylab)
    }


plot.pole.grid.pisf.by.n <- function( res = c(5,10), pisf.res = 7:10,  rda = 4, taus = 0.01, col=3, 
                           ylab="Successful episodes") {

    y <- matrix(0, length(pisf.res), length(res)+1)
    for (i in 1:length(res)) {
      L <- read.table(ps("~/dp/experiments/pole/pole_grid_policy_iteration_", res[i]^4,".res"))
      y[,i] <- rep(apply(L,2,mean)[col], length(pisf.res))
      }
    
    y.ub <- y
    y.lb <- y
    for (i in 1:length(pisf.res)) {
      D <- gen.data("~/dp/experiments/pole/pole_grid_pisf_s", pisf.res[i]^4, taus =taus, m = rda^4*2 , col = col, pick.func = which.max)
      y[i, length(res)+1] <- D$y
      if (!is.na(D$ub)) y.ub[i, length(res)+1] <- D$ub
      else y.ub[i, length(res)+1] <- D$y
      if (!is.na(D$lb)) y.lb[i, length(res)+1] <- D$lb
      else y.lb[i, length(res)+1] <- D$y
      }
      leg <- NULL
      for (i in 1:length(res)) leg <- c(leg, paste("PI", res[i]^4))
      leg <- c(leg, ps("PISF(|S|,", rda^4*2, ")"))
      plot.algs(pisf.res^4, y, Y.ub = y.ub, Y.lb = y.lb, leg=leg, plot.leg=TRUE, leg.pos="bottomright", 
                ylab=ylab, xlab="|S|")
     }


plot.pole.grid.pisf.by.m <- function( res = 15, pisf.res = 20,  rda = 2:4, taus = 0.01, col=3, 
                           ylab="Successful episodes") {

    y <- matrix(0, length(rda), length(res))
    for (i in 1:length(res)) {
      L <- read.table(ps("~/dp/experiments/pole/pole_grid_policy_iteration_", res[i]^4,".res"))
      y[,i] <- rep(apply(L,2,mean)[col], length(pisf.res))
      }
    
    D <- gen.data("~/dp/experiments/pole/pole_grid_pisf_s", pisf.res^4, taus =taus, m = rda^4*2 , col = col, pick.func = which.max)
    y <- cbind(y, D$y)
    y.ub <- cbind(y, D$ub)
    y.lb <- cbind(y, D$lb)
    leg <- NULL
    for (i in 1:length(res)) leg <- c(leg, paste("PI", res[i]^4))
    leg <- c(leg, ps("PISF(", pisf.res^4, ", m)"))
    
    plot.algs(rda^4*2, y, Y.ub = y.ub, Y.lb = y.lb, leg=leg, plot.leg=TRUE, leg.pos="bottomright", 
                ylab=ylab, xlab="m")
    }


plot.double.pole.grid.pi <- function( res = seq(2, 7, by = 1),  col=3,  ylab= "Successful episodes",
                                       transf.func = id) {

    y <- gen.data("~/dp/experiments/double_pole/double_pole_grid_policy_iteration", res^6, col = col, 
                  pick.func=which.max, transf.func = transf.func)
    plot.algs(res^6, y$y, Y.ub = y$ub, Y.lb = y$lb, xlab="|S|", ylab = ylab)
    }


plot.puddle.grid <- function( ns = seq(2, 15, by = 1), col = 1, ylab="Return") {
    y <- gen.data("~/dp/experiments/puddle/puddle_grid_policy_iteration", ns^2, col = col, pick.func=which.max)
    plot.algs(ns^2, y$y, Y.ub = y$ub, Y.lb = y$lb, ylab = ylab, xlab="|S|")
    }

plot.puddle.grid.pisf <- function( n=15,  ms = seq(2,9,by=1)^2*4, taus = c(0.3), col = 1, ylab="Return") {
      y <- matrix(0, length(ms), 1)
      L <- read.table(ps("~/dp/experiments/puddle/puddle_grid_policy_iteration_", n^2,".res"))
      y[,1] <- rep(apply(L,2,mean)[col], length(ms))
      y <- cbind(y, gen.data("~/dp/experiments/puddle/puddle_grid_pisf_s", n^2, taus =taus, m = ms , col = col, pick.func = which.max)$y)
      
      leg <- c(ps("PI(", n^2, ")"), ps("PISF(", n^2, "," , "m)"))
      plot.algs(ms, y, leg = leg, plot.leg = TRUE, leg.pos="bottomright", ylab = ylab, xlab = "m")
     }

    

generate.mountain.vf <- function( res.pi = 100, res.pisf = seq(10,50,by=10) ) {
  cp <- seq(-1.2, 0.5, length = 100)
  cv <- seq(-0.07, 0.07, length = 100)

  V <- as.matrix(read.table(ps("~/dp/experiments/mountain/pi_value_function_",res.pi,".mat"))) 
  image(cp, cv, V, xlab="", ylab="")
  dev.copy2eps(file = ps("~/tex/canada/fig/pi_value_function_", res.pi, ".eps"))

  for (rp in res.pisf) {
    V <- as.matrix(read.table(ps("~/dp/experiments/mountain/pisf_value_function_",rp,".mat"))) 
    image(cp, cv, V, xlab="", ylab="")
    dev.copy2eps(file = ps("~/tex/canada/fig/pisf_value_function_", rp, ".eps"))
    }
  }
  

plot.pole.grid.pi <- function( col=3,  ns = seq(2, 8, by = 1),
                  filename = "~/dp/experiments/pole/easy_pole_grid_policy_iteration" ) {
    D  <- gen.data(filename, ns^4, col = col)
    plot.algs(ns^4, D$y, Y.lb = D$lb, Y.ub = D$ub)
  }
    


plot.pole.grid.pisf.pi.by.m <- function(col=3, n.pi = 8, n.pisf = 8, ms = 2:5, taus = seq(0.1,0.9,by=0.2),
       filename = "~/dp/experiments/pole/easy_pole_grid") {

    ms <- ms^4*2
    B <- gen.data(ps(filename,"_policy_iteration"), n.pi^4, col=col)
    D <- gen.data(ps(filename, "_pisf_s"), n.pisf^4, taus, ms, col = col)
    yy <- cbind(B$y, D$y)
    ub <- cbind(B$ub, D$ub)
    lb <- cbind(B$lb, D$lb)

   plot.algs(ms, yy, Y.lb = lb, Y.ub = ub)
  }


plot.pole.grid.pisf.pi.by.n <- function(col=3, n.pi = 8, ns = 8:15, m = 4, taus = seq(0.1,0.9,by=0.2),
       filename = "~/dp/experiments/pole/easy_pole_grid") {
#     plot.pole.grid.pisf.pi.by.n(4, n.pi=8, ns = 6:10, m=4)
    m <- m^4*2
    B <- gen.data(ps(filename,"_policy_iteration"), n.pi^4, col=col)
    D <- gen.data(ps(filename, "_pisf_s"), ns^4, taus, m, col = col)
    yy <- cbind(B$y, D$y)
    ub <- cbind(B$ub, D$ub)
    lb <- cbind(B$lb, D$lb)

   plot.algs(ns^4, yy, Y.lb = lb, Y.ub = ub)
  }


plot.double.pole.grid.pisf.kbrl <- function(col=3, n.kbrl = 30000, n.pisf = 200000, ms = c(seq(100,400,by=100),500), taus = c(10,1,0.1),
       filename = "~/dp/experiments/double_pole/double_pole", transf.func = id) {

    B <- gen.data(ps(filename,"_kbrl"), n.kbrl, col=col, taus = taus, transf.func = transf.func)
    D <- gen.data(ps(filename, "_pisf"), n.pisf, taus = taus, taus2=taus, ms = ms, col = col, transf.func = transf.func)
    yy <- cbind(B$y, D$y)
    ub <- cbind(B$ub, D$ub)
    lb <- cbind(B$lb, D$lb)

   plot.algs(ms, yy, Y.lb = lb, Y.ub = ub)
  }
 
 
 
 plot.epilepsy.pi <- function(col = 1, ns = 2:12, filename = "~/ep_grid/experiments/before_optimization/epilepsy_grid_policy_iteration",  transf.func = id) {
    B <- gen.data(filename, ns^5, col=col, transf.func = transf.func)
    plot.algs(ns^5, B$y, xlab = "|S|", ylab= "Return")
    }

## KBSF ##

show.mnt.kbrl <- function(sizes=seq(1000, 9000, by = 1000) ) {
   T <- gen.data("~/dp/experiments/mountain/mountain_kbrl", sizes = sizes, taus = c(0.01), col = 2, pick.func=which.min)
   plot.algs(sizes, Y = T$y, Y.lb = T$lb, Y.ub = T$ub, ylim=c(0,300), plot.leg=TRUE,
            ylab = "Steps to escape", xlab = "n",  leg = "KBRL(n)")
   }

   
show.mnt.m <- function(size=5000, ms = seq(10,150,by=20), col = 2, pick.func = which.min, ylim = c(0,300),
                       taus = c(1, 0.1, 0.01), taus2 = c(1, 0.1, 0.01), show.log = FALSE , ylab = "Steps to escape",
                       transf.func = id) {

   T <- gen.data("~/dp/experiments/mountain/mountain_kbrl", sizes = size, taus = 0.01, 
                 col = col, pick.func=pick.func, transf.func = transf.func)
   T$y <- c(T$y, rep(T$y, length(ms) - 1) )
   T$ub <- c( T$ub, rep(T$ub, length(ms) - 1) )
   T$lb <- c( T$lb, rep(T$lb, length(ms) - 1) )
   
     
   K <- gen.data("~/dp/experiments/mountain/mountain_kbsf", sizes = size, taus = taus, 
                taus2 = taus2, col = col, pick.func=pick.func, ms = ms,
                transf.func = transf.func)

   T$y <-  cbind(T$y, K$y)
   T$ub <- cbind( T$ub, K$ub )
   T$lb <- cbind( T$lb, K$lb )

   plot.algs(ms, Y = T$y, Y.lb = T$lb, Y.ub = T$ub, ylim=ylim, 
             alg.ind = c(1,3), plot.leg = TRUE, ylab = ylab, 
             xlab = "m", show.int = c(FALSE, TRUE), pch = c(NA, alg.pch(3)[3]),
             leg = c("KBRL(5000)", "KBSF(5000,m)"), leg.pos = "bottomright")
            
   }


show.mnt.n <- function(sizes= seq(1000, 10000, by = 1000), m = 100, col = 2, pick.func = which.min, ylim = c(0,300),
                       taus = c(1, 0.1, 0.01), taus2 = c(1, 0.1, 0.01), ylab = "Steps to escape",
                       transf.func = id, ...) {
   T <- gen.data("~/dp/experiments/mountain/mountain_kbrl", sizes = sizes, taus = 0.01, col = col, 
                 pick.func=pick.func, transf.func = transf.func)

   L <- gen.data("~/dp/experiments/mountain/mountain_kbsf", sizes = sizes, taus = c(1, 0.1, 0.01),
                taus2 = c(1, 0.1, 0.01), col = col, pick.func=pick.func, ms = m, transf.func = transf.func)

   T$y <-  cbind(T$y, L$y)
   T$ub <- cbind( T$ub, L$ub )
   T$lb <- cbind( T$lb, L$lb )
   
   plot.algs(sizes, Y = T$y, Y.lb = T$lb, Y.ub = T$ub, ylim = ylim, 
             alg.ind = c(1,3), plot.leg = TRUE, leg = c("KBRL(n)", "KBSF(n,100)"),
             ylab = ylab, xlab = "n", ...)
   }


show.pd.m <- function(size=8000, ms = seq(10,150,by=20), col = 1, pick.func = which.max, ylim = c(0.3,3.2),
                       taus = c(1, 0.1, 0.01), taus2 = c(1, 0.1, 0.01), show.log = FALSE, ylab = "Return",
                       transf.func = id, leg.pos = "bottomright", dir = "~/fontes/R/kbsf/puddle/", 
                       nm = NULL, ...) {
   T <- gen.data(ps(dir, "puddle_kbrl"), sizes = size, taus = taus, 
                 col = col, pick.func=pick.func, transf.func = transf.func, nm = nm)
   T$y <- c(T$y, rep(T$y, length(ms) - 1) )
   T$ub <- c( T$ub, rep(T$ub, length(ms) - 1) )
   T$lb <- c( T$lb, rep(T$lb, length(ms) - 1) )

## LSPI
#    L <- gen.data("~/dp/experiments/puddle/puddle_lspi", sizes = size, taus = taus, 
#                 col = col, pick.func=pick.func, ms = ms)
## KBRL(m)
#     R <- gen.data("~/dp/experiments/puddle/puddle_kbrl", sizes = ms, taus = c(1), #0.1 and 0.01 generated NAN
#                  col = col, pick.func=pick.func)
    K <- gen.data(ps(dir,"puddle_kbsf"), sizes = size, taus = taus, 
                 taus2 = taus2, col = col, pick.func=pick.func, ms = ms, transf.func = transf.func, nm = nm)

   T$y <-  cbind(T$y, K$y)
   T$ub <- cbind( T$ub, K$ub )
   T$lb <- cbind( T$lb, K$lb )
   
   if (show.log) {
      T$y <-  log(T$y)
      T$ub <- log(T$ub)
      T$lb <- log(T$lb)
      }

   plot.algs(ms, Y = T$y, Y.lb = T$lb, Y.ub = T$ub, ylim=ylim, 
             alg.ind = c(1,3), plot.leg = TRUE, ylab = ylab, 
             xlab = "m", show.int = c(FALSE, TRUE), pch = c(NA, alg.pch(3)[3]),
             leg = c(ps("KBRL(",size,")"), ps("KBSF(", size,",m)")), leg.pos = leg.pos, ...)
   }


show.pd.n <- function(sizes= seq(1000, 10000, by = 1000), m = 100, col = 1, pick.func = which.max, ylim = c(0.3,3.2),
                       taus = c(1, 0.1, 0.01), taus2 = c(1, 0.1, 0.01), show.log = FALSE, 
                       transf.func = id, ylab = "Return", leg.pos = "bottomright", 
                       dir = "~/fontes/R/kbsf/puddle/",... ) {
   T <- gen.data(ps(dir,"puddle_kbrl"), sizes = sizes, taus = 0.1, 
                 col = col, pick.func=pick.func, transf.func = transf.func)

   L <- gen.data(ps(dir,"puddle_kbsf"), sizes = sizes, taus = c(1, 0.1, 0.01), 
                taus2 = c(1, 0.1, 0.01), col = col, pick.func=pick.func, ms = m, transf.func = transf.func)

   T$y <-  cbind(T$y, L$y)
   T$ub <- cbind( T$ub, L$ub )
   T$lb <- cbind( T$lb, L$lb )

   if (show.log) {
      T$y <-  log(T$y)
      T$ub <- log(T$ub)
      T$lb <- log(T$lb)
      }

   plot.algs(sizes, Y = T$y, Y.lb = T$lb, Y.ub = T$ub, alg.ind = c(1,3), 
            plot.leg = TRUE, xlab = "n", ylab = ylab, leg = c("KBRL(n)", "KBSF(n,100)"), 
            leg.pos = leg.pos, ...)
   }
   


show.pl.m <- function(size=50000, ms = seq(10,150,by=20), col = 3, pick.func = which.max, ylim = c(0,1),
                       taus = c(1, 0.1, 0.01), taus2 = 1, transf.func = id, ylab = "Successful episodes", 
                       leg.pos = "bottomright", dir = "~/fontes/R/kbsf/pole/", 
                       show.tree = FALSE, nt = 30, mnes = c(2000,200,20), ncp = 4, pr = 0.6, 
                      confidence.level = 0.99, ... ) {
## LSPI
    L <- gen.data(ps(dir,"pole_lspi"), sizes = size, taus = taus, 
                 col = col, pick.func=pick.func, ms = ms, transf.func = transf.func)
   
    K <- gen.data(ps(dir,"pole_kbsf"), sizes = size, taus = taus, 
                taus2 = taus2, col = col, pick.func=pick.func, ms = ms, transf.func = transf.func)
    
    L$y <-  cbind(L$y, K$y)
    L$ub <- cbind( L$ub, K$ub )
    L$lb <- cbind( L$lb, K$lb )

    if (show.tree)
    {
       T    <- matrix(0, length(ms), length(mnes))
       T.ub <- T
       T.lb <- T
       for (i in 1:length(mnes))
       {
          fn <- ps(paste(ps(dir,"pole_tree"), size, nt, mnes[i], ncp, pr, sep="_"), ".res")
          X <- read.table(fn)
          mx <- mean(X[,col])
          ci <- qnorm(1-(1-confidence.level)/2)*sd(X[,col])/sqrt(nrow(X))
          T[,i] <- mx
          T.ub[,i] <- mx + ci
          T.lb[,i] <- mx - ci 
       }
      
      L$y <-  cbind(L$y, T)
      L$ub <- cbind( L$ub, T.ub )
      L$lb <- cbind( L$lb, T.lb )
    }
    
    leg <- c(expression("LSPI(5x"*10^4*",m)"), expression("KBSF(5x"*10^4*",m)"))
    alg.ind <- c(2,3)   
    if (show.tree) 
    {
       leg <- c(leg, paste("TREE-", mnes))
       alg.ind <- c(alg.ind, c(1,4,5)) ##FIX
    }
    
    L$y[L$y == 0] <- 1e-5
    L$ub[L$ub == 0] <- 1e-5
    L$lb[L$lb == 0] <- 1e-5
    
    plot.algs(ms, Y = L$y, Y.lb = L$lb, Y.ub = L$ub, ylim=ylim, ylab = ylab, alg.ind = alg.ind, plot.leg = TRUE, 
            leg.pos = leg.pos, xlab = "m", leg = leg, ...)
   }


show.dpl.m <- function(size=1000000, ms = seq(20,200,by=20), col = 3, pick.func = which.max, ylim = NULL,
                       taus = c(1, 0.1, 0.01), taus2 = 1, transf.func = id, ylab= "Successful episodes",
                       leg.pos = "bottomright", y.inter = 1, inset = c(0.0, 0.25), 
                       dir = "~/fontes/R/kbsf/double_pole/", ...) {
    L <- gen.data(ps(dir,"double_pole_lspi"), sizes = size, taus = taus, 
                 col = col, pick.func=pick.func, ms = ms, transf.func = transf.func)
   
    K <- gen.data(ps(dir,"double_pole_kbsf"), sizes = size, taus = taus, 
                taus2 = taus2, col = col, pick.func=pick.func, ms = ms, transf.func = transf.func)

   L$y <-  cbind(L$y, K$y)
   L$ub <- cbind( L$ub, K$ub )       
   L$lb <- cbind( L$lb, K$lb )

 
   plot.algs(ms, Y = L$y, Y.lb = L$lb, Y.ub = L$ub, ylim=ylim, ylab = ylab, alg.ind = c(2,3), plot.leg = TRUE, 
            leg.pos = leg.pos, xlab = "m", leg = c(expression("LSPI("*10^6*",m)"), expression("KBSF("*10^6*",m)")), 
            y.inter = y.inter, inset = inset, ...)
   }
   


show.dpl.times <- function(size1=50000, ms1 = seq(10,100,by=10), size2=1000000, ms2 = seq(20,200,by=20), 
                       col = 4, pick.func = which.max, transf.func = log, 
                       taus = c(1, 0.1, 0.01), taus2 = 1, ylab= "Seconds",
                       leg.pos = "bottomright", dir.pl = "./kbsf/pole/", 
                       dir.dpl = "./kbsf/double_pole/", ...) {
    L1 <- gen.data(ps(dir.pl,"pole_lspi"), sizes = size1, taus = taus, 
                 col = col, pick.func=pick.func, ms = ms1, transf.func = transf.func)
   
    K1 <- gen.data(ps(dir.pl,"pole_kbsf"), sizes = size1, taus = taus, 
                taus2 = taus2, col = col, pick.func=pick.func, ms = ms1, transf.func = transf.func)
   
    L2 <- gen.data(ps(dir.dpl, "double_pole_lspi"), sizes = size2, taus = taus, 
                 col = col, pick.func=pick.func, ms = ms2, transf.func = transf.func)
   
    K2 <- gen.data(ps(dir.dpl,"double_pole_kbsf"), sizes = size2, taus = taus, 
                taus2 = taus2, col = col, pick.func=pick.func, ms = ms2, transf.func = transf.func)

   L1$y <-  cbind(L1$y,K1$y)
   L1$ub <- cbind( L1$ub, K1$ub)       
   L1$lb <- cbind( L1$lb, K1$lb)

   L2$y <-  cbind(L2$y, K2$y)
   L2$ub <- cbind( L2$ub, K2$ub )       
   L2$lb <- cbind( L2$lb, K2$lb )
   ylim <- c(min(L1$lb, L2$lb), max(L1$ub, L2$ub))
 
   plot.algs(ms1, Y = L1$y, Y.lb = L1$lb, Y.ub = L1$ub, ylim=ylim, ylab = ylab, alg.ind = c(2,3), plot.leg = TRUE, 
            leg.pos = leg.pos, xlab = "m", leg = c(expression("LSPI("*10^6*",m)"), expression("KBSF("*10^6*",m)")), 
            y.inter = 1, inset = c(0.0, 0.25), xlim = c(min(ms1),max(ms2)), ...)

   plot.algs(ms2, Y = L2$y, Y.lb = L2$lb, Y.ub = L2$ub, ylim=ylim, ylab = ylab, alg.ind = c(2,3), plot.leg = TRUE, 
            leg.pos = leg.pos, xlab = "m", leg = c(expression("LSPI("*10^6*",m)"), expression("KBSF("*10^6*",m)")), 
            y.inter = 1, inset = c(0.0, 0.25), new = FALSE, ...)

}
   

generate.all <- function(dir = "~/tmp/")  {
   x11(w=6, h=5)

   ## MOUNTAIN ##
#    show.mnt.m()
#    dev.copy2eps(file = "mountain_m.eps", font="Bookman")
#    
#    show.mnt.m(col = 4, transf.func = log, ylab = "Seconds (log)", ylim=NULL)
#    dev.copy2pdf(file = ps(dir, "mountain_m_time.pdf"), font="Bookman")
# 
#    show.mnt.n()
#    dev.copy2pdf(file = ps(dir, "mountain_n.pdf"), font="Bookman")
#    
#    show.mnt.n(col = 4, ylab = "Seconds", ylim=NULL, leg.pos="topleft")
#    dev.copy2pdf(file = ps(dir, "mountain_n_time.pdf"), font="Bookman")


   ## PUDDLE ##
   show.pd.m()
   dev.copy2pdf(file = ps(dir, "puddle_m.pdf"), font="Bookman")
   
   show.pd.m(col = 4, log="y", ylab = "Seconds (log)", ylim=NULL, leg.pos = "left")
   dev.copy2pdf(file = ps(dir, "puddle_m_time_log.pdf"), font="Bookman")

   show.pd.m(col = 4, transf.func = id, ylab = "Seconds", ylim=NULL, leg.pos = "left")
   dev.copy2pdf(file = ps(dir, "puddle_m_time.pdf"), font="Bookman")


   show.pd.n()
   dev.copy2pdf(file = ps(dir, "puddle_n.pdf"), font="Bookman")
   
   show.pd.n(col = 4, ylab = "Seconds (log)", ylim=NULL, leg.pos="right", log="y")
   dev.copy2pdf(file = ps(dir, "puddle_n_time_log.pdf"), font="Bookman")

   show.pd.n(col = 4, ylab = "Seconds", ylim=NULL, leg.pos="topleft")
   dev.copy2pdf(file = ps(dir, "puddle_n_time.pdf"), font="Bookman")


   ## POLE ##
   show.pl.m()
   dev.copy2pdf(file = ps(dir, "pole_m.pdf"), font="Bookman")
   
   show.pl.m(col = 4, ylab = "Seconds", ylim = NULL, leg.pos = "topleft")
   dev.copy2pdf(file = ps(dir, "pole_m_time.pdf"), font="Bookman")
   
   show.pl.m(col = 4, ylab = "Seconds (log)", ylim = NULL, leg.pos = "bottomright", log="y", y.inter = 1, inset = c(0.0, 0.0))
   dev.copy2pdf(file = ps(dir, "pole_m_time_log.pdf"), font="Bookman")

## DOUBLE POLE ##
   show.dpl.m()
   dev.copy2pdf(file = ps(dir, "double_pole_m.pdf"), font="Bookman")
   
   show.dpl.m(col = 4, ylab = "Seconds", ylim = NULL, leg.pos = "topleft")
   dev.copy2pdf(file = ps(dir, "double_pole_m_time.pdf"), font="Bookman")

   show.dpl.m(col = 4, ylab = "Seconds (log)", ylim = NULL, leg.pos = "bottomright", log="y", y.inter = 1, inset = c(0.0, 0.0))
   dev.copy2pdf(file = ps(dir, "double_pole_m_time_log.pdf"), font="Bookman")
   }
   
   
   
## EPILEPSY

show.ep  <- function(col = 2,  boxplot = TRUE,
            freqs = c(seq(0.0, 2, by = 0.5)),
            n = 500000,
            m = 50000,
            taus1 = c(1, 0.1, 0.01),
            taus2 = c(1, 0.1, 0.01),
            ks = 3,
            confidence.level = 0.99,
            ylab = "Fraction of stimulation",
            dir = "~/ep/experiments/") {
            

   T <- list(read.table(ps(dir, "epilepsy_random.res"))[,col])
 
   for (i in 1:length(freqs)) {
      v <- read.table(ps(dir,"epilepsy_fixed_", freqs[i],".res"))[,col]
      T <- c(T, list(v))      
      }

#    minv <- Inf
#    K <- NULL
#    
#    for (t1 in 1:length(taus1)) {
#       for (t2 in 1:length(taus2)) {
#          for (k in 1:length(ks)) {
#             Q <- read.table(ps(dir,"epilepsy_kbsf_",sprintf("%5.0f", n),"_",m, "_", taus1[t1], "_", taus2[t2], "_", ks[k], ".res"))[,col]
#             meanv <- mean(Q) 
#             if (meanv < minv) {
#                K <- Q
#                minv <- meanv
#                }
#             }
#           }
#        }
   
#     K <- read.table(ps(dir,"epilepsy_kbsf_1000000_18000_0.1_1_3.res"))[,col]
   
#     T <- c(T, list(K))

    T <- c(T, list(read.table(ps(dir, "epilepsy_kbsf_", sprintf("%5.0f", n), "_", m, "_0.1_1_3.res"))[,col]))
    T <- c(T, list(read.table(ps(dir, "no_penalty/epilepsy_kbsf_", sprintf("%5.0f", n), "_", m, "_0.1_1_3.res"))[,col]))
    T <- data.frame(T)
    
    names <- c("rand", paste(freqs), "kbsf1", "kbsf2")
    colnames(T) <- names
    if (boxplot) boxplot(T, range = 0)
    else{
#        Y    <- array(0, length(T))
#        Y.ub <- array(0, length(T))
#        Y.lb <- array(0, length(T))
#        
#        for (k in 1:length(T)) {
#           Y[k] <- mean(T[[k]])
#           ci <- qnorm(1-(1-confidence.level)/2) * sd(T[[k]])/sqrt(length(T[[k]]))
#           Y.ub[k] <- Y[k] + ci
#           Y.lb[k] <- Y[k] - ci
#           }
#        plot.algs(c(-0.5, c(freqs), 2.5, 3.0), Y, Y.ub = Y.ub, Y.lb = Y.lb, show.shad = FALSE, 
#                 show.int = TRUE, t = "p" )

       error.bars(T, alpha = 0.01, xlab="Method", ylab=ylab)
       }
     T <- data.frame(T)
     colnames(T) <- names
     T
     }
        


find.best.mne.pareto <- function(cols = c(2,3), 
      n.tree = 200000, nt = 30, mnes = 10, 
      dir, col = NULL) {
         if (is.null(col)) col <- cols 
         B <- NULL
         min.sum <- Inf
         for (ne in mnes) {
            T <- read.table(ps(dir, "epilepsy_tree_rl_", sprintf("%5.0f", n.tree), "_", nt, "_", ne, "_-1_0.6.res"))
            s <- sum(apply(T[,cols], 2, mean))
            if (s < min.sum) {
               min.sum <- s
               B <- T[,col]
               }
           }
      B
      }


find.best.m.pareto <- function(filename, n, taus, m, cols, taus2 = NULL, sf = ".res" , return.cols = NULL) {
  
  if (is.null(return.cols)) return.cols = cols

  T <- load.data(filename, n, taus, m, taus2 = taus2, sf=sf, col=cols[1])
  T <- rbind(T, load.data(filename, n, taus, m,  taus2 = taus2, sf=sf, col=cols[2]))
  
  f <- apply(T, 2, mean)
  idx <- which.min(f)
  
  print(paste(filename, taus[idx]))
  
  T <- NULL
  for (i in 1:length(return.cols)) T <- cbind(T,load.data(filename, n, taus[idx], m,  taus2 = taus2, sf = sf, col = return.cols[i]))
  T
  }
      

show.ep.pareto  <- function(cols = c(2,3), 
            freqs = c(seq(0.0, 1.5, by = 0.5)),
            n = 500000,
            m = 50000,
            m.lspi = 50,
            nt = 30,
            mne =  c(20, 200),
            taus1 = c(1, 0.1, 0.01),
            confidence.level = 0.99,
            dir = "./kbsf/ep/",
            delta = 0.005) {

#      x11(w=6, h=6)     
     set.par()
     par(ps=18)
     
     T <- NULL
    

   for (i in 1:length(freqs)) {
      v <- read.table(ps(dir,"epilepsy_fixed_", freqs[i],".res"))[,cols]
      T <- c(T, list(v))      
      }

      T <- c(T, list(find.best.m.pareto(ps(dir, "small_penalty/epilepsy_lspi"), n, taus1, m.lspi, cols)))
      T <- c(T, list(find.best.m.pareto(ps(dir, "medium_penalty/epilepsy_lspi"), n, taus1, m.lspi, cols)))
      T <- c(T, list(find.best.m.pareto(ps(dir, "large_penalty/epilepsy_lspi"), n, taus1, m.lspi, cols)))

      T <- c(T, list(read.table(ps(dir, "small_penalty/epilepsy_tree_rl_", sprintf("%5.0f", n),"_", nt, "_", mne[1], "_-1_0.6.res"))[,cols]))
      T <- c(T, list(read.table(ps(dir, "medium_penalty/epilepsy_tree_rl_", sprintf("%5.0f", n),"_", nt, "_", mne[1], "_-1_0.6.res"))[,cols]))
      T <- c(T, list(read.table(ps(dir, "large_penalty/epilepsy_tree_rl_", sprintf("%5.0f", n),"_", nt, "_", mne[1], "_-1_0.6.res"))[,cols]))

      T <- c(T, list(read.table(ps(dir, "small_penalty/epilepsy_tree_rl_", sprintf("%5.0f", n),"_", nt, "_", mne[2], "_-1_0.6.res"))[,cols]))
      T <- c(T, list(read.table(ps(dir, "medium_penalty/epilepsy_tree_rl_", sprintf("%5.0f", n),"_", nt, "_", mne[2], "_-1_0.6.res"))[,cols]))
      T <- c(T, list(read.table(ps(dir, "large_penalty/epilepsy_tree_rl_", sprintf("%5.0f", n),"_", nt, "_", mne[2], "_-1_0.6.res"))[,cols]))

      T <- c(T, list(find.best.m.pareto(ps(dir, "small_penalty/epilepsy_kbsf"), n, taus1, m, taus2 = 1, cols, sf="_6.res")))
      T <- c(T, list(find.best.m.pareto(ps(dir, "medium_penalty/epilepsy_kbsf"), n, taus1, m, taus2 = 1,  cols, sf="_6.res")))
      T <- c(T, list(find.best.m.pareto(ps(dir, "large_penalty/epilepsy_kbsf"), n, taus1, m, taus2 = 1, cols, sf="_6.res")))


      x <- array(0, length(T))
      y <- array(0, length(T))
      ubx <- array(0, length(T))
      lbx <- array(0, length(T))
      uby <- array(0, length(T))
      lby <- array(0, length(T))
      
      
      for (i in 1:length(T)) {
         x[i] <- mean(T[[i]][,2])
         y[i] <- mean(T[[i]][,1])
         
         ci <- qnorm(1-(1-confidence.level)/2)*sd(T[[i]][,2])/sqrt(length(T[[i]][,2]))
         ubx[i] <- x[i] + ci
         lbx[i] <- x[i] - ci
         
         ci <- qnorm(1-(1-confidence.level)/2)*sd(T[[i]][,1])/sqrt(length(T[[i]][,1]))
         uby[i] <- y[i] + ci
         lby[i] <- y[i] - ci
         }
     
      col <- c(rep("BLACK", length(freqs)), rep("RED",3), rep("GREEN", 3), rep("LIGHTGREEN", 3), rep("BLUE",3))

      plot(NULL, xlim = c(min(lbx)-0.01, max(ubx)+0.01), ylim = c(min(lby) - 2 * delta, max(uby)), xlab = "Fraction of stimulation", ylab = "Fraction of seizures")
    
#      plot(x,y, ylim = c(0.1,0.25), col = col, pch=pch, , xlab = "Fraction of stimulation", ylab = "Fraction of seizures")

      for (i in 1:length(y)) {
          polygon(c(lbx[i], ubx[i], ubx[i], lbx[i], lbx[i]), c(lby[i], lby[i], uby[i], uby[i], lby[i]), col=col[i], border = col[i])
          }
     
#       ind <- c(length(freqs) + 1, length(freqs) + 3)
#       text(x[-ind], lby[-ind] - delta, c(ps(freqs,"Hz"), 
#       expression("LSPI"[-20]), expression("KBSF"[-40]), expression("KBSF"[-20]),expression("KBSF"[l])), col = col[-ind])
# #      text(x[ind[1]], uby[ind[1]] + delta, expression("LSPI"[-40]), col = col[ind[1]])
#     text(x[ind[2]] + 4 * delta, y[ind[2]], expression("LSPI"[l]), col = col[ind[2]])

      labels <- c(ps(freqs,"Hz"),
      expression("LSPI"[-40]), expression("LSPI"[-20]),expression("LSPI"[-10]), 
      expression("FQIT(20)"[-40]), expression("FQIT(20)"[-20]),expression("FQIT(20)"[-10]), 
      expression("FQIT(200)"[-40]), expression("FQIT(200)"[-20]),expression("FQIT(200)"[-10]), 
      expression("KBSF"[-40]), expression("KBSF"[-20]),expression("KBSF"[-10]))
  
      rm <- c(5, 6, 7, 9, 12)
      cex <- 0.8
      text(x[-rm], lby[-rm] - delta, labels[-rm], cex=cex)

       text(x[rm[1]] + 2 * delta, lby[rm[1]] - delta, labels[rm[1]], cex=cex)
       text(x[rm[2]] - 2 * delta, lby[rm[2]] - delta, labels[rm[2]], cex=cex)
       text(x[rm[3]] + 5 * delta, lby[rm[3]] + 0.5 * delta, labels[rm[3]], cex=cex)
       text(x[rm[4]] + 5 * delta, lby[rm[4]] + 3 * delta, labels[rm[4]], cex=cex)
       text(x[rm[5]] + 5 * delta, lby[rm[5]] + 3 * delta, labels[rm[5]], cex=cex)
       
#       leg <- c(
#       expression("FQI(5x"*10^5*") with "*eta*"=20"), 
#       expression("FQI with "*eta*"=200"), 
#       expression("LSPI(5x"*10^5*",50)"),
#       expression("KBSF(5x"*10^5*",5x"*10^4*")")
#       )
#       legend("bottomleft", leg, col= c("GREEN", "LIGHTGREEN", "RED", "BLUE"), bg="WHITE",
#             lty=1, pch=NULL, lwd = 3, y.intersp = 1, pt.bg=c("GREEN", "LIGHTGREEN", "RED", "BLUE"), inset = 0)      
     }


     
show.ep.times <- function(
      col = 4, 
      n = 500000, 
      m = 50000, 
      m.lspi = 50,
      nt = 30, 
      mne = c(20, 200), 
      taus1 = c(1, 0.1, 0.01),
      dir = "./kbsf/ep/", 
      confidence.level = 0.99, 
      cols.pareto=c(2,3)) {
      
#       x11(h = 8, w=4)
#       set.par()
      par(ps = 8)
      
      T <- NULL

      # Small penalty
      v <- read.table(ps(dir, "small_penalty/epilepsy_tree_rl_", sprintf("%5.0f", n),"_", nt, "_", mne[1], "_-1_0.6.res"))[,col]
      ## remove when experiments are finished
      v <- c(v, rep(mean(v), 500 - length(v)))
      T <- c(T, list(v))
      T <- c(T, list(read.table(ps(dir, "small_penalty/epilepsy_tree_rl_", sprintf("%5.0f", n),"_", nt, "_", mne[2], "_-1_0.6.res"))[,col]))

      T <- c(T, list(find.best.m.pareto(ps(dir, "small_penalty/epilepsy_lspi"), n, taus1, m.lspi, cols.pareto, return.col = col)))

      T <- c(T, list(find.best.m.pareto(ps(dir, "small_penalty/epilepsy_kbsf"), n, taus1, m, cols.pareto, taus2 = 1, sf = "_6.res", return.col = col)))
      
      
      # Medium penalty
      v <- read.table(ps(dir, "medium_penalty/epilepsy_tree_rl_", sprintf("%5.0f", n),"_", nt, "_", mne[1], "_-1_0.6.res"))[,col]
      ## remove when experiments are finished
      v <- c(v, rep(mean(v), 500 - length(v)))
      T <- c(T, list(v))
      
      T <- c(T, list(read.table(ps(dir, "medium_penalty/epilepsy_tree_rl_", sprintf("%5.0f", n),"_", nt, "_", mne[2], "_-1_0.6.res"))[,col]))

      T <- c(T, list(find.best.m.pareto(ps(dir, "medium_penalty/epilepsy_lspi"), n, taus1, m.lspi, cols.pareto, return.col = col)))

      T <- c(T, list(find.best.m.pareto(ps(dir, "medium_penalty/epilepsy_kbsf"), n, taus1, m, cols.pareto, taus2 = 1, sf = "_6.res", return.col = col)))

     
      # Large penalty
      v <- read.table(ps(dir, "large_penalty/epilepsy_tree_rl_", sprintf("%5.0f", n),"_", nt, "_", mne[1], "_-1_0.6.res"))[,col]
      ## remove when experiments are finished
      v <- c(v, rep(mean(v), 500 - length(v)))
      T <- c(T, list(v))

      T <- c(T, list(read.table(ps(dir, "large_penalty/epilepsy_tree_rl_", sprintf("%5.0f", n),"_", nt, "_", mne[2], "_-1_0.6.res"))[,col]))

      T <- c(T, list(find.best.m.pareto(ps(dir, "large_penalty/epilepsy_lspi"), n, taus1, m.lspi, cols.pareto, return.col = col)))

      T <- c(T, list(find.best.m.pareto(ps(dir, "large_penalty/epilepsy_kbsf"), n, taus1, m, cols.pareto, taus2 = 1, sf = "_6.res", return.col = col)))


      T <- data.frame(T)
    
      uby <- array(0, ncol(T))
      lby <- array(0, ncol(T))
      y <- array(0, ncol(T))
      
      for (i in 1:ncol(T)) {
        
         ci <- qnorm(1-(1-confidence.level)/2)*sd(T[,i])/sqrt(nrow(T))
         y[i] <- mean(T[,i])
         uby[i] <- y[i] + ci
         lby[i] <- y[i] - ci
         }

      col <- c(rep(c("GREEN", "LIGHTGREEN", "RED", "BLUE"), 3))

      set.par()
      par(ps=12)
      par(mar = c(5, 7, 2, 2) -1 + 0.1)
      
      
      spc <- rep(0.1, 12)
      spc[c(5,9)] <- 0.5
      
      names <- c(expression("FQIT(20)"[-10]), expression("FQIT(200)"[-10]), expression("LSPI"[-10]), expression("KBSF"[-10]),
        expression("FQIT(20)"[-20]), expression("FQIT(200)"[-20]), expression("LSPI"[-20]), expression("KBSF"[-20]),
        expression("FQIT(20)"[-40]), expression("FQIT(200)"[-40]), expression("LSPI"[-40]), expression("KBSF"[-40]))

      y <- barplot(lby, col = col, ylab = "", xlab = "Seconds (log)", names = names,
        border = "WHITE", log= "x", space = spc, horiz = TRUE, beside = FALSE, las=1)

      
#         barplot(y, col = col, add = TRUE, border = "WHITE", log= "y", space = spc, horiz = TRUE)
#         barplot(uby, col = col, add = TRUE, border = "WHITE", log= "y", space = spc, horiz = TRUE)
     
#       mtext(, 1)
      }
      


show.ep.features <- function(col = 2, 
            freqs = c(seq(0.0, 2, by = 0.5)),
            boxplot = TRUE,
            filename = "~/ep/fixed/epilepsy_fixed",
            confidence.level = 0.99, 
            num.avg = 10) {

   T <- list(read.table("~/ep/random/epilepsy_random_1.res")[,col])
   T <- c(T, list(read.table("~/ep/random/epilepsy_random_4.res")[,col]))
 
   for (i in 1:length(freqs)) {
      v <- read.table(ps(filename,"_", freqs[i],".res"))[,col]
      T <- c(T, list(v))      
      }
   
      T <- c(T, list(read.table("~/ep/features/no_penalty.res")[,col]))      
      T <- c(T, list(read.table("~/ep/features/small_penalty.res")[,col]))      
      T <- c(T, list(read.table("~/ep/features/best.res")[,col]))      
      
     if (boxplot) boxplot(T, names = c("rand1", "rand2", paste(freqs), "kbsf5", "kbsf6", "kbsf7"), range = 0)
     }
        


generate.all.epilepsy <- function() {
   
   x11(w = 8, h = 5)
   show.ep.no.features(col=2)
   dev.copy2eps(file = "~/tex/canada/fig/epilepsy_no_features.eps")
   
   show.ep.no.features(col=3)
   dev.copy2eps(file = "~/tex/canada/fig/epilepsy_no_features_stimulation.eps")

   show.ep.features(col=2)
   dev.copy2eps(file = "~/tex/canada/fig/epilepsy_features.eps")
   
   show.ep.features(col=3)
   dev.copy2eps(file = "~/tex/canada/fig/epilepsy_features_stimulation.eps")
   
   x11(w = 5, h = 5)
   show.ep.features.pareto()
   dev.copy2eps(file = "~/tex/canada/fig/epilepsy_features_pareto.eps")
   }
   
   



show.ep2 <- function(col = 1, freqs = c(0.005, seq(0.5, 2, by = 0.5)),
            boxplot = TRUE,
            confidence.level = 0.99, 
            filename = "~/ep/experiments/epilepsy_fixed") {
   T <- apply(read.table("~/ep/experiments/epilepsy_random.res"), 2, mean)[2:3]
   T <- rbind(T,apply(read.table("~/ep/experiments/epilepsy_random_is4.res"), 2, mean)[2:3])
 
   for (i in 1:length(freqs)) {
      v <- apply(read.table(ps(filename,"_", freqs[i],".res")), 2, mean)[2:3]
      T <- rbind(T,v)
      }
#    T[,1] <- T[,1] / max(T[,1])   
   
   T <- rbind(T, apply(read.table("~/ep/kbsf_results.res"), 2, mean)[2:3])
   plot(T[,2],T[,1])
   text(T[,2]-0.005,T[,1]-0.005, c("random", "rand_fr", paste(freqs), "kbsf"))
   }





print("dp.data.analysis.R loaded")   
