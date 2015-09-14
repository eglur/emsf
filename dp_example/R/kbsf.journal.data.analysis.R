source("util.R")
source("data.plot.R")
source("dp.data.analysis.R")
source("online.rl.data.analysis.R")


generate.all <- function(dir = "~/tmp/")  {
#    x11(w=7, h=7)
   set.par.new()

   ## PUDDLE ##
   show.pd.m(setpar = FALSE)
   dev.copy2pdf(file = ps(dir, "puddle_m.pdf"), font="Bookman")
   
   show.pd.m(col = 4, log="y", ylab = "Seconds (log)", ylim=NULL, leg.pos = "left", setpar = FALSE)
   dev.copy2pdf(file = ps(dir, "puddle_m_time_log.pdf"), font="Bookman")

   show.pd.m(col = 4, transf.func = id, ylab = "Seconds", ylim=NULL, leg.pos = "left", setpar = FALSE)
   dev.copy2pdf(file = ps(dir, "puddle_m_time.pdf"), font="Bookman")


   show.pd.n(setpar = FALSE)
   dev.copy2pdf(file = ps(dir, "puddle_n.pdf"), font="Bookman")
   
   show.pd.n(col = 4, ylab = "Seconds (log)", ylim=NULL, leg.pos="right", log="y", setpar = FALSE)
   dev.copy2pdf(file = ps(dir, "puddle_n_time_log.pdf"), font="Bookman")

   show.pd.n(col = 4, ylab = "Seconds", ylim=NULL, leg.pos="topleft", setpar = FALSE)
   dev.copy2pdf(file = ps(dir, "puddle_n_time.pdf"), font="Bookman")


   ## POLE ##
   show.pl.m(setpar = FALSE)
   dev.copy2pdf(file = ps(dir, "pole_m.pdf"), font="Bookman")
   
   show.pl.m(col = 4, ylab = "Seconds", ylim = NULL, leg.pos = "topleft", setpar = FALSE)
   dev.copy2pdf(file = ps(dir, "pole_m_time.pdf"), font="Bookman")
   
   show.pl.m(col = 4, ylab = "Seconds (log)", ylim = NULL, leg.pos = "bottomright", log="y", y.inter = 1, inset = c(0.0, 0.0), setpar = FALSE)
   dev.copy2pdf(file = ps(dir, "pole_m_time_log.pdf"), font="Bookman")

## DOUBLE POLE ##
   show.dpl.m(setpar = FALSE)
   dev.copy2pdf(file = ps(dir, "double_pole_m.pdf"), font="Bookman")
   
   show.dpl.m(col = 4, ylab = "Seconds", ylim = NULL, leg.pos = "topleft", setpar = FALSE)
   dev.copy2pdf(file = ps(dir, "double_pole_m_time.pdf"), font="Bookman")

   show.dpl.m(col = 4, ylab = "Seconds (log)", ylim = NULL, leg.pos = "bottomright", log="y", y.inter = 1, inset = c(0.0, 0.0), setpar = FALSE)
   dev.copy2pdf(file = ps(dir, "double_pole_m_time_log.pdf"), font="Bookman")


## HIV ##
   plot.hiv.m()
   dev.copy2pdf(file = ps(dir, "hiv_return_m.pdf"), font="Bookman")

   plot.hiv.m(2, ylim=NULL, ylab = "Seconds", delta.y.fqit = 0.02, delta.y.kbsf = 0.8, delta.x.kbsf = 60)
   dev.copy2pdf(file = ps(dir, "hiv_time_m.pdf"), font="Bookman")

   plot.hiv.nag()
   dev.copy2pdf(file = ps(dir, "hiv_return.pdf"), font="Bookman")

   plot.hiv.nag(2, ylab = "Seconds", ylim = NULL, delta.y.fqit = 0.03, delta.x.kbsf = 0.4)
   dev.copy2pdf(file = ps(dir, "hiv_time.pdf"), font="Bookman")

## TRIPLE POLE ##

   show.triple.pole()
   dev.copy2pdf(file = ps(dir, "triple_pole_episodes.pdf"), font="Bookman")

   show.triple.pole(suf = "tim", ylab = "Seconds (log)", log = "y" , ylim = c(40, 100000), leg.pos = "bottomright")
   dev.copy2pdf(file = ps(dir, "triple_pole_times.pdf"), font="Bookman")

   show.triple.pole.batch()
   dev.copy2pdf(file = ps(dir, "triple_pole_episodes_egreedy.pdf"), font="Bookman")

   show.triple.pole.batch("tim", ylab = "Seconds (log)", log = "y", leg.pos = "bottomright", ylim = c(40, 100000))
   dev.copy2pdf(file = ps(dir, "triple_pole_times_egreedy.pdf"), font="Bookman")

## HELICOPTER ##

   plot.helicopter()
   dev.copy2pdf(file = ps(dir, "helicopter_step.pdf"), font="Bookman")

   plot.helicopter("time", ylab = "Seconds", new.window = FALSE)
   dev.copy2pdf(file = ps(dir, "helicopter_time.pdf"), font="Bookman")

   plot.helicopter("time_step", ylab = "Seconds", new.window = FALSE)
   dev.copy2pdf(file = ps(dir, "helicopter_time_step.pdf"), font="Bookman")

   }



set.par.new <- function() {
# set parameters used  in the plots
   par(
      lwd = 1.5,
#       mar = c(5, 6, 2, 2) -1 + 0.1,
      mar = c(3, 5, 2, 2) -1 + 0.1,
      ps=12, cex = 1.5,
      omd = c(0, 1, 0, 1) )
   }


show.pd.rs.methods <- function(size=8000, ms = (4:13)^2, col = 1, 
                       ylim = c(0.3,3.2),
                       taus = c(0.1), taus2 = c(0.1), show.log = FALSE, 
                        ylab = "Return", leg.pos = "bottomright", 
                        base.dir = "~/fontes/R/kbsf/puddle/", 
                        nm = 1, 
                        rs.methods = c("k-means", "k-centers", "grid",  "random"), 
                        leg.text = c("K-means", "K-centers", "Evenly distributed", "Random"), ...) 
{
   
   T <- list(y = NULL, ub = NULL, lb = NULL)
   
   for (i in 1:length(rs.methods))
   { 
      
      print(rs.methods[i])
      K <- gen.data(paste(base.dir, rs.methods[i], "/puddle_kbsf", sep=""), sizes = size, taus = taus, 
                 taus2 = taus2, col = col, pick.func=which.max, ms = ms,
                 transf.func = id, nm = nm)

      T$y <-  cbind(T$y, K$y)
      T$ub <- cbind( T$ub, K$ub )
      T$lb <- cbind( T$lb, K$lb )
   }

   
   if (show.log) {
      T$y <-  log(T$y)
      T$ub <- log(T$ub)
      T$lb <- log(T$lb)
      }

   plot.algs(ms, Y = T$y, Y.lb = T$lb, Y.ub = T$ub, ylim=ylim, 
             plot.leg = TRUE, ylab = ylab, 
             xlab = "m", leg = leg.text, leg.pos = leg.pos, ...)
   }


show.pl.trees.kbsf <- function(col = 3, n = 50000, nt = 30, mne = c(20,200,2000),
                         ncp = 4, pr = 0.6, m = 100, taud = c(1, 0.1, 0.01), 
                         tauk = 1, dir="./kbsf/pole/", confidence.level = 0.99,
                         ylim=c(0,1), ...)
   {
   
   K <- gen.data(ps(dir,"pole_kbsf"), sizes = n, taus = taud, taus2 = tauk, col = col, ms = m)
   
   x    <- mne
   y    <- matrix(0, length(x), 1)
   y.ub <- y
   y.lb <- y
         
   for (i in 1:nrow(y))
   {
      fn <- ps(dir,"pole_tree_",paste(n,nt,mne[i],ncp,pr,sep="_"),".res")
      D <- as.matrix(read.table(fn)[,col])
      y[i,1] <- mean(D)
      ci <- qnorm( 1-(1-confidence.level)/2) * sd(D) / sqrt(length(D))
      y.ub[i,1] <- y[i,1] + ci
      y.lb[i,1] <- y[i,1] - ci
   }
   print(y)
   
   y <- cbind(rep(K$y,nrow(y)), y)
   print(y)
   y.ub <- cbind(rep(K$ub,nrow(y)), y.ub)
   y.lb <- cbind(rep(K$lb,nrow(y)), y.lb)
   
   mp(x,y,y.lb,y.ub, ylim = ylim, ...)
  
   }


show.dpl.trees.kbsf <- function(col = 3, n = 50000, nt = 30, mne = c(20,200,2000),
                         ncp = 4, pr = 0.6, m = 100, taud = c(1, 0.1, 0.01), 
                         tauk = 1, dir="./kbsf/pole/", confidence.level = 0.99)
   {
   
   K <- gen.data(ps(dir,"pole_kbsf"), sizes = n, taus = taud, taus2 = tauk, col = col, ms = m)
   
   x    <- mne
   y    <- matrix(0, length(x), 1)
   y.ub <- y
   y.lb <- y
         
   for (i in 1:nrow(y))
   {
      fn <- ps(dir,"pole_tree_",paste(n,nt,mne[i],ncp,pr,sep="_"),".res")
      D <- as.matrix(read.table(fn)[,col])
      y[i,1] <- mean(D)
      ci <- qnorm( 1-(1-confidence.level)/2) * sd(D) / sqrt(length(D))
      y.ub[i,1] <- y[i,1] + ci
      y.lb[i,1] <- y[i,1] - ci
   }
   print(y)
   
   y <- cbind(rep(K$y,nrow(y)), y)
   print(y)
   y.ub <- cbind(rep(K$ub,nrow(y)), y.ub)
   y.lb <- cbind(rep(K$lb,nrow(y)), y.lb)
   
   mp(x,y,y.lb,y.ub)
  
   }   
                         

show.pd.bagging <- function(col = 1, n = 8000, m = seq(10,100,by=20), nm = c(1,2, 5,10),
                     taud = 0.1, tauk = 0.1, 
                     dir = "./kbsf/puddle/bagging/",
                     ylab = "Return"
                           )
{
   X <- matrix(0, length(m), length(nm))
   for (i in 1:length(m))
   {
      for (j in 1:length(nm))
      {
         fn <- ps(dir, paste("puddle_kbsf", n, nm[j], m[i], taud, tauk, sep="_"), ".res")
         D <- read.table(fn)
         X[i,j] <- mean(D[,col])
      }
   }
   
   set.par()
   matplot(m, X, t="o", ylab=ylab)
}


plot.hiv.m <- function(
   col = 1, # time = col 2 
   tree.filename = "hiv_tree",
   nt = 30,
   mne = c(50, 100, 200),
   ncp = 8,
   it = 50,
   nfi = 10,
   kbsf.filename = "hiv_kbsf_batch",
   m = c(seq(2000,10000,by=2000)),
   td = 1,
   tk = 1,
   mmx = -1,
   tmd = 0,
   nnd = 3,
   nnk = 2,
   nag = 1,
   sdt = 0.1,
   confidence.level = 0.99,
   dir = "./kbsf/hiv/",
   ylab = "Return", 
   ylim=c(0,1e9), # to guarantee that it looks like nag
   delta.y.fqit = 0.04, 
   delta.y.kbsf  = 0.15,
   delta.x.kbsf = 120,
   ...
   )
{
   
   x    <- m
   y    <- matrix(0, length(x), length(mne))
   y.ub <- y
   y.lb <- y
         
   for (i in 1:length(mne))
   {
      fn <- ps(dir,tree.filename,"_",paste(nt,mne[i],ncp,it,nfi,sep="_"),".res")
      D <- as.matrix(read.table(fn)[,col])
      y[,i] <- rep(mean(D), nrow(y))
      ci <- qnorm( 1-(1-confidence.level)/2) * sd(D) / sqrt(length(D))
      y.ub[,i] <- y[,i] + ci
      y.lb[,i] <- y[,i] - ci
   }
   
   z    <- matrix(0, length(x), 1)
   z.ub <- z
   z.lb <- z

   for (i in 1:length(m))
   {
      fn <- ps(dir,kbsf.filename, "_", paste(m[i], td, tk, mmx, tmd, nnd, nnk, nag, sdt, sep = "_"),".res")
      D <-  as.matrix(read.table(fn))[,col]
      j <- 1
      z[i,j] <- mean(D)
      ci <- qnorm( 1-(1-confidence.level)/2) * sd(D) / sqrt(length(D))
      z.ub[i,j] <- z[i,j] + ci
      z.lb[i,j] <- z[i,j] - ci
   }
   
   if (is.null(ylim)) ylim = c(min(min(y.lb),min(z.lb)), max(max(y.ub),max(z.ub)))
#    par(mar = c(6, 6, 3, 2) -1 + 0.1)
   mp(x, y, y.ub, y.lb, t="l", plot.error = TRUE, show.error = FALSE, lty = 1, col="DARKGREEN", lwd=2, xlab="m", ylab = ylab, ylim = ylim, ...)
   mp(x, z, z.ub, z.lb, col="BLUE", lwd=2, add = TRUE, pch = 23, ...)
   text(m[length(m)%%2]-100,y[length(m)%%2,] + delta.y.fqit * y.ub[1,1] ,paste("FQIT(",mne,")", sep=""), pos=4) 
   text(m[length(m)-2] + delta.x.kbsf, z[length(m),1] + delta.y.kbsf  *  z[length(m),1],"KBSF(m)", pos=4)
}



plot.hiv.nag <- function(
   col = 1,
   tree.filename = "hiv_tree",
   nt = 30,
   mne = c(50, 100, 200),
   ncp = 8,
   it = 50,
   nfi = 10,
   kbsf.filename = "hiv_kbsf_batch",
   m = 10000,
   td = 1,
   tk = 1,
   mmx = -1,
   tmd = 0,
   nnd = 3,
   nnk = 2,
   nag = c(1,seq(5,30,by=5)),
   sdt = 0.1,
   confidence.level = 0.99,
   dir = "./kbsf/hiv/",
   ylab = "Return", 
   ylim = c(0,1e9),
   delta.y.fqit = - 0.05,
   delta.x.kbsf = 0.2, 
   ...
   )
{
   
   x    <- nag
   y    <- matrix(0, length(x), length(mne))
   y.ub <- y
   y.lb <- y
         
   for (i in 1:length(mne))
   {
      fn <- ps(dir,tree.filename,"_",paste(nt,mne[i],ncp,it,nfi,sep="_"),".res")
      D <- as.matrix(read.table(fn)[,col])
      y[,i] <- rep(mean(D), nrow(y))
      ci <- qnorm( 1-(1-confidence.level)/2) * sd(D) / sqrt(length(D))
      y.ub[,i] <- y[,i] + ci
      y.lb[,i] <- y[,i] - ci
   }
   
   z    <- matrix(0, length(x), 1)
   z.ub <- z
   z.lb <- z

   for (i in 1:length(nag))
   {
      fn <- ps(dir,kbsf.filename, "_", paste(m, td, tk, mmx, tmd, nnd, nnk, nag[i], sdt, sep = "_"),".res")
print(fn)
      D <-  as.matrix(read.table(fn))[,col]
      j <- 1
      z[i,j] <- mean(D)
      ci <- qnorm( 1-(1-confidence.level)/2) * sd(D) / sqrt(length(D))
      z.ub[i,j] <- z[i,j] + ci
      z.lb[i,j] <- z[i,j] - ci
   }
   
   if (is.null(ylim)) ylim = c(min(min(y.lb),min(z.lb)), max(max(y.ub),max(z.ub)))

#    par(mar = c(5, 5, 3, 2) -1 + 0.1)
   mp(x, y, y.ub, y.lb, t="l", plot.error = TRUE, show.error = FALSE, 
      lty = 1, col="DARKGREEN", lwd=2, xlab="Number of agents", 
      ylab = ylab, ylim = ylim, ...)
   
   if (col == 2)  ## THIS IS A TEMPORARY HACK TO FIX THE RUN TIME OF NAG=25
   {
      z[6,1] <- mean(c(z[5,1], z[7,1]))
      z.ub[6,1] <- mean(c(z.ub[5,1], z.ub[7,1]))
      z.lb[6,1] <- mean(c(z.lb[5,1], z.lb[7,1]))
   }
   
   mp(x, z, z.ub, z.lb, col="BLUE", lwd=2, add = TRUE, pch = 23, ...)
   text(x[length(x)-2]+2,y[1,] + delta.y.fqit * y.ub[1,1] ,paste("FQIT(", mne,")", sep=""), pos=4) 
   text(3,z[3,1] + delta.x.kbsf * z.ub[3,1] ,"KBSF(10000)", pos=4) 
}


helicopter.process.raw.data <- function(process.sarsa = TRUE, process.kbsf = TRUE, num.cols = 5000   )
   {
      if (process.sarsa)
      {
         print("Processing SARSA returns...")
         D <- read.table("./kbsf/helicopter_raw/ST_NCE_Returns.txt")
         D <- D[,1:20000]
         D <- reduce.cols(D, num.cols)
         wt(D, "./kbsf/helicopter/sarsa_returns.res")

         print("Processing SARSA steps...")
         D <- read.table("./kbsf/helicopter_raw/ST_NCE_Episodes.txt")
         D <- D[,1:20000]
         D <- reduce.cols(D, num.cols)
         wt(D, "./kbsf/helicopter/sarsa_steps.res")

         
         print("Processing SARSA times...")
         T <- read.table("./kbsf/helicopter_raw/ST_NCE_Time.txt")
         T <- T[,1:100000]
         T <- reduce.cols(T,num.cols)
         
         print("\t Computing average time per step")
         D2 <- D
         D2[D2==0] <- 1
         TA <- T/D2
         wt(TA, "./kbsf/helicopter/sarsa_time_step.res")

         print("\t Computing accumulated time")
         T <- sum.over.cols(T)
         wt(T, "./kbsf/helicopter/sarsa_time.res")
      }
      
      if (process.kbsf)
      {
         print("Processing KBSF returns...")
         D <- read.table("./kbsf/helicopter_raw/KBSF_NCE_Returns.txt")
         D <- D[,1:20000]
         D <- reduce.cols(D, num.cols)
         wt(D, "./kbsf/helicopter/kbsf_returns.res")

         print("Processing KBSF steps...")
         D <- read.table("./kbsf/helicopter_raw/KBSF_NCE_Episodes.txt")
         D <- D[,1:20000]
         D <- reduce.cols(D, num.cols)
         wt(D, "./kbsf/helicopter/kbsf_steps.res")

         
         print("Processing KBSF times...")
         T <- read.table("./kbsf/helicopter_raw/KBSF_NCE_Time.txt")
         T <- T[,1:100000]
         T <- reduce.cols(T,num.cols)
         
         print("\t Computing average time per step")
         D2 <- D
         D2[D2==0] <- 1
         TA <- T/D2
         wt(TA, "./kbsf/helicopter/kbsf_time_step.res")

         print("\t Computing accumulated time")
         T <- sum.over.cols(T)
         wt(T, "./kbsf/helicopter/kbsf_time.res")
      }

      
   }



plot.helicopter <- function
   (
   suf = "steps",
   filenames = c("sarsa", "kbsf"),
   dir = "./kbsf/helicopter/",
   confidence.level = 0.99,
    ylab = "Steps", 
    num.cols = 5000,
    new.window = TRUE,
    ...
   )
{
#    if (new.window) x11(w=14, h=4)

   D <- read.table(ps(dir,filenames[1],"_", suf,".res"))
#    D <- reduce.cols(D, num.cols)
   Y <- matrix(apply(D, 2, mean), ncol(D), 1)
#    Y <- matrix(D[sample(1:nrow(D),1),],ncol(D), 1)
   ci <- qnorm( 1-(1-confidence.level)/2) * apply(D,2,sd) / sqrt(nrow(D))
   Y.ub <- Y + ci
   Y.lb <- Y - ci

   D <- read.table(ps(dir,filenames[2],"_", suf,".res"))
#    D <- reduce.cols(D, num.cols)
  Y <- cbind(Y, apply(D, 2, mean))
#    Y <- cbind(Y, matrix(D[sample(1:nrow(D),1),], ncol(D), 1))
   ci <- qnorm( 1-(1-confidence.level)/2) * apply(D,2,sd) / sqrt(50) ## FIX WHEN RYAN FINISH RUNS
   Y.ub <- cbind(Y.ub, Y[,2] + ci)
   Y.lb <- cbind(Y.lb, Y[,2] - ci)
   set.par()
   mp(seq(1,100000, length=nrow(Y)), Y, Y.ub, Y.lb, show.error = FALSE, pch = NULL, 
      t="l", col=c("BLACK", "BLUE"), xlab="Episodes", ylab=ylab,...)
   leg(leg = c("SARSA", "iKBSF"), col=c("BLACK", "BLUE"), pch="")
# #    Y <- cbind(D, read.table(ps(dir,filenames[2])))
   
}



make.leg.taus <- function(taus, taus2) {
# makes a legend with labels sigma = sds[1, 2, ...]
   l <- expression()
   for(t in taus) {
      f <- substitute(expression(tau == ta), list(ta=t))[[2]]
      for(t2 in taus2) {
         f2 <- substitute(expression(bar(tau) == ta2), list(ta2=t2))[[2]]
         l <- c(l, bquote(.(f)*" , " *.(f2)))
         print(l)
         }
      }
   l
   }
   
   
show.pd.all.kbsf <- function(n=8000, ms = seq(10,150,by=20), col = 1, ylim = NULL,
                       taus = c(1, 0.1, 0.01), taus2 = c(1, 0.1, 0.01), ylab = "Return",
                       leg.pos = "bottomright", dir = "~/fontes/R/kbsf/puddle/",
                       filename = "puddle_kbsf", 
                       nm = NULL, new.window = FALSE, ...) 
{
   D <- matrix(0, length(ms), length(taus) * length(taus2))
   
   for (i in 1:length(ms))
   {
      for (j in 1:length(taus))
      {
         for (k in 1:length(taus2))
         {
            fn <- ps(dir,paste(filename, n, ms[i], taus[j], taus2[k], sep="_"), ".res")
            T <- as.matrix(read.table(fn))[,col]
            D[i,(j-1)*length(taus) + k] <- mean(T)
         }
      }
   }
   
   K <- matrix(0, length(ms), 1)
   for (i in 1:length(ms))
   {
      fn <- ps(dir,paste(filename, n, ms[i], 0.1, 0.1, sep="_"), ".res")
      K[i,1] <- mean(as.matrix(read.table(fn))[,col])
   }

   if (new.window) x11(w=10, h=7)
   set.par.new()

   par(mar = c(5, 5, 3, 11) - 1 + 0.1)
   
   matplot(ms, D, t="o", ylab="Return", col=rgb(0.5,0.5,0.5), lty=1, pch=1:ncol(D), xlab="m")
   mp(ms, K, add=TRUE, col="blue")
   
   par(xpd=TRUE)
   set.par() # set.par.new increases the dist between lines
   par(cex=1)
   leg <- make.leg.taus(taus, taus2)
   colors <- c(rep(rgb(0.5,0.5,0.5),4), "BLUE", rep(rgb(0.5,0.5,0.5),4))
   legend(160,3, leg=leg, bty = "n", lty=1, 
          pch=c(1:4, 23, 5:ncol(D)), cex=1, y.intersp = 1.5, col=colors, pt.bg = colors)
}   
   


show.pd.kbrl.tau <- function(n=1000, col = 1, ylim = NULL,
                       taus = round(10^seq(-1.3333333,1,l=8), d=2), ylab = "Return",
                       leg.pos = "bottomright", dir = "~/dp_trans/files/",
                       filename = "puddle_kbrl", 
                       nm = 1, new.window = FALSE, 
                       confidence.level = 0.99, ...) 
{
   D <- matrix(0, length(taus), 1)
   D.lb <- D
   D.ub <- D

   for (i in 1:length(taus))
   {
         fn <- ps(dir,paste(filename, n, nm, taus[i], sep="_"), ".res")
         T <- as.matrix(read.table(fn))[,col]
         D[i,1] <- mean(T)
         ci <- qnorm( 1-(1-confidence.level)/2) * sd(T) / sqrt(length(T)) 
         D.ub[i,1] <- D[i,1] + ci
         D.lb[i,1] <- D[i,1] - ci
   }
   

   mp(taus, D, D.ub, D.lb, t="o", ylab="Return", xlab=expression(tau* "(log)") , log="x", xaxt = "n", ylim=c(-1, 2.5), ...)
   x <- c(0.1, 1, 10)
   axis(1, at = x, labels = paste(round(x, d = 2)))
   leg("topright", "KBRL(1000)")
}   
   




print("kbsf.journal.data.analysis.R loaded")   
