source("data.plot.R")

## PUDDLE WORLD
show.puddle.n.new <- function(
                     suf = "ret", 
                     nk = 2000,
                     n = 8000, 
                     m = 100, 
                     td = 0.01, 
                     tk = 0.01, 
                     nnd = -1,
                     nnk = -1,
                     mmx = -1,
                     thr = 0,
                     umi = c(1000, 2000, 4000, 8000), 
                     uvi = c(1000, 2000, 4000, 8000), 
                     ep = 1, 
                     filename = "~/fontes/c++/online_rl/files/puddle_kbsf",
                     ylab = "Return",
                     confidence.level = 0.99,
                     plot.error = TRUE, 
                     ylim = NULL, 
                     num.avg = 50, 
                     log = "")

{
  
  set.par()
  Y <- NULL
  Y.ui <- NULL
  Y.li <- NULL
  
  for (i in 1:length(umi))
  {
    fn <- paste(filename, fo(nk), fo(n), m, td, tk, nnd, nnk, mmx, thr, fo(umi[i]), fo(uvi[i]), ep, suf, sep="_")
    fn <- ps(fn, ".res")
    R <- as.matrix(read.table(fn))

    ci <- qnorm(1-(1-confidence.level)/2) * apply(R,2,sd) / sqrt(nrow(R))
    M <- apply(R, 2, mean)
    Y    <- cbind(Y, M)
    Y.ui <- cbind(Y.ui, M + ci)
    Y.li <- cbind(Y.li, M - ci)
    print(i)
  }
 
  mp(seq(1, n, length=nrow(Y)), Y, Y.li = Y.li, Y.ui = Y.ui, lwd = 3, ylab = ylab, xlab = "Number of sample transitions",
     num.avg = num.avg, plot.error = plot.error, log = log, ylim = ylim, show.error = rep(FALSE,ncol(Y)))
  leg("topleft", leg=make.leg.iota(umi), lwd = 3)
}


show.puddle.n <- function(
                    suf = "ret", 
                    n = 8000, 
                    m = 100, 
                    td = 0.1, 
                    tk = 0.01,
                    umi = c(1000, 2000, 4000, 8000),
                    uvi = c(1000, 2000, 4000, 8000), 
                    ep = 1, 
                    si = 500, 
                    ylab = "Return", 
                    filename = "~/fontes/c++/online_rl/files/puddle_kbsf",
                    confidence.level = 0.99,
                     plot.error = TRUE, ylim = NULL, num.avg = 50, log = "")

{
  
  set.par()
  Y <- NULL
  Y.ui <- NULL
  Y.li <- NULL
  
  for (i in 1:length(umi))
  {
    fn <- paste(filename, n, m, td, tk, umi[i], uvi[i], ep, suf, sep="_")
    fn <- ps(fn, ".res")
    R <- as.matrix(read.table(fn))

    ci <- qnorm(1-(1-confidence.level)/2) * apply(R,2,sd) / sqrt(nrow(R))
    M <- apply(R, 2, mean)
    Y    <- cbind(Y, M)
    Y.ui <- cbind(Y.ui, M + ci)
    Y.li <- cbind(Y.li, M - ci)
    print(i)
  }
 
  mp(seq(1, n, length=nrow(Y)), Y, Y.li = Y.li, Y.ui = Y.ui, lwd = 3, ylab = ylab, xlab = "Number of sample transitions",
     num.avg = num.avg, plot.error = plot.error, log = log, ylim = ylim, show.error = rep(FALSE,ncol(Y)))
  leg("topleft", leg=make.leg.iota(umi), lwd = 3)
}


show.puddle.ep <- function(
                     suf = "ret", 
                     n = 8000, 
                     m = 100, 
                     td = 0.1, 
                     tk = 0.01, 
                     umi = 1000, 
                     uvi = 1000, 
                     ep = c(1, 0.75, 0.5, 0.25,0), 
                     filename = "~/fontes/c++/online_rl/files/puddle_kbsf",
                     ylab = "Return",
                     confidence.level = 0.99,
                     plot.error = TRUE, 
                     ylim = NULL, 
                     num.avg = 50, 
                     log = "")

{
  
  set.par()
  Y <- NULL
  Y.ui <- NULL
  Y.li <- NULL
  for (i in 1:length(ep))
  {
    fn <- paste(filename, n, m, td, tk, umi, uvi, ep[i], suf, sep="_")
    fn <- ps(fn, ".res")
    R <- as.matrix(read.table(fn))

    ci <- qnorm(1-(1-confidence.level)/2) * apply(R,2,sd) / sqrt(nrow(R))
    M <- apply(R, 2, mean)
    Y    <- cbind(Y, M)
    Y.ui <- cbind(Y.ui, M + ci)
    Y.li <- cbind(Y.li, M - ci)
    print(i)
  }
 
  mp(seq(1, n, length=nrow(Y)), Y, Y.li = Y.li, Y.ui = Y.ui, lwd = 3, ylab = ylab, xlab = "Number of sample transitions",
     num.avg = num.avg, plot.error = plot.error, log = log, ylim = ylim, show.error = rep(FALSE,ncol(Y)))
  leg("topleft", leg=make.leg.epsilon(ep), lwd = 3)
}





show.puddle.thr <- function(suf = "ret", 
                 nk = 500,
                 n = 7500, 
                 m = 30, 
                 td = 0.1, 
                 tk = 0.1,
                 mmx = 100,
                 thr = c(0.1, 0.3, 0.5, 0.7),
                 umi = 1000,
                 uvi = 1000, 
                 ep = 0.25, 
                 si = 1000, 
                 ylab = "Return", 
                 leg.position = "topright",
                 filename = "~/fontes/c++/online_rl/files/puddle_kbsf",
                 confidence.level = 0.99,
                  plot.error = TRUE, ylim = NULL, num.avg = 50, log = "")
{
  
  set.par()
  Y <- NULL
  Y.li <- NULL
  Y.ui <- NULL
  for (i in 1:length(thr))
  {
    fn <- paste(filename, nk, n, m, td, tk, mmx, thr[i], umi, uvi, ep, suf, sep="_")
    fn <- ps(fn, ".res")
    R <- as.matrix(read.table(fn))

    ci <- qnorm(1-(1-confidence.level)/2) * apply(R,2,sd) / sqrt(nrow(R))
    M <- apply(R, 2, mean)
    Y    <- cbind(Y, M)
    Y.ui <- cbind(Y.ui, M + ci)
    Y.li <- cbind(Y.li, M - ci)
    print(i)
  }
 
  mp(seq(1, n, length=nrow(Y)), Y, Y.li = Y.li, Y.ui = Y.ui, lwd = 3, ylab = ylab, xlab = "Number of sample transitions",
     num.avg = num.avg, plot.error = plot.error, log = log, ylim = ylim, show.error = rep(FALSE,ncol(Y)))
  leg("topleft", leg=make.leg.epsilon(ep), lwd = 3)
}

fo <- function(n) 
{
   n.sx <- n
   if (n >= 1e5 && n < 1e6) n.sx <- sprintf("%5.0f", n)
    else if (n >= 1e6 && n < 1e7) n.sx <- sprintf("%6.0f", n)
    else if (n >= 1e7 && n < 1e8) n.sx <- sprintf("%7.0f", n)
    n.sx
}

show.double.pole.m <- function(suf = "eps", n = 500000, ms = seq(200,500,by=100),
		    td = 1, tk = 10,
		    umi = 100000, uvi = 100000, 
		    ep = 1, 
		    filename = "~/fontes/c++/online_rl/files/double_pole_kbsf",
		    ylab = "Successful episodes")
{
  
  set.par()
  T <- NULL
  for (i in 1:length(ms))
  {
 
    fn <- paste(filename, fo(n), ms[i], td, tk, fo(umi), fo(uvi), ep, suf, sep="_")
    fn <- ps(fn, ".res")
    R <- as.matrix(read.table(fn))
    T <- cbind(T, apply(as.matrix(read.table(fn)), 2, mean))
  }
  mp(seq(1,n + 500000,length=nrow(T)), T, lwd = 3, ylab = ylab, xlab = "n")
  leg("right", leg=paste("m=", ms), lwd = 3)
}


show.double.pole <- function(suf = "stp", nk = 500000, n = 500000, ms = seq(200,500,by=100),
              td = 1, tk = 10,
              umi = 100000, uvi = 100000, 
              ep = 1, 
              confidence.level = 0.99, 
              filename = "~/fontes/c++/online_rl/files/double_pole_kbsf",
              ylab = "Successful episodes",
              log = "",
              leg.position = "topright",
              plot.error = TRUE,
              ylim = NULL
              )
{
  
  set.par()
  Y    <- NULL
  Y.ui <- NULL
  Y.li <- NULL
  for (i in 1:length(ms))
  {
 
    fn <- paste(filename, fo(nk), fo(n), ms[i], td, tk, fo(umi), fo(uvi), ep, suf, sep="_")
    fn <- ps(fn, ".res")
    R  <- as.matrix(read.table(fn))
    ci <- qnorm(1-(1-confidence.level)/2) * apply(R,2,sd) / sqrt(nrow(R))
    M <- apply(R, 2, mean)
    Y    <- cbind(Y, M)
    Y.ui <- cbind(Y.ui, M + ci)
    Y.li <- cbind(Y.li, M - ci)
  }
  
  mp(seq(nk, n + nk, length=nrow(Y)), Y, Y.li = Y.li, Y.ui = Y.ui, lwd = 3, ylab = ylab, xlab = "n",
     num.avg = num.avg, plot.error = plot.error, log = log, ylim = ylim)
  leg(leg.position, leg=paste("m=", ms), lwd = 3)
}


rename.double.pole <- function(suf = "eps", n = 500000, ms = seq(200,500,by=100),
              td = 1, tk = 10,
              umi = 100000, uvi = 100000, 
              ep = 1, 
              filename = "~/fontes/c++/online_rl/files/double_pole_kbsf",
              ylab = "Successful episodes")
{
  
  for (i in 1:length(ms))
   {
 
    fn <- paste(filename, fo(n), ms[i], td, tk, fo(umi), fo(uvi), ep, suf, sep="_")
    fn <- ps(fn, ".res")
    R <- as.matrix(read.table(fn))
    wt(R, ps(paste(filename, fo(n), fo(n), ms[i], td, tk, fo(umi), fo(uvi), ep, suf, sep="_"), ".res"))
  }
}



show.triple.pole <- function(suf = "eps", 
                             nk = 1000000, 
                             n = 10000000, 
                             m = 200,
                             td = 1, 
                             tk = 100,
                             nnd = 10,
                             nnk = 50 ,
                             mmx = -1,
                             thr = 0.01,
                             umi = 1000000, 
                             uvi = 1000000, 
                             si = 1000000, 
                             ep = 0.3, 
                             confidence.level = 0.99, 
                             filename.kbsf=  "~/fontes/R/kbsf/triple_pole/triple_pole_kbsf_1000000_9000000_200_1_100_10_50_-1_0.01_1000000_1000000_0.3",
                             filename.tree = c(
                                               "~/fontes/c++/dp_trans/files/triple_pole_tree_rl_1000000_30_10000_8_0.6.res",
                                               "~/fontes/c++/dp_trans/files/triple_pole_tree_rl_1000000_30_1000_8_0.6.res",
                                               "~/fontes/c++/dp_trans/files/triple_pole_tree_rl_1000000_30_100_8_0.6_merged.res"),
                             ylab = "Successful episodes",
                             log = "",
                             leg.position = "topleft",
                             show.leg = TRUE,
                             plot.error = TRUE,
                              ylim = c(0.1,1), # to make sure it matches the batch version of the plot
                             plot.tree = TRUE, 
                              ...
                              )
{
  
  Y    <- matrix(0, n/uvi, 0)
  Y.ui <- Y
  Y.li <- Y
 
  
  if (plot.tree)
  {
      col.tree <- NULL
      if (suf == "ret") col.tree <- 1
      if (suf == "stp") col.tree <- 2
      if (suf == "eps") col.tree <- 3
      if (suf == "tim") col.tree <- 4
      
      
      for (i in 1:length(filename.tree))
      {
         R <- read.table(filename.tree[i])
         print(paste(filename.tree[i], nrow(R)))
         ci <- qnorm(1-(1-confidence.level)/2) * sd(R[,col.tree]) / sqrt(nrow(R))
         m <- mean(R[,col.tree])
         Y <- cbind(Y,m)
#          if (i == 3) ci <- 0
         Y.ui <- cbind(Y.ui,m + ci) 
         Y.li <- cbind(Y.li, m-ci)
      }

  }
  
   
    fn <- ps(filename.kbsf,"_",suf,".res")
    R  <- as.matrix(read.table(fn))
    print(paste(fn, nrow(R)))
    ci <- qnorm(1-(1-confidence.level)/2) * apply(R,2,sd) / sqrt(nrow(R))
    M <- apply(R, 2, mean)
    Y    <- cbind(Y, M)
    Y.ui <- cbind(Y.ui, M + ci)
    Y.li <- cbind(Y.li, M - ci)


  show.error <- rep(FALSE, 1 + length(filename.tree))
  pch <- c(rep(NA, length(filename.tree)), 23)
# This is for log
#   Y[1,1] <- 0.1
#   Y.li[1,1] <- 0.1
#   Y.ui[1,1] <- 0.1
  mp(seq(nk, n, length=nrow(Y)), Y, Y.li = Y.li, Y.ui = Y.ui, lwd = 3, ylab = ylab, xlab = "Number of sample transitions",
     num.avg = num.avg, plot.error = plot.error, log = log, ylim = ylim, show.error = show.error, pch=pch, 
    col = c("BLACK", "RED", "BLUE"), ...)
  
  if (show.leg)
  {
#    p <- locator(1) #list(x = 1504297, y = 0.6391959)  #
    if (suf == "tim") p <- list(x = 3124415 + 50, y= 0.1153789 + 50) # TIME
    else if (suf == "eps") p <- list(x = 3 * nk , y= 0.323995) # EPISODES
    else p <- locator(1) 
    text(p, "Batch KBSF")
#   print(p)
  
#   p <- locator(1) # list(x = 1539517, y = 0.6125268) # 
#   text(p, "KBSF")
#   print(p)
  
#   p <- list(x = 1592347, y = 0.6039546) #locator(1)
  
  arrows(p$x - 1.2e6, p$y, nk , Y[1,3] - 0.01, col = "BLUE")
  print(p)
  
  
  leg(leg.position, leg=c("FQIT(1000)", "FQIT(100)", "iKBSF"), lwd = 3, pch=c(NA,NA,23))
  }
}

triple.pole.generate.all <- function(dir = "~/tmp/")
{
   # As in the NIPS 2012 paper
   x11(w=7, h=7)
   show.triple.pole()
   dev.copy2pdf(file = ps(dir,"/triple_pole_episodes.pdf"))
   show.triple.pole("tim", log="y", ylab="Seconds (log)", leg.pos="bottomright")
   dev.copy2pdf(file = ps(dir,"/triple_pole_times.pdf"))
   
   x11(w=4.5, h=7)
   show.triple.pole("nrs", plot.tree=FALSE, show.leg=FALSE, ylab="Number of representative states", col="BLUE")
   dev.copy2pdf(file = ps(dir,"/triple_pole_rep_states.pdf"))
}



show.triple.pole.batch <- function(
suf = "eps", 
filenames.tree =
c( 
"triple_pole_tree_10000000_30_10000_8_50_10_1000000_0.3",
"triple_pole_tree_10000000_30_1000_8_50_10_1000000_0.3"
),
filenames.kbsf.batch=
c(
"triple_pole_kbsf_10000000_200_1_100_10_50_-1_0_1000000_0.3",
"triple_pole_kbsf_10000000_3000_1_100_10_50_-1_0_1000000_0.3",
"triple_pole_kbsf_10000000_4000_1_100_10_50_-1_0_1000000_0.3",
"triple_pole_kbsf_10000000_5000_1_100_10_50_-1_0_1000000_0.3"
),
filenames.kbsf.inc =
c(
"triple_pole_kbsf_1000000_9000000_200_1_100_10_50_-1_0.01_1000000_1000000_0.3"
),
names = c("FQIT(10000)", "FQIT(1000)", "KBSF-200", "KBSF-3000", "KBSF-4000", "KBSF-5000", "iKBSF"),
plot.file = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE),
n = 10^7,
n.step = 10^6,
confidence.level = 0.99, 
ylab = "Successful episodes",
log = "",
leg.position = "bottomleft",
plot.error = TRUE,
ylim = c(0.1, 1), # to make sure it looks like the non-batch version
cols = NULL,
plot.tree = TRUE,
dir = "~/fontes/R/kbsf/triple_pole/"
)
{
  
  Y    <- NULL
  Y.ui <- NULL
  Y.li <- NULL
  
   # Tree batch #
   for (i in 1:length(filenames.tree))
   {
      fn <- ps(dir,filenames.tree[i],"_", suf, ".res")
      R  <- as.matrix(read.table(fn))
      ci <- qnorm(1-(1-confidence.level)/2) * apply(R,2,sd) / sqrt(nrow(R))
      M <- apply(R, 2, mean)
      Y    <- cbind(Y, M)
      Y.ui <- cbind(Y.ui, M + ci)
      Y.li <- cbind(Y.li, M - ci)
   }

   # KBSF batch #
   for (i in 1:length(filenames.kbsf.batch))
   {
      fn <- ps(dir,filenames.kbsf.batch[i],"_", suf, ".res")
      R  <- as.matrix(read.table(fn))
      ci <- qnorm(1-(1-confidence.level)/2) * apply(R,2,sd) / sqrt(nrow(R))
      M <- apply(R, 2, mean)
      Y    <- cbind(Y, M)
      Y.ui <- cbind(Y.ui, M + ci)
      Y.li <- cbind(Y.li, M - ci)
   }
   
   # KBSF incremental #
   for (i in 1:length(filenames.kbsf.inc))
   {
      fn <- ps(dir,filenames.kbsf.inc[i],"_", suf, ".res")
      R  <- as.matrix(read.table(fn))
      ci <- qnorm(1-(1-confidence.level)/2) * apply(R,2,sd) / sqrt(nrow(R))
      M <- apply(R, 2, mean)
      Y    <- cbind(Y, M)
      Y.ui <- cbind(Y.ui, M + ci)
      Y.li <- cbind(Y.li, M - ci)
   }
   
   # This is for log
   Y[Y == 0] <- 1e-5
   Y.li[Y.li == 0] <- 1e-5
   Y.ui[Y.ui == 0] <- 1e-5
   
   plot.cols <- (1:ncol(Y))[plot.file]
   Y <- Y[, plot.cols]
   Y.li <- Y.li[,plot.cols]
   Y.ui <- Y.ui[,plot.cols]
   total <- sum(plot.file)
   
  mp(seq(n.step, n, length=nrow(Y)), Y, Y.li = Y.li, Y.ui = Y.ui, lwd = 3, ylab = ylab, xlab = "Number of sample transitions",
     num.avg = num.avg, plot.error = plot.error, log = log, ylim = ylim, show.error=rep(FALSE, total),
     col = c("DARKGREEN", "BLACK", "BLUE"))
  
print(names[plot.file])  
leg(leg.position, leg=names[plot.file], lwd = 3, col = c("DARKGREEN", "BLACK", "BLUE"), pt.bg = c("DARKGREEN", "BLACK", "BLUE"))
}










triple.pole.merge.results.kbsf <- function(
   dir = "~/online_rl/files/",
   filename = "triple_pole_kbsf_1000000_10000000_200_1_100_10_50_-1_0.01_1000000_1000000_0.3_"
   )
{
   for (sf in c("eps", "tim", "nrs", "ort", "ret", "stp"))
   {
      T <- NULL
    
      fn <- ps(dir,filename, sf, ".res")
      T <- rbind(T,read.table(fn))

      fn <- ps(dir, "s2/", filename, sf,".res")
      T <- rbind(T,read.table(fn))
 
      fn <- ps(dir,"s3/", filename, sf, ".res")
      T <- rbind(T,read.table(fn))

      fn <- ps(dir,"s4/", filename, sf, ".res")
      T <- rbind(T,read.table(fn))

      if (nrow(T) > 100) T <- T[sample(1:nrow(T), 100), ]
 
      fn <- ps(dir, filename, sf, "_merged.res")
      
      wt(T, fn)
      print(paste("Saved", fn))
      pd(T)
   }
}
   
   
triple.pole.merge.results.tree <- function(
   dir = "~/fontes/c++/dp_trans/files/",
   filename = "triple_pole_tree_rl_1000000_30_100_8_0.6"
   )
{
      T <- NULL
    
      fn <- ps(dir,filename,".res")
      T <- rbind(T,read.table(fn))

      fn <- ps(dir, "s2/", filename, ".res")
      T <- rbind(T,read.table(fn))
 
#       fn <- ps(dir,"s3/", filename, ".res")
#       T <- rbind(T,read.table(fn))

      fn <- ps(dir, filename, "_merged.res")
      wt(T, fn)
      print(paste("Saved", fn))
      pd(T)
}
      
triple.pole.merge.results.tree2 <- function(
#    dir = "~/fontes/c++/dp_trans/files/",
#    filename = "triple_pole_tree_rl_1000000_30_10000_8_0.6"
    dir = "~/batch_online_rl/files/",
    filename = "triple_pole_tree_10000000_30_1000_8_50_10_1000000_0.3_tim"
   )
{
      fn <- ps(dir,filename,".res")
      T <- read.table(fn)

      for (i in 3)
      {   
         fn <- ps(dir, "s", i, "/", filename, ".res")
         T <- rbind(T,read.table(fn))
      }
      
      # correct negative times (problem with the timing function)
#       ind <- T[,4] > 0
#       T[!ind, 4] <- mean( T[ind,4] )

      if (nrow(T) > 100) T <- T[sample(1:nrow(T), 100),]

      fn <- ps(dir, filename, "_merged.res")
      wt(T, fn)
      print(paste("Saved", fn))
      pd(T)
}
      

triple.pole.merge.results.tree3 <- function(
    dir = "~/batch_online_rl/files/",
    filename = "triple_pole_tree_10000000_30_10000_8_50_10_1000000_0.3",
    suf = c("ret", "stp", "eps", "tim"),
    target.dir = "~/fontes/R/kbsf/inc/"                                          
   )
{
      for (i in 1:length(suf))
      {
         fn <- ps(dir,filename,"_", suf[i], ".res")
         T <- read.table(fn)
         
         fn <- ps(dir,"incomplete/", filename,"_", suf[i], ".res")
         T <- rbind(T, read.table(fn))

         fn <- ps(target.dir, filename, "_", suf[i], ".res")
         wt(T, fn)
         print(paste("Saved", fn))
         if (nrow(T) < 50) print("warning: less than 50 runs")
      }
}
      

      
compare.kbsf.triple.pole <- function(suf = "eps", 
files = c(
"triple_pole_kbsf_test_1000000_9000000_200_1_100_10_50_10000_0.01_1000000_1000000_0.3",  
"triple_pole_kbsf_test_1000000_9000000_200_1_100_9_9_10000_0.01_1000000_1000000_0.3",
"triple_pole_kbsf_test_1000000_9000000_200_1_100_20_20_10000_0.01_1000000_1000000_0.3",  
"triple_pole_kbsf_test_1000000_9000000_200_1_10_10_50_10000_0.01_1000000_1000000_0.3",
"triple_pole_kbsf_test_1000000_9000000_200_1_100_50_50_10000_0.01_1000000_1000000_0.3",
"./old/triple_pole_kbsf_1000000_10000000_200_1_100_10_50_-1_0.01_1000000_1000000_0.3"
),
dir = "~/online_rl/files/",
confidence.level = 0.99
)
{
   
   Y <- NULL
   Y.ui <- NULL
   Y.li <- NULL
   
   for (i in 1:length(files))
   {
      fn <- ps(dir,files[i],"_",suf,".res")
      R <- read.table(fn)
      if (ncol(R) > 10) R <- R[,1:10]
      if (nrow(R) > 50) R <- R[sample(1:nrow(R), 50), ]
      ci <- qnorm(1-(1-confidence.level)/2) * apply(R,2,sd) / sqrt(nrow(R))
      M <- apply(R, 2, mean)
      Y    <- cbind(Y, M)
      Y.ui <- cbind(Y.ui, M + ci)
      Y.li <- cbind(Y.li, M - ci)
      
   }
   
  mp(seq(10^6, 10^7, length=nrow(Y)), Y, Y.li = Y.li, Y.ui = Y.ui, lwd = 3, ylab = suf, xlab = "Number of sample transitions",
     num.avg = 50, plot.error = TRUE)
  leg("topleft", paste(seq(1:length(files))))
   
}
   

   
fix.kbsf.triple.pole <- function(suf = c("ret", "stp", "eps", "tim", "ort", "nrs"),
   file.source = "triple_pole_kbsf_1000000_10000000_200_1_100_10_50_-1_0.01_1000000_1000000_0.3",
   file.target = "triple_pole_kbsf_1000000_9000000_200_1_100_10_50_-1_0.01_1000000_1000000_0.3",
   dir = "~/online_rl/files/",
   dir2 = "old/"
   )

{
   for (i in 1:length(suf))
   { 
      fn  <- ps(dir,dir2,file.source,"_",suf[i],".res")
      T <- read.table(fn)
      T <- T[,1:10]
      fn2 <- ps(dir,file.target,"_",suf[i],".res")
      wt(T,fn2)
      print("From...")
      print(paste(fn))
      print("To...")
      print(paste(fn2))
   }
   
}

print("online.data.analysis.R loaded")