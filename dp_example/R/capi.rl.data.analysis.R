source("data.plot.R")


## POLE BALANCING
show.pole <- function(
                     suf = "eps", 
                     n = 20000,
                     nt = 30,
                     mne = 100, 
                     mnea = 1000,
                     ncp = 4,
                     max_it = 50,
                     nfi = 10,
                     umi = 4000,
                     ep = 0.5,
                     filenames = c("tree", "capi_tree_nag", "capi_tree"),
                     dir = "~/batch_online_rl/capi/pole_",
                     ylab = "Successful Episodes",
                     confidence.level = 0.95,
                     plot.error = TRUE, 
                     ylim = c(0,1), 
                     num.avg = 50, 
                     log = "")

{
  
  set.par()
  Y <- NULL
  Y.ui <- NULL
  Y.li <- NULL
  
  for (i in 1:length(filenames))
  {
    if (filenames[i] == "tree") fn <- ps(dir, paste(filenames[i],n, nt, mne, ncp, max_it, nfi, umi, ep, suf, sep="_"))
    else fn <- ps(dir, paste(filenames[i],n, nt, mne, mnea, ncp, max_it, nfi, umi, ep, suf, sep="_"))
       
    fn <- ps(fn, ".res")
    R <- as.matrix(read.table(fn))

    ci <- qnorm(1-(1-confidence.level)/2) * apply(R,2,sd) / sqrt(nrow(R))
    M <- apply(R, 2, mean)
    Y    <- cbind(Y, M)
    Y.ui <- cbind(Y.ui, M + ci)
    Y.li <- cbind(Y.li, M - ci)
  }
 
  mp(seq(n / nrow(Y), n, length=nrow(Y)), Y, Y.li = Y.li, Y.ui = Y.ui, lwd = 3, ylab = ylab, xlab = "Number of sample transitions",
     num.avg = num.avg, plot.error = plot.error, log = log, ylim = ylim, show.error = rep(FALSE,ncol(Y)))
  leg("bottomright", leg=filenames, lwd = 3)
}


show.pole.mnea <- function(
                     suf = "eps", 
                     n = 30000,
                     nt = 30,
                     mne = 20, 
                     mnea = c(20, 100, 500),
                     ncp = 4,
                     max_it = 50,
                     nfi = 10,
                     umi = 5000,
                     ep = 0.5,
                     dir = "~/batch_online_rl/capi/pole_",
                     ylab = "Successful Episodes",
                     confidence.level = 0.95,
                     plot.error = TRUE, 
                     show.shadow = TRUE,
                     show.error = FALSE,
                     ylim = NULL, 
                     num.avg = 50, 
                     log = "",
                     algs = c("capi_tree", "capi_tree_nag"),
                     leg.pos = "bottomright",
                     ...)

{
  
  set.par()
  Y <- NULL
  Y.ui <- NULL
  Y.li <- NULL
  
   # First the tree
   fn <- ps(dir, paste("tree",n, nt, mne, ncp, max_it, nfi, umi, ep, suf, sep="_"))
   fn <- ps(fn, ".res")
   R <- as.matrix(read.table(fn))

   ci <- apply(R,2,sd) / sqrt(nrow(R)) # qnorm(1-(1-confidence.level)/2) * apply(R,2,sd) / sqrt(nrow(R))
   M <- apply(R, 2, mean)
   Y    <- cbind(Y, M)
   Y.ui <- cbind(Y.ui, M + ci)
   Y.li <- cbind(Y.li, M - ci)

  # Now CAPI
  for (alg in algs)
    {
   for (i in 1:length(mnea))
    {
      fn <- ps(dir, paste(alg,n, nt, mne, mnea[i], ncp, max_it, nfi, umi, ep, suf, sep="_"))
      fn <- ps(fn, ".res")
      R <- as.matrix(read.table(fn))

#       ci <- qnorm(1-(1-confidence.level)/2) * apply(R,2,sd) / sqrt(nrow(R))
      ci <- apply(R,2,sd) / sqrt(nrow(R))
      M <- apply(R, 2, mean)
      Y    <- cbind(Y, M)
      Y.ui <- cbind(Y.ui, M + ci)
      Y.li <- cbind(Y.li, M - ci)
    }
 }
 
 leg="TREE-FQI"
 leg <- c(leg, paste("TREE-CAPI", mnea))
#  for (i in 1:length(algs)) leg <- c(leg, paste(algs[i], mnea))

  mp(seq(n / nrow(Y), n, length=nrow(Y)), Y, Y.li = Y.li, Y.ui = Y.ui, lwd = 3, ylab = ylab, xlab = "Number of sample transitions",
     num.avg = num.avg, plot.error = plot.error, log = log, ylim = ylim, 
     show.shadow=show.shadow, show.error = show.error, ...)
   # plot FQI again to make error bars show up
  Z <- matrix(Y[,1], nrow(Y), 1)
  Z.ui <- matrix(Y.ui[,1], nrow(Y), 1)
  Z.li <- matrix(Y.li[,1], nrow(Y), 1)
  
  mp(seq(n / nrow(Y), n, length=nrow(Y)), Z, Y.li = Z.li, Y.ui = Z.ui, lwd = 3, ylab = ylab, xlab = "Number of sample transitions",
     num.avg = num.avg, plot.error = plot.error, log = log, ylim = ylim, 
     show.shadow=show.shadow, show.error = show.error, add = TRUE, ...)

  
  leg(leg.pos, leg, lwd = 3)
}



plot.dpi <- function(
                     suf = "eps", 
                     nis = c(60, 190, 600, 1897, 6000, 18974, 60000),
                     hor = 50,
                     nar = 1,
                     max_it = 5,
                     nt = 30,
                     mne = c(20, 100, 500),
                     ncp = 4,
                     epsilon = 1,
                     dir = "~/batch_online_rl/capi/pole_",
                     ylab = "Successful Episodes",
                     confidence.level = 0.95,
                     plot.error = TRUE, 
                     show.shadow = TRUE,
                     show.error = FALSE,
                     ylim = NULL, 
                     num.avg = 50, 
                     log = "",
                     den = 1,
                     leg.pos = "topleft",
                     mne.capi = c(20, 500),
                     capi.leg = "TREE-CAPI 20 / 500",
                     ...
                    )
{
  set.par()
  
  Y <- matrix(0, length(nis), length(mne))
  Y.ui <- Y
  Y.li <- Y
  # DPI
  for (k in 1:length(mne))
  {
   for (i in 1:length(nis))
      {
         fn <- ps(dir, paste("dpi_tree" ,nis[i], hor, nar, max_it, nt, mne[k], ncp, epsilon, suf, sep="_"))
         fn <- ps(fn, ".res")
         R <- as.matrix(read.table(fn)) / den
         
#          ci <- qnorm(1-(1-confidence.level)/2) * apply(R,2,sd) / sqrt(nrow(R))
         ci <- apply(R,2,sd) / sqrt(nrow(R))
         M <- apply(R, 2, mean)
         Y[i,k]    <- M[length(M)]
         Y.ui[i,k] <- M[length(M)] + ci
         Y.li[i,k] <- M[length(M)] - ci
         
   #       y[i] <- apply(R,2,mean)[ncol(R)]
      }
  }


     # Now, CAPI
   Z <- matrix(0, 1, length(mne.capi))
   Z.ui <- Z
   Z.li <- Z
   
   for (i in 1:length(mne.capi))
   {
      fn <- ps(dir, paste("capi_tree", 30000, 30, 20, mne.capi[i], 4, 50, 10, 5000, 0.5, suf, sep="_"))
      fn <- ps(fn, ".res")
      R <- as.matrix(read.table(fn))
      ci <- apply(R,2,sd) / sqrt(nrow(R))
      M <- apply(R, 2, mean)
      Z[1,i]    <- M[length(M)]
      Z.ui[1,i] <- M[length(M)] + ci[length(ci)]
      Z.li[1,i] <- M[length(M)] - ci[length(ci)]
   }

  
   if (is.null(ylim))
   {
      ylim <- c(min(min(Y.li),min(Z.li)), max(max(Y.ui), max(Z.ui)))
   }
   x <- (nis * 2 * hor * nar * max_it) / 30000
#    print(x)
#    plot(x, y / 500, xlab = "Number of sample transitions (x 30,000)", lwd = 3, ylab = ylab, log = log, t="o", ylim=ylim)
   
   mp(x, Y, Y.li = Y.li, Y.ui = Y.ui, lwd = 3, ylab = ylab, xlab = "Number of sample transitions (x 30,000)",
     num.avg = num.avg, plot.error = plot.error, log = log, ylim = ylim, 
     show.shadow=show.shadow, show.error = show.error, xaxt = "n")
   
   leg <- paste("TREE-DPI", mne)
   leg(leg.pos, leg, lwd = 3)   
   
   mp(1, Z, Y.li = Z.li, Y.ui = Z.ui, lwd = 3, ylab = ylab, xlab = "Number of sample transitions (x 30,000)",
     num.avg = num.avg, plot.error = plot.error, log = log, ylim = ylim, 
     show.shadow=show.shadow, show.error = show.error, add = TRUE,
     col = get.col(6)[6], pch = get.pch(4)[4], xaxt = "n")
   x <- c(1, 10, 100, 1000)
   axis(1, at = x, labels = paste(round(x, d = 0)))
   
   x <- c(4, 4.5)
   if (length(capi.leg) == 1) x <- 6
   text(x, Z[,1:length(capi.leg)], capi.leg)
  
}
 

generate.all <- function()
{
  
   show.pole.mnea(mnea = c(20, 500), algs = c("capi_tree"), show.shadow = FALSE, show.error = TRUE)
   dev.copy2pdf(file = "~/tex/dev/reports/fig/capi_fqi_episodes.pdf")
   
   show.pole.mnea("stp", ylab="Number of steps", mnea = c(20, 500), algs = c("capi_tree"), show.shadow = FALSE, show.error = TRUE)
   dev.copy2pdf(file = "~/tex/dev/reports/fig/capi_fqi_steps.pdf")

   show.pole.mnea("ret", ylab="Return", mnea = c(20, 500), algs = c("capi_tree"), show.shadow = FALSE, show.error = TRUE)
   dev.copy2pdf(file = "~/tex/dev/reports/fig/capi_fqi_return.pdf")

   show.pole.mnea("tim", ylab="Seconds", mnea = c(20, 500), algs = c("capi_tree"), show.shadow = FALSE, show.error = TRUE, leg.pos = "topleft")
   dev.copy2pdf(file = "~/tex/dev/reports/fig/capi_fqi_time.pdf")

   plot.dpi("eps", ylab = "Successful episodes", log="x", show.shadow = FALSE, show.error = TRUE, leg.pos = "left")
   dev.copy2pdf(file = "~/tex/dev/reports/fig/capi_dpi_episodes.pdf")

   plot.dpi("stp", ylab = "Number of steps", log="x", show.shadow = FALSE, show.error = TRUE, leg.pos = "bottomright")
   dev.copy2pdf(file = "~/tex/dev/reports/fig/capi_dpi_steps.pdf")

   plot.dpi("ret", ylab = "Return", log="x", show.shadow = FALSE, show.error = TRUE, leg.pos = "bottomright")
   dev.copy2pdf(file = "~/tex/dev/reports/fig/capi_dpi_return.pdf")

   
   plot.dpi("tim", ylab = "Seconds", log="x", show.shadow = FALSE, show.error = TRUE, leg.pos = "topright", capi.le = c("TREE-CAPI 20", "TREE-CAPI 500"))
   dev.copy2pdf(file = "~/tex/dev/reports/fig/capi_dpi_time.pdf")

   
}


# ALL: (2, 10, 20, 50, 100, 200, 300, 500, 1000, 2000)
#  running 2, 20, 10, 50, 50, 200, 1000


show.hiv <- function(
		     col = 1,
                     nt = 30,
                     mne = 50, 
		     mne.capi = NULL,
                     mnea = c(2, 10, 20, 50, 100, 200, 300, 500, 1000, 2000),
                     ncp = 8,
                     max_it = 50,
                     nfi = 10,
                     dir = "~/hiv/files/hiv_",
                     ylab = "Return",
                     confidence.level = 0.95,
                     plot.error = TRUE, 
                     num.avg = 50,
		     log = "x",
		     ylim = NULL
                    )
{
    if (is.null(mne.capi)) mne.capi <- mne

    set.par()
    Y <- matrix(0, length(mnea), 2) 
    Y.ui <- Y
    Y.li <- Y
    
    # first the tree
    fn <- ps(dir, ps(paste("tree", nt, mne,  ncp, max_it, nfi, sep="_"),".res"))
    T <- read.table(fn)[,col]
    Y[,1] <- mean(T)
    ci <- sd(T) / sqrt(length(T))
    Y.ui[,1] <- Y[,1] + ci 
    Y.li[,1] <- Y[,1] - ci

    
    for (i in 1:length(mnea))
    {
	
      fn <- ps(dir, ps(paste("capi_tree", nt, mne.capi, mnea[i], ncp, max_it, nfi, sep="_"),".res"))
      T <- read.table(fn)[,col]
      ci <- sd(T) / sqrt(length(T))
      
      print(paste(fn, mean(T)))
      
      Y[i,2] <- mean(T)
      Y.ui[i,2] <- Y[i,2] + ci 
      Y.li[i,2] <- Y[i,2] - ci

      print((Y[i,2] - Y[i,1]) / Y[i,1])
    }
  
    mp(mnea, Y, Y.li = Y.li, Y.ui = Y.ui, lwd = 3,  ylab = ylab, xlab = expression(eta[pi]*" (log)"),  num.avg = 
           num.avg, plot.error = plot.error, log = log, ylim = ylim,  show.error = c(FALSE,TRUE), show.shadow = c(TRUE,FALSE),
           pch = c(NA, get.pch(1)))
    leg("bottomleft", leg=c("Tree-FQI", "Tree-CAPI"), lwd = 3, pch = c(NA,get.pch(1)) )
}


print("capi.data.analysis.R loaded")