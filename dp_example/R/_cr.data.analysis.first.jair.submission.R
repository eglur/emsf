source("util.R")


set.par <- function() 
{
# set parameters used  in the plots
   par(
      lwd = 1.5,
      mar = c(5, 6, 2, 2) -1 + 0.1,
      ps  = 20)
}

get.lty <- function(n) 
{
   types <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
   types[(1:n - 1) %% length(types) + 1]
}


get.col <- function(n, alpha = 1) 
{
   col <- c(rgb(0,0,0, alpha), rgb(1,0,0, alpha), rgb(0,0,1, alpha), rgb(0,1,0,alpha), rgb(0,1,1,alpha), 
      rgb(1,0,1,alpha), rgb(1,1,0,alpha), rgb(0.5,0.5,0.5,alpha))
   col[(1:n - 1) %% length(col) + 1]
}

get.pch <- function(n) 
{
   if (n <= 5) seq(21, 21 + n -1, by = 1)
   else c(21:25, 1:(n-6)) 
}

mp <- function(x, Y, ...)
{
   p <- ncol(Y)
   lty= get.lty(p)
   matplot(x, Y, t="o", bg = get.col(p), col=get.col(p), lty= get.lty(p), pch= get.pch(p), ...)
}

leg <- function(pos = "topleft", leg, y.intersp = 1.5, ...)
{
   p <- length(leg)
   legend(pos, leg, pt.bg = get.col(p), col=get.col(p), lty= get.lty(p), pch= get.pch(p), y.intersp  = y.intersp ,...)
}


show.sizes <- function(nc = 2:6, lt = 10, res = c(2, 10, 30),
   dir = "~/cr/experiments/", red = TRUE, ...)
{
   
   if (red) suffix <- "_red"
   else suffix <- ""

   p <- length(nc)
   R <- matrix(0, p,  2 + length(res))
   
   for (i in 1:p)
   {
      R[i,1] <- read.table(ps(dir, "pi_sizes_", nc[i], "_", lt, ".txt"))[1,1]
      R[i,2] <- read.table(ps(dir, "pi_red_sizes_", nc[i], "_", lt, ".txt"))[1,1]
      
      for (j in 1:length(res))
      {
         R[i, 2 + j] <- read.table(ps(dir, "pisf", suffix, "_sizes_", nc[i], "_", lt, "_", res[j], ".txt"))[1,2]
      }
   }
   
   set.par()
   mp(nc, R, ylab="|S|", xlab="Number of components", ...)
   leg("topleft", c("PI", "PI-RED", ps("PISF-", res)))

   # Since they were not actually run...
   p <- length(nc)
   points(c(nc[p], nc[p]), R[p, 1:2], pch=get.pch(2), bg="WHITE", col=get.col(2))
   
}


<<<<<<< cr.data.analysis.R
show.times <- function(nc = 2:6, lt = 10, res = c(2, 10, 30),
   dir = "~/cr/experiments/", red = TRUE, nf = FALSE, add.build = FALSE, ...)
=======
show.times <- function(nc = 2:6, lt = 10, res = c(100,1000),
   dir = "~/cr/experiments/", red = TRUE, add.build = FALSE, ...)
>>>>>>> 1.6
{
   
   if (red) suffix <- "_red"
   else suffix <- ""

   if (nf) nfs <- "_nf"
   else nfs <- ""

   p <- length(nc)
   R <- matrix(0, p, 2 + length(res))
   
   for (i in 1:p)
   {
      T <- read.table(ps(dir, "pi",nfs,"_time_", nc[i], "_", lt, ".txt"))
      if (add.build) R[i,1] <- T[1,1] + T[2,1]
      else R[i,1] <- T[2,1]
      
      T <- read.table(ps(dir, "pi_red",nfs,"_time_", nc[i], "_", lt, ".txt"))
      if (add.build) R[i,2] <- T[1,1] + T[2,1]
      else R[i,2] <- T[2,1]
      
       for (j in 1:length(res))
       {
          T <- read.table(ps(dir, "pisf", suffix, nfs, "_time_", nc[i], "_", lt, "_", res[j], ".txt"))
          if (add.build) R[i, 2 + j] <- T[1,1] + T[2,1]
          else R[i, 2 + j] <- T[2,1]
       }
   }
   
   set.par()
   mp(nc, R, ylab="Seconds", xlab="Number of components", ...)

   leg("topleft", c("PI", "PI-RED", ps("PISF-", res)))
   
    # Since they were not actually run...
#    p <- length(nc)
#    points(c(nc[p], nc[p]), R[p, 1:2], pch=get.pch(2), bg="WHITE", col=get.col(2))

}


<<<<<<< cr.data.analysis.R
show.error.v <- function(nc = 2:5, lt = 10, res = c(2, 10, 30),
   dir = "~/cr/experiments/", red = TRUE, nf = FALSE, load.last = FALSE, ...)
=======
show.error.v <- function(nc = 2:5, lt = 10, res = c(100, 1000),
   dir = "~/cr/experiments/", red = TRUE, ...)
>>>>>>> 1.6
{
   
   if (red) suffix <- "_red"
   else suffix <- ""

<<<<<<< cr.data.analysis.R
   if (nf) snf <- "_nf"
   else snf <- ""

   R <- NULL
   if (load.last) R <- read.table("~/tmp/last_error_v.txt")
   else
=======
   p <- length(nc)
   R <- matrix(0, p, lt  + length(res))
   
   for (i in 1:p)
>>>>>>> 1.6
   {
      R[i,1] <- mean(read.table(ps(dir, "pi_red_error_", nc[i], "_", lt, ".txt")))
    
      for (ma in 1:(lt-1))
      {
         vt <- read.table(ps(dir, "std_error_", nc[i], "_", lt, "_", ma, ".txt"))
         R[i,1 + ma] <- mean(vt)
      }
    
      
#       R[i,2] <- -Inf
#       best <- 1
#       for (ma in 1:(lt-1))
#       {
#          vt <- read.table(ps(dir, "std_error_", nc[i], "_", lt, "_", ma, ".txt"))
#          e <- mean(vt)
#          if (e > R[i,2]) 
#          {
#             R[i,2] <- e
#             best <- ma
#          }
#       }
      
#       print(paste("best", best))
      
      for (j in 1:length(res))
      {
<<<<<<< cr.data.analysis.R
	 print(ps(dir, "pi", snf, "_v_", nc[i], "_", lt, ".txt"))
         vo <- read.table(ps(dir, "pi", snf, "_v_", nc[i], "_", lt, ".txt"))

         
         R[i,1] <- -Inf
         for (ma in 0:(lt-1))
         {
	    print(ps(dir, "std", snf, "_v_", nc[i], "_", lt, "_", ma,".txt"))
            vt <- read.table(ps(dir, "std", snf, "_v_", nc[i], "_", lt, "_", ma,".txt"))
            e <- sum((vt - vo) / abs(vo)  * 100) / nrow(vo)
            if (e > R[i,1]) R[i,1] <- e
         }
      
         for (j in 1:length(res))
         {
	    print(ps(dir, "pisf", suffix, snf, "_v_", nc[i], "_", lt, "_", res[j], ".txt"))
            vt <- read.table(ps(dir, "pisf", suffix, snf, "_v_", nc[i], "_", lt, "_", res[j], ".txt"))
    
            R[i,1 + j] <- sum((vt - vo) / abs(vo)  * 100) / nrow(vo)
         }
=======
         vt <- read.table(ps(dir, "pisf", suffix, "_error_", nc[i], "_", lt, "_", res[j], ".txt"))
         R[i,lt + j] <- mean(vt)
>>>>>>> 1.6
      }
  }
     
   set.par()
<<<<<<< cr.data.analysis.R
   mp(nc, R, ylab="Expected loss (%)", xlab="Number of components", ...)
   leg("left", c("STD", ps("PISF-", res)))
   wt(R, "~/tmp/last_error_v.txt")
=======
   mp(nc, R, ylab="Expected gain", xlab="Number of components", ...)
   leg("bottomleft", c("PI", paste("STD",1:(lt-1)), ps("PISF-", res)))
#    wt(R, "~/tmp/last_error_v.txt")
>>>>>>> 1.6
}


show.diff.v <- function(nc = 6, lt = 10, res = c(2, 10, 30),
   dir = "~/cr/experiments/", red = TRUE, nf = FALSE, load.last = FALSE, ...)
{
   
   if (red) suffix <- "_red"
   else suffix <- ""

   if (nf) snf <- "_nf"
   else snf <- ""

   R <- NULL
   if (load.last) R <- read.table("~/tmp/last_diff_v.txt")
   else
   {
      R <- array(0, length(res))
      
      std.avg <- -Inf
      vs <- NULL
      for (ma in 0:(lt-1))
      {
         vt <- read.table(ps(dir, "std", snf, "_v_", nc, "_", lt, "_", ma, ".txt"))
         avg.ret <- mean(vt)
         if (avg.ret > std.avg) 
         { 
            std.avg <- avg.ret
            vs <- vt
         }
      }
      
      for (j in 1:length(res))
      {
         vt <- read.table(ps(dir, "pisf", suffix, snf, "_v_", nc, "_", lt, "_", res[j], ".txt"))
         R[j] <- sum((vt - vs) / abs(vs)  * 100) / nrow(vs)
      }
   }
   
   set.par()

   barplot(R, ylab="Expected gain (%)", names = ps("PISF-", res) , ...)
   wt(R, "~/tmp/last_diff_v.txt")
}


generate.tables <- function(nc = c(2,3,4,5,6), lt = 10, res = c(100,1000),
                            dir = "~/cr/experiments/")
{
for (suf  in c("sizes", "time","error"))
{  

   for (i in 1:length(nc))
   {
      cat(ps(i+1, " & "))
      if (i < length(nc))
      {
         T <- read.table(ps(dir,"pi_", suf, "_", nc[i], "_",lt,".txt"))
         d <- 1
         if (suf == "time") T[,1] <- apply(T,1,sum)
         else if (suf == "error") T <- T * 1000
         else if (suf == "sizes") d <- 0
            
         cat(paste("$ ", round(mean(T[,1]), d=d), "$  $\\ci{", round(qnorm(1-(1-0.95)/2)*sd(T[,1])/sqrt(nrow(T)), d=d),"} $ & "))
      }
      else cat("-- & ")
      
      if (i < length(nc))
      {
         T <- read.table(ps(dir,"pi_red_", suf, "_", nc[i], "_",lt,".txt")) 
         d <- 1
         if (suf == "time") T[,1] <- apply(T,1,sum)
         else if (suf == "error") T <- T * 1000
         else if (suf == "sizes") d <- 0
         cat(paste("$ ", round(mean(T[,1]), d=d), "$  $\\ci{", round(qnorm(1-(1-0.95)/2)*sd(T[,1])/sqrt(nrow(T)), d=d),"} $ & "))
      }
      else cat("-- & ")

      if (suf == "sizes") col <- 2
      else col <- 1
         
      for (r in 1:length(res))
      {   
         T <- read.table(ps(dir,"pisf_red_",suf,"_", nc[i], "_",lt,"_", res[r], ".txt")) 
         d <- 1
         if (suf == "time") T[,col] <- apply(T,1,sum)
         else if (suf == "error") T <- T * 1000
         else if (suf == "sizes") d <- 0
         cat(paste("$", round(mean(T[,col]), d=d), "$  $\\ci{", round(qnorm(1-(1-0.95)/2)*sd(T[,col])/sqrt(nrow(T)), d=d)), "} $")
         if (r < length(res)) cat(" & ")
      }

      cat("\\\\ \n")
   }
      cat("\n")
   
}

}



generate.tables.transpose <- function(nc = c(2,3,4,5,6), lt = 10, res = c(100,1000),
                            dir = "~/cr/experiments/")
{

# SIZES
for (suf  in c("sizes", "time"))
{  
   print(suf)
   print(" ")
   # PI
   cat("PI & ")
   for (i in 1:length(nc))
      {
         if (i < length(nc))
         {
            T <- read.table(ps(dir,"pi_", suf, "_", nc[i], "_",lt,".txt"))
            if (suf == "time") T[,1] <- apply(T,1,sum)
            cat(paste("$ ", round(mean(T[,1]), d=1), "\\pm", round(qnorm(1-(1-0.95)/2)*sd(T[,1])/sqrt(nrow(T)), d=1)," $ & "))
         }
         else cat("-- \\\\ \n")
      }

   # PI-RED
   cat("PI-RED & ")
   for (i in 1:length(nc))
      {
         if (i < length(nc))
         {
            T <- read.table(ps(dir,"pi_red_", suf, "_", nc[i], "_",lt,".txt")) 
            if (suf == "time") T[,1] <- apply(T,1,sum)
            cat(paste("$ ", round(mean(T[,1]), d=1), "\\pm", round(qnorm(1-(1-0.95)/2)*sd(T[,1])/sqrt(nrow(T)), d=1)," $ & "))
         }
         else cat("-- \\\\ \n")
      }
      
   # PISF
   col <- 2
   for (r in res)
   {   
      cat(ps("PISF-",r," & "))

      for (i in 1:length(nc))
         {
               T <- read.table(ps(dir,"pisf_red_",suf,"_", nc[i], "_",lt,"_", r, ".txt")) 
               if (suf == "time") T[,col] <- apply(T,1,sum)
               cat(paste("$", round(mean(T[,col]), d=1), "\\pm", round(qnorm(1-(1-0.95)/2)*sd(T[,col])/sqrt(nrow(T)), d=1)), "$")
               if (i < length(nc)) cat(" & ")
         }
         cat("\\\\ \n")
   }
   }

}


print("cr.data.analysis.R")
