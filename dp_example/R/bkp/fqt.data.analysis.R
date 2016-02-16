source("data.plot.R")

show.all <- function(show.times = TRUE,  ...)
{
   show.comp(TRUE, "ret", "puddle", ...)
   dev.copy2pdf(file = "~/tmp/puddle_nt.pdf")
   
   x11()
   show.comp(FALSE, "ret", "puddle", ...)
   dev.copy2pdf(file = "~/tmp/puddle_n.pdf")

   x11()
   show.comp(TRUE, "ret", "mountain", ...)
   dev.copy2pdf(file = "~/tmp/mountain_nt.pdf")

   x11()
   show.comp(FALSE, "ret", "mountain", ...)
   dev.copy2pdf(file = "~/tmp/mountain_n.pdf")

   if (show.times)
   {
      x11()
      show.comp(TRUE, "tim", "puddle", show.leg = TRUE, ...)
      dev.copy2pdf(file = "~/tmp/puddle_nt_time.pdf")

      x11()
      show.comp(FALSE, "tim", "puddle", ...)
      dev.copy2pdf(file = "~/tmp/puddle_n_time.pdf")

      x11()
      show.comp(TRUE, "tim", "mountain", show.leg = TRUE, ...)
      dev.copy2pdf(file = "~/tmp/mountain_nt_time.pdf")

      
      x11()
      show.comp(FALSE, "tim", "mountain", ...)
      dev.copy2pdf(file = "~/tmp/mountain_n_time.pdf")

      }
      
}

show.comp <- function(
             vary.nt = TRUE,
             crit = "ret",
             env = "puddle",
             n = c(2000, 4000, 6000, 8000, 10000, 12000),
             nt = c(1, 5, 10, 15, 20, 25, 30),
             mne = 30,
             ncp = 10,
             nfi = 150, 
             it = 200,
             epsilon = 1,
             dir = "~/fontes/c++/batch_online_rl/fqt/",
             ylab = "Return",
             log = "",
             plot = c(TRUE, TRUE, TRUE, FALSE),
             show.leg = FALSE
             )


{
   set.par()
   
   if (vary.nt) inds <- 1:length(nt)
   else inds <- 1:length(n)

   Y <- matrix(0, length(inds), sum(plot))
   Y.li <- Y
   Y.ui <- Y
   
   
   nt.value <- 30
   n.value  <- 10000
   for (i in inds)
   {
      if (vary.nt) nt.value <- nt[i]
      else n.value <- n[i]
      
      col <- 1
      
      if (plot[1])
      {
         # TREE-RL NON-FIXED
         fn <- paste(dir, paste(env, "tree", n.value, nt.value, mne, ncp, it, nfi, n.value, epsilon, crit, sep="_"), ".res", sep="")
         D <- as.matrix(read.table(fn))
         if (nrow(D) < 30) print(fn)
         y <- D[,ncol(D)]
         Y[i,col] <- mean(y)
         ci <- sd(y) / sqrt(length(y))
         Y.li[i,col] <- Y[i,col] - ci
         Y.ui[i,col] <- Y[i,col] + ci
       
         col <- col + 1
      }
      
      if (plot[2])
      {
         # TREE-RL FIXED
         fn <- paste(dir, paste(env, "tree", n.value, nt.value, mne, ncp, it, 1, n.value, epsilon, crit, sep="_"), ".res", sep="")
         D <- as.matrix(read.table(fn))
         if (nrow(D) < 30) print(fn)
         y <- D[,ncol(D)]
         Y[i,col] <- mean(y)
         ci <- sd(y) / sqrt(length(y))
         Y.li[i,col] <- Y[i,col] - ci
         Y.ui[i,col] <- Y[i,col] + ci 
       
         col <- col + 1
      }
      
      if (plot[3])
      {
         # ON-LINE TREE-RL FIXED (FQT)
         fn <- paste(dir, paste(env, "fqt", n.value, nt.value, mne, ncp, it, n.value, epsilon, crit, sep="_"), ".res", sep="")
         D <- as.matrix(read.table(fn))
         if (nrow(D) < 30) print(fn)
         y <- D[,ncol(D)]
         Y[i,col] <- mean(y)
         ci <- sd(y) / sqrt(length(y))
         Y.li[i,col] <- Y[i,col] - ci
         Y.ui[i,col] <- Y[i,col] + ci 
       
         col <- col + 1
      }
      

      if (plot[4])
      {
            # ON-LINE TREE-RL ADAPTIVE (FQT-ADP)
            fn <- paste(dir, paste(env, "fqt_adp", n.value, nt.value, mne, ncp, it, n.value, epsilon, crit, sep="_"), ".res", sep="")
            D <- as.matrix(read.table(fn))
            if (nrow(D) < 30) print(fn)
            y <- D[,ncol(D)]
            Y[i,col] <- mean(y)
            ci <- sd(y) / sqrt(length(y))
            Y.li[i,col] <- Y[i,col] - ci
            Y.ui[i,col] <- Y[i,col] + ci 
       
         col <- col + 1
      }
      
   }
 
 
   if (crit == "ret") ylab <- "Return"
   else if (crit == "tim") ylab <- "Seconds"
   
   x <- n
   if (vary.nt) x <- nt
      
     
   xlab <- "n"
   if (vary.nt) xlab <- "|T|"
   
   mp(x, Y, Y.li = Y.li, Y.ui = Y.ui, lwd = 3,  ylab = ylab, xlab = xlab, log = log)
   txt <- c("FQTI", "FQTI-F", "TBSF", "TBSF-ADP")[plot]
   if (show.leg) leg("topleft", txt)
}




show.inc <- function(
             crit = "ret",
             env = "mountain",
             n = 10000,
             umi = 1000,
             nt = 30,
             mne = 30,
             ncp = 10,
             nfi = 150, 
             it = 200,
             epsilon = 1,
             dir = "~/fontes/c++/batch_online_rl/fqt/",
             ylab = "Return",
             log = "",
             plot = c(TRUE, TRUE, TRUE)
             )


{
   set.par()
   
   Y <- NULL
   Y.li <- Y
   Y.ui <- Y
   
      col <- 1
      
      if (plot[1])
      {
         # TREE-RL NON-FIXED
         fn <- paste(dir, paste(env, "tree", n, nt, mne, ncp, it, nfi, umi, epsilon, crit, sep="_"), ".res", sep="")
         D <- as.matrix(read.table(fn))
         if (nrow(D) < 30) print(fn)
         y <- apply(D,2,mean)
         Y <- cbind(Y, y)
         ci <- apply(D,2,sd) / sqrt(nrow(D))
         Y.li <- cbind(Y.li, y - ci)
         Y.ui <- cbind(Y.ui, y + ci)
       
         col <- col + 1
      }
      
     
      if (plot[2])
      {
         # ON-LINE TREE-RL FIXED (FQT)
         fn <- paste(dir, paste(env, "fqt", n, nt, mne, ncp, it, umi, epsilon, crit, sep="_"), ".res", sep="")
         D <- as.matrix(read.table(fn))
         if (nrow(D) < 30) print(fn)
         y <- apply(D,2,mean)
         Y <- cbind(Y, y)
         ci <- apply(D,2,sd) / sqrt(nrow(D))
         Y.li <- cbind(Y.li, y - ci)
         Y.ui <- cbind(Y.ui, y + ci)
       
         col <- col + 1
      }
      

      if (plot[3])
      {
         # ON-LINE TREE-RL ADAPTIVE (FQT-ADP)
         fn <- paste(dir, paste(env, "fqt_adp", n, nt, mne, ncp, it, umi, epsilon, crit, sep="_"), ".res", sep="")
         D <- as.matrix(read.table(fn))
         if (nrow(D) < 30) print(fn)
         y <- apply(D,2,mean)
         Y <- cbind(Y, y)
         ci <- apply(D,2,sd) / sqrt(nrow(D))
         Y.li <- cbind(Y.li, y - ci)
         Y.ui <- cbind(Y.ui, y + ci)
       
         col <- col + 1
      }
      
 
 
   if (crit == "ret") ylab <- "Return"
   else if (crit == "tim") ylab <- "Seconds"
   
   x <- seq(umi, n, by = umi)
 
     
   xlab <- "n"
   
   mp(x, Y, Y.li = Y.li, Y.ui = Y.ui, lwd = 3,  ylab = ylab, xlab = xlab, log = log)
}




show.hiv <- function(
             nt = 30,
             mne = c(200,100, 50, 10),
             ncp = 8,
             it = 50,
             env = "hiv",
             dir = "~/hiv/fqt/",
             ylab = "Return",
             log = "",
             plot = c(TRUE, TRUE, TRUE,FALSE, TRUE, TRUE),
             crit = 1,
             leg.position = "right"
             )
{
   set.par()
   
   Y <- matrix(0, length(mne), sum(plot))
   Y.li <- Y
   Y.ui <- Y
   
     
      for (i in 1:length(mne))
      {
         col <- 1
         
         if (plot[1])
         {
            fn <- paste(dir, "hiv_tree_30_200_8_200_150.res", sep="")
            D <- as.matrix(read.table(fn))
            if (nrow(D) < 30) print(fn)
            Y[i,col] <- mean(D[,crit])
            ci <- sd(D[,crit]) / sqrt(nrow(D))
            Y.li[i,col] <- Y[i,col] - ci
            Y.ui[i,col] <- Y[i,col] + ci
            col <- col + 1
            
         }
         
         if (plot[2])
         {
            fn <- paste(dir, "hiv_tree_30_100_8_200_150.res", sep="")
            D <- as.matrix(read.table(fn))
            if (nrow(D) < 30) print(fn)
            Y[i,col] <- mean(D[,crit])
            ci <- sd(D[,crit]) / sqrt(nrow(D))
            Y.li[i,col] <- Y[i,col] - ci
            Y.ui[i,col] <- Y[i,col] + ci
            col <- col + 1
            
         }      

         if (plot[3])
         {
            # FQT FIXED
            fn <- paste(dir, paste(env, "fqt", nt, mne[i], ncp, it, sep="_"), ".res", sep="")
            D <- as.matrix(read.table(fn))
            if (nrow(D) < 30) print(fn)
            Y[i,col] <- mean(D[,crit])
            ci <- sd(D[,crit]) / sqrt(nrow(D))
            Y.li[i,col] <- Y[i,col] - ci
            Y.ui[i,col] <- Y[i,col] + ci
            col <- col + 1
         }
         
         if (plot[4])
         {
            # FQT ADP
            fn <- paste(dir, paste(env, "fqt_adp", nt, mne[i], ncp, it, sep="_"), ".res", sep="")
            D <- as.matrix(read.table(fn))
            if (nrow(D) < 30) print(fn)
            Y[i,col] <- mean(D[,crit])
            ci <- sd(D[,crit]) / sqrt(nrow(D))
            Y.li[i,col] <- Y[i,col] - ci
            Y.ui[i,col] <- Y[i,col] + ci
            col <- col + 1
         }

         if (plot[5])
         {
            # FQT FIXED with 300 patients
            fn <- paste(dir, "p300/", paste(env, "fqt", nt, mne[i], ncp, it, sep="_"), ".res", sep="")
            D <- as.matrix(read.table(fn))
            if (nrow(D) < 30) print(fn)
            Y[i,col] <- mean(D[,crit])
            ci <- sd(D[,crit]) / sqrt(nrow(D))
            Y.li[i,col] <- Y[i,col] - ci
            Y.ui[i,col] <- Y[i,col] + ci
            col <- col + 1
         }

         if (plot[6])
         {
            # FQT FIXED with 600 patients
            fn <- paste(dir, "p600/", paste(env, "fqt", nt, mne[i], ncp, it, sep="_"), ".res", sep="")
            D <- as.matrix(read.table(fn))
            if (nrow(D) < 30) print(fn)
            Y[i,col] <- mean(D[,crit])
            ci <- sd(D[,crit]) / sqrt(nrow(D))
            Y.li[i,col] <- Y[i,col] - ci
            Y.ui[i,col] <- Y[i,col] + ci
            col <- col + 1
         }         
         
         
      }
      
     
   xlab <- "mne"
   
   mp(mne, Y, Y.li = Y.li, Y.ui = Y.ui, lwd = 3,  ylab = ylab, xlab = xlab, log = log)
   leg(leg.position, c("TREE-200", "TREE-100", "FQT", "FQT-ADP", "FQT 300")[plot])
}



show.hiv.box <- function(
             nt = 30,
             mne = c(200,100, 50, 10),
             ncp = 8,
             it = 50,
             env = "hiv",
             dir = "~/hiv/fqt/",
             ylab = "Return",
             log = "",
             plot = c(TRUE, TRUE, TRUE,FALSE, TRUE, TRUE, FALSE),
             crit = 1,
             leg.position = "right"
             )
{
   set.par()
   
   if (crit == 2) plot[7] <- FALSE
      
   Y <- matrix(0, 50, length(mne) * sum(plot[3:6]) + sum(plot[c(1,2,7)]))
           
         col <- 1
        
        if (plot[1])
         {
               # FQT FIXED
               fn <- paste(dir, "hiv_tree_30_200_8_200_150.res", sep="")
               print(fn)
                D <- as.matrix(read.table(fn))
               if (nrow(D) < 30) print(fn)
               Y[,col] <- D[,crit]
               col <- col + 1
         }
         
        if (plot[2])
         {
               # FQT FIXED
               fn <- paste(dir, "hiv_tree_30_100_8_200_150.res", sep="")
               print(fn)

               D <- as.matrix(read.table(fn))
               if (nrow(D) < 30) print(fn)
               Y[,col] <- D[,crit]
               col <- col + 1
         }
         

         if (plot[3])
         {
            for (i in 1:length(mne))
            {
               # FQT FIXED
               fn <- paste(dir, paste(env, "fqt", nt, mne[i], ncp, it, sep="_"), ".res", sep="")
               print(fn)

               D <- as.matrix(read.table(fn))
               if (nrow(D) < 30) print(fn)
               Y[,col] <- D[,crit]
               col <- col + 1
              }
         }
         
         if (plot[4])
         {
            for (i in 1:length(mne))
            {
               # FQT ADP
               fn <- paste(dir, paste(env, "fqt_adp", nt, mne[i], ncp, it, sep="_"), ".res", sep="")
               print(fn)

               D <- as.matrix(read.table(fn))
               if (nrow(D) < 30) print(fn)
               Y[,col] <- D[,crit]
               col <- col + 1
            }
         }

         if (plot[5])
         {
           for (i in 1:length(mne))
            {
               # FQT FIXED with 300 patients
               fn <- paste(dir, "p300/", paste(env, "fqt", nt, mne[i], ncp, it, sep="_"), ".res", sep="")
               print(fn)

               D <- as.matrix(read.table(fn))
               if (nrow(D) < 30) print(fn)
               Y[,col] <- D[,crit]
               col <- col + 1
            }
         }
         
         
         if (plot[6])
         {
           for (i in 1:length(mne))
            {
               # FQT FIXED with 600 patients
               fn <- paste(dir, "p600/", paste(env, "fqt", nt, mne[i], ncp, it, sep="_"), ".res", sep="")
               print(fn)

               D <- as.matrix(read.table(fn))
               if (nrow(D) < 30) print(fn)
               Y[,col] <- D[,crit]
               col <- col + 1
            }
         }
         
         if (plot[7])
         {
               # FQT FIXED with 300 patients
               fn <-  paste(dir, "hiv_random.res", sep="")
               print(fn)

               D <- as.matrix(read.table(fn))
               if (nrow(D) < 30) print(fn)
               Y[,col] <- D[,crit]
               col <- col + 1
         }
      
      
      boxplot(Y, log = log, range = 0)
}


show.hiv.pareto <- function(
             nt = 30,
             mne = c(200, 100, 50, 10),
             ncp = 8,
             it = 50,
             env = "hiv",
             dir = "~/hiv/fqt/",
             ylab = "Return",
             plot = c(TRUE, TRUE, TRUE,FALSE, TRUE, TRUE, FALSE),
             leg.position = "right", 
             deltax = 2.8e3,
             deltay = 1e6,
             log = "y",
             scale.ci = 1,
             ...
             )
{
      par(
      lwd = 2,
      mar = c(5, 6, 3, 2) -1 + 0.1,
      ps  = 16)
   
      
   num.algs <- length(mne) * sum(plot[3:6]) + sum(plot[c(1,2,7)])

   x <- array(0, num.algs)
   y <- array(0, num.algs)
   ubx <- array(0, num.algs)
   lbx <- array(0, num.algs)
   uby <- array(0, num.algs)
   lby <- array(0, num.algs)
      
        alg <- 1
        
        if (plot[1])
         {
               # FQT FIXED
               fn <- paste(dir, "hiv_tree_30_200_8_200_150.res", sep="")
               D <- as.matrix(read.table(fn))
               if (nrow(D) < 30) print(fn)
               
               x[alg] <- mean(D[,2])
               y[alg] <- mean(D[,1])
      
               ci <- scale.ci * sd(D[,2]) / sqrt(nrow(D))
               ubx[alg] <- x[alg] + ci
               lbx[alg] <- x[alg] - ci
               
               ci <- scale.ci * sd(D[,1]) / sqrt(nrow(D))
               uby[alg] <- y[alg] + ci
               lby[alg] <- y[alg] - ci

               alg <- alg + 1
         }
         
        if (plot[2])
         {
               # FQT FIXED
               fn <- paste(dir, "hiv_tree_30_100_8_200_150.res", sep="")

               D <- as.matrix(read.table(fn))
               if (nrow(D) < 30) print(fn)
               
               x[alg] <- mean(D[,2])
               y[alg] <- mean(D[,1])
      
               ci <- scale.ci * sd(D[,2]) / sqrt(nrow(D))
               ubx[alg] <- x[alg] + ci
               lbx[alg] <- x[alg] - ci
               
               ci <- scale.ci * sd(D[,1]) / sqrt(nrow(D))
               uby[alg] <- y[alg] + ci
               lby[alg] <- y[alg] - ci

               alg <- alg + 1         
         }
         

         if (plot[3])
         {
            for (i in 1:length(mne))
            {
               # FQT FIXED
               fn <- paste(dir, paste(env, "fqt", nt, mne[i], ncp, it, sep="_"), ".res", sep="")

               D <- as.matrix(read.table(fn))
               if (nrow(D) < 30) print(fn)
               
               x[alg] <- mean(D[,2])
               y[alg] <- mean(D[,1])
      
               ci <- scale.ci * sd(D[,2]) / sqrt(nrow(D))
               ubx[alg] <- x[alg] + ci
               lbx[alg] <- x[alg] - ci
               
               ci <- scale.ci * sd(D[,1]) / sqrt(nrow(D))
               uby[alg] <- y[alg] + ci
               lby[alg] <- y[alg] - ci

               alg <- alg + 1
            }
         }
         
         if (plot[4])
         {
            for (i in 1:length(mne))
            {
               # FQT ADP
               fn <- paste(dir, paste(env, "fqt_adp", nt, mne[i], ncp, it, sep="_"), ".res", sep="")

               D <- as.matrix(read.table(fn))
               if (nrow(D) < 30) print(fn)
               
               x[alg] <- mean(D[,2])
               y[alg] <- mean(D[,1])
      
               ci <- scale.ci * sd(D[,2]) / sqrt(nrow(D))
               ubx[alg] <- x[alg] + ci
               lbx[alg] <- x[alg] - ci
               
               ci <- scale.ci * sd(D[,1]) / sqrt(nrow(D))
               uby[alg] <- y[alg] + ci
               lby[alg] <- y[alg] - ci

               alg <- alg + 1
            }
         }

         if (plot[5])
         {
           for (i in 1:length(mne))
            {
               # FQT FIXED with 300 patients
               fn <- paste(dir, "p300/", paste(env, "fqt", nt, mne[i], ncp, it, sep="_"), ".res", sep="")

               D <- as.matrix(read.table(fn))
               if (nrow(D) < 30) print(fn)
               
               x[alg] <- mean(D[,2])
               y[alg] <- mean(D[,1])
      
               ci <- scale.ci * sd(D[,2]) / sqrt(nrow(D))
               ubx[alg] <- x[alg] + ci
               lbx[alg] <- x[alg] - ci
               
               ci <- scale.ci * sd(D[,1]) / sqrt(nrow(D))
               uby[alg] <- y[alg] + ci
               lby[alg] <- y[alg] - ci

               alg <- alg + 1
           }
         }
         
         
         if (plot[6])
         {
           for (i in 1:length(mne))
            {
               # FQT FIXED with 600 patients
               fn <- paste(dir, "p600/", paste(env, "fqt", nt, mne[i], ncp, it, sep="_"), ".res", sep="")

               D <- as.matrix(read.table(fn))
               if (nrow(D) < 30) print(fn)
               
               x[alg] <- mean(D[,2])
               y[alg] <- mean(D[,1])
      
               ci <- scale.ci * sd(D[,2]) / sqrt(nrow(D))
               ubx[alg] <- x[alg] + ci
               lbx[alg] <- x[alg] - ci
               
               ci <- scale.ci * sd(D[,1]) / sqrt(nrow(D))
               uby[alg] <- y[alg] + ci
               lby[alg] <- y[alg] - ci

               alg <- alg + 1            }
         }
         
         if (plot[7])
         {
               # random
               fn <-  paste(dir, "hiv_random.res", sep="")

               D <- as.matrix(read.table(fn))
               if (nrow(D) < 30) print(fn)
               
               x[alg] <- mean(D[,2])
               y[alg] <- mean(D[,1])
      
               ci <- scale.ci * sd(D[,2]) / sqrt(nrow(D))
               ubx[alg] <- x[alg] + ci
               lbx[alg] <- x[alg] - ci
               
               ci <- scale.ci * sd(D[,1]) / sqrt(nrow(D))
               uby[alg] <- y[alg] + ci
               lby[alg] <- y[alg] - ci

               alg <- alg + 1
         }
      

       fn <-  paste(dir, "hiv_random.res", sep="")
       D <- as.matrix(read.table(fn))
       baseline <- mean(D[,1])
      
      x <- x / 60^2
      lbx <- lbx / 60^2
      ubx <- ubx / 60^2
      deltax <- deltax / 60^2
      
      col <- c("BLACK", "BLACK", rep("RED", length(mne)), rep("BLUE", length(mne)), rep(rgb(0,0.5,0), length(mne)), "VIOLET")

      xlim <- c(min(lbx) - deltax, max(ubx) + 1.5 * deltax)
      
      plot(NULL, xlim = xlim, ylim = c(baseline -deltay * 0.1, max(uby) + deltay), xlab = "Hours", ylab = "Return (log)", log=log, ...)

      for (i in 1:length(y)) {
          polygon(c(lbx[i], ubx[i], ubx[i], lbx[i], lbx[i]), c(lby[i], lby[i], uby[i], uby[i], lby[i]), col=col[i], border = col[i])
          }
      
      xx <- seq(xlim[1], xlim[2], length = 100)
      yy <- rep(baseline, 100)
      lines(xx, yy, lty = "dashed")
      text(deltax * 3, baseline + deltay*0.2, "Random Policy")
      
      prb <- c(1,2, 3, 4, 7, 8,9, 11) # c(1, 2, 6,7)
      txt <- c("200", "100", rep(mne, 3))
      text(x[-prb] + deltax, y[-prb] + deltay, txt[-prb], col = col[-prb])
      
      text(x[1] + deltax * 1.3, y[1], txt[1])
      text(x[2] + deltax * 1.3, y[2], txt[2])
      text(x[3] + deltax, y[3] - deltay, txt[3], col = "RED")
      text(x[4] + deltax, y[4] + 0.7 * deltay, txt[4], col = "RED")
      text(x[7], y[7] - 1.8 * deltay, txt[7], col = "BLUE")
      text(x[8] +  0 * deltax, y[8] + 1.8 * deltay, txt[8], col = "BLUE")
      text(x[9] + 0.8 * deltax, y[9] - 0.5 * deltay, txt[9], , col = "BLUE")
      text(x[11] + 0.8 * deltax, y[11] - 1.5 * deltay, txt[11], , col = rgb(0,0.5,0))
      
      txt <- c("FQIT(30 patients)", "TBSF(30 patients)", "TBSF(300 patients)", "TBSF(600 patients)")
      col <- unique(col)
      legend("right", txt, col = col, fill = col, pt.bg = col)
}


merge.results.fqt.hiv <- function
            (
             alg = "fqt",
             nt = 30,
             mne = c(200,100, 50, 10, 1),
             ncp = 8,
             it = 50,
             env = "hiv",
             dir = "~/hiv/fqt/",
             avg = 1:50
             )
   {
      
      
      for (i in 1:length(mne))
      {
         D <- NULL
     
         for (j in avg)
         {
            # FQT FIXED
            fn <- paste(dir, paste(env, alg, nt, mne[i], ncp, it, j, sep="_"), ".res", sep="")
            if (file.exists(fn)) D <- rbind(D,as.matrix(read.table(fn)))
            else print(fn)
         }
      
         fn <- paste(dir, paste(env, alg, nt, mne[i], ncp, it, sep="_"), ".res", sep="")
         
         pd(D)
         if (file.exists(fn)) print("FILE EXISTS! DID NOT SAVE")
         else wt(D, fn)
      }
   }
      
      
merge.results.tree.hiv <- function
            (
             alg = "tree",
             nt = 30,
             mne = c(200,100),
             ncp = 8,
             it = 50,
             nfi = 40,
             env = "hiv",
             dir = "~/hiv/fqt/",
             avg = 1:50
             )
   {
      
      
      for (i in 1:length(mne))
      {
         D <- NULL
     
         for (j in avg) 
         {
            
            fn <- paste(dir, paste(env, alg, nt, mne[i], ncp, it, nfi, j, sep="_"), ".res", sep="")
            D <- rbind(D,as.matrix(read.table(fn)))
         }
      
         fn <- paste(dir, paste(env, alg, nt, mne[i], ncp, it, nfi,  sep="_"), ".res", sep="")
         
         pd(D)
         if (file.exists(fn)) print("FILE EXISTS!")
         else wt(D, fn)
      }
   }
      

print("ftq.data.analysis.R loaded")