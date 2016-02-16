# # Experiments with archetypes in dynamic programming (NOT RL)
# The functions are roughly in the same order as used in the thesis

source("tese.plot.R")
source("mdp.R")
source("nmf.R")
source("ci.95.R")
source("dp.R")
source("arch.R")

  

gen.mp.simplex <- function(size, m, r.min = -1, r.max = 1, sd = 10, ...) {
# generates a Markov process in a (m-1)-dimensional simplex
   K  <- normal.transition.matrix(nrow = m, ncol = size, sd=sd, ...)
   rb <- runif(m,r.min, r.max)
   D <- matrix(runif(m*size), size, m)
   D <- D / apply(D, 1, sum)
   P <- D %*% K
   r <- D %*% rb
   list(r = r, P = P, K = K, D = D, rb = rb)
   }
   

ap.mp.simplex <- function(size = 100, p = 0.5, num.avg = 50,
                 nmf.function = nmf.mult, max.iter = 100, ...) {
# uses Lee and Seung's algorithm to approximate an artificial Markov process
# of size 'size' contained in a (m-1)-dimensional simplex
# returns the history of the approximation (so it's possible to plot it)
   his <- matrix(0, num.avg, max.iter)
   m <- round(p * size)
   for (i in 1: num.avg) {
      MP <- gen.mp.simplex(size, m, ...)
      #M <- cbind(MP$r,MP$P)
      F <- nmf.function(MP$P, m, max.iter, hist = TRUE)
      his[i,] <- F$hist
      }
   his <- his / size
   as.data.frame(his)
   }
   

ap.mp.simplex.script <- function(perc = c(0.1, 0.3, 0.5), sd=1,
                              max.iter=100, size = 100, num.avg = 50) {
   xh <- matrix(0, max.iter, length(perc))
   for (p in 1:length(perc)) {
      xh[,p] <- mean(ap.mp.simplex(size = size, p = perc[p], num.avg = num.avg,
                 nmf.function = nmf.mult, max.iter = max.iter, sd=sd))
      }
   xh
   }                 
   
   
plot.hist <- function(xh = NULL, m = c(10,30,50),
               eps.filename = "./fig_tese/xi_hist.eps") {
# plots the sse series of nmf.mult 
   if (is.null(xh)) xh <- read.table("./res_tese/dp/xh.txt")
   leg <- paste("m=", m)
   mp(xh, ylim= c(0,max(xh)), ylab = expression(xi[m](P^pi, DK)), 
      xlab = "Iterações", leg = leg, t="l", pch=NULL)
   dev.copy2eps(file = eps.filename)
   }

   
plot.normal <- function(sd = c(1,5,10,50), size = 100) {
# plots the theoretical normal distributions used in the thesis to generate the
# transition matrices 
   D <- matrix(0, size, length(sd))
   for (i in 1:length(sd)) {
      D[,i] <- normal(1:size, size %/% 2, sd[i])
      }
   mp(D)
   }
      
   
plot.dist.simplex <- function(sds = c(0.1,0.75,3), size = 20, noise = 0,
                        delta = 0.4, ncenter = 3,
                        dir = "./fig_tese/", 
                        save.to.eps = TRUE, show.leg = FALSE){
# plots points distributed in a 2-dimensional simplex according to a normal
# distribution
   cex = 2
   
   c <- cos(pi/3)
   s <- sin(pi/3)
   
   v <- matrix(c(0,0,c,s,1,0), 3, 2, byrow = TRUE)
   margins <- v + c(-delta, 0, delta, -delta, delta, 0)
   set.par.tese()
   par(mar=c(0,0,0,0))
   plot(margins, bty = "n", xlab = "", ylab = "", xaxt="n", yaxt="n", t="n",
         fig=c(0,1,0,1), new = TRUE)
   polygon(x = c(0,c,1), y = c(0,s,0))
   pos <- c(2,3,4)
   if (show.leg) pos <- c(2,4,4)
   text(x = c(0,c,1), y = c(0,s,0), labels = c("[1,0,0]", "[0,1,0]",
         "[0,0,1]"), pos=pos, offset = 0.5, cex = cex)
   
   col <- col.tese(length(sds)+1)[-1]
   for (s in 1:length(sds)) {
      D <- normal.transition.matrix(size, 3, ncenter, sds[s], noise = noise, 
      										lim = 1:ncenter)
      points(D %*% v, pch = s, col=col[s], cex=2.5)
      }
      
   if (show.leg) {
      leg(make.leg.widths(sds), pch =1:length(sds),col=col,
            cex=cex, pt.cex = 3, inset = 0, bty="n", x="topleft", 
	    lty="blank",lwd=2.5, y.intersp=1.2)
      }
      
   if (save.to.eps) {
      filename <- paste(dir, "simplex2d_n",noise,"_nc",ncenter, ".eps", sep="")
      dev.copy2eps(file = filename)
      }
   }
   
    
plot.dist.simplex.script <- function() {
# script to generate the plots just like in the thesis
   plot.dist.simplex(noise = 0.01, show.leg = TRUE)
   plot.dist.simplex(noise = 0.1, show.leg = FALSE)
   plot.dist.simplex(noise = 0.3, show.leg = FALSE)
   }
   
  

plot.dist.simplex.no.labels <- function(sds = c(0.1,0.75,3), size = 20, noise = 0,
                        delta = 0.4, ncenter = 3,
                        dir = "./fig_tese/", 
                        save.to.eps = TRUE, show.leg = FALSE){
# plots points distributed in a 2-dimensional simplex according to a normal
# distribution
   cex = 2
   
   c <- cos(pi/3)
   s <- sin(pi/3)
   
   v <- matrix(c(0,0,c,s,1,0), 3, 2, byrow = TRUE)
   margins <- v + c(-delta, 0, delta, -delta, delta, 0)
   set.par.tese()
   par(mar=c(0,0,0,0))
   plot(margins, bty = "n", xlab = "", ylab = "", xaxt="n", yaxt="n", t="n",
         fig=c(0,1,0,1), new = TRUE)
   polygon(x = c(0,c,1), y = c(0,s,0))
   pos <- c(2,3,4)
   if (show.leg) pos <- c(2,4,4)
#    text(x = c(0,c,1), y = c(0,s,0), labels = c("[1,0,0]", "[0,1,0]",
#          "[0,0,1]"), pos=pos, offset = 0.5, cex = cex)
   
   col <- col.tese(length(sds)+1)[-1]
   for (s in 1:length(sds)) {
      D <- normal.transition.matrix(size, 3, ncenter, sds[s], noise = noise, 
      										lim = 1:ncenter)
      points(D %*% v, pch = s, col=col[s], cex=3.5)
      }
      
   if (show.leg) {
      leg(make.leg.widths(sds), pch =1:length(sds),col=col,
            cex=cex, pt.cex = 3, inset = 0, bty="n", x="topleft", 
	    lty="blank",lwd=2.5, y.intersp=1.2)
      }
      
   if (save.to.eps) {
      filename <- paste(dir, "simplex2d_n",noise,"_nc",ncenter, ".eps", sep="")
      dev.copy2eps(file = filename)
      }
   }
   

plot.dist.simplex.script.paper <- function(noise =0.2,
	sds=c(0.15,0.75,3.00),size=15, delta = 0.1) {
# script to generate the plots just like in the paper
   plot.dist.simplex.no.labels(noise = noise, sds = sds, size=size, ncenter = 1, show.leg = TRUE, delta = delta)
   plot.dist.simplex.no.labels(noise = noise, sds = sds, size=size, ncenter = 2, delta = delta)
   plot.dist.simplex.no.labels(noise = noise, sds = sds, size=size, ncenter = 3, delta = delta)
   }
  

factor.P.mult <- function(max.iter = 100, size = 100, 
                          perc=0.2, 
                          p.centers = c(0.1,0.3,0.5,0.7,0.9),
                          sds = c(1,5,10,50), 
                          noises =c(0, 1e-4, 1e-3),
                          num.avg = 50, 
                          load.p=TRUE, 
                          dir ="./res_tese/dp/",
                          dir.res ="./res_tese/dp/", ... ) {
# factorate only the transition matrix P using "nmf.mult"
# the MPs were generated first, which may cause surprise
   ms <- perc * size
   nc <- p.centers * size
   for (c in 1:length(nc)) {
       print(paste("nc", c, "/", length(nc)))
       for (s in 1:length(sds)) {
         print(paste("sd", s, "/", length(sds)))
         for (n in 1:length(noises)) {
            print(paste("noise", n, "/", length(noises)))
            res <- matrix(0, num.avg, length(ms))
            for (i in 1:num.avg) {
               cat(paste(i, " ", sep=""))
               P <- NULL
               R <- NULL
               M <- NULL
               filename <- paste(dir, "P_norm_s", size, "_nc",nc[c], 
                                    "_sd", sds[s], "_n", noises[n], "_", i,
                                    ".txt", sep = "")
               if (!load.p) {
                  P <- normal.transition.matrix(size, ncenter = nc[c], 
                                             sd = sds[s], noise = noises[n])
                  wt(P, filename)
                  }
               else {
                  P <- as.matrix(read.table(filename))
                  }
                              
               for (m in 1:length(ms)) {
                  D <- nmf.mult(P, ms[m], max.iter = max.iter,...)
                  res[i,m] <- sum((P - D$D %*% D$K)^2)/nrow(P) 
                  }
               }
            
            for (m in 1:length(ms)) {
                  filename <- paste(dir.res, "mult_P_s", size, "_nc", nc[c],
                                    "_sd", sds[s], "_n", noises[n], "_p",
                                    perc[m], ".txt", sep="")
                  res.mat <- matrix(res[,m], num.avg,1)
                  colnames(res.mat) <- "P-DK"
                  write.table(res.mat, filename, quote = FALSE, row.names=FALSE,
                            col.names=TRUE)
                  }  
            print(" ")
            }
         }
      }
   print("All done!!")
   }
             
                                                                              
inf.norm <- function(x) {
   x <- abs(x)
   if (is.null(dim(x))) max(x)
   else {
      x <- apply(x, 1, sum)
      max(x)
      }
   }



## SE FOR RODAR ISSO DE NOVO TEM QUE CARREGAR 'M' DE ARQUIVO ##
factor.M <- function(factor.functions=c(nmf.mult, nmf.kmeans),
                        max.iters = c(100,10),
                        names = c("mult", "kmeans"),
                        sizes = c(100), 
                        p.centers = c(0.1,0.3,0.5,0.7,0.9), 
                        sds = c(1,5,10,50), 
                        noises = 0,
                        perc= 0.2,
                        df = c(0.1,0.3,0.5,0.7,0.9),
                        num.avg = 50, 
                        load.p=TRUE, 
                        save.m = TRUE,
                        dp.function= solve.mp,
                        min.reward = 0,
                        max.reward= 1, 
                        dir ="./res_tese/dp/",
                        dir.res = "./res_tese/dp/", ... ) {
# 'factor.functions' is a list with the functions to do the fatorization
# 'max.iter' is an array with the maximum number of iterations per function
# 'sizes' is the number of states |S| of each markov process (MP)
# 'sds' is the standard deviation to generate the transition matrix
# 'df' are the discount factors
# 'perc' is the percentual of reduction of the MP, i.e., k = perc * |S|
# 'num.avg' is the number of times each MP will be reduced
# 'names' is the strings with the names of the method (to write to file)
# 'load.mp' defines whether the MP should be created or loaded from file
# 'dp.function' is the function used to find the value function 
#    (currently, either "value.iteration" or "solve.mc")
# 'min.reward' is the minimum reward value possible in the MP
# 'max.reward' is the maximum reward value possible in the MP
# 'dir' is the directory where results should be saved

   for (s in 1:length(sizes)) {
     ncenters <- round(p.centers * sizes[s]) 
     for (nc in 1:length(ncenters)) {
      print(paste("nc", nc, "/", length(ncenters)))
      for (i in 1:length(sds)) {
            print(paste("sd", i, "/", length(sds)))
            for (n in 1:length(noises)) {
            print(paste("noise", n, "/", length(noises)))
            res <- array(0,c(length(factor.functions), length(perc), 
                        num.avg, 6 + 5 * length(df)))
                        
            for (j in 1:num.avg) {
               cat(paste(j, " ", sep=""))
               # create (or load) MP
               P <- NULL
               R <- NULL
               M <- NULL
               filename <- paste(dir, "P_norm_s", sizes[s], "_nc",ncenters[nc],
                                 "_sd", sds[i], "_n", noises[n], "_", j, ".txt",
                                 sep = "")
               if (!load.p) {
                  P <- normal.transition.matrix(sizes[s], ncenter =
                         ncenters[nc], sd = sds[i], noise = noises[n])
                  }
               else {
                  P <- as.matrix(read.table(filename))
                  }
               
               R <- runif(sizes[s], min.reward, max.reward)
               M <- cbind(R, P)
               
               if (save.m) {
                  filename <- paste(dir, "M_norm_s", sizes[s], "_nc",
                                    ncenters[nc], "_sd", sds[i], "_n",
                                    noises[n], "_", j, ".txt", sep = "")
                                    
                  wt(M,filename)
                  }                                    
                  
               
               Q <- matrix(0, length(R), length(df))
               Qc <- matrix(0, length(R), length(df))
               
               #print("Solving the MP exactly for...")
               for (d in 1:length(df)) {
                  # "solve" the MC exactly for each discount factor
                  #print(paste("Discount factor =", df[d]))
                  Q[,d] <- dp.function(R, P, df = df[d])
                  Qc[,d] <- Q[,d] - mean(Q[,d]) # center the data
                  }
            
               for (f in 1:length(factor.functions)) {
                  #print(paste("Running", names[f]))
                  for (p in 1:length(perc)) {
                     m <- round(perc[p] * sizes[s])
                    # print(paste("Perc.: ", perc[p], " (m =", m,")", sep=""))
                     
                     D<-factor.functions[[f]](M, m, max.iter = max.iters[f],
                                              P.only = FALSE, ...)
                     
                     PA <- D$D %*% D$K
                     RA <- D$D %*% D$r
                     
                     #Euclidean norm
                     res[f,p,j,1]<-sum((P - PA)^2)/nrow(P) 
                     res[f,p,j,2]<-sum((R - RA)^2)/length(R) 
                     # Max norm
                     res[f,p,j,3]<-inf.norm(P - PA)
                     res[f,p,j,4]<-inf.norm(R)
                     res[f,p,j,5]<-inf.norm(RA)
                     res[f,p,j,6]<-inf.norm(R - RA)
                     
                     P2 <- D$K %*% D$D
                     for (d in 1:length(df)) {
                        Q2 <- D$D %*% dp.function(D$r, P2, df = df[d])
                        Q2c <- Q2 - mean(Q2)
                        b <- 7 + (d - 1) * 5
                        # Euclidean norm
                        res[f,p,j,b]   <- sum((Q[,d] - Q2)^2) / length(Q2)
                        res[f,p,j,b+1] <- sum((Qc[,d] - Q2c)^2) / length(Q2c)
                        # Max norm
                        res[f,p,j,b+2] <- inf.norm(Q[,d])
                        res[f,p,j,b+3] <- inf.norm(Q2)
                        res[f,p,j,b+4] <- inf.norm(Q[,d]-Q2)
                        }
                     }
                  }
               }
            
            res.mat <- matrix(0, num.avg, 8 + 5 * length(df))
            for (f in 1:length(factor.functions)) {
               for (p in 1:length(perc)) {
                  filename <- paste(dir.res, names[f],"_M_s", sizes[s], 
                                    "_nc", ncenters[nc], "_sd", sds[i],
                                    "_n", noises[n], "_p", perc[p], ".txt",
                                    sep="")
                  res.mat <- res[f,p,,]
                  df.colnames <- outer(c("mse_V", "mse_mV",
                                           "max_V", "max_Va",
                                           "max_V_Va)"),df,paste,sep="_")
                  dim(df.colnames) <- NULL
                  colnames(res.mat) <- c("mse_P","mse_r",
                                         "max_P_Pa", "max_r",
                                         "max_ra", "max_r_ra",
                                         df.colnames)
                  write.table(res.mat, filename, quote = FALSE, row.names=FALSE,
                        col.names=TRUE)  
                  }
               }
            }
         print(" ")
         }
        }
      }
   print("All done!!")
   }
    


res.cube <- function(cols, 
                     stat.function = mean, 
                     names = c("mult", "kmeans"), 
                     sizes = 100, 
                     p.centers = c(0.1,0.3,0.5,0.7,0.9),
                     sds = c(1,5,10,50), 
                     noises=c(0,1e-4, 1e-3), 
                     perc=seq(0.05,0.5,length=10), 
                     dir = "./res_tese/dp/", 
                     type = "P") {

# creates an "hypercube" with the results            

   res <- array(0, c(length(names), length(sizes), length(perc),                
                    length(p.centers), length(sds),
                     length(noises), length(cols)))
   n.centers <- p.centers * sizes                   
   names(res) <- c("alg", "size", "m", "nc", "sd", "noise", "cols")
   for (n in 1:length(names)) 
      for (s in 1:length(sizes)) 
         for (p in 1:length(perc)) 
            for (c in 1:length(n.centers))   
               for (std in 1:length(sds))
                  for (ns in 1:length(noises)) {
                     filename <- paste(dir, names[n],"_", type,"_s", sizes[s], 
                                    "_nc", n.centers[c], "_sd",
                                    sds[std], "_n", noises[ns], "_p",
                                    perc[p], ".txt", sep="")
                     res[n,s,p,c,std,ns,] <- stat.function(read.table(filename,
                                             head=TRUE)[,cols])              
                  }
                  
   res   
   }                        
                        

  

plot.cube.2d <- function(res.cube, mask, plot.function = mp, ...) {
   ld <- NULL
   nrow <- NULL
   ncol <- NULL
   xind <- NULL
   yind <- NULL
   for (i in 1:7) {
      if (mask[i] > 0) ld <- c(ld, mask[i])
      else {
         ld <- c(ld,list(1:dim(res.cube)[i]))
         if (mask[i] == 0) { # key for x axis
            nrow <- dim(res.cube)[i]
            xind <- i
            }
         else { # key for the curve colors (< 0)
            ncol <- dim(res.cube)[i]
            yind <- i
            }
         }
      }
  byrow <- FALSE
  if (xind > yind) byrow <- TRUE
  table <- matrix(res.cube[ld[[1]], ld[[2]], ld[[3]], ld[[4]],
                           ld[[5]],ld[[6]], ld[[7]]],
                  nrow, ncol, byrow = byrow)
  plot.function(data.frame(table), ...)
  }
      
               
plot.P.script <- function(
                        size = 100,
                        p.centers = c(0.1,0.3,0.5,0.7,0.9),
                        sds=c(1,5,10,50),
                        noises = c(0,1e-4, 1e-3), 
                        type="Ponly",
                        cube = res.cube(1, perc=0.2, type=type,
                                       names="mult", noises=noises,
                                       sds=sds, p.centers = p.centers), 
                        dir= "./fig_tese/",
                        leg.pos =c("topleft", "bottomright", "topright",
                                   "right")
                        ){
# plots the results of factor.experiments() stored in "cube"
   
   #ylim <- c(min(cube), max(cube))
   n.centers = size * p.centers
   l <- make.leg.noises(noises)
   for (n in 1:length(sds)) {
      plot.cube.2d(cube, c(1,1,1,0,n,-1,1), x=n.centers, t="o", 
                     ylab=expression(xi[m](P^pi, DK)),
							xlab=expression(vartheta), 
                     leg=l, lwd=1.5, cex=2, cex.lab=1.2, leg.pos=leg.pos[n],
                     text.width = NULL)
      filename <- paste(dir, "mult_", type, "_sd",sds[n],".eps", sep="")
      dev.copy2eps(file = filename)
      }
   }


# This is to recover the total square error
mean.n <- function(x,n=100) mean(x*n)


plot.P.script.paper <- function(
                        size = 100,
                        p.centers = c(0.1,0.3,0.5,0.7,0.9),
                        sds=c(1,3,5),
                        noises = 1e-4, 
                        type="P",
                        cube = res.cube(1, perc=0.2, type=type,
                                       names="mult", noises=noises,
                                       sds=sds, p.centers = p.centers,
                                       stat.function = mean.n,
                                       dir="./res_paper/"), 
                        dir= "./fig_tese/",
                        leg.pos = "topleft"
                        ){
# plots the results of factor.experiments() stored in "cube"
   
   #ylim <- c(min(cube), max(cube))
   n.centers = size * p.centers
   l <- make.leg.widths(sds)
   plot.cube.2d(cube, c(1,1,1,0,-1,1,1), x=n.centers, t="o", 
               ylab=expression("SSE(P,DK)"),
	       xlab=expression(tilde(rho)), 
               leg=l, lwd=1.5, cex=2, cex.lab=1.2, leg.pos=leg.pos,
               text.width = NULL, y.intersp = 1.2)
   filename <- paste(dir, "res_sf_P.eps", sep="")
   dev.copy2eps(file = filename)
   }



plot.M.script <- function(
                        size = 100,
                        p.centers = c(0.1,0.3,0.5,0.7,0.9),
                        sds=c(1,5,10,50),
                        noises = c(0), 
                        type="M",
                        cube1 = res.cube(1, perc=0.2, type=type,
                                       names="mult", noises=noises,
                                       sds=sds, p.centers = p.centers), 
                        cube2 = res.cube(2, perc=0.2, type=type,
                                       names="mult", noises=noises,
                                       sds=sds, p.centers = p.centers), 
                        dir= "./fig_tese/",
                        noise.col=1
                        ){

# plots the difference of approximating a Markov process with and without the
# reward function.
# 'cube1' is the cube with results of the first column, 'cube2' has the values
# of the second column
   #ylim <- c(min(min(cube1),min(cube2)), max(max(cube1),max(cube2)))
   n.centers <- p.centers * size
   for (s in 1:length(sds)) {
      D <- cbind(cube1[1,1,1,,s,noise.col,1], cube2[1,1,1,,s,noise.col,1])
      D <- cbind(D[,1], D[,2], D[,1] + D[,2])
      mp(D, x=n.centers, t="o", ylab=expression(xi[m]),
         xlab=expression(vartheta), 
         leg=c(expression(P^pi), expression(r^pi), expression(M^pi)), 
         lwd=1.5, cex=2, cex.lab=1.2, leg.pos = "topleft")
      filename <- paste(dir, "mult_M_sd",sds[s],".eps", sep="")
      dev.copy2eps(file = filename)
      }
   }
         
         

plot.M.script.paper <- function(
                        size = 100,
                        p.centers = c(0.1,0.3,0.5,0.7,0.9),
                        sds=c(1,3,5),
                        noises = 1e-4, 
                        type="M",
                        cube1 = res.cube(1, perc=0.2, type=type,
                                       names="mult", noises=noises,
                                       sds=sds, p.centers = p.centers,
                                       stat.function = mean.n,
                                       dir="./res_paper/"), 
                        cube2 = res.cube(2, perc=0.2, type=type,
                                       names="mult", noises=noises,
                                       sds=sds, p.centers = p.centers,
                                       stat.function = mean.n,
                                       dir="./res_paper/"), 
                        dir= "./fig_tese/",
                        noise.col=1
                        ){

# plots the difference of approximating a Markov process with and without the
# reward function.
# 'cube1' is the cube with results of the first column, 'cube2' has the values
# of the second column
   #ylim <- c(min(min(cube1),min(cube2)), max(max(cube1),max(cube2)))
   n.centers <- p.centers * size
   for (s in 1:length(sds)) {
      D <- cbind(cube1[1,1,1,,s,noise.col,1], cube2[1,1,1,,s,noise.col,1])
      D <- cbind(D[,1] + D[,2], D[,1], D[,2])
      mp(D, x=n.centers, t="o", ylab="SSE",
         xlab=expression(tilde(rho)), 
         leg=c(expression(M), expression(P), expression(r)), 
         lwd=1.5, cex=2, cex.lab=1.2, leg.pos = "topleft", y.intersp=1.2)
      filename <- paste(dir, "res_sf_M_sd",sds[s],".eps", sep="")
      dev.copy2eps(file = filename)
      }
   }
         
         
plot.corr <- function(cols=c(17,22,27), 
                      dfs = c(0.5,0.7,0.9), 
                      ncenter = 50,
                      sd = 10, 
                      noise = 0, 
                      perc= 0.2, 
                      dir = "./res_tese/dp/", 
                      eps.dir = "./fig_tese/") {

# plots the correlation between the Mpi approximation error and the error on the
# computation of vpi  
   filename <- paste(dir, "mult_M_s100_nc", ncenter, "_sd", sd, "_n", noise,
                      "_p", perc, ".txt",sep="")
   D <- read.table(filename, h = TRUE)
   mse.v <- D[,cols]
   ylim <- c(min(mse.v), max(mse.v))
   
   mse.M <- D[,1] + D[,2]
   xlim <- c(min(mse.M), max(mse.M))
   x11(w=4,h=4)
   for (c in 1:length(cols)) {
      par(mar = c(5, 5, 2, 2) + 0.1)
      plot(mse.M, D[,cols[c]], ylim = ylim, xlim = xlim,
         xlab=expression(xi[m](M^pi, tilde(M)^pi)),
         ylab=expression(xi[m](v^pi, tilde(v)^pi)),
         #xlab="",
         #ylab="",
         cex=1.5, cex.axis=1.2, cex.lab=1.5)
      lines(xlim,xlim)
      cr <- round(cor(mse.M, D[,cols[c]]), d=2)
      eps.name <- paste(eps.dir,"corr_df", dfs[c], "_cor",cr, ".eps",sep="")
      dev.copy2eps(file = eps.name)
      }
  }
  
         

plot.error.bound<-function(
                  filename="./res_tese/dp/mult_M_s100_nc50_sd10_n0_p0.2.txt",
                  cols=c(11,16,21,26,31), 
                  dfs = c(0.1,0.3,0.5,0.7,0.9),
                  p.pa = 3, 
                  r.ra = 6, 
                  r = 4, 
                  dir="./fig_tese/") {
# plots the theoretical bounds and the actual errors  

   D <- read.table(filename, head=TRUE)
   
   x11(w=4,h=4)
   for (d in 1:length(dfs)) {
      coef <- 1/(1-dfs[d])
      bound <- coef * (D[,r.ra] + coef * dfs[d] * D[,p.pa] * D[,r])
      ylim = c(0, max(bound))
      print(paste(dfs[d], mean(bound)))
      print(paste(dfs[d], mean(D[,cols[d]])))
      #par(mar = c(5, 5, 2, 2) + 0.1)
      plot(D[,cols[d]], bound, ylim = ylim, 
           xlab=expression("||" * v^pi - tilde(v)^pi * "||"), 
           ylab="Cota superior", cex=1, cex.axis=1.2, cex.lab=1.5)
      x <- c(min(D[,cols[d]]),max(D[,cols[d]]))
      lines(x,x)
      dev.copy2eps(file=paste(dir,"lim_df",dfs[d],".eps", sep=""))
      }
   }
         
   
         

plot.v.script <- function( size = 100,
                           p.centers = c(0.1,0.5,0.9),
                           sds=c(1,5,10,50),
                           dfs=seq(0.1,0.9,l=5),
                           cube = res.cube(
                                    seq(7,27,by=5), 
                                    p.centers = p.centers,
                                    sds=sds,
                                    noise = 0,
                                    perc =0.2, 
                                    type = "M"),
                          dir="./fig_tese/") {
    n.centers = size * p.centers
    l <- make.leg.ncenters(n.centers)
    
    ylim <- c(min(cube[1,1,1, , ,1, ]), max(cube[1,1,1, , ,1, ]))
    
    for (sd in 1:length(sds)) {
       plot.cube.2d(cube, c(1,1,1,-1,sd,1,0), x=dfs, t="o", 
                      ylab=expression(xi[m](v^pi, tilde(v)^pi)),
                      xlab=expression(gamma), 
                      leg=l, lwd=1.5, cex=2, cex.lab=1.2, 
                      leg.pos="topleft")
       filename <- paste(dir, "v_sd",sds[sd],".eps", sep="")
       dev.copy2eps(file = filename)
       }
    }


         
plot.v.script.paper <- function( size = 100,
                           p.centers = c(0.9,0.5,0.1),
                           sds=c(1,3,5),
                           dfs=seq(0.1,0.9,l=5),
                           cube = res.cube(
                           			names = "mult",
                                    seq(7,27,by=5), 
                                    p.centers = p.centers,
                                    sds=sds,
                                    noise = 1e-4,
                                    perc =0.2, 
                                    type = "M",
                                    stat.function = mean.n,
                                    dir = "./res_paper/"),
                          dir="./fig_tese/") {
   n.centers = size * p.centers
   l <- expression()
   for(nc in n.centers) {
      l <- c(l, substitute(expression(tilde(rho) == nce), list(nce=nc))[[2]]) 
      }

 
    ylim <- c(min(cube[1,1,1, , ,1, ]), max(cube[1,1,1, , ,1, ]))
    
    for (sd in 1:length(sds)) {
       plot.cube.2d(cube, c(1,1,1,-1,sd,1,0), x=dfs, t="o", 
       					 ylab=expression("SSE(v,"*tilde(v)*")"),
                      xlab=expression(gamma), 
                      leg=l, lwd=1.5, cex=2, cex.lab=1.2, 
                      leg.pos="topleft", y.intersp=1.2)
       filename <- paste(dir, "res_sf_v_sd",sds[sd],".eps", sep="")
       dev.copy2eps(file = filename)
       }
	    
	 l <- make.leg.widths(sds)  
    for (pc in 1:length(p.centers)) {
       plot.cube.2d(cube, c(1,1,1,pc,-1,1,0), x=dfs, t="o", 
       					 ylab=expression("SSE(v,"*tilde(v)*")"),
                      xlab=expression(gamma), 
                      leg=l, lwd=1.5, cex=2, cex.lab=1.2, 
                      leg.pos="topleft",y.intersp=1.2)
       filename <- paste(dir, "res_sf_v_sr",n.centers[pc],".eps", sep="")
       dev.copy2eps(file = filename)
       }
    }
         
         
plot.nmf.kmeans.script <- function(size = 100,
                           p.centers = 0.5,
                           sds=1,
                           dfs=0.9,
                           perc = c(0.1,0.3,0.5,0.7,0.9),
                           cols = c(1,2,27), ## ALWAYS CHECK ##
                           cube = res.cube(
                                    cols = cols, 
                                    names = c("mult", "kmeans"),
                                    p.centers = p.centers,
                                    sds=sds,
                                    noise = 0,
                                    perc = perc, 
                                    type = "M",
                                    stat.function = mean),
                          dir="./fig_tese/") {
    n.centers = size * p.centers
    n.arq = size * perc
    l <- c("Lee e Seung", "K-means")
    
    D <- cbind(cube[1,1,,1,1,1,1]+ cube[1,1,,1,1,1,2], 
               cube[2,1,,1,1,1,1]+ cube[2,1,,1,1,1,2])
    ylim <- c(0, max(D))
    
    
    plot.cube.2d(cube, c(-1,1,0,1,1,1,1), x=n.arq, t="o", 
                       ylab=expression(xi[m](P^pi, tilde(P)^pi)),
                       xlab="m", ylim = ylim,
                       leg=l, lwd=1.5, cex=2, cex.lab=1.2, 
                       leg.pos="topright")
    filename <- paste(dir, "mult_kmeans_P.eps", sep="")
    dev.copy2eps(file = filename)
    
    mp(D, x=n.arq, t="o", ylab=expression(xi[m](M^pi, tilde(M)^pi)),
         xlab="m",  leg=l, lwd=1.5, cex=2, cex.lab=1.2,
         leg.pos = "topright", ylim = ylim)
    filename <- paste(dir, "mult_kmeans_M.eps", sep="")
    dev.copy2eps(file = filename)
    
    
    plot.cube.2d(cube, c(-1,1,0,1,1,1,3), x=n.arq, t="o", 
                       ylab=expression(xi[m](v^pi, tilde(v)^pi)),
                       xlab="m", 
                       leg=l, lwd=1.5, cex=2, cex.lab=1.2, 
                       leg.pos="topright")
    filename <- paste(dir, "mult_kmeans_v.eps", sep="")
    dev.copy2eps(file = filename)
    }


save.MDP <- function(P, R, filename) {
	num.actions <- dim(P)[3]
	BP <- matrix(0, nrow(P) * num.actions, ncol(P))
	BR <- array(0, nrow(P) * num.actions)
	
	for (a in 1:num.actions) {
		b <- (a - 1) * nrow(P) + 1
		e <- b + nrow(P) - 1
		BP[b:e,] <- P[,,a]
		BR[b:e] <- R[,a]
		}
	wt(cbind(BP,BR), filename)
	}
	

load.MDP <- function(num.states, num.actions, filename) {
	P <- array(0, c(num.states, num.states, num.actions))
	R <- matrix(0, num.states, num.actions)
	BP <- as.matrix(read.table(filename))
	BR <- BP[,ncol(BP)]
	BP <- BP[,-ncol(BP)]
	
	for (a in 1:num.actions) {
		b <- (a - 1) * nrow(P) + 1
		e <- b + nrow(P) - 1
		P[,,a] <- BP[b:e,]
		R[,a] <- BR[b:e]
		}
	list(P = P, R = R)
	}
	

#  factor.MDP(nmf.mult, 100, "mult", sd=3, noise=1e-4, perc=seq(0.1,0.9,by=0.1),
# num.actions=5:9, dir="./res_paper/")

factor.MDP <- function(factor.function= nmf.kmeans,
                        max.iter.ff = 10,
                        name = "kmeans",
                        size = 100, 
                        p.center = 0.5, 
                        sd = 5, 
                        noise = 0,
                        perc= c(0.1, 0.3, 0.5),
                        dfs = c(0.3, 0.5, 0.7, 0.9),
                        num.actions = 2:10,
                        max.iter.pi = 20,
                        num.avg = 50, 
                        solve.mp.function= solve.mp,
                        min.reward = 0,
                        max.reward= 1, 
                        save.model = FALSE,
                        load.model = TRUE,
                        dir ="./res_tese/dp/",... ) {
# 'factor.functions' is a list with the functions to do the fatorization
# 'max.iter' is an array with the maximum number of iterations per function
# 'sizes' is the number of states |S| of each markov process (MP)
# 'sds' is the standard deviation to generate the transition matrix
# 'df' are the discount factors
# 'perc' is the percentual of reduction of the MP, i.e., k = perc * |S|
# 'num.avg' is the number of times each MP will be reduced
# 'names' is the strings with the names of the method (to write to file)
# 'load.mp' defines whether the MP should be created or loaded from file
# 'dp.function' is the function used to find the value function 
#    (currently, either "value.iteration" or "solve.mc")
# 'min.reward' is the minimum reward value possible in the MP
# 'max.reward' is the maximum reward value possible in the MP
# 'dir' is the directory where results should be saved
   ncenters <- round(p.center * size)
   P.pi <- matrix(0, size, size)
   R.pi <- array(0, size)
   
   for (a in num.actions) {
      print(paste("action", a))
      res <- array(0, c(num.avg, length(perc), length(dfs), 3))
      for (i in 1:num.avg) {
         cat(paste(i, " ", sep=""))
         
         P <- NULL
         R <- NULL
         
         # Generate the MDP
         if (load.model) {
         	filename <- paste(dir, "MDP_s", size, "_a", a, "_nc", ncenters, 
         						"_sd", sd, "_n", noise, "_", i, ".txt", sep="")
         	M <- load.MDP(size, a, filename)
         	R <- M$R
         	P <- M$P
         	}
         else {
            P <-mdp.transition.matrix(size, a, ncenters = ncenters, 
                                 sd= sd, noise =noise)
         	R <- matrix(runif(size * a, min.reward, max.reward), 
                          size, a)
				}           

         # Save the MDP to use it later
         if (save.model) {
         	filename <- paste(dir, "MDP_s", size, "_a", a, "_nc", ncenters, 
         						"_sd", sd, "_n", noise, "_", i, ".txt", sep="")
         	save.MDP(P, R, filename)
         	}
         
         # Solve the MDP
         tr <- NULL
         for (df in dfs) {
            tr <- c(tr,list(policy.iteration(R, P, df, max.iter = max.iter.pi,
                     solve.mp.function = solve.mp.function)))
                     
            }
    
         for (p in 1:length(perc)) {
            m <- round(perc[p] * size)
            
            for (d in 1:length(dfs)) {
               
               ap <- arch.policy.iteration(R, P, dfs[d], m, max.iter =
                        max.iter.pi, solve.mp.function = solve.mp.function, 
                        factor.function = factor.function, max.iter.ff =
                        max.iter.ff)
                        
               res[i,p,d,1] <- sum(tr[[d]]$pi != ap$pi) / size
               
              
              # compute the value function of the true and approximate policy
              v.true <- apply(tr[[d]]$Q,1,max)
              for (s in 1:size) {
                  R.pi[s] <- R[s,ap$pi[s]]
                  P.pi[s, ] <- P[s,, ap$pi[s]]
                  }
               v.ap <- solve.mp.function(R.pi, P.pi, dfs[d], ...)
               
               res[i,p,d,2] <- sum((v.true - v.ap)^2) / size
               res[i,p,d,3] <- ap$mse ##wasting memory##
               }
            }
         }
      print(" ")
      
#        print(apply(res[,1,1,],2,mean))
#        print(apply(res[,2,1,],2,mean))
#        print(apply(res[,3,1,],2,mean))
       
       for (p in 1:length(perc)) {
          filename <- paste(dir, name, "_MDP_s", size, "_nc", ncenters,
                                     "_sd", sd, "_n", noise, "_p",
                                     perc[p], "_a", a, ".txt",
                                     sep="")
          res.mat <- cbind(matrix(res[,p,,1], num.avg, length(dfs)),
                           matrix(res[,p,,2], num.avg, length(dfs)),
                           res[,p,1,3]
                           )
          colnames(res.mat) <- c(paste("pi_err_", dfs, sep=""), 
                                 paste("r_err_",  dfs, sep=""),
                                 "M_err")
          write.table(res.mat, filename, quote = FALSE, row.names=FALSE,
                             col.names=TRUE)
          }
      }
   print("All done!!!")
   }
 
 

solve.MDP.straits <- function(  size = 100, 
                        num.actions = 5,
                        prob.pass = 0.5, 
                        df = 0.7,
                        max.iter.pi = size+10,
                        solve.mp.function= solve.mp,
                        num.avg = 50, 
                        save.model = FALSE,
                        load.model = TRUE,
                        dir ="./res_paper/tmp/",... ) {
   
   times <- NULL
   
   for (a in num.actions) {
   	
      print(paste("action", a))
      
      policies <- matrix(0,num.avg, size)
      times    <- matrix(0, num.avg, 1) 
      
      for (i in 1:num.avg) {
         cat(paste(i, " ", sep=""))
        	
        	R <- NULL
        	P <- NULL
         
         # Generate  (or load) the MDP
         if (load.model) {
         	filename <- paste(dir, "MDP_straits_s", size, "_a", a, 
         						"_pp", prob.pass, "_", i,".txt", sep="")
         	M <- load.MDP(size, a, filename)
         	R <- M$R
         	P <- M$P
         	}
         else {
  				MDP <- mdp.straits.no.goal(size, a, prob.pass = prob.pass)
  				P <- MDP$P                     
  				R <- MDP$R
     			}           

         # Save the MDP to use it later
         if (save.model) {
         	filename <- paste(dir, "MDP_straits_s", size, "_a", a, 
         						"_pp", prob.pass, "_", i,".txt", sep="")
         	save.MDP(P, R, filename)
         	}
         
         # Solve the MDP and save the time spent
         pi <- array(1,size) 
         times[i,] <- system.time(policies[i,] <- policy.iteration(R, P, df, 
          max.iter= max.iter.pi,solve.mp.function = solve.mp.function,
			 pi=pi)$pi, 
          TRUE)[1]
         }
         
   	filename <- paste(dir, "policy_straits_s", size, "_a", a, 
         						"_pp", prob.pass,".txt", sep="")
      wt(policies, filename)
      
      filename <- paste(dir, "times_straits_s", size, "_a", a, 
         						"_pp", prob.pass,".txt", sep="")
      wt(times, filename)
		}
	}



solve.MDP.straits.sf <- function(factor.function= nmf.kmeans,
                        max.iter.ff = 10,
                        name = "kmeans",
                        size = 100, 
                        num.actions = 5,
                        prob.pass = 0.5,
                        df = 0.7,
			perc= 0.2,
                        num.avg = 50, 
                        solve.mp.function= solve.mp,
                        dir ="./res_paper/tmp/",... ) {
   
   res <- NULL
   for (a in num.actions) {
      print(paste("action", a))
      
      res <- array(0, c(length(perc), num.avg, 3))
      
      # load the decision policies
   	filename <- paste(dir, "policy_straits_s", size, "_a", a, 
         						"_pp", prob.pass,".txt", sep="")
      policies <- read.table(filename)
      
      for (i in 1:num.avg) {
         cat(paste(i, " ", sep=""))
         
         # load the MDP 
        	filename <- paste(dir, "MDP_straits_s", size, "_a", a, 
         						"_pp", prob.pass, "_", i,".txt", sep="")
         P <- load.MDP(size, a, filename)
        	R <- P$R
        	P <- P$P
         
         # Solve the MDP using PISF
  
         for (p in 1:length(perc)) {
            m <- round(perc[p] * size)
            
            pi <- array(1, size)
            t <- system.time(
            Q <- PISF(R, P, df, m, pi = pi, max.iter = m+10, solve.mp.function
				=solve.mp.function, factor.function =factor.function,	
				max.iter.ff = max.iter.ff), TRUE)[1]
				pi <- Q$pi
				tf <- Q$factor.time                    
				                    
            res[p,i,1] <- sum(pi != policies[i,]) / size
            res[p,i,2] <- t   
            res[p,i,3] <- tf
            }
      }

      # save results
      for (p in 1:length(perc)) {
         	filename <- paste(dir, name,"_straits_s", size, "_a", a, 
         						"_pp", prob.pass, "_p" ,perc[p], ".txt", sep="")
         wt(res[p,,], filename)
          }
      print(" ")
      }
   res   
   }
   
   
script.straits <- function(sizes = seq(100,500,by=100),
                        	num.actions = 5,
                        	prob.pass = 0.5,
                        	df = 0.7,
                        	factor.functions = list(nmf.kmeans),
                        	names = c("kmeans"),
                        	max.iter.ffs = 10,
				perc= 0.2,
                        	num.avg = 50, 
                        	run.PI = TRUE,
                        	dir = "./res_paper/tmp/"){
									   
	for (s in sizes) {
		if (run.PI) {
         solve.MDP.straits(s, num.actions, prob.pass, df,
                        num.avg = num.avg, save.model = TRUE,
                        load.model = FALSE, dir = dir)
         }
		
		for (f in 1:length(factor.functions)) {
			solve.MDP.straits.sf(factor.functions[[f]], max.iter.ffs[f],
                        names[f], s, num.actions, prob.pass,
                        df, perc, num.avg, dir = dir)
         }
      }
  print("All done!")
  }
                        

plot.straits <- function(sizes = seq(100,500,by=100),
                        	num.actions = 5,
                        	prob.pass = 0.5,
                        	df = 0.7,
                        	names = c("kmeans"),
                        	names.leg = c("Policy iteration","PISF + k-means"),
                        	perc= 0.2,
			        num.avg = 50,
                        	dir = "./res_paper/tmp/",
                        	dir.eps = "./fig_tese/",
                        	eps.filename = c("times.eps", "errors.eps")){
	
	T <- matrix(0, length(sizes), length(names)+1)
	R <- matrix(0, num.avg, length(sizes))
	
	for (i in 1:length(sizes)) {
	    filename <- paste(dir, "times_straits_s", sizes[i], "_a", num.actions, 
         						"_pp", prob.pass,".txt", sep="")
	    T[i,1] <- mean(read.table(filename))
      
	    for (j in 1:length(names)) {
        	filename <- paste(dir, names[j],"_straits_s", sizes[i], "_a",
		  num.actions,  "_pp", prob.pass,"_p" ,perc, ".txt",sep="")
		A <- read.table(filename)       
	        T[i,1+j] <- mean(A[,2])
		R[,i] <- A[,1]     ##huumm, not good...
		}
	 }
   
   mp(T, ylab = "s", xlab = "|S|", leg = names.leg, t="b", x=sizes,
      leg.pos = "topleft",y.intersp=1.2, lwd=1.5, cex=2, cex.lab=1.2)
   dev.copy2eps(file = paste(dir.eps,eps.filename[1],sep=""))
   
   x11()
   # revert order of R's columns
#    R2 <- R[,ncol(R)]
#    for (c in (ncol(R)-1):1) {
#    	R2 <- cbind(R2,R[,c])
#    	}
#    	
#    mp(R2, ylab = expression(chi), xlab = "|S|", leg = rev(names.leg[-1]), t="b",
# 		x=sizes, col = rev(col.tese(ncol(R)+1)[-1]), lty=rev(2:(ncol(R)+1)),
# 		leg.pos = "bottomright",y.intersp=1.2, lwd=1.5, cex=2, cex.lab=1.2)
      R <- data.frame(R)
      colnames(R) <- paste(sizes)
      means <- apply(R,2,mean)
      lim.sup <- apply(R,2,max) - means
      lim.inf <- means - apply(R,2,min)
      set.par.tese()
      col <- col.tese(2)[2]
      plotCI(means, uiw = lim.sup, liw = lim.inf, xlab="|S|", ylab=expression(chi), 		ps=20, cex=2, pch = 22, pt.bg = col, slty="dashed", x=sizes,
	      col = col)
#       
#       stripchart(R,  xlab="|S|", ylab=expression(chi), ps=20, cex=2, vertical=TRUE)
     dev.copy2eps(file = paste(dir.eps,eps.filename[2],sep=""))
   }
	      

plot.straits.box <- function(sizes = seq(100,500,by=100),
                        	num.actions = 5,
                        	prob.pass = 0.5,
                        	df = 0.7,
                        	names = c("kmeans","heuristic"),
                        	names.leg = c("PI", "K-means","Heuristic"),
                        	perc= 0.2,
                        	dir = "./res_paper/tmp/"){
	
	T <- matrix(0, 2, (length(names))*length(sizes))
	R <- matrix(0, length(sizes), length(names))
	
	for (i in 1:length(sizes)) {
		b <- (i-1) * (length(names))          						
      for (j in 1:length(names)) {
        	filename <- paste(dir, names[j],"_straits_s", sizes[i], "_a",
								num.actions,  "_pp", prob.pass,"_p" ,perc, 
								".txt",sep="")
			A <- read.table(filename)       
      	T[1,b+j] <- mean(A[,2])
      	T[2,b+j] <- mean(A[,3])
      	R[i,j] <- mean(A[,1])
      	}
      }
   barplot(T)
   #mp(T, ylab = expression(ms), xlab = "|S|", leg = names.leg, t="b", x=sizes)
   #x11()
   #mp(R, ylab = expression(chi), xlab = "|S|", leg = names.leg[-1], t="b",
	#	x=sizes)
   }

		                        	

factor.MDP.heuristic <- function(
                        num.bits = 10, 
                        dif.min = 1e-1,
                        df = 0.999,
                        noise = 0,
                        flip.bit.function = flip.bit.li
                        ) {
         
     M <-mdp.bitflip(num.bits, flip.bit.function, noise = noise)
     P <- M$P
     R <- M$R

     time.pi <- system.time(
     tr <- policy.iteration(R, P, df, max.iter = nrow(R)), TRUE)[1]
          
     time.pisf <- system.time( 
          ap <- PISF(R, P, df, 0, max.iter = nrow(R), 
               factor.function = nmf.fast.heuristic2, max.iter.ff = 0, 
               n = nrow(R), dif.min = dif.min), TRUE)[1]
      list(e = sum(tr$pi != ap$pi) / nrow(R),
           m = ap$m,
           t = time.pisf / time.pi)
    }

script.heuristic <- function(num.bits = 2:10, dif.min = seq(0,0.5, by=0.1),
                              dir = "./res_paper/", name = "heuristic2",
                              verbose = FALSE) {
     res.m <- matrix(0, length(num.bits), length(dif.min))
     res.e <- matrix(0, length(num.bits), length(dif.min))

     for (i in 1:length(num.bits)) {
          for (j in 1:length(dif.min)) {
               if (verbose) print(paste("Bits: ", num.bits[i],
                                   " Dif: ", dif.min[j]))
               L <- factor.MDP.heuristic(num.bits[i], dif.min[j])
               res.m[i,j] <- L$m
               res.e[i,j] <- L$e
               print(paste(2^num.bits[i],L$m,L$e))
               }
           }
     filename <- paste(dir,name,"_m_bm", num.bits[1], "_bM",
          num.bits[length(num.bits)], "_dm", dif.min[1],
          "_dM", dif.min[length(dif.min)],".txt", sep="")
     wt(res.m, filename)

     filename <- paste(dir,name,"_e_bm", num.bits[1], "_bM",
          num.bits[length(num.bits)], "_dm", dif.min[1],
          "_dM", dif.min[length(dif.min)],".txt", sep="")
     wt(res.e, filename)

     list(M = res.m, E = res.e)
     }
									

plot.heuristic <- function(res, err, num.bits = 2:10, dif.min = seq(0,0.5,
by=0.1), delta1 = 10, delta2 = -50) {
# To reproduce the figure in the paper:
#    x11(h=6.5, w=10)
#        res <- read.table("./res_paper/heuristic2_m_bm2_bM10_dm0_dM0.5.txt")
#         err <-read.table("./res_paper/heuristic2_e_bm2_bM10_dm0_dM0.5.txt")
# plot.heuristic(res[5:9 ,c(1,4,6)],err[5:9 ,c(1,4,6)], num.bits=6:10, dif.min
# =c(0,0.3,0.5), delta2=-40,delta1=0)
   
 res <- cbind(2^(num.bits), res)
     mp(res, x=num.bits, t="o", ylab="Number of states",
         xlab=expression("Number of bits"), 
         leg=c("Original MDP", make.leg.difs(dif.min)) , 
         lwd=1.5, cex=2, cex.lab=1.2, leg.pos = "topleft", y.intersp=1.2,
          ylim = c(min(res)-25, max(res)))
 for (i in 1:nrow(err)) {
     for (j in 1:ncol(err)) {
          if (err[i,j] > 0) {
                text(num.bits[i] + delta1, res[i,j+1] + delta2, 
                paste(round(err[i,j],2)), col=col.tese(nrow(err))[j+1],cex=0.9)
#            print(paste(i,j))
#            text(locator(1), 
#                 paste(round(err[i,j],2)), col=col.tese(nrow(err))[j+1],cex=1)
            }

          }
     }
}

factor.MDP.sep <- function(factor.function= nmf.kmeans,
                        max.iter.ff = 10,
                        name = "kmeans",
                        size = 100, 
                        p.center = 0.5, 
                        sd = 5, 
                        noise = 0,
                        perc= c(0.1, 0.3, 0.5),
                        dfs = c(0.3, 0.5, 0.7, 0.9),
                        num.actions = 2:10,
                        max.iter.pi = 20,
                        num.avg = 50, 
                        solve.mp.function= solve.mp,
                        min.reward = 0,
                        max.reward= 1, 
                        load.model = TRUE,
                        dir ="./res_tese/dp/",... ) {
# 'factor.functions' is a list with the functions to do the fatorization
# 'max.iter' is an array with the maximum number of iterations per function
# 'sizes' is the number of states |S| of each markov process (MP)
# 'sds' is the standard deviation to generate the transition matrix
# 'df' are the discount factors
# 'perc' is the percentual of reduction of the MP, i.e., k = perc * |S|
# 'num.avg' is the number of times each MP will be reduced
# 'names' is the strings with the names of the method (to write to file)
# 'load.mp' defines whether the MP should be created or loaded from file
# 'dp.function' is the function used to find the value function 
#    (currently, either "value.iteration" or "solve.mc")
# 'min.reward' is the minimum reward value possible in the MP
# 'max.reward' is the maximum reward value possible in the MP
# 'dir' is the directory where results should be saved
   ncenters <- round(p.center * size)
   P.pi <- matrix(0, size, size)
   R.pi <- array(0, size)
   
   for (a in num.actions) {
      print(paste("action", a))
      res <- array(0, c(num.avg, length(perc), length(dfs), 3))
      for (i in 1:num.avg) {
         cat(paste(i, " ", sep=""))
         
         P <- NULL
         R <- NULL
         if (load.model) {
         	filename <- paste(dir, "MDP_s", size, "_a", a, "_nc", ncenters, 
         						"_sd", sd, "_n", noise, "_", i, ".txt", sep="")
         	M <- load.MDP(size, a, filename)
         	R <- M$R
         	P <- M$P
         	}
         else {
            P <-mdp.transition.matrix(size, a, ncenters = ncenters, 
                                 sd= sd, noise =noise)
         	R <- matrix(runif(size * a, min.reward, max.reward), 
                           size, a)
				}           
         
         tr <- NULL
         for (df in dfs) {
            tr <- c(tr,list(policy.iteration(R, P, df, max.iter = max.iter.pi,
                     solve.mp.function = solve.mp.function)))
            }
    
         for (p in 1:length(perc)) {
            m <- round(perc[p] * size)
            
            # Factorate the MDP
            Pb <- array(0,c(m,m,a))
            Rb <- matrix(0, m, a)
	         D <- NULL
	         
            M.mse <- 0         
            for (u in 1:a) {
             	M <- cbind(R[,u], P[,,u])
            	D <- factor.function(M, m, max.iter = max.iter.ff,
                                  P.only = FALSE, ...)
													 
            	Pb[,,u] <- D$K %*% D$D
            	Rb[,u]  <- D$r
            	
            	M.mse <- M.mse + sum((D$D %*% D$K - P[,,u])^2)
            	M.mse <- M.mse + sum((D$D %*% D$r - R[,u])^2)
            	}
            	
				M.mse <- M.mse / size
            
            for (d in 1:length(dfs)) {
               
               ap <- policy.iteration(Rb, Pb, df, max.iter = max.iter.pi,
                     solve.mp.function = solve.mp.function)
               
               ap$Q <- D$D %*% ap$Q
               ap$pi <- apply(ap$Q, 1, which.max)
               
               res[i,p,d,1] <- sum(tr[[d]]$pi != ap$pi) / size
               
              
              # compute the value function of the true and approximate policy
              v.true <- apply(tr[[d]]$Q,1,max)
              for (s in 1:size) {
                  R.pi[s] <- R[s,ap$pi[s]]
                  P.pi[s, ] <- P[s,, ap$pi[s]]
                  }
               v.ap <- solve.mp.function(R.pi, P.pi, dfs[d], ...)
               
               res[i,p,d,2] <- sum((v.true - v.ap)^2) / size
               res[i,p,d,3] <- M.mse ##wasting memory##
               }
            }
         }
      print(" ")
      
#        print(apply(res[,1,1,],2,mean))
#        print(apply(res[,2,1,],2,mean))
#        print(apply(res[,3,1,],2,mean))
       
       for (p in 1:length(perc)) {
          filename <- paste(dir, name, "_sep_MDP_s", size, "_nc", ncenters,
                                     "_sd", sd, "_n", noise, "_p",
                                     perc[p], "_a", a, ".txt",
                                     sep="")
          res.mat <- cbind(matrix(res[,p,,1], num.avg, length(dfs)),
                           matrix(res[,p,,2], num.avg, length(dfs)),
                           res[,p,1,3]
                           )
          colnames(res.mat) <- c(paste("pi_err_", dfs, sep=""), 
                                 paste("r_err_",  dfs, sep=""),
                                 "M_err")
          write.table(res.mat, filename, quote = FALSE, row.names=FALSE,
                             col.names=TRUE)
          }
      }
   print("All done!!!")
   }


plot.MDP <- function(delta = 0, #delta = 0 for pi error, delta = 1 for q error
                     stat.function = mean,
                     name = "kmeans",
                     size = 100, 
                     p.center = 0.5, 
                     sd = 5, 
                     noise = 0,
                     perc= c(0.1,0.3,0.5),
                     dfs = c(0.3, 0.5, 0.7, 0.9),
                     num.actions = 2:10,
                     dir ="./res_tese/dp/",
                     dir.eps = "./fig_tese/",
                     suffix = "") {   
   
   ncenters <- round(size * p.center)
   m = round(perc * size)
   
   D <- array(0, c(length(dfs), length(num.actions), length(perc)))
   
   for (d in 1:length(dfs)) {
      for (a in 1:length(num.actions)) {
         for (p in 1:length(perc)) {
            filename <- paste(dir, name, "_MDP_s", size, "_nc",ncenters,
                                       "_sd", sd, "_n", noise, "_p",
                                       perc[p], "_a", num.actions[a], ".txt",
                                       sep="")
            b <- delta * length(dfs) + d
            D[d,a,p] <- stat.function(read.table(filename, head = TRUE)[,b])
            }
         }
      }
      
    for (d in 1:length(dfs)) {
      #ylim = c(min(D), max(D))
       
       mp(D[d,,], x=num.actions, t="o", ylab=expression(chi[m]),
          xlab="|A|",leg=paste("m =",m), 
          lwd=1.5, cex=2, cex.lab=1.2, leg.pos = "topleft")
       
      filename <- paste(dir.eps, name, "_MDP", suffix,"_df",dfs[d],".eps",
                        sep="")
      dev.copy2eps(file = filename)
      }
     D
  }

                           
plot.MDP.paper <- function(delta = 0, #delta = 0 for pi error, delta = 1 for
                     stat.function = mean,
                     name = "mult",
                     size = 100, 
                     p.center = 0.5, 
                     sd = 3, 
                     noise = 1e-4,
                     perc= seq(0.1,0.9,by=0.1),
                     num.actions = 2:10,
                     a.ind = 4,
                     m.ind = 2,
                     dir ="./res_paper/",
                     dir.eps = "./fig_tese/",
                     suffix = "") {   
   
   ncenters <- round(size * p.center)
   m <- round(perc * size)
   dfs <- c(0.9,0.7,0.5) ## this can't be changed
   
   D <- array(0, c(length(dfs), length(num.actions), length(perc)))
   
   for (d in 1:length(dfs)) {
      for (a in 1:length(num.actions)) {
         for (p in 1:length(perc)) {
            filename <- paste(dir, name, "_MDP_s", size, "_nc",ncenters,
                                       "_sd", sd, "_n", noise, "_p",
                                       perc[p], "_a", num.actions[a], ".txt",
                                       sep="")
            b <- 5-d
            D[d,a,p] <- stat.function(read.table(filename, head = TRUE)[,b])
            }
         }
      }
      
     
      mp(t(D[,,m.ind]), x=num.actions, t="o", ylab=expression(chi),
          xlab="|A|",leg=make.leg.dfs(dfs), 
          lwd=1.5, cex=2, cex.lab=1.2, leg.pos = "topleft")
       
     filename <- paste(dir.eps, "res_sf_MDP", suffix,"_m",m[m.ind],".eps",
                        sep="")
      
      dev.copy2eps(file = filename)
      
      
      
      mp(t(D[,a.ind,]), x=m, t="o", ylab=expression(chi),
          xlab="m",leg=make.leg.dfs(dfs), 
          lwd=1.5, cex=2, cex.lab=1.2, leg.pos = "topright")
       
     filename <- paste(dir.eps, "res_sf_MDP", suffix,"_a",num.actions[a.ind],
								".eps",sep="")
      dev.copy2eps(file = filename)
   D        
  }



plot.MDP.paper.sep <- function(delta = 0, #delta = 0 for pi error, delta = 1 for
                     stat.function = mean,
                     name = "mult",
                     size = 100, 
                     p.center = 0.5, 
                     sd = 3, 
                     noise = 1e-4,
                     perc= seq(0.1,0.9,by=0.1),
                     num.actions = 2:10,
                     df = 0.9,
                     fix.perc = 0.2,
                     fix.a = 5,
                     col1 = 4,
                     col2 = 1,
                     dir ="./res_paper/",
                     dir.eps = "./fig_tese/",
                     suffix = "") {   
   
   ncenters <- round(size * p.center)
   m <- round(perc * size)
     
   D <- matrix(0, length(num.actions), 2)
   
   for (a in 1:length(num.actions)) {
      filename <- paste(dir, name, "_MDP_s", size, "_nc",ncenters,
                        "_sd", sd, "_n", noise, "_p", fix.perc, "_a", 
                        num.actions[a], ".txt", sep="")
      D[a,1] <- stat.function(read.table(filename, head = TRUE)[,col1])
      
      filename <- paste(dir, name, "_sep_MDP_s", size, "_nc",ncenters,
                        "_sd", sd, "_n", noise, "_p", fix.perc, "_a", 
                        num.actions[a], ".txt", sep="")
      D[a,2] <- stat.function(read.table(filename, head = TRUE)[,col2])
      }
      
     
      mp(D, x=num.actions, t="o", ylab=expression(chi),
          xlab="|A|",leg=c("PI", "PISF"), 
          lwd=1.5, cex=2, cex.lab=1.2, leg.pos = "topleft")
       
  		}




load.MDP.table <- function(delta = 1, 
                     name = "kmeans",
                     size = 100, 
                     p.center = 0.5, 
                     sd = 5, 
                     noise = 0,
                     perc= c(0.1,0.3,0.5),
                     dfs = c(0.9),
                     pos.df = 8, # position of df in the file (i.e., col number)
                     num.actions = 2:10,
                     dir ="./res_tese/dp/",
                     dir.eps = "./fig_tese/",
                     suffix = "") {
   
   ncenters <- round(size * p.center)
   T <- array(0, c(length(num.actions), length(perc), 3))
   
   for (a in 1:length(num.actions)) {
      for (p in 1:length(perc)) {
         filename <- paste(dir, "kmeans_MDP_s", size, "_nc",ncenters,
                                    "_sd", sd, "_n", noise, "_p",
                                    perc[p], "_a", num.actions[a], ".txt",
                                    sep="")
         A <- read.table(filename, head = TRUE)[,pos.df]
         T[a,p,1] <- mean(A)
         T[a,p,2] <- max(A)
         T[a,p,3] <- sd(A)
         }
      }
   T
   }                           
                        
print.MDP.table <- function(T,
                            perc= c(0.1,0.3,0.5),
                            num.actions = 2:10,
                            d = 2) {
   
   T <- round(T, d=d)
   for (a in 1:length(num.actions)) {
      for (p in 1:length(perc)) {
         for (i in 1:3) {
            if (p == length(perc) && i ==3 ) 
                  cat(paste(T[a,p,i], "\\\\",'\n')) 
            else cat(paste(T[a,p,i], " & ")) 
            }
         }
      }
   }
                            
   


                           
plot.MDP.sd1 <- function(name = "kmeans",
                     size = 100, 
                     p.center = 0.5, 
                     sd = 1, 
                     noise = 0,
                     perc= c(0.1,0.3,0.5, 0.7, 0.9),
                     dfs = c(0.5, 0.7, 0.9),
                     num.actions = 10,
                     dir ="./res_tese/dp/",
                     dir.eps = "./fig_tese/") {   
   
   ncenters <- round(size * p.center)
   m <- round(perc * size)
   
   D <- matrix(0, length(perc), length(dfs))
   E <- matrix(0, length(perc), 1)
   
   for (p in 1:length(perc)) {
         filename <- paste(dir, "kmeans_MDP_s", size, "_nc",ncenters,
                          "_sd", sd, "_n", noise, "_p",
                          perc[p], "_a", num.actions, ".txt",
                          sep="")
         T <- read.table(filename, head = TRUE)                          
         D[p,] <- mean(T[,1:length(dfs)])
         E[p,1] <- mean(T[,ncol(T)])
         }
      
     mp(D, x=m, t="o", ylab=expression(chi[m]),
         xlab="m",leg=make.leg.dfs(dfs), 
         lwd=1.5, cex=2, cex.lab=1.2, leg.pos = "topleft")
     print(D)
      
     filename <- paste(dir.eps, "kmeans_MDP_sd1_chi.eps", sep="")
     dev.copy2eps(file = filename)
     
   
     mp(E, x=m, t="o", ylab=expression(xi[m](M, tilde(M))),
         xlab="m", lwd=1.5, cex=2, cex.lab=1.2, leg.pos = "topleft")
      
     filename <- paste(dir.eps, "kmeans_MDP_sd1_xi.eps", sep="")
     dev.copy2eps(file = filename)
    }



plot.MDP.dist <- function(name = "kmeans",
                     size = 100, 
                     p.center = 0.5, 
                     sd = 1, 
                     noise = 0,
                     perc= c(0.1,0.3,0.5, 0.7, 0.9),
                     dfs = 0.9,
                     num.actions = 10,
                     num.avg = 50,
                     dir ="./res_tese/dp/",
                     dir.eps = "./fig_tese/") {   
   
   ncenters <- round(size * p.center)
   m = round(perc * size)
   
   D <- data.frame(matrix(0, num.avg, length(perc)))
   for (p in 1:length(perc)) {
         filename <- paste(dir, "kmeans_MDP_s", size, "_nc",ncenters,
                          "_sd", sd, "_n", noise, "_p",
                          perc[p], "_a", num.actions, ".txt",
                          sep="")
#        b <- length(dfs) + d                                       
         D[,p] <- read.table(filename, head = TRUE)[,3]
         }
      
     boxplot(D, range=0) 
     #filename <- paste(dir.eps, "kmeans_MDP_df",dfs[d],".eps", sep="")
     #dev.copy2eps(file = filename)
    }

                           
                     
                     
                     

print("tese.experiments.dp.R loaded")                 
