# graphical functions used to generate the thesis' figures
source("ci.95.R")

col.tese <- function(n) {
   col <- c("BLACK", "RED", "BLUE", "GREEN3", "VIOLET")
   col[1:n]
   }

mp <- function(D, x = NULL, leg = NULL, t="l", pch = 21:(20+ncol(D)),
                  lwd=2.5, leg.pos = "topright", text.width = NULL, 
                  bty = "n", ps = 20, inset=0.05, y.intersp = 1,
                  col = col.tese(ncol(D)), bg = col, ...) { 
# alias for 'matplot'
   if (is.null(x)) x <- 1:nrow(D)
   set.par.tese()
   par(ps=ps)
   m <- ncol(D)
   matplot(x, D, lwd = lwd, col = col, bg = bg, t=t, pch = pch, ...)
   if (!is.null(leg)) leg(leg, lty = 1:m, lwd = lwd, col = col,
                           pt.bg=col,pch=pch, x=leg.pos, text.width =
                           text.width, bty = bty, inset = inset, 
                           y.intersp = y.intersp)
   }
   
pl <- function(...) { 
# alias for 'plot'
   set.par.tese()
   plot(...)
   #if (!is.null(legend)) leg(legend)
   }
   
bp <- function(...) {
   set.par.tese()
   boxplot(...)
   #if (!is.null(legend)) leg(legend)
   }
      
   
leg <- function(legend, x="topright", inset = 0.05, lty =1:length(legend),...){
   legend(legend=legend, x= x, inset = inset, bg="WHITE", lty = lty, ...)
   }

set.par.tese <- function() {
# set parameters used throughtout the thesis
   par(
      lty = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"),
      lwd = 1.5,
      mar = c(5, 6, 4, 2) + 0.1,
      ps  = 20)
   }
   
lines.ci <- function(D, col=rgb(0.5,0.5,0.5), plot = TRUE, ...) {
# add the lines corresponding to the mean of the data with the associated
# confidence interval
   means <- apply(D,2,mean)
   sds <- apply(D,2,sd)
   ci <- ci.95(means, sds, ncol(D))
   y2 <- rev(ci[,2])
   if (plot) pl(means, t="l",...)
   else lines(means, ...)
   x <- c(1:ncol(D), ncol(D):1)   
   polygon(x, c(ci[,1],y2), col=col, border=NA)
   }


make.leg.widths <- function(sds) {
# makes a legend with labels sigma = sds[1, 2, ...]
   l <- expression()
   for(s in sds) {
      l <- c(l, substitute(expression(sigma == sig), list(sig=s))[[2]]) 
      }
   l
   }



make.leg.noises <- function(noises) {
# makes a legend with labels eta = noises[1, 2, ...]
   l <- expression()
   for(n in noises) {
      l <- c(l, substitute(expression(eta == noi), list(noi=n))[[2]]) 
      }
   l
   }


make.leg.ncenters <- function(n.centers) {
# makes a legend with labels lambda = n.centers[1, 2, ...]
   l <- expression()
   for(nc in n.centers) {
      l <- c(l, substitute(expression(vartheta == nce), list(nce=nc))[[2]]) 
      }
   l
   }


make.leg.dfs <- function(dfs) {
# makes a legend with labels gamma = dfs[1, 2, ...]
   l <- expression()
   for(d in dfs) {
      l <- c(l, substitute(expression(gamma == d2), list(d2=d))[[2]]) 
      }
   l
   }


make.leg.difs <- function(difs) {
# makes a legend with labels epsilon = dfs[1, 2, ...]
   l <- expression()
   for(d in difs) {
      l <- c(l, substitute(expression(epsilon == d2), list(d2=d))[[2]]) 
      }
   l
   }



lines.int <- function(data, pch = 25, color = "BLACK", add = FALSE, ...) {
## not finished
# plots a sequence of values with an associated interval (which can be a
# confidence interval, standard deviation, etc.)
# ncol(data) must be 3; the first the lower bound of the interval, the second
# the mean and the third the upper bound
  bg <- c("WHITE",color, "WHITE")
  matplot(data, pch=pch, xaxt="n", t="o", lty=c(0,1,0), col = color, bg=bg,
          add=add, ...)
  segments(1:nrow(data), data[,1], 1:nrow(data), data[,3], col = color)
  }

mp.int <- function(data,...) {
## not finished
# matplot with intervals
   set.par.tese()
   ylim <- c(min(data), max(data))
   ncol <- ncol(data) %/% 3
   color <- col.tese(ncol)
   pch <- 21:(21+ncol)
   add <- FALSE
   b <- seq(1,ncol(data),by=3)
   print(b)
   for (i in 1:ncol) {
      print(b[i])
      lines.int(data[,b[i]:(b[i]+2)], color = color[i], pch=pch[i],
               add=add, ylim=ylim, ...)
      add <- TRUE
      }
   }
   
print("tese.plot.R loaded")