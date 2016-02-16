source("util.R")
library(plotrix)

set.par <- function() 
{
                                        # set parameters used  in the plots
    par(
        lwd = 2,
        mar = c(5, 6, 3, 2) -1 + 0.1,
        ps  = 20)
}

get.lty <- function(n) 
{
    types <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
    types[(1:n - 1) %% length(types) + 1]
}

                                        # original
get.col <- function(n, alpha = 1) 
{
    col <- c(rgb(0,0,0, alpha), rgb(1,0,0, alpha), rgb(0,0,1, alpha), rgb(0,1,0,alpha), rgb(0,1,1,alpha), 
             rgb(1,0,1,alpha), rgb(1,1,0,alpha), rgb(0.5,0.5,0.5,alpha))
    col[(1:n - 1) %% length(col) + 1]
}

                                        # ## cr paper
                                        # get.col <- function(n, alpha = 1) 
                                        # {
                                        #    col <- c(rgb(0,0,0, alpha), rgb(1,0,0, alpha), rgb(0,0,1, alpha), rgb(0,0.5,0,alpha), rgb(0,0.3,0,alpha), 
                                        #       rgb(0,0.1,0,alpha), rgb(1,1,0,alpha), rgb(0.5,0.5,0.5,alpha))
                                        #    col[(1:n - 1) %% length(col) + 1]
                                        # }


get.pch <- function(n) 
{
    if (n <= 5) seq(21, 21 + n -1, by = 1)
    else c(21:25, 1:(n-6)) 
}

mp <- function(x, Y, 
               Y.ui = NULL, 
               Y.li = NULL,
               plot.error = TRUE, 
               show.error = NULL,
               show.shadow = NULL, 
               alpha = 0.1,
               transparency = TRUE, # transparency doesn't work with EPS files
               num.avg = 50, 
               ylim = NULL,
               pch = NULL,
               t = "o",
               lty = NULL,
               col = NULL,
               inds = NULL, # this is for plotting only a subset of the columns of Y
               ...)
{
    p <- ncol(Y)

    if (is.null(inds)) inds <- 1:ncol(Y)

    Y <- matrix(Y[,inds], nrow(Y), length(inds))
    Y.ui <- matrix(Y.ui[,inds], nrow(Y), length(inds))
    Y.li <- matrix(Y.li[,inds], nrow(Y), length(inds))

    
    if (is.null(col)) col <- get.col(p)[inds]
    if (is.null(lty)) lty <- get.lty(p)[inds]
    if (is.null(pch)) pch <- get.pch(p)[inds]
    
    if (is.null(ylim))
    {
        ylim <- c(min(Y), max(Y))
        if (!is.null(Y.ui) && !is.null(Y.li) && plot.error) ylim <- c(min(Y.li), max(Y.ui))
    }
    
    matplot(x, Y, t=t, bg = col, col= col, lty= lty, pch= pch, ylim = ylim, ...)
    
    
                                        # Error bars
    if (!is.null(Y.ui) && !is.null(Y.li) && plot.error) {
        
        if (is.null(show.error))  show.error <- rep(TRUE, p)
       else if (length(show.error)==1) show.error <- rep(show.error, p)
        if (is.null(show.shadow)) show.shadow <- rep(TRUE, p)
       else if (length(show.shadow)==1) show.shadow <- rep(show.shadow, p)

        for (i in 1:ncol(Y)) {
            if (show.error[i])
            {
                
                plotCI(x, Y[,i], ui=Y.ui[,i], li = Y.li[,i], add=TRUE,
                       col=col[i], scol= col[i], lty=lty[i], pch=pch[i], 
                       pt.bg=col[i], sfrac=0)
            }
            
            if (show.shadow[i]) {
                px <- c(x,rev(x))
                py <- c(Y.li[,i], rev(Y.ui[,i]))
                
                if (transparency) polygon(px,py, col = rgb(0,0,0, alpha), border = NA)
                else
                {
                    polygon(px,py, col = rgb(1 - alpha, 1 - alpha, 1 - alpha), border = NA)
                    lines(x, Y[,i], t=t, col= col[i], pch= pch[i])
                }

            }
            
        }
    }
    
}



leg <- function(pos = "topleft", leg, y.intersp = 1.5, 
                pch = NULL, 
                col = NULL,
                lty = NULL, 
                pt.bg = NULL, 
                ...)
{
    p <- length(leg)
    if (is.null(pch)) pch <- get.pch(p)
    if (is.null(col)) col <- get.col(p)
    if (is.null(lty)) lty <- get.lty(p)
    if (is.null(pt.bg)) pt.bg <- get.col(p)

    legend(pos, leg, pt.bg = pt.bg, col=col, lty= lty, y.intersp  = y.intersp, pch = pch, ...)
}



make.leg.epsilon <- function(eps) {
                                        # makes a legend with labels epsilon = dfs[1, 2, ...]
    l <- expression()
    for(d in eps) {
        l <- c(l, substitute(expression(epsilon == d2), list(d2=d))[[2]]) 
    }
    l
}


make.leg.iota<- function(iot) {
                                        # makes a legend with labels epsilon = dfs[1, 2, ...]
    l <- expression()
    for(d in iot) {
        l <- c(l, substitute(expression(iota == d2), list(d2=d))[[2]]) 
    }
    l
}

print("data.plot.R loaded")
