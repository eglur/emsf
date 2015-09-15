library(gtools)
source("emsf.R")
source("util.R")
source("data.plot.R")


emsf.experiment.generate.P <- function(n, alpha.dir, na = 1)
{
    P <- array(0, c(n, n, na))
    for (a in 1:na) P[,,a] <- rdirichlet(n,rep(alpha.dir, n))
    P 
}

emsf.experiment.fixed.tc <- function(n, alpha.dir, m, tc, alphas, max.it, na = 1, num.avg = 50, dir = "./files/", idx = "")
{
    R <- NULL
    D <- NULL
    
    for (i in 1:length(alphas))
    {
        L <- NULL
        
        for (j in 1:num.avg)
        {
            print(paste("Running alpha = ", alphas[i], "  run", j))
            P <- emsf.experiment.generate.P(n, alpha.dir, na)
            L <- cbind(L, emsf(P, m, tc, tc, alphas[i], max.it))
                                        #       print(L)
        }
        
        R <- cbind(R, apply(L, 1, mean))
        colnames(R)[ncol(R)] <- paste(alphas[i])
        
        D <- cbind(D, apply(L, 1, sd) / sqrt(num.avg))
        colnames(D)[ncol(D)] <- paste(alphas[i])
        
        prefix <- ps(dir, paste("emsf", n, m, tc, max.it, alpha.dir, sep="_"))
        write.table(R, ps(prefix, idx, "_kld.txt"), row.names=FALSE, col.names=TRUE, quote=FALSE)
        write.table(D, ps(prefix, idx, "_ser.txt"), row.names=FALSE, col.names=TRUE, quote=FALSE)
    }
    
    list(R = R, D = D)
}


make.leg.tcs.alphas <- function(tcs, alphas) {
                                        # makes a legend with labels epsilon = dfs[1, 2, ...]
    l <- expression()
    for (tc in tcs)
    {
        for(alpha in alphas) {
            l <- c(l, substitute(expression(t[c] == TC~alpha == ALPHA), list(TC=tc, ALPHA=alpha))[[2]]) 
        }
    }
                                        #    print(l)
    l
}


emsf.experiment.plot.results <- function(n = 100, m = 10, tcs = seq(100, 500, by = 100), max.it = 25*10^3, 
                                         alpha.dir = 0.5, alphas = c(0.1, 0.3, 0.5, 0.7, 1), inds = 1:5, num.points = 30, dir = "./files/", idx = "")
{
    
    R <- NULL
    D <- NULL
    for (i in 1:length(tcs))
    {
        prefix <- ps(dir, paste("emsf", n, m, tcs[i], max.it, alpha.dir, sep="_"))
        TMP <- as.matrix(read.table(ps(prefix, idx, "_kld.txt"), header = TRUE))
        R <- cbind(R, TMP[seq(1, nrow(TMP), length=num.points),inds])
        
        prefix <- ps(dir, paste("emsf", n, m, tcs[i], max.it, alpha.dir, sep="_"))
        TMP <- as.matrix(read.table(ps(prefix, idx, "_ser.txt"), header = TRUE))
        D <- cbind(D, TMP[seq(1, nrow(TMP), length=num.points),inds])
    }
    
    print(R)
    mp(seq(1, max.it, l = num.points), R, R + D , R - D, xlab = expression(tau), ylab = expression("KL"[rho]*"(P, DK)"))
    
    l <- make.leg.tcs.alphas(tcs,alphas[inds])
    leg("topright",l)
}
