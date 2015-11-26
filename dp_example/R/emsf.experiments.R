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


emsf.experiment.plot.results.by.tc <- function(save=FALSE)
{
    tcs <- seq(500, 1000, by = 100)

    for (tc in tcs) {
        emsf.experiment.plot.results(tcs=c(tc), inds=1:4, save=TRUE)
    }
}

emsf.experiment.plot.results.aaai <- function(save=FALSE,
                                              ylim=NULL,
                                              pos=NULL,
                                              cex=NULL,
                                              width=width,
                                              height=height)
{
    tcs <- c(600, 900)
    pch1 <- pch2 <- c(18, 17, 16, 15, 25)

    col1 <- c("darkgreen", "darkred", "darkorange", "darkblue", "darkgreen")
    col2 <- c("green", "red", "black", "blue", "black")

    col1 <- c("seagreen3", "deepskyblue1", "black", "firebrick1", "black")
    col2 <- c("seagreen4", "deepskyblue4", "black", "firebrick4", "black")

    inds <- c(1,2,4)
    emsf.experiment.plot.results(save=TRUE,
                                 tcs=tcs,
                                 alphas=c(0.1, 0.3, 0.5, 0.7),
                                 inds=inds,
                                 pch=c(pch1[inds], pch2[inds]),
                                 col=c(col1[inds], col2[inds]),
                                 ylim=ylim,
                                 pos=pos,
                                 cex=cex,
                                 width=width,
                                 height=height)

    ## comando utilizado para plotar:
    ## source("emsf.experiments.R"); emsf.experiment.plot.results.aaai(save=TRUE, ylim=NULL, cex=1.5, width=8.75, height=7)
}


emsf.experiment.plot.results <- function(n = 100, m = 10,
                                         tcs = seq(100, 1000, by = 100),
                                         max.it = 25*10^3,
                                         alpha.dir = 0.5,
                                         alphas = c(0.1, 0.3, 0.5, 0.7, 1),
                                         inds = 1:5,
                                         num.points = 25, dir = "./files/", idx = "", save=FALSE,
                                         pch=NULL,
                                         col=NULL,
                                         ylim=NULL,
                                         pos=NULL,
                                         cex=NULL,
                                         width=width,
                                         height=height,
                                         ...)
{
    R <- NULL
    D <- NULL
    for (i in 1:length(tcs))
    {
        prefix <- ps(dir, paste("emsf", n, m, tcs[i], max.it, alpha.dir, sep="_"))
        path <- ps(prefix, idx, "_kld.txt")
        print(path)
        TMP <- as.matrix(read.table(path, header = TRUE))
        R <- cbind(R, TMP[seq(1, nrow(TMP), length=num.points),inds])

        prefix <- ps(dir, paste("emsf", n, m, tcs[i], max.it, alpha.dir, sep="_"))
        TMP <- as.matrix(read.table(ps(prefix, idx, "_ser.txt"), header = TRUE))
        D <- cbind(D, TMP[seq(1, nrow(TMP), length=num.points),inds])
    }

    par(cex=1.65, mai=c(1.3, 1.4, 0.1, 0.25)) # Margens em polegadas (down, left, top, right)
    cex <- 1
    mp(seq(1, max.it, l = num.points), R, R + D , R - D, xlab = expression(tau), ylab = expression("KL"[rho]*"(P, DK)", main=paste(tcs)), ylim=ylim, pch=pch, col=col, cex=cex, show.shadow=FALSE,
       cex.lab=cex,
       cex.axis=cex,
       lty=c(1,1,1,1,1,1,1,1,1,1),
       lwd=3)

    l <- make.leg.tcs.alphas(tcs,alphas[inds])
    leg(pos="topright",l[1:3], pch=pch[1:3], col=col[1:3], cex=cex, border=NULL, box.lwd=0, bty="n", lty=c(1,1,1,1,1,1,1,1,1,1), lwd=3)
    leg(pos="bottomleft",l[4:6], pch=pch[4:6], col=col[4:6], cex=cex, border=NULL, box.lwd=0, bty="n", lty=c(1,1,1,1,1,1,1,1,1,1), lwd=3)

    grid(lwd=2)
    if (save) {
        dev.copy2pdf(file=paste(sep="", "~/online_em_sf/fig/emsf_tc_alpha.pdf"), width=width, height=height)
    }
}
