source("util.R")


pp <- function(A, tau.max = max(A)) {
# each column is an algorithm; the rows are the problems
	A <- apply(A / apply(A,1,min), 2, sort)
	np <- nrow(A)
	A <- apply(rbind(A,A), 2, sort)
	d <- seq(0, (np-1)/np, by = 1/np)
	d <- sort(c(d,d+1/np))
	plot(1,ylim=c(0,1), xlim=c(1,tau.max), pch=NA_integer_,
		xlab="x", ylab=expression(rho(x)))
	cl <- rainbow(ncol(A))
	for (i in 1:ncol(A)) lines(A[,i],d, lwd=2, col=cl[i])
	legend("bottomright", paste("A", 1:ncol(A),sep=""), col=cl, lwd=2)
	}


pp.st <- function(A, S = matrix(0, nrow(A),ncol(A)), 
			    max.tau = max(A)+2*max(S),
			    length = 1000, lim.inf = 0,
			    stochastic = TRUE, print.legend = TRUE,
             print.cols = 1:ncol(A), ylim = NULL,
             leg.pos = "bottomright",
             labels = FALSE) {
# each column is an algorithm; the rows are the problems
# S is a matrix with the same dimensions as A with the associated standard dvts.
	mins <- apply(A,1,min)
   mins[mins == Inf] <- 1
 	A <- A / mins
	S <- S / mins
	x <- seq(lim.inf,max.tau,l=length)
	
	PP <- matrix(0, length(x), length(print.cols))
	
	for (i in 1:length(x)) {
		for (j in 1:ncol(PP)) {
			for (k in 1:nrow(A)) {
               l <- print.cols[j]
               if (A[k,l] != Inf) {
                    PP[i,j] <- PP[i,j] + pnorm(x[i],A[k,l], S[k,l])
                    }
				}
			PP[i,j] <- PP[i,j] / nrow(A)
        if (round(x[i],2) == 1) print(paste(j,round(PP[i,j],3)))
			}
		}   

     rainbow.cols <- colors.plot(ncol(A))

#      rainbow.cols <- c("RED", "GREEN3","BLUE")

	matplot(x,PP, t="l", col=rainbow.cols[print.cols], lwd=2,
			lty=print.cols, xlab="", ylab="", ylim = ylim)
   if (print.legend) legend(leg.pos, make.leg.rhos(ncol(PP),
                         stochastic), col=rainbow.cols[print.cols], lwd=2,
                         lty=print.cols, bg="WHITE")
   if (labels) {                                        
     labels <- make.leg.rhos(ncol(PP), stochastic)
# Stochastic
        x <- c(3.380543, 8.759058, 3.380543, 8.759058, 8.759058, 8.759058)
        y <- c(0.31, 0.3534203, 0.2081985, 0.31, 0.2187986, 0.158)
# Standard
#        x <- c(3.96, 3.96, 9.35, 9.35, 9.35, 9.35)
#        y <- c(0.3449402,0.2622591, 0.3046597,0.265, 0.221,0.1795781)
     for (i in 1:length(labels)) {
#            pos <- locator(1)
#            print(pos)
          text(x[i],y[i], labels[i])
          }
     }
	}


make.leg.rhos <- function(n, stochastic) {
# makes a legend with labels sigma = sds[1, 2, ...]
   l <- expression()
   if (stochastic) {
	  for(i in 1:n) {
			 l <- c(l, substitute(expression(bar(rho)[ind](x)), list(ind=i))[[2]]) 
			 }
		}
    else {
     for(i in 1:n) {
          l <- c(l, substitute(expression(rho[ind](x)), list(ind=i))[[2]]) 
          }
      }
   l
   }


generate.artificial.examples <- function(
     A = matrix(c(1,1,1,5,3,Inf,5.5,5.5,1,1,2,4,4,6.5,8),5,3)) {

     S <- matrix(0,5,3)
     
     x11(w=5,h=5)
     par(mar=c(2,2,0,0)+0.1)

     # pp base
      pp.st(A,S,10, stochastic = FALSE)
      dev.copy2eps(file = "~/tmp/pp_base.eps")
 
      #pps with the same sd
      S[,] <- 0.1
      pp.st(A,S,10, stochastic = TRUE)
      dev.copy2eps(file = "~/tmp/spp_0_1.eps")
 
      S[,] <- 0.5
      pp.st(A,S,10, stochastic = TRUE)
      dev.copy2eps(file = "~/tmp/spp_0_5.eps")
 
      S[,] <- 1
      pp.st(A,S,10, stochastic = TRUE)
      dev.copy2eps(file = "~/tmp/spp_1.eps")

     #pps in which s1 has a larger sd (equal for all problems)
      S[,]  <- 0.1
      S[,1] <- 5
      pp.st(A,S,10, stochastic = TRUE)
      dev.copy2eps(file = "~/tmp/ppp_a1.eps")
      
      S[,]  <- 0.1
      S[,2] <- 5
      pp.st(A,S,10, stochastic = TRUE)
      dev.copy2eps(file = "~/tmp/ppp_a2.eps")

      S[,]  <- 0.1
      S[,3] <- 5
      pp.st(A,S,10, stochastic = TRUE)
      dev.copy2eps(file = "~/tmp/ppp_a3.eps")

     #pps in which s1 has a larger sd (equal for all problems)
#      S[,]  <- 0.1
#      S[,1] <- 1
#      pp.st(A,S,30, stochastic = TRUE)
#      dev.copy2eps(file = "~/tmp/spp_s1_1.eps")
# 
#      S[,1] <- 5
#      pp.st(A,S,30, stochastic = TRUE)
#      dev.copy2eps(file = "~/tmp/spp_s1_5.eps")
# 
# 
#      S[,1] <- 10
#      pp.st(A,S,30, stochastic = TRUE)
#      dev.copy2eps(file = "~/tmp/spp_s1_10.eps")
# 
#      #pps in which s3 has a larger sd (different across problems)
#      S[,]  <- 0.1
#      S[2,3] <- 3
#      pp.st(A,S,15, stochastic = TRUE)
#      dev.copy2eps(file = "~/tmp/spp_s3_1.eps")
# 
#      S[1:3,3] <- 3
#      pp.st(A,S,15, stochastic = TRUE)
#      dev.copy2eps(file = "~/tmp/spp_s3_5.eps")
# 
#      S[1:3,3] <- 5
#      pp.st(A,S,15, stochastic = TRUE)
#      dev.copy2eps(file = "~/tmp/spp_s3_10.eps")
     

#        dev.off()

     }


generate.gecco.examples <- function(dir = "./files_pp/", length = 1000) {
     z <- function(A) matrix(0, nrow(A), ncol(A))
     
      x11(w=5,h=5)
     par(mar=c(2,2,0,0)+0.1)

     A <- as.matrix(t(read.table(paste(dir,"M-3.dat",sep=""))))
     S <- as.matrix(t(read.table(paste(dir,"S-3.dat",sep=""))))
     o <- c(2,4,5,15,12,10)
     A <- cbind(A[,o], A[,-o])
 
    S <- cbind(S[,o], S[,-o])

         pp.st(A,z(S),10,length=length, print.leg = FALSE, ylim=c(0,0.7),
                   stochastic = FALSE)
          dev.copy2eps(file = "~/tmp/pp_gecco_10.eps")

          pp.st(A,S,10,length=length, print.leg = FALSE,ylim=c(0,0.7))
          dev.copy2eps(file = "~/tmp/ppp_gecco_10.eps")
 
#          pp.st(A,z(S),10,length=length,print.leg = TRUE,
#                    print.cols=1:6, ylim=c(0,0.4),
#                    leg.pos="topleft", stochastic = FALSE, labels = TRUE)
#         dev.copy2eps(file = "~/tmp/pp_gecco_10_red.eps")

#           pp.st(A,S,10,length=length, print.leg =TRUE,
#                     print.cols=1:6, ylim=c(0,0.4),
#                     leg.pos="topleft", labels = TRUE)
#          dev.copy2eps(file = "~/tmp/ppp_gecco_10_red.eps")
     
#      pp.st(A,z(S),5,length=length,print.leg = TRUE,print.cols=c(5,15,10,12),
#           ylim=c(0,0.23))
#      dev.copy2eps(file = "~/tmp/pp_gecco_5_red.eps")
#      pp.st(A,S,5,length=length, print.leg =TRUE,print.cols=c(5,15,10,12),
#           ylim=c(0,0.23))
#      dev.copy2eps(file = "~/tmp/ppp_gecco_5_red.eps")

#       dev.off()

     }
    


print("performance.profiles.R loaded")

