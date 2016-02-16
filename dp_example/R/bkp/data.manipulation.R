# These functions were first develop for the LSPI algorithm, but
# later I realized the KBRL algorithm uses a very similar structure.
# So I moved these functions here and both modules (LSPI and KBRL)
# now share the same datafile structure, which is a list (s,a,r,s2)

source("util.R") 


gp <- function(x, y) {
#generate two-dimensional points from two intervals (like a grid)
	pts <- matrix(0, length(x)*length(y), 2)
	k <- 1
	for (ix in x) {
		for (iy in y) {
			pts[k,1] <- ix
			pts[k,2] <- iy
			k <- k +1
			}
		}
	pts
	}

normalize <- function(X, means=apply(X,2,mean),stdevs=apply(X,2,sd),
                      return.params=FALSE) {
   stdevs[stdevs==0] <- 1
  N <- nrow(X)
  p <- ncol(X)
  X <- (X - matrix(rep(means,N),N,p,byrow=TRUE))/
  matrix(rep(stdevs,N),N,p,byrow=TRUE)
  # another way:
  #min.x <- apply(X, 2,min)
  #max.x <- apply(X, 2,max)
  #min.x <- matrix(rep(min.x, nrow(X)), nrow(X), ncol(X), byrow = TRUE)
  #max.x <- matrix(rep(max.x, nrow(X)), nrow(X), ncol(X), byrow = TRUE)
  #X <- (X - min.x) / (max.x - min.x)
  if (return.params) list(data=X,means=means,stdevs=stdevs)
  else X
  }


unnormalize <- function(X, means, stdevs) {
  N <- nrow(X)
  p <- ncol(X)
  X <- X * matrix(rep(stdevs,N),N,p,byrow=TRUE) + 
		matrix(rep(means,N),N,p,byrow=TRUE)
  X
  }	

  
make.empty.sars <- function() {
	list(s = NULL, a = NULL, r = NULL, s2 = NULL, g= NULL, 
                  means = NULL, stdevs = NULL)
	}
 
collect.transitions <- function(sa, transition.function, A, 
                        normalize = TRUE, ...) {
# sa is a list l = (s, a)
# returns a list with entries s, a, r, s2
   illegals <- NULL
   
   s2 <- matrix(0, nrow(sa$s), ncol(sa$s))
   r <- array(0, nrow(sa$s))
   g <- array(FALSE, nrow(sa$r))
	
   for (i in 1:length(sa$a)) {
      t <- transition.function(sa$s[i,], A[sa$a[i]], ...)
      
      if (is.null(t$s)) illegals <- c(illegals, i)
      else {
         s2[i, ] <- t$s
         r[i] <- t$r
         g[i] <- t$g
         }
      }
	
   if (!is.null(illegals)) {
		#remove illegal moves
      sa$s <- sa$s[-illegals, ] 
      sa$a <- sa$a[-illegals] 
      r	 <- r[-illegals] 
      s2 <- s2[-illegals, ] 
      }
	
   sars <- list(s = sa$s, a = sa$a, r = r, s2 = s2, g = g, means = 0, stdevs=1)
   if (normalize) { ## is this the best way to normalize?
          S <- rbind(sars$s, sars$s2)
          tmp <- normalize(S, return.params = TRUE)  
          sars$s <- normalize(sars$s, means=tmp$means, stdevs= tmp$stdevs)
          sars$s2 <- normalize(sars$s2,means=tmp$means,stdevs= tmp$stdevs)
          sars <- list(s = sars$s, a = sars$a, r = sars$r, s2 =
          sars$s2, g = sars$g, means = tmp$means, stdevs = tmp$stdevs)
          }

   sars
   }



collect.episodes <- function(transition.function, num.samples, steps, A, min.s,
                      max.s, goal.only = TRUE, normalize = TRUE, 
                     filename = NULL, ...) {
     samples <- list(s = NULL, a = NULL, r = NULL, s2 = NULL, g = NULL)
     s <- array(0, length(min.s))
     while (length(samples$a) < num.samples) {
          for (i in 1:length(min.s)) s[i] <- runif(1, min.s[i], max.s[i])
          a <- runif.int(1,1,length(A))
          goal <- FALSE
          stp <- 0
          sas <- list(s = NULL, a = NULL, r = NULL, s2 = NULL, g = NULL)
          while (!goal && stp < steps) {
               stp <- stp + 1
               t <- transition.function(s, A[a], ...)
               sas$s <- rbind(sas$s, s)
               sas$s2 <- rbind(sas$s2, t$s)
               sas$r <- c(sas$r, t$r)
               sas$a <- c(sas$a, a)
               sas$g <- c(sas$g, t$g)
               
               goal <- t$g
               s <- t$s
               a <- runif.int(1,1,length(A))               
               }
          if (!(goal.only && !goal)) {
               remaining <- num.samples- length(samples$a)
               b <- 1
               if (length(sas$a) > remaining) b <- length(sas$a) - remaining + 1
               e <- length(sas$a)
               samples$s <- rbind(samples$s, sas$s[b:e,])
               samples$s2 <- rbind(samples$s2, sas$s2[b:e,])
               samples$r <- c(samples$r, sas$r[b:e])
               samples$a <- c(samples$a, sas$a[b:e])
               samples$g <- c(samples$g, sas$g[b:e])            
               }
               
          }
     if (normalize) {
          S <- rbind(samples$s, samples$s2)
          tmp <- normalize(S, return.params = TRUE)  
          samples$s <- normalize(samples$s, means=tmp$means, stdevs= tmp$stdevs)
          samples$s2 <- normalize(samples$s2,means=tmp$means,stdevs= tmp$stdevs)
          samples <- list(s = samples$s, a = samples$a, r = samples$r, s2 =
             samples$s2, g = samples$g, means = tmp$means, stdevs = tmp$stdevs)
          }
   if (!is.null(filename)) save.sars(samples, filename)
   samples
   }


normalize.sars <- function(sars, ...) {
	tmp <- normalize(sars$s, return.params = TRUE, ...)  
	sars$s <- tmp$data
	sars$s2 <- normalize(sars$s2, means = tmp$means, stdevs = tmp$stdevs)
	sars <- list(s = sars$s, a = sars$a, r = sars$r, s2 = sars$s2, 
					 g = sars$g, means = tmp$means, stdevs = tmp$stdevs)
	sars
	}

unnormalize.sars <- function(sars) {
	sars$s <- unnormalize(sars$s, sars$means, sars$stdevs)
	sars$s2 <- unnormalize(sars$s2, sars$means, sars$stdevs)	
	sars$means[ ] <- 0
	sars$stdevs[ ] <- 1
	sars
	}


make.sa <- function(S, num.actions) {
	sa <- list(s = NULL, a  = NULL)
	for (a in 1:num.actions) {
		sa$s <- rbind(sa$s,S)
		sa$a <- c(sa$a, rep(a, nrow(S)))
		}
	sa
	}


make.lsars <- function(sars) {
	lsars <- NULL
	num.actions <- max(sars$a)
	for (a in 1:num.actions) {
		# select the data referring to the current action
		ind <- sars$a == a
		s <- sars$s[ind,]
		r <- sars$r[ind]
		s2 <- sars$s2[ind,]
		g <- sars$g[ind]
		
		lsars <- c(lsars, list(list(s=s, r = r, s2 = s2, g = g)))
		}
	lsars <- c(lsars, list(means = sars$means, stdevs = sars$stdevs,
               num.actions = num.actions))
	lsars
	}


make.sars <- function(lsars) {
   sars <- list(s = NULL, a = NULL, r = NULL, s2 = NULL, g = NULL)
   for (a in 1:lsars$num.actions) {
      # select the data referring to the current action
      sars$s <- rbind(sars$s, lsars[[a]]$s)
      sars$a <- c(sars$a, rep(a, length(lsars[[a]]$r)))
      sars$r <- c(sars$r, lsars[[a]]$r)
      sars$s2 <- rbind(sars$s2, lsars[[a]]$s2)
      sars$g <- c(sars$g, lsars[[a]]$g)
      }
  c(sars, list(means = lsars$means, stdevs = lsars$stdevs))
   }

	
save.sars <- function(sars, filename) {
	wt(sars$s,paste(filename,".s",sep=""))
	wt(sars$a,paste(filename,".a",sep=""))
	wt(sars$r,paste(filename,".r",sep=""))
	wt(sars$s2,paste(filename,".s2",sep=""))
	if (!is.null(sars$g))      wt(sars$g,paste(filename,".g",sep=""))
	if (!is.null(sars$means))  wt(sars$means,paste(filename,".means",sep=""))
	if (!is.null(sars$stdevs)) wt(sars$stdevs,paste(filename,".stdevs",sep=""))
	}
	
load.sars <- function(filename) {	
	sars <- list(s = NULL, a = NULL, r = NULL, s2 = NULL, g = NULL, means =
					 NULL, stdevs = NULL)
	sars$s <- as.matrix(read.table(paste(filename,".s",sep="")))
	sars$a <- drop(as.matrix(read.table(paste(filename,".a",sep=""))))
	sars$r <- drop(as.matrix(read.table(paste(filename,".r",sep=""))))
	sars$s2 <- as.matrix(read.table(paste(filename,".s2",sep="")))
	
	if (file.exists(paste(filename,".g",sep=""))) {
		sars$g <- drop(as.matrix(read.table(paste(filename,".g",sep=""))))
		}
	if (file.exists(paste(filename,".means",sep=""))) {
		sars$means <- drop(as.matrix(read.table(paste(filename,".means",sep=""))))
		}
	if (file.exists(paste(filename,".stdevs",sep=""))) {
		sars$stdevs<-drop(as.matrix(read.table(paste(filename,".stdevs",sep=""))))
		}
	sars
	}

			
merge.sars <- function(sars1, sars2, unnormalize = c(FALSE,FALSE),
								normalize=FALSE) {
	if (unnormalize[1]) sars1 <- unnormalize.sars(sars1)
	if (unnormalize[2]) sars2 <- unnormalize.sars(sars2)

	sars1$s <- rbind(sars1$s, sars2$s)
	sars1$a <- c(sars1$a, sars2$a)
	sars1$r <- c(sars1$r, sars2$r)
	sars1$g <- c(sars1$g, sars2$g)
	sars1$s2 <- rbind(sars1$s2, sars2$s2)
	
	sars1$means <- 0
	sars1$stdevs <- 1

	if (normalize) sars1 <- normalize.sars(sars1)

	sars1
	}

print("data.manipulation.R loaded")