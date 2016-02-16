
source("pole.double.R")
source("data.manipulation.R")
source("dp.R")
source("arch.kernel.R")
source("kbrl.R")
source("lspi.rbfn.optimized.R")
source("tese.plot.R")

pole.data <- function (num.episodes, 
							  normalize = FALSE, 
							  lsars = FALSE, 
                       double = FALSE, 
                       p = 0.75,
                       p2 = 0,
                       max.pos        = 2.4, 
                       max.vel        = 2.4,
                       max.angle      = (pi/5), 
                       max.angle.vel  = (pi/5), 
                       max.angle2      = (pi/5), 
                       max.angle.vel2  = (pi/5),
                       A = c(-10,0,10), 
                       verbose = FALSE, ...){
   
   max.pos        <-  p * max.pos 
   max.vel        <-  p * max.vel
   max.angle      <-  p * max.angle 
   max.angle.vel  <-  p * max.angle.vel
   max.angle2     <-  p2* max.angle2
   max.angle.vel2 <-  p2* max.angle.vel2
	pos <- runif(num.episodes, -max.pos, max.pos)
	vel <- runif(num.episodes, -max.vel, max.vel)
	angle1     <- runif(num.episodes, -max.angle, max.angle)
	angle1.vel <- runif(num.episodes, -max.angle.vel, max.angle.vel)
	angle2     <- NULL
	angle2.vel <- NULL
	if (double) {
		angle2     <- runif(num.episodes, -max.angle2, max.angle2)
		angle2.vel <- runif(num.episodes, -max.angle.vel2, max.angle.vel2)
		}
		
	s.initial <- cbind(pos,vel,angle1,angle1.vel)
	if(double) s.initial <- cbind(s.initial, angle2, angle2.vel)

	sars <- collect.episodes(s.initial, pole.double.transition, A,
									max.steps = 10^3, normalize = normalize, 
									verbose = verbose, double = double, ...)
		
	if (lsars) make.lsars(sars)
	else sars
	}


pole.make.ST <- function(n.grid,
                       p = 0.75,
							  max.pos = c(2.4, 2.4, (pi/5), (pi/5)),
                       min.pos = -max.pos,
                       d = 4) {
   max.pos <- p*max.pos
   min.pos <- p*min.pos
   grid <- matrix(0, n.grid^d, d)
 	uniform.grid(grid, max.pos, min.pos)
	}
	
	
info.sars <- function(sars) {
	goals <- sars$g
	ep <- (1:length(goals))[sars$g]
	
	l <- array(0, length(ep))
	last <- 0
	for (i in 1:length(ep)) {
		l[i] <- ep[i] - last
		last <- ep[i]
		}
		
	print(paste("Number of transitions: ", length(sars$g)))
	print(paste("Number os episodes: ", length(l)))
	print(paste("Mean length: ", mean(l)))
	print(paste("Max. length: ", max(l)))
	print(paste("Min. length: ", min(l)))
	print(paste("Sd. length: ", sd(l)))
	}
	

	
convert.index <- function(inds, n) {
	d <- ncol(inds)
	id <- inds[,d]
	block <- 1
	for (j in (d-1):1) {
		block <- block * n
		id <- id + (inds[,j] - 1) * block
		}
	id
	}

 
pole.build.mdp.grid <- function(n.grid, sars, d=4, num.actions=3) {
# sars must not be normalized! (in order to be possible to test the policy)	
 	s.max <- apply(abs(sars$s), 2, max)
 	s.max <- apply(rbind(s.max, apply(abs(sars$s2), 2, max)),2,max)
 	
 	s.delta <- array(0, d)
	for (i in 1:d) {
		seq <- seq(-s.max[i], s.max[i], l=n.grid) ## not elegant...
		s.delta[i] <- seq[2]-seq[1]
		}
	
	n <- n.grid^d
 	P <- array(0, c(n,n,num.actions))
 	R <- array(0, c(n,n,num.actions))
	 	
	ind1 <- matrix(0,nrow(sars$s),d)
	ind2 <- matrix(0,nrow(sars$s),d)
	# calculate the index of states
	for (j in 1:d) {
		ind1[,j] <- round((sars$s[,j]  + s.max[j]) / s.delta[j]) + 1
		ind2[,j] <- round((sars$s2[,j] + s.max[j]) / s.delta[j]) + 1
		}
		
	l <- convert.index(ind1, n.grid)
	c <- convert.index(ind2, n.grid)
	
	N <- matrix(0, n, num.actions) # stores the number of transitions 
	for (i in 1:length(l)) {
		P[l[i],c[i],sars$a[i]] <- 	P[l[i],c[i],sars$a[i]] + 1
		R[l[i],c[i],sars$a[i]] <- 	R[l[i],c[i],sars$a[i]] + sars$r[i]
      N[l[i],sars$a[i]] <- N[l[i],sars$a[i]] + 1
		}
		
	r <- matrix(0, n, num.actions)
	for (a in 1:num.actions) {
		den <- apply(P[,,a], 1, sum)
		non.zero <- den != 0
		P[non.zero,,a] <- P[non.zero,,a] /den[non.zero] 
		for (i in 1:n) r[i,a] <- P[i,,a] %*% R[i,,a]
		}
	list(P = P, R = r, N = N, n.grid = n.grid, s.max = s.max, s.delta = s.delta)
	}
	

pole.test.policy.grid <- function(pi, s, s.max, s.delta, n.grid,
								max.steps=180000,
								A = c(-10,0,10), 
								double = FALSE, ...) {
	g <- FALSE
	it <- 0
	while (!g && it < max.steps) {
		ind <- round((s + s.max) / s.delta) + 1
		ind <- pmax(pmin(ind, length(pi)),1)
		id <- convert.index(matrix(ind,1,length(ind)), n.grid)
		a <- pi[id]
		t <- pole.double.transition(s,A[a], double = double)
		g <- t$g
		s <- t$s
		
		it <- it + 1
		}
	it
	}
	

find.mdp.nonzero <- function(R, P, S, 
										tau = 1) {
   
   num.actions <- ncol(R)
   num.states <- nrow(R)
   
 	non.zero <- array(FALSE, num.states)
 	for (a in 1:num.actions) {
 		non.zero <- non.zero | (apply(P[,,a],1 ,sum) != 0)
 		}
 	print(sum(!non.zero) / nrow(P) * 100)
 	
	m <- sum(non.zero)	
	Pb <- array(0, c(m,m,num.actions))
	Rb <- matrix(0,m, num.actions)
	
	SN <- normalize(S)
	rbfn <-  make.rbfn.kbrl(SN[non.zero,], tau) 
	D <- rbfn.norm.design.matrix(rbfn, SN)
			
	for (a in 1:num.actions) {
		K <- P[non.zero,,a]
		Pb[,,a] <- K %*% D
		Rb[,a] <- R[non.zero,a]
		}
   list (P = Pb, R = Rb, D = D)	
	}
	
	
pole.grid <- function(n.grid, ST,
								sars = NULL, 
								exclude.zeros = FALSE,
								tau = 1,
								num.episodes = 10^4, 
							 	double =FALSE,
							 	df = 0.99,
							 	max.iter.pi = 15,
							 	verbose = TRUE,
							 	max.steps.test = 3000,
							 	M = NULL,
							 	pi = NULL,
							 	...){
							 
	if (is.null(sars)) {
		if (verbose) print("Generating data...")
		sars <- pole.data(num.episodes, normalize = FALSE, lsars = FALSE,
								double = double, verbose = verbose, ...)
		}
	
	if (is.null(M)) {
		if (verbose) print("Building the MDP...")
		M <-  pole.build.mdp.grid(n.grid,sars)  	
		}
	
	D <- NULL
	if (exclude.zeros) {
		if (verbose) print("Reducing the MDP...")
		S <- uniform.grid(matrix(0,nrow(M$P),ncol(ST)),M$s.max,-M$s.max)
		MZ <- find.mdp.nonzero(R = M$R, P = M$P, S = S, tau = tau)
		M$R <- MZ$R
		M$P <- MZ$P
		D <- MZ$D
		if (verbose) print(paste("New size of the MDP:", nrow(M$P)))
		}
		
		
	if (is.null(pi)) {
		if (verbose) print("Running policy iterarion...")
		PIR <- policy.iteration(M$R, M$P, df=df, max.iter=max.iter.pi)
		pi <- PIR$pi
		if (exclude.zeros) {
			Q <- MZ$D %*% PIR$Q
			pi <- apply(Q,1,which.max)
			}
		}
	
	if (verbose) print("Testing the policy...")
	
	ns <- array(0, nrow(ST))
	
	for (i in 1:nrow(ST)) {
		ns[i] <- pole.test.policy.grid(pi, ST[i,], M$s.max, M$s.delta,
													M$n.grid, max.steps = max.steps.test,...)
		if (verbose) print(paste("s",i,ns[i]))
		}
	ns
	}


pole.experiment.grid <- function(n.grid=c(3,4,5,6,7), 
											sars = NULL,
											ST = pole.make.ST(3),
				 							filename = "pole.grid",
											dir = "./res_tese/pole/",
			  								verbose = TRUE,
											...){
	for (i in 1:length(n.grid)) {
		if (verbose) print(paste("n =", n.grid[i],"(", i,"/",length(n.grid),")"))
		ns <- matrix(0, 1, nrow(ST))
		ns[1,]<- pole.grid(n.grid[i],ST, sars = sars, exclude.zeros =FALSE,
								 verbose = verbose, ...)
      if (verbose) print(paste("n =", n.grid[i], "mean =", mean(ns)))
      wt(ns, paste(dir,filename,"_n",n.grid[i],".txt",sep=""))
      }
	if (verbose) print("All done!")
	}
													

	
select.pole.kmeans <- function(S, k, metric=NULL, iter.max=20,
								s.max, s.delta, n.grid) {
 	d <- ncol(S)

 	D <- normalize(S, return.params = TRUE)
 	
 	S <- D$data
 	C <- kmeans(S, k, iter.max = iter.max)$centers
	C <- unnormalize(C, means=D$means, stdevs = D$stdevs)
	
	# calculate the index of states
	ind <- matrix(0,nrow(C),d)
	for (j in 1:d) {
		ind[,j] <- round((C[,j]  + s.max[j]) / s.delta[j]) + 1
		}
	convert.index(ind, n.grid)
	}
	

select.pole.grid <- function(S, k, metric=NULL, iter.max=10,
								s.max, s.delta, n.grid) {
 	d <- ncol(S)
 	grid <- matrix(0,k,d)
 	C <- uniform.grid(grid, s.max, -s.max)
 	
 	# calculate the index of states
	ind <- matrix(0,nrow(C),d)
	for (j in 1:d) {
		ind[,j] <- round((C[,j]  + s.max[j]) / s.delta[j]) + 1
		}
	 convert.index(ind, n.grid)
	}
	 	
	
 select.pole.clara <- function(S, k, metric, stand = TRUE,
										s.max, s.delta, n.grid) { 
# needs library("cluster")
    C <- S[clara(S, k, metric = metric, stand = stand, medoids.x =
				FALSE)$i.med,]
 	# calculate the index of states
	ind <- matrix(0,nrow(C),ncol(C))
	for (j in 1:ncol(C)) {
		ind[,j] <- round((C[,j]  + s.max[j]) / s.delta[j]) + 1
		}
	 convert.index(ind, n.grid)
	}
	 
	 
select.pole.pam <- function(S, k, metric, stand = TRUE,
										s.max, s.delta, n.grid) {
# needs library("cluster")
    C <- S[pam(S, k, metric = metric, stand = stand)$id.med,]
 	# calculate the index of states
	ind <- matrix(0,nrow(C),ncol(C))
	for (j in 1:ncol(C)) {
		ind[,j] <- round((C[,j]  + s.max[j]) / s.delta[j]) + 1
		}
	 convert.index(ind, n.grid)
	}


arch.kernel.policy.iteration.alt <- function(R, P, df, m, S, N, 
                                        keep.all = FALSE,
													 tau = 1,
                                        pi = NULL, 
                                        max.iter = Inf,
                                        solve.mp.function = solve.mp,
                                        selection.function= select.pole.clara,
                                        metric = "euclidean", 
                                        verbose = FALSE, 
                                        ...) {
# Policy iteration using archetype algorithm with kernels 
# This version excludes the lines P[i,] == 0
# R is a |S| x |A| matrix with the rewards associated with the (s,a) transitions
# P is a |S| x |S| x |A| matrix with the transition probabilities P(s_i| s_j, a)
# df is the discount factor
# m is the number of archetypes
# S is the "observational" representation of states
# pi is the initial policy
# max.iter is the maximum number of iterations 
# solve.mp.function is the function used to compute the value of a policy
#    (this can be used to implement the "modified policy iteration algorihtm",
#     which is the method recommended by Puterman (p.186))
# method is the dissimilarity function used (see function "dist")
   
   num.actions <- ncol(R)
   num.states <- nrow(R)
   
   if (verbose) print("Selecting archetypes...")
   k <- m %/% num.actions
   m <- k * num.actions
   if (keep.all) {
      k <- -Inf
      for (a in 1:num.actions) k <- max(sum(N[,a] != 0), k)
      m <- k * num.actions
      }
   if (verbose) print(paste("Using", m, "archetypes"))
   
# 	ind <- matrix(0, k, num.actions)
# 	for (a in 1:num.actions) {
# 		ind[,a] <- selection.function(S[apply(P[,,a],1,sum)!=0,], k, metric,...)
# 		}

  ind <- matrix(0, k, num.actions)
  for (a in 1:num.actions) {
     ind[,a] <- order(N[,a], decreasing=TRUE)[1:k]
     }

   K <- matrix(0, m, num.states)
   rb <- array(0, m)
   for (a in 1:num.actions) {
      b <- (a-1) * k + 1
      e <- b + k - 1
      K[b:e,] <- P[ind[,a],,a]   
      rb[b:e] <- R[ind[,a],a]
      }
  
   # determine matrix D
   if (verbose) print("Allocating states...")
   SN <- normalize(S)
	Da <- array(0, c(num.states, k, num.actions))
	for (a in 1:num.actions) {
		rbfn <-  make.rbfn.kbrl(SN[ind[,a],], tau)
		Da[,,a] <- rbfn.norm.design.matrix(rbfn,SN)
		}
	
	if (is.null(pi)) pi <- sample(1:num.actions, num.states, TRUE)

   pi.old <- array(0, num.states)
   Q <- matrix(0, num.states, num.actions)
      
   if (verbose) print("Running policy iteration. This may take a while...")
   if (verbose) print("Iteration: ")
   it <- 0
   while (it < max.iter && !identical(pi, pi.old)) {
      if (verbose) cat(paste(it+1," "))
      # policy evaluation
      D <- matrix(0, num.states, m)
      for (i in 1:nrow(D)) {
         b <- (pi[i] - 1) * k + 1
         e <- b + k - 1 
         D[i, b:e] <- Da[i,,pi[i]]
         }
      Ppi <- K %*% D
		V <- solve.mp.function(rb, Ppi, df)

      # policy improvement
      V <- D %*% V

      for (a in 1:num.actions) {
         Q[,a] <- R[,a] + df * P[,,a] %*% V
         }
      pi.old <- pi
      pi <- apply(Q,1,which.max)
      it <- it + 1
      }
   if (verbose) {
      print(" ")
      print("Done.")
      }
   list(pi = pi, Q = Q, mse = mse)
   }
   

pole.pisf <- function(n.grid, ST, m, tau,
								sars = NULL, 
								num.episodes = 10^4, 
							 	double =FALSE,
							 	df = 0.99,
							 	max.iter.pi = 15,
							 	max.steps.test = 3000,
							 	M = NULL,
							 	pi = NULL,
							 	it.kmeans = 20,
								selection.function = select.pole.clara,
                        keep.all = FALSE,                
							 	verbose = FALSE,
							 	...){
							 
	if (is.null(sars)) {
		if (verbose) print("Generating data...")
		sars <- pole.data(num.episodes, normalize = FALSE, lsars = FALSE,
								double = double, verbose = verbose, ...)
		}
	
	if (is.null(M)) {
		if (verbose) print("Building the MDP...")
		M <-  pole.build.mdp.grid(n.grid,sars)  	
		}
	
	
	S <- uniform.grid(matrix(0,nrow(M$P),ncol(ST)),M$s.max,-M$s.max)
	
	
	if (is.null(pi)) {
		if (verbose) print("Running PISF algorithm...")
		pi <- arch.kernel.policy.iteration.alt(M$R, M$P, df=df, m = m, S = S, 
                            N = M$N, keep.all = keep.all,
									 max.iter = max.iter.pi,
									 tau = tau, selection.function = selection.function,
									 s.max = M$s.max, s.delta = M$s.delta, n.grid =
									 M$n.grid)$pi
		}
	
	if (verbose) print("Testing the policy...")
	
	ns <- array(0, nrow(ST))
	
	for (i in 1:nrow(ST)) {
		ns[i] <- pole.test.policy.grid(pi, ST[i,], M$s.max, M$s.delta,
													M$n.grid, max.steps = max.steps.test,...)
		if (verbose) print(paste("s",i,ns[i]))
		}
	ns
	}


pole.experiment.pisf2 <- function(n.grid, M = NULL, sars= NULL,
											ST = pole.make.ST(3),
											m = seq(100,600,by=100), 
											tau = c(0.01),
			  								num.avg = 1,
				 							filename = "pole.pisf",
											dir = "./res_tese/pole/",
											selection.function = select.pole.clara,
			  								verbose = TRUE
											){
# See script below...											
	if (is.null(M)) M <- pole.build.mdp.grid(n.grid, sars)
	for (i in 1:length(m)) {
		if (verbose) print(paste("m =", m[i],"(", i,"/",length(m),")"))
		for (t in 1:length(tau)) {
			if (verbose) print(paste("tau =", tau[t],"(", t,"/",length(tau),")"))
			ns <- matrix(0, num.avg, nrow(ST))
			for (j in 1:num.avg) {
				ns[j,] <- pole.pisf(Inf, ST, m[i], tau[t], 0, verbose=FALSE,
												M=M, selection.function=selection.function)
												
				if (verbose) print(paste("m =",m[i],"tau =", tau[t], "run =",j, 
								"mean =", mean(ns[j,]))) 
				}
			if (verbose) print(paste("m =", m[i], "tau =", tau[t], "mean =",
									mean(ns)))
			wt(ns, paste(dir,filename,"_n",n.grid,
                        "_m",m[i],"_t",tau[t],".txt",sep=""))
			}
		}
	if (verbose) print("All done!")
	}
													
pole.experiment.pisf <- function(sars, n.grid = c(3,4,5,6,7), pm =
										c(0.2,0.4,0.6,0.8), tau=0.01) {
	for (ng in n.grid) {
		m <- round(pm * ng^4)
		pole.experiment.pisf2(ng, sars=sars, m=m, tau=tau)
		}
	}
		

		
print.line <- function(filename, max.steps=3000) {
	t <- as.matrix(read.table(filename))
	dim(t) <- length(t)
	print(paste(round(sum(t==max.steps)/length(t)*100,d=2),
               round(mean(t),d=2),
					round(max(t), d=2),
	            round(min(t), d=2),
				   round(sd(t), d=2),
 				   sep=" & "))
	}
			
   
print.table.grid <- function(n.grid = c(3,4,5,6,7),
				 				filename = "pole.grid",
								dir = "./res_tese/pole/",
								max.steps=3000) {
	for (i in 1:length(n.grid)) {
		print.line(paste(dir,filename,"_n",n.grid[i],".txt",sep=""), max.steps)
		}   
	} 
		  				

print.table.pisf <- function(n.grid = c(3,4,5,6,7),
								pm = c(0.2,0.4,0.6,0.8),
								tau = c(0.01),
				 				filename = "pole.pisf",
								dir = "./res_tese/pole/",
								max.steps=3000) {
	for (i in 1:length(n.grid)) {
		for (j in 1:length(pm)) {
			m <- round(pm[j] * n.grid[i]^4)
			print.line(paste(dir,filename,"_n",n.grid[i],"_m", m,"_t",tau,
							".txt",sep=""), max.steps)
		   }   
	   } 
  }
			  				
  
	
pole.arch.clara <- function(S, num.archs, goals, p.goals = 0.1, stand = TRUE,
									metric = "euclidean") {
   g <- (1:length(goals))[goals]
   num.goals <- round(p.goals * num.archs)
   if (length(g) > num.goals) g <- sample(g,num.goals)
   
   num.archs <- num.archs - length(g)
	id <- clara(S, num.archs, metric = metric, stand = stand, medoids.x =
				FALSE)$i.med
	c(g,id)
	}

  

pole.arch.random <- function(S, num.archs, goals, p.goals = 0.1, stand = TRUE,
										metric = "euclidean") {
   g <- (1:length(goals))[goals]
   num.goals <- round(p.goals * num.archs)
   if (length(g) > num.goals) g <- sample(g,num.goals)
   num.archs <- num.archs - length(g)
   c(g,sample(1:nrow(S), num.archs))
	}
	
	
   
pole.kbrl <- function(ST, m, tau,
								lsars = NULL, 
								num.episodes = 10^4, 
							 	double =FALSE,
							 	df = 0.99,
							 	max.iter.pi = 15,
							 	verbose = FALSE,
							 	max.steps.test = 3000,
							 	pi = NULL,
								selection.function = pole.arch.clara,
        						p.goals = 0.1,
							 	...){
							 	
# lsars should be normalized
			
   
	if (is.null(lsars)) {
		if (verbose) print("Generating data...")
		lsars <- pole.data(num.episodes, normalize = TRUE, lsars = TRUE,
							double = double, verbose = verbose, ...)
		}
	
 	k <- m %/% lsars$num.actions
	m <- k * lsars$num.actions
   
   if (verbose) print("Selecting states...")
   id <- NULL
	for (a in 1:lsars$num.actions) {
      id <- c(id,
       list(selection.function(lsars[[a]]$s, num.archs=k, stand =FALSE,
				goals = lsars[[a]]$g, p.goals = p.goals)))
   	}
  	
   for (a in 1:lsars$num.actions) {
		lsars[[a]]$s  <- lsars[[a]]$s[id[[a]],]
		lsars[[a]]$a  <- lsars[[a]]$a[id[[a]]]
		lsars[[a]]$r  <- lsars[[a]]$r[id[[a]]]
		lsars[[a]]$s2 <- lsars[[a]]$s2[id[[a]],]
		lsars[[a]]$g  <- lsars[[a]]$a[id[[a]]]
		}
		
	
   if (verbose) print("Generating the MDP...")
	M <- kbrl.mdp(lsars, tau = tau, goal.reward = -1) 
	
	if (verbose) print("Running policy iteration...")
   Q <- policy.iteration(M$R, M$P, df = df, pi = NULL, max.iter = max.iter.pi,
                     solve.mp.function = solve.mp)$Q
   V <- apply(Q,1,max)
   
   rbfns <- NULL
   for (a in 1:lsars$num.actions) {
      rbfns <- c(rbfns,
         		list(make.rbfn.kbrl(lsars[[a]]$s, tau, num.actions=1)))
      b <- (a-1) * k + 1
      e <- b + k - 1
      rbfns[[a]]$w <- lsars[[a]]$r + df * V[b:e]
  	   }
   
   if (verbose) print("Testing the policy...")
	
	ns <- array(0, nrow(ST))
	
	for (i in 1:nrow(ST)) {
		ns[i] <- pole.test.policy.multiple.rbfns(rbfns, ST[i,], means =
						lsars$means, stdevs = lsars$stdevs, 
      				max.steps = max.steps.test)
		if (verbose) print(paste("s",i,ns[i]))
		}
	ns
	}


pole.experiment.kbrl <- function(lsars,  
											ST = pole.make.ST(3),
											m = seq(600,1500,by=300), 
											tau = c(0.3),
			  								num.avg = 20,
				 							filename = "pole.kbrl",
											dir = "./res_tese/pole/",
			  								verbose = TRUE,
                                 max.iter.pi = 15,
                                 max.steps.test = 3000,
                                 selection.function = pole.arch.clara,
                                 p.goals =  0.1
                                 ){
# lsars must be normalized!!
	for (i in 1:length(m)) {
		if (verbose) print(paste("m =", m[i],"(", i,"/",length(m),")"))
		for (t in 1:length(tau)) {
			if (verbose) print(paste("tau =", tau[t],"(", t,"/",length(tau),")"))
			ns <- matrix(0, num.avg, nrow(ST))
			for (j in 1:num.avg) {
			
			
				ns[j,] <- pole.kbrl(ST, m[i], tau[t], lsars, 
                     max.iter.pi = max.iter.pi, max.steps.test =
							max.steps.test, selection.function=selection.function,
							p.goals = p.goals, verbose=FALSE)
				
				if (verbose) print(paste("m =",m[i],"tau =", tau[t],"run =",j, 
								"mean =", mean(ns[j,]))) 
				}
			if (verbose) print(paste("m =", m[i], "tau =", tau[t], "mean =",
									mean(ns)))
			wt(ns, paste(dir,filename,"_m",m[i],"_t",tau[t],".txt",sep=""))
			}
		}
	if (verbose) print("All done!")
	}
													
	
											
print.line.kbrl <- function(filename, max.steps=3000) {
	t <- as.matrix(read.table(filename))
	tt <- apply(t, 1, mean)
	dim(t) <- length(t)
	print(paste(round(sum(t==max.steps)/length(t)*100,d=2),
               round(mean(t),d=2),
					round(max(t), d=2),
	            round(min(t), d=2),
				   round(sd(t), d=2),
				   round(sum(tt==max.steps)/length(tt)*100,d=2),
 				   sep=" & "))
	}
			

print.table.kbrl <- function(m = seq(600,1500,by=300), 
								tau = c(0.3),
				 				filename = "pole.kbrl",
								dir = "./res_tese/pole/",
								max.steps=3000) {
	for (i in 1:length(m)) {
		print.line.kbrl(paste(dir,filename,"_m",
								m[i],"_t",tau,".txt",sep=""), max.steps)
		}   
	}
	

                  

pole.kbsf <- function(ST, m, tau.a, tau.q,
								lsars = NULL, 
								p.data = 0.1,
								p.goals = NULL,
								num.episodes = 2000, 
							 	double =FALSE,
							 	df = 0.99,
							 	max.iter.pi = 15,
							 	max.steps.test = 3000,
							 	pi = NULL,
								selection.function = arch.grid,
								p.centers = 0.3,
							 	verbose = FALSE,
							 	...){
							 	
# lsars should be normalized
	if (is.null(lsars)) {
		if (verbose) print("Generating data...")
		lsars <- pole.data(num.episodes, normalize = TRUE, lsars = TRUE,
							double = double, verbose = verbose, ...)
		}
	
	# reduce the dataset keeping the proportion of 'p.goals' goals
	calc.prop <- FALSE
	if (is.null(p.goals)) calc.prop <- TRUE
	if (p.data < 1) {
      for (a in 1:lsars$num.actions) {
         n <- nrow(lsars[[a]]$s)
         p <- round(p.data * n)
         if (calc.prop) {
         	p.goals <- sum(lsars[[a]]$g) / nrow(lsars[[a]]$s)
         	}
         ig <- (1:n)[lsars[[a]]$g]
         num.goals <- min(max(round(p.goals * p),1),p)
         ind <- ig[sample(1:length(ig),num.goals)]
         
         p <- p - num.goals
         ind <- c(ind, sample(1:n,p))
         
         lsars[[a]]$s <- lsars[[a]]$s[ind,]
         lsars[[a]]$a <- lsars[[a]]$a[ind]
         lsars[[a]]$r <- lsars[[a]]$r[ind]
         lsars[[a]]$s2 <- lsars[[a]]$s2[ind,]
         lsars[[a]]$g <- lsars[[a]]$g[ind]
         } 
      }
	
	
	rbfn <- kbrl.arch(lsars = lsars, tau.a = tau.a, tau.q = tau.q, df = df,
							num.archs = m, C = NULL, arch.function=selection.function, 
                     iter.pi = max.iter.pi, 
                  	goal.reward = -1, run.value.iteration=FALSE, 
                  	p.centers = p.centers,
                  	verbose = verbose)
   
   
   if (verbose) print("Testing the policy...")
	
	ns <- array(0, nrow(ST))
	for (i in 1:nrow(ST)) {
		ns[i] <- pole.test.policy.rbfn(rbfn, ST[i,], means =
						lsars$means, stdevs = lsars$stdevs, 
      				max.steps = max.steps.test)
		if (verbose) print(paste("s",i,ns[i]))
		}
	ns
	}


pole.experiment.kbsf <- function(ST, 
									lsars,
									num.archs = c(25,50,75,100,125),
									p.data = c(1),
									p.goals = NULL,
									p.centers = 0.02,
									tau.a = 0.3,
									tau.q = 0.01,
						   	 	max.iter.pi = 15,
									selection.function = arch.clara,
									num.avg = 20,
									dir = "./res_tese/pole/",
									filename = "pole.kbsf",
									verbose=TRUE) {
	if (length(p.centers) == 1) p.centers <- rep(p.centers, length(p.data))
	for (m in 1:length(num.archs)) {
		for (pd in 1:length(p.data)) {
			res <- matrix(0, num.avg, nrow(ST))
			for (i in 1:num.avg) {
            res[i,] <- pole.kbsf(ST, num.archs[m], tau.a, tau.q,
            							lsars = lsars, p.data=p.data[pd], 
            							p.centers = p.centers[pd],
											selection.function = selection.function,
											p.goals = p.goals, max.iter.pi = max.iter.pi)
            if (verbose) print(paste(i, "Num. archs", num.archs[m], "Prop.data",
									p.data[pd], ":",mean(res[i,])))
				}
			wt(res,paste(dir,filename,"_m",num.archs[m],"_pd",
								p.data[pd],".txt",sep=""))
			}
		}
	if (verbose) print("All done.")
	} 
				
											

print.table.kbsf <- function(num.archs = c(25,50,75,100),
								p.data = 1,
								filename = "pole.kbsf",
								dir = "./res_tese/pole/",
								max.steps=3000) {
	for (i in 1:length(num.archs)) {
		for (j in 1:length(p.data)) {
         print.line.kbrl(paste(dir,filename,"_m", num.archs[i],
								"_pd",p.data[j],".txt",sep=""), max.steps)
			}
		}   
	}
	

pi.op <- function(n, m, num.actions) {
	30 * n +                    # compute indexes
	num.actions * (2 * m) + 	 #  compute Pa
	num.actions * (m ^ 2)       # compute Ra   
	}
	
pisf.op <- function(n, m, k, num.actions) {
	pi.op(n,m,num.actions) +    # generates the initial MDP
	m * log(m) + 					 # order rows
	2 * num.actions * m * k     # builds Da and K	
	}

clara.op <- function(n,k,it,c=40) {
 	it * (k * (c + k)^2 + k * (n - k))	
	}
	
kbrl.op <- function(n, m, num.actions) {
	num.actions * clara.op(n,round(m/num.actions),20) + 
 	num.actions * m^2  # compute Pa 
	}

kbsf.op <- function(n, m, num.actions) {
	clara.op(n,m,20) +            # cluster
	num.actions * (2 * n * m) + # compute Da and Ka
	num.actions * (m^2 * n)    # multiply Da and Ka
	}
	


generate.bar.plot.pi <- function(names=c("PI","PISF"),
										n = 94787,
										m = c(1296,778),
										num.actions = 3,
										col = c("LIGHTGREY", "DARKGREY"),
										filename.eps="./fig_tese/pole_bar_pi") {
	 # do not change the order...
	 t <- matrix(0, 2, 2, byrow = TRUE)
	 rownames(t) <- c("Definição", "Solução")
	 colnames(t) <- names
    t[1,1] <- pi.op(n,m[1], num.actions)
	 t[1,2] <- pisf.op(n,m[1],m[2], num.actions)
    t[2,] <- 15 * m^3
    t[2,2] <- 15 * (m[2]^2 * m[1] + m[1] * m[2]) 
    
 	 print(t)
 	 x11(w=4,h=7)
    barplot(t[1,], ylab="Número de operações", col=c(col[1],col[1]),
				legend=FALSE)
    dev.copy2eps(file=paste(filename.eps, "_def.eps",sep=""))
	 print(paste("Definição:", t[1,2]/t[1,1]))
    
    barplot(t, ylab="Número de operações", col=col, legend=TRUE)
#     legend(locator(1), legend=rownames(t), fill=col)
    dev.copy2eps(file=paste(filename.eps, ".eps",sep=""))
	 print(paste("Total:", sum(t[,2]) / sum(t[,1])))
	 } 			
	
	
generate.bar.plot.kbrl <- function(names=c("KBRL","KBSF"),
										n = 94787,
										m = c(1500,75),
										num.actions = 3,
										col = c("LIGHTGREY", "DARKGREY"),
										filename.eps="./fig_tese/pole_bar_kbrl") {
	 # do not change the order...
	 t <- matrix(0, 2, 2, byrow = TRUE)
	 rownames(t) <- c("Definição", "Solução")
	 colnames(t) <- names
	 t[1,1] <- kbrl.op(n,m[1], num.actions)	 
    t[1,2] <- kbsf.op(n,m[2], num.actions)
    t[2,] <- 15 * m^3
    
    print(t)
    x11(w=4,h=7)
    barplot(t[1,], ylab="Número de operações", col=c(col[1],col[1]),
				legend=FALSE)
    dev.copy2eps(file=paste(filename.eps, "_def.eps",sep=""))
	 print(paste("Definição:", t[1,2]/t[1,1]))
    print(paste("Solução:", t[2,2]/t[2,1]))
    
    barplot(t, ylab="Número de operações", col=col, legend=TRUE)
  #  legend(locator(1), legend=rownames(t), fill=col)
    dev.copy2eps(file=paste(filename.eps, ".eps",sep=""))
	 print(paste("Total:", sum(t[,2]) / sum(t[,1])))
	 } 			
				

pole.double.experiment <- function( sars = NULL, lsars = NULL, num.episodes=500,
												ST = NULL, n.grid = 3,	
												p.trans = 0.3,
												num.archs = c(250,200,150,100,50),
												tau.a = 0.3,
												tau.q = 0.01,
												tau.lspi = 0.95,
												selection.function = arch.kmeans,
												p.centers = 0.2,
												max.iter.pi = 15,
												num.avg = 20,
												df = 0.999,
												max.steps.test = 3000,
												A = c(-10,10),
												run = c(TRUE,TRUE),
												max.iter.selection = 20, 
												dir = "./res_tese/pole/double/",
												filename = c("pole.double.lspi",
																"pole.double.kbsf"),
												verbose=TRUE) {
# the samples (including sars) must be normalized

	if(is.null(sars)) {
		if (verbose) print("Generating data...")
		sars <- pole.data(num.episodes, normalize = TRUE, lsars = FALSE, 
                       double = TRUE, p=p.trans, p2=0, A=A)
   	}
   
   if (is.null(lsars)) lsars <- make.lsars(sars)
   
#    num.episodes <- sum(sars$g)
      	
	if (is.null(ST)) {
		if (verbose) print("Generating test set...")
		max.pos <- p.trans * c(2.4, 2.4, pi/5,pi/5)
		min.pos <- -max.pos
		grid <- matrix(0, n.grid^4, 4)
 		ST <- uniform.grid(grid, max.pos, min.pos)
		ST <- cbind(ST, 0, 0)
      }
  	
	for (m in 1:length(num.archs)) {
		res.lspi <- array(0, c(num.avg, nrow(ST), length(tau.lspi)))
		res.kbsf <- matrix(0, num.avg, nrow(ST))
		
		for (i in 1:num.avg) {
			if (verbose) print(paste("Selecting archetypes for m = ",num.archs[m]))
         C <-selection.function(lsars, num.archs[m],max.iter=max.iter.selection)
         if (run[1]) {
            for (t in 1:length(tau.lspi)) {
            	rbfn <- make.rbfn.kbrl(C, tau.lspi[t], 
                 							num.actions = lsars$num.actions, 
                                    num.neighbors = ncol(C)+1,
                                    same.width = TRUE, p.centers = 1)
               rbfn <- lspi.rbfn(sars, rbfn, df = df, 
               						num.iterations = max.iter.pi,
                                 precision = 1e-6, verbose = FALSE, 
                                 rbfn.dm = rbfn.norm.design.matrix)
               
               for (j in 1:nrow(ST)) {
                  res.lspi[i,j,t] <- pole.test.policy.rbfn(rbfn, ST[j,], means =
                              lsars$means, stdevs = lsars$stdevs, 
                              max.steps = max.steps.test, double = TRUE, A=A)
                  }
               if (verbose) print(paste(i,"LSPI m", num.archs[m],"t",
										tau.lspi[t],":",mean(res.lspi[i,,t])))
			    	}
			    }
                           
         if (run[2]) {
            rbfn <- kbrl.arch(lsars = lsars, tau.a = tau.a, tau.q = tau.q, 
            				df = df, num.archs = num.archs[m], C = C,              
        						arch.function=NULL, 
                        iter.pi = max.iter.pi, num.neighbors=ncol(C)+1, 
                        goal.reward = -1, run.value.iteration=FALSE, 
                        p.centers = p.centers,
                        verbose = FALSE)
         
            	for (j in 1:nrow(ST)) {
                  res.kbsf[i,j] <- pole.test.policy.rbfn(rbfn, ST[j,], means =
                              lsars$means, stdevs = lsars$stdevs, 
                              max.steps = max.steps.test, double = TRUE, A=A)
                  }
            if (verbose) print(paste(i,"KBSF m", num.archs[m],
											mean(res.kbsf[i,])))
            }
         }
      
      if (run[1]) {
         for (t in 1:length(tau.lspi)) {
            wt(res.lspi[,,t],paste(dir,filename[1],"_n",num.episodes,
                           "_m",num.archs[m],"_t", tau.lspi[t],".txt",sep=""))
            }
         }
      if (run[2]) wt(res.kbsf,paste(dir,filename[2],"_n",num.episodes,
                        "_m",num.archs[m],"_ta",tau.a, 
                        "_tq", tau.q, ".txt",sep=""))
		}
	}
												


print.table.double <- function(num.episodes= 1000,
												num.archs = c(50,100,150,200,250),
												tau.a = 0.3,
												tau.q = 0.01,
												tau.lspi = c(0.95),
												max.steps.test = 3000,
												run = c(TRUE,TRUE),
												dir = "./res_tese/pole/double/",
												filename =c("pole.double.lspi",
                                                "pole.double.kbsf")) {

	for (n in num.episodes) {
		for (m in num.archs) {
			print.line.kbrl(paste(dir,filename[1],"_n",n,"_m", m,"_t",
										 tau.lspi,".txt", sep=""), max.steps.test)
			}
		}
	
	for (n in num.episodes) {
		for (m in num.archs) {
			print.line.kbrl(paste(dir,filename[2],"_n",n,"_m", m,".txt",
								sep=""), max.steps.test)
			}
		}

	}


print.table.double.single <- function(num.episodes= 1000,
												num.archs = c(50,100,150,200),
												tau.a = 0.3,
												tau.q = 0.01,
												tau.lspi =0.95,
												max.steps.test = 3000,
												dir = "./res_tese/pole/double/",
												filename =c("pole.double.lspi.single",
                                                "pole.double.kbsf.single"),
                                    run = c(TRUE, TRUE)) {

	if (run[1]) {
      for (n in num.episodes) {
         for (m in num.archs) {
         	for (t in tau.lspi) {
               print.line(paste(dir,filename[1],"_n",n,"_m", m,"_t",
                  t,".txt", sep=""), max.steps.test)
               }
            }
         }
      }
	
	if (run[2]) {
      for (n in num.episodes) {
         for (m in num.archs) {
            print.line(paste(dir,filename[2],"_n",n,"_m",
						m,".txt",sep=""), max.steps.test)
            }
         }
		}
	}
	
	
print.table.lspi <- function(num.episodes= 1000,
										n = 17319,
      								num.archs = c(50,100,150,200),
										tau.lspi =c(0.01, seq(0.1,0.9,by=0.1),0.95,0.99),
										max.steps.test = 3000,
										dir = "./res_tese/pole/double/lspi/test_set/",
										filename =c("pole.double.lspi")) {
   for (m in num.archs) {
      for (t in tau.lspi) {
         cat(paste("LSPI(" ,n, "," ,m, ") &", t))
         print.line.kbrl(paste(dir,filename[1],"_n",num.episodes,"_m", m,"_t",
            t,".txt", sep=""), max.steps.test)
         }
      }
	
	}

# LSPI: m^2 * n + (m * |A|)^3 (per iteration) + policy improvement(?)
# KBSF: |A| * (m * n + m^2 * n) (initial)
#              + m^3(per iteration) + policy improvement


split.files.merged <- function(num.episodes= 1800,
										num.archs = c(50,100,150,200,250),
										tau.a = 0.3,tau.q = 0.01,
										tau.lspi =0.95,
										dir = "./res_tese/pole/double/",
										filename = c("pole.double.lspi.merged",
                                          "pole.double.kbsf.merged"),
                              suffix = "single" ,run = c(TRUE, TRUE) ) {
# The experiments with the merged samples were carried out together. 
# This function splits the data referring to the single state and test set cases
   if (run[1]) {
         for (n in num.episodes) {
            for (m in num.archs) {
               for (t in tau.lspi) {
                  test.name <- paste(dir,filename[1],"_n",n,"_m", m,"_t",
                     t,".txt", sep="")
                  T <- read.table(test.name)
                  wt(T[,1:(ncol(T)-1)], test.name)
                  single.name <- paste(dir,filename[1],".single_n",n,"_m",
											m,"_t",t,".txt", sep="")
                  wt(T[,ncol(T)], single.name)
                  }
               }
            }
         }
      
      if (run[2]) {
         for (n in num.episodes) {
            for (m in num.archs) {
               test.name <- paste(dir,filename[2],"_n",n,"_m",
                     				m,".txt",sep="")
                T <- read.table(test.name)
                wt(T[,1:(ncol(T)-1)], test.name)
                single.name <- paste(dir,filename[2],".single_n",n,"_m",
                     				m,".txt",sep="")
                wt(T[,ncol(T)], single.name)
               }
            }
         }
   }
      
	

pole.double.data.uniform <- function(n=17319, p.goals= NULL, ng = 0, p = 1,
									  normalize = TRUE, 
									  max.pos=c(2.4, 2.4, pi/5,pi/5, pi/5,pi/5),
									  min.pos = -max.pos, A = c(-10,10)) {
	
	max.pos <-  p * max.pos
	min.pos <-  p * min.pos
	
	if (is.null(ng)) ng <- round(p.goals*n)
	
	sars <- make.empty.sars()
	
	if (ng > 0) {
      s <- matrix(0, ng, length(max.pos))
      for (i in 1:length(max.pos)) {
         s[,i] <- runif(ng, -max.pos[i], max.pos[i])
         }
      sars <- collect.episodes(s, pole.double.transition, c(-10,10),
                                 normalize=FALSE, double = TRUE)
      sars$s <- sars$s[sars$g,]
      sars$a <- sars$a[sars$g]
      sars$r <- sars$r[sars$g]
      sars$s2 <- sars$s2[sars$g,]
      sars$g <- sars$g[sars$g]
      }
	
	
	n <- n - length(sars$g)
	s <- matrix(0, n, length(max.pos))
	for (i in 1:length(max.pos)) {
		s[,i] <- runif(n, -max.pos[i], max.pos[i])
		}
		
	sa <- list(s = s, a = sample(1:2, nrow(s), TRUE))
	sars2 <- collect.transitions(sa, pole.double.transition, c(-10,10),
										 normalize=FALSE, sd=0, double=TRUE)
	
	merge.sars(sars, sars2, unnormalize=c(FALSE, FALSE), normalize=normalize)
	}
	

pole.double.data.grid <- function(n=17319, ng = 1000, n.grid=5, p = 1,
									  normalize = TRUE, 
									  max.pos=c(2.4, 2.4, pi/5,pi/5, pi/5,pi/5),
									  min.pos = -max.pos, A = c(-10,10)) {
	S <- uniform.grid(matrix(0, n.grid^6, 6), p * max.pos, p * min.pos, p=0)
	
	list(s=S, a= sample(1:2, nrow(S), replace=TRUE))
	sars <- collect.transitions(sa, pole.double.transition, c(-10,10),
										normalize=FALSE, sd=0, double=TRUE)
 	sars2 <- pole.double.data.uniform(n-length(sars$r), ng = 0,
 										p=p, normalize=FALSE, max.pos=max.pos,
 										min.pos=min.pos, A=A)
 	merge.sars(sars, sars2, unnormalize=c(FALSE, FALSE), normalize=TRUE)
	}
	
	

pole.double.data.terminal <- function(udsars, n=200, pos=c(2.4,-2.4),
										ang=c(pi/5,-pi/5), delta=1e-6) {
## "udsars" can not be normalized!!!
	max.others <- apply(udsars$s2,2,max)
	min.others <- apply(udsars$s2,2,min)
	
	S <- matrix(0, 4 * n, 6)
	
	for (i in 1:ncol(S)) {
		S[,i] <- runif(4 * n, min.others[i], max.others[i])
		}
		
	S[,5] <- 0
	S[,6] <- 0
	
	S[1:n,1] <- pos[1] - delta # delta is to make sure the transition will happen
	S[(n+1):(2*n),1] <- pos[2] + delta
	S[(2*n+1):(3*n),3] <- ang[1] - delta
	S[(3*n+1):(4*n),3] <- ang[2] + delta
		
   sars <- collect.episodes(S, pole.double.transition, c(-10,10),
							  normalize=FALSE, double=TRUE)
   sars$s <- sars$s[sars$g,]
   sars$a <- sars$a[sars$g]
   sars$r <- sars$r[sars$g]
   sars$s2 <- sars$s2[sars$g,]
   sars$g <- sars$g[sars$g]

	merge.sars(sars, udsars, unnormalize=c(FALSE, FALSE), normalize=FALSE)
	}
		

generate.samples <- function(load.first.from.file = TRUE,
									  num.episodes = seq(1000,2000, by=200),
									  p.trans=0.3, p.trans2 = 0, A=c(-10,10),
									  dir = "./files/",
									  filename= "sars.pole.double") {
	sars <- NULL
	b <- 2
	generated <- 0
	if (load.first.from.file) {
		sars <- load.sars(paste(dir,filename,".",num.episodes[1],sep=""))
		generated <- num.episodes[1]
		}
	else b <- 1
	
	for (i in b:length(num.episodes)) {
		tsars <- pole.data(num.episodes[i] - generated, normalize = FALSE, 
								   lsars = FALSE, double = TRUE, 
                       		p=p.trans, p2=p.trans2, A=A)
		sars <- merge.sars(sars, tsars, c(TRUE, FALSE), TRUE)
		save.sars(sars, paste(dir,filename,".",num.episodes[i],sep=""))
		generated <- num.episodes[i]
		}
	}
	  

arch.double.fixed <- function(lsars, num.archs = NULL, max.iter=NULL,
											  np=c(2,2,2,2,3,3)) {
	S <- lsars[[1]]$s2
	for (a in 2:lsars$num.actions) {
		S <- rbind(S, lsars[[a]]$s2)
		}
	max.pos <- apply(S, 2, max) 
	min.pos <- apply(S, 2, min)
	delta <- (max.pos - min.pos) / (np+1)
	
	C <- matrix(0, prod(np), 6)
	r <- 1
	for (i1 in 1:np[1]) {
	   for (i2 in 1:np[2]) {
         for (i3 in 1:np[3]) {
            for (i4 in 1:np[4]) {
               for (i5 in 1:np[5]) {
                  for (i6 in 1:np[6]) {
                  	i <- c(i1,i2,i3,i4,i5,i6)
                     C[r,] <- min.pos + i * delta
                     r <- r + 1
                     }
                  }
               }
            }
         }
      }
   C
   }
	
			
									  


script.sample <- function(num.episodes=seq(1000,2000, by=200),
									ST,
									num.archs = 144, 
									p.centers = 0.2 / c(1,1.2,1.4,1.6,1.8,2)^2, 
									dir="./files/",
									filename="sars.pole.double"){

	for (i in 1:length(num.episodes)) {
		sars <- load.sars(paste(dir,filename,".",num.episodes[i],sep=""))
		
		pole.double.experiment(sars, num.episodes=num.episodes[i], ST=ST,
         num.archs=num.archs, num.avg = 1, 
         selection.function = arch.double.fixed,		
         filename=c("pole.double.lspi.fixed","pole.double.kbsf.fixed"),	
         p.center=p.centers[i])
      }
   }



generate.samples.uniform <- function(load.first.from.file = FALSE,				
                                 num.trans=seq(10000,50000,by=10000),
											dir="./files/",
											filename="sars.pole.double.uniform"){
   if (!load.first.from.file) {
      usars <- pole.double.data.uniform(n=num.trans[1])
      save.sars(usars,paste(dir,filename,".",num.trans[1],sep=""))
      }
  else usars <- load.sars(paste(dir,filename,".",num.trans[1],sep=""))
   
   for (n in 2:length(num.trans)) {
   	num <- num.trans[n] - num.trans[n-1]
      usars <- merge.sars(usars, pole.double.data.uniform(n=num),
      						c(TRUE, TRUE), TRUE)
      save.sars(usars,paste(dir,filename,".",num.trans[n],sep=""))
      }
   }

												

script.sample.uniform <- function(
											num.trans=seq(10000,100000,by=10000),
											ST,
										   p.centers = 0.8 / (1:10)^2,
											dir="./files/",
											filename="sars.pole.double.uniform",
											tau.q = 0.6,
											tau.a=0.3,
											tau.lspi = 0.95, 
											num.avg = 20,
											run=c(TRUE, TRUE)){
	for (i in 1:length(num.trans)) {
		sars <- load.sars(paste(dir,filename,".",num.trans[i],sep=""))
		pole.double.experiment(sars, num.episodes=num.trans[i], ST=ST,
         num.archs=100,p.center=p.centers[i], tau.a=tau.a, tau.q=tau.q,
         tau.lspi = tau.lspi, run=run, num.avg = num.avg,
         filename=c("pole.double.lspi.uniform",
							"pole.double.kbsf.uniform"))
      }
   }
   
   
   
   
script.sample.uniform.preliminar <- function(ST, 
								tau.lspi=c(0.01, seq(0.1,0.9,by=0.1), 0.9, 0.95,0.99),
								tau.q = c(0.01, seq(0.1,0.9,by=0.1), 0.9, 0.95,0.99),	
								tau.a=0.2, 
								num.avg=5,
								run=c(TRUE, TRUE)) {
	if (run[1]) {
      for (t in tau.lspi) {
         script.sample.uniform(num.trans=20000, ST=ST, tau.lspi=t,
										 num.avg=num.avg, run=c(TRUE, FALSE))
         }
      }
	
	if (run[2]) {
      for (t in tau.q) {
         script.sample.uniform(num.trans=20000, ST=ST, tau.q=t, num.avg =
							num.avg, p.centers = 0.2,tau.a=tau.a, run=c(FALSE, TRUE))
         }
      }
	}
	
print.table.uniform.preliminar <- function(							
							   tau.lspi=c(0.01,seq(0.1,0.9,by=0.1), 0.9, 0.95,0.99),
								tau.q = c(0.01, seq(0.1,0.9,by=0.1), 0.9, 0.95,0.99),	
								tau.a=0.2, 
								run=c(TRUE, TRUE),
								max.steps=3000,
								dir="./res_tese/pole/double/uniform_pre/",
								filename=c("doublepole.double.lspi.uniform",
                                 "doublepole.double.kbsf.uniform")) {
	if (run[1]) {
      for (t in tau.lspi) {
      	cat(t)
      	print.line.kbrl(paste(dir,filename[1],"_n20000_m100_t",t,
											".txt",sep=""))
         }
      }
	
	if (run[2]) {
      for (t in tau.q) {
      	cat(t)
      	print.line.kbrl(paste(dir,filename[2],"_n20000_m100_ta",tau.a, 
      								 "_tq", t,".txt", sep=""))
         }
      }
	}


print.table.samples <- function(	num.episodes = c(1800,1000,1000),
											num.archs = 150,
										   tau.lspi = 0.95,
											dir ="./res_tese/pole/double/",
                                 filename = c("pole.double.lspi",
                                             "pole.double.kbsf"),
                                 suffixes = c("merged", "uniform","grid"),
                                 max.steps = 3000) {
	for (s in 1:length(suffixes)) {
		fn <- paste(dir,filename[1],".",suffixes[s],sep="")
		fn <- paste(fn,"_n",num.episodes[s],"_m",num.archs,"_t",tau.lspi,".txt",
						sep="")
		print.line.kbrl(fn, max.steps = max.steps)
		}
		
	for (s in 1:length(suffixes)) {
		fn <- paste(dir,filename[2],".", suffixes[s],sep="")
		fn <- paste(fn,"_n",num.episodes[s],"_m",num.archs,".txt",sep="")
		print.line.kbrl(fn, max.steps = max.steps)
		}
	}
                                 

plot.sample.size.series <- function(
                                    num.episodes = seq(1000,4000,by=1000),
 												num.archs = 100,
 												tau.lspi =0.95,
 												tau.a=0.3,
 												tau.q = 0.01,
												dir="./res_tese/pole/double/sample_size/",
												filename = c("pole.double.lspi",
            			                          "pole.double.kbsf"),
								filename.eps=c("./fig_tese/pole_sample.eps"),
							   					xlab="Número de episódios"){
	D <- matrix(0, length(num.episodes), 2)
	m <- num.archs
	for (i in 1:length(num.episodes)) {
		n <- num.episodes[i]
       t <- as.matrix(read.table(paste(dir,filename[1],"_n",n,"_m", m,"_t",
                      tau.lspi,".txt", sep="")))
 		D[i,2] <- mean(t) 
		
		t <- as.matrix(read.table(paste(dir,filename[2],"_n",n,"_m", m,
							".txt",sep="")))
		D[i,1] <- mean(t)      			                          
		}
 	
 	mp(D, t="o",ps=17,x=num.episodes, 
 	   xlab=xlab, ylab="Passos", 
 	   leg=c("KBSF","LSPI"), leg.pos=locator(1))
  	dev.copy2eps(file=filename.eps[1])
	}


plot.sample.size.series.uniform <- function(
                                    num.episodes = seq(10000,60000,by=10000),
 												num.archs = 100,
 												tau.lspi =0.95,
 												tau.a=0.3,
 												tau.q = 0.6,
												dir="./res_tese/pole/double/uniform/",
												filename = c("pole.double.lspi.uniform",
            			                          "pole.double.kbsf.uniform"),
								filename.eps=c("./fig_tese/pole_sample_uniform.eps"),
							   					xlab="Número de transições"){
	D <- matrix(0, length(num.episodes), 2)
	m <- num.archs
	for (i in 1:length(num.episodes)) {
		n <- num.episodes[i]
       t <- as.matrix(read.table(paste(dir,filename[1],"_n",n,"_m", m,"_t",
                      tau.lspi,".txt", sep="")))
 		D[i,2] <- mean(t) 
 		
		t <- as.matrix(read.table(paste(dir,filename[2],"_n",n,"_m", m,
							"_ta", tau.a, "_tq", tau.q, ".txt",sep="")))
		D[i,1] <- mean(t)      			                          
		}
 	
 	mp(D, t="o",ps=17,x=num.episodes, 
 	   xlab=xlab, ylab="Passos", 
 	   leg=c("KBSF","LSPI"), leg.pos="bottomright")
  	dev.copy2eps(file=filename.eps[1])
	}



plot.kbsf.double.pole <- function(num.episodes = 1000,
										num.archs = c(50,250),
										dir = "./res_tese/pole/double/",
										filename="pole.double.kbsf") {
	D <- NULL
	for (m in num.archs) {
		T <- as.matrix(read.table(paste(dir,filename,"_n",num.episodes,"_m", m,
						          ".txt", sep="")))
		D <- cbind(D,apply(T,2,mean))
		D <- cbind(D, mean(mean(T)))
		}
	mp(D, t="o", x=num.archs)
	}
												 	
					
plot.kbsf.double.pole.sd <- function(num.episodes = 1000,
										num.archs = seq(50,250, by=50),
										dir = "./res_tese/pole/double/",
										filename="pole.double.kbsf",
                     filename.eps="./fig_tese/pole_kbsf_sd.eps") {
	D <- matrix(0, length(num.archs),2)
	for (m in 1:length(num.archs)) {
		T <- as.matrix(read.table(paste(dir,filename,"_n",num.episodes,"_m", 
										num.archs[m], ".txt", sep="")))
		D[m,1] <- sd(apply(T, 2, mean))
		D[m,2] <- sd(apply(T, 1, mean))
		}
	mp(D, t="o", x=num.archs, xlab="m", ylab="Desvio-padrão",
		leg = c("Intra-execução", "Inter-execução"), ps=17)
	dev.copy2eps(file=filename.eps)
	}
										
										
										

print("tese.experiments.pole.R loaded")