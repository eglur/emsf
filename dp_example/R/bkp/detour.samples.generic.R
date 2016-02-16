# these are generic versions of the functions implemented on "adp.R", in which
# the inital states of each action are not necessarily the same. I haven't
# worked on this in awhile, so there can be inconsistencies

source("data.manipulation.R")
source("rbfn.gaussian.R")
source("rbfn.rl.R") # only needed if 'kbrl.mountain' is present


make.rbfns.kbrl <- function(lsars, tau = 0.1, old.rbfns = NULL, num.neighbors =
                     3) {
# 'a' in 'sars' is the action index, not the action itself (not zero-based,
# i.e., starts at '1')
	rbfns <- NULL
	num.actions <- length(lsars)
	D <- NULL
	for (a in 1:num.actions) {
		# add a new rbfn to the list
		rbfns <- c(rbfns, list(make.null.rbfn()))
		
		# set the centers
		rbfns[[a]]$c <- lsars[[a]]$s
		
		# set the widths; the width of the i-th RBF will be 'tau' at its nearest neighbor
		num.rbfs <- nrow(lsars[[a]]$s)
				
		dists <- matrix(0, num.rbfs, num.rbfs)
		for (i in 1:(num.rbfs-1)) {
			for (j in (i+1):num.rbfs) {
				dists[i,j] <- rbfns[[a]]$norm(rbfns[[a]]$c[i,] - rbfns[[a]]$c[j,])
				dists[j,i] <- dists[i,j]
				}
			}
		for (i in 1:num.rbfs) {
			dists[i,i] <- Inf
			rbfns[[a]]$s[i] <- -mean(sort(dists[i,])[1:(num.neighbors)]) / log(tau) 
			}
		
	
		# set the output weights
		if (is.null(old.rbfns)) rbfns[[a]]$w <- matrix(lsars[[a]]$r, num.rbfs, 1) ##think
		else rbfns[[a]]$w <- rbfn.norm.ouput(old.rbfns[[a]], lsars[[a]]$s)
		 }
	
		
	rbfns 
	}
	
		
kbrl.rbfn <- function(rbfns, lsars, df, epsilon = 1e-5, verbose = FALSE, goal.reward = 1, la = FALSE) {
# This implementation uses the analogy with an RBF network pointed out 
# by me in the RGD analysis.
# 'rbfns' is a list with |A| networks; the number and centers of the RBFs might be different
# 'rs2' is a list with |A| entries; each rs2[[i]] is again a list with the rewards 'r' and final states 's''
# 'goal.reward' is the reward associated with the goal. In infinite-horizon tasks it should be NULL
	
	# Compute the design matrices; there is one for each action and each rbfn
	H <- NULL
	for (i in 1:length(lsars)) {
		H1 <-  NULL
		for (j in 1:length(lsars)) {
			H1 <- c(H1, list(rbfn.norm.design.matrix(rbfns[[j]], lsars[[i]]$s2)))
			}
		H <- c(H, list(H1))
		}
		
	A <- NULL
	if (la) {
		for (i in 1:length(lsars)) {
			A <- c(A, list(t(rbfn.norm.design.matrix(rbfns[[i]], lsars[[i]]$s))))
			A[[i]] <- A[[i]] / apply(A[[i]], 1, sum) ### REALLY think... why?!?
			}
		}
		
	goals <- NULL
	if (!is.null(goal.reward)) {
		for (i in 1:length(lsars)) goals <- c(goals, list(lsars[[i]]$r == goal.reward))
		}
	
	num.actions <- length(lsars)
	max.dif <- -Inf
	count <- 0
	while (abs(max.dif) > epsilon) {
		max.dif <- -Inf
		
		for (i in 1:num.actions) {
			Q <- matrix(0, nrow(lsars[[i]]$s2), num.actions)
			for (a in 1:num.actions) {
				Q[,a] <- H[[i]][[a]] %*% rbfns[[a]]$w
				Q[goals[[i]],a] <- 0
				}
				
			old.w <-  rbfns[[i]]$w
			#the new 'w' will be used in the next iteration
			B <- as.matrix(lsars[[i]]$r + df * apply(Q, 1, max)) 
			
			if (la) rbfns[[i]]$w <-  A[[i]] %*% B
			else rbfns[[i]]$w <- B
			max.dif <- max(max.dif, max(abs(old.w - rbfns[[i]]$w)))
			}
			
		if (verbose) print(paste("Max. difference:",max.dif))
		count <- count + 1
		}
		print(count)
	rbfns
	}



kbrl.rbfn2 <- function(rbfns, lsars, df, epsilon = 1e-5, verbose = FALSE, goal.reward = 1, la = FALSE) {
# This implementation uses the analogy with an RBF network pointed out 
# by me in the RGD analysis.
# 'rbfns' is a list with |A| networks; the number and centers of the RBFs might be different
# 'rs2' is a list with |A| entries; each rs2[[i]] is again a list with the rewards 'r' and final states 's''
# 'goal.reward' is the reward associated with the goal. In infinite-horizon tasks it should be NULL
	
	# Compute the design matrices; there is one for each action and each rbfn
	H <- NULL
	R <- NULL
	goals <- NULL
	for (i in 1:length(lsars)) {
		H1 <-  NULL
		A <-t(rbfn.norm.design.matrix(rbfns[[i]], lsars[[i]]$s))
		A <- A / apply(A, 1, sum)
		R <- c(R, list(A %*% lsars[[i]]$r))
		
		if (!is.null(goal.reward)) {
			G <- A %*%(lsars[[i]]$r == goal.reward)
	 		goals <- c(goals, list(G > 0.5))
			}
			
		for (j in 1:length(lsars)) {
			L <- A %*% rbfn.norm.design.matrix(rbfns[[j]], lsars[[i]]$s2)
			H1 <- c(H1, list(L))
			}
		H <- c(H, list(H1))
		}
	
	goals <- NULL
	
	print("Started!")
	num.actions <- length(lsars)
	max.dif <- -Inf
	Q <- matrix(0, nrow(rbfns[[1]]$c), length(rbfns))
	count <- 0
	while (abs(max.dif) > epsilon) {
		max.dif <- -Inf
		for (i in 1:num.actions) {
			old.Q <-  Q
			Q[,i] <- R[[i]] + df * H[[i]][[1]] %*% apply(Q, 1, max)
			Q[goals[[i]],i] <- 0
			}
		max.dif <- max(abs(Q - old.Q))
		if (verbose) print(paste("Max. difference:",max.dif))
		count <- count + 1
		}
	print(count)
	for (i in 1:num.actions) rbfns[[i]]$w <- Q[,i]
	rbfns
	}

			
kbrl.mountain.grid <- function(np = c(5,7,9,11,13,15), size.dataset = NULL, v.real, tau = 0.01, df = 0.995, uniform = TRUE, ...) {
	
	errors <- matrix(0, length(np), 2)
	rbfns <- NULL
	cp <- seq(-1.2, 0.5, length = sqrt(length(v.real)))
	cv <- seq(-0.07, 0.07, length = sqrt(length(v.real)))
	original.S <-  gp(cp, cv)
	
	lsars <- NULL
	means <- NULL
	stdevs <- NULL
	if (!is.null(size.dataset)) {
		cp <- seq(-1.2, 0.5, length = size.dataset)
		cv <- seq(-0.07, 0.07, length = size.dataset)
		sa <- make.sa(gp(cp,cv),3)
		sars <- collect.transitions(sa, mountain.car.transition, c(-1,0,1), normalize = TRUE)
		means <- sars$means
		stdevs <- sars$stdevs
		lsars <- make.lsars(sars)
		}
		
		
	s<- NULL
	for (p in 1:length(np)) {
		# generate rbf network
		print("Generating rbfn...")
		cp <- seq(-1.2, 0.5, length = np[p])
		cv <- seq(-0.07, 0.07, length = np[p])
		s <- gp(cp, cv)
		sa <- make.sa(s, 3)
		sars.rbfn <- collect.transitions(sa, mountain.car.transition, c(-1,0,1), normalize = FALSE)
		if (is.null(size.dataset)) {
			sars.rbfn <- normalize.sars(sars.rbfn)
			means <- sars.rbfn$means
			stdevs <- sars.rbfn$stdevs
			}			
		else sars.rbfn <- normalize.sars(sars.rbfn, means = sars$means, stdevs = sars$stdevs)
		lsars.rbfn <- make.lsars(sars.rbfn)
		rbfns <- make.rbfns.kbrl(lsars.rbfn, tau)
		
		#prepare and run KBRL
		print("Running KBRL...")
		if (is.null(size.dataset)) lsars <- lsars.rbfn
		rbfns <- kbrl.rbfn(rbfns, lsars, df, la = !is.null(size.dataset), ...)
		
		S <- NULL
		S <- normalize(original.S, means = means, stdevs = stdevs)
		
		#compare V generated by KBRL and the "real" one
		print("Computing error...")
		o <- matrix(0, nrow(S), length(rbfns))
		for (i in 1:length(rbfns)) {
			o[,i] <- rbfn.norm.output(rbfns[[i]], S)
			}
		v <- apply(o, 1, max)
		# SSE with respect to the true V
		vm <- matrix(v, sqrt(length(v)), sqrt(length(v)), byrow = TRUE)
		sx <- seq(min(rbfns[[1]]$c[,1]), max(rbfns[[1]]$c[,1]), length = 100)
		sy <- seq(min(rbfns[[1]]$c[,2]), max(rbfns[[1]]$c[,2]), length = 100)
		image(sx,sy, -vm)		
		#text(x = rbfns[[1]]$c[,1], y = rbfns[[1]]$c[,2], labels =  paste(round(rbfns[[1]]$w[,1],d=2)))
		sse <- (v.real - vm)^2
		errors[p,1] <- sum(sse)
		
		# Bellman error
		actions <- apply(o, 1, which.max)
		sars2 <- collect.transitions(list(s = original.S, a = actions), mountain.car.transition, c(-1,0,1), normalize = FALSE)
		sars2$s2 <- normalize(sars2$s2, means = means, stdevs = stdevs)
		for (i in 1:length(rbfns)) {
			o[,i] <- rbfn.norm.output(rbfns[[i]], sars2$s2)
			}
		v2 <- apply(o, 1, max)
		v2[sars2$r == 1] <- 0
		b <- (sars2$r + df * v2 - v)^2
		bm <- matrix(b, sqrt(length(b)), sqrt(length(b)), byrow = TRUE)
		
		
		errors[p,2] <- sum(b)
		print(paste("SSE:", errors[p,1], "Bellman:", errors[p,2]))
		}
		
	errors = errors
	}
	

		
print("kbrl.R loaded")		
		




kbrl <- function(rbfns, rs2, df, epsilon = 1e-5, verbose = FALSE) {
# 	OUT-OF-DATE
# KBRL as described by Ormoneit and Sen; notice that the action-value function of
# the states s' are stores in the network (I used a network to improve efficiency, since
# there's a C function implementing it).
# 'rbfns' is a list with |A| networks; the number and centers of the RBFs migh be different
# rs2 is a list with |A| entries; each rs2[[i]] is again a list with the rewards 'r' and final states 's''
	H <- NULL
	for (i in 1:length(rs2)) {
		H1 <-  NULL
		for (j in 1:length(rs2)) {
			H1 <- c(H1, list(rbfn.norm.design.matrix(rbfns[[j]], rs2[[i]]$s2)))
			}
		H <- c(H, list(H1))
		}
		
	num.actions <- length(rs2)
	max.dif <- -Inf
	while (abs(max.dif) > epsilon) {
		max.dif <- -Inf
		
		for (i in 1:num.actions) {
			Q <- matrix(0, nrow(rs2[[i]]$s2), num.actions)
			R <- matrix(0, nrow(rs2[[i]]$s2), num.actions)
			for (a in 1:num.actions) {
				Q[,a] <- H[[i]][[a]] %*% rbfns[[a]]$w
				R[,a] <- H[[i]][[a]] %*% rs2[[a]]$r
				}
			U <- R + df * Q							
			old.w <-  rbfns[[i]]$w
			#the new 'w' will be used in the next iteration
			rbfns[[i]]$w <- as.matrix(apply(U, 1, max))
			max.dif <- max(max.dif, max(abs(old.w - rbfns[[i]]$w)))
			}
			
		if (verbose) print(paste("Max. difference:",max.dif))
		}
	rbfn <- make.null.rbfn()
	rbfn$c <- rbfns[[1]]$c
	rbfn$s <-  rbfns[[1]]$s
	rbfn$w <-  rbfns[[1]]$w
	for (i in 2:length(rbfns)) rbfn$w <- cbind(rbfn$w, rbfns[[i]]$w)
	rbfn
	}
		
print("detour.samples.generic.R loaded")