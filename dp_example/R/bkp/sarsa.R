 source("backpropagation.R")
source("rbfn.R")

e.greedy <- function(s, Q, epsilon) {
	if (runif(1) > epsilon) a <- which.max	(Q)
	else sample(1:length(Q), 1)
	}



sarsa.2D <- function(x, y, A, transition.function,  lambda = 0, 
							action.selection = e.greedy, epsilon = 0.15, 
							 learning.rate.global = 0.01, df = 0.9, 
							 max.episodes = 2000, graphic = FALSE, 
							 show.cost.to.go = FALSE, show.persp = TRUE,
							 save.Q = TRUE, load.Q = FALSE,
							 file.name.Q = "./files/sarsa.2D.Q.txt", 
							 save.V = TRUE, file.name.V = "./files/sarsa.2D.V.txt", ...) {

	dim.x <- length(x)
	dim.y <- length(y)
	delta.x <- x[2] - x[1]  # considering the grid is uniform
     delta.y <-y[2] - y[1]  # that is, the distance between two points is constant (like an AP)
	lim.x <- c(x[1], x[length(x)])
	lim.y <- c(y[1], y[length(y)])
	
	
	V <- NULL
     goals <- matrix(FALSE, dim.x, dim.y)
     illegals <- matrix(FALSE, dim.x, dim.y)

	Q <- NULL
	if  (load.Q) {
		Q <-as.matrix(read.table(file.name.Q)) 
		dim(Q) <- c(dim.x, dim.y,length(A))
		}
	else 	Q <- array(runif(dim.x * dim.y * length(A), -1e-5, 1e-5), c(dim.x, dim.y,length(A)))
	traces <- array(0, c(dim.x, dim.y,length(A)))

	episode <- 1
	while (episode <= max.episodes) {
		print(paste("Episode",episode))
		valid <- FALSE
		while (!valid) {
			s <- c(runif(1, 1, lim.x[2]) , runif(1, 1, lim.y[2]))
			s[1] <- min(lim.x[2], round((s[1] - lim.x[1]) / delta.x) +1)
			s[2] <- min(lim.y[2], round((s[2] - lim.y[1]) / delta.y) +1)
			valid <- !illegals[s[1],s[2]] && !goals[s[1],s[2]]
			}
	
		a <- action.selection(s, Q[s[1],s[2],], epsilon)
		goal <- FALSE
		illegal <- FALSE

		while (!goal && !illegal) {
			t <- transition.function(s,A[a], ...)
	
			if (!is.null(t$s))  {
				a2 <- action.selection(t$s, Q[t$s[1],t$s[2],], epsilon)
				# TD update
				td.error <-  (t$r + df * Q[t$s[1], t$s[2], a2] - Q[s[1],s[2],a])
				traces[s[1],s[2],a] <- traces[s[1],s[2],a] + 1

				Q <- Q + learning.rate.global * td.error * traces
				traces <- traces * lambda * df
				
				s <- t$s
				s[1] <- min(lim.x[2], round((s[1] - lim.x[1]) / delta.x) +1)
				s[2] <- min(lim.y[2], round((s[2] - lim.y[1]) / delta.y) +1)
		
				a <- a2
				goal <- t$g
				}
			else {
				illegal <- TRUE
				illegals[s[1],s[2]] <- TRUE
				}
			}
		
		if (!illegal) {
			goals[s[1],s[2]] <- TRUE
			traces[,,] <- 0
			
			if (save.Q) {
				Q.save <- matrix(Q,dim(Q)[1], dim(Q)[2]*dim(Q)[3])
				write.table(Q.save, file.name.Q, quote = FALSE, col.names = FALSE, row.names = FALSE)
				}
			
			V <- as.matrix(Q[,,1])
			for (a in 2:length(A)) V <- pmax(V, Q[,,a])
	
			V[goals] <- max(Q) + t$r
			V[illegals] <- min(Q) - t$r
	
			if (graphic) {
				if (show.cost.to.go) V = - V
				if (show.persp) persp(V, theta = 35, phi = 35, shade= 0.5)
				else image(V)
				}

			if (save.V) write.table(V, file.name.V, quote = FALSE, col.names = FALSE, row.names = FALSE)
					
			episode <- episode + 1
			}
		}
      V
	}


sarsa.ap <- function(S, A, transition.function,  num.hidden, td= FALSE,
							lambda = 0, action.selection = e.greedy, 
							lr.hidden = 0.1, lr.output = 0.01, mom = 0.8, 
							epsilon = 0.15, df = 0.9, max.episodes = 2000, 
							graphic = FALSE, show.cost.to.go = FALSE, 
							show.persp = TRUE, save.Q = TRUE, load.Q = FALSE,
							file.name.Q = "./files/sarsa.ap.Q.txt", 
							save.V = TRUE, file.name.V = "./files/sarsa.ap.V.txt", 
							save.at = max.episodes %/% 10,  max.transitions = 1000, ...) {
# if 'action.selection' is a policy, this algorithm becomes the TD algorithm
# in this case, one might want to use only one approximator instead of
# one for each action. This is done by setting up td = TRUE

	dim.S <- ncol(S)
	means.S <- apply(S, 2, mean)
	sds.S <- apply(S, 2, sd)
	sds.S[sds.S == 0] <- 1
	S.normal <- (S - matrix(rep(means.S,nrow(S)),nrow(S),dim.S,byrow=TRUE)) /
	  matrix(rep(sds.S,nrow(S)),nrow(S),dim.S,byrow=TRUE)

	dim.out <- 1
	if (!td) dim.out <- length(A)

     Q <- matrix(0,nrow(S), dim.out)
	
	dim.V <- NULL
	if (graphic) dim.V <- sqrt(nrow(S))
	
	nn <- vector("list", length = dim.out)
	for (a in 1:dim.out) {
		nn[[a]] <- make.nn(S, matrix(0,1,1), num.hidden)
		}

	
	s <- matrix(0,1,dim.S) # it must be a matrix because of the neural net
	q <- array(0,dim.out)
     q.next <- array(0,dim.out)
	
	V <- array(0, nrow(S))
	
	episode <- 1
	while (episode < max.episodes) {
		print(paste("Episode",episode))
		
		ind <- round(runif(1,1,nrow(S)))
		s <- S[ind, , drop = FALSE]	 # it has chosen from S, so that sarsa can be applied to discrete environments
		s.normal <- S.normal[ind, , drop = FALSE]
		for (i in 1:dim.out) q[i] <- nn.output(nn[[i]], s.normal)
      	
		a <- action.selection(s, q, epsilon)

		goal <- FALSE
		illegal <- FALSE

		transitions <- 0
		if (td) max.transitions <- Inf
		while (!goal && !illegal && transitions < max.transitions) {
			t <- transition.function(s, A[a], ...)
			transitions <- transitions + 1
			
			if (!is.null(t$s))  {
                    a2 <- NULL
				target <- NULL
				ts.normal <- NULL                                                                                 

				if (t$g) target <- t$r
				else {
					ts.normal <- (t$s[, , drop = FALSE] - means.S) / sds.S # it has to be normalized because it might not be in S.normal
			
					for (i in 1:dim.out) q.next[i] <- nn.output(nn[[i]], ts.normal) 
					a2 <- action.selection(t$s, q.next, epsilon)
					
					ind <- 1
					if (!td) ind <- a2
					target <- t$r + df * q.next[ind]
					}

				ind <- 1
				if (!td) ind <- a

				nn[[ind]] <- bp.step(nn[[ind]], cbind(1, s.normal), target, lr.hidden = lr.hidden, lr.output = lr.output, mom = mom)$nn

				s <- t$s
				s.normal <- ts.normal
				a <- a2                                                                     
				goal <- t$g
				}
			else 	illegal <- TRUE
			}

		if (!illegal) { # if it was a valid transition
			if (episode %% save.at == 0) {
				for (i in 1:dim.out) Q[,i] <- nn.output(nn[[i]], S.normal)
	
				V <- as.matrix(Q[,1])
				if (!td) {
					for (i in 2:dim.out) V <- pmax(V, Q[,i])
				}
				
				if (graphic) {
					V.graphic <- matrix(V, dim.V, dim.V)
					if (show.cost.to.go) V.graphic = - V.graphic
					if (show.persp) persp(V.graphic, theta = 35, phi = 35, shade= 0.5, ticktype="detailed")
					else image(V.graphic)
					}
	
				if (save.V) write.table(V, file.name.V, quote = FALSE, col.names = FALSE, row.names = FALSE)
				}
					
			episode <- episode + 1
			}
		}
      V
	}


sarsa.pl <- function(S, A, transition.function,  num.hidden, diameter,
							td= FALSE, lambda = 0, action.selection = e.greedy, 
							lr.hidden = 0.1, lr.output = 0.01, mom = 0.8, 
							epsilon = 0.15, df = 0.9, max.episodes = 2000, 
							grid = TRUE,	deltas = NULL, graphic = FALSE, 
							show.cost.to.go = FALSE, show.persp = TRUE,  save.nn = TRUE,
							load.nn  = FALSE, nn.file.name = "./files/sarsa.pl.nn.txt",
							save.V = TRUE, V.file.name = "./files/sarsa.pl.V.txt",
							save.at = max.episodes %/% 10, max.transitions = 1000, ...) {
# if 'action.selection' is a policy, this algorithm becomes the TD algorithm
# in this case, one might want to use only one approximator instead of
# one for each action. This is done by setting up td = TRUE
# grid defines if the states should be fixed in a grid
# deltas is the increment in each dimension of the grid
	dim.S <- ncol(S)
	means.S <- apply(S, 2, mean)
	sds.S <- apply(S, 2, sd)
	sds.S[sds.S == 0] <- 1
	max.S <- apply(S, 2, max)
	min.S <- apply(S, 2, min)

	S.normal <- (S - matrix(rep(means.S,nrow(S)),nrow(S),dim.S,byrow=TRUE)) /
	  matrix(rep(sds.S,nrow(S)),nrow(S),dim.S,byrow=TRUE)

	

	dim.out <- 1
	if (!td) dim.out <- length(A)

     Q <- matrix(0,nrow(S), dim.out)
	
	dim.V <- NULL
	if (graphic) dim.V <- sqrt(nrow(S))
	
	nn <- NULL
	if (load.nn) load(nn.file.name)
	else 	{
		nn <- vector("list", length = diameter)
		for (d in 1:diameter) {
			nn[[d]] <- vector("list", dim.out)
			for (a in 1:dim.out) {
				nn[[d]][[a]] <- make.nn(S, matrix(0,1,1), num.hidden)
				}
			}
		}

	
	s <- matrix(0,1,dim.S) # it must be a matrix because of the neural net
	q <- array(0,dim.out)
     q.next <- array(0,dim.out)
	
	V <- array(0, nrow(S))
	
	episode <- 1
	while (episode <= max.episodes) {
		print(paste("Episode",episode))
		
		ind <- round(runif(1,1,nrow(S)))
		s <- S[ind, , drop = FALSE]	 # it has to be chosen from S, so that sarsa can be applied to discrete environments
		s.normal <- S.normal[ind, , drop = FALSE]

		if (!td) for (i in 1:dim.out) q[i] <- nn.output(nn[[diameter]][[i]], s.normal)
      	
		a <- action.selection(s, q, epsilon)

		goal <- FALSE
		illegal <- FALSE
		transitions   <- 0
		if (td) max.transitions <- Inf
		while (!goal && !illegal && transitions < max.transitions) {
			t <- transition.function(s, A[a], ...)
			transitions <- transitions + 1
			
			if (!is.null(t$s))  {
                    a2 <- NULL
				ts.normal <- NULL                                                                                 

				target <- t$r
				ind <- 1
				if (!td) ind <- a
				nn[[1]][[ind]] <- bp.step(nn[[1]][[ind]], cbind(1, s.normal), target, lr.hidden = lr.hidden, lr.output = lr.output, mom = mom)$nn
				ts.normal <- (t$s[, , drop = FALSE] - means.S) / sds.S # it has to be normalized because it might not be in S.normal
				if (!td) for (i in 1:dim.out) q.next[i] <- nn.output(nn[[diameter]][[i]], ts.normal) 
				a2 <- action.selection(t$s, q.next, epsilon)
				ind2 <- 1
				if (!td) ind2 <- a2
				for (d in 2:diameter) {
					if (!t$g) target <- t$r + df *  nn.output(nn[[d-1]][[ind2]], ts.normal)
					nn[[d]][[ind]] <- bp.step(nn[[d]][[ind]], cbind(1, s.normal), target, lr.hidden = lr.hidden, lr.output = lr.output, mom = mom)$nn
					}

				s <- t$s
				s.normal <- ts.normal
				a <- a2
				goal <- t$g
				}
			else 	illegal <- TRUE
			}
          
		if (!illegal ) { # if it was a valid transition
			if (episode %% save.at == 0) {
				
				if (save.V || graphic) {

					for (i in 1:dim.out) Q[,i] <- nn.output(nn[[diameter]][[i]], S.normal)
		
					V <- as.matrix(Q[,1])
					if (!td) for (i in 2:dim.out) V <- pmax(V, Q[,i])
					
					if (graphic) {
						V.graphic <- matrix(V, dim.V, dim.V)
						if (show.cost.to.go) V.graphic = - V.graphic
						if (show.persp) persp(V.graphic, theta = 35, phi = 35, shade= 0.5, ticktype="detailed")
						else image(V.graphic)
						}
		
					if (save.V) write.table(V, V.file.name, quote = FALSE, col.names = FALSE, row.names = FALSE)
					}
				
				if (save.nn) save(nn, file =  nn.file.name)
				}
				
			episode <- episode + 1
			}
		}


	for (i in 1:dim.out) Q[,i] <- nn.output(nn[[diameter]][[i]], S.normal)
	V <- as.matrix(Q[,1])
	if (!td) for (i in 2:dim.out) V <- pmax(V, Q[,i])
      V
	}


sarsa.plus <- function(S, A, transition.function,  num.hidden, diameter,
							td= FALSE, lambda = 0, action.selection = e.greedy, 
							lr.hidden = 0.1, lr.output = 0.01, mom = 0.8, 
							epsilon = 0.15, df = 0.9, max.episodes = 2000, 
							grid = TRUE,	deltas = NULL, graphic = FALSE, 
							show.cost.to.go = FALSE, show.persp = TRUE,  save.nn = TRUE,
							load.nn  = FALSE, nn.file.name = "./files/sarsa.pl.nn.txt",
							save.V = TRUE, V.file.name = "./files/sarsa.pl.V.txt",
							save.at = max.episodes %/% 10, ...) {
# if 'action.selection' is a policy, this algorithm becomes the TD algorithm
# in this case, one might want to use only one approximator instead of
# one for each action. This is done by setting up td = TRUE
# grid defines if the states should be fixed in a grid
# deltas is the increment in each dimension of the grid
	dim.S <- ncol(S)
	means.S <- apply(S, 2, mean)
	sds.S <- apply(S, 2, sd)
	sds.S[sds.S == 0] <- 1
	max.S <- apply(S, 2, max)
	min.S <- apply(S, 2, min)

	S.normal <- (S - matrix(rep(means.S,nrow(S)),nrow(S),dim.S,byrow=TRUE)) /
	  matrix(rep(sds.S,nrow(S)),nrow(S),dim.S,byrow=TRUE)

	

	dim.out <- 1
#	if (!td) dim.out <- length(A)

     Q <- matrix(0,nrow(S), dim.out)
	
	dim.V <- NULL
	if (graphic) dim.V <- sqrt(nrow(S))
	
	nn <- NULL
	if (load.nn) load(nn.file.name)
	else 	{
		nn <- vector("list", length = diameter)
		for (d in 1:diameter) {
			nn[[d]] <- vector("list", dim.out)
			for (a in 1:dim.out) {
				nn[[d]][[a]] <- make.nn(S, matrix(0,1,1), num.hidden)
				}
			}
		}

	
	s <- matrix(0,1,dim.S) # it must be a matrix because of the neural net
	q <- array(0,dim.out)
     q.next <- array(0,dim.out)
	
	V <- array(0, nrow(S))
	
	episode <- 1
	while (episode <= max.episodes) {
		print(paste("Episode",episode))
		
		ind <- round(runif(1,1,nrow(S)))
		s <- S[ind, , drop = FALSE]	 # it has to be chosen from S, so that sarsa can be applied to discrete environments
		s.normal <- S.normal[ind, , drop = FALSE]

#		if (!td) for (i in 1:dim.out) q[i] <- nn.output(nn[[diameter]][[i]], s.normal)
      	
		a <- action.selection(s, q, epsilon)

		goal <- FALSE
		illegal <- FALSE

		while (!goal && !illegal) {
			t <- transition.function(s, A[a], ...)
		#	if (grid) {
		#		ind <- round((t$s - min.S) / deltas) + 1
		#		ind <- 70 * (ind[1]-1) + ind[2]
		#		t$s <- S[ind,, drop = FALSE]  # REVER; ISSO AQUI TÁ BIDIMENSIONAL...
		#		}
			
			if (!is.null(t$s))  {

				target <- t$r
				ind <- 1
#				if (!td) ind <- a
				bpo <- bp.step(nn[[1]][[ind]], cbind(1, s.normal), target, lr.hidden = lr.hidden, lr.output = lr.output, mom = mom)
				nn[[1]][[ind]] <- bpo$nn
				o1 <- bpo$output
								
				ts.normal <- (t$s[, , drop = FALSE] - means.S) / sds.S # it has to be normalized because it might not be in S.normal
#				if (!td) for (i in 1:dim.out) q.next[i] <- nn.output(nn[[diameter]][[i]], ts.normal) 
				a2 <- action.selection(t$s, q.next, epsilon)
				ind2 <- 1
#				if (!td) ind2 <- a2

				o2 <- 0
				for (d in 2:diameter) {
					o1 <- o1 + nn.output(nn[[d]][[ind]], s.normal)
					if (!t$g) {
						o2 <- o2 + nn.output(nn[[d-1]][[ind2]], ts.normal)
						target <- t$r + df *  o2
						}
					nn[[d]][[ind]] <- bp.step(nn[[d]][[ind]], cbind(1, s.normal), target, lr.hidden = lr.hidden, lr.output = lr.output, mom = mom, o = o1)$nn
					}

				s <- t$s
				s.normal <- ts.normal
				a <- a2
				goal <- t$g
				}
			else 	illegal <- TRUE
			}
          
		if (!illegal ) { # if it was a valid transition
			if (episode %% save.at == 0) {
				
				if (save.V || graphic) {
					Q[,] <- 0
					
					for (d in 1:diameter) {
						for (i in 1:dim.out) Q[,i] <- Q[,i] + nn.output(nn[[d]][[i]], S.normal)
						}
		
					V <- as.matrix(Q[,1])
#					if (!td) for (i in 2:dim.out) V <- pmax(V, Q[,i])
					
					if (graphic) {
						V.graphic <- matrix(V, dim.V, dim.V)
						if (show.cost.to.go) V.graphic = - V.graphic
						if (show.persp) persp(V.graphic, theta = 35, phi = 35, shade= 0.5, ticktype="detailed")
						else image(V.graphic)
						}
		
					if (save.V) write.table(V, V.file.name, quote = FALSE, col.names = FALSE, row.names = FALSE)
					}
				
				if (save.nn) save(nn, file =  nn.file.name)
				}
				
			episode <- episode + 1
			}
		}

	Q[,] <- 0
	
	for (d in 1:diameter) {
		for (i in 1:dim.out) Q[,i] <- Q[,i] + nn.output(nn[[d]][[i]], S.normal)
		}

	V <- as.matrix(Q[,1])
#	if (!td) for (i in 2:dim.out) V <- pmax(V, Q[,i])
      V
	}



sarsa.pl.plus <- function(S, A, transition.function,  num.hidden, diameter,
							td= FALSE, lambda = 0, action.selection = e.greedy, 
							lr.hidden = 0.1, lr.output = 0.01, mom = 0.8, 
							epsilon = 0.15, df = 0.9, max.episodes = 2000, 
							grid = TRUE,	deltas = NULL, graphic = FALSE, 
							show.cost.to.go = FALSE, show.persp = TRUE,  save.nn = TRUE,
							load.nn  = FALSE, nn.file.name = "./files/sarsa.pl.nn.txt",
							save.V = TRUE, V.file.name = "./files/sarsa.pl.V.txt",
							save.at = max.episodes %/% 10, ...) {
# if 'action.selection' is a policy, this algorithm becomes the TD algorithm
# in this case, one might want to use only one approximator instead of
# one for each action. This is done by setting up td = TRUE
# grid defines if the states should be fixed in a grid
# deltas is the increment in each dimension of the grid
	dim.S <- ncol(S)
	means.S <- apply(S, 2, mean)
	sds.S <- apply(S, 2, sd)
	sds.S[sds.S == 0] <- 1
	max.S <- apply(S, 2, max)
	min.S <- apply(S, 2, min)

	S.normal <- (S - matrix(rep(means.S,nrow(S)),nrow(S),dim.S,byrow=TRUE)) /
	  matrix(rep(sds.S,nrow(S)),nrow(S),dim.S,byrow=TRUE)

	

	dim.out <- 1
#	if (!td) dim.out <- length(A)

     Q <- matrix(0,nrow(S), dim.out)
	
	dim.V <- NULL
	if (graphic) dim.V <- sqrt(nrow(S))
	
	nn <- NULL
	if (load.nn) load(nn.file.name)
	else 	{
		nn <- vector("list", length = diameter)
		for (d in 1:diameter) {
			nn[[d]] <- vector("list", dim.out)
			for (a in 1:dim.out) {
				nn[[d]][[a]] <- make.nn(S, matrix(0,1,1), num.hidden)
				}
			}
		}

	
	s <- matrix(0,1,dim.S) # it must be a matrix because of the neural net
	q <- array(0,dim.out)
     q.next <- array(0,dim.out)
	
	V <- array(0, nrow(S))
	
	episode <- 1
	while (episode <= max.episodes) {
		print(paste("Episode",episode))
		
		ind <- round(runif(1,1,nrow(S)))
		s <- S[ind, , drop = FALSE]	 # it has to be chosen from S, so that sarsa can be applied to discrete environments
		s.normal <- S.normal[ind, , drop = FALSE]

#		if (!td) for (i in 1:dim.out) q[i] <- nn.output(nn[[diameter]][[i]], s.normal)
      	
		a <- action.selection(s, q, epsilon)

		goal <- FALSE
		illegal <- FALSE

		while (!goal && !illegal) {
			t <- transition.function(s, A[a], ...)
		#	if (grid) {
		#		ind <- round((t$s - min.S) / deltas) + 1
		#		ind <- 70 * (ind[1]-1) + ind[2]
		#		t$s <- S[ind,, drop = FALSE]  # REVER; ISSO AQUI TÁ BIDIMENSIONAL...
		#		}
			
			if (!is.null(t$s))  {
                    a2 <- NULL
				ts.normal <- NULL                                                                                 

				target <- t$r
				ind <- 1
#				if (!td) ind <- a
				nn[[1]][[ind]] <- bp.step(nn[[1]][[ind]], cbind(1, s.normal), target, lr.hidden = lr.hidden, lr.output = lr.output, mom = mom)$nn
								
				ts.normal <- (t$s[, , drop = FALSE] - means.S) / sds.S # it has to be normalized because it might not be in S.normal
#				if (!td) for (i in 1:dim.out) q.next[i] <- nn.output(nn[[diameter]][[i]], ts.normal) 
				a2 <- action.selection(t$s, q.next, epsilon)
				ind2 <- 1
#				if (!td) ind2 <- a2
				for (d in 2:diameter) {
					o1 <- 0
					for (i in 1:d) o1 <- o1 + nn.output(nn[[i]][[ind]], s.normal)
					if (!t$g) {
						o2 <- 0
						for (i in 1:(d-1)) o2 <- o2 + nn.output(nn[[i]][[ind2]], ts.normal)
						target <- t$r + df *  o2
						}
					nn[[d]][[ind]] <- bp.step(nn[[d]][[ind]], cbind(1, s.normal), target, lr.hidden = lr.hidden, lr.output = lr.output, mom = mom, o = o1)$nn
					}

				s <- t$s
				s.normal <- ts.normal
				a <- a2
				goal <- t$g
				}
			else 	illegal <- TRUE
			}
          
		if (!illegal ) { # if it was a valid transition
			if (episode %% save.at == 0) {
				
				if (save.V || graphic) {
					Q[,] <- 0
					
					for (d in 1:diameter) {
						for (i in 1:dim.out) Q[,i] <- Q[,i] + nn.output(nn[[d]][[i]], S.normal)
						}
		
					V <- as.matrix(Q[,1])
#					if (!td) for (i in 2:dim.out) V <- pmax(V, Q[,i])
					
					if (graphic) {
						V.graphic <- matrix(V, dim.V, dim.V)
						if (show.cost.to.go) V.graphic = - V.graphic
						if (show.persp) persp(V.graphic, theta = 35, phi = 35, shade= 0.5, ticktype="detailed")
						else image(V.graphic)
						}
		
					if (save.V) write.table(V, V.file.name, quote = FALSE, col.names = FALSE, row.names = FALSE)
					}
				
				if (save.nn) save(nn, file =  nn.file.name)
				}
				
			episode <- episode + 1
			}
		}

	Q[,] <- 0
	
	for (d in 1:diameter) {
		for (i in 1:dim.out) Q[,i] <- Q[,i] + nn.output(nn[[d]][[i]], S.normal)
		}

	V <- as.matrix(Q[,1])
#	if (!td) for (i in 2:dim.out) V <- pmax(V, Q[,i])
      V
	}



sarsa.pl.rbfn <- function(S, A, transition.function,  num.hidden, diameter,
							td= FALSE, lambda = 0, action.selection = e.greedy, 
							lr.hidden = 0.1, lr.output = 0.01, mom = 0.8, 
							epsilon = 0.15, df = 0.9, max.episodes = 2000, 
							graphic = FALSE, show.cost.to.go = FALSE, 
							show.persp = TRUE, save.Q = TRUE, load.Q = FALSE,
							file.name.Q = "./files/sarsa.ap.Q.txt", 
							save.V = TRUE, file.name.V = "./files/sarsa.ap.V.txt", ...) {
# if 'action.selection' is a policy, this algorithm becomes the TD algorithm
# in this case, one might want to use only one approximator instead of
# one for each action. This is done by setting up td = TRUE

	dim.S <- ncol(S)
	means.S <- apply(S, 2, mean)
	sds.S <- apply(S, 2, sd)
	sds.S[sds.S == 0] <- 1
	S.normal <- (S - matrix(rep(means.S,nrow(S)),nrow(S),dim.S,byrow=TRUE)) /
	  matrix(rep(sds.S,nrow(S)),nrow(S),dim.S,byrow=TRUE)

	dim.out <- 1
	if (!td) dim.out <- length(A)

     Q <- matrix(0,nrow(S), dim.out)
	
	dim.V <- NULL
	if (graphic) dim.V <- sqrt(nrow(S))
	
	nn <- vector("list", length = diameter)
	for (d in 1:diameter) {
		nn[[d]] <- vector("list", dim.out)
		for (a in 1:dim.out) {
			nn[[d]][[a]] <- make.rbfn(S.normal, matrix(0,1,1), num.hidden)
			}
		}

	
	s <- matrix(0,1,dim.S) # it must be a matrix because of the neural net
	q <- array(0,dim.out)
     q.next <- array(0,dim.out)
	
	V <- array(0, nrow(S))
	
	episode <- 1
	while (episode < max.episodes) {
		print(paste("Episode",episode))
		
		ind <- round(runif(1,1,nrow(S)))
		s <- S[ind, , drop = FALSE]	 # it has chosen from S, so that sarsa can be applied to discrete environments
		s.normal <- S.normal[ind, , drop = FALSE]
		for (i in 1:dim.out) q[i] <- rbfn.output(nn[[diameter]][[i]], s.normal)
		
		a <- action.selection(s, q, epsilon)

		goal <- FALSE
		illegal <- FALSE

		while (!goal && !illegal) {
			t <- transition.function(s, A[a], ...)
			
			if (!is.null(t$s))  {
                    a2 <- NULL
				target <- NULL
				ts.normal <- NULL                                                                                 

				target <- t$r
				ind <- 1
				if (!td) ind <- a
				 nn[[1]][[ind]] <- step.bp.gaussian(nn[[1]][[ind]], s.normal, target)
				
				ts.normal <- (t$s[, , drop = FALSE] - means.S) / sds.S # it has to be normalized because it might not be in S.normal
				for (i in 1:dim.out) q.next[i] <- rbfn.output(nn[[diameter]][[i]], ts.normal) 
				a2 <- action.selection(t$s, q.next, epsilon)
				ind2 <- 1
				if (!td) ind2 <- a2

				for (d in 2:diameter) {
					target <- t$r + df * rbfn.output(nn[[d-1]][[ind2]], ts.normal) 
                    	nn[[d]][[ind]] <- step.bp.gaussian(nn[[d]][[ind]], s.normal, target)
					}

				s <- t$s
				s.normal <- ts.normal
				a <- a2
				goal <- t$g
				}
			else 	illegal <- TRUE
			}

		if (!illegal) { # if it was a valid transition
			for (i in 1:dim.out) Q[,i] <- rbfn.output(nn[[1]][[i]], S.normal)

			V <- as.matrix(Q[,1])
			if (!td) {
				for (i in 2:dim.out) V <- pmax(V, Q[,i])
			}
			
			if (graphic) {
				V.graphic <- matrix(V, dim.V, dim.V)
				if (show.cost.to.go) V.graphic = - V.graphic
 				if (show.persp) persp(V.graphic, theta = 35, phi = 35, shade= 0.5, ticktype="detailed")
				else image(V.graphic)
				}

			if (save.V) write.table(V, file.name.V, quote = FALSE, col.names = FALSE, row.names = FALSE)
					
			episode <- episode + 1
			}
		}
      V
	}

print("sarsa.R loaded")
