
#library("AMORE")

value.iteration.2D <- function(x, y,  A, transition.function,  df
                        = 0.9, epsilon = 1e-6, graphic = FALSE,
                        show.cost.to.go = FALSE, show.persp = TRUE,
                        save.function = TRUE, load.function = FALSE, file.name =
                        "./files/value.iteration.2D.output.txt", 
                        max.iterations = 10000, ...) {
# Only works in deterministic 2D environments
# x is the x axis
# y is the y axis
# A is the action space
# epsilon is the maximum difference on the value function between
# two consecutive iterations
# df is the discount factor
	dim.x <- length(x)
	dim.y <- length(y)
	delta.x <- x[2] - x[1]  # considering the grid is uniform
   delta.y <- y[2] - y[1]  # that is, the distance between two points is 
                             # constant (like an AP)

	V <- matrix(0, dim.x, dim.y)
	if  (load.function) V <-as.matrix(read.table(file.name))
	
	# build the matrix with next states
	S.next <- array(0, c(dim.x, dim.y, length(A), 2)) # keep the INDICES of the next states, not the states themselves
	r.next <- array(0, c(dim.x, dim.y, length(A)))
	g.next <- array(FALSE , c(dim.x, dim.y, length(A)))
	
	# keep track of goal and illegal states
	goal <- matrix(FALSE, dim.x, dim.y)
	illegal <- matrix(FALSE, dim.x, dim.y) #like walls in mazes
	
	nst <- 0
	for (i in 1:dim.x) {
		for (j in 1:dim.y) {
				for (a in 1:length(A)) {
					t <- transition.function(c(x[i],y[j]),A[a], ...)
					if (is.null(t$s)) illegal[i,j] <- TRUE #redundant, but that's ok; it only happens once
					else { 
						S.next[i,j,a,1] <- round((t$s[1] - x[1]) / delta.x) +1
						S.next[i,j,a,2] <- round((t$s[2] - y[1]) / delta.y) +1
						g.next[i,j,a] <- t$g
						if (S.next[i,j,a,1] == i && S.next[i,j,a,2] == j && !t$g) nst <- nst + 1
						r.next[i,j,a] <- t$r
						if (g.next[i,j,a]) goal[S.next[i,j,a,1], S.next[i,j,a,2]] <- TRUE #redundant again...
						}
					# ADD PROBABILITY LATER
					}
				}
		}
	
	print(paste("Number of suspicious self-transitions:", nst))
	
	avg.r <- abs(mean(r.next)) # keep the average reward
	
	tolerance <- epsilon * (1 - df) / (2 * df)
	max.dif <- Inf
	iteration <- 0
	while (max.dif > tolerance && iteration < max.iterations) {
		print(paste("Iteration ", iteration + 1))
	
		max.dif <- -Inf
		for (i in 1:dim.x) {
			for (j in 1:dim.y) {
				if (!illegal[i,j]) {
					if (!goal[i,j]) {
						v.max <- -Inf
						for (a in 1:length(A)) {
							v <- r.next[i,j,a] + df * V[S.next[i,j,a,1], S.next[i,j,a,2]]
							if (v >  v.max) v.max <- v
							}
						if (abs(V[i,j] - v.max) > max.dif) max.dif <- abs(V[i,j] - v.max)
						V[i,j] <- v.max
						}
					else V[i,j] <- 0
					}
			}
		}

		V.complete <- V
		V.complete[goal] <- max(V.complete) + avg.r
		V.complete[illegal] <- min(V.complete) - avg.r

		
		if (graphic) {
			if (show.cost.to.go) V.complete = - V.complete
			if (show.persp) persp(V.complete, theta = 35, phi = 35, shade= 0.5, ticktype = "detailed")
			else image(V.complete)
			}
		if (save.function) write.table(V.complete, file.name, col.names = FALSE, 
							     row.names = FALSE, quote = FALSE)
		iteration <- iteration + 1
		}
	
		V[goal] <- max(V) + avg.r
		V[illegal] <- min(V) - avg.r
		V
	}
	

value.iteration.2D.Q <- function(x, y,  A, transition.function,  df
                        = 0.9, epsilon = 1e-6, max.iterations = 10000, 
                        border.penalty = 0, ...) {
# The way it is, it only works in deterministic 2D environments
# returns Q instead of V
# x is the x axis
# y is the y axis
# A is the action space
# epsilon is the maximum difference on the value function between
# two consecutive iterations
# df is the discount factor
	dim.x <- length(x)
	dim.y <- length(y)
	delta.x <- x[2] - x[1]  # considering the grid is uniform
     delta.y <- y[2] - y[1]  # that is, the distance between two points is 
                             # constant (like an AP)

	Q <- array(0, c(dim.x, dim.y, length(A)))
	#if  (load.function) V <-as.matrix(read.table(file.name))
	
	# build the matrix with next states
     S.next <- array(0, c(dim.x, dim.y, length(A), 2)) # keep the INDICES of
                                   # the next states (not the states themselves)
	r.next <- array(0, c(dim.x, dim.y, length(A)))
	g.next <- array(FALSE , c(dim.x, dim.y, length(A)))
	
	# keep track of goal and illegal states
	goal <- matrix(FALSE, dim.x, dim.y)
	illegal <- matrix(FALSE, dim.x, dim.y) #like walls in mazes
	
     
     nst <- 0
	for (i in 1:dim.x) {
		for (j in 1:dim.y) {
				for (a in 1:length(A)) {
					t <- transition.function(c(x[i],y[j]),A[a], ...)
                         if (is.null(t$s)) illegal[i,j] <- TRUE #redundant
					else { 
						S.next[i,j,a,1] <-round((t$s[1]-x[1])/delta.x)+1
                              border <- FALSE                        
                              if (S.next[i,j,a,1] < 1) {
                                 S.next[i,j,a,1] <- 1
                                 border <- TRUE
                                 }
                              else if (S.next[i,j,a,1] > dim.x) {
                                 S.next[i,j,a,1] <- dim.x
                                 border <- TRUE
                                 }
                                 
                          	S.next[i,j,a,2] <-round((t$s[2]-y[1])/delta.y)+1
                              if (S.next[i,j,a,2] < 1) {
                                 S.next[i,j,a,2] <- 1
                                 border <- TRUE
                                 }
                              else if (S.next[i,j,a,2] > dim.y) {
                                 S.next[i,j,a,2] <- dim.y
                                 border <- TRUE
                                 }
                              
                         	g.next[i,j,a] <- t$g
                              if (border) g.next[i,j,a] <- g.next[i,j,a] +
                                                            border.penalty
                              
						r.next[i,j,a] <- t$r
                              if (S.next[i,j,a,1] == i && S.next[i,j,a,2] ==
                                  j && !t$g) nst <- nst + 1
						if (g.next[i,j,a]) goal[S.next[i,j,a,1],
                                    S.next[i,j,a,2]] <- TRUE #redundant again...
						}
					# ADD PROBABILITY LATER
					}
				}
		}
	
	print(paste("Number of suspicious self-transitions:", nst))
	
	avg.r <- abs(mean(r.next)) # keep the average reward
	
	
	tolerance <- epsilon * (1 - df) / (2 * df)
	max.dif <- Inf
	iteration <- 0
	while (max.dif > tolerance && iteration < max.iterations) {
		print(paste("Iteration ", iteration + 1))
		max.dif <- -Inf
		for (i in 1:dim.x) {
			for (j in 1:dim.y) {
				if (!illegal[i,j]) {
					if (!goal[i,j]) {
						for (a in 1:length(A)) {
                                   V.next <- -Inf
                                   for (a2 in 1:length(A)) {
                                     V.next <- max(V.next, Q[S.next[i,j,a,1],
                                                            S.next[i,j,a,2],a2])
                                     }
							Q.old <- Q[i,j,a]
							Q[i,j,a] <- r.next[i,j,a] + df * V.next
							dif <- abs(Q[i,j,a] - Q.old)
                                   if (dif > max.dif) max.dif <- dif
							}
						}
					else Q[i,j,] <- 0
					}
                  }
            }

		iteration <- iteration + 1
          image(Q.to.V(Q))        
		}
	Q
	}
	

Q.to.V <- function(Q) {
   V <- Q[,,1]
   for (a in 2:dim(Q)[3]) {
      V <- pmax(V, Q[,,a])
      }
   V
   }

Q.to.p <- function(Q) {
   p <- matrix(0,dim(Q)[1], dim(Q)[2])
   for (i in 1:dim(Q)[1]) {
      ## for pole
      p[i,] <- apply(Q[i,,], 1, which.max)
      p[i,Q[i,,1] == Q[i,,2]] <- 2
      p[i,Q[i,,2] == Q[i,,3]] <- 2
      }
   p
   }
      

Q2D.to.Q1D <- function(Q2) {
   Q1 <- matrix(0, dim(Q2)[1] * dim(Q2)[2], dim(Q2)[3])
   for (a in 1:dim(Q2)[3]) {
      Q1[,a] <- array(t(Q2[,,a]))
      }
   Q1
   }  
    



## CHECK FROM THIS POINT DOWN
value.iteration.ap <- function(S, A, transition.function, num.hidden, 
										num.epochs = 2000, tolerance = 0.01, df = 0.9, 
										graphic = FALSE, show.cost.to.go = FALSE,
										show.persp = TRUE,  save.function = TRUE, load.function = FALSE,
										file.name = "./files/ap.value.iteration.output.txt",
										learning.rate.global = 0.01, momentum.global = 0.8,
										hidden.layer = "sigmoid", method="ADAPTgdwm", ...) {
# S is the state space
# A is the action space
# num.hidden is the number of hidden units to be used in the approximator
# tolerance is the maximum difference on the value function between two consecutive iterations
# df is the discount factor

# NORMALIZAR A ENTRADA!!!

	V <- NULL
	if  (load.function) V <-as.matrix( read.table(file.name))
	else V <- matrix(0, nrow(S), 1)
	dim <- sqrt(length(V))
			

	nn <- newff(n.neurons = c(ncol(S), num.hidden, 1), learning.rate.global = learning.rate.global, 
			momentum.global = momentum.global ,  error.criterium="LMS",  Stao=NA, 
			hidden.layer = hidden.layer,  method =method)
	
	# build the matrix with next states
	S.next <- array(0, c(nrow(S), length(A), ncol(S)))
	V.next <- matrix(0, nrow(S), length(A))
	r.next <- matrix(0, nrow(S), length(A))
	g.next <- matrix(FALSE , nrow(S), length(A))
	goal <- array(FALSE, nrow(S))
	illegal <- array(FALSE, nrow(S))	

	for (s in 1:nrow(S)) { # populate the matrices with transitions
			for (a in 1:length(A)) {
				t <- transition.function(S[s,, drop= FALSE], A[a], ...)
				if (is.null(t$s)) illegal[s] <- TRUE # redundant
				else {
					S.next[s,a,] <- t$s
					g.next[s,a] <- t$g
					r.next[s,a] <- t$r

					if (identical(t$s, S[s,,drop = FALSE]) && t$g) goal[s] <- TRUE # absorbing state
					}
				# ADD PROBABILITY LATER
				}
			}

	# remove the illegal states, like walls in mazes	
	V.original <- V
	
	V <- V[!illegal,, drop = FALSE]
	S <- S[!illegal, , drop = FALSE]
	S.next <- S.next[!illegal , , ,  drop = FALSE]
     V.next <- V.next[!illegal , , drop = FALSE]
	r.next <- r.next[!illegal , , drop = FALSE]
	g.next <- g.next[!illegal , , drop = FALSE]
     goal <- goal[!illegal]
	avg.r <- abs(mean(r.next))

	# normalize
	no <- normalize(S, return.params = TRUE)
	S.normal <- no$data
	S.next.normal <- S.next
	for (a in 1:length(A)) S.next.normal[,a,] <- normalize(S.next[,a,], no$means, no$stdevs)

	max.dif <- tolerance + 1
	iteration <- 1
	while (max.dif > tolerance) {
		print(paste("Iteration ", iteration))
		nn <- train(nn, S.normal, V, error.criterium="LMS", report=FALSE, show.step=num.epochs / 10, n.shows=10)$net
     	
		for (a in 1:length(A)) {
			V.next[, a] <- sim.MLPnet(nn, S.next.normal[,a,])
			V.next[goal,a] <- 0				
			}

		V.tmp <- r.next + df * V.next
	
		V.new <- V.tmp[,1]
		for (a in 2:length(A)) V.new <- pmax(V.new, V.tmp[,a])
	
		max.dif <- max(abs(V.new - V))
		
		V <- V.new

		V.tmp <- V # I could use V directly, but it is dangerous, because it would depend on a posterior assignment 
		V.tmp[goal] <- max(V.tmp) + avg.r 
		V.complete <- V.original
		V.complete[!illegal,] <- V.tmp
		V.complete[illegal,] <- min(V.complete) - avg.r
         
		if (graphic) {
			c <- 1
			 if (show.cost.to.go) c <- -1
			M <- matrix(c * V.complete, dim, dim, byrow = TRUE)
			if (show.persp) persp(M, phi = 35, theta = 35, main = paste("Iteration ", iteration), shade = 0.5, ticktype = "detailed")
			else image(M)
			}
		if (save.function) write.table(V.complete, file.name, col.names = FALSE, 
							     row.names = FALSE, quote = FALSE)
		iteration <- iteration + 1
		}
	V[goal] <- max(V) + avg.r
	V.original[!illegal,] <- V
	V.original[illegal,] <- min(V) - avg.r
	V.original
	}

			
generate.policy.2D <- function(x, y, V, A, transition.function, delta.x = x[2] - x[1], delta.y = y[2] - y[1], ...) {
	
	policy <- matrix(0, length(x), length(y))
	
	delta.x <- x[2] - x[1]  # considering the grid is uniform
     delta.y <-y[2] - y[1]  # that is, the distance between two points is constant (like an AP)
	lim.x <- c(x[1], x[length(x)])
	lim.y <- c(y[1], y[length(y)])
	for (i in 1:length(x)) {
		for (j in 1:length(y)) {
			s <- c(x[i],y[j])
			best.a <- 0 #Any
			V.max <- -Inf
			ind <- c(0,0)
			problem <- TRUE
			for (a in 1:length(A)) {
				s.next <- transition.function(s, A[a], ...)$s
			#	print(paste("[",i,",",j,"]", "[",ind[1],",",ind[2],"]"))
				if (!is.null(s.next)) {
					ind[1] <- round((s.next[1] - lim.x[1]) / delta.x) +1
					ind[2] <- round((s.next[2] - lim.y[1]) / delta.y) +1
					if (ind[1] != i || ind[2] != j) problem <- FALSE
					if (V[ind[1], ind[2]] >= V.max) {
						V.max <- V[ind[1], ind[2]]
						best.a <- a
						}
					}
				}
			if (problem) {
				print("Warning: self-transition. If nota a terminal state, the discretization is probably too coarse.")
				print(s)
				}
			policy[i,j] <- best.a
			}
		}
	policy
	}
		    





print("dp.2D.R  loaded")