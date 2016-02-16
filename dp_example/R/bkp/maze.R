

maze.transition <- function(s, a, G, W = matrix(FALSE, nrow(G), ncol(G)), rr = 0, gr = 1, wr = 0) {
# G is a matrix with the goal states
# W is a matrix with the walls
# rr is the regular reward, that is, the reward associated with a normal transition
# gr is the goal reward
#wr is the reward delivered to the agent when it runs into a wall
	g <- FALSE
	r <- wr
	p <- 1
	if (W[s[1],s[2]]) { # if s is invalid
		s <- NULL
		r <- 0
		p <- 0
		}
	else {
		if (G[s[1],s[2]]) { #if the current state is a goal
			r <- 0
			g <- TRUE
			p <- 1
			}
		else { # only here a real transition will happen
			if (a == "N" && s[2] < ncol(G) && !W[s[1], s[2]+1]) {
				s[2] <- s[2] +1
				r <- rr
				}
			else if (a == "S" && s[2] > 1 && !W[s[1], s[2]-1]) {
					s[2] <- s[2] -1
					r <- rr
					}
				else if (a == "E" && s[1] < ncol(G) && !W[s[1]+1, s[2]]) {
							s[1] <- s[1] +1
							r <- rr
							}
						else if (a == "W" && s[1] > 1 && !W[s[1]-1, s[2]]) {
								s[1] <- s[1] -1
								r <- rr
								}
			if (G[s[1], s[2]]) { # here s is already the next state
				r <- gr
				g <- TRUE
				}
			}
		}
		
	list(s = s, r = r, g = g,  p =p) 
	}


continuous.maze.transition <- function(s, a, G, W = NULL, rr = 0, gr = 1, wr = -1, step = 0.02) {
# G is a matrix with the goal locations
# W is a matrix with the penalty regions; it has both their location and radius (ncol(W) = 3)
	dist <- function(s, g) sqrt(sum((s - g)^2))
	
	g <- FALSE
	r <- rr
	p <- 1

	i <- 1
	while (i <= nrow(G) && !g) { # check if the current state is one of the goals
		if (dist(s, G[i,1:2]) < step) {
			r <- 0
			g <- TRUE
			p <- 1 
			}
		i <- i +1
		}

	if (!g) { # only here a real transition will happen
		if (a == "N") s[2] <- s[2] +step
		else if (a == "S") 	s[2] <- s[2] -step
			else if (a == "E") 	s[1] <- s[1] +step
				else if (a == "W") s[1] <- s[1] -step
		
		
		i <- 1
		while (i <= nrow(G) && !g) { # check if the current state is one of the goals
			if (dist(s, G[i,1:2]) < step) {
				r <- gr
				g <- TRUE
				p <- 1 
				}
			i <- i +1
			}
	
		if (!g) { 
			i <- 1
			penalty <- FALSE
			while (i <= nrow(W) && !penalty) { # check if the current state is in one of penalty regions
				if (dist(s, W[i,1:2]) < W[i,3] * step) {
					r <- wr
					penalty <- TRUE
					p <- 1 
					}
				i <- i +1
				}
			}

		s[1] <- min(max(s[1], 0), 1)
		s[2] <- min(max(s[2], 0), 1)
		}

	list(s = s, r = r, g = g,  p =p) 
	}



puddle.maze.transition <- function(s, a, G, P = NULL, rr = -1, gr = 0, pr = -400, step = 0.05) {
# G is a matrix with the goal locations and radius (dim(G) = num.goals x 3)
# P is a matrix with the "puddles",  each one represented by its two extremes and its radius (dim(P) = num.puddles x 5)
	
	dist <- function(s, g) sqrt(sum((s - g)^2))
	coef <- function(p1, p2, p3)  { ((p3[1] - p1[1]) * (p2[1] - p1[1]) + (p3[2] - p1[2]) * (p2[2] - p1[2])) / sum((p2 - p1)^2) }
	
	g <- FALSE
	r <- rr
	p <- 1

	i <- 1
	while (i <= nrow(G) && !g) { # check if the current state is one of the goals
		if (dist(s, G[i,1:2]) < G[i,3]) {
			r <- 0
			g <- TRUE
			p <- 1 
			}
		i <- i +1
		}

	if (!g) { # only here a real transition will happen
		if (a == "N") s[2] <- s[2] +step
		else if (a == "S") 	s[2] <- s[2] -step
			else if (a == "E") 	s[1] <- s[1] +step
				else if (a == "W") s[1] <- s[1] -step
		
		
		i <- 1
		while (i <= nrow(G) && !g) { # check if the current state is one of the goals
			if (dist(s, G[i,1:2]) < G[i,3]) {
				r <- 0
				g <- TRUE
				p <- 1 
				}
			i <- i +1
			}
	
		max.dist <- -Inf
		if (!g && !is.null(P)) { # check the puddles
			for (i in 1:nrow(P)) {
				p1 <- P[i,1:2]
				p2 <- P[i,3:4]
				u <- coef(p1, p2, s)
				if (u >= 0 && u <= 1) { # the projection is in between the extremes
					p <- c(p1[1] + u * (p2[1] - p1[1]), p1[2] + u *(p2[2] - p1[2])) 
					d <- dist(s, p)
					if (d < P[i, 5]) {
						d <- P[i,5] - d # the distance to the edge
						if (d > max.dist) max.dist <- d
						}
					}
				}
			}
		
		
		if (max.dist > - Inf) r <- r + pr * max.dist

		s[1] <- min(max(s[1], 0), 1)
		s[2] <- min(max(s[2], 0), 1)
		}

	list(s = s, r = r, g = g,  p =p) 
	}




maze.test <- function(G, W = matrix(FALSE, nrow(G), ncol(G)), rr = 0, gr = 1, wr = 0) {
	s <- c(sample(1:nrow(G),1), sample(1:ncol(G),1))
	a <- c("N", "S", "E", "W")
	goal <- FALSE
	while(!goal) {
		print(s)
		a <- sample(a)
		t <- maze.transition(s, a[1], G, W, rr, gr, wr)
		s <- t$s
		goal <- t$g
		}
	print(s)
	}

continuous.maze.test <- function(G, P, rr = -1, gr = 0, pr = -400, step = 0.05) {
	s <- array(runif(2))
	a <- c("N", "S", "E", "W")
	goal <- FALSE
	while(!goal) {
		print(s)
		a <- sample(a)
		t <- continuous.maze.transition(s, a[1], G, P, rr, gr, pr, step)
		s <- t$s
		goal <- t$g
		}
	print(s)
	}


maze.generate.policy <- function(V, A=c("N","S","E","W"),
								transition.function=maze.transition, ...) {
	dim.x <- nrow(V)
	dim.y <- ncol(V)
	policy <- matrix(0, dim.x, dim.y)
	for (i in 1:dim.x) {
		for (j in 1:dim.y) {
			s <- c(i,j)
			best.a <- 0 #Any
			V.max <- -Inf
			for (a in 1:length(A)) {
				s.next <- transition.function(s, A[a], ...)$s
				if (!is.null(s.next)) {
					if (V[s.next[1], s.next[2]] > V.max) {
						V.max <- V[s.next[1], s.next[2]]
						best.a <- a
						}
					}
				}
			policy[i,j] <- best.a
			}
		}
	policy
	}

maze.random.start <- function(G, ...) {
	s <-  NULL
	while (is.null(s)) {
		s <- array(0,2)
		s[1] <- round(runif(1, 0.51, nrow(G) + 0.49))
		s[2] <- round(runif(1, 0.51, ncol(G) + 0.49))
		}
	s
	}


continuous.maze.op <- function(s, Q = NULL, epsilon = NULL,  delta =0.05, lim.max = 21, policy = ps.puddle) {
	ind1 <- max(1,min(round(s[1] / delta) + 1, lim.max))
     ind2 <- max(1,min(round(s[2] / delta)  + 1, lim.max))
	
	policy[ind1, ind2]
	}

print("maze.R loaded")

                                 
