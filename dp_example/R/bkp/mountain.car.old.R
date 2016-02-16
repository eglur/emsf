

mountain.car.transition <- function(s, a,  min.vel = -0.07, max.vel = 0.07,
   min.pos = -1.2 , max.pos = 0.5, regular.reward = 0, goal.reward = 1,
   sd = 0) {
# s[1] is the car position                                                       
# s[2] is the car velocity
# a is the action: -1, 1 or 0 
     r <- regular.reward # standard reward
     g <- FALSE # in principle, the next state is not a goal
     p <- 1

	if (s[1] > max.pos) { # if the current s is a goal no changes will be made to s
		s[1] <- max.pos
		r <- 0 
		g <- TRUE 
		p <- 1
		}
	else {
      a <- a + rnorm(1,0,sd)
      
		s[2] <- s[2] +  0.001 * a - 0.0025 * cos(3 * s[1])
		s[1] <- s[1] + s[2]
		
		if (s[2] > max.vel) s[2] <-  max.vel
		else if (s[2] < min.vel) s[2] <-  min.vel
			
		if (s[1] < min.pos) {
			s[1] <-  min.pos;
			s[2] <-  0;
			}
		else if (s[1] > max.pos) { # s here is already the next state
			s[1] <- max.pos
			r <- goal.reward;
			g <- TRUE
			p <- 1
			}
		}
	
	list(s = s, r = r, g = g, p =p) # p is the probability of the transition; g indicates if the state is one of the goals
	}


mountain.test.policy <- function(Q, S, df, 
                               max.steps=300, 
                               min.vel = -0.07, 
                               max.vel = 0.07,
                               min.pos = -1.2 , 
                               max.pos = 0.5,...) {
# Q is n x n x |A|, where n = sqrt(|S|)
   A <- c(-1,0,1)
   deltax <- (max.pos-min.pos)/ (nrow(Q) - 1)
   deltay <- (max.vel-min.vel)/ (ncol(Q) - 1)
   
   res.return <- array(0, nrow(S))
   res.steps  <- array(0, nrow(S))  
   
   for (i in 1:nrow(S)) {
      goal <- FALSE
      tr <- 0
      num.steps <- 0
      s <- S[i,]
      ddf <- 1
      
#       plot(0,0,xlim=c(0,1),ylim=c(0,1))      
 
      while(!goal && num.steps < max.steps) {
         indx <- min(max((s[1] - min.pos) %/% deltax + 1, 1), nrow(Q))
         indy <- min(max((s[2] - min.vel) %/% deltay + 1, 1), ncol(Q))
         a <- which.max(Q[indx, indy,])
         t <- mountain.car.transition(s, A[a], min.vel = min.vel,
                     max.vel = max.vel, min.pos = min.pos, max.pos=max.pos,...)
         s <- t$s
         goal <- t$g
         num.steps <- num.steps + 1
         tr <- tr + ddf * t$r
         ddf <- ddf * df
#          points(s[1],s[2])
         }
      res.return[i] <- tr
      res.steps[i]  <- num.steps
      }
   list(st = res.steps, tr = res.return)
   }

     

mountain.car.transition.all <- function(S, A = array(c(-1,0,1)),  
		min.vel = -0.07, max.vel = 0.07, min.pos = -1.2 , max.pos = 0.5) {
# S[,1] are the car positions
# S[,2] are the car velocities
# A is the action space {-1, 1, 0}
# the return value is a c(nrow(S),2*length(A)), in which each row corresponds 
# to the possible successor states of each original one
	
	vels <- matrix(rep(S[,2], length(A)), nrow(S), length(A))
	pos <- matrix(rep(S[,1], length(A)), nrow(S), length(A))
	
	vels <- vels + outer(-0.0025 * cos(3 * S[,1]),  0.001 * A, "+")
	pos <- pos + vels

	vels[vels > max.vel] <- max.vel
	vels[vels < min.vel] <- min.vel

	ind <- pos < min.pos
     pos[ind] <- min.pos
	vels[ind] <- 0
	
	goals <- matrix(pos > max.pos, nrow(S), length(A))
	
	rewards <- matrix(-1, nrow(S), length(A))
	rewards[goals] <- 0
	
	S <- cbind(pos[,1], vels[,1])
	for (i in 2:length(A)) S <- cbind(S, pos[,i], vels[,i])

	list( S = S, rewards = rewards, goals = goals)
	}

 

mountain.car.test <- function(min.vel = -0.07, max.vel = 0.07, min.pos = -1.2 , max.pos = 0.5) {
	s <- c(runif(1,min.vel, max.vel), runif(1,min.pos, max.pos))
	a <- c(-1,0,1)
	goal <- FALSE
	while(!goal) {
		print(s)
		a <- sample(a)
		t <- mountain.car.transition(s,a[1],min.vel, max.vel, min.pos,max.pos)
		s <- t$s
		goal <- t$g
		}
	}


mountain.car.op <- function(s, Q = NULL, epsilon = NULL,  policy = p.mountain, 
                     min.pos = -1.2, min.vel = -0.07, delta.x, delta.y) {
	policy[round((s[1] - min.pos) / delta.x) + 1, round((s[2] - min.vel) / 
                              delta.y) + 1]
	}
	
  
mountain.car.test.policy <- function(policy = p.mountain, S = NULL, A = c(-1,   
                              0, 1), max.steps = 200, min.vel = -0.07, max.vel =
                              0.07, min.pos = -1.2 , max.pos = 0.5, 
                              v.pos = seq(min.pos, max.pos, l=nrow(policy)), 
                              v.vel = seq(min.vel, max.vel, l=ncol(policy)), 
                              delta.x =v.pos[2] - v.pos[1],
                              delta.y = v.vel[2] -v.vel[1],
                              sd = 0) {

	if (is.null(S)) {
		S <- matrix(0,1,2)
		S[1,1] <- v.pos[round(runif(1,1,length(v.pos)))]
		S[1,2] <- v.vel[round(runif(1,1,length(v.vel)))]
		}

        res <- array(0, nrow(S))  
        for (i in 1:nrow(S)) {
          goal <- FALSE
          num.steps <- 0
          s <- S[i,]
          while(!goal && num.steps < max.steps) {
                         a <- mountain.car.op(s, policy = policy, delta.x =
                                        delta.x, delta.y = delta.y)
                         a <- rnorm(1, A[a], sd)
                         t <- mountain.car.transition(s, a, min.vel, max.vel,
                                             min.pos, max.pos)
                         s <- t$s
                         goal <- t$g
                         ind1 <- round((t$s[1] - min.pos) / delta.x) + 1
                         ind2 <- round((t$s[2] - min.vel) / delta.y) + 1
                         S[1] <- v.pos[ind1]
                         S[2] <- v.vel[ind2]
          #		x <- as.matrix(v.mountain)
          #		x[ind1,ind2] <- mv + 1
          #		image(pos, vel, x, main = paste(A[a]))	
                         
                         num.steps <- num.steps + 1
                         }
               res[i] <- num.steps
               }
	res
	}


mountain.car.test.policy.rbfn <- function(rbfn, S = NULL,
                                          means, stdevs, 
                                          max.steps = 200, 
                                          num.avg = 20,
                                          min.vel = -0.07, max.vel = 0.07,
                                          min.pos = -1.2 , max.pos = 0.5, 
                                          sd = 0,
                                          rbfn.out = rbfn.norm.output, 
                                          A = c(-1, 0, 1), ...) {
   res <- matrix(0, num.avg, nrow(S))
   for (n in 1:num.avg) {
      for (i in 1:nrow(S)) {
         goal <- FALSE
         num.steps <- 0 
         s <- S[i,]
         sn <- (s - means) / stdevs
         while (!goal && num.steps < max.steps) {
            a <- which.max(rbfn.out(rbfn, matrix(sn,1,length(sn))))
            t <- mountain.car.transition(s, A[a] + rnorm(1,0,sd), ...)
#          total.reward <- total.reward + t$r
            num.steps <- num.steps + 1
            s <- t$s
            sn <- (s - means) / stdevs
            goal <- t$g
            }
         res[n,i] <- num.steps
         }
      }
   res
   }



fix.mountain <- function(vm, value = 0) {
# correct the wrong values of state due to self-transitions	
	ind <- vm <= value
	for (i in 1:length(ind)) {
		if (ind[i]) {
			vm[i] <- 0.5 * vm[i-1] + 0.5*vm[i+1]
			}
		}
		
	vm
	}


make.S.mountain <- function(size, min.pos = -1.2, max.pos = 0.5, 
                            min.vel = -0.07, max.vel = 0.07) {
   cp <- seq(min.pos, max.pos, l=size)
   cv <- seq(min.vel, max.vel, l = size)
   gp(cp,cv)
   }



print("mountain.car.R loaded")