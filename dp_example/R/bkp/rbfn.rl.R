compute.v.rbfn <- function(rbfn, S, rbfn.out = rbfn.norm.output) {
	# each output unit of the rbfn is an action
	apply(rbfn.out(rbfn, S), 1, max)
	}
	
compute.policy.rbfn <- function(rbfn, S, rbfn.out = rbfn.norm.output) {
	apply(rbfn.out(rbfn, S), 1, which.max)
	}
	
compute.pv.rbfn <- function(rbfn, S, rbfn.out = rbfn.norm.output) {
	o <- rbfn.out(rbfn, S)
	v <- apply(o,1,max)
	p <- apply(o,1,which.max)
	list(p = p, v = v)
	}
	
compute.v.rbfns <- function(rbfns, S, rbfn.out = rbfn.norm.output) {
	# each RBF network is an action
	o <- matrix(0, nrow(S), length(rbfns))
	for (i in 1:length(rbfns)) {
		o[,i] <- rbfn.out(rbfns[[i]], S)
		}
	apply(o, 1, max)
	}


control.rbfn <- function(rbfn, s, transition.function, A, sample.mean,
                  sample.sd, iter.max = 300, rbfn.out = 
                  rbfn.norm.output, df = 1, epsilon = 0, ...) {
	total.reward <- 0
	num.steps <- 0
	goal <- FALSE                                    
	ddf <- 1
# 	plot(0,0,xlim=c(0,1),ylim=c(0,1))
	while (!goal && num.steps < iter.max) {
		s.n <- (s - sample.mean) / sample.sd
      a <- which.max(rbfn.out(rbfn, matrix(s.n,1,length(s.n))))
      if (runif(1) < epsilon) a <- sample(1:length(A),1)
		t <- transition.function(s,A[a], ...)
		total.reward <- total.reward + ddf * t$r
		num.steps <- num.steps + 1
		s <- t$s
		goal <- t$g
		ddf <- ddf * df
# 		points(s[1],s[2])
      }
	list(ns = num.steps, tr = total.reward)
	}


control.rbfns <- function(rbfns, s, transition.function, A, sample.mean,
                  sample.sd, iter.max = 300, rbfn.out = 
                  rbfn.norm.output, df = 1, epsilon = 0, ...) {
 # one RBFN for each action in A
   total.reward <- 0
   num.steps <- 0
   goal <- FALSE  
   ddf <- 1
   q <-  array(0, length(A))
   while (!goal && num.steps < iter.max) {
      s.n <- (s - sample.mean) / sample.sd
      for (i in 1:length(q)) {
          q[i] <- rbfn.out(rbfns[[i]], matrix(s.n,1,length(s.n)))
          }
      a <- which.max(q)
      if (runif(1) < epsilon) a <- sample(1:length(A),1)
      
      t <- transition.function(s,A[a], ...)
      total.reward <- total.reward + ddf * t$r
      num.steps <- num.steps + 1
      s <- t$s
      goal <- t$g
      ddf <- ddf * df
      }
   list(ns = num.steps, tr = total.reward)
   }


compute.policy.rbfn <- function(S, rbfn, rbfn.out = rbfn.norm.output) {
	apply(rbfn.out(rbfn, S), 1, which.max)
	}


print("rbfn.rl.R loaded")