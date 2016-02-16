

mc.transition <- function(s, a,  
  min.vel = -0.07, max.vel = 0.07,
  min.pos = -1.2 , max.pos = 0.5, 
  regular.reward = 0, goal.reward = 1,
  sd = 0) {
# s[1] is the car position                                                       
# s[2] is the car velocity
# a is an index for A = c(-1, 0, 1)
  A <- c(-1, 0, 1)   
  r <- regular.reward # standard reward
  g <- FALSE # in principle, the next state is not a goal
  p <- 1

	if (s[1] > max.pos) { # if the current s is a goal no changes will be made to s
# 		s[1] <- max.pos
		r <- 0 
		g <- TRUE 
		p <- 1
		}
	else {
    A[a] <- A[a] + rnorm(1,0,sd)
      
		s[2] <- s[2] +  0.001 * A[a] - 0.0025 * cos(3 * s[1])
		s[1] <- s[1] + s[2]
		
		if (s[2] > max.vel) s[2] <-  max.vel
		else if (s[2] < min.vel) s[2] <-  min.vel
			
		if (s[1] < min.pos) {
			s[1] <-  min.pos;
			s[2] <-  0;
			}
		else if (s[1] > max.pos) { # s here is already the next state
# 			s[1] <- max.pos
			r <- goal.reward;
			g <- TRUE
			p <- 1
			}
		}
	
	list(s = s, r = r, g = g, p =p) # p is the probability of the transition; g indicates if the state is one of the goals
	}



mc.state.index <- function(s, cells,
  min.vel = -0.07, max.vel = 0.07,
  min.pos = -1.2 , max.pos = 0.5 
  )
{
# cells = array(2) gives the number of cells in each dimension
  if (s[1] > max.pos) cells[1] * cells[2] + 1
  else
  {
    delta1 <- (max.pos - min.pos) / cells[1]
    delta2 <- (max.vel - min.vel) / cells[2]
    ind1 <- min((s[1] - min.pos) %/% delta1 + 1, cells[1])
    ind2 <- min((s[2] - min.vel) %/% delta2, cells[2] - 1)

    ind2 * cells[1] + ind1 
  }
}


mc.evaluate.policy <- function(pi, cells,                               
  test.set = rbind(gp(seq(-1,0.15,l=5), seq(-0.07,0.02,l=5)), c(0,0)),
  max.steps = 300)
{
  
  R <- array(0, nrow(test.set))

  for (i in 1:length(R))
  {
    
    s <- test.set[i,]
    g <- FALSE
    st <- 1
    while (!g && st <= max.steps)
    {
   
      a <- pi[mc.state.index(s, cells)]

      t <- mc.transition(s, a)
      
      g <- t$g
      s <- t$s
      st <- st + 1
      R[i] <- R[i] + t$r
    }
    
  }
  mean(R)
} 


print("mountain.car.R loaded")