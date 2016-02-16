source("util.R")


## Guestrin's SysAdmin 

pg.ring <- function(i, s, a) {  
# network topology is a ring
     prob <- 0
     if (a == i) prob <- 1
     else {
          j <- i - 1
          if (j == 0) j <- length(s)
          if (s[i] == 1 && s[j] == 1) prob <-  0.95
          if (s[i] == 1 && s[j] == 0) prob <-  0.475
          if (s[i] == 0 && s[j] == 1) prob <-  0.0475
          if (s[i] == 0 && s[j] == 0) prob <- 0.0238
          }
     prob
     }


pg.star <- function(i, s, a) { 
# network topology is a star (with a central server)
     prob <- 0
     if (a == i) prob <- 1
     else {
          j <- 1 # (server)
          if (j == 0) j <- length(s)
          if (s[i] == 1 && s[j] == 1) prob <-  0.95
          if (s[i] == 1 && s[j] == 0) prob <-  0.475
          if (s[i] == 0 && s[j] == 1) prob <-  0.0475
          if (s[i] == 0 && s[j] == 0) prob <- 0.0238
          }
     prob
     }


prob.sysadmin <- function(o,o2,a, pg.function = pg.ring) {
     prob <- 1
     for (i in 1:length(o)) {
          prob.1 <- pg.function(i, o, a)
          if (o2[i] == 1) prob <- prob * prob.1
          else prob <- prob * (1 - prob.1)
          }
     prob
     }

P.sysadmin <- function(size,a, pg.function = pg.ring) {
     n <- 2^size
     P <- matrix(0,n,n)
     o <- array(0, size)
     for (i in 1:n) {
          o2 <- array(0, size)
          for (j in 1:n) {
               P[i,j] <- prob.sysadmin(o,o2,a, pg.function)
               o2 <- inc.odometer(o2,1,0)               
               }
          o <- inc.odometer(o,1,0)
          }
     P
     }

mdp.sysadmin <- function(size, pg.function = pg.ring) {
     n <- 2^size
     M <- NULL
     for (a in 0:size) {
        M <- rbind(M, P.sysadmin(size,a, pg.function))  
     }
     M
}

## end(Guestrin)

## Elevator

bin2dec <- function(v) {
  num <- 0
  for (i in length(v):1) if (v[i]) num <- num + 2^(length(v)-i)
  num
}

ind.state <- function(elev, build, num.floors) bin2dec(c(elev[2:length(elev)], build)) * num.floors + elev[1] 

people.disappeared <- function(build1, build2) {
    flipped <- FALSE
    i <- 1
    while (!flipped && i <= length(build1)) {
      if (build1[i] == 1 && build2[i] == 0) flipped <- TRUE
      i <- i + 1
      }
    flipped
    }


prob.trans.build <- function(build1, build2, prob) {
    if (people.disappeared(build1, build2)) 0
    else {
      prob.total <- 1
      for (i in 1:length(build1)) {
	if (build1[i] == 0) {
	  if (build2[i] == 1) prob.total <- prob.total * prob
	  else prob.total <- prob.total * (1 - prob)
	}
      }
      prob.total
      }
    }



prob.trans <- function(elev, build, a, num.floors, prob){
# a = 1 keep current dir
# a = 2 stop and down
# a = 3 stop and up

  # first, change the part of the state vector describing the elevator
  n.elev <- elev # next elevator state
  if (a == 1) { # keep going
    if (n.elev[2] == 0 && n.elev[1] > 1)          n.elev[1] <- n.elev[1] - 1 # one floor down
    if (n.elev[2] == 1 && n.elev[1] < num.floors) n.elev[1] <- n.elev[1] + 1 # one floor up
    }
  else {
    if (n.elev[1] == 1) n.elev[4] <- 0 # people going down leave the elevator
    if (n.elev[1] == num.floors) n.elev[3] <- 0 # people going up leave the elevator

    if (a == 2) n.elev[2] <- 0  # stop and down
    if (a == 3) n.elev[2] <- 1  # stop and up
    }
  
   # Now, check every possible combination of elev and build
   p <- array(0, num.floors * 2^(3 + length(build))) 
   n.build <- array(0, length(build))
   for (i in 1:(2^(length(build)))) {
	 j <- ind.state(n.elev, n.build, num.floors)
         p[j] <- prob.trans.build(build, n.build, prob)
         n.build <- inc.odometer(n.build,1,0)               
         }
  p
  }
    
P.elev <- function(a, num.floors = 3, prob = 0.2) {
   comb <- 2^(3 + 2 * num.floors)
   size <- num.floors * comb
   
   P <- matrix(0, size, size)
   ind <- 1
   s <- array(0, 3 + 2 * num.floors)
   for (i in 1:comb) {
    for (floor in 1:num.floors) {
	elev  <- c(floor, s[1:3])
	build <- s[4:length(s)]
	P[ind,] <- prob.trans(elev, build, a, num.floors, prob)
	ind <- ind + 1
	}
    s <- inc.odometer(s, 1, 0)
    }
  P
  }


mdp.elev <- function(num.floors = 3, prob = 0.2) {
  M <- NULL
  for (a in 1:3) M <- rbind(M, P.elev(a, num.floors, prob))
  M
  }



