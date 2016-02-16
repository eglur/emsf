

puddle.world.transition <- function(s, a, 
   # Parameters set as in Sutton, 96, "Generalization..."
   regular.reward = 0, goal.reward = 5, puddle.reward = -10,
   step = 0.05, sd = 0.01, goal.radius = 0.05, 
   # vertical and horizontal puddles
   ph = list(x = c(0.10, 0.45), y = 0.75, r = 0.1), # y1 == y2, x1<x2 
   pv = list(x = 0.45, y = c(0.4, 0.8), r = 0.1)   # x1 == x2, y1<y2 
   ) {
   
   if (a == "N") s[2] <- min(1, s[2] + step)
   else if (a == "S") s[2] <- max(0, s[2] - step)
        else if (a == "E") s[1] <- min(1, s[1] + step)
             else if (a == "W") s[1] <- max(0, s[1] - step)
                  else print(paste("Puddle world warning: invalid action",a))
   s <- s + rnorm(2,0,sd)
   
   r <- regular.reward
   g <- FALSE
   if (s[1] > 1 - goal.radius && s[2] > 1 - goal.radius) {
      r <- goal.reward
      g <- TRUE
      }
   else { # Must check the puddles...
      ins <- 0
      # Horizontal puddle
         if (s[1] >= ph$x[1] && s[1] <= ph$x[2] && 
             s[2] >= ph$y - ph$r && s[2] <= ph$y + ph$r) { 
            # s is inside ph's box
            if (s[2] > ph$y) ins <- ph$y + ph$r - s[2]
            else ins <- s[2] - ph$y + ph$r
            }
         else {
            if (s[1] < ph$x[1]) {
            # check left circunference
               d <- (s[1]-ph$x[1])^2  + (s[2]-ph$y)^2 - ph$r^2
               if (d < 0) ins <- -d
               }
            else if (s[1] > ph$x[2]) {
               # check right circunference
               d <- (s[1]-ph$x[2])^2  + (s[2]-ph$y)^2 - ph$r^2
               if (d < 0) ins <- -d
               }
            }
            
      # Vertical puddle
         if (s[2] >= pv$y[1] && s[2] <= pv$y[2] && 
            s[1] >= pv$x - pv$r && s[1] <= pv$x + pv$r) { 
            # s is inside pv's box
            if (s[1] > pv$x) ins <- max(ins, pv$x + pv$r - s[1])
            else ins <- max(ins, s[1] - pv$x + pv$r)
            }
         else {
            if (s[2] < pv$y[1]) {
            # check bottom circunference
               d <- (s[1]-pv$x)^2  + (s[2]-pv$y[1])^2 - pv$r^2
               if (d < 0) ins <- max(ins, -d)
               }
            else if (s[2] > pv$y[2]) {
               # check top circunference
               d <- (s[1]-pv$x)^2  + (s[2]-pv$y[2])^2 - pv$r^2
               if (d < 0) ins <- max(ins, -d)
               }
            }
       r <- r + puddle.reward * ins
       }

   list(s = s, r = r, g = g)
	}



puddle.test.policy <- function(Q, S, df, max.steps=300, ...) {
# Q is n xn x |A|, where n = sqrt(|S|)
   A <- c("N", "S", "E", "W")
   delta <- 1/ (nrow(Q) - 1)
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
         indx <- min(max(s[1] %/% delta + 1, 1), nrow(Q))
         indy <- min(max(s[2] %/% delta + 1, 1), ncol(Q))
         a <- which.max(Q[indx, indy,])
         t <- puddle.world.transition(s, A[a], ...)
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


print("puddle.world.R loaded")