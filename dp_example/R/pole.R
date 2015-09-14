dyn.load(paste("./c/pole", .Platform$dynlib.ext, sep = "")) #melhorar

pole.transition <- function(s, a, noise.sd = 0) {
   r <- 0
   g <- 0
   a <- rnorm(1, a, noise.sd)
   R <- .C("pole_transition", as.double(t(s)), as.double(a), as.double(r),
            as.integer(g))
   g <- FALSE
   if (R[[4]] == 1) g <- TRUE
   list(s = R[[1]], r = R[[3]], g = g)
   }


pole.test.policy <- function(Q, S, df, 
                           max.steps=3000, 
                           A = c(-10,10),
                           max.pos = c(2.4, 4.8, pi/5, 4.4),
                           min.pos = -c(2.4, 4.8, pi/5, 4.4)) {
# Q is n x n x n x n x |A|, where n os number of partitions
   n <- dim(Q)[1]

   delta <- (max.pos-min.pos)/ (n - 1)
   
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
         ind <- pmin(pmax((s - min.pos) %/% delta + 1, 1), n)
         a <- which.max(Q[ind[1], ind[2],ind[3],ind[4],])
         t <- pole.transition(s, A[a])
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
 
 

print("pole.R loaded")