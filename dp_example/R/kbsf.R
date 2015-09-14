kbsf  <- function(lsars, rbfn.q, tau.a, df, 
                  epsilon.vi = 0.1, iter.vi = Inf, span.vi = TRUE, 
                  iter.pi = 15, goal.reward = 1, 
                  run.value.iteration = FALSE, p.centers = 1,
                  rbfn.function = make.rbfn.simple, 
                  verbose = FALSE, 
                  use.c.code = TRUE, ...) {
    
# The transitions in lsars associated with action a should all begin at the same
#   set of states  
     
     if (verbose) print("Creating the MDP. This may take a while...")
     num.actions <- lsars$num.actions
     num.archs <- nrow(rbfn.q$c)
 
     R <- matrix(0, num.archs, num.actions)
     P <- array(0, c(num.archs, num.archs, num.actions))
     G <- matrix(1, num.archs, num.actions)
      
   
     for (a in 1:num.actions) {
         # compute the detour matrix
         if (verbose) print("Computing rbfn.q output")
         D <- rbfn.norm.design.matrix(rbfn.q, lsars[[a]]$s2)
         
         # compute the return matrix
         if (verbose) print("Making rbfn.a")
         rbfn.a <- rbfn.function(lsars[[a]]$s,tau.a, num.actions,
                                  p.centers = p.centers, 
                                  use.c.code = use.c.code, ...)
   
      if (verbose) print("Computing rbfn.a output")
         K <- rbfn.norm.design.matrix(rbfn.a, rbfn.q$c)
                  
         if (verbose) print("Multiplying...")
         P[,,a] <- K %*% D
         for (i in 1:nrow(K)) R[i,a] <- K[i,] %*% lsars[[a]]$r
         G[,a] <- goal.reward - R[,a]
         }
      if (run.value.iteration) {
         if (verbose) print("Running value iteration")   
         rbfn.q$w <- value.iteration(R, P, G = G, epsilon = epsilon.vi, df =
                                   df, iter.max=iter.vi, span=span.vi)$Q
         }
      else {
         if (verbose) print("Running policy iteration")   
         rbfn.q$w <- policy.iteration(R, P, df = df, max.iter = iter.pi)$Q
         }
      rbfn.q
#       Maneira alternativa de se gerar a RBF, incluindo os pontos s2
#      (embora teoricamente pareça melhor, na prática gera resultados piores)
#       rbfn <- make.rbfn.kbrl(lsars[[1]]$s,tau.a,
#                                         num.actions,num.neighbors,...)
#       for (a in 1:num.actions) {
#           vy <- apply(rbfn.norm.output(rbfn_q, lsars[[a]]$s2),1,max)
#           rbfn$w[,a] <-  vy + lsars[[a]]$r 
#           }
#       rbfn
    }

print("kbsf.R loaded")