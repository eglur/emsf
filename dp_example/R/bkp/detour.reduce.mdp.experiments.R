source("dp.R") 
source("util.R") # for wt()


kl <- function(P,Q) {
# Kullback_leibler divergence
   sum(P * log(P/Q) - P + Q)
   }
		

reduce.mc.experiments <- function(factor.functions, sizes = c(100), df = c(0.9, 
                           0.99,0.999), perc=c(0.2, 0.15, 0.1, 0.05), num.mcs =
                           5, num.avg = 20, perc.iter = c(0.3, 0.5,0.7), names
                           = paste(1:length(factor.functions)), load.mc = TRUE,
                           dp.function = solve.mc, min.reward = 0, max.reward =
                           1, dir = "./res/", ... ) {
# 'sizes' is the number of states |S| of the markov chain (MC)
# 'perc' is the percentual of reduction of the MC, i.e., k = perc * |S|
# 'df' are the discount factors
# 'num.mcs' is the number of MCs generated for each size
# 'num.avg' is the number of times each MC will be reduced
# 'dp.function' is the function used to find the value function 
#    (currently, either "value.iteration" or "solve.mc")
   perc.iter <- sort(perc.iter)
   for (s in 1:length(sizes)) {
      for (i in 1:num.mcs) {
         # create (or load) MC
         P <- NULL
         R <- NULL
         W <- NULL
         filename <- paste(dir, "random_mc_RP_s", sizes[s], "_", i, ".txt",
                        sep = "")
         if (!load.mc) {
            P <- normal.transition.matrix(sizes[s], sd = 0.5)
            R <- runif(sizes[s], min.reward, max.reward)
            W <- cbind(R, P)
            wt(W, filename)
            }
         else {
            W <- as.matrix(read.table(filename))
            R <- W[,1]
            P <- W[, 2:ncol(W)]
            }
         
         Q <- matrix(0, length(R), length(df))
         print(paste("Solving the MC", i,"exactly for..."))
         for (d in 1:length(df)) {
            # "solve" the MC exactly for each discount factor
            print(paste("Discount factor =", df[d]))
            Q[,d] <- dp.function(R, P, df = df[d])
            Q[,d] <- Q[,d] - mean(Q[,d]) # center the data
            }
         
         for (f in 1:length(factor.functions)) {
            print(paste("Running", names[f]))
            for (p in 1:length(perc)) {
               k <- round(perc[p] * sizes[s])
               print(paste("Perc.:", perc[p], " (k =", k,")"))
               
               iter <- floor(perc.iter * (nrow(P)^3 - k^3) / (nrow(P)^2 * k))
				      
               print(c("Running", iter,"iterations"))
               res <- array(0, c(num.avg, 3 + length(df), length(perc.iter)))
                     
               for (j in 1:num.avg) {
                  print(paste("Run #", j))
                  D <- NULL
                  iter.total <- 0
                  for (t in 1:length(iter)){
                     iter.max <-  iter[t] - iter.total
                     iter.total <- iter.total + iter.max
                     D <- factor.functions[[f]](W, k, iter.max = iter.max, K =
                           D$K, L = D$L, R = D$R, ...)
                     P2 <- D$K %*% D$L
                     res[j,1,t] <- mse(R, D$L %*% D$R)  
                     LK <- D$L %*% D$K
                     res[j,2,t] <- mse(P,LK )
                     res[j,3,t] <- kl(P, LK)
                        
                     for (d in 1:length(df)) {
                        Q2 <- D$L %*% dp.function(D$R, P2, df = df[d])
                        Q2 <- Q2 - mean(Q2)
                        res[j, 3 + d, t] <-   mse(Q[,d], Q2)
                        }
                     }
                  }
                     
               # 'tmp' is used to keep the colnames
               tmp <- matrix(0, dim(res)[1], dim(res)[2])
               cn <-  c("MSE R", "MSE P", "KL P") # column names
               colnames(tmp) <- c(cn, paste("MSE V", df ))
               
               for (t in 1:length(iter)) {
                  filename <- paste(dir, "random_mc_s", sizes[s], "_", names[f],
                              "_p", perc[p], "_it", perc.iter[t], "_", i,
                              ".txt", sep="")
                  tmp[,] <- res[,,t] # just to keep the colnames
                  write.table(tmp, filename, quote = FALSE, row.names=FALSE,
                     col.names=TRUE)
                  }
               }
            }
         }
      }
   print("All done!!")
   } 
   



concatenate.results <- function(dir = "./res/", names = c("nmf", "kmeans",
                        "fuzzy_kmeans_m1.25"), columns = 1, sizes = 100, df =
                        c(0.9, 0.99,0.999), perc = c(0.5, 0.3, 0.1, 0.05),
                        num.mcs = 5) {
   res <- NULL
   for (n in names) { 
      for (s in sizes) {
         for (p in perc) {
            t <- NULL
            for (i in 1:num.mcs) {
               filename <- paste(dir, "random_mc_s", s, "_", n, "_p", p, "_", i,
                              ".txt", sep="")
               t <- rbind(t, read.table(filename, header = TRUE)[,columns])
               }
            colnames(t) <- paste(n, colnames(t))
            res <- cbind(res, as.matrix(t))
            }
         }
      }
   as.data.frame(res)
   }
				

	
reduce.mdp.experiments <- function(size = 100, num.actions = 3, 
                           perc = 0.5, factor.function = factor.kmeans, num.avg
                           = 20, sd = 1, min.reward = 0, max.reward = 1, df =
                           0.99, reduce.mdp.function = reduce.mdp,...){
   res <- matrix(0, num.avg, 2)
   for (i in 1:num.avg) {
      P <-mdp.transition.matrix(size,num.actions,normal.transition.matrix,
            sd=sd)
      R <- matrix(runif(size * num.actions, min.reward, max.reward), size,
            num.actions)
      ## remove
      #P[1,,] <- 0
      #P[1,1,] <- 1
      #R[1,] <- 1
      			
      t <- system.time(
         Q <- value.iteration(R, P, df = df, span = TRUE)$Q,
         TRUE)[1]
		
      t2 <- system.time(
         D <- reduce.mdp.function(R, P, perc,factor.function =
               factor.function,...),
         TRUE)[1]
      t2 <- t2 + system.time(
         Q2 <- value.iteration(D$R, D$P, df = df, span = TRUE)$Q, 
         TRUE)[1]
      t2 <- t2 + system.time(
         Q2 <- expand.Q(D,Q2), 
         TRUE)[1]
			
      p <- apply(Q,1,which.max)
      p2 <-  apply(Q2,1,which.max)
      res[i,1] <- sum(p != p2) / length(p)
      res[i,2] <- t / t2
      }
   res
   }
		
		
	
print("detour.reduce.mdp.experiments.R loaded")						
				
												  
