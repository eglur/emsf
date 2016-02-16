
lspi.size.to.detour.size <- function(num.rbfs, num.samples, num.actions) {
# Given the number of RBFS, sample transitions and actions, returns the number
# of states/RBFS that should be used by detour so that both this algorithm and
# LSPI would have roughly the same computational cost (this is conservative, in
# the sense that detour will execute less operations than LSPI) 
   c <- -(num.rbfs^2 * num.samples + (num.rbfs * num.actions)^3 
            + num.samples * num.rbfs * num.actions)
   floor(Re(polyroot(c(c, 0, num.actions, 1))[1]))
   }
   


lspi.detour <- function(num.avg, sars, S.test, lspi.function,
                           detour.function,                        
                           rbfn.sizes.lspi = c(3^2, 4^2, 5^2), 
                           taus.lspi = c(0.7, 0.5, 0.3, 0.1), 
                           taus.detour =c(1e-1, 1e-3, 1e-5, 1e-7), df = 0.995,
                           it = 10, max.it.control = 3000, verbose = FALSE,
                           rbfn.out = rbfn.norm.output, dir = "./res/", 
                           cluster = FALSE, ...){
   num.samples <- nrow(sars$s) / 3
   
   res.lspi <- array(0, c(length(rbfn.sizes.lspi), length(taus.lspi), num.avg))
   res.detour <-array(0,c(length(rbfn.sizes.lspi), length(taus.detour),num.avg))
   
   center <- NULL
   if (cluster) center.function <- center.kmeans
   else center.function <- center.grid
   
   
   for (s in 1:length(rbfn.sizes.lspi)) {
      if (verbose) {
         print(paste("Running size", s, "of", length(rbfn.sizes.lspi)))
         }
         
      rbfn.size.detour <- lspi.size.to.detour.size(rbfn.sizes.lspi[s],
                              num.samples, 3)
      print(rbfn.size.detour)                              
                              
      if (verbose) print("Running LSPI")      
      for (t in 1:length(taus.lspi)) {
         if (verbose) {
            print(paste("Running tau", t, "of", length(taus.lspi)))
            }
         
         res <- lspi.function(num.avg, sars, S.test,rbfn.sizes.lspi[s], df,
                  taus.lspi[t], it, max.it.control,verbose = verbose, 
                  rbfn.out = rbfn.out, cluster = cluster)
         res.lspi[s,t,] <- apply(res, 1, mean)
         
         filename <- paste(dir, "lspi_s", rbfn.sizes.lspi[s], "_t",
                     taus.lspi[t],".txt", sep = "")
         wt(res.lspi[s,t, ], filename)
         if (verbose) print(paste("Average result:", mean(res.lspi[s,t,])))
         }
      
      
      if (verbose) print("Running detour")      
      for (t in 1:length(taus.detour)) {
         if (verbose) {
            print(paste("Running tau", t, "of", length(taus.detour)))
            }
         res <- detour.function(num.avg, sars, S.test, rbfn.size.detour, df,
                  taus.detour[t], it, max.it.control, verbose = verbose,
                  center.function = center.function)
         res.detour[s,t,] <- apply(res, 1, mean)
         filename <- paste(dir, "mountain_detour_s", rbfn.sizes.lspi[s], "_t",
                     taus.lspi[t],".txt", sep = "")
         wt(res.detour[s,t, ], filename)
         if (verbose) print(paste("Average result:", mean(res.detour[s,t,])))
         }
         
      }
   
   list(lspi = res.lspi, detour = res.detour)
   }     
   

compile.results <- function(res, rownames, colnames, stat.function = mean) {
   df <- matrix(0, dim(res)[1], dim(res)[2])
   colnames(df) <- colnames
   rownames(df) <- rownames
   for (i in 1:nrow(df)) {
      for (j in 1:ncol(df)) {
         df[i,j] <- stat.function(res[i,j,])
         }
       }
   df
   }

   
   
print("detour.samples.experiments.R loaded")
   
