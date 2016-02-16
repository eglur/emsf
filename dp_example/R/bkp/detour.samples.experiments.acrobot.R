source("acrobot.R")
source("rbfn.rl.R")
source("detour.samples.R")
source("util.R")
source("lspi.rbfn.R")

acrobot.data <- function(n = 10, normalize = TRUE, max.values = c(3.12349, 
                  3.14018, 6.15194, 12.97330), min.values = c(-3.06010,
                  -3.14110, -4.88625, -10.30510), lsars = FALSE, ...) {
# 'n' is the number of points in the dataset
# The values for 'max.values' and 'min.values' were obtained in the experiments
# with the simulator written in C++ (see RGD paper: Barreto & Anderson, 2007)
   S <- matrix(0, n, 4)
   #S <- uniform.grid(S, max.values, min.values)
   for (i in 1:4) {
      S[,i] <- runif(n, min.values[i], max.values[i])
      }         
   sars <- collect.transitions(make.sa(S, 3), acrobot.transition, c(-1,0,1),
            TRUE,...)
   if (lsars) make.lsars(sars)
   else sars
   }
   


acrobot.lspi <- function(num.avg, sars, rbfn.size = 3^4, df = 1, tau = 0.05, it
                  = 10, max.it = 1000, verbose = FALSE,...) {
# ainda não testei depois que retirei o teste de dentro do LSPI	
	res <- array(0, num.avg)
	for (i in 1:num.avg) {
		print(paste("Run #", i))
		rbfn <- make.rbfn.tau(sars$s, 3, rbfn.size, FALSE, tau = tau)
		rbfn <-  lspi.rbfn(sars, rbfn, df, it, verbose = verbose, ...)
          res[i] <- acrobot.control(graphic = TRUE, policy = rbfn.policy, rbfn
                        = rbfn, sample.mean = sars$means, sample.stdev =
                        sars$stdevs, max.it = max.it)

		#wt(res, filename)
		}
	res
	}
	
	

acrobot.detour <- function(num.avg, lsars, rbfn.size = 3^4, df = 1, tau = 0.05,
                     it = 10, verbose = FALSE, graphic = FALSE, max.it =       
                     1000, filename, ...) {
	
   res <- array(0, num.avg)
   prop.centers <- rbfn.size / nrow(lsars[[1]]$s) 
   for (i in 1:num.avg) {
      rbfn <- detour.rbfn(lsars, tau, df, prop.centers, ...)
      res[it] <- acrobot.control(graphic = graphic, policy = rbfn.policy, rbfn
                        = rbfn, sample.mean = lsars$means, sample.stdev =
                        lsars$stdevs, max.it = max.it)
      print(res[it])
      }
   res
   }
   
   
print("detour.samples.experiments.acrobot.R loaded")             