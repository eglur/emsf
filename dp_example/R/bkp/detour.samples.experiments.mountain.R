source("mountain.car.R")
source("rbfn.rl.R")
source("detour.samples.R")
source("lspi.rbfn.R")
source("util.R")
source("detour.samples.experiments.R")


mountain.data <- function (size, normalize = TRUE, p = 0, random = TRUE,
                  lsars = FALSE, ...){
	# make sure there is at least one goal in the dataset
	num.goals <- max(round(p*size), 1)
	# The state s = (0.5,0.01) will reach the goal even if a = -1
	s <- gp(0.5, seq(0.07, 0.01, length = num.goals))  
	remaining <- size - num.goals
		
     if (random) {
       # generate the dataset with size-num.goals points
       cp <- runif(remaining, -1.2, 0.5)
       cv <- runif(remaining,-0.07, 0.07)
       s <- rbind(s, cbind(cp,cv))
       }
     else {
       # uniform grid
       grid.side <- floor(sqrt(remaining))
	  cp <- seq(-1.2, 0.5, length = grid.side)
	  cv <- seq(-0.07, 0.07,length = grid.side)
	  s <- rbind(s, gp(cp,cv))
	  # complete the dataset with random points
	  remaining <- size - nrow(s)
       if (remaining > 0) {
         cp <- runif(remaining, -1.2, 0.5)
         cv <- runif(remaining,-0.07, 0.07)
         s <- rbind(s, cbind(cp,cv))
         }
       }
	
	sa <- make.sa(s, 3)
	sars <- collect.transitions(sa, mountain.car.transition, c(-1,0,1),
               normalize = normalize, ...)
	if (lsars) make.lsars(sars)
	else sars
	}
	

mountain.test.policy <- function(S, rbfn, sample.mean, sample.sd, max.it = 
                           300, rbfn.out = rbfn.norm.output) {
#    cp <- seq(-1.2, 0.5, length = 100)
#    cv <- seq(-0.07, 0.07, length = 100)
#    S2 <- gp(cp, cv)
#    S2 <- normalize(S2, means = sample.mean, stdevs = sample.sd)
#    o <- rbfn.out(rbfn, S2)
#    p <- apply(o, 1, max)
#    image(-t(matrix(p,100,100)))
   ns <- array(0, nrow(S))
   for (s in 1:nrow(S)) {
      ns[s] <- control.rbfn(rbfn, S[s,], mountain.car.transition, c(-1,0,1), 
               sample.mean, sample.sd, max.it, rbfn.out)$ns
      }
   ns
   }
   
   
               
             
mountain.lspi <- function(num.avg, sars, S.test, rbfn.size = 3^2, df = 1, tau =
                   1e-1, it = 10, max.it.control = 300, verbose = FALSE,
                   rbfn.out = rbfn.norm.output, cluster = FALSE, ...) {
	res <- matrix(0, num.avg, nrow(S.test))
	for (i in 1:num.avg) {
          if (verbose) print(paste("Run #", i))
          
          C <- NULL
          if (cluster) {
            C <- kmeans(sars$s, rbfn.size)$centers
            }
          else {
            C <- matrix(0, rbfn.size, ncol(S.test))
            max.X <- apply(sars$s, 2, max)
            min.X <- apply(sars$s, 2, min)
            C <- uniform.grid(C, max.X, min.X)
            }
          
          rbfn <- make.rbfn.kbrl(C, 3, tau)       
          rbfn <- lspi.rbfn(sars, rbfn, df, it, verbose = FALSE, ...)
          
	 if (verbose) print("Testing the policy...")
	 res[i,] <- mountain.test.policy(S.test, rbfn, sars$means,
                        sars$stdevs, max.it = max.it.control, rbfn.out =
                        rbfn.out)
          if (verbose) print(paste("Average number of steps:", mean(res[i,]))) 
	  }
      res
      }
	
	
mountain.detour <- function(num.avg, sars, S.test, rbfn.size = 3^2, df = 1, tau
                   = 1e-1, it = 10, max.it.control = 300, verbose = FALSE,...){
   lsars <- make.lsars(sars)
   res <- matrix(0, num.avg, nrow(S.test))
   prop.centers <- rbfn.size / nrow(lsars[[1]]$s) 
   for (i in 1:num.avg) {
      if (verbose) print(paste("Run #", i))
	 rbfn <- detour.rbfn(lsars, tau, df, prop.centers, iter.pi = it,
               run.value.iteration = FALSE, ...)
      if (verbose) print("Testing the policy...")
      res[i,] <- mountain.test.policy(S.test, rbfn, lsars$means,
                        lsars$stdevs, max.it = max.it.control, rbfn.out =
                        rbfn.norm.output)
      if (verbose) print(paste("Average number of steps:", mean(res[i,]))) 
	}
   res
   }
   
   


## consider to remove, since I have a generic version
mountain.lspi.detour <- function(num.avg, sars, S.test, 
                           rbfn.sizes.lspi = c(3^2, 4^2, 5^2), 
                           taus.lspi = c(0.7, 0.5, 0.3, 0.1), 
                           taus.detour =c(1e-1, 1e-3, 1e-5, 1e-7), df = 0.995,
                           it = 10, max.it.control = 300, verbose = FALSE,
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
         
         res <- mountain.lspi(num.avg, sars, S.test,rbfn.sizes.lspi[s], df,
                  taus.lspi[t], it, max.it.control,verbose = FALSE, rbfn.out,
                  cluster)
         res.lspi[s,t,] <- apply(res, 1, mean)
         
         filename <- paste(dir, "mountain_lspi_s", rbfn.sizes.lspi[s], "_t",
                     taus.lspi[t],".txt", sep = "")
         wt(res.lspi[s,t, ], filename)
         if (verbose) print(paste("Average result:", mean(res.lspi[s,t,])))
         }
      
      
      if (verbose) print("Running detour")      
      for (t in 1:length(taus.detour)) {
         if (verbose) {
            print(paste("Running tau", t, "of", length(taus.detour)))
            }
         res <- mountain.detour(num.avg, sars, S.test, rbfn.size.detour, df,
                  taus.detour[t], it, max.it.control, verbose = FALSE,
                  center.function = center.function)
         res.detour[s,t,] <- apply(res, 1, mean)
         filename <- paste(dir, "mountain_detour_s", rbfn.sizes.lspi[s], "_t",
                     taus.lspi[t],".txt", sep = "")
         wt(res.detour[s,t, ], filename)
         if (verbose) print(paste("Average result:", mean(res.detour[s,t,])))
         }
         
      }
   
   list(res.lspi = res.lspi, res.detour = res.detour)
   }     
                                 


# CONCLUSÕES ATÉ AGORA:
# ver ./res/persistent/mountain_lspi_detour.txt


mountain.detour.compare.center <- function(lsars, tau, df, prop.centers,
                              center.functions, Q, num.avg = 20, 
                              verbose = FALSE, sqrt.test.size = 7, ...) {
   cp <- seq(-1.2, 0.5, length = sqrt.test.size)
   cv <- seq(-0.07, 0.07, length = sqrt.test.size)
   S.test <- gp(cp,cv)
   
   V <- apply(Q, 1, max)
   V <- V - mean(V)
   op <- apply(Q, 1, which.max)
   ## think
   op[op == 2] <- 3
   
   grid.side <- sqrt(length(V))
   cp <- seq(-1.2, 0.5, length = grid.side)
   cv <- seq(-0.07, 0.07, length = grid.side)
   S <- gp(cp,cv)
   S <- normalize(S, lsars$means, lsars$stdevs)
   
   V.tmp <- t(matrix(op, grid.side, grid.side))
   image(-V.tmp)
   
   time <- matrix(0, num.avg, length(center.functions))
   sse  <- matrix(0, num.avg, length(center.functions))
   pd  <- matrix(0, num.avg, length(center.functions))
   pt  <- matrix(0, num.avg, length(center.functions))
   
   # first, solve using all the points
   time.all <- system.time(rbfn <- detour.rbfn(lsars, tau, df,1,...),TRUE)[1]
   R <- compute.pv.rbfn(rbfn, S, rbfn.norm.output)
   
   VD <- R$v
   VD <- VD - mean(VD)
   sse.all <- sse(V, VD)
   
   p <- R$p
   pd.all <- sum(p != op) / length(p)
   
   pt.all <- mountain.test.policy(S.test, rbfn, lsars$mean, lsars$stdev)
   
   V.tmp <- t(matrix(p, grid.side, grid.side))
   image(-V.tmp)
     
     
   if (verbose) print(paste("SSE all:", sse.all, "pd all:", pd.all, "pt all:", 
                 pt.all, "time all:", time.all))
   
   for (i in 1:num.avg) {
      for (cf in 1:length(center.functions)) {
         time[i,cf] <- system.time(rbfn <- detour.rbfn(lsars, tau, df,
                        prop.centers, center.functions[[cf]], ...),TRUE)[1] 
         R <- compute.pv.rbfn(rbfn, S, rbfn.norm.output)
         VD <- R$v
         VD <- VD - mean(VD)
         sse[i,cf] <- sse(V, VD) 
         
         p <- R$p
         pd[i,cf] <- sum(p != op) / length(p)
         pt[i,cf] <- mountain.test.policy(S.test, rbfn, lsars$mean, lsars$stdev)

         V.tmp <- t(matrix(p, grid.side, grid.side))
         image(-V.tmp)
         
         if (verbose) print(paste(i,cf, ":: SSE:", sse[i,cf],"pd:", pd[i,cf], 
                        "pt", pt[i,cf], "time:", time[i,cf]))         
         }
      }
   sse <- cbind(sse.all, sse)
   time <- cbind(time.all, time)
   pd <- cbind(pd.all, pd)
   pt <- cbind(pt.all, pt)
    
    list(sse = sse, time = time, pd = pd, pt = pt) 
    }

   
#    CONCLUSÕES ATÉ AGORA:
#    Nos experimentos que eu rodei eu só obtive resultados razoáveis quando
# nrow(lsars[[1]]$s) >= 400. Com menos pontos do que isso, a solução dos MDPs
# reduzidos são triviais (ex., todas as ações iguais a '3'). Não fiz muitos
# testes com o parâmetro 'prop.center', mas pelo que percebi '0.3' é o menor
# valor que gera soluções competitivas. Em relação ao parâmetro 'tau', o melhor
# valor varia. Se os estados são agrupados pela matriz K, então o melhor valor
# para 'tau' fica em torno de 1e-2; no caso em que K é a matriz identidade (ou
# seja, no KBRL padrão), o melhor valor parece ser em torno de 3e-1. 
#    Esse foi o melhor resultado do KBRL usando 400 pontos e tau = 0.3:
# SSE all: 135.76 pd all: 0.3996 pt all: 44.69 time all: 67.153
# Consegui resultados semelhantes com os MDPs reduzidos, usando tau = 1e-2,
# prop.center = 0.3 e os mesmos 400 pontos. Nos experimentos que eu fiz eu
# comparei os resultados dos MDPs reduzidos com o do KBRL usando todos os
# pontos. Seria interessante comparar os MDPs reduzidos com o KBRL usando o
# mesmo número de pontos (ou seja, usando a função 'center.random' para
# configurar os centros, mas usando K = I, onde I é a matriz identidade).
# Tentar rodar os experimentos com df = 1, o que me parece mais adequado nesse
# caso

   
   
   
   
   