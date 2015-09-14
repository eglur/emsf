source("tese.experiments.mountain.R")

mountain.kbrl.kbsf<- function(
                              num.archs = 100,
                              num.transitions = c(200,300,400,500,1000),
                              tau.k = 1e-3, tau.q = 1e-3, tau.a = 0.01,
                              random = FALSE,
                              sd = 0,
                              gamma = 0.95, epsilon = 0.01, iter.max = Inf,
                              names = c("mountain_or", "kbrl_mountain","kbsf_mountain"),
                              vm.file = "./files/mountain.value.function.200_by_200.txt",
                              dir = "./fig/",
                              run = c(TRUE, TRUE),
                              p = 0,
                              col=grey.colors(10),
                              verbose = FALSE) {


vm <- as.matrix(read.table(vm.file))
n <- sqrt(length(vm))
cp <-  seq(-1.2, 0.5, l=n)
cv <-  seq(-0.07, 0.07, l=n)
S <- gp(cp,cv)
                  

image(cp, cv, vm,xlab=expression(x), ylab=expression(dot(x)),
      col=col)
dev.copy2eps(file=paste(dir,names[1],".eps",sep=""))

if (verbose) print("Collecting data...")
sars <- mountain.data(num.archs, normalize = TRUE, p = p, random =
                        random, lsars = FALSE, sd = sd)
lsars <- make.lsars(sars)
                    
if (run[1]) {
	if (verbose) print("Running KBRL")
   # run KBRL
   rbfn <- kbrl(lsars, tau.k, gamma, epsilon, iter.max)
	SN <- normalize(S, means = lsars$means, stdevs = lsars$stdevs)
   qm <- rbfn.gaussian.output(rbfn, SN)
   vm <- apply(qm, 1, max)
   vm <- matrix(vm, n, n, byrow=TRUE)
	image(cp, cv, vm,xlab=expression(x), ylab=expression(dot(x)),
         col=col)
   dev.copy2eps(file=paste(dir,names[2],".eps",sep=""))
  	}
               

if (run[2]) {
	if (verbose) print("Running KBSF")
   # run KBSF
   for (a in 1:length(num.transitions)) {
   	if (verbose) print(paste("Number of transitions = ", num.transitions[a]))
   	sars2 <- mountain.data(num.transitions[a], 
   					normalize = TRUE, p = p, random = random, 
   					lsars = FALSE, sd = sd)
      sars <- sars2 #merge.sars(sars, sars2, unnormalize=c(TRUE,TRUE), normalize=TRUE)
      lsars <- make.lsars(sars)
                        
      rbfn <- kbrl.arch(lsars, tau.q= tau.q, tau.a = tau.a,
         num.archs = num.archs, df =gamma, epsilon.vi =
         epsilon, iter.vi =iter.max)
         
		  SN <- normalize(S, means = lsars$means, stdevs = lsars$stdevs)
        qm <- rbfn.gaussian.output(rbfn, SN)
        vm <- apply(qm, 1, max)
        vm <- matrix(vm, n, n, byrow=TRUE)
       image(cp, cv, vm,xlab=expression(x), ylab=expression(dot(x)),
              col=col)
       dev.copy2eps(file=paste(dir,names[3],num.transitions[a],".eps",sep=""))
  		}      
   }
}

