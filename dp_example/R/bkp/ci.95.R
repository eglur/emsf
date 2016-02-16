ci.95 <- function(mean, sd, n) {
#computes the 95% confidence interval given the mean and sd
 	se.95 <-1.96 *  sd / sqrt(n)
	cbind(mean - se.95, mean + se.95)
	}
     
ci.95.sample <- function(data, p.max = 0.05) {
## CONFERIR O VALOR DE p.max!!!    
# computes the 95% confidence interval of a population, but check for 
#     normality before
# 'data' is a matrix with columns representing variables and rows representing
#     cases
   normal <- TRUE
   for (i in 1:ncol(data)) {
      p <- shapiro.test(data[,i])$p.value
      if (p > p.max) normal <- FALSE
      }
   if (normal) {
      means <- apply(data, 2, mean)
      sds   <- apply(data, 2, sd)
      n     <- nrow(data)
      ci <- ci.95(means,sds,n)
      cbind(ci[,1], means, ci[,2])
      }
   else NULL
   }
   
print("ci.95.R loaded")   
