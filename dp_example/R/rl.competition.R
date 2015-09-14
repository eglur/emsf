source("util.R")
source("data.plot.R")

show.heli <- function(base = "steps_fixed", avg = 50:79, number.points = 100, 
		      dir = "~/helicopter/results/steps/", ...)
{
   
   set.par()
  
  batch <- 600000 / (number.points + 1)
  D <- matrix(0, number.points, length(avg))
  
  for (i in 1:length(avg))
  {
    fn <- paste(dir, base, "_", avg[i], ".txt", sep="")
    T <- as.matrix(read.table(fn))
    
    k <- 1
    d <- 0
    sum <- 0
    for (j in 1:nrow(T))
    {
       sum <- sum + T[j,1]
       d <- d + 1
       if (sum >= batch)
       {
         D[k,i] <- sum / d
         d <- 0
         sum <- 0
         k <- k + 1
       }
    }
  
  }
  
  matplot(D, t="l", lwd = 2,  ylab = "Steps flying", xlab = paste("Transitions (x",  600000 / number.points, ")"), ...)
 
}


fix.steps <- function(base = "steps", avg = 50:79,
		      dir = "~/fontes/rlglue/rl_competition_2013/helicopter/results/steps/")
{
  
  for (i in 1:length(avg))
  {
    fn <- paste(dir, base, "_", avg[i], ".txt", sep="")
    T <- as.matrix(t(read.table(fn)))
    wt(T, paste(dir, "steps_fixed", "_", avg[i], ".txt", sep=""))
  }
  D
  
}


print("rl.competition.R loaded")