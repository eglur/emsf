read.matrix <- function(filename) {
   as.matrix(read.table(filename))
   }

ps <- function(...){
   paste(..., sep="")
   }


wt <- function(T, filename) {
	write.table(T, filename, row.names=FALSE, col.names=FALSE, quote=FALSE)
	}

rt <- function(filename) {
	as.matrix(read.table(filename))
	}

pp <- function(X) {
	persp(as.matrix(X), theta = 35, phi = 35, ticktype = "detailed", shade = 0.5)
	}

pd <- function(X) {
	if (is.null(dim(X))) print(length(X))
     else print(dim(X))
	}

sy <- function(I) {
	system(I)
	}
     
plot.both <- function(X,Y) {
    matplot(cbind(X,Y), t="l")
    }     

runif.int <- function(n, min=0, max=1) {
# returns a random integer between min and max
	round((min - 0.4999) + runif(n) * (max - min + 0.9998))
	}
	
sse <- function(X, Y) {
	sum((X - Y)^2)
	}
	
mse <- function(X, Y) {
	sse(X, Y) / length(X)		
	}
	
euclidean.norm <- function(v) {
   sqrt(crossprod(v,v))
   }
	
norm.inf <- function(A) {
   if (!is.matrix(A)) A <- matrix(A, length(A), 1)
   max(apply(abs(A), 1, sum))
   }

norm.1 <- function(A) {
   if (!is.matrix(A)) A <- matrix(A, length(A), 1)
   max(apply(abs(A), 2, sum))
   }


norm.frobenius <- function(A) {
   sqrt(sum(A^2))
   }



uniform.grid <- function(grid, max.X, min.X, p=1) {
## FIX: no need to pass "grid"; just nrow(grid)
	n <- round(nrow(grid)^(1 / ncol(grid)))
     
      pos <- matrix(0, n, ncol(grid))

	for (i in 1:ncol(grid)) {
		delta <- p * (max.X[i] - min.X[i]) / (2 * n)
		pos[,i] <- seq(min.X[i] + delta, max.X[i] - delta, length = n)
		}

	od <- array(1, ncol(grid))
     
     nc <- n^ncol(grid)
	for (i in 1:nc) {
		for (j in 1:ncol(grid)) grid[i,j] <- pos[od[j], j]
		od <- inc.odometer(od, n, 1)
		}
	# spread the remaining points randomly
	if (nc < nrow(grid)) {
		for (j in 1:ncol(grid)) 	grid[(nc+1):nrow(grid),j] <-
                                    runif(nrow(grid) - nc, min.X[j], max.X[j])
		}
	grid
	}
	

inc.odometer <- function(o, max=1, min=0, i = length(o)) {
	if (o[i] + 1 > max) {
		o[i] <- min
		if (i != 1) {
			i <- i - 1
			o <- inc.odometer(o, max, min, i)
			}
		}
	else o[i] <- o[i] + 1
	o
	}

colors.plot <- function(n) {
     cc <- c("DARKGREY", "RED", "BLUE", "GREEN3", "VIOLET","ORANGE",
             "DARKORCHID", "CORAL4", "DARKCYAN","DARKBLUE")
     ind <- 1:n %% length(cc)
     ind[ind==0] <- length(cc)
     cc[ind]
     }

reduce.cols <- function(D, ncols) {
   X <- matrix(0, nrow(D), ncols)
   ind <- seq(1, ncol(D), length=ncols+1)
   for (i in 1:(length(ind)-1))
   {
      X[,i] <- apply(D[,ind[i]:(ind[i+1]-1)], 1, mean)
   }
   X
}


sum.over.cols <- function(D)
{
   for (i in 2:ncol(D))
   {
      D[,i] <- apply(as.matrix(D[,(i-1):i]),1,sum)
   }
   D
}


print("util.R loaded")