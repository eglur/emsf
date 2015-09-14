dyn.load(paste("./c/proj_sum_one", .Platform$dynlib.ext, sep = "")) #melhorar

# 'proj.sum.one' projects the vectors in the rows of X onto the ncol(X)-dimensional space
# composed by the vectors whose elements sum up to 1.

proj.sum.one <- function(X) {
	if (is.null(dim(X))) X <- t(as.matrix(X))
	U <- matrix(0, nrow(X), ncol(X))
	S <- apply(X, 1, sum)
	
	# sort the elements and keep the indices
	O <- matrix(apply(X, 1, order, decreasing = TRUE), nrow(X), ncol(X), byrow = TRUE)
	for (i in 1:nrow(X)) X[i,] <- X[i, O[i,]]
	
	U <- .C("proj_sum_one", as.double(t(X)), as.integer(nrow(X)), 
				  as.integer(ncol(X)), as.double(S), as.double(t(U)))
	U <- matrix(U[[length(U)]], nrow(X), ncol(X), byrow = TRUE)
	
	# put the elements in the previous order
	for (i in 1:nrow(U)) U[i, O[i,]] <- U[i,]
	U
	}

print("proj.sum.one.R loaded")	
	
