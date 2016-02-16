dyn.load(paste("./c/dist_mat", .Platform$dynlib.ext, sep = "")) #melhorar

dist.mat <- function(M1, M2) {
# computes the square Euclidean distance between the rows of two matrices
   if (is.null(nrow(M1))) M1 <- matrix(M1, 1, length(M1))
   if (is.null(nrow(M2))) M2 <- matrix(M2, 1, length(M2))
   R <- matrix(0,nrow(M1),nrow(M2)) 
	R <- .C("dist_mat", as.double(t(M1)), as.double(t(M2)),
            as.integer(nrow(M1)), as.integer(nrow(M2)), as.integer(ncol(M1)),R)
	matrix(R[[length(R)]], nrow(M1), nrow(M2), byrow = TRUE)
	}


print("dist.mat.R loaded")	
	
