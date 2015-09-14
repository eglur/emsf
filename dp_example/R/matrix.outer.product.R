dyn.load(paste("./c/matrix_outer_product", .Platform$dynlib.ext, sep = "")) #melhorar

matrix.outer.product <- function(A, B) {
# computes the outer product of two matrices summing the outer products of the individual vectors
	n.col <- ncol(A)
	O <- matrix(0, n.col, n.col)
	O <- .C("matrix_outer_product", as.double(t(A)), as.double(t(B)),  
				as.integer(n.col),  as.integer(nrow(A)), as.double(O))
	
	matrix(O[[length(O)]], n.col, n.col, byrow = TRUE)
	}

print("matrix.outer.product.R loaded")