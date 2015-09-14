dyn.load(paste("./c/matrix_outer_product_bellman", .Platform$dynlib.ext, sep = "")) #melhorar

matrix.outer.product.bellman <- function(A, B, num.rbfs, a1, a2) {
# computes the outer product of two matrices summing the outer products of the individual vectors
# it's faster than the original version, but only works for matrices of the LSPI shape
# only performs the operations with no zeros involved, but for that it needs num.rbfs, a1 and a2
	n.col <- ncol(A)
	O <- matrix(0, n.col, n.col)
	O <- .C("matrix_outer_product_bellman", as.double(t(A)), as.double(t(B)),  
				as.integer(n.col),  as.integer(nrow(A)), as.integer(num.rbfs), 
				as.integer(a1), as.integer(a2), as.double(O))
	
	matrix(O[[length(O)]], n.col, n.col, byrow = TRUE)
	}

print("matrix.outer.product.bellman.R loaded")
