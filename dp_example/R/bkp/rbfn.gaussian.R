source("rbfn.R")

dyn.load(paste("./c/rbfn_gaussian", .Platform$dynlib.ext, sep = "")) #melhorar

rbfn.gaussian.design.matrix <- function(rbfn, X) {
     H <- matrix(0, nrow(X), nrow(rbfn$c)) 
	H <- .C("rbfn_gaussian_design_matrix", as.double(t(X)), 
            as.integer(nrow(X)), as.integer(ncol(X)), as.double(t(rbfn$c)),
            as.double(rbfn$s), as.integer(nrow(rbfn$c)),as.double(t(H)))
	H <-	matrix(H[[length(H)]], nrow(X), nrow(rbfn$c), byrow = TRUE)
	H
	}

rbfn.gaussian.output <- function(rbfn, X)  {
   rbfn.gaussian.design.matrix(rbfn, X) %*% rbfn$w
   }

rbfn.norm.design.matrix <- function(rbfn, X) {
	# I didn't call the previous function to improve efficiency
     H <- matrix(0, nrow(X), nrow(rbfn$c)) 
	H <- .C("rbfn_gaussian_design_matrix", as.double(t(X)), 
            as.integer(nrow(X)), as.integer(ncol(X)), as.double(t(rbfn$c)),
            as.double(rbfn$s), as.integer(nrow(rbfn$c)), as.double(t(H)))
	H <-	matrix(H[[length(H)]], nrow(X), nrow(rbfn$c), byrow = TRUE)
	H / apply(H, 1, sum)
	}
	


rbfn.norm.output <- function(rbfn, X) {
   H <- rbfn.norm.design.matrix(rbfn, X) 
   H[is.na(H)] <- 0 ## Isso foi feito para lidar com altas dimensões. REPENSAR
   H %*% rbfn$w
   }


rbfn.bias.design.matrix <- function(rbfn, X) {
     H <- matrix(0, nrow(X), nrow(rbfn$c) +1) # make room for bias
	H <- .C("rbfn_bias_design_matrix", as.double(t(X)), as.integer(nrow(X)), 
            as.integer(ncol(X)), as.double(t(rbfn$c)), as.double(rbfn$s), 
            as.integer(nrow(rbfn$c)), as.double(t(H)))
     H <- matrix(H[[length(H)]], nrow(X), nrow(rbfn$c) + 1, byrow = TRUE)
	H
	}


rbfn.bias.output <- function(rbfn, X) {
   rbfn.bias.design.matrix(rbfn, X) %*% rbfn$w
   }

print("rbfn.gaussian.R loaded")	
	
