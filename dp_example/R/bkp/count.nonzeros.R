dyn.load(paste("./c/count_nonzeros", .Platform$dynlib.ext, sep = "")) #melhorar

count.nonzeros <- function(D, K) {
   res <- 0
   a <- array(0,1)
   if (is.null(dim(D))) dim(D) <- c(1,length(D))
   if (is.null(dim(K))) dim(K) <- c(length(K),1)
   if (ncol(D) != nrow(K)) print("ncol(D) must be equal to nrow(K)")
   else {
     R <- .C("count_nonzeros", as.double(t(D)), as.double(K),
            as.integer(nrow(D)), as.integer(ncol(D)), out = as.integer(a))
     res <- R$out
     }
   res  
	}


print("count.nonzeros.R loaded")	
	
