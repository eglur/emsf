dyn.load(paste("./c/pole_double_center", .Platform$dynlib.ext, sep = ""))

pole.double.center.transition <- function(s, a, noise.sd = 0) {
   r <- 0
   g <- 0
   a <- rnorm(1, a, noise.sd)
   R <- .C("pole_double_center_transition", as.double(t(s)), as.double(a),
         as.double(r),as.integer(g))
   g <- FALSE
   if (R[[4]] == 1) g <- TRUE
   list(s = R[[1]], r = R[[3]], g = g)
   }


print("pole.double.center.R loaded")