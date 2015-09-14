
least.squares.params <- function(X, y) {
	X <- as.matrix(cbind(1,X))
	y <- drop(as.matrix(y))
	solve(t(X) %*% X) %*% t(X) %*% y
	}


least.squares.batch <- function(X, theta) {
	as.matrix(cbind(1,X)) %*% theta
	}


lwr <- function(X, y, x, tau = 0.1) {
	X <- cbind(1,X)
	x <- c(1,x)
	W <- matrix(rep(x, nrow(X)), nrow(X), ncol(X), byrow = TRUE)
	W <- apply((X-W)^2, 1, sum)
	W <- 0.5 * exp(-W/(2*tau^2))
	W <- diag(W)
	theta <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y
	t(theta) %*% x
	}                                                                          
                             


lwr.batch <- function(X, y, X2, tau) {
	X <- as.matrix(X)
	y <- drop(as.matrix(y))
	X2 <- as.matrix(X2)
     h <- array(0, nrow(X2))

	for (i in 1:nrow(X2)) h[i] <- lwr(X,y,X2[i,],tau)
	h
	}


lwr.best.tau <- function(X,y,Vx,vy,tau=seq(0.1,0.9,by=0.1)) {
# X are the input points of the training set
# y are the output points of the training set
# Vx are the input points of the validation set
# vy are the output points of the validation set
	errors <- array(0,length(tau))
	for (t in 1:length(tau)) {
		h   <- lwr.batch(X,y,Vx,tau[t])
		errors[t] <- mse(h,vy)
		}
	tau[which.min(errors)]
	}
		
		
		
lwr.plot <- function(X, y, tau, length = 100) {
	X <- as.matrix(X)
	y <- drop(as.matrix(y))

	X.plot <- matrix(seq(min(X), max(X), l = length), length, 1)
     y.plot <- array(0, nrow(X.plot))
	
	for (i in 1:nrow(X.plot)) y.plot[i] <- lwr(X,y,X.plot[i,],tau)

	plot(cbind(X.plot,y.plot), col = "RED", t="l")
	points(cbind(X,y))
	}
	


logistic <- function(x) {
	1/(1+exp(-x))
	}


logistic.regression.grad <- function(X, y, alpha = 1e-4, 
                                     epsilon = 1e-10, plot = FALSE) {
# X = input variables
# y = outpu variable (0 or 1)
# alpha = learning rate
# epsilon = maximum diference between to successive estimates of theta
# plot = logical; should the data be plotted?
	
	X <- as.matrix(cbind(1,X))
	y <- as.matrix(y)
	
	m <- nrow(X)
	n <- ncol(X)

	theta <- runif(n,-1e-10,1e-10)
 
	dif.theta <- Inf
	while (dif.theta > epsilon) {
		theta.prime <- theta
		grad <- array(0,n)
		error <- 0

		h <- logistic(theta%*%t(X))
		for (i in 1:m) {
			grad <- grad + (y[i,] - h[i]) * X[i,]
			}
		theta <- theta + alpha * grad
		dif.theta <- sum((theta-theta.prime)^2)
		
# 		Computes the log-likelihood
# 		ll <- 0
# 		for (i in 1:m) {
# 			ll <- ll + y[i] * log(h[i]) + (1-y[i]) * log(1 - h[i])
# 			}
#           print(ll)
          if (plot) {
			plot(X[y==0,2:3], pch=1, col="RED")
			points(X[y==1,2:3], pch=2, col = "BLUE")
			y1 <- X[1,2]
			y2 <- (-y1 * theta[2] - theta[1]) / theta[3]
			z1 <- X[2,2]
			z2 <- (-z1 * theta[2] - theta[1]) / theta[3]
			ca <- (z2-y2) / (z1-y1)
			r <- ca * (X[,2] - y1) +  y2
			lines(cbind(X[,2],r), col = "BLACK")
			save <- readline()
			if (save == "y") dev.copy2eps(file = "~/class.eps")
			}
		}
	theta
	}


logistic.regression.newton <- function(X, y, alpha = 1e-4, epsilon = 1e-10,
                                       plot = FALSE) {
# X = input variables
# y = outpu variable (0 or 1)
# alpha = learning rate
# epsilon = maximum diference between to successive estimates of theta
# plot = logical; should the data be plotted?
	
	X <- as.matrix(cbind(1,X))
	y <- as.matrix(y)
	
	m <- nrow(X)
	n <- ncol(X)

	theta <- runif(n,-1e-10,1e-10)
	dif.theta <- Inf
 
	dif.theta <- Inf
	while (dif.theta > epsilon) {
		theta.prime <- theta
		grad <- array(0,n)
		H <- matrix(0,n,n)
		error <- 0

		h <- logistic(theta%*%t(X))
		for (i in 1:m) {
			grad <- grad + (y[i,] - h[i]) * X[i,]
			H <- H - (X[i,] %o% X[i,]) * h[i] * (1 - h[i])
			}
		theta <- drop(theta - solve(H) %*% grad)
	
		dif.theta <- sum((theta-theta.prime)^2)

# 		Computes the log-likelihood
# 		ll <- 0
# 		for (i in 1:m) {
# 			ll <- ll + y[i] * log(h[i]) + (1-y[i]) * log(1 - h[i])
# 			}
#           print(ll)

		}

          if (plot) {
			max.x <- apply(X,2,max)
			min.x <- apply(X,2,min)
			plot(X[y==0,2:3], pch=1, col="RED", xlim=c(min.x[1],max.x[2]),
				ylim = c(min.x[2], max.x[2]))
			points(X[y==1,2:3], pch=2, col = "BLUE")
			y1 <- X[1,2]
			y2 <- (-y1 * theta[2] - theta[1]) / theta[3]
			z1 <- X[2,2]
			z2 <- (-z1 * theta[2] - theta[1]) / theta[3]
			ca <- (z2-y2) / (z1-y1)
			r <- ca * (X[,2] - y1) +  y2
			lines(cbind(X[,2],r), col = "BLACK")
			save <- readline()
			if (save =="y") dev.copy2eps(file = "~/class.eps")
			}

	theta
	}


logistic.regression.params <- logistic.regression.newton

logistic.regression.class <- function(x, theta) {
	o <- logistic(c(1,x) %*% theta)
	c <- 0
	if (o > 0.5) c <- 1
	list(c = c, o = o)
	}

logistic.regression.batch <- function(X, theta) {
	sum <- 0
	X <- as.matrix(X)
	h <- array(0,nrow(X))
	for (i in 1:nrow(X)) {
		h[i] <- logistic.regression.class(X[i,], theta)$c
		}
	h
	}


gda.params <- function(X,y) {
	X <- as.matrix(X)
	y <- as.matrix(y)

	m <- nrow(X)
	n <- ncol(X)
	
	phi.y <- sum(y==1) / length(y)
	
	X0 <- X[y==0,]
	mu0 <- apply(X0,2,sum) / nrow(X0)

	X1 <- X[y==1,]
	mu1 <- apply(X1,2,sum) / nrow(X1)
	
	Sigma <- matrix(0, n, n)
	D <- X0 - matrix(mu0, nrow(X0), n, byrow = TRUE)
	for (i in 1:nrow(X0)) Sigma <- Sigma + D[i,] %o% D[i,]
	
	D <- X1 - matrix(mu1, nrow(X1), n, byrow = TRUE)
	for (i in 1:nrow(X1)) Sigma <- Sigma + D[i,] %o% D[i,]
	
	Sigma <- Sigma / m
	Sigma.inv <- solve(Sigma)	
	list(phi.y = phi.y, mu0 = mu0, mu1 = mu1, Sigma.inv = Sigma.inv)
	}
	

gda.class <- function(x, P) {
	p0 <- exp(-(x-P$mu0)%*%P$Sigma.inv%*%(x-P$mu0)) * (1-P$phi.y)
	p1 <- exp(-(x-P$mu1)%*%P$Sigma.inv%*%(x-P$mu1)) * P$phi.y
	c <- 0
	if (p1 > p0) c <- 1
	list(c = c, p0 = p0, p1 = p1)
	}


gda.batch <- function(X, P) {
	X <- as.matrix(X)
	h <- array(nrow(X))
	for (i in 1:nrow(X)) {
		h[i] = gda.class(X[i,], P)$c
		}
	h
	}
	
	
	
naive.bayes.params <- function(X,y) {
	X <- as.matrix(X)
	y <- as.matrix(y)

	m <- nrow(X)
	n <- ncol(X)
	
	phi.y <- sum(y==1) / length(y)
	
	d <- -Inf # maximum number of distinct values
	for (i in 1:n) {
		v <- length(unique(X[,i]))
		if (v > d) d <- v
		}

	Phi <- array(0,c(2,n,d))
	
	for (i in 1:2){
		X2 <- X[y==(i-1),]
		for (j in 1:n) {
			for (k in 1:d) {
				Phi[i,j,k] <- sum(X2[,j] == k)
				}
			}
		Phi[i,,] <- (Phi[i,,] + 1)/ (nrow(X2) + d)
		}
	
	list(phi.y = phi.y, Phi = Phi)
	}


naive.bayes.class <- function(x,P) {
	p <- c(1,1)
	for (c in 1:2) {
		for (i in 1:length(x)) {
			p[c] <- p[c] * P$Phi[c,i,x[i]]
			}
		}
		
     p[1] <- p[1] * (1 - P$phi.y)
     p[2] <- p[2] * P$phi.y
	
	c <- 0
	if (p[2] > p[1]) c <- 1
	list (c = c, p0 = p[1], p1 = p[2])
	}


naive.bayes.batch <- function(X, P) {
	X <- as.matrix(X)
	h <- array(0, nrow(X))
	for (i in 1:nrow(X)) {
		h[i] <- naive.bayes.class(X[i,], P)$c
		}
	h
	}


naive.bayes.params.latex <- function(P) {
	for (j in 1:dim(P$Phi)[2]) {
		for (i in 1:dim(P$Phi)[1]) {
			for (k in 1:dim(P$Phi)[3]) {
				cat(paste(round(P$Phi[i,j,k],d=3), "& "))
				}
			}
		cat(" \\ ", fill=TRUE)
	}
}



back.prop.params <- function(X, y, p = 10, alpha = c(1e-4, 1e-3),
                    epsilon = 1e-10, max.iter = 1e4, plot = FALSE, 
				verbose = FALSE) {
# X = input variables
# y = outpu variable (0 or 1)
# p = number of hidden units
# alpha = learning rate
# epsilon = maximum diference between to successive estimates of theta
# max.iter = total maximum number of iterations
# plot = logical; should the data be plotted?
	
	X <- as.matrix(cbind(1,X))
	y <- as.matrix(y)
	
	m <- nrow(X)
	n <- ncol(X)

	theta1 <- matrix(runif(n*p,-10,10), n, p)
	theta2 <- runif(p,-10,10)

	dif.theta <- Inf
	it <- 0
	while (dif.theta > epsilon && it < max.iter) {

		theta1.prime <- theta1
		theta2.prime <- theta2

		# Propagation
		phi <- tanh(X %*% theta1)
		h <- phi %*% theta2

		# Backpropagation
		e <- h - y

		U <- t(matrix(rep(e, n), m, n) * X) %*% (1 - phi^2)
		U <- matrix(rep(theta2, n), n, p, byrow = TRUE) * U
 		theta1 <- theta1 - alpha[1] * U
		theta2 <- theta2 - alpha[2] * t(phi) %*% e
		
		dif.theta <- sum((theta1 - theta1.prime)^2) + 
					sum((theta2 -theta2.prime)^2)
		if (plot) {
			plot(cbind(X[,-1],y))
			o <- order(X[,-1])
			lines(cbind(X[o,-1],h[o]),col="RED")
			}
		if (verbose) print(dif.theta)
		it <- it + 1
		}
	list(theta1 = theta1, theta2 = theta2)     
	}
		

back.prop.early.stop.params <- function(X, y, p = 10, alpha = c(1e-4, 1e-3),
                    epsilon = 1e-10, perc.validation = 0.3,
				iter.validation = 500, max.iter = 1e4,
				plot = FALSE, verbose = FALSE) {
# X = input variables
# y = outpu variable (0 or 1)
# p = number of hidden units
# alpha = learning rate
# epsilon = maximum diference between to successive estimates of theta
# perc.validation = percentage of the data used as a validation set
# iter.validation = number of iterations with no improvement on the validation set
# max.iter = maximum total number of iterarions
# plot = logical; should the data be plotted?
	
	X <- as.matrix(cbind(1,X))
	y <- as.matrix(y)



	nv <- round(perc.validation * nrow(X))
	ind <- sample(1:nrow(X), nv, replace = FALSE)
	Xv <- as.matrix(X[ind,])
	yv <- as.matrix(y[ind,])
	
	X <- as.matrix(X[-ind,])
	y <- as.matrix(y[-ind])

	m <- nrow(X)
	n <- ncol(X)

	theta1 <- matrix(runif(n*p,-1,1), n, p)
	theta2 <- runif(p,-1,1)

	min.error <- Inf
	best.theta1 <- NULL
	best.theta1 <- NULL
	it <- 0
	no.improvements <- 0
	hist.error <- NULL
	hist.error.val <- NULL
	while (no.improvements < iter.validation && it < max.iter) {
		# Propagation
		phi <- tanh(X %*% theta1)
		h <- phi %*% theta2

		# Backpropagation
		e <- h - y
		hist.error <- c(hist.error,sum(e^2))
		
		U <- t(matrix(rep(e, n), m, n) * X) %*% (1 - phi^2)
		U <- matrix(rep(theta2, n), n, p, byrow = TRUE) * U
 		theta1 <- theta1 - alpha[1] * U
		theta2 <- theta2 - alpha[2] * t(phi) %*% e

		# Propagation validation set 
		phi <- tanh(Xv %*% theta1)
		h <- phi %*% theta2
		ev <- sum((h - yv)^2)
		hist.error.val <- c(hist.error.val,ev)

		if (ev < min.error) {
			min.error <- ev
			best.theta1 <- theta1
			best.theta2 <- theta2
			no.improvements <- 0
			}
		else no.improvements <- no.improvements + 1
		it <- it + 1
		
		if (verbose) print(ev)

		if (plot) {
			matplot(cbind(hist.error, hist.error.val),t="l")
			}
		}
	list(theta1 = best.theta1, theta2 = best.theta2)     
	}
		
back.prop.out <- function(x,P) {
	x <- as.matrix(c(1,x))
	phi <- tanh(t(x) %*% P$theta1)
	h <- phi %*% P$theta2
	}


back.prop.select.model <- function(X, y, ps, k = 10, verbose=TRUE,
		training.algorithm = back.prop.params, 
		out.function = back.prop.out, class = FALSE, ...) {
	p <- NULL
	min.error <- Inf
	for (i in 1:length(ps)) {
		error <- k.fold.cv(k,X,y,training.algorithm,
				out.function = back.prop.out,class=class,p=ps[i],...)
		if (error < min.error) {
			p <- ps[i]
			min.error <- error
			}
		if (verbose) print(paste("Hidden units:", ps[i], error))
		}
	list(p=p, error = min.error)
	}


back.prop.class.params <- function(X, y, p = 10, alpha = c(1e-4, 1e-3),
                    epsilon = 1e-10, compute.class.err = FALSE) {
# X = input variables
# y = outpu variable (0 or 1)
# p = number of hidden units
# alpha = learning rate
# epsilon = maximum diference between to successive estimates of theta
# plot = logical; should the data be plotted?
	
	X <- as.matrix(cbind(1,X))
	y <- as.matrix(y)
	
	m <- nrow(X)
	n <- ncol(X)

	theta1 <- matrix(runif(n*p,-1,1), n, p)
	theta2 <- runif(p,-1,1)

	dif.theta <- Inf
	while (dif.theta > epsilon) {
		theta1.prime <- theta1
		theta2.prime <- theta2

		# Propagation
		phi <- logistic(X %*% theta1)
		h <- logistic(phi %*% theta2)

		# Backpropagation
		e <- h - y

		U <- t(matrix(rep(e, n), m, n) * X) %*% (phi*(1-phi))
		U <- matrix(rep(theta2, n), n, p, byrow = TRUE) * U
 		theta1 <- theta1 - alpha[1] * U
		theta2 <- theta2 - alpha[2] * t(phi) %*% e
		
		dif.theta <- sum((theta1 - theta1.prime)^2) + 
					sum((theta2 -theta2.prime)^2)
		if (compute.class.err) {
			print(paste("Class. error:", sum(abs((h > 0.5) - y))))
			}
		print(dif.theta)
		}
	list(theta1 = theta1, theta2 = theta2)     
		}


dot.prod <- function(v1,v2, params) v1%*%v2

gaussian <- function(v1,v2, params) exp(-sum((v1-v2)^2) / params$sigma) 

poli <- function(v1,v2,params) (v1%*%v2+params$c)^params$d

svm.params <- function(X, y, C = 5, ker = dot.prod, ker.params = NULL, sigf =7){
# y[i] must be -1 or 1
# EX. ker.params = list(sigma=0.1) para o kernel gaussiano

 	X <- as.matrix(X)
#   	X <- scale(X)
 	y <- array(as.matrix(y))
	y[y==0] <- -1 # just to make sure
 
 	m <- nrow(X)

 	K <- matrix(0,m,m)

 	for (i in 1:m) {
 		for (j in i:m) {
 			K[i,j] <- ker(X[i,],X[j,], ker.params) * y[i] * y[j]
 			K[j,i] <- K[i,j]
 			}
 		}

	c <- matrix(rep(-1,m))
	A <- t(y)
	b <- 0
     l <- matrix(rep(0,m))
	u <- matrix(rep(C,m))
	r <- 0

	alpha <- ipop(c, K, A, b, l, u, r, sigf)
	alpha <- primal(alpha)

 	zero <- 10^(-sigf+2) ## this seems to give the same number of SV as ksvm
 	ind <- alpha > zero

 	sv <- as.matrix(X[ind,])
 	y.alpha <- alpha[ind] * y[ind]
 
  	ind <- ind & (alpha <  (C - zero)) 
 	lv <- as.matrix(X[ind,])
  	ly <- y[ind]
  
  	p <- length(ly)
  	b <- 0							
  	for (i in 1:p) {
  		b <- b + ly[i]
  		for (j in 1:nrow(sv)) {
  			b <- b - y.alpha[j] * ker(lv[i,], sv[j,], ker.params)
  			}
  		}
  	b <- b / p
   	list(sv = sv, y.alpha = y.alpha,  b=b, ker = ker, ker.params=ker.params)
	## So far I know that sv, alpha, and b are idential to ksvm's
	}

svm.out <- function(x,P) {
	o <- 0
	for (i in 1:nrow(P$sv)) {
           o <- o + P$y.alpha[i] * P$ker(P$sv[i,],x, P$ker.params)
		}
  	o + P$b
	}


svm.class <- function(x, P) {
	o <- 0
	for (i in 1:nrow(P$sv)) {
           o <- o + P$y.alpha[i] * P$ker(P$sv[i,],x, P$ker.params)
		}
  	o <- o + P$b
	c <- -1                       
	if (o > 0) c <- 1
	list(c=c)
	}
	
svm.batch <- function(X, P) {  ## can be one function for all the models
	X <- as.matrix(X)
#    	X <- scale(X)

	h <- array(0, nrow(X))
	for (i in 1:nrow(X)) {
		h[i] <- svm.class(X[i,], P)$c
		}
	h
	}

svm.select.model.dot.prod <- function(X, y, Cs , k = 10, verbose=TRUE) {
	C <- NULL
	min.error <- Inf
	for (i in 1:length(Cs)) {
		error <- k.fold.cv(k,X,y,svm.params,svm.class,class=TRUE,
			ker=dot.prod, ker.params = NULL)
		if (error < min.error) {
			C <- Cs[i]
			min.error <- error
			}
		if (verbose) print(paste("C", Cs[i], error))
		}
	list(C=C, error = min.error)
	}


svm.select.model.gaussian <- function(X, y, Cs , sigmas, k = 10, verbose=TRUE) {
	C <- NULL
	sigma <- NULL
	min.error <- Inf
	for (i in 1:length(Cs)) {
		for (j in 1:length(sigmas)) {
			error <- k.fold.cv(k,X,y,svm.params,svm.class,class=TRUE,
				ker=gaussian, ker.params = list(sigma=sigmas[j]))
			if (error < min.error) {
				C <- Cs[i]
				sigma <- sigmas[j]
				min.error <- error
				}
			if (verbose) print(paste("C", Cs[i], "sigma", sigmas[j], error))
			}
		}
	list(C=C, sigma = sigma, error = min.error)
	}
			

svm.select.model.poli <- function(X, y, Cs , cs, ds, k = 10, verbose=TRUE) {
	C <- NULL
	c <- NULL
	d <- NULL
	min.error <- Inf
	for (i in 1:length(Cs)) {
		for (j in 1:length(ds)) {
			for (l in 1:length(cs)) {
				error <- k.fold.cv(k,X,y,svm.params,svm.class,class=TRUE,
					ker=poli, ker.params = list(d=ds[j],c=cs[l]))
				if (error < min.error) {
					C <- Cs[i]
					d <- ds[j]
					c <- cs[l]
					min.error <- error
					}
				if (verbose) print(paste("C", Cs[i], "d", ds[j],"c", cs[l],
										error))
				}
			}
		}
	list(C=C, d = d, c = c, error = min.error)
	}
			
		

class.batch <- function(X, P, class.function) {
	X <- as.matrix(X)
	h <- array(0, nrow(X))
	for (i in 1:nrow(X)) {
		h[i] <- class.function(X[i,], P)$c
		}
	h
	}

regression.batch <- function(X, P, out.function) {
	X <- as.matrix(X)
	h <- array(0, nrow(X))
	for (i in 1:nrow(X)) {
		h[i] <- out.function(X[i,], P)
		}
	h
	}

k.fold.cv <- function(k,X,y,training.algorithm, out.function,class= TRUE, ...){ 
	p <- nrow(X) %/% k
	Ind <- matrix(0, k, p)

	av <- 1:nrow(X) ## indices still available
	for (i in 1:nrow(Ind)) {
		j <- sample(1:length(av), ncol(Ind), replace = FALSE)
		Ind[i,] <- av[j]
		av <- av[-j]
		}

	error <- 0
	for (i in 1:nrow(Ind)) {
		Xts <- X[Ind[i,],]
		Xtr <- X[-Ind[i,],]
		yts <- y[Ind[i,]]
		ytr <- y[-Ind[i,]]
		P <- training.algorithm(X=Xtr, y=ytr, ...)
		if (class) {
			h <- class.batch(X=Xts, P=P, class.function = out.function)
			error <- error + error.class(h,yts)
			}
		else {
			h <- regression.batch(X=Xts, P=P, out.function = out.function)
			error <- error + mse(h,yts)
			}

		}

	error <- error / nrow(Ind)

	error
	}
	

mutual.information <- function(X,y) {
	mi <- array(0, ncol(X))

	vals.y <- unique(y)
	py <- array(0, length(vals.y))
	for (i in 1:length(py)) py[i] <- sum(y == vals.y[i]) / length(y)

	for (i in 1:ncol(X)) {
		vals.x <- unique(X[,i])
		for (vx in vals.x) {
			px <- sum(X[,i]==vx) / nrow(X)
			for (j in 1:length(vals.y)) {
				pxy <- sum(X[y==vals.y[j],i] == vx) / nrow(X)
				if (pxy != 0) mi[i] <- mi[i] + pxy * log(pxy/(px * py[j]))
				}
			}
		}
	mi
	}
				
		

select.features.mi <- function(X,y,training.algorithm, out.function, 
					k = 5, class=TRUE, verbose = TRUE, latex = FALSE, ...){
	mi <- mutual.information(X,y)
	o <- order(mi, decreasing = TRUE)
	min.error <- Inf

     atr <- NULL
	best.atr <- NULL
	for (i in 1:length(o)) {
		atr <- c(atr, o[i])
		Xtr <- as.matrix(X[,atr])
		error <- k.fold.cv(k,Xtr,y,training.algorithm, out.function,class,...)
		if (error < min.error) {
			min.error <- error
			best.atr <- atr
			}
		if (verbose) {
			if (!latex) print(c(atr,round(error,d=3)))
			else {
				for (i in 1:length(atr)) {
					cat(atr[i])
					if (i < length(atr)) cat(", ")
					}
				cat(" & ")
				cat(round(error,d=3))
				cat(" \\") 
				cat("\\ \n")
				}
			}

		}
	atr
	}


select.features.forward <- function(X,y,training.algorithm, out.function, 
					k = 5, class=TRUE, verbose = TRUE, latex = FALSE, ...){
	min.error <- Inf
	remaining.atr <- 1:ncol(X)
     selected.atr <- NULL
	best.atr <- NULL
	while (length(remaining.atr) > 0) {
		min.e <- Inf
		min.i <- Inf
		for (i in 1:length(remaining.atr)) {
			sa <- c(selected.atr, remaining.atr[i])
			Xtr <- as.matrix(X[,sa])
			error <- k.fold.cv(k,Xtr,y,training.algorithm,
					out.function,class,...)
			if (error < min.e) {
				min.e <- error
				min.i <- i
				}
			if (verbose) {
				if (!latex) print(c(sa,round(error,d=2)))
				else {
					for (i in 1:length(sa)) {
						cat(sa[i])
						if (i < length(sa)) cat(", ")
						}
					cat(" & ")
					cat(round(error,2))
					cat(" \\") 
					cat("\\ \n")
					}
				}
			}
		selected.atr <- c(selected.atr,remaining.atr[min.i])
		remaining.atr <- remaining.atr[-min.i]
		if (min.e < min.error) {
			min.error <- min.e
			best.atr <- selected.atr
			}
		}
			
	best.atr
	}



bagging.params <- function(k, X, y, training.algorithm,...) {
	lP <- list()
	for (i in 1:k) {
		ind <- sample(1:nrow(X), nrow(X), replace = TRUE)
		Xtr <- X[ind,]
		ytr <- y[ind,]
		P <- training.algorithm(Xtr,ytr,...)
		lP <- c(lP, list(P))
		}
	lP
	}


bagging.out <- function(X, lP, out.function, class = FALSE) {
	hc <- matrix(0, nrow(X), length(lP))
	for (i in 1:length(lP)) {
		if (class) hc[,i] <- class.batch(X, lP[[i]],out.function)
		else hc[,i] <- regression.batch(X, lP[[i]],out.function)
		}
	h <- NULL
	if (class) {
		# find the classes of the problem
		classes <- NULL
		g <- apply(hc, 1, unique)
		for (i in 1:length(g)) classes <-unique(c(classes,g[[i]]))
		G <- matrix(0, nrow(X), length(classes))
		for (i in 1:ncol(G)) {
			G[,i] <- apply(hc==classes[i], 1, sum)
			}
		G <- apply(G,1,which.max)
		h <- array(0,length(G))
		for (i in 1:length(G)) h[i] <- classes[G[i]]
		}
	else h <- apply(hc,1,mean)
	h
	}
	
	
	



mse <- function(h, y) {
	h <- as.array(as.matrix(h))
	y <- as.array(as.matrix(y))
	sum((h-y)^2) / length(h)
	}

error.class <- function(h,y) {
	h <- as.array(as.matrix(h))
	y <- as.array(as.matrix(y))
	sum(h != y) / length(h)
	}



		
print("ml.R loaded")

