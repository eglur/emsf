source("ml.R")


## PRIMEIRA LISTA DE EXERCÍCIOS ##

script.gen.points.reg.1 <- function() {
	X <- seq(-3, 1, l = 50)
 	X <- sample(X, length(X))
	y <- -12 * X + 44
 	y <- y + rnorm(length(y), 0, 5)
	
	ind <- sample(1:length(X), 20, replace = FALSE)
	X.test <- X[ind]
	y.test <- y[ind]

	X <- X[-ind]
	y <- y[-ind]


	plot(cbind(X,y))
	points(cbind(X.test, y.test), col = "RED")

	X <- matrix(X,length(X), 1)
	theta <- least.squares.params(X, y)
	lines(cbind(X, cbind(1,X) %*% theta, col = "RED"))
	print(paste("Training error: ", least.squares.batch(X,y,theta)))
	print(paste("Test error: ", least.squares.batch(X.test,y.test,theta)))

	wt(X, "./ml_files/reg_1_tr_X.dat")
	wt(y, "./ml_files/reg_1_tr_Y.dat")
	wt(X.test, "./ml_files/reg_1_ts_X.dat")
	wt(y.test, "./ml_files/reg_1_ts_Y.dat")
	}


script.gen.points.reg.2 <- function() {
	X <- as.matrix(read.table("./ml_files/q2x.dat"))
	y <- drop(as.matrix(read.table("./ml_files/q2y.dat")))

	m <- round(0.75 * nrow(X))
	ind <- sample(1:nrow(X), m, replace = FALSE)

	X.test <- X[-ind,]
	y.test <- y[-ind]
	
	X <- X[ind,]
	y <- y[ind]

	for (tau in c(0.01, 0.05, 0.1, 0.3, 0.8, 2, 10)) {
		print(paste("Tau =", tau, ":", lwr.batch(X,y,X.test,y.test,tau)))
		x11(w=4, h=4)
		lwr.plot(X,y,tau)
		}

	wt(X, "./ml_files/reg_2_tr_X.dat")
	wt(y, "./ml_files/reg_2_tr_Y.dat")
	wt(X.test, "./ml_files/reg_2_ts_X.dat")
	wt(y.test, "./ml_files/reg_2_ts_Y.dat")
	}

gen.points.gaussian <- function(m = c(50,50), mu0 = c(1,1), mu1 = c(3,3),
						  Sigma1 =  matrix(c(2,-0.8,-0.8,1), 2, 2),
						  Sigma2 = Sigma1, 
						  plot = FALSE) {
	X1 <- rmnorm(m[1], mu0, Sigma1)
	X2 <- rmnorm(m[2], mu1, Sigma2)


	if (plot) {
		max.X <- apply(rbind(X1,X2), 2, max)
		min.X <- apply(rbind(X1,X2), 2, min)

		plot(X1, col = "RED", pch = 1, xlim = c(min.X[1], max.X[1]), ylim =
				c(min.X[2], max.X[2]))
		points(X2, col = "GREEN", pch = 2)
		}

	y1 <- array(0, m[1])
	y2 <- array(1, m[2])

	list (X = rbind(X1,X2), y = c(y1,y2))
	}
	

script.gen.points.class.1 <- function() {
	# Classification Problem #1: 2 Gaussians with the same covariance matrix
	D <- gen.points.gaussian(m = c(50,50), mu0 = c(1,1), mu1 = c(3,3),
						  Sigma1 =  matrix(c(2,-0.8,-0.8,2), 2, 2),
						  Sigma2 = matrix(c(2,-0.8,-0.8,2), 2, 2),
						  plot = TRUE)
	wt(D$X, "./ml_files/class_1_tr_X.dat")
	wt(D$y, "./ml_files/class_1_tr_Y.dat")

     t <- logistic.regression.newton(D$X, D$y)
     G <- gda.params(D$X, D$y)

	D <- gen.points.gaussian(m = c(25,15), mu0 = c(1,1), mu1 = c(3,3),
						  Sigma1 =  matrix(c(2,-0.8,-0.8,2), 2, 2),
						  Sigma2 = matrix(c(2,-0.8,-0.8,2), 2, 2),
						  plot = TRUE)
	wt(D$X, "./ml_files/class_1_ts_X.dat")
	wt(D$y, "./ml_files/class_1_ts_Y.dat")

     print(paste("LG:",logistic.regression.batch(D$X, D$y, t)))
     print(paste("GDA:", gda.batch(D$X, D$y, G)))
	}


script.gen.points.class.2 <- function() {
	# Classification Problem #2: 2 Gaussians with different covariance matrices
	D <- gen.points.gaussian(m = c(50,50), mu0 = c(1,1), mu1 = c(3,3),
						  Sigma1 =  matrix(c(2,-0.8,-0.8,2), 2, 2),
						  Sigma2 = matrix(c(2,0.8,0.8,2), 2, 2),
						  plot = TRUE)
	wt(D$X, "./ml_files/class_2_tr_X.dat")
	wt(D$y, "./ml_files/class_2_tr_Y.dat")

     t <- logistic.regression.newton(D$X, D$y)
     G <- gda.params(D$X, D$y)

	D <- gen.points.gaussian(m = c(10,30), mu0 = c(1,1), mu1 = c(3,3),
						  Sigma1 =  matrix(c(2,-0.8,-0.8,2), 2, 2),
						  Sigma2 = matrix(c(2,0.8,0.8,2), 2, 2),
						  plot = TRUE)
	wt(D$X, "./ml_files/class_2_ts_X.dat")
	wt(D$y, "./ml_files/class_2_ts_Y.dat")

     print(paste("LG:",logistic.regression.batch(D$X, D$y, t)))
     print(paste("GDA:", gda.batch(D$X, D$y, G)))
	}




script.gen.points.class.3 <- function() {
	# Classification Problem #3: discrete attributes
	D <- read.table("./ml_files/kr-vs-kp.txt")

	X <- D[,1:(ncol(D)-1)]
	y <- D[,ncol(D)]

	m <- round(0.75 * nrow(X))
	
	ind <- sample(1:nrow(X), m, replace = FALSE)
	
	X.tr <- X[ind,]
	y.tr <- y[ind]

	X.ts <- X[-ind,]
	y.ts <- y[-ind]

	wt(X.tr, "./ml_files/class_3_tr_X.dat")
	wt(y.tr, "./ml_files/class_3_tr_Y.dat")
	wt(X.ts, "./ml_files/class_3_ts_X.dat")
	wt(y.ts, "./ml_files/class_3_ts_Y.dat")

	N <- naive.bayes.params(X.tr, y.tr)
	print(paste("NB training:", naive.bayes.batch(X.tr, y.tr, N)))
	print(paste("NB test:", naive.bayes.batch(X.ts, y.ts, N)))
	}


## ERRO DE TESTE NA LISTA DE EXERCÍCIOS 1 ##

erros.lista1 <- function(sobrenome.aluno, 
                  suffix = "_ts_Y.dat", 
                  prefix.reg = "reg_", 
                  prefix.class = "class_", 
                  num.reg = 2, 
                  num.class = 3, 
                  dir.arqs = "~/docs/aulas/machine_learning/lncc_barreto/dados/", 
                  dir.alunos = "~/docs/aulas/machine_learning/lncc_barreto/2014/lista_exercicios/",
                  meus.reg = c(43.53, 0.065), 
                  meus.class = c(0.05,0.15,0.12)) 
{
	for (i in 1:num.reg) {
		y <- read.table(paste(dir.arqs,prefix.reg,i,suffix, sep=""))
		erro <- mse(y,read.table(paste(dir.alunos,sobrenome.aluno,"_",
prefix.reg,i,".dat", sep="")))
		print(paste("Regressão ",i,": ", round(erro,d=3), " (",
meus.reg[i],")", sep = ""))
		}
	for (i in 1:num.class) {
		y <- read.table(paste(dir.arqs,prefix.class,i,suffix, sep=""))
		ya <- read.table(paste(dir.alunos,sobrenome.aluno,"_",
prefix.class,i,".dat", sep=""))
		ya[ya > 0.5]  <- 1
     	ya[ya <= 0.5] <- 0
		erro <- error.class(y,ya)
		print(paste("Classificacao ",i,": ", round(erro,d=3), " (",
meus.class[i],")", sep = ""))
		}
}


## ERRO DE TESTE NA LISTA DE EXERCÍCIOS 2 ##

erros.lista2 <- function(sobrenome.aluno, 
                          dir.arqs = "~/docs/aulas/machine_learning/lncc_barreto/dados/", 
                          dir.alunos = "~/docs/aulas/machine_learning/lncc_barreto/2014/lista_exercicios/",
                          meus.casa = c(41.60,31.97), 
                          meus.diag = c(0.0429,0.0429),
                          entregou = c(TRUE, TRUE, TRUE, TRUE)) 
{
	y <- read.table(paste(dir.arqs,"casas_ts_Y.dat", sep=""))

	if (entregou[1]) {
	erro <- mse(y,read.table(paste(dir.alunos,sobrenome.aluno,
				"_casas_simples.dat", sep="")))
	print(paste("Casas simples: ", round(erro,d=3), " (",
meus.casa[1],")", sep = ""))
	}

	if (entregou[2]) {
	erro <- mse(y,read.table(paste(dir.alunos,sobrenome.aluno,
				"_casas_comite.dat", sep="")))

	print(paste("Casas comite: ", round(erro,d=3), " (",
meus.casa[2],")", sep = ""))
     }


	y <- read.table(paste(dir.arqs,"diagnostico_ts_Y.dat", sep=""))

	if (entregou[3]) {
	erro <- error.class(y,read.table(paste(dir.alunos,sobrenome.aluno,
				"_diagnostico_simples.dat", sep="")))
	print(paste("Diagnostico simples: ", round(erro,d=4), " (",
meus.diag[1],")", sep = ""))
	}


	if (entregou[4]) {
	erro <- error.class(y,read.table(paste(dir.alunos,sobrenome.aluno,
				"_diagnostico_comite.dat", sep="")))

	print(paste("Diagnostico comite: ", round(erro,d=4), " (",
meus.diag[2],")", sep = ""))
	}
}


## PRIMEIRA PROVA	 ##

script.gen.points.test <- function() {
	# Classification Problem #2: 2 Gaussians with different covariance matrices
	D <- gen.points.gaussian(m = c(100,100), mu0 = c(1,1), mu1 = c(3,3),
						  Sigma1 =  matrix(c(2,-0.8,-0.8,2), 2, 2),
						  Sigma2 = matrix(c(1.5,1,1,1.5), 2, 2),
						  plot = TRUE)
	x11(h=5,w=5)
	plot(D$X[D$y==1,], xlim = c(min(D$X[,1]),max(D$X[,1])),ylim =
          c(min(D$X[,2]),max(D$X[,2])), xlab = expression(x[1]), ylab =
		expression(x[2]) )
	points(D$X[D$y==0,], col = "RED", pch = 2)
     dev.copy2eps(file =
"~/machine_learning/lncc_barreto/exercicios/fig/dists.eps")
	}


## SEGUNDA PROVA ##

questao.svm <- function(C=1e10) {
	C1 <- matrix(c(3,3.5,3.5,4,7,5,4.5,5.5,5,5), 5, 2)
	C2 <- matrix(c(5,5.5,6,6,7,7,2.5,2,3,3.5,2.5,3.5), 6, 2)
	max.x <- apply(rbind(apply(C1,2,max), apply(C2,2,max)),2,max)
	min.x <- apply(rbind(apply(C1,2,min), apply(C2,2,min)),2,min)

	X <- rbind(C1,C2)
	y <- array(1,nrow(X))
	y[1:nrow(C1)] <- -1

	x1 <- min.x[1]-2
	x2 <- max.x[1]+2

	for (c in C) {
		P <- svm.params(X,y,C=c,ker=dot.prod)

          x11(w=5,h=5)
		plot(C1, col = "BLUE", xlim=c(min.x[1]-1,max.x[1]+1),
			ylim=c(min.x[2]-1,max.x[2]+1), xlab=expression(x[1]), 
			ylab=expression(x[2]), pch=15)
		points(C2, col = "RED", pch=17)
# 		text(X+0.2, paste(1:nrow(X)))
# 	 	points(P$sv, col="BLACK", pch=8)

		w <- array(0,2)
		for (i in 1:nrow(P$sv)) {
			w <- w + P$y.alpha[i] * P$sv[i,]
			}
	
          print(w)
		print(P$b)

		y1 <- (-w[1] * x1 - P$b) / w[2]
		y2 <- (-w[1] * x2 - P$b) / w[2]
		lines(cbind(c(x1,x2),c(y1,y2)))
		
		dev.copy2eps(file=paste("questao_svm_C_Inf_solucao.eps",sep=""))
		}
	}

questao.kmeans <- function() {
# 	g1 <- matrix(c(1,2,1,4,1,1,2,4),4,2)
# 	g2 <- matrix(c(5,7,4,6),2,2)

 	g1 <- matrix(c(0,1,2,5,1,2,0,5),4,2)
 	g2 <- matrix(c(6,7,7,6),2,2)

 	c2 <- matrix(c(4,5),1,2)
 	c1 <- matrix(c(6,6),1,2)
	lim.x1 <- c(0,8)
	lim.x2 <- c(0,8)
     x11(w=5,h=5)
	plot(NULL,xlab=expression(theta[1]),ylab=expression(theta[2]),pch=0,
		ylim=lim.x2, xlim=lim.x1,cex=2, 
		xaxp=c(lim.x1[1],lim.x1[2],lim.x1[2]-lim.x1[1]),
		yaxp=c(lim.x2[1],lim.x2[2],lim.x2[2]-lim.x2[1]))
	s1 <- seq(lim.x1[1],lim.x1[2],l=lim.x1[2]-lim.x1[1]+1)
	abline(h=s1,lty="dotted",col="darkgrey")
	s2 <- seq(lim.x2[1],lim.x2[2],l=lim.x2[2]-lim.x2[1]+1)
	abline(v=s2,lty="dotted",col="darkgrey")
	
	points(g1, pch=0, cex=2)
	points(g2, pch=15, cex=2)
	points(c1, pch=17, cex=2)
	points(c2, pch=2, cex=2)


	text(c1+0.3,expression(mu[1]))
	text(c2+0.3,expression(mu[2]))
	dev.copy2eps(file="kmeans.eps")

	plot(NULL,xlab=expression(theta[1]),ylab=expression(theta[2]),pch=0,
		ylim=lim.x2, xlim=lim.x1,cex=2, 
		xaxp=c(lim.x1[1],lim.x1[2],lim.x1[2]-lim.x1[1]),
		yaxp=c(lim.x2[1],lim.x2[2],lim.x2[2]-lim.x2[1]))
	s1 <- seq(lim.x1[1],lim.x1[2],l=lim.x1[2]-lim.x1[1]+1)
	abline(h=s1,lty="dotted",col="grey")
	s2 <- seq(lim.x2[1],lim.x2[2],l=lim.x2[2]-lim.x2[1]+1)
	abline(v=s2,lty="dotted",col="grey")
 	dev.copy2eps(file="kmeans_empty.eps")
	}


questao.kmeans.solucao <- function() {
# 	g1 <- matrix(c(1,2,1,4,1,1,2,4),4,2)
# 	g2 <- matrix(c(5,7,4,6),2,2)

 	g1 <- matrix(c(0,1,2,1,2,0),3,2)
 	g2 <- matrix(c(6,7,5,7,6,5),3,2)

 	c2 <- matrix(c(2,2),1,2)
 	c1 <- matrix(c(6.5,6.5),1,2)
	lim.x1 <- c(0,8)
	lim.x2 <- c(0,8)
     x11(w=5,h=5)
	plot(NULL,xlab=expression(theta[1]),ylab=expression(theta[2]),pch=0,
		ylim=lim.x2, xlim=lim.x1,cex=2, 
		xaxp=c(lim.x1[1],lim.x1[2],lim.x1[2]-lim.x1[1]),
		yaxp=c(lim.x2[1],lim.x2[2],lim.x2[2]-lim.x2[1]))
	s1 <- seq(lim.x1[1],lim.x1[2],l=lim.x1[2]-lim.x1[1]+1)
	abline(h=s1,lty="dotted",col="darkgrey")
	s2 <- seq(lim.x2[1],lim.x2[2],l=lim.x2[2]-lim.x2[1]+1)
	abline(v=s2,lty="dotted",col="darkgrey")
	
	points(g1, pch=0, cex=2)
	points(g2, pch=15, cex=2)
	points(c1, pch=17, cex=2)
	points(c2, pch=2, cex=2)


	text(c1+0.3,expression(mu[1]))
	text(c2+0.3,expression(mu[2]))
	dev.copy2eps(file="kmeans_solucao_1.eps")

 	c2 <- matrix(c(1,1),1,2)
 	c1 <- matrix(c(6,6),1,2)
	lim.x1 <- c(0,8)
	lim.x2 <- c(0,8)
     x11(w=5,h=5)
	plot(NULL,xlab=expression(theta[1]),ylab=expression(theta[2]),pch=0,
		ylim=lim.x2, xlim=lim.x1,cex=2, 
		xaxp=c(lim.x1[1],lim.x1[2],lim.x1[2]-lim.x1[1]),
		yaxp=c(lim.x2[1],lim.x2[2],lim.x2[2]-lim.x2[1]))
	s1 <- seq(lim.x1[1],lim.x1[2],l=lim.x1[2]-lim.x1[1]+1)
	abline(h=s1,lty="dotted",col="darkgrey")
	s2 <- seq(lim.x2[1],lim.x2[2],l=lim.x2[2]-lim.x2[1]+1)
	abline(v=s2,lty="dotted",col="darkgrey")
	
	points(g1, pch=0, cex=2)
	points(g2, pch=15, cex=2)
	points(c1, pch=17, cex=2)
	points(c2, pch=2, cex=2)


	text(c1+0.3,expression(mu[1]))
	text(c2+0.3,expression(mu[2]))
	dev.copy2eps(file="kmeans_solucao_2.eps")
	}



print("ml.assignments.R loaded")