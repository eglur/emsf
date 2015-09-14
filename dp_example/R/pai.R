prob.acerto <- function(n,k, p = 0.5) {
# veja entrada "Binomial distribution" na Wikipedia
	choose(n,k) * p^k * (1-p)^(n-k)
	}
	

calc.probs <- function(num.questoes, num.opcoes=2) {
	n <- num.questoes
	p <- 1/num.opcoes
	probs <- array(0,n)
	for (i in 1:n) {
		probs[i] <- prob.acerto(n,i,p)
		}
	for (i in (n-1):1) {
		probs[i] <- probs[i] + probs[i+1]
		}
	probs
	}
	
calc.notas <- function(num.questoes, num.opcoes=2, nota.max=num.questoes) {
	probs <- calc.probs(num.questoes, num.opcoes)
	delta <- probs[length(probs)]
	(1-probs+delta) * 1:nota.max
	}
	

print("pai.R loaded")