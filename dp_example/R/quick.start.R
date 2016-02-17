source("blackjack.experiments.R")
source("emsf.experiments.R")

## Figure 3
emsf.experiment.fixed.tc(n=100, alpha.dir=0.5, m=10, tc=600,  alphas=c(0.1, 0.3, 0.7), max.it=25*10^3, na=1, num.avg=50)
emsf.experiment.fixed.tc(n=100, alpha.dir=0.5, m=10, tc=900,  alphas=c(0.1, 0.3, 0.7), max.it=25*10^3, na=1, num.avg=50)

emsf.experiment.plot.results.aaai(save=TRUE, ylim=NULL, cex=1.5, width=8.75, height=7)

## Figure 4
bj.ml.experiment(num.episodes=1e4, epsilon=0.15, tc=100, num.avg=100, ne=1e6)

bj.qlearning.experiment(alpha=0.10, num.episodes=1e4, epsilon=0.15, tc=100, num.avg=100, ne=1e6)

bj.emsf.comp.experiment(m=10, alpha=1.0, tcc=100, num.episodes=1e4, epsilon=0.15, tc=100, num.avg=100, ne=1e6)
bj.emsf.comp.experiment(m=10, alpha=0.5, tcc=100, num.episodes=1e4, epsilon=0.15, tc=100, num.avg=100, ne=1e6)
bj.emsf.comp.experiment(m=10, alpha=0.1, tcc=100, num.episodes=1e4, epsilon=0.15, tc=100, num.avg=100, ne=1e6)

plot.results.aaai()
