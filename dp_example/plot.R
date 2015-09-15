epsilons=c("0.01", "0.05", "0.10", "0.20", "0.3", "1")
pref="/home/rafaelbeirigo/emsf/dp_example/v_cnt_bj_20_2_30000_1000000_0_100_"
suf="_01.log"

plot(read.table("/home/rafaelbeirigo/emsf/dp_example/v_cnt_bj_20_2_30000_1000000_0_100_1_01.log"))

## D <- NULL
## plot()
## for (epsilon in epsilons) {
##     lines(read.table(paste(sep="", pref, epsilon, suf)))
##     tmp <- read.table(paste(sep="", pref, epsilon, suf))
##     D <- cbind(D, tmp)
## }


