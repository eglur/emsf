D <- NULL

tmp <- read.table("v_emsf_bj_100000_1000000_0_100_100_30_1_300.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

matplot(D, t="l", main="EMSF-SK + PISF utilizando 100000 lotes (avaliando em 1000000 jogos)")
grid()

