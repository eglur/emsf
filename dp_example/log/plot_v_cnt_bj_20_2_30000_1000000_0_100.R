D <- NULL

tmp <- read.table("v_cnt_bj_20_2_30000_1000000_0_100.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

matplot(D, t="o", main="Counting + PI utilizando 30000 lotes (avaliando em 1000000 jogos)")
grid()
