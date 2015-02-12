D <- NULL

tmp <- read.table("v_cnt_bj_20_2_30000_1000000_0_100.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

tmp <- read.table("v_emsf_bj_30000_1000000_0_100_100_30_1_300_10.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

tmp <- read.table("v_emsf_bj_30000_1000000_0_100_100_30_1_300_50.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

tmp <- read.table("v_emsf_bj_30000_1000000_0_100_100_30_1_300_100.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

matplot(D, t="l", main="v_emsf_bj_30000_1000000_0_100_100_30_1_300_10")
grid()

