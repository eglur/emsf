D <- NULL

tmp <- read.table("~/emsf/dp_example/log/v_cnt_bj_20_2_100000_1000000_0_100.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

tmp <- read.table("v_emsf_bj_100000_1000000_0_100_100_30_1_300.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

matplot(D, t="l")
