D <- NULL

tmp <- read.table("v_emsf_bj_30000_1000000_0_100_100_30_1_300_75.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

matplot(D, t="l", main="v_emsf_bj_30000_1000000_0_100_100_30_1_300_75")
legend("right", c("test2", "test2"), t="l")

grid()

