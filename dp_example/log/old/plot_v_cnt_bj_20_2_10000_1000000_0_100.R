D <- NULL

tmp <- read.table("v_cnt_bj_20_2_10000_1000000_0_100.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

matplot(D, t="o", log="y", main="n=203; sr=20; maxit=; D, K aleat.")

