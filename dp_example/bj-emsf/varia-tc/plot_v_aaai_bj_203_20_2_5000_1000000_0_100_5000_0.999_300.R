D <- NULL
legends=list()

tmp <- read.table("v_aaai_bj_203_20_2_5000_1000000_0_100_5000_0.999_300_5_0.15_0.2.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.15")

tmp <- read.table("v_aaai_bj_203_20_2_5000_1000000_0_100_5000_0.999_300_5_0.15_0.4.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.15")

tmp <- read.table("v_aaai_bj_203_20_2_5000_1000000_0_100_5000_0.999_300_5_0.15_0.6.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.15")

tmp <- read.table("v_aaai_bj_203_20_2_5000_1000000_0_100_5000_0.999_300_5_0.15_0.8.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.15")

tmp <- read.table("v_aaai_bj_203_20_2_5000_1000000_0_100_5000_0.999_300_5_0.15_1.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.15")

tmp <- read.table("v_aaai_bj_203_20_2_5000_1000000_0_100_5000_0.999_300_10_0.15_0.2.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.15")

tmp <- read.table("v_aaai_bj_203_20_2_5000_1000000_0_100_5000_0.999_300_10_0.15_0.4.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.15")

tmp <- read.table("v_aaai_bj_203_20_2_5000_1000000_0_100_5000_0.999_300_10_0.15_0.6.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.15")

tmp <- read.table("v_aaai_bj_203_20_2_5000_1000000_0_100_5000_0.999_300_10_0.15_0.8.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.15")

tmp <- read.table("v_aaai_bj_203_20_2_5000_1000000_0_100_5000_0.999_300_10_0.15_1.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.15")

matplot(D, t="o")
grid()
