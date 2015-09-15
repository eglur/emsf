D <- NULL
legends=list()

tmp <- read.table("v_aaai_bj_203_20_2_5000_1000000_0_100_5000_1_300_10_0.1_0.1.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.1")

tmp <- read.table("v_aaai_bj_203_20_2_5000_1000000_0_100_5000_1_300_10_0.4_0.1.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.4")

tmp <- read.table("v_aaai_bj_203_20_2_5000_1000000_0_100_5000_1_300_10_1_0.1.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("1")

tmp <- read.table("v_aaai_bj_203_20_2_5000_1000000_0_100_5000_1_300_10_0.1_0.4.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.1")

tmp <- read.table("v_aaai_bj_203_20_2_5000_1000000_0_100_5000_1_300_10_0.4_0.4.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.4")

tmp <- read.table("v_aaai_bj_203_20_2_5000_1000000_0_100_5000_1_300_10_1_0.4.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("1")

tmp <- read.table("v_aaai_bj_203_20_2_5000_1000000_0_100_5000_1_300_10_0.1_1.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.1")

tmp <- read.table("v_aaai_bj_203_20_2_5000_1000000_0_100_5000_1_300_10_0.4_1.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.4")

tmp <- read.table("v_aaai_bj_203_20_2_5000_1000000_0_100_5000_1_300_10_1_1.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("1")

tmp <- read.table("v_aaai_bj_203_20_2_5000_1000000_0_100_5000_1_300_100_0.1_0.1.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.1")

tmp <- read.table("v_aaai_bj_203_20_2_5000_1000000_0_100_5000_1_300_100_0.4_0.1.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.4")

tmp <- read.table("v_aaai_bj_203_20_2_5000_1000000_0_100_5000_1_300_100_1_0.1.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("1")

tmp <- read.table("v_aaai_bj_203_20_2_5000_1000000_0_100_5000_1_300_100_0.1_0.4.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.1")

tmp <- read.table("v_aaai_bj_203_20_2_5000_1000000_0_100_5000_1_300_100_0.4_0.4.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.4")

tmp <- read.table("v_aaai_bj_203_20_2_5000_1000000_0_100_5000_1_300_100_1_0.4.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("1")

tmp <- read.table("v_aaai_bj_203_20_2_5000_1000000_0_100_5000_1_300_100_0.1_1.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.1")

tmp <- read.table("v_aaai_bj_203_20_2_5000_1000000_0_100_5000_1_300_100_0.4_1.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.4")

tmp <- read.table("v_aaai_bj_203_20_2_5000_1000000_0_100_5000_1_300_100_1_1.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("1")

matplot(D, t="o")
grid()
