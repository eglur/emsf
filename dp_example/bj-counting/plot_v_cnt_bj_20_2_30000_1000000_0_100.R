D <- NULL
legends=list()

tmp <- read.table("v_cnt_bj_20_2_30000_1000000_0_100_0.01.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.01")

tmp <- read.table("v_cnt_bj_20_2_30000_1000000_0_100_0.05.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.05")

tmp <- read.table("v_cnt_bj_20_2_30000_1000000_0_100_0.1.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.1")

tmp <- read.table("v_cnt_bj_20_2_30000_1000000_0_100_0.2.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.2")

tmp <- read.table("v_cnt_bj_20_2_30000_1000000_0_100_0.3.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.3")

tmp <- read.table("v_cnt_bj_20_2_30000_1000000_0_100_0.4.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.4")

tmp <- read.table("v_cnt_bj_20_2_30000_1000000_0_100_0.5.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.5")

## tmp <- read.table("v_cnt_bj_20_2_30000_1000000_0_100_0.6.log")
## tmp <- apply(tmp, 2, mean)
## D <- cbind(D, tmp)
## legends[[length(legends)+1]]=paste("0.6")

tmp <- read.table("v_cnt_bj_20_2_30000_1000000_0_100_0.7.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.7")

tmp <- read.table("v_cnt_bj_20_2_30000_1000000_0_100_0.8.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.8")

tmp <- read.table("v_cnt_bj_20_2_30000_1000000_0_100_0.9.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("0.9")

tmp <- read.table("v_cnt_bj_20_2_30000_1000000_0_100_1.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)
legends[[length(legends)+1]]=paste("1")

matplot(D, type="l")
grid()
legend("bottomright", legend=legends)
