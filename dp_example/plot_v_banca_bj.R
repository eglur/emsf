D <- NULL

tmp <- read.table("v_banca_bj_1000_10_01.log")
D <- cbind(D, t(tmp))
tmp <- read.table("v_banca_bj_10000_100_01.log")
D <- cbind(D, t(tmp))
tmp <- read.table("v_banca_bj_100000_1000_01.log")
D <- cbind(D, t(tmp))
tmp <- read.table("v_banca_bj_1000000_10000_01.log")
D <- cbind(D, t(tmp))

## tmp <- read.table("v_banca_bj_10000000_100000_01.log")
## D <- cbind(D, t(tmp))
## tmp <- read.table("v_banca_bj_100000000_1000000_01.log")
## D <- cbind(D, t(tmp))
## tmp <- read.table("v_banca_bj_1000000000_10000000_01.log")
## D <- cbind(D, t(tmp))

matplot(D, t="l", main="Valor de aplicar a política da banca para diferentes tam. de aval.")
