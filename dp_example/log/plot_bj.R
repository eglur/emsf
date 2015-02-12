D <- NULL

tmp <- read.table("v_emsf_bj_30000_1000000_0_100_100_30_1_300_10.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

tmp <- read.table("v_emsf_bj_30000_1000000_0_100_100_30_1_300_50.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

tmp <- read.table("v_emsf_bj_30000_1000000_0_100_100_30_1_300_100.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

tmp <- read.table("v_cnt_bj_20_2_30000_1000000_0_100.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

D <- cbind(D, -0.079234) # Estrategia da banca

par(mai=c(0.825, 0.85, 0.02,0.015)) # Margens em polegadas (down, left, top, right)

matplot(seq(1, 30000, by=1500), D[seq(1, 100, by=5),], xlab="Number of games", ylab="Expected return", cex.lab=1.5, type="b", col=c(4, 3, 2, 4, 1), pch=1:5, lwd=1, cex=0.65)

legend("bottomright", c("EMSF+PISF (m=10)", "EMSF+PISF (m=50)", "EMSF+PISF (m=100)", "CNT+PI", "Dealer's strategy"), col=c(4, 3, 2, 4, 1), pch=1:5, lwd=1)

grid()
