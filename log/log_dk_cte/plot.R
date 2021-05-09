D <- NULL

tmp <- read.table("/home/rafaelbeirigo/emsf/dp_example/log/e_emsf_100_1_5000_10_1e+20_30_20.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

tmp <- read.table("/home/rafaelbeirigo/emsf/dp_example/log/e_emsf_100_1_5000_10_1e+20_30_30.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

tmp <- read.table("/home/rafaelbeirigo/emsf/dp_example/log/e_emsf_100_1_5000_10_1e+20_30_40.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)


tmp <- read.table("/home/rafaelbeirigo/emsf/dp_example/log/e_emsf_100_1_5000_10_1e+20_30_50.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

matplot(D, t="o" )
