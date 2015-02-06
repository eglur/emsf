D <- NULL

## tmp <- read.table("/home/rafaelbeirigo/emsf/dp_example/log/e_cnt.log")
## tmp <- apply(tmp, 2, mean)
## D <- cbind(D, tmp)

## tmp <- read.table("/home/rafaelbeirigo/emsf/dp_example/log/e_emsf_a.log")
## tmp <- apply(tmp, 2, mean)
## D <- cbind(D, tmp)

tmp <- read.table("/home/rafaelbeirigo/emsf/dp_example/log/e_emsf_b.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

tmp <- read.table("/home/rafaelbeirigo/emsf/dp_example/log/e_emsf_c.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

tmp <- read.table("/home/rafaelbeirigo/emsf/dp_example/log/e_emsf_d.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

tmp <- read.table("/home/rafaelbeirigo/emsf/dp_example/log/e_emsf_e.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

tmp <- read.table("/home/rafaelbeirigo/emsf/dp_example/log/e_emsf_f.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

tmp <- read.table("/home/rafaelbeirigo/emsf/dp_example/log/e_emsf_g.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

matplot(D, t="o")
