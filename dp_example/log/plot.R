D <- NULL

## tmp <- read.table("/home/rafaelbeirigo/emsf/dp_example/log/e_cnt.log")
## tmp <- apply(tmp, 2, mean)
## D <- cbind(D, tmp)

## tmp <- read.table("/home/rafaelbeirigo/emsf/dp_example/log/e_emsf_a.log")
## tmp <- apply(tmp, 2, mean)
## D <- cbind(D, tmp)

## tmp <- read.table("/home/rafaelbeirigo/emsf/dp_example/log/e_emsf_b.log")
## tmp <- apply(tmp, 2, mean)
## D <- cbind(D, tmp)

## tmp <- read.table("/home/rafaelbeirigo/emsf/dp_example/log/e_emsf_c.log")
## tmp <- apply(tmp, 2, mean)
## D <- cbind(D, tmp)

tmp <- read.table("/home/rafaelbeirigo/log/e_emsf_h.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

tmp <- read.table("/home/rafaelbeirigo/log/e_emsf_i.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

tmp <- read.table("/home/rafaelbeirigo/log/e_emsf_j.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

tmp <- read.table("/home/rafaelbeirigo/log/e_emsf_k.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

matplot(D, t="o" )
