D <- NULL

## tmp <- read.table("/home/rafaelbeirigo/emsf/dp_example/emsf/dp_example/n_100_T_5000/t_cnt.log")
## tmp <- apply(tmp, 2, mean)
## D <- cbind(D, tmp)

## tmp <- read.table("/home/rafaelbeirigo/emsf/dp_example/emsf/dp_example/n_100_T_5000/t_emsf_a.log")
## tmp <- apply(tmp, 2, mean)
## D <- cbind(D, tmp)

## tmp <- read.table("/home/rafaelbeirigo/emsf/dp_example/emsf/dp_example/n_100_T_5000/t_emsf_b.log")
## tmp <- apply(tmp, 2, mean)
## D <- cbind(D, tmp)

## tmp <- read.table("/home/rafaelbeirigo/emsf/dp_example/emsf/dp_example/n_100_T_5000/t_emsf_c.log")
## tmp <- apply(tmp, 2, mean)
## D <- cbind(D, tmp)

tmp <- read.table("/home/rafaelbeirigo/emsf/dp_example/n_100_T_5000/t_emsf_h.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

tmp <- read.table("/home/rafaelbeirigo/emsf/dp_example/n_100_T_5000/t_emsf_i.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

tmp <- read.table("/home/rafaelbeirigo/emsf/dp_example/n_100_T_5000/t_emsf_j.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

tmp <- read.table("/home/rafaelbeirigo/emsf/dp_example/n_100_T_5000/t_emsf_k.log")
tmp <- apply(tmp, 2, mean)
D <- cbind(D, tmp)

matplot(D, t="o" )
