set term png enhanced size 1366,768 font ",20"
set encoding utf8

set output "grafico.png"

titulo = ""

set rmargin 4
set title titulo
set xlabel ""
set ylabel ""

plot "v_cnt_bj_20_2_30000_1000000_0_100_1_01.log.dat" w l lw 2, \
     "v_cnt_bj_20_2_30000_1000000_0_100_0.9_01.log.dat" w l lw 2, \
     "v_cnt_bj_20_2_30000_1000000_0_100_0.8_01.log.dat" w l lw 2, \
     "v_cnt_bj_20_2_30000_1000000_0_100_0.7_01.log.dat" w l lw 2, \
     "v_cnt_bj_20_2_30000_1000000_0_100_0.6_01.log.dat" w l lw 2, \
     "v_cnt_bj_20_2_30000_1000000_0_100_0.5_01.log.dat" w l lw 2, \
     "v_cnt_bj_20_2_30000_1000000_0_100_0.4_01.log.dat" w l lw 2, \
     "v_cnt_bj_20_2_30000_1000000_0_100_0.3_01.log.dat" w l lw 2, \
     "v_cnt_bj_20_2_30000_1000000_0_100_0.2_01.log.dat" w l lw 2, \
     "v_cnt_bj_20_2_30000_1000000_0_100_0.1_01.log.dat" w l lw 2, \
     "v_cnt_bj_20_2_30000_1000000_0_100_0.05_01.log.dat" w l lw 2, \
     "v_cnt_bj_20_2_30000_1000000_0_100_0.01_01.log.dat" w l lw 2
