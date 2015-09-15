set term png enhanced size 1366,768 font ",20"
set encoding utf8

set output "grafico.png"

titulo = ""

set rmargin 4
set title titulo
set xlabel ""
set ylabel ""

plot "v_cnt_bj_20_2_30000_1000000_0_100_15.log" w l lw 2, \
"v_cnt_bj_20_2_30000_1000000_0_100_14.log" w l lw 2, \
"v_cnt_bj_20_2_30000_1000000_0_100_13.log" w l lw 2, \
"v_cnt_bj_20_2_30000_1000000_0_100_12.log" w l lw 2, \
"v_cnt_bj_20_2_30000_1000000_0_100_11.log" w l lw 2, \
"v_cnt_bj_20_2_30000_1000000_0_100_10.log" w l lw 2, \
"v_cnt_bj_20_2_30000_1000000_0_100_09.log" w l lw 2, \
"v_cnt_bj_20_2_30000_1000000_0_100_08.log" w l lw 2, \
"v_cnt_bj_20_2_30000_1000000_0_100_07.log" w l lw 2, \
"v_cnt_bj_20_2_30000_1000000_0_100_06.log" w l lw 2, \
"v_cnt_bj_20_2_30000_1000000_0_100_05.log" w l lw 2, \
"v_cnt_bj_20_2_30000_1000000_0_100_04.log" w l lw 2, \
"v_cnt_bj_20_2_30000_1000000_0_100_03.log" w l lw 2, \
"v_cnt_bj_20_2_30000_1000000_0_100_02.log" w l lw 2, \
"v_cnt_bj_20_2_30000_1000000_0_100_01.log" w l lw 2
