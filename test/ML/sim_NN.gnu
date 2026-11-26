# Sample written by sim_NN (x1, x2, network, Ackley): gnuplot -p -e 'file="ackley2d_KAN-TITO-h5-g8x8.out"' sim_NN.gnu
set xlabel "x1"; set ylabel "x2"; set zlabel "z"
set ticslevel 0
splot file using 1:2:4 with points pt 7 ps 0.4 lc rgb "gray" title "Ackley", \
      file using 1:2:3 with points pt 7 ps 0.6 palette title file noenhanced
