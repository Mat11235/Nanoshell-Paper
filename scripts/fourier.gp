reset
set arrow nohead from first -0.0011666698, graph 0 to  first -0.0011666698, graph 1

set ylabel "Magnitude"
set xlabel "ℏ{/Symbol W} (eV)"
plot[-0.004:0.002] "../data/output/fourier.dat" u ($1):(abs($2)) w p pt 7 ps 3 lc rgb "red" t "", "" u ($1):(abs($2)) w p pt 6 ps 3 lc rgb "black" t ""
