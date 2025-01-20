set xlabel "ℏ{/Symbol w}"
plot[0:0.04] "../data/output/fourier.dat" w lp pt 6 t "Real part", "" u 1:3 w p pt 7 t "Imaginary part"
