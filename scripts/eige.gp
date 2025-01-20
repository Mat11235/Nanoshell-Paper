reset
# set multiplot
# set size 1, 0.5
# plot "../data/output/eigenvalues.dat" w l t "real(k_1)", "" u 1:4  w l t "real(k_2)", "" u 1:6  w l t "real(k_3)"
# set origin 0, 0.5
# plot "../data/output/eigenvalues.dat" u 1:3 w l t "imag(k_1)", "" u 1:5  w l t "imag(k_2)", "" u 1:7  w l t "imag(k_3)"
# unset multiplot

set xzeroaxis
set multiplot
set key left
set xrange [2.2:3.2]
set size 1, 0.5
set lmargin 10
plot [][-0.001:0.0002] "../data/output/eigenvalues.dat"  u 1:6  w l t "real(k_3)"
set origin 0, 0.5
plot "../data/output/eigenvalues.dat"  u 1:7  w l t "imag(k_3)"
unset multiplot
