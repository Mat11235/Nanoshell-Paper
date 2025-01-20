#!/bin/bash
reset 
# DEFINING COMMAND TO READ FROM files
at(file, row, col) = system( sprintf("awk -v row=%d -v col=%d 'NR == row {print $col}' %s", row, col, file) )
# DONE

# READING SPECTRUM RANGE (omin omax) FROM "../data/input/nanosphere_eV.dat"
file="../data/input/nanosphere_eV.dat" ;
row=1 ; col=5
omin = at(file,row,col)
row=1 ; col=6
omax = at(file,row,col)
# DONE

# READING GAIN QUANTITY (G) FROM "../data/input/nanosphere_eV.dat"
row=1 ; col=4
G = abs(at(file,row,col))
# DONE

# READING GAIN THRESHOLD (Gth) FROM "../data/output/frohlich.dat"
file="../data/output/frohlich.dat" ;
row=1 ; col=2
Gth = abs(at(file,row,col))
# DONE

# READING FREQUENCY (ome) FROM "../data/output/omega.dat"
file="../data/output/omega.dat" ;
row=1 ; col=1
omega = abs(at(file,row,col))
# DONE

# READING TOTAL TIME AND TIME_PUMP_ON (ome T tpump) FROM "../data/input/time.dat"
file="../data/input/time.dat" ;
row=1 ; col=1
T = abs(at(file,row,col))
row=1 ; col=2
tpump = at(file,row,col)
# DONE

# READING STATIONARY QUASI STATIC ALPHA (Re(alph_QS), Im(alph_QS)) FROM "../data/output/omega.dat"
file="../data/output/omega.dat" ;
row=1 ; col=2
reaQS = at(file,row,col)
row=1 ; col=3
imaQS = at(file,row,col)
# DONE

set term unknown
if (G<=Gth){
    plot "../data/output/stationary.dat" u 1:2, "" u 1:3
    if (GPVAL_Y_MAX>300) {
    set yrange [-300:300];
    } else {
            set yrange [GPVAL_Y_MIN:GPVAL_Y_MAX];
            }
    } else {
    plot "../data/output/numtime.dat" u 1:2, "" u 1:3
    set yrange [GPVAL_Y_MIN:GPVAL_Y_MAX];
    }

# set xrange [0:20];
set xrange [0:T];
set term pdf color enhanced size 20cm, 8cm font ",16";
set output '../img/intime.pdf';
set style line 6 lt 1 lc 0 lw 4;
set style line 7 lt 1 lc 8 lw 4;
set style line 8 dt (18,6) lc rgb "purple" lw 3;
set style line 9 dt (18,6) lc rgb "dark-green" lw 3;
set size 2,1.;
set title sprintf("G = %.2f",G);
set multiplot
set size 0.56,1.;

set ylabel "Normalized Polarizability {/Symbol a}/(4{/Symbol p}a^3)"

# set grid ytics

if (G!=0) { 
    set label 33 at first 1.1*tpump, graph 0.95 "{/=14 PUMP ON}";
    set object 66 rect from first tpump, graph 0 to first T, graph 1  fc rgb "light-blue" fillstyle transparent solid;
    }

set xlabel "Time (ps)"
set key bottom
plot "../data/output/numtime.dat" u 1:2 w l ls 6 t "numerical {/Symbol a}'",\
     "../data/output/numtime.dat" u 1:3 w l ls 7 t "numerical {/Symbol a}''",\
     "../data/output/anltime.dat" u 1:2 w l ls 8 t "analytical {/Symbol a}'",\
     "../data/output/anltime.dat" u 1:3 w l ls 9 t "analytical {/Symbol a}''";

set xrange [omin:omax];
unset label 33
set title sprintf("ℏ{/Symbol w} = %.2f",omega) tc rgb "blue";
unset ylabel

unset ytics
set y2label "Normalized Polarizability {/Symbol a}/(4{/Symbol p}a^3)"
set xlabel "Energy (eV)"

set origin 0.55,0
set rmargin 17

set object 66 rect from graph 0, graph 0 to graph 1, graph 1
set arrow 22 nohead from first omega, graph 0 to first omega, graph 1 lc rgb "blue" lw 2 front
set arrow 23 nohead from graph 0, first reaQS to first omega, first reaQS lc 0 lw 2 dt 2
set arrow 24 nohead from graph 0, first imaQS to first omega, first imaQS ls 7 lw 2 dt 2


plot[omin:omax] "../data/output/stationary.dat"  w l ls 6 t "{/Symbol a}'", \
                "" u 1:3 w l ls 7 t "{/Symbol a}''";

unset multiplot
unset term
! convert -density 300 ../img/intime.pdf ../img/intime.png
