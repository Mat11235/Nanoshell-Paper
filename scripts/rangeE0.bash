#!/bin/bash
echo "> Compiling codes..."
g++ -Wall -I/usr/local/include -L/usr/local/lib ../src/frohlich.cxx -o ../bin/fro -lgsl -lgslcblas -lm -larmadillo
g++ -Wall -I/usr/include/ -L/usr/local/lib ../src/nanoshell_num_anl.cxx -o ../bin/nsh -lgsl -lgslcblas -lm -larmadillo
echo "> done!"
# File path
file_path="../data/input/nanosphere_eV.dat"

# Read the variables from the file
read r Dome ome0 G omemi omema mtl mdl atv hst E0 rho sol< <(grep -v '^#' "$file_path")

# Backup the original file
cp "$file_path" "${file_path}.bak"

fro_output=($(./../bin/fro))
ome_sp=${fro_output[0]}
new_G=`echo "1.01*${fro_output[1]}"| bc -l`

echo "> Calculating for different E0..."
for i in 8 9 10 11; do 
    newE0="1.e-"$i
    echo $r $Dome $ome_sp $new_G $omemi $omema $mtl $mdl $atv $hst $newE0 $rho $sol> $file_path
    ./../bin/nsh 2.811 # $ome_sp
    cp "../data/output/numtime.dat" "../data/output/numtime$(printf "%02d" $i).dat"
done
echo "> done!"

# Reinstate the original content of the file
echo "> Reinstating the original content of $file_path"
{
  echo "$r $Dome $ome0 $G $omemi $omema $mtl $mdl $atv $hst $E0 $rho $sol"
  echo "# r Dome ome0 G omemi omema mtl mdl atv hst E0 rho sol"
} > "$file_path"
echo "> done!"
