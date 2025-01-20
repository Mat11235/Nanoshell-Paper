#!/bin/bash
export LC_NUMERIC="en_US.UTF-8"
alifol="../img/sgl_alp_emi/"
echo
echo "> Compiling codes..."
g++ -Wall -I/usr/local/include -L/usr/local/lib ../src/sfrohlich.cxx -o ../bin/sfr -lgsl -lgslcblas -lm -larmadillo
g++ -Wall -I/usr/include/ -L/usr/local/lib ../src/single_ome_al_p3.cxx -o ../bin/oap -lgsl -lgslcblas -lm -larmadillo
g++ -Wall -I/usr/include/ -L/usr/local/lib ../src/single_omeG_p3.cxx -o ../bin/oGp -lgsl -lgslcblas -lm -larmadillo
echo "> ...Done!"
echo
none=1.e-30
echo "> Setting the gain at the frohlich frequency"
echo "> and setting the probe field to "$none"..."
read  a  Dome  ome21  G  omemi  omema  metal  model  gain_model  solvent E0 rap host< ../data/input/nanosphere_eV.dat
fro=(`./../bin/sfr`)
omeG=${fro[0]}
Gth=${fro[1]}
echo $a $Dome $omeG $G $omemi $omema $metal $model $gain_model $solvent $none $rap $host> ../data/input/nanosphere_eV.dat
focurr=$metal"-in-"$solvent"/"
echo "> ...Done!"
echo
echo "> the emission threshold is Gth = "$Gth
echo
echo "> Producing BOUNDARY MAPS"
echo $a $Dome $omeG $G 2. 4.5 $metal $model $gain_model $solvent $E0 $rap $host> ../data/input/nanosphere_eV.dat
./../bin/oGp
ogpfol="../img/oGp/"
mnewdir=$ogpfol$focurr
echo "> Removing old files from the image folder: "$ogpfol
rm -fr $mnewdir
echo "> ...Done!"
echo "> Creating image folder: "$mnewdir
mkdir $mnewdir
echo $a $Dome $omeG $G $omemi $omema $metal $model $gain_model $solvent $none $rap $host> ../data/input/nanosphere_eV.dat
gnuplot emi_bndr.gp
mv $ogpfol*.png $mnewdir

echo
echo
echo " Producing SPECTRA"
echo
anewdir=$alifol$focurr
echo "> Removing old files from the image folder: "$alifol
rm -fr $anewdir
echo "> ...Done!"
echo "> Creating image folder: "$anewdir
mkdir $anewdir
echo "> Creating image folders: map_bdr"
mkdir $mnewdir"map_bdr"
echo "> ...Done!"
echo 
focurr=$metal"-in-"$solvent"/"
echo "> Removing old files from the image folder: GIF"
gifdir="../img/GIF/"$focurr
rm -fr $gifdir
echo "> Creating image folders: GIF/"$focurr
echo $gifdir
mkdir $gifdir
mkdir $gifdir"alpha/"
for ii in {0..20..1}; do
    iG=`echo $ii*0.1|bc -l`
    iG=`printf '%.3f\n' $iG`
    iname=`printf "%03d" $ii`
    mG=`echo $iG*$Gth|bc -l`
    echo -e "   > Calculating the emission spectrum for G = "$iG"*Gth\t=\t"$mG
    echo $a $Dome $omeG $mG 2. 4.5 $metal $model $gain_model $solvent $E0 $rap $host> ../data/input/nanosphere_eV.dat
    ./../bin/sap
    echo $a $Dome $omeG $G $omemi $omema $metal $model $gain_model $solvent $none $rap $host> ../data/input/nanosphere_eV.dat
    gnuplot emi_al.gp
    mv "../img/emi_al.pdf" $anewdir"/emi_"$iG"Gth.pdf"
    convert -density 300 $anewdir"/emi_"$iG"Gth.pdf" $anewdir"/emi_"$iname".png"
    rm $anewdir"/emi_"$iG"Gth.pdf"
    echo $a $Dome $omeG $mG $omemi $omema $metal $model $gain_model $solvent $E0 $rap $host> ../data/input/nanosphere_eV.dat
    gnuplot -c emi_bndr_step.gp $mG
    mv $ogpfol"iso_alka_step.png" $mnewdir"map_bdr/emi_"$iname".png"
    convert -resize x940 $anewdir"/emi_"$iname".png" "../img/temp0.png"
    convert -resize x940 $mnewdir"map_bdr/emi_"$iname".png" "../img/temp1.png"
    convert +append "../img/temp0.png" "../img/temp1.png" $gifdir"alpha/emi_"$iname".png"
    rm "../img/temp0.png" 
    rm "../img/temp1.png"
done
convert -delay 20 -loop 0 -dispose previous $gifdir"alpha/*.png"  $gifdir"alpha.gif"
echo "> Resetting the original input file..."
echo $a $Dome $ome21 $G $omemi $omema $metal $model $gain_model $solvent $E0 $rap $host> ../data/input/nanosphere_eV.dat
echo "> ...Done!"
