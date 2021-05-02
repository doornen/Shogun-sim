#!/bin/bash

INPUTFILE="INPUTFILE ../../../SimulationResults/Builder/Shogun/32mg/Case"
FILEEXT0="mg"
FILEEXT1="ps.root"
FILEEXT2="_32mg32mg200mev"
OUTPUTFILE="OUTPUTFILE ../../../SimulationResults/Reconstructor/Shogun/32mg/Case"
SPECTRABINANDRANGE="SPECTRABINANDRANGE 250 0. 2000."
BETADOPPLERAVERAGE="BETADOPPLERAVERAGE "
BETATOFAVERAGE="BETATOFAVERAGE 0.5677 "
DECAYPOSITIONZ="DECAYPOSITIONZ 0.0"
END="END"

FILE="./input/Reconstructor.in"

beta="0"
BETA=(0.5676 0.5628 0.5576 0.5522 0.5465 0.5404)
a=0
for arg1 in 1 500 1000 1500 2000 2500
do
    
    for arg2 in 0 10 20 30 40 50 60 70 80 90 100
    
    do
        echo $beta
        echo $arg1
        echo $arg2
        beta=${BETA[$a]}
        echo $beta
        echo ""$INPUTFILE$1$FILEEXT2$arg1$FILEEXT0$arg2$FILEEXT1"
"$OUTPUTFILE$1$FILEEXT2$arg1$FILEEXT0$arg2$FILEEXT1"
"$SPECTRABINANDRANGE"
"$BETADOPPLERAVERAGE$beta"
"$BETATOFAVERAGE"
"$DECAYPOSITIONZ"
"$END" 
" > $FILE
        
        time ./ShogunReconstructor
    done
    let a=a+1
done

exit 0




