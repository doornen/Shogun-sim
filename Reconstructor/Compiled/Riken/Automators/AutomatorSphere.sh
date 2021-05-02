#!/bin/bash

INPUTFILE="INPUTFILE ../../../SimulationResults/Builder/Shogun/Efficiency/Case"
FILEEXT0="_"
FILEEXT1="MeV"
FILEEXT2="keV.root"
OUTPUTFILE="OUTPUTFILE ../../../SimulationResults/Reconstructor/Shogun/Efficiency/Case"
SPECTRABINANDRANGE="SPECTRABINANDRANGE 1000 0. 5000."
BETADOPPLERAVERAGE="BETADOPPLERAVERAGE "
BETATOFAVERAGE="BETATOFAVERAGE "
DECAYPOSITIONZ="DECAYPOSITIONZ 0.0"
END="END"


FILE="./input/ShogunReconstructor.in"

beta="0"
BETA1="0.4295"
BETA2="0.5081"
BETA3="0.5677"
#BETA3="0.6152"
#BETA3="0.6541"


for arg1 in 100 150 200  
do
    a=0
    for arg2 in 250 500 750 1000 1250 1500 1750 2000 2250 2500 2750 3000 
    
    do
        echo $beta
        echo $arg1
        echo $arg2
        
        if [ $arg1 -eq 100 ] 
        then
            beta=$BETA1
            echo $beta
        fi
        if [ $arg1 -eq 150 ]
        then 
            beta=$BETA2
            echo $beta
        fi
        if [ $arg1 -eq 200 ]
        then 
            beta=$BETA3
            echo $beta
        fi
        if [ $arg1 -eq 250 ]
        then 
            beta=$BETA4
            echo $beta
        fi
        if [ $arg1 -eq 300 ]
        then 
            beta=$BETA4
            echo $beta
        fi
        
        echo ""$INPUTFILE$1$FILEEXT0$arg1$FILEEXT1$arg2$FILEEXT2"
"$OUTPUTFILE$1$FILEEXT0$arg1$FILEEXT1$arg2$FILEEXT2"
"$SPECTRABINANDRANGE"
"$BETADOPPLERAVERAGE$beta"
"$BETATOFAVERAGE$beta"
"$DECAYPOSITIONZ"
"$END" 
" > $FILE
        
        time ./ShogunReconstructor
    done
done

exit 0




