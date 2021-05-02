#!/bin/bash
cd input
ln -sf Reconstructor.in.dummy ShogunReconstructor.in
cd ..
INPUTFILE="INPUTFILE ../../../SimulationResults/Builder/Shogun/Efficiency_"
FILEEXT1="mev"
FILEEXT2="kev_gagg_conf1_no_res.root"
OUTPUTFILE="OUTPUTFILE ../../../SimulationResults/Reconstructor/Shogun/Efficiency_"
SPECTRABINANDRANGE="SPECTRABINANDRANGE 2000 0. 20000."
BETADOPPLERAVERAGE="BETADOPPLERAVERAGE "
BETATOFAVERAGE="BETATOFAVERAGE "
DECAYPOSITIONZ="DECAYPOSITIONZ 0.0"
FIFIND="FIFIND 2"
ADDBACK="ADDBACK 1 15"
END="END"

FILE="./input/ShogunReconstructor.in"

beta="0"
BETA1="0.4295"
BETA2="0.5081" 
BETA3="0.5677"

for arg1 in 100 150 200 
do
    a=0
    for arg2 in 250 500 1000 2000 5000 10000  
    
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
        
        echo "
"$INPUTFILE$arg1$FILEEXT1$arg2$FILEEXT2"
"$OUTPUTFILE$arg1$FILEEXT1$arg2$FILEEXT2"
"$SPECTRABINANDRANGE"
"$BETADOPPLERAVERAGE$beta"
"$BETATOFAVERAGE$beta"
"$DECAYPOSITIONZ"
"$FIFIND"
"$ADDBACK"
"$END" 
" > $FILE
        
        time ./run.sh
    done
done

exit 0
