#!/bin/bash

INPUTFILE="INPUTFILE ../../../SimulationResults/Builder/RisingFast/efficiency/"
FILEEXT1="mev_500kev_"
FILEEXT2="ps.root"
OUTPUTFILE="OUTPUTFILE ../../../SimulationResults/Reconstructor/RisingFast/efficiency/"
SPECTRABINANDRANGE="SPECTRABINANDRANGE 500 0. 1000."
BETADOPPLERAVERAGE="BETADOPPLERAVERAGE "
BETATOFAVERAGE="BETATOFAVERAGE "
DECAYPOSITIONZ="DECAYPOSITIONZ "
CLUSTERINCLUDE="CLUSTERINCLUDE 1"
MINIBALLINCLUDE="MINIBALLINCLUDE 1"
HECTORINCLUDE="HECTORINCLUDE 1"
END="END"
FILE="./input/RisingReconstructor.in"

beta="0"
BETA1="0.4295"
BETA2="0.5081"
BETA3="0.5677"


POSITIONZ1=("0.0" "0.2065" "0.413" "0.6195" "0.826" "1.0325" "1.239" "1.4455" "1.652" "1.8585" "2.065");
POSITIONZ2=("0.0" "0.25854" "0.5171" "0.77562" "1.03416" "1.2927" "1.55124" "1.90978" "2.06832" "2.3269" "2.5854");
POSITIONZ3=( "0.0" "0.301" "0.602" "0.903" "1.204" "1.505" "1.806" "2.107" "2.408" "2.709" "3.012");



for arg1 in 100 150 200 
do
    a=0
    for arg2 in 000 010 020 030 040 050 060 070 080 090 100 
    
    do
        echo $beta
        echo $arg1
        echo $arg2
        
        if [ $arg1 -eq 100 ] 
        then
            beta=$BETA1
            echo $beta
            positionz=${POSITIONZ1[$a]}
            echo $positionz
        fi
        if [ $arg1 -eq 150 ]
        then 
            beta=$BETA2
            echo $beta
            positionz=${POSITIONZ2[$a]}
            echo $positionz
        fi
        if [ $arg1 -eq 200 ]
        then 
            beta=$BETA3
            echo $beta
            positionz=${POSITIONZ3[$a]}
            echo $positionz
        fi
        let a=$a+1    
        
        echo ""$INPUTFILE$arg1$FILEEXT1$arg2$FILEEXT2"
"$OUTPUTFILE$arg1$FILEEXT1$arg2$FILEEXT2"
"$CLUSTERINCLUDE"
"$MINIBALLINCLUDE"
"$HECTORINCLUDE"
"$SPECTRABINANDRANGE"
"$BETADOPPLERAVERAGE$beta"
"$BETATOFAVERAGE$beta"
"$DECAYPOSITIONZ$positionz"
"$END" 
" > $FILE
        
        time ./RisingReconstructor
    done
done

exit 0
