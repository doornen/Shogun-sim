#!/bin/bash

INPUTFILE="INPUTFILE ../../../SimulationResults/Builder/Rising/LifetimePaper/34Cl/34cl461kev"
FILEEXT1="ps.root"
OUTPUTFILE="OUTPUTFILE ../../../SimulationResults/Reconstructor/Rising/LifetimePaper/34Cl/34cl461kev"
FILEEXT2="psBetaBigger.root"
SPECTRABINANDRANGE="SPECTRABINANDRANGE 250 0. 1000."
BETADOPPLERAVERAGE="BETADOPPLERAVERAGE "
BETATOFAVERAGE="BETATOFAVERAGE  0.5631"
BETAAFTERTARGETAVERAGE="BETAAFTERTARGETAVERAGE 0.53227"
DECAYPOSITIONZ="DECAYPOSITIONZ "
CLUSTERINCLUDE="CLUSTERINCLUDE 1"
MINIBALLINCLUDE="MINIBALLINCLUDE 1"
CLUSTERFIFIND="CLUSTERFIFIND "
MINIBALLFIFIND="MINIBALLFIFIND "
HECTORINCLUDE="HECTORINCLUDE 1"
STATISTICSREDUCTIONFACTOR="STATISTICSREDUCTIONFACTOR 1"
END="END"
FILE="./input/RisingReconstructor.in"

Beta=("0.5455"   "0.5436" "0.5420" "0.5406" "0.5396" "0.5387" "0.5380" "0.5375" "0.5369" "0.5365" "0.5362" "0.5359" "0.5356" "0.5354" "0.5352" "0.5349" "0.5348" "0.5346" "0.5345" "0.5343" "0.5342");
ZPosition=("0.0" "0.029"  "0.057"  "0.085"  "0.113"  "0.1407" "0.168"  "0.196"  "0.2239" "0.2520" "0.2791" "0.3073" "0.3351" "0.3621" "0.3894" "0.4179" "0.4456" "0.4734" "0.5006" "0.5285" "0.5560");

a=0
LIMIT=20
for ((arg1=0; arg1 <= LIMIT ; arg1++))
do
    BETA=${Beta[$a]}
    ZPOSITION=${ZPosition[$a]}
    for arg2 in 1  2  
    do   
        echo ""$INPUTFILE$arg1$FILEEXT1"
"$OUTPUTFILE$arg1$FILEEXT2"
"$CLUSTERINCLUDE"
"$MINIBALLINCLUDE"
"$HECTORINCLUDE"
"$MINIBALLFIFIND$arg2"
"$CLUSTERFIFIND$arg2"
"$SPECTRABINANDRANGE"
"$BETADOPPLERAVERAGE$BETA"
"$BETATOFAVERAGE"
"$BETAAFTERTARGETAVERAGE"
"$DECAYPOSITIONZ$ZPOSITION"
"$STATISTICSREDUCTIONFACTOR"
"$END" 
" > $FILE
        
        time ./RisingReconstructor
    done
    let a=$a+1    
done

exit 0
