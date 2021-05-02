#!/bin/bash

INPUTFILE="INPUTFILE ../../../SimulationResults/Builder/Rising/36kPaper/36K"
FILEEXT1="keV"
FILEEXT2="ps.root"
OUTPUTFILE="OUTPUTFILE ../../../SimulationResults/Reconstructor/Rising/36kPaper/36K"
SPECTRABINANDRANGE="SPECTRABINANDRANGE 250 0. 2000."
BETADOPPLERAVERAGE="BETADOPPLERAVERAGE "
BETATOFAVERAGE="BETATOFAVERAGE  0.5631"
BETAAFTERTARGETAVERAGE="BETAAFTERTARGETAVERAGE 0.5297"
DECAYPOSITIONZ="DECAYPOSITIONZ "
CLUSTERINCLUDE="CLUSTERINCLUDE 1"
MINIBALLINCLUDE="MINIBALLINCLUDE 1"
CLUSTERFIFIND="CLUSTERFIFIND "
MINIBALLFIFIND="MINIBALLFIFIND "
HECTORINCLUDE="HECTORINCLUDE 1"
STATISTICSREDUCTIONFACTOR="STATISTICSREDUCTIONFACTOR 1"
END="END"
FILE="./input/RisingReconstructor.in"

Beta=("0.5463" "0.544" "0.5421" "0.5405" "0.5393" "0.5383" "0.5374" "0.5367" "0.5361" "0.5356" "0.5352" "0.5348" "0.5345" "0.5342" "0.5340" "0.5337" "0.5335" "0.5333" "0.5332" "0.5330" "0.5329");
ZPosition=("0.0" "0.029"  "0.057"  "0.085"  "0.113"  "0.1407" "0.168"  "0.196"  "0.2239" "0.2520" "0.2791" "0.3073" "0.3351" "0.3621" "0.3894" "0.4179" "0.4456" "0.4734" "0.5006" "0.5285" "0.5560");

a=0
LIMIT1=822
LIMIT2=30

for ((arg2=0; arg2 <= LIMIT2 ; arg2=arg2+2))
do
for ((arg1=802; arg1 <= LIMIT1 ; arg1=arg1+2))
do
    BETA=${Beta[$a]}
    ZPOSITION=${ZPosition[$a]}
    for arg3 in 1  2  
    do   
        echo ""$INPUTFILE$arg1$FILEEXT1$arg2$FILEEXT2"
"$OUTPUTFILE$arg1$FILEEXT1$arg2$FILEEXT2"
"$CLUSTERINCLUDE"
"$MINIBALLINCLUDE"
"$HECTORINCLUDE"
"$MINIBALLFIFIND$arg3"
"$CLUSTERFIFIND$arg3"
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
    done
    let a=$a+1    
done
exit 0
