#!/bin/bash
ln -sf Dali2Reconstructor.in.dummy ./input/Dali2Reconstructor.in

INPUTFILE="INPUTFILE ../../../SimulationResults/Builder/Dali2/Efficiency/"
FILEEXT1="MeV"
FILEEXT2="keV_beampipe_1mm_pb_6p_res.root"
OUTPUTFILE="OUTPUTFILE ../../../SimulationResults/Reconstructor/Dali2/Efficiency/"
SPECTRABINANDRANGE="SPECTRABINANDRANGE 1000 0. 5000."
BETADOPPLERAVERAGE="BETADOPPLERAVERAGE "
BETATOFAVERAGE="BETATOFAVERAGE "
DECAYPOSITIONZ="DECAYPOSITIONZ 0.0"
ADDBACK="ADDBACK 1 15"
TRIGGER="TRIGGER 1"
FIFIND="FIFIND "
THRESHOLD="ENERGYTHRESHOLD 50"
END="END"

FILE="./input/Dali2Reconstructor.in"

beta="0"
BETA1="0.4295"
BETA2="0.5081"
BETA3="0.5677"
BETA4="0.6152"
BETA5="0.6541"
find="0"


for arg1 in 100 150 200 250  
do
    a=0
    for arg2 in 250 500 750 1000 1250 1500 1750 2000 2250 2500 2750 3000     
    do
        for arg3 in 1 2 
        do
            echo $beta
            echo $arg1
            echo $arg2
            echo $arg3
            echo $find
            
            find=$arg3
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
#        if [ $arg1 -eq 300 ]
#        then 
#            beta=$BETA5
#            echo $beta
#       fi
            
            echo ""$INPUTFILE$arg1$FILEEXT1$arg2$FILEEXT2"
"$OUTPUTFILE$arg1$FILEEXT1$arg2$FILEEXT2"
"$SPECTRABINANDRANGE"
"$BETADOPPLERAVERAGE$beta"
"$BETATOFAVERAGE$beta"
"$DECAYPOSITIONZ"
"$ADDBACK"
"$TRIGGER"
"$FIFIND$find"
"$THRESHOLD"
"$END" 
" > $FILE
            
            time ./Dali2Reconstructor
        done
    done
done
exit 0




