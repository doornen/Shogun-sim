#!/bin/bash

BEAMISOTOPE="BEAMISOTOPE 57 29 29"
BEAMENERGY="BEAMENERGY 150 5"
BEAMPOSITION="BEAMPOSITION 0.0 0.5 0.0 0.5"
BEAMANGLE="BEAMANGLE 0.0 0.5 0.0 360.0"
TARGET="TARGET 4 7.0 7.0 500"
TARGETANGULARBROADENING="TARGETANGULARBROADENING 1 1."
MASSCHANGE="MASSCHANGE 1 1"
BORREL="BORREL 1 8."
GOLDHABER="GOLDHABER 1 90."
GAMMAINPUT="GAMMAINPUT ./input/74ni.in"
THETARANGE="THETARANGE 0.0 180.0"
NUMBEROFEVENTS="NUMBEROFEVENTS 200000"
DEFAULTCUTVALUE="DEFAULTCUTVALUE 0.001"
OUTPUTFILE1="OUTPUTFILE ../SimulationResults/Generator/74ni/75cu_74ni_150mev_500mgfe_1024kev"
OUTPUTFILE2="ps.root"
DEDXTABLE="DEDXTABLE 1 ./dEdXTables/CuOnFedEdX.in  ./dEdXTables/NiOnFedEdX.in "
END="END"


FILE="./input/EventGenerator.in"
FILE2="./input/74ni.in"


halflife1="0.5"
halflife2="0.6"
halflife3="0.7"
halflife4="0.8"
halflife5="0.9"
halflife6="1.0"
halflife7="1.1"
halflife8="1.2"
halflife9="1.3"
halflife10="1.4"
halflife11="1.5"



for arg1 in "$halflife1" "$halflife2" "$halflife3" "$halflife4" "$halflife5" "$halflife6" "$halflife7" "$halflife8" "$halflife9" "$halflife10" "$halflife11" 
do

echo "LEVEL 0   0.    0.  0.
LEVEL 1 100. 1024. $arg1
DECAY 1 0 100.
END
" > $FILE2

echo ""$BEAMISOTOPE"
"$BEAMENERGY"
"$BEAMPOSITION"
"$BEAMANGLE"
"$TARGET"
"$TARGETANGULARBROADENING"
"$MASSCHANGE"
"$BORREL"
"$GOLDHABER"
"$GAMMAINPUT"
"$THETARANGE"
"$NUMBEROFEVENTS"
"$DEFAULTCUTVALUE"
"$OUTPUTFILE1$arg1$OUTPUTFILE2"
"$DEDXTABLE"
"$END" 
" > $FILE

time EventGenerator run_nothing.mac
done

exit 0
