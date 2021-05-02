#!/bin/bash

BEAMISOTOPE="BEAMISOTOPE 37 20 20"
BEAMENERGY1="BEAMENERGY "  
BEAMENERGY2=" 0.1"
BEAMPOSITION="BEAMPOSITION 0.0 0.01 0.0 0.01"
BEAMANGLE="BEAMANGLE 0.0 0.01 0.0 360.0"
TARGET="TARGET 2 7.0 7.0 0.1"  
TARGETANGULARBROADENING="TARGETANGULARBROADENING 0 0.00"
MASSCHANGE="MASSCHANGE 0 0"
BORREL="BORREL 0 8."
GOLDHABER="GOLDHABER 0 90."
GAMMAINPUT="GAMMAINPUT ./input/efficiency.in"
THETARANGE="THETARANGE 0.0 180.0"
NUMBEROFEVENTS="NUMBEROFEVENTS 200000"
DEFAULTCUTVALUE="DEFAULTCUTVALUE 0.01"
OUTPUTFILE1="OUTPUTFILE ../SimulationResults/Generator/efficiency/"
OUTPUTFILE2="mev_500kev_"
OUTPUTFILE3="ps.root"
DEDXTABLE="DEDXTABLE 0 ./dEdXTables/CaOnBedEdX.in"
END="END"


FILE="./input/EventGenerator.in"
FILE2="./input/efficiency.in"

energy1="100"
energy2="150"
energy3="200"


halflife0="000"
halflife1="010"
halflife2="020"
halflife3="030"
halflife4="040"
halflife5="050"
halflife6="060"
halflife7="070"
halflife8="080"
halflife9="090"
halflife10="100"

for arg1 in "$energy1" "$energy2" "$energy3" 
do
for arg2 in "$halflife0" "$halflife1" "$halflife2" "$halflife3" "$halflife4" "$halflife5" "$halflife6" "$halflife7" "$halflife8" "$halflife9" "$halflife10" 

do
echo "LEVEL 0   0.    0.  0.
LEVEL 1 100. 500. $arg2.
DECAY 1 0 100.
END
" > $FILE2

echo ""$BEAMISOTOPE"
"$BEAMENERGY1$arg1$BEAMENERGY2"
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
"$OUTPUTFILE1$arg1$OUTPUTFILE2$arg2$OUTPUTFILE3"
"$DEDXTABLE"
"$END" 
" > $FILE

time EventGenerator run_nothing.mac
done
done

exit 0
