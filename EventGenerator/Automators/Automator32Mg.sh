#!/bin/bash

BEAMISOTOPE="BEAMISOTOPE 32 12 12"
BEAMENERGY="BEAMENERGY 200.0 2.000"
BEAMPOSITION="BEAMPOSITION 0.0 1.00 0.0 1.00"
BEAMANGLE="BEAMANGLE 0.0 0.658 0.0 360.0"
TARGET="TARGET 3 7. 7. "
TARGETANGULARBROADENING="TARGETANGULARBROADENING 0 1.563"
MASSCHANGE="MASSCHANGE 0 0"
BORREL="BORREL 0 8."
GOLDHABER="GOLDHABER 0 90."
GAMMAINPUT="GAMMAINPUT ./input/32Mg2+.in"
THETARANGE="THETARANGE 0.0 180.0"
NUMBEROFEVENTS="NUMBEROFEVENTS 100000"
DEFAULTCUTVALUE="DEFAULTCUTVALUE 0.01"
OUTPUTFILE1="OUTPUTFILE ../SimulationResults/Generator/32mg/32mg32mg200mev"
OUTPUTFILE2="mg"
OUTPUTFILE3="ps.root"
DEDXTABLE="DEDXTABLE 1 ./dEdXTables/MgOnCLowEnergies.in ./dEdXTables/MgOnCLowEnergies.in"
END="END"


FILE="./input/EventGenerator.in"
FILE2="./input/32Mg2+.in"

rm ./input/EventGenerator.in

for arg1 in 0 10 20 30 40 50 60 70 80 90 100
do
echo "LEVEL 0   0.    0.  0.
LEVEL 1 100. 1000 $arg1.
DECAY 1 0 100.
END
" > $FILE2

for arg2 in 1 500 1000 1500 2000 2500 
do
echo ""$BEAMISOTOPE"
"$BEAMENERGY"
"$BEAMPOSITION"
"$BEAMANGLE"
"$TARGET$arg2"
"$TARGETANGULARBROADENING"
"$MASSCHANGE"
"$BORREL"
"$GOLDHABER"
"$GAMMAINPUT"
"$THETARANGE"
"$NUMBEROFEVENTS"
"$DEFAULTCUTVALUE"
"$OUTPUTFILE1$arg2$OUTPUTFILE2$arg1$OUTPUTFILE3"
"$DEDXTABLE"
"$END" 
" > $FILE

time EventGenerator run_nothing.mac
done
done
exit 0
