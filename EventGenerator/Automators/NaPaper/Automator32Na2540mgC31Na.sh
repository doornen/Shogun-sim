#!/bin/bash

BEAMISOTOPE="BEAMISOTOPE 32 11 11"
BEAMENERGY="BEAMENERGY 263.7 1.83"  
BEAMPOSITION="BEAMPOSITION -0.07 .53 -0.19 .68"
BEAMANGLE="BEAMANGLE 0.0 0.0 0.0 360.0"
TARGET="TARGET 3 3.0 3.0 2540"  
TARGETANGULARBROADENING="TARGETANGULARBROADENING 0 0.00"
MASSCHANGE="MASSCHANGE 1 0"
BORREL="BORREL 1 8."
GOLDHABER="GOLDHABER 1 90."
GAMMAINPUT="GAMMAINPUT ./input/31Na.in"
THETARANGE="THETARANGE 0.0 180.0"
NUMBEROFEVENTS="NUMBEROFEVENTS 100000"
DEFAULTCUTVALUE="DEFAULTCUTVALUE 0.001"
OUTPUTFILE1="OUTPUTFILE ../SimulationResults/Generator/31na/32na31na262mev2540mgC365kev"
OUTPUTFILE2="ps.root"
DEDXTABLE="DEDXTABLE 1 ./dEdXTables/NaOnC.in ./dEdXTables/NaOnC.in"
END="END"


FILE="./input/EventGenerator.in"
FILE2="./input/31Na.in"

for arg1 in 0 5 10 15 20 25 30 35 40 45 50
do
echo "LEVEL 0   0.    0.  0.
LEVEL 1 100. 365.  $arg1.
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