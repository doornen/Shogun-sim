#!/bin/bash

BEAMISOTOPE="BEAMISOTOPE 112 50 50"
BEAMENERGY1="BEAMENERGY "  
BEAMENERGY2=" 0.001"
BEAMPOSITION="BEAMPOSITION 0.0 0.01 0.0 0.01"
BEAMANGLE="BEAMANGLE 0.0 0.01 0.0 360.0"
TARGET="TARGET 1 7.0 7.0 500"  
TARGETANGULARBROADENING="TARGETANGULARBROADENING 0 0.00"
MASSCHANGE="MASSCHANGE 0 0"
BORREL="BORREL 0 8."
GOLDHABER="GOLDHABER 0 90."
GAMMAINPUT="GAMMAINPUT ./input/112Sn.in"
THETARANGE="THETARANGE 0.0 180.0"
NUMBEROFEVENTS="NUMBEROFEVENTS 1000000"
DEFAULTCUTVALUE="DEFAULTCUTVALUE 0.0001"
OUTPUTFILE1="OUTPUTFILE ../SimulationResults/Generator/LifetimePaper/112Sn/"
OUTPUTFILE2="mev"
OUTPUTFILE3="fs.root"
DEDXTABLE="DEDXTABLE 1 ./dEdXTables/SnOnAu.in ./dEdXTables/SnOnAu.in"
END="END"


FILE="./input/EventGenerator.in"
FILE2="./input/112Sn.in"

for arg1 in 150
do
#for arg2 in 185 370 740    
for arg2 in 0 250 500 750 999
do
for arg3 in 0 
do
echo "LEVEL 0   0.    0.  0.
LEVEL 1 100. 1256.85  $arg3.$arg2  
DECAY 1 0 100.
END
" > $FILE2
arg4=0
let arg4=$arg2+$arg3*1000

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
"$OUTPUTFILE1$arg1$OUTPUTFILE2$arg4$OUTPUTFILE3"
"$DEDXTABLE"
"$END" 
" > $FILE

time EventGenerator run_nothing.mac
done
done
done
exit 0
