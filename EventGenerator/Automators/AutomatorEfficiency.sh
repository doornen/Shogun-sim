#!/bin/bash
cd input
ln -sf EventGenerator.in.dummy EventGenerator.in
cd ..

BEAMISOTOPE="BEAMISOTOPE 37 20 20"
BEAMENERGY1="BEAMENERGY "  
BEAMENERGY2=" 0.001"
BEAMPOSITION="BEAMPOSITION 0.0 0.01 0.0 0.01"
BEAMANGLE="BEAMANGLE 0.0 0.01 0.0 360.0"
TARGET="TARGET 1 7.0 7.0 0.01"  
TARGETANGULARBROADENING="TARGETANGULARBROADENING 0 0.00"
MASSCHANGE="MASSCHANGE 0 0"
BORREL="BORREL 0 8."
GOLDHABER="GOLDHABER 0 90."
GAMMAINPUT="GAMMAINPUT ./input/Efficiency.in"
THETARANGE="THETARANGE 0.0 180.0"
NUMBEROFEVENTS="NUMBEROFEVENTS 1000000"
DEFAULTCUTVALUE="DEFAULTCUTVALUE 0.01"
OUTPUTFILE1="OUTPUTFILE ../SimulationResults/Generator/Efficiency/"
OUTPUTFILE2="mev"
OUTPUTFILE3="kev.root"
DEDXTABLE="DEDXTABLE 0 ./dEdXTables/CaOnBedEdX.in ./dEdXTables/CaOnBedEdX.in"
END="END"


FILE="./input/EventGenerator.in"
FILE2="./input/Efficiency.in"

for arg1 in 150  
do
for arg2 in 500 1000 1500 2000 2500 3000 3500 4000 4500 5000   
do

echo "LEVEL 0   0.    0.  0.
LEVEL 1 100. $arg2.  0.
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
