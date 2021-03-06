#!/bin/bash

INPUTFILE="INPUTFILE ../SimulationResults/Generator/32ne/33na_32ne260mev_700kev"
FILEEXT1="ps.root"
OUTPUTFILE="OUTPUTFILE ../SimulationResults/Builder/dali2/32ne/33na_32ne260mev_700kev_nores"
DALI2INCLUDE="DALI2INCLUDE 1"
DALI2ENERGYRESOLUTION="DALI2ENERGYRESOLUTION 2 0.0 0.5"
POSDETECTORONTARGETRESOLUTION="POSDETECTORONTARGETRESOLUTION .3"
ENERGYDETECTORAFTERTARGETINCLUDE="ENERGYDETECTORAFTERTARGETINCLUDE 0"
POSDETECTORAFTERTARGETDISTANCE="POSDETECTORAFTERTARGETDISTANCE 100."
POSDETECTORAFTERTARGETRESOLUTION="POSDETECTORAFTERTARGETRESOLUTION .3"
BETARESOLUTION="BETARESOLUTION 0.001"
BEAMPIPEINCLUDE="BEAMPIPEINCLUDE 1"
TARGETHOLDERINCLUDE="TARGETHOLDERINCLUDE 1"
END="END"

FILE="./input/EventBuilder.in"

for arg1 in 0 5 10 15 20 25 30 35 40 45 50
do

echo ""$INPUTFILE$arg1$FILEEXT1"
"$OUTPUTFILE$arg1$FILEEXT1"
"$DALI2INCLUDE"
"$DALI2ENERGYRESOLUTION"
"$POSDETECTORONTARGETRESOLUTION"
"$ENERGYDETECTORAFTERTARGETINCLUDE"
"$POSDETECTORAFTERTARGETDISTANCE"
"$POSDETECTORAFTERTARGETRESOLUTION"
"$BETARESOLUTION"
"$BEAMPIPEINCLUDE"
"$TARGETHOLDERINCLUDE"
"$END" 
" > $FILE

time EventBuilder run_nothing.mac
done

exit 0
