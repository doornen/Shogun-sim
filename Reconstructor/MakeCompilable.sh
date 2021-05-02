#!/bin/bash
dummy1="int main"
dummy2="int Dali2Reconstructor"
if [ $1 -eq 1 ] 
then
    cat  Dali2Reconstructor.C | awk -v rep="$dummy1" -v str="$dummy2" '{ gsub(str,rep) ; print }' > temp.txt
    mv temp.txt Dali2Reconstructor.C
fi
if [ $1 -eq 0 ]
then 
    cat  Dali2Reconstructor.C | awk -v rep="$dummy2" -v str="$dummy1" '{ gsub(str,rep) ; print }' > temp.txt
    mv temp.txt Dali2Reconstructor.C
fi
exit 0

