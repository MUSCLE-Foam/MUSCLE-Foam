#!/bin/bash

echo "Cleaning $PWD case"

zeros=""
while [ ${#zeros} -lt 8 ]
do
    timeDir="0.${zeros}[1-9]*"
    rm -rf ./${timeDir} ./-${timeDir} > /dev/null 2>&1
    zeros="0$zeros"
done

rm -rf ./[1-9]* ./-[1-9]* ./log* ./log.* ./logSummary.* \
    ./*.OpenFOAM ./*.blockMesh ./*.foam

rm -rf 0
rm -rf constant/polyMesh
rm -rf postProcessing*
rm -f *.pdf
rm -f *.png
rm -f *.txt
rm -f *.obj
rm -f *.msh
rm -rf processor*
rm -fr __pycache__
rm -fr plots/__pycache__
rm -fr plots/*.png
rm -f slurm*
rm -f yplus*
rm -f wallShearStress