#!/bin/bash

source $FOAM_SRC/../bin/tools/RunFunctions
source $FOAM_SRC/../bin/tools/CleanFunctions

cleanCase
rm -rf 0
rm -rf postProcessing/*
rm -f *.pdf
rm -f *.png
rm -f *.txt
rm -f *.obj
rm -f log*
rm -f *.msh
rm -rf processor*
rm -fr __pycache__
rm -f slurm*
rm -f constant/fvModels
rm -f system/fvSolution
rm -f system/controlDict