#!/bin/bash

# This case corresponds to the CNRS benchmark on multi-physics codes
# as described in Tiberga et al. (2020)
# https://doi.org/10.1016/j.anucene.2020.107428

# Usage: ./cavityBenchmark.sh <case> <neutronicsSolver>
# Example: ./cavityBenchmark.sh 0.2 SP3
# neutronicsSolver defaults to diffusion if no second argument is provided

source $FOAM_SRC/../bin/tools/RunFunctions

CASE="$1"
SOLVER="${2:-diffusion}"

cp -r 0.orig 0

runApplication blockMesh

case $CASE in

    0.1)
        echo "Case 0.1: Velocity field"
        sed -i 's/VARVEL/0.5/g' 0/U
        cat system/fvSolution.orig > system/fvSolution
        sed -i 's/VAROUTERCORR/1/g' system/fvSolution
        sed -i 's/VARUMAXITER/1000/g' system/fvSolution
        sed -i 's/VARHMAXITER/0/g' system/fvSolution
        sed -i 's/VARNOUTERCORR/0/g' system/fvSolution
        ;;
    0.2)
        echo "Case 0.2: Neutronics"
        sed -i 's/VARVEL/0/g' 0/U
        cat system/fvSolution.orig > system/fvSolution
        sed -i 's/VAROUTERCORR/0/g' system/fvSolution
        sed -i 's/VARNOUTERCORR/1/g' system/fvSolution
        ;;
    1.1)
        echo "Case 1.1: Circulating fuel"
        sed -i 's/VARVEL/0.5/g' 0/U
        cat system/fvSolution.orig > system/fvSolution
        sed -i 's/VAROUTERCORR/1/g' system/fvSolution
        sed -i 's/VARUMAXITER/1000/g' system/fvSolution
        sed -i 's/VARHMAXITER/0/g' system/fvSolution
        sed -i 's/VARNOUTERCORR/1/g' system/fvSolution
        ;;
    1.3)
        echo "Case 1.3: Buoyancy"
        sed -i 's/VARVEL/0/g' 0/U
        cat system/fvSolution.orig > system/fvSolution
        sed -i 's/VAROUTERCORR/1/g' system/fvSolution
        sed -i 's/VARUMAXITER/1000/g' system/fvSolution
        sed -i 's/VARHMAXITER/1000/g' system/fvSolution
        sed -i 's/VARNOUTERCORR/1/g' system/fvSolution
        cat constant/neutronicsHeatSource > constant/fvModels
        sed -i 's/VARGAMMA/1e6/g' constant/fvModels
        sed -i 's/VARTREF/900/g' constant/fvModels
        ;;
    1.4)
        echo "Case 1.4: Full coupling"
        sed -i 's/VARVEL/0.5/g' 0/U
        cat system/fvSolution.orig > system/fvSolution
        sed -i 's/VAROUTERCORR/1/g' system/fvSolution
        sed -i 's/VARUMAXITER/1000/g' system/fvSolution
        sed -i 's/VARHMAXITER/1000/g' system/fvSolution
        sed -i 's/VARNOUTERCORR/1/g' system/fvSolution
        cat constant/neutronicsHeatSource > constant/fvModels
        sed -i 's/VARGAMMA/1e6/g' constant/fvModels
        sed -i 's/VARTREF/900/g' constant/fvModels
        ;;
    *)
        echo "Unknown or not implemented case"
        exit
        ;;
esac

cat system/controlDict.orig > system/controlDict
sed -i "s/VARSOLVER/${SOLVER}/g" system/controlDict

muscleFoamRun > log &