#!/bin/bash

# This case corresponds to the CNRS benchmark on multi-physics codes
# as described in Tiberga et al. (2020)
# https://doi.org/10.1016/j.anucene.2020.107428

ORD="${1:-S6}"

if [[ "$ORD" != "S2" && "$ORD" != "S6" ]]; then
    echo "Error: Invalid SN discretization. Only S2 and S6 are valid" >&2
    exit 1
fi

source $FOAM_SRC/../bin/tools/RunFunctions

cp -r 0.orig 0

cat "constant/neutronFlightDirectionsDict_${ORD}" > constant/neutronFlightDirectionsDict

runApplication blockMesh