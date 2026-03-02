#!/bin/bash

ORD="${1:-S6}"

if [[ "$ORD" != "S2" && "$ORD" != "S6" ]]; then
    echo "Error: Invalid SN discretization. Only S2 and S6 are valid" >&2
    exit 1
fi

source $FOAM_SRC/../bin/tools/RunFunctions

cp -r 0.orig 0
cat "constant/neutronFlightDirectionsDict_${ORD}" > constant/neutronFlightDirectionsDict

runApplication blockMesh