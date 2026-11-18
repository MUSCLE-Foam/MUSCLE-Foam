
#!/bin/bash

# Create 0 folder
cp -r 0.orig 0

# Mesh
blockMesh > log.blockMesh
createZones > log.createZones
setFields > log.setFields