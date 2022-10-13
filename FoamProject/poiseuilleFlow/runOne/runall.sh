#!/bin/bash

# Copy mesh file to run directory
cp ./../../../../Meshing/flat_plate_in_channel.msh .

# Generate mesh
blockMesh

# Run case
simpleFoam
