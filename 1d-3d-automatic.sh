#!/bin/bash

#important paths
UGSHELL=/home/lreinhardt/Code/ug4/bin/ugshell
ANAMORPH=/home/lreinhardt/Code/anamorph/bin/am_cellgen
APPS=/home/lreinhardt/Code/ug4/apps/cable_neuron_app

#NeuGen, Grids and Geometries
NEUGENBASE=$1
NI=$2

#Parameters

#1. NeuGen Netzwerk erstellen 
#(cli Version?)

#2. UGX-Datei aus NeuGen-TXT erstellen - "import_txt" Methode
$UGSHELL -ex $APPS/neti_import.lua -name $NEUGENBASE -method txt

#3. 3d Neuron auswählen/extrahieren 
#TODO automatisiert nach einem funktionierenden Neuron suchen
#until [ $? -eq 0 ]; do
#$ANAMORPH -i $NEUGENBASE-3d.swc -force-meshing -meshing-cansurf-angularsegments 6 -meshing-#triangle-height 1.0 -preserve-crease-edges -mesh-pp-gec 1.5 0.125 0.5 5 -no-mesh-pp-hc
#done
$UGSHELL -ex $APPS/conversion_ugx2swc.lua -grid $NEUGENBASE.ugx -out $NEUGENBASE-3d.swc

#(3.a SWC-Datei zu UGX konvertieren - "import_geometry_and_generate_grid" Methode
$UGSHELL -ex $APPS/neti_import.lua -name $NEUGENBASE-3d -method swc

#4. 3d Geometrie mit AnaMorph erstellen (Plasmamembran)
$ANAMORPH -i $NEUGENBASE-3d.swc -force-meshing -meshing-cansurf-angularsegments 6 -meshing-triangle-height 1.0 -preserve-crease-edges -mesh-pp-gec 1.5 0.125 0.5 5 -no-mesh-pp-hc

#5. Kleinere 3d Geometrie mit AnaMorph erstellen (ER)
#TODO Wie macht man da nochmal eine kleinere Geometrie?

#6. Plasmamembran und ER mergen, Tetrahedisieren
#TODO ProMesh cli?


