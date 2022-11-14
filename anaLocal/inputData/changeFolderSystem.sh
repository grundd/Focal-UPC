#!/bin/bash
# bash script to rename files in a given name in all subdirectories

for i in $(seq -f "%03g" 1 11)
do
    echo "folder aliDPG_01/kCohJpsiToElRad_${i}_1000ev/"
    mkdir -p aliDPG_01/kCohJpsiToElRad/${i}/
    cd aliDPG_01/kCohJpsiToElRad_${i}_1000ev/
    mv FOCAL.Hits.root galice.root Kinematics.root sim.log ~/alice/analyses/FocalUpc/anaLocal/inputData/aliDPG_01/kCohJpsiToElRad/${i}/
    ls
    cd ..
    cd ..
    echo ""
done