#!/bin/bash
# bash script to rename files in a given name in all subdirectories

for dir in sim02/*/; do
    echo "$dir"
    mv "$dir/focalClusters_v02.root" "$dir/focalClusters_g02_p01.root"
done