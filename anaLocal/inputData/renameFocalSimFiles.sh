#!/bin/bash
# bash script to rename files in a given name in all subdirectories

for dir in sim01/*/; do
    echo "$dir"
    mv "$dir/focalClusters_v1.root" "$dir/focalClusters_v01.root"
done