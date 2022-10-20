#!/bin/bash

# Copies FOCAL simulation outputs from grid
# The user has to specify the input directory wrt alienFocalSimDir and the output directory in his/her own local filesystem
# In the output directory, the same structure as in the simulation inputDir will be replicated
#
# Usage example:
# ./copyFromGrid.sh /alice/cern.ch/user/f/focal/sim/v1.1_root6/boxPhotonG4 /home/myuser/data/ root_archive.zip 20

inputDir=$1
outDir=$2
filename=$3
nthreads=$4

mkdir -p $outDir/$inputDir
alien_ls alien://$inputDir/ > $outDir/chunkList.txt
nchunks=`cat $outDir/chunkList.txt | wc -l`
echo "Copying $nchunks files to $outDir from grid"
chunkList=`cat $outDir/chunkList.txt`

for chunk in $chunkList; do
  
   while [ $(ps -ef | grep alien_cp | wc -l) -gt $nthreads ]; do
     sleep 5s;
   done
    
   mkdir -p $outDir/$chunk
   echo "Copy chunk $chunk"
   if [ -s $outDir/$chunk/$filename ]
      then
         echo "    -->  Chunk exists already"
      else
         alien_cp alien://$inputDir/$chunk/$filename file://$outDir/$chunk/  &
      fi
   done
done
   
