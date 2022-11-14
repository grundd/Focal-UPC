# shell script must be first allowed: chmod +x runJpsiSim.sh
#!/bin/bash
# to run it do (inside ali shell):
# ./runGeneratorOnly.sh

sSim="kCohPsi2sToElPi"
nSim=1
nEvPerSim=100000

echo "******************************************"
echo "Process ${kCohPsi2sToElPi}:"
echo "${nSim} simulation(s) will be performed"
echo "${nEvPerSim} event per each simulation"
echo "******************************************"

#for i in {1..$nSim}
# https://stackoverflow.com/questions/8789729/how-to-zero-pad-a-sequence-of-integers-in-bash-so-that-all-have-the-same-width
for i in $(seq -f "%03g" 1 $nSim)
do
    echo "Running simmulation ${i}"
    # get random seed:
    ranSeed=$((($RANDOM*$RANDOM*$RANDOM*$RANDOM) % 1000000000 + 1))
    echo "Rnd seed: ${ranSeed}"
    # run the simulation:
    $ALIDPG_ROOT/bin/aliroot_dpgsim.sh --run 294925 --system Pb-Pb --energy 5500.0 --mode sim --uid $ranSeed \
    --nevents $nEvPerSim --generator Starlight --process kCohPsi2sToElPi --simulation GeneratorOnly --ymin 3.4 --ymax 6.0
    # create a folder to store the results:
    folderName="${sSim}_${nEvPerSim}ev/${i}"
    mkdir -p $folderName
    mv galice.root Kinematics.root sim.log $folderName
    # delete unwanted files:
    rm -r grpdump.sh
    rm -r QA.root
    rm -r simwatch.log
    rm -r slight.txt
done