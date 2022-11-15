# shell script must be first allowed: chmod +x runJpsiSim.sh
#!/bin/bash
# to run it do (inside ali shell):
# ./runStarlightGenOnly.sh

# options to set:
first=1
simulations=1
last=$(($first+$simulations-1))
evPerSim=100000
process="kCohPsi2sToElPi"

# print info:
echo ""
echo "******************************************"
echo "Process ${process}:"
echo " -> ${simulations} simulation(s) will be performed"
echo " -> ${evPerSim} event per each simulation"
echo "******************************************"
echo ""

#for i in {1..$simulations}
# https://stackoverflow.com/questions/8789729/how-to-zero-pad-a-sequence-of-integers-in-bash-so-that-all-have-the-same-width
for i in $(seq -f "%03g" $first $last)
do
    # get random seed:
    ranSeed=$((($RANDOM*$RANDOM*$RANDOM*$RANDOM) % 1000000000 + 1))
    folderName="${process}/${i}"
    # print info:
    echo ""
    echo "******************************************"
    echo "Running simmulation ${i}:"
    echo " -> random seed: ${ranSeed}"
    echo " -> results stored in ${folderName}"
    echo "******************************************"
    echo ""
    # run the simulation:
    $ALIDPG_ROOT/bin/aliroot_dpgsim.sh --run 294925 --system Pb-Pb --energy 5500.0 --mode sim --uid $ranSeed \
    --nevents $evPerSim --generator Starlight --process kCohPsi2sToElPi --simulation GeneratorOnly --ymin 3.4 --ymax 6.0
    # create a folder to store the results:
    mkdir -p $folderName
    mv galice.root Kinematics.root sim.log $folderName
    # delete unwanted files:
    rm -r grpdump.sh
    rm -r QA.root
    rm -r simwatch.log
    rm -r slight.txt
done