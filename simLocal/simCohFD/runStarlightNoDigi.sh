# shell script must be first allowed: chmod +x runJpsiSim.sh
#!/bin/bash
# to run it do (inside ali shell):
# ./runStarlightNoDigi.sh

# options to set:
first=13
simulations=8
last=$(($first+$simulations-1))
evPerSim=1000
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
    $ALIDPG_ROOT/bin/aliroot_dpgsim.sh --run 294925 --system Pb-Pb --energy 5500.0 --mode sim --detector FOCAL --uid $ranSeed \
    --nevents $evPerSim --generator Starlight --process $process --simulation NoDigitization --focalGeometryFile geometry_02.txt \
    --ymin 3.4 --ymax 6.0
    # create a folder to store the results:
    mkdir -p $folderName
    mv FOCAL.Hits.root galice.root Kinematics.root sim.log $folderName
    # delete unwanted files:
    rm -r GRP/
    rm -r ACORDE.Hits.root
    rm -r EMCAL.Hits.root
    rm -r FIT.Hits.root
    rm -r geometry.root
    rm -r gphysi.dat
    rm -r grpdump.sh
    rm -r HMPID.Hits.root
    rm -r ITS.Hits.root
    rm -r MCStepLoggerVolMap.dat
    rm -r MUON.Hits.root
    rm -r PHOS.HITS.root
    rm -r PMD.Hits.root
    rm -r QA.root
    rm -r simwatch.log
    rm -r slight.txt
    rm -r TOF.Hits.root
    rm -r TPC.Hits.root
    rm -r TRD.Hits.root
    rm -r Trigger.root
done