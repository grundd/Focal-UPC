# shell script must be first allowed: chmod +x runBoxSim.sh
#!/bin/bash
# to run it do (inside ali shell):
# ./runBoxSim.sh

nSim=10
nEvPerSim=1000
pdgCode=11 # electrons 11, photons 22

echo "${nSim} simulation(s) will be performed with ${nEvPerSim} event per each"

#for i in {1..$nSim}
# https://stackoverflow.com/questions/8789729/how-to-zero-pad-a-sequence-of-integers-in-bash-so-that-all-have-the-same-width
for i in $(seq -f "%03g" 1 $nSim)
do
    echo "Running simmulation ${i}"
    # run the simulation:
    ranSeed=$((($RANDOM*$RANDOM*$RANDOM*$RANDOM) % 1000000000 + 1))
    echo "Random seed: ${ranSeed}"
    $ALIDPG_ROOT/bin/aliroot_dpgsim.sh --mode sim --run 294925 --uid $ranSeed --generator Upgrade:FOCAL_Generators:box \
    --nevents $nEvPerSim --simulation NoDigitization --detector FOCAL --focalGeometryFile geometry_02.txt --etamin 3.4 --etamax 5.8 \
    --ptmin 0.0 --ptmax 2.0 --pdg $pdgCode
    # create the folder to store the results:
    #folderName="BoxPhotons_${i}_${nEvPerSim}ev"
    folderName="BoxElectrons_${i}_${nEvPerSim}ev"
    mkdir -p $folderName
    mv FOCAL.Hits.root galice.root Kinematics.root sim.log $folderName
    # delete unwanted files:    
    rm -r GRP/
    #rm -r ACORDE.Digits.root
    rm -r ACORDE.Hits.root
    #rm -r EMCAL.Digits.root
    rm -r EMCAL.Hits.root
    #rm -r EMCAL.SDigits.root
    #rm -r FIT.Digits.root
    rm -r FIT.Hits.root
    rm -r geometry.root
    rm -r gphysi.dat
    rm -r grpdump.sh
    #rm -r HMPID.Digits.root
    rm -r HMPID.Hits.root
    #rm -r HMPID.SDigits.root
    #rm -r ITS.Digits.root
    rm -r ITS.Hits.root
    rm -r MCStepLoggerVolMap.dat
    #rm -r MUON.Digits.root
    rm -r MUON.Hits.root
    #rm -r MUON.SDigits.root
    #rm -r PHOS.Digits.root
    rm -r PHOS.HITS.root
    #rm -r PHOS.SDigits.root
    #rm -r PMD.Digits.root
    rm -r PMD.Hits.root
    #rm -r PMD.SDigits.root
    rm -r QA.root
    rm -r simwatch.log
    #rm -r TOF.Digits.root
    rm -r TOF.Hits.root
    #rm -r TOF.SDigits.root
    #rm -r TPC.Digits.root
    rm -r TPC.Hits.root
    #rm -r TRD.Digits.root
    rm -r TRD.Hits.root
    #rm -r TRD.SDigits.root
    rm -r Trigger.root
done