# shell script must be first allowed: chmod +x run.sh
#!/bin/bash
# to run it do (inside ali shell):
# ./run.sh

# version as of Nov 09, 2022

sim="cohJpsiNoFIT"

# calculate:
# - rapidity dependence of Starlight cross sections
# - expected yields of J/psi, psi' and Y(1S) in Run 4
# - acceptance in feed-down processes
if false; then 
    aliroot -q 'StarlightRapDep.C'
fi

# run primary (grid) and secondary (main) analysis
# over all input data for a selected process
if false; then 
    aliroot -q 'FocalUpcGrid_RunAnalysis.C(kTRUE,"'$sim'")'
    aliroot -q 'AnaMain.C("'$sim'")'
fi

# run clusterizer over a specific HITS file
if true; then
    aliroot -q 'ClusterizeGrid.C(kTRUE,"geometry_02.txt","parameters_02.txt","inputData/aliDPG_v02/kCohJpsiToElRad/062/")'
fi

# plot event displays: cohJpsi
if false; then 
    aliroot -q 'EventDisplay.C("cohJpsi")'
fi

# plot event displays: box electrons
if false; then 
    aliroot -q 'EventDisplay.C("boxEle")'
fi

# run primary analysis over a selected input: cohJpsi
if false; then 
    aliroot -q 'FocalUpcGrid.C(kTRUE,kFALSE,"inputData/sim02/kCohJpsiToElRad_001_1000ev/","results/sim02_g02_p02/cohJpsi/001/")'
fi

# run primary analysis over a selected input: box electrons
if false; then 
    aliroot -q 'FocalUpcGrid.C(kTRUE,kTRUE,"inputData/sim02/BoxElectrons_001_1000ev/","results/sim02_g02_p02/boxEle/001/")'
fi