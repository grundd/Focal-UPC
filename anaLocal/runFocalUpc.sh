# shell script must be first allowed: chmod +x run.sh
#!/bin/bash
# to run it do (inside ali shell):
# ./run.sh

# version as of Nov 07, 2022

# run primary analysis over all input data: box electrons
sim="cohJpsi"
if true; then 
    aliroot -q 'FocalUpcGrid_RunAnalysis.C("'$sim'")'
fi

# plot event displays: cohJpsi
if false; then 
    aliroot -q 'FocalUpcEventDisplay.C("cohJpsi")'
fi

# plot event displays: box electrons
if false; then 
    aliroot -q 'FocalUpcEventDisplay.C("boxEle")'
fi

# run primary analysis over a selected input: cohJpsi
if false; then 
    aliroot -q 'FocalUpcGrid.C(kTRUE,kFALSE,"inputData/sim02/kCohJpsiToElRad_001_1000ev/","results/sim02_g02_p02/cohJpsi/001/")'
fi

# run primary analysis over a selected input: box electrons
if false; then 
    aliroot -q 'FocalUpcGrid.C(kTRUE,kTRUE,"inputData/sim02/BoxElectrons_001_1000ev/","results/sim02_g02_p02/boxEle/001/")'
fi