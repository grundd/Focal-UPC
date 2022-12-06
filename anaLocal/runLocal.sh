# shell script must be first allowed: chmod +x run.sh
#!/bin/bash
# to run it do (inside ali shell):
# ./run.sh

# version as of Nov 24, 2022

process="incPsi2s"
overwrite=kFALSE

# calculate:
# - rapidity dependence of Starlight cross sections
# - expected yields of J/psi, psi' and Y(1S) in Run 4
# - acceptance in feed-down processes
if false; then 
    aliroot -q 'StarlightRapDep.C'
fi

# run primary (grid) and secondary (main) analysis over input data for a selected process
if false; then 
    aliroot -q 'FocalUpcGrid_RunAnalysis.C(kTRUE,"'$process'",'$overwrite')'
    aliroot -q 'AnaMain.C("'$process'")'
fi

# run primary (grid) and secondary (main) analysis over all processes
if false; then 
    for pcs in "cohJpsi" "incJpsi" "cohFD" "incFD" "cohPsi2s" "incPsi2s"
    do
        aliroot -q 'FocalUpcGrid_RunAnalysis.C(kTRUE,"'$pcs'",'$overwrite')'
        aliroot -q 'AnaMain.C("'$pcs'")'
    done
    # make the combined invariant mass fit 
    aliroot -q AnaMain_SignalExtraction.C
fi

# do invariant mass of the combined sample
if true; then
    aliroot -q AnaMain_SignalExtraction.C
fi

# run clusterizer over a specific HITS file
if false; then
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
    aliroot -q 'FocalUpcGrid.C(kTRUE,"cohJpsi",'$overwrite',"inputData/aliDPG_v02/kCohJpsiToElRad/001/","results/sim02_g02_p02/cohJpsi/001/")'
fi

# run primary analysis over a selected input: box electrons
if false; then 
    aliroot -q 'FocalUpcGrid.C(kTRUE,"boxEle",'$overwrite',"inputData/aliDPG_v02/BoxElectrons/001/","results/sim02_g02_p02/boxEle/001/")'
fi