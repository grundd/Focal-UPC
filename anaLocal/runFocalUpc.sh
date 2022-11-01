# shell script must be first allowed: chmod +x run.sh
#!/bin/bash
# to run it do (inside ali shell):
# ./run.sh

option="cohJpsi"

#aliroot -q 'FocalUpc_SimulateClusters.C("sim01/kCohJpsiToElRad_001_1000ev/")'

aliroot -q 'FocalUpc_AnalysisPrimary.C("'$option'")'
aliroot -q 'FocalUpc_AnalysisSecondary.C("'$option'")'
aliroot -q 'FocalUpc_InvMassFit.C("'$option'")'