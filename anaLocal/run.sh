# shell script must be first allowed: chmod +x run.sh
#!/bin/bash
# to run it do (inside ali shell):
# ./run.sh

#aliroot -q 'AnalysisFocalUPC_Jpsi_v2.C("09-26-2022_kCohJpsiToElRad_1000ev/", 10)'
#aliroot -q 'AnalysisFocalUPC_Jpsi_v2.C("10-05-2022_kCohJpsiToElRad_3000ev/", 10)'
#aliroot -q 'AnalysisFocalUPC_Jpsi_v2.C+("10-12-2022_BoxElectrons_1000ev/",kTRUE)'
aliroot -q 'FocalUpcAnalysis_SimulateClusters.C("sim01/kCohJpsiToElRad_001_1000ev/")'
#aliroot -q 'FocalUpcAnalysis.C+(kTRUE)'
