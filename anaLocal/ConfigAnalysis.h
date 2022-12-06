// ConfigAnalysis.h
// David Grund, Nov 09, 2022

const Int_t nSim = 3;
// which version of MC simulations:
// (matters only for local analysis)
Int_t simFiles = 3; 
// 1 -> produced before the update on Oct 19, 2022; geometry_01.txt
// 2 -> produced before the update on Nov 30, 2022; geometry_02.txt
// 3 -> produced after the update on Nov 30; geometry_03.txt
// which geometry and parameters files will be used to run the clusterizer:
// geometry.txt
Int_t fileGeom = 3;
// 1 -> geometry_01.txt (before the update on Oct 19, 2022)
// 2 -> geometry_02.txt (before the update on Nov 30, 2022)
// 3 -> geometry_03.txt (after the update on Nov 30)
Int_t filePara = 3;
// parameters.txt
// 1 -> parameters_01.txt: Calib.Const1 set to 12900. (current official version)
// 2 -> parameters_02.txt: Calib.Const1 set to 11220. (old value)
// 3 -> parameters_03.txt: Calib.Const1 set to 6000.

TString sGeomFile = "";
TString sParaFile = "";

// number of available input files for each simulation version (1000 events each):
const Int_t nBoxEle[nSim] = {16, 16, 0};
const Int_t nBoxPho[nSim] = {6, 0, 0};
const Int_t nCohJpsi[nSim] = {11, 97, 70};
const Int_t nIncJpsi[nSim] = {0,  72, 50};
const Int_t nCohFD[nSim] = {0, 18, 10};
const Int_t nIncFD[nSim] = {0, 18, 10};
const Int_t nBkg[nSim]   = {0, 13, 0};
const Int_t nCohPsi2s[nSim] = {0, 24, 17};
const Int_t nIncPsi2s[nSim] = {0, 20, 15};
const Int_t nCohJpsiNoFIT[nSim] = {0, 10, 0};
Int_t nFiles(0.);
TString inDir = "";
TString outDir = "";

void ConfigLocalAnalysis(TString sim)
{
    sGeomFile = Form("geometry_%02i.txt",fileGeom);
    sParaFile = Form("parameters_%02i.txt",filePara);

    // set the input directory:
    inDir = Form("inputData/aliDPG_v%02i/",simFiles);

    // box electrons
    if(sim == "boxEle") {
        nFiles = nBoxEle[simFiles-1];
        inDir += "BoxElectrons/";
    }
    // box photons
    else if(sim == "boxPho") {
        nFiles = nBoxPho[simFiles-1];
        inDir += "BoxPhotons/";
    }
    // kCohJpsiToElRad
    else if(sim == "cohJpsi") {
        nFiles = nCohJpsi[simFiles-1];
        inDir += "kCohJpsiToElRad/";
    }
    // kIncohJpsiToElRad
    else if(sim == "incJpsi") {
        nFiles = nIncJpsi[simFiles-1];
        inDir += "kIncohJpsiToElRad/";
    }
    // kCohPsi2sToElPi
    else if(sim == "cohFD") {
        nFiles = nCohFD[simFiles-1];
        inDir += "kCohPsi2sToElPi/";
    }
    // kIncohPsi2sToElPi
    else if(sim == "incFD") {
        nFiles = nIncFD[simFiles-1];
        inDir += "kIncohPsi2sToElPi/";
    }
    // kTwoGammaToElMedium
    else if(sim == "bkg") {
        nFiles = nBkg[simFiles-1];
        inDir += "kTwoGammaToElMedium_noEtaCut/";
    }
    // kCohPsi2sToEl
    else if(sim == "cohPsi2s") {
        nFiles = nCohPsi2s[simFiles-1];
        inDir += "kCohPsi2sToEl/";
    }
    // kIncohPsi2sToEl
    else if(sim == "incPsi2s") {
        nFiles = nIncPsi2s[simFiles-1];
        inDir += "kIncohPsi2sToEl/";
    }
    // kCohJpsiToElRad without FIT
    else if(sim == "cohJpsiNoFIT") {
        nFiles = nCohJpsiNoFIT[simFiles-1];
        inDir += "kCohJpsiToElRadNoFIT/";
    }
    else {
        cout << " ERROR: Configuration not supported. Choose between:\n"
             << "\t\"boxEle\",\n"
             << "\t\"boxPho\",\n"
             << "\t\"cohJpsi\",\n"
             << "\t\"incJpsi\",\n"
             << "\t\"cohFD\",\n"
             << "\t\"incFD\",\n"
             << "\t\"bkg\",\n"
             << "\t\"cohPsi2s\",\n"
             << "\t\"incPsi2s\".\n"
             << "\t\"cohJpsiNoFIT\",\n"
             << "\tTerminating...\n";
        return;
    }
    // set the output folder
    outDir = Form("results/sim%02i_g%02i_p%02i/%s/",simFiles,fileGeom,filePara,sim.Data());

    return;
}