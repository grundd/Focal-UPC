// ConfigAnalysis.h
// David Grund, Nov 09, 2022

// which version of MC simulations:
// (matters only for local analysis)
Int_t simFiles = 2; 
// 1 -> produced before the update on Oct 19, 2022; with geometry_01.txt
// 2 -> produced after the update; with geometry_02.txt
// which geometry and parameters files will be used to run the clusterizer:
// geometry.txt
Int_t fileGeom = 2;
// 1 -> geometry_01.txt (used before the update on Oct 19, 2022)
// 2 -> geometry_02.txt (used after the update)
Int_t filePara = 2;
// parameters.txt
// 1 -> parameters_01.txt: Calib.Const1 set to 12900. (current official version)
// 2 -> parameters_02.txt: Calib.Const1 set to 11220. (old value)

TString sGeomFile = "";
TString sParaFile = "";

// number of available input files for each simulation version (1000 events each):
const Int_t nBoxEle[2] = {16,16};
const Int_t nBoxPho[2] = {6, 0};
const Int_t nCohJpsi[2] = {11,21};
const Int_t nIncJpsi[2] = {0, 1};
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
    else {
        cout << " ERROR: Configuration not supported. Choose between \"boxEle\",\"boxPho\",\"cohJpsi\",\"incJpsi\". Terminating..." << endl;
        return;
    }
    // set the output folder
    outDir = Form("results/sim%02i_g%02i_p%02i/%s/",simFiles,fileGeom,filePara,sim.Data());

    return;
}