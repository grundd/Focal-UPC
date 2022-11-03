// FocalUpc_ConfigAnalysis.h
// David Grund, Oct 31, 2022

#include "TString.h"
#include "TSystem.h"

// *****************************************************************************
// options to set:
Bool_t plotEvInfo = kFALSE;
Bool_t printEvInfo = kFALSE;
Bool_t matchDirectly = kFALSE;
// superclusterizer:
Bool_t doSupercls = kTRUE; // if true, create superclusters using:
const Float_t minSeedE = 5; // [GeV]
const Float_t radius = 15; // [cm]
// selections:
const Float_t cutM = 0.0; // [GeV]; if > 0, filter out all (super)cluster pairs having mass below cutM
const Float_t cutE = 30.0; // [GeV]; if > 0, filter out all (super) clusters having energy below cutE
// *****************************************************************************
// input files and configuration:
Int_t vSim = 2;
const Int_t nInputBoxEle[2] = {16,16};
const Int_t nInputBoxPho[2] = {6, 0};
const Int_t nInputCohJpsi[2] = {11,11};
const Int_t nInputIncJpsi[2] = {0, 1};
// version of input simulations:
// 1 -> files produced before Oct 19, 2022
// 2 -> files produced after Oct 19, 2022, after an update of FoCal&AliRoot
Int_t vGeomFile = 2;
// version of the FoCal software and the "geometry.txt" file using which the clusters were produced 
// 1 -> "geometry_01.txt" (used before the update on Oct 19, 2022)
// 2 -> "geometry_02.txt" (used after the update on Oct 19, 2022)
Int_t vParaFile = 2;
// version of the "parameters.txt" file using which the clusters were produced
// 1 -> "parameters_01.txt": Calib.Const1 set to 12900. (currently the official version)
// 2 -> "parameters_02.txt": Calib.Const1 set to 11220. (old value)
// *****************************************************************************

Int_t nInputFiles(0.);
TString sInputFiles = "";
Bool_t isBoxSim(kFALSE);
TString sOut = "";

Bool_t ConfigAnalysis(TString sSim)
{
    gSystem->Load("libpythia6_4_28.so");

    // box electrons
    if(sSim == "boxEle") {
        nInputFiles = nInputBoxEle[vSim-1];
        sInputFiles = "BoxElectrons_";
    }
    // box photons
    else if(sSim == "boxPho") {
        nInputFiles = nInputBoxPho[vSim-1];
        sInputFiles = "BoxPhotons_";
    }
    // kCohJpsiToElRad
    else if(sSim == "cohJpsi") {
        nInputFiles = nInputCohJpsi[vSim-1];
        sInputFiles = "kCohJpsiToElRad_";
    }
    // kIncohJpsiToElRad
    else if(sSim == "incJpsi") {
        nInputFiles = nInputIncJpsi[vSim-1];
        sInputFiles = "kIncohJpsiToElRad_";
    }
    else {
        cout << "Configuration not supported. Choose between sSim = \"boxEle\",\"boxPho\",\"cohJpsi\",\"incJpsi\". Terminating..." << endl;
        return kFALSE;
    }
    
    if(sSim == "boxEle" || sSim == "boxPho") isBoxSim = kTRUE;

    // mass filtering only for J/psi simulations
    if(cutM > 0 && isBoxSim) {
        cout << "Cannot do mass cleaning for box simulations of electrons/photons. Terminating... " << endl;
        return kFALSE;
    }

    // create output folder
    sOut = Form("results/sim%02i_g%02i_p%02i/",vSim,vGeomFile,vParaFile);
    sOut += "_";
    sOut += sSim;
    if(doSupercls) sOut += Form("_supCl");
    if(!isBoxSim && matchDirectly) sOut += Form("_dirMtch");
    if(cutM > 0)   sOut += Form("_cutM%.1f",cutM);
    if(cutE > 0)   sOut += Form("_cutE%.1f",cutE);
    sOut += "/";
    gSystem->Exec("mkdir -p " + sOut);

    return kTRUE;
}