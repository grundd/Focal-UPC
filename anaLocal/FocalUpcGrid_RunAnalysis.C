// FocalUpcGrid_RunAnalysis.C
// David Grund, Nov 07, 2022

// root headers
#include "TSystem.h"
#include "TROOT.h"

Bool_t isLocal = kTRUE;
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

void FocalUpcGrid_RunAnalysis(TString sim = "")
{
    gSystem->Load("libpythia6_4_28.so");
    TString sGeomFile = Form("geometry_%02i.txt",fileGeom);
    TString sParaFile = Form("parameters_%02i.txt",filePara);

    // ******************************************************************************************************************
    // local analysis:
    // ******************************************************************************************************************
    if(isLocal)
    {
        // number of available input files for each simulation version (1000 events each):
        const Int_t nBoxEle[2] = {16,16};
        const Int_t nBoxPho[2] = {6, 0};
        const Int_t nCohJpsi[2] = {11,11};
        const Int_t nIncJpsi[2] = {0, 1};
        Int_t nFiles(0.);
        // input and output folders:
        TString inDir = "";
        TString outDir = "";
        TString outSubDir = "";
        // set the input directory:
        inDir = Form("inputData/sim%02i/",simFiles);
        // box electrons
        if(sim == "boxEle") {
            nFiles = nBoxEle[simFiles-1];
            inDir += "BoxElectrons_";
        }
        // box photons
        else if(sim == "boxPho") {
            nFiles = nBoxPho[simFiles-1];
            inDir += "BoxPhotons_";
        }
        // kCohJpsiToElRad
        else if(sim == "cohJpsi") {
            nFiles = nCohJpsi[simFiles-1];
            inDir += "kCohJpsiToElRad_";
        }
        // kIncohJpsiToElRad
        else if(sim == "incJpsi") {
            nFiles = nIncJpsi[simFiles-1];
            inDir += "kIncohJpsiToElRad_";
        }
        else {
            cout << " ERROR: Configuration not supported. Choose between \"boxEle\",\"boxPho\",\"cohJpsi\",\"incJpsi\". Terminating..." << endl;
            return;
        }
        // set the output folder
        outDir = Form("results/sim%02i_g%02i_p%02i/%s/",simFiles,fileGeom,filePara,sim.Data());
        // is it box simulation?
        Bool_t isBoxSim(kFALSE);
        if(sim == "boxEle" || sim == "boxPho") isBoxSim = kTRUE;
        // compile the macros
        gROOT->ProcessLine(".L ClusterizeGrid.C");
        gROOT->ProcessLine(".L FocalUpcGrid.C");
        // run the analysis
        for(Int_t i = 0; i < nFiles; i++) 
        {
            // input folder ending:
            TString sIn = Form("%s%03i_1000ev/",inDir.Data(),i+1);
            TString sOut = Form("%s%03i/",outDir.Data(),i+1);
            gSystem->Exec("mkdir -p " + sOut);
            // run the clusterizer:
            TString sClFile = Form("%sfocalClusters.root",sOut.Data());
            if(gSystem->AccessPathName(sClFile.Data()))
            {
                cout << " MESSAGE: cluster file not found! Runing the clusterizer now." << endl;
                TString sCmd = Form("ClusterizeGrid(\"%s\",\"%s\",\"%s\",\"%s\")",
                    sGeomFile.Data(),sParaFile.Data(),sIn.Data(),sOut.Data());
                gROOT->ProcessLine(sCmd.Data());
            }
            // run the analysis:
            TString sCmd = Form("FocalUpcGrid(%i,%i,\"%s\",\"%s\")",isLocal,isBoxSim,sIn.Data(),sOut.Data());
            gROOT->ProcessLine(sCmd.Data());
        }
        cout << endl << " FINISHED! " << endl;
        return;
    }

    // ******************************************************************************************************************
    // analysis on Grid:
    // ******************************************************************************************************************
    else 
    {
        // run the clusterizer:
        gROOT->ProcessLine(Form(".x ClusterizeGrid.C(\"%s\",\"%s\")",sGeomFile.Data(),sParaFile.Data()));
        // run the analysis:
        gROOT->ProcessLine(".x FocalUpcGrid_Analysis.C(kFALSE,kFALSE)");
        return;
    }
}