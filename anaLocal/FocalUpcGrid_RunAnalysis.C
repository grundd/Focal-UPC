// FocalUpcGrid_RunAnalysis.C
// David Grund, Nov 07, 2022

// root headers
#include "TSystem.h"
#include "TROOT.h"
// my headers
#include "ConfigAnalysis.h"

void FocalUpcGrid_RunAnalysis(Bool_t isLocal, TString sim)
{
    gSystem->Load("libpythia6_4_28.so");

    // ******************************************************************************************************************
    // local analysis:
    // ******************************************************************************************************************
    if(isLocal)
    {
        ConfigLocalAnalysis(sim);
        // compile the macros
        gROOT->ProcessLine(".L ClusterizeGrid.C");
        gROOT->ProcessLine(".L FocalUpcGrid.C");
        // run the analysis
        for(Int_t i = 0; i < nFiles; i++) 
        {
            // ending of input and output directories:
            TString sIn = Form("%s%03i/",inDir.Data(),i+1);
            TString sOut = Form("%s%03i/",outDir.Data(),i+1);
            gSystem->Exec("mkdir -p " + sOut);
            // run the clusterizer:
            TString sClFile = Form("%sfocalClusters.root",sOut.Data());
            if(gSystem->AccessPathName(sClFile.Data()))
            {
                cout << " MESSAGE: cluster file not found! Running the clusterizer now." << endl;
                TString sCmd = Form("ClusterizeGrid(kTRUE,\"%s\",\"%s\",\"%s\",\"%s\")",
                    sGeomFile.Data(),sParaFile.Data(),sIn.Data(),sOut.Data());
                gROOT->ProcessLine(sCmd.Data());
            }
            // run the analysis:
            TString sCmd = Form("FocalUpcGrid(kTRUE,\"%s\",kFALSE,\"%s\",\"%s\")",sim.Data(),sIn.Data(),sOut.Data());
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
        gROOT->ProcessLine(Form(".x ClusterizeGrid.C(kFALSE,\"%s\",\"%s\")",sGeomFile.Data(),sParaFile.Data()));
        // run the analysis:
        gROOT->ProcessLine(Form(".x FocalUpcGrid.C(kFALSE,\"%s\")",sim.Data()));
        return;
    }
}