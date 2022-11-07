// FocalUpcEventDisplay.C
// David Grund, Nov 07, 2022

// root headers
#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
// aliroot headers
#include "AliRunLoader.h"
// focal headers
#include "AliFOCALClusterizerv2.h"
// my headers
#include "FocalUpcEventDisplay.h"

// which version of MC simulations:
// (matters only for local analysis)
Int_t simFiles = 2; 
// 1 -> produced before the update on Oct 19, 2022; with geometry_01.txt
// 2 -> produced after the update; with geometry_02.txt
// dataset number:
Int_t nSim = 1;
// which geometry and parameters files will be used to run the clusterizer:
// geometry.txt
Int_t fileGeom = 2;
// 1 -> geometry_01.txt (used before the update on Oct 19, 2022)
// 2 -> geometry_02.txt (used after the update)
Int_t filePara = 2;
// parameters.txt
// 1 -> parameters_01.txt: Calib.Const1 set to 12900. (current official version)
// 2 -> parameters_02.txt: Calib.Const1 set to 11220. (old value)

void PlotEventDisplays(TString, TString, TString, Bool_t);

void FocalUpcEventDisplay(TString sim)
{
    // set the input folder:
    TString inDir = Form("inputData/sim%02i/",simFiles);
    // box electrons
    if(sim == "boxEle") inDir += "BoxElectrons_";
    // box photons
    else if(sim == "boxPho") inDir += "BoxPhotons_";
    // kCohJpsiToElRad
    else if(sim == "cohJpsi") inDir += "kCohJpsiToElRad_";
    // kIncohJpsiToElRad
    else if(sim == "incJpsi") inDir += "kIncohJpsiToElRad_";
    else {
        cout << " ERROR: Configuration not supported. Choose between \"boxEle\",\"boxPho\",\"cohJpsi\",\"incJpsi\". Terminating..." << endl;
        return;
    }
    inDir += Form("%03i_1000ev/",nSim);
    // set the output folder and subfolder:
    TString outDir = Form("results/sim%02i_g%02i_p%02i/%s/%03i/",simFiles,fileGeom,filePara,sim.Data(),nSim);
    TString outSubDir = "eventDisplay/";  
    gSystem->Exec("mkdir -p " + outDir + outSubDir);
    // plot event displays:
    PlotEventDisplays(inDir,outDir,outSubDir,kTRUE);
}

void PlotEventDisplays(TString inDir = "", TString outDir = "", TString outSubDir = "", Bool_t test = kFALSE)
{
    gSystem->Load("libpythia6_4_28.so");
    TString sGeomFile = Form("geometry_%02i.txt",fileGeom);
    TString sParaFile = Form("parameters_%02i.txt",filePara);

    // if clusters not yet produced, run the clusterizer
    TString sClFile = Form("%sfocalClusters.root",outDir.Data());
    if(gSystem->AccessPathName(sClFile.Data()))
    {
        cout << " MESSAGE: cluster file not found! Runng the clusterizer now." << endl;
        TString sCmd = Form(".x ClusterizeGrid.C(\"%s\",\"%s\",\"%s\",\"%s\")",
            sGeomFile.Data(),sParaFile.Data(),inDir.Data(),outDir.Data());
        gROOT->ProcessLine(sCmd.Data());
    }
    // check if focalClusters.root were produced properly
    if(gSystem->AccessPathName(sClFile.Data()))
    {
        cout << " ERROR: cluster file produced but not found! Terminating." << endl;
        return;
    } 
    // open the file with clusters
    TFile* fCls = new TFile(sClFile.Data());
    cout << " MESSAGE: Loading clusters from: " << fCls->GetName() << endl;

    // load the clusterizer
    AliFOCALClusterizerv2* clusterizer = new AliFOCALClusterizerv2();
    clusterizer->InitGeometry(sGeomFile.Data());
    AliFOCALGeometry* geometry = clusterizer->GetGeometry();

    // define ALICE run loader: open galice.root
    AliRunLoader* runLoader = AliRunLoader::Open(inDir + "galice.root");
    if(!runLoader) 
    {
        cout << " ERROR: AliRunLoader not good! Terminating." << endl;
        return;   
    }
    if(!runLoader->GetAliRun()) runLoader->LoadgAlice();
    if(!runLoader->TreeE()) runLoader->LoadHeader();
    if(!runLoader->TreeK()) runLoader->LoadKinematics();

    Int_t nEv(0.);
    if(test) nEv = 10;
    else nEv = runLoader->GetNumberOfEvents();
    cout << "Plotting event displays of " << nEv << " events:" << endl;
    // loop over MC events
    Float_t progress = 0.; // perc
    for(Int_t iEv = 0; iEv < nEv; iEv++) 
    {
        cout << "Ev " << iEv+1 << ":" << endl;
        // update the progress bar
        if((iEv+1) % (Int_t)(nEv/10.) == 0) {
            progress += 10.;
            cout << "[" << progress << "%] done." << endl;
        }

        // get current MC event
        Int_t isEventOk = runLoader->GetEvent(iEv);
        // the method GetEvent(i) returns zero if event is loaded succesfully, if not we continue
        if (isEventOk != 0) {
            cout << "Ev " << iEv+1 << " not OK, skipping." << endl;
            continue;
        }

        // get the stack of MC tracks for this event
        AliStack* stack = runLoader->Stack();

        // get a tree with clusters for this event 
        // (separate tree for each event in subfolders "Event0, Event1, ...")
        TTree* tCls = NULL;
        if(fCls->GetDirectory(Form("Event%i",iEv))) fCls->GetDirectory(Form("Event%i",iEv))->GetObject("fTreeR",tCls);
        else {
            cout << "  (!) cannot find Ev " << iEv+1 << " in a cluster file " << fCls->GetName() << ". Skipping..." << endl;
            fCls->ls();
            continue;
        }

        // get the branch with segment-by-segment clusters
        TBranch* bClsSeg = NULL;
        bClsSeg = tCls->GetBranch("AliFOCALClusterItr");
        TClonesArray* arrClsSeg = NULL;
        bClsSeg->SetAddress(&arrClsSeg);
        bClsSeg->GetEvent(0);
        Int_t nClsSeg = arrClsSeg->GetEntries();

        // get the branch with final summed clusters
        TBranch* bClsSum = NULL;
        bClsSum = tCls->GetBranch("AliFOCALCluster");
        TClonesArray* arrClsSum = NULL;
        bClsSum->SetAddress(&arrClsSum);
        bClsSum->GetEvent(0);
        Int_t nClsSum = arrClsSum->GetEntries();
        
        // canvas showing a 2d view of an event (Z-XY space)
        TCanvas* cMain = NULL;
        cMain = PrepareCanvas();
        DrawFOCAL(geometry,cMain);
        // plot the tracks of particles anticipating straight tracks in the direction of the initial momentum vector
        // (primary particles and optionally their direct daughters)
        DrawTracksMC(stack,cMain);
        // plot the positions of all segment-by-segment clusters
        DrawClusters(arrClsSeg,cMain,kTRUE);
        // plot the positions of all summed clusters
        DrawClusters(arrClsSum,cMain,kFALSE);
        // draw prefiltered clusters
        // (...)
        // save the canvas
        cout << " * ";
        cMain->SaveAs(Form("%s%sEv%03i.pdf",outDir.Data(),outSubDir.Data(),iEv+1));
        delete cMain;
    }
    delete runLoader;
    delete clusterizer;
    return;
}