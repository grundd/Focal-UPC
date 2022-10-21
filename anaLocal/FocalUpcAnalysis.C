// FocalUpcAnalysisJpsi.C
// David Grund, Oct 19, 2022

#include <iostream>
#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()
// root headers
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TObjArray.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TROOT.h"
// aliroot headers
#include "AliRunLoader.h"
#include "AliFOCALClusterizerv2.h"
// my headers
#include "FocalUpcAnalysis_Utilities.h"
#include "FocalUpcAnalysis_DefineOutput.h"

// *****************************************************************************
// options to set:
Bool_t plotEvInfo = kFALSE;
Bool_t printEvInfo = kFALSE;
// superclusterizer:
Bool_t doSupercls = kFALSE; // if true, create superclusters using:
const Float_t minSeedE = 5; // [GeV]
const Float_t radius = 8; // [cm]
// selections:
const Float_t cutM = 0.0; // [GeV]; if > 0, filter out all (super)cluster pairs having mass below cutM
const Float_t cutE = 10.0; // [GeV]; if > 0, filter out all (super) clusters having energy below cutE 
// cuts used during matching with MC particles:
const Float_t cutdEta = 0.4; // [-] difference in eta of a MC particle and a (super)cluster
const Float_t cutdPhi = 0.4; // [-] difference in phi angles of a MC particle and a (super)cluster
// *****************************************************************************
// input files and configuration:
Int_t vSim = 2;
const Int_t nInputBoxEl[2] = {16,16};
const Int_t nInputBoxPh[2] = {6, 0};
const Int_t nInputJpsi[2] = {11,11};
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

TObjArray* arrTH1F = NULL;
TObjArray* arrTH2F = NULL;
TObjArray* arrTPrf = NULL;

void DoFocalAnalysis(TString dataset, Int_t pdgSim, Int_t nEv = -1);
// dataset -> name of a folder where the input files are stored
// nEv -> number of events from each input file to be analyzed 
// (if -1, then analyze all available events)

void FocalUpcAnalysis(Int_t pdgSim, Bool_t testOnly = kFALSE)
// pdgSim = 11  -> box simulations of electrons
//        = 22  -> box simulations of photons
//        = 443 -> kCohJpsiToElRad
// testOnly = false -> will run over all events in all available input files
//          = true  -> will run only over 100 events from the first file
{
    gSystem->Load("libpythia6_4_28.so");

    Int_t nInputFiles(0.);
    TString sInputFiles = "";
    // box electrons
    if(pdgSim == 11) {
        nInputFiles = nInputBoxEl[vSim-1];
        sInputFiles = "BoxElectrons_";
    }
    // box photons
    else if(pdgSim == 22) {
        nInputFiles = nInputBoxPh[vSim-1];
        sInputFiles = "BoxPhotons_";
    }
    // kCohJpsiToElRad
    else if(pdgSim == 443) {
        nInputFiles = nInputJpsi[vSim-1];
        sInputFiles = "kCohJpsiToElRad_";
    }
    else {
        cout << "Configuration not supported. Choose between pdgSim = 11, 22, 443. Terminating... " << endl;
        return;
    }

    // define output
    Int_t nTH1F(-1), nTPrf(-1), nTH2F(-1);
    // if box simulations
    if(pdgSim == 11 || pdgSim == 22) 
    {
        arrTH1F = new TObjArray(kBx_TH1F_all); DefineHisto_Bx_TH1F(arrTH1F); nTH1F = kBx_TH1F_all;
        arrTPrf = new TObjArray(kBx_TPrf_all); DefineHisto_Bx_TPrf(arrTPrf); nTPrf = kBx_TPrf_all;
        arrTH2F = new TObjArray(kBx_TH2F_all); DefineHisto_Bx_TH2F(arrTH2F); nTH2F = kBx_TH2F_all;
    } 
    // if starlight simulations
    else if(pdgSim == 443)
    {
        arrTH1F = new TObjArray(kJp_TH1F_all); DefineHisto_Jp_TH1F(arrTH1F); nTH1F = kJp_TH1F_all;
        arrTH2F = new TObjArray(kJp_TH2F_all); DefineHisto_Jp_TH2F(arrTH2F); nTPrf = kJp_TPrf_all;
        arrTPrf = new TObjArray(kJp_TPrf_all); DefineHisto_Jp_TPrf(arrTPrf); nTH2F = kJp_TH2F_all;
    }

    // mass filtering only for J/psi simulations
    if(cutM > 0 && (pdgSim == 11 || pdgSim == 22)) {
        cout << "Cannot do mass cleaning for box simulations of electrons/photons. Terminating... " << endl;
        return;
    }

    // run the analysis over selected input data
    if(nInputFiles == 0) { 
        cout << "No input files to analyze. Terminating... " << endl;
        return;
    } 
    if(testOnly) DoFocalAnalysis(sInputFiles + "001_1000ev/",pdgSim,100);
    else for(Int_t i = 0; i < nInputFiles; i++) DoFocalAnalysis(Form("%s%03i_1000ev/",sInputFiles.Data(),i+1),pdgSim);
    
    // post-analysis of histograms
    // if box simulations
    if(pdgSim == 11 || pdgSim == 22)
    {
        // integrate the histograms at kBx_totE_mcE and kBx_totEFromSegCls_mcE
        cout << "integral of hBx_totE_mcE: " << ((TH2F*)arrTH2F->At(kBx_totE_mcE))->Integral(1,nBinsE,1,nBinsE) << endl;
        cout << "integral of hBx_totEFromSegCls_mcE: " << ((TH2F*)arrTH2F->At(kBx_totEFromSegCls_mcE))->Integral(1,2*nBinsE,1,2*nBinsE) << endl;
    }
    // if starlight simulations
    else if(pdgSim == 443)
    {
        // calculate the fraction of coh J/psi -> ee events that have energy of both electron larger than 20 GeV
        Double_t nAll = ((TH2F*)arrTH2F->At(kJp_mcJE1E_mcJE2E))->Integral(1,nBinsE,1,nBinsE);
        Double_t nAbove20 = ((TH2F*)arrTH2F->At(kJp_mcJE1E_mcJE2E))->Integral(11,nBinsE,11,nBinsE);
        cout << "total no. of coh J/psi events: " << nAll << endl;
        cout << "no. of J/psi events with energy of both electron higher than 20 GeV: " << nAbove20 << endl;
        /*
        hMCJpsiEnRatio->Sumw2();
        hMCJpsiEnRatio->Divide(hMCJpsiEnAll);
        */
    }
    // draw output histograms
    TString outputFolder = Form("results/sim%02i_g%02i_p%02i/",vSim,vGeomFile,vParaFile);
    if(pdgSim == 11)       outputFolder += "_boxEle";
    else if(pdgSim == 22)  outputFolder += "_boxPho";
    else if(pdgSim == 443) outputFolder += "_Jpsi";
    if(doSupercls) outputFolder += Form("_supCl");
    if(cutM > 0)   outputFolder += Form("_cutM%.1f",cutM);
    if(cutE > 0)   outputFolder += Form("_cutE%.1f",cutE);
    outputFolder += "/";    
    gSystem->Exec("mkdir -p " + outputFolder);    
    for(Int_t i = 0; i < nTH1F; i++) if((TH1F*)arrTH1F->At(i)) DrawHisto<TH1F>(((TH1F*)arrTH1F->At(i)),outputFolder);
    for(Int_t i = 0; i < nTPrf; i++) if((TProfile*)arrTPrf->At(i)) DrawHisto<TProfile>(((TProfile*)arrTPrf->At(i)),outputFolder);
    for(Int_t i = 0; i < nTH2F; i++) if((TH2F*)arrTH2F->At(i)) DrawHistoCOLZ(((TH2F*)arrTH2F->At(i)),outputFolder);

    return;
}

void DoFocalAnalysis(TString dataset, Int_t pdgSim, Int_t nEv)
{
    if(plotEvInfo || printEvInfo) gSystem->Exec(Form("mkdir -p results/sim%02i_g%02i_p%02i/%seventInfo/",vSim,vGeomFile,vParaFile,dataset.Data()));

    // if clusters not yet produced, run the clusterizer
    TString sClFile = Form("inputData/sim%02i/%sfocalClusters_g%02i_p%02i.root",vSim,dataset.Data(),vGeomFile,vParaFile);
    if(gSystem->AccessPathName(sClFile.Data()))
    {
        cout << " MESSAGE: cluster file not found! Running the clusterizer now." << endl;
        // protection against accidentally overwriting the old cluster files:
        if(vGeomFile == 1) return;
        TString sCmd = Form(".x ClusterizeJpsi.C(\"inputData/sim%02i/%s\",%i,%i)",vSim,dataset.Data(),vGeomFile,vParaFile);
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
    TString geomFile = Form("geometry_%02i.txt",vGeomFile);
    clusterizer->InitGeometry(geomFile.Data());
    AliFOCALGeometry* geometry = clusterizer->GetGeometry();

    // define ALICE run loader: open galice.root
    AliRunLoader* runLoader = NULL;
    runLoader = AliRunLoader::Open(Form("inputData/sim%02i/%sgalice.root",vSim,dataset.Data()));
    if(!runLoader) 
    {
        cout << " ERROR: AliRunLoader not good! Terminating." << endl;
        return;   
    }
    if(!runLoader->GetAliRun()) runLoader->LoadgAlice();
    if(!runLoader->TreeE()) runLoader->LoadHeader();
    if(!runLoader->TreeK()) runLoader->LoadKinematics();

    Int_t nEvents(nEv);
    if(nEv == -1) nEvents = runLoader->GetNumberOfEvents();
    cout << "Analyzing " << nEvents << " events:" << endl;
    // loop over MC events contained within Kinematics.root
    Double_t progress = 0.; // perc
    for(Int_t iEv = 0; iEv < nEvents; iEv++) 
    {
        // update the progress bar
        if((iEv+1) % (Int_t)(nEvents/10.) == 0){
            progress += 10.;
            cout << "[" << progress << "%] done." << endl;
        }

        // get current MC event
        Int_t isEventOk = runLoader->GetEvent(iEv);
        // the method GetEvent(i) returns zero if event is loaded succesfully, if not we continue
        if (isEventOk != 0)
        {
            cout << "Event " << iEv << " not OK, skipping." << endl;
            continue;
        }

        // output text file
        TString sEvFile = Form("results/sim%02i_g%02i_p%02i/%seventInfo/Ev%i",vSim,vGeomFile,vParaFile,dataset.Data(),iEv);
        ofstream of;
        if(printEvInfo) {
            of.open((sEvFile + ".txt").Data());
            of << std::fixed << std::setprecision(0);
        }

        // get the stack of MC tracks for this event
        AliStack* stack = runLoader->Stack();

        // get a tree with clusters for this event 
        // (separate tree for each event in subfolders "Event0, Event1, ...")
        TTree *tCls = NULL;
        if(fCls->GetDirectory(Form("Event%i",iEv))) fCls->GetDirectory(Form("Event%i",iEv))->GetObject("fTreeR",tCls);
        else 
        {
            cout << " * MESSAGE: Cannot find event " << iEv << " in a cluster file " << fCls->GetName() << ". Skipping this event." << endl;
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
        if(printEvInfo) of << "seg. clusters: " << nClsSeg << endl;

        // find the total energy deposited in FoCal in this event from cls per segment
        Float_t ETotFromSegCls(0.);
        for(Int_t iCl = 0; iCl < nClsSeg; iCl++) 
        {
            AliFOCALCluster* clust = (AliFOCALCluster*) arrClsSeg->At(iCl);
            // count only the energy deposited in EMCal
            if(clust->Segment() < 6) ETotFromSegCls += clust->E();
        }

        // get the branch with final summed clusters
        TBranch* bClsSum = NULL;
        bClsSum = tCls->GetBranch("AliFOCALCluster");
        TClonesArray* arrClsSum = NULL;
        bClsSum->SetAddress(&arrClsSum);
        bClsSum->GetEvent(0);
        Int_t nClsSum = arrClsSum->GetEntries();
        if(printEvInfo) of << "sum. clusters: " << nClsSum << endl;

        // prepare a list of prefiltered clusters 
        // for the moment, add all clusters (later we will apply selections)
        TList listClsPref;
        listClsPref.SetOwner(kFALSE);
        // find the total and maximum cluster energy in this event
        Float_t ETot(0.), EClMax(-1.);
        Float_t EHCalTot(0.), EIsoR2Tot(0.), EIsoR4Tot(0.);
        Float_t ETotWithEHcalTot(0.), ETotWithEIsoR2Tot(0.), ETotWithEIsoR4Tot(0.);
        Int_t iMaxClE(-1); // index of a cluster with maximum energy
        for(Int_t iCl = 0; iCl < nClsSum; iCl++) 
        {
            AliFOCALCluster* clust = (AliFOCALCluster*) arrClsSum->At(iCl);
            listClsPref.Add(clust);
            // get the total energy summing over summed clusters
            ETot += clust->E();
            EHCalTot += clust->GetHCALEnergy();
            EIsoR2Tot += clust->GetIsoEnergyR2();
            EIsoR4Tot += clust->GetIsoEnergyR4();
            ETotWithEHcalTot += clust->E() + EHCalTot;
            ETotWithEIsoR2Tot += clust->E() + EIsoR2Tot;
            ETotWithEIsoR4Tot += clust->E() + EIsoR4Tot;
            // get the maximum cluster energy
            if(clust->E() > EClMax) {
                EClMax = clust->E();
                iMaxClE = iCl;
            } 
        }
        Int_t nClsPref = listClsPref.GetEntries();

        // prepare superclusters
        if(doSupercls)
        {
            // create new list: listClsPrefSort
            TList listClsPrefSort;
            listClsPrefSort.SetOwner(kFALSE);
            // fill it with the clusters from listClsPref, sorting them according to their energy
            while(listClsPrefSort.GetEntries() < nClsPref)
            {
                Float_t EClMaxFound = -1.;
                AliFOCALCluster* clustMaxE = NULL;
                for(Int_t iCl = 0; iCl < listClsPref.GetEntries(); iCl++) 
                {
                    AliFOCALCluster* clust = (AliFOCALCluster*)listClsPref.At(iCl);
                    Float_t ECl = clust->E();
                    if(ECl > EClMaxFound) {
                        EClMaxFound = ECl;
                        clustMaxE = clust;
                    }
                }
                listClsPrefSort.Add(clustMaxE);
                listClsPref.Remove(clustMaxE);
            }
            // while the energy of the most energetic cluster in listClsPrefSort is above minSeedE, 
            // take this cluster as a seed for a new supercluster
            // then combine it with all other clusters lying within the predefined radius
            AliFOCALCluster* clustTop;
            Float_t EClTop(-1.);
            if(nClsPref > 0) {
                clustTop = (AliFOCALCluster*)listClsPrefSort.At(0);
                EClTop = clustTop->E();
            }
            while(EClTop > minSeedE)
            {
                Float_t xSeed = clustTop->X();
                Float_t ySeed = clustTop->Y();
                Float_t xSupCl = clustTop->X() * EClTop;
                Float_t ySupCl = clustTop->Y() * EClTop;
                Float_t zSupCl = clustTop->Z() * EClTop;
                Float_t ESupCl = EClTop;
                listClsPrefSort.Remove(clustTop);
                for(Int_t iCl = listClsPrefSort.GetEntries()-1; iCl >= 0; iCl--) 
                {
                    AliFOCALCluster* clust = (AliFOCALCluster*)listClsPrefSort.At(iCl);
                    Float_t xCl = clust->X();
                    Float_t yCl = clust->Y();
                    Float_t zCl = clust->Z();
                    Float_t ECl = clust->E();
                    Float_t distXY = TMath::Sqrt(TMath::Power(xSeed-xCl,2) + TMath::Power(ySeed-yCl,2)); 
                    // if this cluster is close enough to the supercluster seed, add it to the supercluster
                    if(distXY < radius) {
                        xSupCl += xCl * ECl;
                        ySupCl += yCl * ECl;
                        zSupCl += zCl * ECl;
                        ESupCl += ECl;
                        listClsPrefSort.Remove(clust);
                    }
                }
                if(listClsPrefSort.GetEntries() > 0) {
                    clustTop = (AliFOCALCluster*)listClsPrefSort.At(0);
                    EClTop = clustTop->E();
                } 
                else EClTop = 0;
                // add the new supercluster to listClsPref
                AliFOCALCluster* supClNew = new AliFOCALCluster(xSupCl/ESupCl,ySupCl/ESupCl,zSupCl/ESupCl,ESupCl,-1);
                listClsPref.Add(supClNew);
            }
            nClsPref = listClsPref.GetEntries();
            if(printEvInfo) of << "identified superclusters: " << nClsPref << endl;
        }

        // if J/psi analysis, do "mass cleaning"
        // i.e., remove all clusters that have a mass below cutM when pared with any other cluster
        if(cutM > 0)
        {
            std::vector<Bool_t> foundLowMass;
            for(Int_t i = 0; i < nClsPref; i++) foundLowMass.push_back(kFALSE);
            // go over the clusters
            for(Int_t iCl1 = 0; iCl1 < nClsPref; iCl1++) 
            {            
                AliFOCALCluster* clust1 = (AliFOCALCluster*)listClsPref.At(iCl1);
                if(!clust1) continue;
                // get energy and coordinates of this cluster
                Float_t xCl1 = clust1->X();
                Float_t yCl1 = clust1->Y();
                Float_t zCl1 = clust1->Z();
                Float_t ECl1 = clust1->E();
                TLorentzVector cl1 = ConvertXYZEtoLorVec(xCl1,yCl1,zCl1,ECl1);
                // go over all possible pairs and combine clusters
                for(Int_t iCl2 = iCl1+1; iCl2 < nClsPref; iCl2++) 
                {
                    AliFOCALCluster* clust2 = (AliFOCALCluster*)listClsPref.At(iCl2);
                    if(!clust2) continue;
                    // get energy and coordinates of this cluster
                    Float_t xCl2 = clust2->X();
                    Float_t yCl2 = clust2->Y();
                    Float_t zCl2 = clust2->Z();
                    Float_t ECl2 = clust2->E();
                    TLorentzVector cl2 = ConvertXYZEtoLorVec(xCl2,yCl2,zCl2,ECl2);
                    // calculate the mass of the two current clusters
                    TLorentzVector cl12 = cl1 + cl2;
                    Double_t mass = cl12.M();
                    if (mass < cutM) 
                    {
                        foundLowMass[iCl1] = kTRUE;
                        foundLowMass[iCl2] = kTRUE;
                        // no need to look for any more pairs for cluster at iCl1
                        break;
                    }
                }
            }
            // remove all clusters giving low mass when paired
            for(Int_t iCl = nClsPref-1; iCl >= 0; iCl--)
            {
                AliFOCALCluster* clust = (AliFOCALCluster*)listClsPref.At(iCl);
                if(foundLowMass[iCl]) listClsPref.Remove(clust);
            }
            // how many clusters there are after prefiltering
            nClsPref = listClsPref.GetEntries();
            if(printEvInfo) of << "pref. (super)clusters after mass cleaning: " << nClsPref << endl;
        }

        // cut on minimum cluster energy
        // remove all remaining clusters (if any) with energy below threshold
        if(cutE > 0)
        {
            for(Int_t iCl = nClsPref-1; iCl >= 0; iCl--) 
            {
                AliFOCALCluster* clust = (AliFOCALCluster*)listClsPref.At(iCl);
                if(clust->E() < cutE) listClsPref.Remove(clust);
            }
            nClsPref = listClsPref.GetEntries();
            if(printEvInfo) of << "pref. (super)clusters with E > " << cutE << " GeV: " << nClsPref << endl;
        }

        // canvas showing a 2d view of an event (z-r space)
        TCanvas* cMain = NULL;
        if(plotEvInfo)
        {
            cMain = PrepareCanvas();
            DrawFOCAL(geometry,cMain);
            // plot the tracks of particles anticipating straight tracks in the direction of the initial momentum vector
            // plot primary particles and their direct daughters 
            DrawTracksMC(stack,cMain);
            // plot the positions of all segment-by-segment clusters
            DrawClusters(arrClsSeg,cMain,kTRUE);
            // plot the positions of all summed clusters
            DrawClusters(arrClsSum,cMain,kFALSE);
            // plot the positions of prefiltered (super)clusters
            DrawPrefClusters(&listClsPref,cMain);
            // save the canvas
            cout << " * ";
            cMain->SaveAs((sEvFile + ".pdf").Data());
        }
        delete cMain;

        // ******************************************************************************************************************
        // analysis of box electron or photon simulations
        // ******************************************************************************************************************
        if(pdgSim == 11 || pdgSim == 22)
        {
            // general histograms
            ((TH1F*)arrTH1F->At(kBx_mcE))->Fill(stack->Particle(0)->Energy());
            ((TH1F*)arrTH1F->At(kBx_mcPt))->Fill(stack->Particle(0)->Pt());
            ((TH2F*)arrTH2F->At(kBx_mcE_nCls))->Fill(stack->Particle(0)->Energy(), nClsPref);
            // MC energy vs total energies
            ((TH2F*)arrTH2F->At(kBx_mcE_totE))->Fill(stack->Particle(0)->Energy(), ETot);
            ((TH2F*)arrTH2F->At(kBx_totE_mcE))->Fill(ETot, stack->Particle(0)->Energy());
            ((TProfile*)arrTPrf->At(kBx_mcE_totE_prof))->Fill(stack->Particle(0)->Energy(), ETot);
            ((TProfile*)arrTPrf->At(kBx_totE_mcE_prof))->Fill(ETot, stack->Particle(0)->Energy());
            // total energy vs HCal energies
            ((TH2F*)arrTH2F->At(kBx_totEwHCalE_mcE))->Fill(ETotWithEHcalTot, stack->Particle(0)->Energy());
            ((TH2F*)arrTH2F->At(kBx_totEwIsoR2E_mcE))->Fill(ETotWithEIsoR2Tot, stack->Particle(0)->Energy());
            ((TH2F*)arrTH2F->At(kBx_totEwIsoR4E_mcE))->Fill(ETotWithEIsoR4Tot, stack->Particle(0)->Energy());
            ((TH2F*)arrTH2F->At(kBx_totE_totHCalE))->Fill(ETot, EHCalTot);
            ((TH2F*)arrTH2F->At(kBx_totE_totIsoR2E))->Fill(ETot, EIsoR2Tot);
            ((TH2F*)arrTH2F->At(kBx_totE_totIsoR4E))->Fill(ETot, EIsoR4Tot);
            ((TH2F*)arrTH2F->At(kBx_totEFromSegCls_mcE))->Fill(ETotFromSegCls*.001, stack->Particle(0)->Energy());
            // MC energy vs maximum cluster energy
            ((TH2F*)arrTH2F->At(kBx_mcE_maxClE))->Fill(stack->Particle(0)->Energy(), EClMax);
            ((TH2F*)arrTH2F->At(kBx_maxClE_mcE))->Fill(EClMax, stack->Particle(0)->Energy());
            ((TProfile*)arrTPrf->At(kBx_mcE_maxClE_prof))->Fill(stack->Particle(0)->Energy(), EClMax);
            ((TProfile*)arrTPrf->At(kBx_maxClE_mcE_prof))->Fill(EClMax, stack->Particle(0)->Energy());
            // correlation between MC energy and energies of prefiltered clusters
            for(Int_t iCl = 0; iCl < nClsPref; iCl++)
            {
                AliFOCALCluster *clust = (AliFOCALCluster*) listClsPref.At(iCl);
                ((TH2F*)arrTH2F->At(kBx_clE_mcE))->Fill(clust->E(), stack->Particle(0)->Energy());
            }
        }
        // ******************************************************************************************************************
        // analysis of J/psi simulations
        // ******************************************************************************************************************
        else 
        {
            // fill histograms with MC kinematics
            TObjArray ppElectrons; // array of physical primary electrons
            for(Int_t iTrk = 0; iTrk < stack->GetNtrack(); iTrk++)
            {
                TParticle *part = stack->Particle(iTrk);
                // if J/psi
                if(part->GetPdgCode() == 443) 
                {
                    ((TH1F*)arrTH1F->At(kJp_mcJPt))->Fill(part->Pt());
                    ((TH2F*)arrTH2F->At(kJp_mcJRap_mcJPt))->Fill(part->Y(), part->Pt());
                }
                // if physical primary (pp) electron (i.e., J/psi direct decay product)
                if(isEleOrPos(part) && stack->IsPhysicalPrimary(iTrk))
                {
                    TParticle* mother = stack->Particle(part->GetMother(0));
                    if(mother->GetPdgCode() != 443) cout << "Unexpected mother of pp electron, terminating." << endl;
                    ((TH2F*)arrTH2F->At(kJp_mcJEEta_mcJEPt))->Fill(part->Eta(), part->Pt());
                    ppElectrons.AddLast(part);
                }
            }
            if(ppElectrons.GetEntries() != 2) {
                cout << "Unexpected number of pp electrons: " << ppElectrons.GetEntries() << ", terminating..." << endl;
                return;
            } else {
                // fill the histogram showing the correlation between energies of the two pp electrons
                TParticle *ppe1 = (TParticle*) ppElectrons[0];
                TParticle *ppe2 = (TParticle*) ppElectrons[1];
                ((TH2F*)arrTH2F->At(kJp_mcJE1E_mcJE2E))->Fill(ppe1->Energy(), ppe2->Energy());
                ((TH2F*)arrTH2F->At(kJp_mcJE1Pt_mcJE2Pt))->Fill(ppe1->Pt(), ppe2->Pt());
            }

            // match each cluster to the closest (in XY distance) MC track, then go by mothers and find its primary physical particle
            // match each cluster to the closest primary physical MC track
            std::vector<Int_t> idxMtchP_prim; // vector of indices of physical primary particles which are mothers of particles to which clusters were matched
            TObjArray arrMtchP_prim; arrMtchP_prim.Clear(); // vector of the above primary particles
            std::vector<Int_t> idxMtchP_primDir; // vector of indices of physical primary particles to which the clusters were matched directly
            TObjArray arrMtchP_primDir; arrMtchP_primDir.Clear(); // vector of the above primary particles
            // initialize the arrays
            for(Int_t iCl = 0; iCl < nClsPref; iCl++)
            {
                idxMtchP_prim.push_back(-1);
                arrMtchP_prim.AddLast(NULL);
                idxMtchP_primDir.push_back(-1);
                arrMtchP_primDir.AddLast(NULL);
            }
            // now do the actual matching
            for(Int_t iCl = 0; iCl < nClsPref; iCl++)
            {
                AliFOCALCluster *clust = (AliFOCALCluster*) listClsPref.At(iCl);
                if(!clust) continue;
                // get energy and coordinates of this cluster
                Float_t xCl = clust->X();
                Float_t yCl = clust->Y();
                Float_t zCl = clust->Z();
                Float_t ECl = clust->E();
                TLorentzVector cl = ConvertXYZEtoLorVec(xCl,yCl,zCl,ECl);

                // fill some histograms with cluster kinematics
                ((TH2F*)arrTH2F->At(kJp_clEta_clPhi))->Fill(cl.Eta(), cl.Phi());
                ((TH2F*)arrTH2F->At(kJp_clEta_clPt ))->Fill(cl.Eta(), cl.Pt());

                // do the matching
                Float_t mtchXY(1e3);
                Float_t mtchXY_dir(1e3);
                TParticle* mtchP_any = NULL;
                TParticle* mtchP_prim = NULL;
                TParticle* mtchP_primDir = NULL;
                Int_t iMtchP_any = -1;
                Int_t iMtchP_prim = -1;
                Int_t iMtchP_primDir = -1;
                for(Int_t iTrk = 0; iTrk < stack->GetNtrack(); iTrk++)
                {
                    TParticle *part = stack->Particle(iTrk);
                    // during matching, skip photons with very low energy (< 0.2 MeV)
                    // trajectories of these often overlap with those of pp electron (when projected as straight lines)
                    // and clusters could very likely be matched with them instead of electrons
                    if(part->GetPdgCode() == 22 && part->Energy() < 0.2) continue;
                    // eta cut
                    Float_t dEta = TMath::Abs(part->Eta() - cl.Eta());
                    if(dEta > cutdEta) continue;
                    // phi cut
                    Float_t dPhi = TMath::Abs(part->Phi() - cl.Phi());
                    if(dPhi > TMath::Pi()) dPhi = 2*TMath::Pi() - dPhi;
                    if(dPhi > cutdPhi) continue;
                    // do matching
                    Float_t x(0.), y(0.);
                    TrackCoordinatesAtZ(part,zCl,x,y);
                    Float_t distXY = TMath::Sqrt(TMath::Power(x-xCl,2) + TMath::Power(y-yCl,2));
                    if(distXY < mtchXY) {
                        // a new closest MC particle found
                        mtchP_any = part;
                        iMtchP_any = iTrk;
                        mtchXY = distXY;
                    }
                    if((distXY < mtchXY_dir) && stack->IsPhysicalPrimary(iTrk)) {
                        // a new closest physical primary MC particle found
                        mtchP_primDir = part;
                        iMtchP_primDir = iTrk;
                        mtchXY_dir = distXY;
                    }
                }
                // if a matching particle found
                if(mtchP_any) {
                    // find the physical primary particle from which the matched particle originates
                    mtchP_prim = mtchP_any;
                    Bool_t physPrimFound = stack->IsPhysicalPrimary(iMtchP_any);
                    while(!physPrimFound) {
                        Int_t iNewMother = mtchP_prim->GetMother(0);
                        mtchP_prim = stack->Particle(iNewMother);
                        physPrimFound = stack->IsPhysicalPrimary(iNewMother);
                        iMtchP_any = iNewMother;
                    }
                    iMtchP_prim = iMtchP_any;
                }
                idxMtchP_prim[iCl] = iMtchP_prim;
                arrMtchP_prim.AddAt(mtchP_prim, iCl);
                idxMtchP_primDir[iCl] = iMtchP_primDir;
                arrMtchP_primDir.AddAt(mtchP_primDir, iCl);
            }
            // fill the histograms
            for(Int_t iCl = 0; iCl < nClsPref; iCl++)
            {
                AliFOCALCluster *clust = (AliFOCALCluster*) listClsPref.At(iCl);
                if(!clust) continue;
                // get energy and coordinates of this cluster
                Float_t xCl = clust->X();
                Float_t yCl = clust->Y();
                Float_t zCl = clust->Z();
                Float_t ECl = clust->E();
                TLorentzVector cl = ConvertXYZEtoLorVec(xCl,yCl,zCl,ECl);
                TParticle *mtchP_prim = (TParticle*)arrMtchP_prim[iCl];
                TParticle *mtchP_primDir = (TParticle*)arrMtchP_primDir[iCl];
                // if matched to a ppp going by mothers
                if(mtchP_prim) {
                    ((TH2F*)arrTH2F->At(kJp_clE_mtchE))->Fill(ECl, mtchP_prim->Energy());
                    ((TH2F*)arrTH2F->At(kJp_clEta_mtchEta))->Fill(cl.Eta(), mtchP_prim->Eta());
                    ((TH2F*)arrTH2F->At(kJp_clPt_mtchPt))->Fill(cl.Pt(), mtchP_prim->Pt());
                    // if electron (ppe)
                    if(isEleOrPos(mtchP_prim)) {
                        ((TH2F*)arrTH2F->At(kJp_primElClE_mtchE))->Fill(ECl, mtchP_prim->Energy());
                        ((TH2F*)arrTH2F->At(kJp_primElClEta_mtchEta))->Fill(cl.Eta(), mtchP_prim->Eta());
                        ((TH2F*)arrTH2F->At(kJp_primElClPt_mtchPt))->Fill(cl.Pt(), mtchP_prim->Pt());
                    }
                }
                // if matched to a ppp directly
                if(mtchP_primDir) {
                    ((TH2F*)arrTH2F->At(kJp_clE_mtchDirE))->Fill(ECl, mtchP_primDir->Energy());
                    ((TH2F*)arrTH2F->At(kJp_clEta_mtchDirEta))->Fill(cl.Eta(), mtchP_primDir->Eta());
                    ((TH2F*)arrTH2F->At(kJp_clPt_mtchDirPt))->Fill(cl.Pt(), mtchP_primDir->Pt());
                    // if electron (ppe)
                    if(isEleOrPos(mtchP_primDir)) {
                        ((TH2F*)arrTH2F->At(kJp_primElClE_mtchDirE))->Fill(ECl, mtchP_primDir->Energy());
                        ((TH2F*)arrTH2F->At(kJp_primElClEta_mtchDirEta))->Fill(cl.Eta(), mtchP_primDir->Eta());
                        ((TH2F*)arrTH2F->At(kJp_primElClPt_mtchDirPt))->Fill(cl.Pt(), mtchP_primDir->Pt());
                    }
                }
                // if matched to a ppp in both ways
                if(mtchP_prim && mtchP_primDir) 
                {
                    ((TH2F*)arrTH2F->At(kJp_mtchE_mtchDirE))->Fill(mtchP_prim->Energy(), mtchP_primDir->Energy());
                    // if both matched to ppe
                    if(isEleOrPos(mtchP_prim) && isEleOrPos(mtchP_primDir)) {
                        ((TH2F*)arrTH2F->At(kJp_mtchE_mtchDirE_primEl))->Fill(mtchP_prim->Energy(), mtchP_primDir->Energy());
                        ((TH2F*)arrTH2F->At(kJp_mtchEta_mtchDirEta_primEl))->Fill(mtchP_prim->Eta(), mtchP_primDir->Eta());
                        ((TH2F*)arrTH2F->At(kJp_mtchPt_mtchDirPt_primEl))->Fill(mtchP_prim->Pt(), mtchP_primDir->Pt());
                    }
                } 
            }
            // combine clusters into pairs
            for(Int_t iCl1 = 0; iCl1 < nClsPref; iCl1++) 
            {
                AliFOCALCluster *clust1 = (AliFOCALCluster*) listClsPref.At(iCl1);
                if(!clust1) continue;
                // get energy and coordinates of this cluster
                Float_t xCl1 = clust1->X();
                Float_t yCl1 = clust1->Y();
                Float_t zCl1 = clust1->Z();
                Float_t ECl1 = clust1->E();
                TLorentzVector cl1 = ConvertXYZEtoLorVec(xCl1,yCl1,zCl1,ECl1);
                // go over all possible pairs of clusters
                for(Int_t iCl2 = iCl1+1; iCl2 < nClsPref; iCl2++) 
                {
                    AliFOCALCluster *clust2 = (AliFOCALCluster*) listClsPref.At(iCl2);
                    if(!clust2) continue;
                    // get energy and coordinates of this cluster
                    Float_t xCl2 = clust2->X();
                    Float_t yCl2 = clust2->Y();
                    Float_t zCl2 = clust2->Z();
                    Float_t ECl2 = clust2->E();
                    TLorentzVector cl2 = ConvertXYZEtoLorVec(xCl2,yCl2,zCl2,ECl2);

                    // add the momentum vectors of the two clusters to get the total momentum
                    TLorentzVector cl12 = cl1 + cl2;
                    Double_t mass = cl12.M();
                    // fill some kinematic histograms
                    ((TH1F*)arrTH1F->At(kJp_clPairM))->Fill(mass);
                    ((TH1F*)arrTH1F->At(kJp_clPairEta))->Fill(cl12.Eta());
                    ((TH1F*)arrTH1F->At(kJp_clPairPt))->Fill(cl12.Pt());
                    if(mass > 2.5) ((TH1F*)arrTH1F->At(kJp_clPairPt_massCut))->Fill(cl12.Pt());
                    // if both paired to a different pp electron:
                    // matching by mothers
                    if((idxMtchP_prim[iCl1] != idxMtchP_prim[iCl2]) && idxMtchP_prim[iCl1] != -1 && idxMtchP_prim[iCl2] != -1)
                    {
                        ((TH1F*)arrTH1F->At(kJp_primElClPairM))->Fill(mass);
                    }
                    // direct matching
                    if((idxMtchP_primDir[iCl1] != idxMtchP_primDir[iCl2]) && idxMtchP_primDir[iCl1] != -1 && idxMtchP_primDir[iCl2] != -1)
                    {
                        ((TH1F*)arrTH1F->At(kJp_primElClPairM_dir))->Fill(mass);
                    }
                } // end of for over iCl2
            } // end of for over iCl1
        }
        of.close();

        /*
        // go over prefiltered clusters and match them with MC particles
        for(Int_t iCl1 = 0; iCl1 < listClsPref.GetEntries(); iCl1++) 
        {
            // fill the histograms
            if(isMatchingParticleJpsiEle1) 
            {
                hMismtchXY->Fill(matchedE1,matchedXY1); 
                hMCJpsiEnMtch->Fill(matchedE1);
                hMCJpsiEnRatio->Fill(matchedE1);
            }

            // go over all possible pairs of prefiltered clusters and combine them
            for(Int_t iCl2 = iCl1+1; iCl2 < listClsPref.GetEntries(); iCl2++) 
            {
                // fill the histograms
                if(isMatchingParticleJpsiEle1 && isMatchingParticleJpsiEle2)
                {
                    hClPairJpsiEleEta->Fill(cl12.Eta());
                    hClPairJpsiElePt->Fill(cl12.Pt());
                    hClPairJpsiEleVsMCJpsiEta->Fill(cl12.Eta(),stack->Particle(0)->Eta());
                    hClPairJpsiEleVsMCJpsiPt->Fill(cl12.Pt(),stack->Particle(0)->Pt());
                }
            }  
        }  
        */ 
    }  // end of for over events in AliRunLoader

    runLoader->Delete();
    fCls->Close();

    delete clusterizer;

    return;
}