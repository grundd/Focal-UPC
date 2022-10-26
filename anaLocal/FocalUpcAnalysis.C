// FocalUpcAnalysisJpsi.C
// David Grund, Oct 19, 2022

#include <iostream>
#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()
// root headers
#include "TROOT.h"
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
Bool_t matchDirectly = kFALSE;
// superclusterizer:
Bool_t doSupercls = kFALSE; // if true, create superclusters using:
const Float_t minSeedE = 5; // [GeV]
const Float_t radius = 15; // [cm]
// selections:
const Float_t cutM = 0.0; // [GeV]; if > 0, filter out all (super)cluster pairs having mass below cutM
const Float_t cutE = 0.0; // [GeV]; if > 0, filter out all (super) clusters having energy below cutE 
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
TString sOut = "";

void DoFocalAnalysis(TString sDataset, Int_t pdgSim, Int_t nEv = -1);
// sDataset -> name of a folder where the input files are stored
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
        cout << "Configuration not supported. Choose between pdgSim = 11, 22, 443. Terminating..." << endl;
        return;
    }

    // define output
    Int_t nTH1F(-1), nTPrf(-1), nTH2F(-1);
    // if box simulations
    if(pdgSim == 11 || pdgSim == 22) 
    {
        arrTH1F = new TObjArray(kB1_all); DefineHisto_BxH1(arrTH1F); nTH1F = kB1_all;
        arrTH2F = new TObjArray(kB2_all); DefineHisto_BxH2(arrTH2F); nTH2F = kB2_all;
        arrTPrf = new TObjArray(kBP_all); DefineHisto_BxPr(arrTPrf); nTPrf = kBP_all;
    } 
    // if starlight simulations
    else if(pdgSim == 443)
    {
        arrTH1F = new TObjArray(kJ1_all); DefineHisto_JpH1(arrTH1F); nTH1F = kJ1_all;
        arrTH2F = new TObjArray(kJ2_all); DefineHisto_JpH2(arrTH2F); nTH2F = kJ2_all; 
        arrTPrf = new TObjArray(kpJp_all); DefineHisto_JpPr(arrTPrf); nTPrf = kpJp_all;
    }

    // mass filtering only for J/psi simulations
    if(cutM > 0 && (pdgSim == 11 || pdgSim == 22)) {
        cout << "Cannot do mass cleaning for box simulations of electrons/photons. Terminating... " << endl;
        return;
    }

    // create output folder
    sOut = Form("results/sim%02i_g%02i_p%02i/",vSim,vGeomFile,vParaFile);
    if(pdgSim == 11)       sOut += "_boxEle";
    else if(pdgSim == 22)  sOut += "_boxPho";
    else if(pdgSim == 443) sOut += "_Jpsi";
    if(doSupercls) sOut += Form("_supCl");
    if(pdgSim == 443 && matchDirectly) sOut += Form("_dirMtch");
    if(cutM > 0)   sOut += Form("_cutM%.1f",cutM);
    if(cutE > 0)   sOut += Form("_cutE%.1f",cutE);
    sOut += "/";    
    gSystem->Exec("mkdir -p " + sOut);

    // run the analysis over selected input data
    if(nInputFiles == 0) { 
        cout << "No input files to analyze. Terminating... " << endl;
        return;
    } 
    if(testOnly) DoFocalAnalysis(sInputFiles + "001_1000ev/",pdgSim,100);
    else for(Int_t i = 0; i < nInputFiles; i++) DoFocalAnalysis(Form("%s%03i_1000ev/",sInputFiles.Data(),i+1),pdgSim);

    // post-analysis of histograms
    ofstream of((sOut + "log.txt").Data());
    // if box simulations
    if(pdgSim == 11 || pdgSim == 22)
    {
        // ...
    }
    // if starlight simulations
    else if(pdgSim == 443)
    {
        // calculate the fraction of coh J/psi -> ee events that have energy of both electron larger than 20 GeV
        Float_t nAll = ((TH2F*)arrTH2F->At(kJ2_mcJEl1En_mcJEl2En))->Integral(1,nBinsEn,1,nBinsEn);
        Float_t nAbove20 = ((TH2F*)arrTH2F->At(kJ2_mcJEl1En_mcJEl2En))->Integral(11,nBinsEn,11,nBinsEn);
        of << "total no. of coh J/psi events: " << nAll << endl;
        of << "no. of J/psi events with energy of both electrons higher than 20 GeV: " << nAbove20 << endl;
    }
    of.close();

    // draw output histograms
    for(Int_t i = 0; i < nTH1F; i++) if((TH1F*)arrTH1F->At(i)) DrawHisto<TH1F>(((TH1F*)arrTH1F->At(i)),sOut);
    for(Int_t i = 0; i < nTH2F; i++) if((TH2F*)arrTH2F->At(i)) DrawHistoCOLZ(((TH2F*)arrTH2F->At(i)),sOut);
    for(Int_t i = 0; i < nTPrf; i++) if((TProfile*)arrTPrf->At(i)) DrawHisto<TProfile>(((TProfile*)arrTPrf->At(i)),sOut);

    return;
}

void DoFocalAnalysis(TString sDataset, Int_t pdgSim, Int_t nEv)
{
    gSystem->Exec(Form("mkdir -p %s%s",sOut.Data(),sDataset.Data()));

    // if clusters not yet produced, run the clusterizer
    TString sClFile = Form("inputData/sim%02i/%sfocalClusters_g%02i_p%02i.root",vSim,sDataset.Data(),vGeomFile,vParaFile);
    if(gSystem->AccessPathName(sClFile.Data()))
    {
        cout << " MESSAGE: cluster file not found! Running the clusterizer now." << endl;
        // protection against accidentally overwriting old cluster files:
        if(vGeomFile == 1) return;
        TString sCmd = Form(".x ClusterizeJpsi.C(\"inputData/sim%02i/%s\",%i,%i)",vSim,sDataset.Data(),vGeomFile,vParaFile);
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
    runLoader = AliRunLoader::Open(Form("inputData/sim%02i/%sgalice.root",vSim,sDataset.Data()));
    if(!runLoader) 
    {
        cout << " ERROR: AliRunLoader not good! Terminating." << endl;
        return;   
    }
    if(!runLoader->GetAliRun()) runLoader->LoadgAlice();
    if(!runLoader->TreeE()) runLoader->LoadHeader();
    if(!runLoader->TreeK()) runLoader->LoadKinematics();

    // define output tree
    TString sTree = Form("%s%s_analysisTree.root", sOut.Data(), sDataset.Data());
    TFile* fOut = new TFile(sTree.Data(),"RECREATE");
    TTree* tOut = new TTree("tOut", "output tree containing cluster pairs");
    // pairs of summed clusters/superclusters
    Float_t fEnClPair, fPtClPair, fEtaClPair, fPhiClPair;
    tOut->Branch("fEnClPair", &fEnClPair, "fEnClPair/F");
    tOut->Branch("fPtClPair", &fPtClPair, "fPtClPair/F");
    tOut->Branch("fEtaClPair", &fEtaClPair, "fEtaClPair/F");
    tOut->Branch("fPhiClPair", &fPhiClPair, "fPhiClPair/F");
    // pairs of J/psi electrons (if cluster pairs was matched with it)
    Float_t fEnJElPair, fPtJElPair, fEtaJElPair, fPhiJElPair;
    tOut->Branch("fEnJElPair", &fEnJElPair, "fEnJElPair/F");
    tOut->Branch("fPtJElPair", &fPtJElPair, "fPtJElPair/F");
    tOut->Branch("fEtaJElPair", &fEtaJElPair, "fEtaJElPair/F");
    tOut->Branch("fPhiJElPair", &fPhiJElPair, "fPhiJElPair/F");
    gROOT->cd();

    Int_t nEvents(nEv);
    if(nEv == -1) nEvents = runLoader->GetNumberOfEvents();
    cout << "Analyzing " << nEvents << " events:" << endl;
    // loop over MC events contained within Kinematics.root
    Float_t progress = 0.; // perc
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
        TString sEvInfo = Form("%s%sEv%i",sOut.Data(),sDataset.Data(),iEv);
        ofstream of;
        if(printEvInfo) {
            of.open((sEvInfo + ".txt").Data());
            of << std::fixed << std::setprecision(0);
        }

        // get the stack of MC tracks for this event
        AliStack* stack = runLoader->Stack();

        // get a tree with clusters for this event 
        // (separate tree for each event in subfolders "Event0, Event1, ...")
        TTree* tCls = NULL;
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
        // for the moment, add all clusters (later we will apply the selections)
        TList* listClsPref = new TList();
        listClsPref->SetOwner(kFALSE);
        // find the total and maximum cluster energy in this event
        Float_t ETot(0.), EClMax(-1.);
        Float_t EHCalTot(0.), EIsoR2Tot(0.), EIsoR4Tot(0.);
        Float_t ETotWithEHcalTot(0.), ETotWithEIsoR2Tot(0.), ETotWithEIsoR4Tot(0.);
        Int_t iClMaxE(-1); // index of a cluster with maximum energy
        for(Int_t iCl = 0; iCl < nClsSum; iCl++) 
        {
            AliFOCALCluster* clust = (AliFOCALCluster*) arrClsSum->At(iCl);
            listClsPref->Add(clust);
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
                iClMaxE = iCl;
            } 
        }
        Int_t nClsPref = listClsPref->GetEntries();

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
                for(Int_t iCl = 0; iCl < listClsPref->GetEntries(); iCl++) 
                {
                    AliFOCALCluster* clust = (AliFOCALCluster*)listClsPref->At(iCl);
                    Float_t ECl = clust->E();
                    if(ECl > EClMaxFound) {
                        EClMaxFound = ECl;
                        clustMaxE = clust;
                    }
                }
                listClsPrefSort.Add(clustMaxE);
                listClsPref->Remove(clustMaxE);
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
                listClsPref->Add(supClNew);
            }
            nClsPref = listClsPref->GetEntries();
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
                AliFOCALCluster* clust1 = (AliFOCALCluster*)listClsPref->At(iCl1);
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
                    AliFOCALCluster* clust2 = (AliFOCALCluster*)listClsPref->At(iCl2);
                    if(!clust2) continue;
                    // get energy and coordinates of this cluster
                    Float_t xCl2 = clust2->X();
                    Float_t yCl2 = clust2->Y();
                    Float_t zCl2 = clust2->Z();
                    Float_t ECl2 = clust2->E();
                    TLorentzVector cl2 = ConvertXYZEtoLorVec(xCl2,yCl2,zCl2,ECl2);
                    // calculate the mass of the two current clusters
                    TLorentzVector cl12 = cl1 + cl2;
                    Float_t mass = cl12.M();
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
                AliFOCALCluster* clust = (AliFOCALCluster*)listClsPref->At(iCl);
                if(foundLowMass[iCl]) listClsPref->Remove(clust);
            }
            // how many clusters there are after prefiltering
            nClsPref = listClsPref->GetEntries();
            if(printEvInfo) of << "pref. (super)clusters after mass cleaning: " << nClsPref << endl;
        }

        // cut on minimum cluster energy
        // remove all remaining clusters (if any) with energy below threshold
        if(cutE > 0)
        {
            for(Int_t iCl = nClsPref-1; iCl >= 0; iCl--) 
            {
                AliFOCALCluster* clust = (AliFOCALCluster*)listClsPref->At(iCl);
                if(clust->E() < cutE) listClsPref->Remove(clust);
            }
            nClsPref = listClsPref->GetEntries();
            if(printEvInfo) of << "pref. (super)clusters with E > " << cutE << " GeV: " << nClsPref << endl;
        }

        // canvas showing a 2d view of an event (Z-XY space)
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
            DrawPrefClusters(listClsPref,cMain);
            // save the canvas
            cout << " * ";
            cMain->SaveAs((sEvInfo + ".pdf").Data());
        }
        delete cMain;

        // ******************************************************************************************************************
        // analysis of box electron or photon simulations
        // ******************************************************************************************************************
        if(pdgSim == 11 || pdgSim == 22)
        {
            // general histograms
            ((TH1F*)arrTH1F->At(kB1_mcE))->Fill(stack->Particle(0)->Energy());
            ((TH1F*)arrTH1F->At(kB1_mcPt))->Fill(stack->Particle(0)->Pt());
            ((TH2F*)arrTH2F->At(kB2_mcE_nCls))->Fill(stack->Particle(0)->Energy(), nClsPref);
            // total energy vs MC energy
            ((TH2F*)arrTH2F->At(kB2_totE_mcE))->Fill(ETot, stack->Particle(0)->Energy());
            ((TProfile*)arrTPrf->At(kBP_totE_mcE))->Fill(ETot, stack->Particle(0)->Energy());
            // total energy vs HCal energies
            ((TH2F*)arrTH2F->At(kB2_totEwHCalE_mcE))->Fill(ETotWithEHcalTot, stack->Particle(0)->Energy());
            ((TH2F*)arrTH2F->At(kB2_totEwIsoR2E_mcE))->Fill(ETotWithEIsoR2Tot, stack->Particle(0)->Energy());
            ((TH2F*)arrTH2F->At(kB2_totEwIsoR4E_mcE))->Fill(ETotWithEIsoR4Tot, stack->Particle(0)->Energy());
            ((TH2F*)arrTH2F->At(kB2_totE_totHCalE))->Fill(ETot, EHCalTot);
            ((TH2F*)arrTH2F->At(kB2_totE_totIsoR2E))->Fill(ETot, EIsoR2Tot);
            ((TH2F*)arrTH2F->At(kB2_totE_totIsoR4E))->Fill(ETot, EIsoR4Tot);
            ((TH2F*)arrTH2F->At(kB2_totEFromSegCls_mcE))->Fill(ETotFromSegCls*.001, stack->Particle(0)->Energy());
            // maximum cluster energy vs MC energy
            ((TH2F*)arrTH2F->At(kB2_maxClE_mcE))->Fill(EClMax, stack->Particle(0)->Energy());
            ((TProfile*)arrTPrf->At(kBP_maxClE_mcE))->Fill(EClMax, stack->Particle(0)->Energy());
            // X and Y of the cluster with the highest energy
            Float_t xClMaxE(-1.), yClMaxE(-1.);
            Bool_t isClWithMaxE = kFALSE;
            if(iClMaxE != -1 && cutM == 0 && cutE == 0 && !doSupercls) {
                isClWithMaxE = kTRUE;
                AliFOCALCluster* clustMaxE = (AliFOCALCluster*) listClsPref->At(iClMaxE);
                xClMaxE = clustMaxE->X();
                yClMaxE = clustMaxE->Y();
            }
            // info per prefiltered (super)clusters
            for(Int_t iCl = 0; iCl < nClsPref; iCl++)
            {
                AliFOCALCluster* clust = (AliFOCALCluster*) listClsPref->At(iCl);
                Float_t xCl = clust->X();
                Float_t yCl = clust->Y();
                Float_t zCl = clust->Z();
                // physical primary electron:
                TParticle* part = stack->Particle(0);
                Float_t x(0.), y(0.);
                TrackCoordinatesAtZ(part,zCl,x,y);
                // correlation between MC energy and energies of prefiltered clusters
                ((TH2F*)arrTH2F->At(kB2_clE_mcE))->Fill(clust->E(), part->Energy());
                // radial distance between the cluster and the MC track
                Float_t DeltaX = xCl - x;
                Float_t DeltaY = yCl - y;
                Float_t DeltaR = TMath::Sqrt(TMath::Power(DeltaX,2) + TMath::Power(DeltaY,2));
                ((TH2F*)arrTH2F->At(kB2_clMcDX_clMcDY))->Fill(DeltaX, DeltaY);
                ((TH2F*)arrTH2F->At(kB2_mcE_clMcSep))->Fill(part->Energy(), DeltaR);
                // radial distance between the cluster and the cluster with maximum energy
                if(isClWithMaxE && iCl != iClMaxE) {
                    Float_t DeltaX_maxE = xCl - xClMaxE;
                    Float_t DeltaY_maxE = yCl - yClMaxE;
                    Float_t DeltaR_maxE = TMath::Sqrt(TMath::Power(DeltaX_maxE,2) + TMath::Power(DeltaY_maxE,2));
                    ((TH2F*)arrTH2F->At(kB2_clMaxClDX_clMaxClDY))->Fill(DeltaX_maxE, DeltaY_maxE);
                    ((TH2F*)arrTH2F->At(kB2_mcE_clMaxClSep))->Fill(part->Energy(), DeltaR_maxE);
                }
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
                TParticle* part = stack->Particle(iTrk);
                // if J/psi
                if(part->GetPdgCode() == 443) {
                    ((TH1F*)arrTH1F->At(kJ1_mcJEn))->Fill(part->Energy());
                    ((TH1F*)arrTH1F->At(kJ1_mcJPt))->Fill(part->Pt());
                    ((TH1F*)arrTH1F->At(kJ1_mcJRap))->Fill(part->Y());
                    ((TH1F*)arrTH1F->At(kJ1_mcJM))->Fill(part->GetCalcMass());
                    ((TH2F*)arrTH2F->At(kJ2_mcJRap_mcJPt))->Fill(part->Y(), part->Pt());
                }
                // if physical primary (pp) electron (i.e., J/psi direct decay product)
                if(isEleOrPos(part) && stack->IsPhysicalPrimary(iTrk)) {
                    TParticle* mother = stack->Particle(part->GetMother(0));
                    if(mother->GetPdgCode() != 443) cout << " * MESSAGE: Unexpected mother of a pp electron." << endl;
                    ((TH2F*)arrTH2F->At(kJ2_mcJElEta_mcJElPt))->Fill(part->Eta(), part->Pt());
                    ppElectrons.AddLast(part);
                }
            }
            // kinematics of the two physical primary electrons
            TParticle* ppEl1 = NULL;
            TParticle* ppEl2 = NULL;
            TLorentzVector* lvJElPair = new TLorentzVector();
            if(ppElectrons.GetEntries() != 2) {
                cout << " * ERROR: Unexpected number of pp electrons: " << ppElectrons.GetEntries() << ", terminating..." << endl;
                runLoader->Delete();
                fCls->Close();
                return;
            } else {
                // fill the histogram showing the correlation between energies and transverse momenta of the two pp electrons
                ppEl1 = (TParticle*) ppElectrons[0];
                ppEl2 = (TParticle*) ppElectrons[1];
                ((TH2F*)arrTH2F->At(kJ2_mcJEl1En_mcJEl2En))->Fill(ppEl1->Energy(), ppEl2->Energy());
                ((TH2F*)arrTH2F->At(kJ2_mcJEl1Pt_mcJEl2Pt))->Fill(ppEl1->Pt(), ppEl2->Pt());
                // fill the histograms showing kinematics of J/psi reconstructed from the two pp electrons
                TLorentzVector lvPPEl1;
                lvPPEl1.SetPxPyPzE(ppEl1->Px(),ppEl1->Py(),ppEl1->Pz(),ppEl1->Energy());
                TLorentzVector lvPPEl2;
                lvPPEl2.SetPxPyPzE(ppEl2->Px(),ppEl2->Py(),ppEl2->Pz(),ppEl2->Energy());
                TLorentzVector Jpsi = lvPPEl1 + lvPPEl2;
                lvJElPair->SetPxPyPzE(Jpsi.Px(),Jpsi.Py(),Jpsi.Pz(),Jpsi.Energy());
                ((TH1F*)arrTH1F->At(kJ1_mcJElPairEn))->Fill(lvJElPair->Energy());
                ((TH1F*)arrTH1F->At(kJ1_mcJElPairPt))->Fill(lvJElPair->Pt());
                ((TH1F*)arrTH1F->At(kJ1_mcJElPairRap))->Fill(lvJElPair->Rapidity());
                ((TH1F*)arrTH1F->At(kJ1_mcJElPairM))->Fill(lvJElPair->M());
            }

            // match clusters to MC tracks
            std::vector<Int_t> idxMtchPhysPrimParts; 
            TObjArray arrMtchPhysPrimParts;
            MatchClsToPhysPrimP(stack,listClsPref,idxMtchPhysPrimParts,&arrMtchPhysPrimParts,matchDirectly);

            // fill the histograms
            for(Int_t iCl = 0; iCl < nClsPref; iCl++)
            {
                AliFOCALCluster* clust = (AliFOCALCluster*) listClsPref->At(iCl);
                if(!clust) continue;
                // get energy and coordinates of this cluster
                Float_t xCl = clust->X();
                Float_t yCl = clust->Y();
                Float_t zCl = clust->Z();
                Float_t ECl = clust->E();
                TLorentzVector cl = ConvertXYZEtoLorVec(xCl,yCl,zCl,ECl);

                ((TH1F*)arrTH1F->At(kJ1_clZ))->Fill(zCl);
                // fill some histograms with cluster kinematics
                ((TH2F*)arrTH2F->At(kJ2_clEta_clPhi))->Fill(cl.Eta(), cl.Phi());
                ((TH2F*)arrTH2F->At(kJ2_clEta_clPt ))->Fill(cl.Eta(), cl.Pt());

                TParticle* mtchPhysPrimPart = (TParticle*)arrMtchPhysPrimParts[iCl];
                // if matched to a ppp going by mothers
                if(mtchPhysPrimPart) {
                    ((TH2F*)arrTH2F->At(kJ2_pppClEn_mtchEn))->Fill(ECl, mtchPhysPrimPart->Energy());
                    ((TH2F*)arrTH2F->At(kJ2_pppClEta_mtchEta))->Fill(cl.Eta(), mtchPhysPrimPart->Eta());
                    ((TH2F*)arrTH2F->At(kJ2_pppClPt_mtchPt))->Fill(cl.Pt(), mtchPhysPrimPart->Pt());
                    // if electron (ppe)
                    if(isEleOrPos(mtchPhysPrimPart)) {
                        ((TH2F*)arrTH2F->At(kJ2_ppeClEn_mtchEn))->Fill(ECl, mtchPhysPrimPart->Energy());
                        ((TH2F*)arrTH2F->At(kJ2_ppeClEta_mtchEta))->Fill(cl.Eta(), mtchPhysPrimPart->Eta());
                        ((TH2F*)arrTH2F->At(kJ2_ppeClPt_mtchPt))->Fill(cl.Pt(), mtchPhysPrimPart->Pt());
                    }
                }
            }
            // combine clusters into pairs
            for(Int_t iCl1 = 0; iCl1 < nClsPref; iCl1++) 
            {
                AliFOCALCluster* clust1 = (AliFOCALCluster*) listClsPref->At(iCl1);
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
                    AliFOCALCluster* clust2 = (AliFOCALCluster*) listClsPref->At(iCl2);
                    if(!clust2) continue;
                    // get energy and coordinates of this cluster
                    Float_t xCl2 = clust2->X();
                    Float_t yCl2 = clust2->Y();
                    Float_t zCl2 = clust2->Z();
                    Float_t ECl2 = clust2->E();
                    TLorentzVector cl2 = ConvertXYZEtoLorVec(xCl2,yCl2,zCl2,ECl2);

                    // add the momentum vectors of the two clusters to get the total momentum
                    TLorentzVector cl12 = cl1 + cl2;
                    // cluster pair kinematics -> tree
                    fEnClPair = cl12.Energy();
                    fPtClPair = cl12.Pt();
                    fEtaClPair = cl12.Eta();
                    fPhiClPair = cl12.Phi();
                    Float_t mass = cl12.M();
                    // pp electron pair kinematics -> tree
                    fEnJElPair = -1e3;
                    fPtJElPair = -1e3;
                    fEtaJElPair = -1e3;
                    fPhiJElPair = -1e3;
                    // fill some kinematic histograms
                    ((TH1F*)arrTH1F->At(kJ1_clPairEn))->Fill(fEnClPair);
                    ((TH1F*)arrTH1F->At(kJ1_clPairPt))->Fill(fPtClPair);
                    if(mass > 2.8) ((TH1F*)arrTH1F->At(kJ1_clPairPt_massCut))->Fill(fPtClPair);
                    ((TH1F*)arrTH1F->At(kJ1_clPairRap))->Fill(cl12.Rapidity());
                    ((TH1F*)arrTH1F->At(kJ1_clPairM))->Fill(mass);
                    // radial separation between the pairs of clusters
                    Float_t sepCl = TMath::Sqrt(TMath::Power(xCl1-xCl2,2) + TMath::Power(yCl1-yCl2,2));
                    ((TH1F*)arrTH1F->At(kJ1_clPairSep))->Fill(sepCl);
                    // radial separation between cluster pair vs between the pair of ppe
                    Float_t x1(0.), y1(0.);
                    TrackCoordinatesAtZ(ppEl1,zCl1,x1,y1);
                    Float_t x2(0.), y2(0.);
                    TrackCoordinatesAtZ(ppEl2,zCl2,x2,y2);
                    Float_t sepMC = TMath::Sqrt(TMath::Power(x1-x2,2) + TMath::Power(y1-y2,2));
                    ((TH2F*)arrTH2F->At(kJ2_clPairSep_mcJElSep))->Fill(sepCl,sepMC);
                    // if both paired to a different pp electron:
                    if((idxMtchPhysPrimParts[iCl1] != idxMtchPhysPrimParts[iCl2]) && idxMtchPhysPrimParts[iCl1] != -1 && idxMtchPhysPrimParts[iCl2] != -1)
                    {
                        fEnJElPair = lvJElPair->Energy();
                        fPtJElPair = lvJElPair->Pt();
                        fEtaJElPair = lvJElPair->Eta();
                        fPhiJElPair = lvJElPair->Phi();
                        ((TH1F*)arrTH1F->At(kJ1_ppeClPairM))->Fill(mass);
                        ((TH1F*)arrTH1F->At(kJ1_ppeClPairSep))->Fill(sepCl); 
                        ((TH2F*)arrTH2F->At(kJ2_ppeClPairEn_mtchEn))->Fill(fEnClPair,fEnJElPair);
                        ((TH2F*)arrTH2F->At(kJ2_ppeClPairRap_mtchRap))->Fill(cl12.Rapidity(),lvJElPair->Rapidity());
                        ((TH2F*)arrTH2F->At(kJ2_ppeClPairPt_mtchPt))->Fill(fPtClPair,fPtJElPair);
                        ((TH2F*)arrTH2F->At(kJ2_ppeClPairM_mtchM))->Fill(mass,lvJElPair->M());
                    } else {
                        ((TH1F*)arrTH1F->At(kJ1_sameppeClPairSep))->Fill(sepCl); 
                    }
                    tOut->Fill();
                } // end of for over iCl2
            } // end of for over iCl1
            delete lvJElPair;
        }
        of.close();
        delete listClsPref;
    }  // end of for over events in AliRunLoader

    runLoader->Delete();
    fCls->Close();
    delete clusterizer;

    fOut->Write("",TObject::kWriteDelete);
    delete fOut;

    return;
}