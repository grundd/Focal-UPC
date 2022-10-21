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
Bool_t makePlots = kFALSE;
Bool_t massFiltering = kTRUE;
Int_t vSim = 2;
// version of the simulated input files:
// 1 -> files produced on Oct 17, 2022
// 2 -> files produced on Oct 19, 2022, after an update of FoCal&AliRoot
Int_t vGeomFile = 2;
// version of the FoCal software and the "geometry.txt" file using which the clusters were produced 
// 1 -> "geometry_01.txt" (used before the update on Oct 19, 2022)
// 2 -> "geometry_02.txt" (used after the update on Oct 19, 2022)
Int_t vParaFile = 2;
// version of the "parameters.txt" file using which the clusters were produced
// 1 -> "parameters_01.txt": Calib.Const1 set to 12900. (current official file)
// 2 -> "parameters_02.txt": Calib.Const1 set to 11220. (old value)
// values of the cuts:
const Float_t cutE = 0.0;  // [GeV] cut on minimal cluster energy
const Float_t cutM = 0.2;  // [GeV] cut on mass of cluster pair
const Float_t cutdEta = 0.4; // [-] cut on the difference in eta for primary MC particle and prefiltered cluster
const Float_t cutdPhi = 0.4; // [-] cut on the difference in phi angles for primary MC particle and prefiltered cluster
// input box simulations of electrons
const Int_t nInputBoxEl = 16;
TString sInputBoxEl[nInputBoxEl] = {
    "BoxElectrons_001_1000ev/",
    "BoxElectrons_002_1000ev/",
    "BoxElectrons_003_1000ev/",
    "BoxElectrons_004_1000ev/",
    "BoxElectrons_005_1000ev/",
    "BoxElectrons_006_1000ev/",
    "BoxElectrons_007_1000ev/",
    "BoxElectrons_008_1000ev/",
    "BoxElectrons_009_1000ev/",
    "BoxElectrons_010_1000ev/",
    "BoxElectrons_011_1000ev/",
    "BoxElectrons_012_1000ev/",
    "BoxElectrons_013_1000ev/",
    "BoxElectrons_014_1000ev/",
    "BoxElectrons_015_1000ev/",
    "BoxElectrons_016_1000ev/"
};
// input box simulations of photons
const Int_t nInputBoxPh = 6;
TString sInputBoxPh[nInputBoxPh] = {
    "BoxPhotons_001_1000ev/",
    "BoxPhotons_002_1000ev/",
    "BoxPhotons_003_1000ev/",
    "BoxPhotons_004_1000ev/",
    "BoxPhotons_005_1000ev/",
    "BoxPhotons_006_1000ev/"
};
// input starlight simulations of kCohJpsiToElRad
const Int_t nInputJpsi = 11;
TString sInputJpsi[nInputJpsi] = {
    "kCohJpsiToElRad_001_1000ev/",
    "kCohJpsiToElRad_002_1000ev/",
    "kCohJpsiToElRad_003_1000ev/",
    "kCohJpsiToElRad_004_1000ev/",
    "kCohJpsiToElRad_005_1000ev/",
    "kCohJpsiToElRad_006_1000ev/",
    "kCohJpsiToElRad_007_1000ev/",
    "kCohJpsiToElRad_008_1000ev/",
    "kCohJpsiToElRad_009_1000ev/",
    "kCohJpsiToElRad_010_1000ev/",
    "kCohJpsiToElRad_011_1000ev/"
};
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

    // define output
    if(pdgSim == 11 || pdgSim == 22) {
        arrTH1F = new TObjArray(kBx_TH1F_all); DefineHisto_Bx_TH1F(arrTH1F);
        arrTH2F = new TObjArray(kBx_TH2F_all); DefineHisto_Bx_TH2F(arrTH2F);
        arrTPrf = new TObjArray(kBx_TPrf_all); DefineHisto_Bx_TPrf(arrTPrf);
    } else if(pdgSim == 443) {
        arrTH1F = new TObjArray(kJp_TH1F_all); DefineHisto_Jp_TH1F(arrTH1F);
        arrTH2F = new TObjArray(kJp_TH2F_all); DefineHisto_Jp_TH2F(arrTH2F);
        arrTPrf = new TObjArray(kJp_TPrf_all); DefineHisto_Jp_TPrf(arrTPrf);
    } else return;
    // run the analysis over selected input data
    if(pdgSim == 11) 
    {
        if(testOnly) DoFocalAnalysis(sInputBoxEl[0],11,10);
        else for(Int_t i = 0; i < nInputBoxEl; i++) DoFocalAnalysis(sInputBoxEl[i],11);
    } 
    else if(pdgSim == 22) 
    {
        if(testOnly) DoFocalAnalysis(sInputBoxPh[0],22,10);
        else for(Int_t i = 0; i < nInputBoxPh; i++) DoFocalAnalysis(sInputBoxPh[i],22);
    } 
    else 
    {
        if(testOnly) DoFocalAnalysis(sInputJpsi[0],443,10);
        else for(Int_t i = 0; i < nInputJpsi; i++) DoFocalAnalysis(sInputJpsi[i],443);
    }
    // post-analysis of histograms
    if(pdgSim == 11 || pdgSim == 22)
    {
        // integrate the histograms at kBx_totE_mcE and kBx_totEFromSegCls_mcE
        cout << "integral of hBx_totE_mcE: " << ((TH2F*)arrTH2F->At(kBx_totE_mcE))->Integral(1,nBinsE,1,nBinsE) << endl;
        cout << "integral of hBx_totEFromSegCls_mcE: " << ((TH2F*)arrTH2F->At(kBx_totEFromSegCls_mcE))->Integral(1,2*nBinsE,1,2*nBinsE) << endl;
    }
    else 
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
    Int_t nTH1F(-1), nTPrf(-1), nTH2F(-1);
    TString outputFolder = Form("results/sim%02i_g%02i_p%02i/",vSim,vGeomFile,vParaFile);
    if(pdgSim == 11) {
        nTH1F = kBx_TH1F_all; nTPrf = kBx_TPrf_all; nTH2F = kBx_TH2F_all;
        outputFolder += "_boxEle";
    } else if(pdgSim == 22) {
        nTH1F = kBx_TH1F_all; nTPrf = kBx_TPrf_all; nTH2F = kBx_TH2F_all;
        outputFolder += "_boxPho";
    } else {
        nTH1F = kJp_TH1F_all; nTPrf = kJp_TPrf_all; nTH2F = kJp_TH2F_all;
        outputFolder += "_Jpsi";
    }
    if(pdgSim == 443) outputFolder+= Form("_cutM%.1f",cutM);
    outputFolder += Form("_cutE%.1f/",cutE);
    gSystem->Exec("mkdir -p " + outputFolder);    
    for(Int_t i = 0; i < nTH1F; i++) if((TH1F*)arrTH1F->At(i)) DrawHisto<TH1F>(((TH1F*)arrTH1F->At(i)),outputFolder);
    for(Int_t i = 0; i < nTPrf; i++) if((TProfile*)arrTPrf->At(i)) DrawHisto<TProfile>(((TProfile*)arrTPrf->At(i)),outputFolder);
    for(Int_t i = 0; i < nTH2F; i++) if((TH2F*)arrTH2F->At(i)) DrawHistoCOLZ(((TH2F*)arrTH2F->At(i)),outputFolder);

    return;
}

void DoFocalAnalysis(TString dataset, Int_t pdgSim, Int_t nEv)
{
    gSystem->Exec(Form("mkdir -p results/sim%02i_g%02i_p%02i/%seventInfo/",vSim,vGeomFile,vParaFile,dataset.Data()));

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
        ofstream of((sEvFile + ".txt").Data());
        of << std::fixed << std::setprecision(0);

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
        of << "seg. clusters: " << nClsSeg << endl;

        // find the total energy deposited in FoCal in this event from cls per segment
        Float_t totEFromSegCls(0.);
        for(Int_t iCl = 0; iCl < nClsSeg; iCl++) 
        {
            AliFOCALCluster* clust = (AliFOCALCluster*) arrClsSeg->At(iCl);
            // count only the energy deposited in EMCal
            //if(clust->Segment() < 6) 
            totEFromSegCls += clust->E();
        }

        // get the branch with final summed clusters
        TBranch* bClsSum = NULL;
        bClsSum = tCls->GetBranch("AliFOCALCluster");
        TClonesArray* arrClsSum = NULL;
        bClsSum->SetAddress(&arrClsSum);
        bClsSum->GetEvent(0);
        Int_t nClsSum = arrClsSum->GetEntries();
        of << "sum. clusters: " << nClsSum << endl;

        // prepare a list of prefiltered clusters 
        // for the moment, add all clusters (later we will apply mass cleaning and cut on minimum energy)
        TList listClsPref;
        listClsPref.SetOwner(kFALSE);
        // find the total and maximum cluster energy in this event
        Float_t totE(0.), totEwHCalE(0.), totEwIsoER2(0.), totEwIsoER4(0.), maxClE(-1.);
        Int_t iMaxClE(-1); // index of a cluster with maximum energy
        for(Int_t iCl = 0; iCl < nClsSum; iCl++) 
        {
            AliFOCALCluster* clust = (AliFOCALCluster*) arrClsSum->At(iCl);
            listClsPref.Add(clust);
            // get the total energy summing over summed clusters
            totE += clust->E();
            totEwHCalE += clust->E() + clust->GetHCALEnergy();
            totEwIsoER2 += clust->E() + clust->GetIsoEnergyR2();
            totEwIsoER4 += clust->E() + clust->GetIsoEnergyR4();
            // get the maximum cluster energy
            if(clust->E() > maxClE) {
                maxClE = clust->E();
                iMaxClE = iCl;
            } 
        }
        Int_t nClsPref = listClsPref.GetEntries();

        // if J/psi analysis, do "mass cleaning"
        // i.e., remove all clusters that have a mass below cutM when pared with any other cluster
        if(pdgSim == 443 && massFiltering)
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
            of << "pref. clusters after mass cleaning: " << nClsPref << endl;
        }

        // cut on minimum cluster energy
        // remove all remaining clusters (if any) with energy below threshold
        for(Int_t iCl = 0; iCl < nClsSum; iCl++) 
        {
            AliFOCALCluster* clust = (AliFOCALCluster*) arrClsSum->At(iCl);
            if(clust->E() < cutE) listClsPref.Remove(clust);
        }
        nClsPref = listClsPref.GetEntries();
        of << "pref. clusters with E > " << cutE << " GeV: " << nClsPref << endl;

        // canvas showing a 2d view of an event (z-r space)
        TCanvas* cMain = NULL;
        if(makePlots)
        {
            cMain = PrepareCanvas();
            DrawFOCAL(geometry,cMain);
            // plot the tracks of particles anticipating straight tracks in the direction of the initial momentum vector
            // plot primary particles and their direct daughters 
            DrawTracksMC(stack,cMain);
            // plot the positions of segment-by-segment clusters
            DrawSegClusters(arrClsSeg,cMain);
            // plot the positions of prefiltered summed clusters
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
            // MC energy vs total energy over all clusters
            ((TH2F*)arrTH2F->At(kBx_mcE_totE))->Fill(stack->Particle(0)->Energy(), totE);
            ((TH2F*)arrTH2F->At(kBx_totE_mcE))->Fill(totE, stack->Particle(0)->Energy());
            ((TH2F*)arrTH2F->At(kBx_totEwHCalE_mcE))->Fill(totEwHCalE, stack->Particle(0)->Energy());
            ((TH2F*)arrTH2F->At(kBx_totEwIsoER2_mcE))->Fill(totEwIsoER2, stack->Particle(0)->Energy());
            ((TH2F*)arrTH2F->At(kBx_totEwIsoER4_mcE))->Fill(totEwIsoER4, stack->Particle(0)->Energy());
            ((TH2F*)arrTH2F->At(kBx_totEFromSegCls_mcE))->Fill(totEFromSegCls*.001, stack->Particle(0)->Energy());
            ((TProfile*)arrTPrf->At(kBx_mcE_totE_prof))->Fill(stack->Particle(0)->Energy(), totE);
            ((TProfile*)arrTPrf->At(kBx_totE_mcE_prof))->Fill(totE, stack->Particle(0)->Energy());
            // MC energy vs maximum cluster energy
            ((TH2F*)arrTH2F->At(kBx_mcE_maxClE))->Fill(stack->Particle(0)->Energy(), maxClE);
            ((TH2F*)arrTH2F->At(kBx_maxClE_mcE))->Fill(maxClE, stack->Particle(0)->Energy());
            ((TProfile*)arrTPrf->At(kBx_mcE_maxClE_prof))->Fill(stack->Particle(0)->Energy(), maxClE);
            ((TProfile*)arrTPrf->At(kBx_maxClE_mcE_prof))->Fill(maxClE, stack->Particle(0)->Energy());
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

            // find superclusters
            // ... to do

            // calibrate the two groups of clusters?
            // ... to do

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
                    if(mass > 2.2) ((TH1F*)arrTH1F->At(kJp_clPairPt_massCut))->Fill(cl12.Pt());
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
                hRelDiffVsMtchEn->Fill(matchedE1,(matchedE1-e1)/matchedE1);
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