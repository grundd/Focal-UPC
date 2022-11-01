// FocalUpc_AnalysisSecondary.C
// David Grund, Oct 31, 2022

// root headers
#include "TLorentzVector.h"
// my headers
#include "FocalUpc_Utilities.h"
#include "FocalUpc_ConfigAnalysis.h"

void MergeTrees()
{
    gSystem->Exec("mkdir -p " + sOut + "merged/");

    TList* lMergedTreesCl = new TList;
    TList* lMergedTreesClPairs = new TList;
    // go over the dataset trees and add them to the list
    for(Int_t i = 0; i < nInputFiles; i++)
    {
        TString sTrees = Form("%s%s%03i_1000ev/_analysisTrees.root", sOut.Data(), sInputFiles.Data(), i+1);
        TFile* fTrees = TFile::Open(sTrees.Data(), "read");
        if(fTrees) Printf("File %s loaded.", fTrees->GetName());
        // tree with clusters
        TTree* tCls = dynamic_cast<TTree*>(fTrees->Get("tOutCls"));
        if(tCls) Printf("Tree %s loaded.", tCls->GetName());
        lMergedTreesCl->Add(tCls);
        // tree with cluster pairs
        TTree* tClPairs = dynamic_cast<TTree*>(fTrees->Get("tOutClPairs"));
        if(tClPairs) Printf("Tree %s loaded.", tClPairs->GetName());
        lMergedTreesClPairs->Add(tClPairs);
    }
    // merge the individual trees into one
    TFile* fMergedTrees = new TFile((sOut + "merged/_mergedTrees.root").Data(),"RECREATE");
    TTree* tMergedCls = TTree::MergeTrees(lMergedTreesCl);
    tMergedCls->SetName("tMergedCls");
    TTree* tMergedClPairs = TTree::MergeTrees(lMergedTreesClPairs);
    tMergedClPairs->SetName("tMergedClPairs");
    fMergedTrees->Write("",TObject::kWriteDelete);
    fMergedTrees->Close();

    return;
}

void FocalUpc_AnalysisSecondary(TString sSim)
{
    if(!ConfigAnalysis(sSim)) {
        cout << "Wrong configuration. Terminating..." << endl;
        return;
    }

    // merge the trees from the individual datasets into a single tree that will be analyzed
    MergeTrees();

    // load the merged tree
    TFile* f = TFile::Open((sOut + "merged/_mergedTrees.root").Data(), "read");
    if(f) Printf("File %s loaded.", f->GetName());
    // tree with clusters
    TTree* tCls = dynamic_cast<TTree*>(f->Get("tMergedCls"));
    if(tCls) Printf("Tree %s loaded.", tCls->GetName());
    SetBranchAddresses_tCls(tCls);
    // tree with cluster pairs
    TTree* tClPairs = dynamic_cast<TTree*>(f->Get("tMergedClPairs"));
    if(tClPairs) Printf("Tree %s loaded.", tClPairs->GetName());
    SetBranchAddresses_tClPairs(tClPairs);

    // go over events in the trees
    // tree with clusters
    Int_t nEntries = tCls->GetEntries();
    Printf("%i entries found in %s.", nEntries, tCls->GetName());
    Float_t progress = 0.; // perc

    TH2F* hJ2_ppeClEta_mtchEta = new TH2F("hJ2_ppeClEta_mtchEta","",100,2.,7.,100,2.,7.);
    hJ2_ppeClEta_mtchEta->SetTitle("#eta of a cluster matched with a pp electron vs #eta of ppe;#eta_{cl} [-];#eta_{matched ppe^{#pm}} [-]");

    for(Int_t iEntry = 0; iEntry < nEntries; iEntry++)
    {
        tCls->GetEntry(iEntry);
        // update the progress bar
        if((iEntry+1) % (Int_t)(nEntries/10.) == 0) {
            progress += 10.;
            cout << "[" << progress << "%] done." << endl;
        }
        if(!areSame(fEtaJEl,-1e3)) hJ2_ppeClEta_mtchEta->Fill(fEtaCl,fEtaJEl);
    }

    // tree with cluster pairs
    nEntries = tClPairs->GetEntries();
    Printf("%i entries found in %s.", nEntries, tClPairs->GetName());
    progress = 0.; // perc

    TH2F* hJ2_ppeClPairRap_mtchRap = new TH2F("hJ2_ppeClPairRap_mtchRap","",100,2.,7.,100,2.,7.);
    hJ2_ppeClPairRap_mtchRap->SetTitle("#it{y} of a cluster pair matched with a pair of pp electrons vs #it{y} of ppe pair;#it{y}_{cl. pair} [-];#it{y}_{matched ppe pair} [-]");
    TH2F* hJ2_ppeClPairEn_mtchEn = new TH2F("hJ2_ppeClPairEn_mtchEn","",100,0.,200.,100,0.,200.);
    hJ2_ppeClPairEn_mtchEn->SetTitle("energy of a cluster pair matched with a pair of pp electrons vs energy of ppe pair;#it{E}_{cl. pair} [GeV];#it{E}_{matched ppe pair} [GeV]");
    
    for(Int_t iEntry = 0; iEntry < nEntries; iEntry++)
    {
        tClPairs->GetEntry(iEntry);
        // update the progress bar
        if((iEntry+1) % (Int_t)(nEntries/10.) == 0) {
            progress += 10.;
            cout << "[" << progress << "%] done." << endl;
        }
        TLorentzVector lvClPair;
        lvClPair.SetPtEtaPhiE(fPtClPair,fEtaClPair,fPhiClPair,fEnClPair);
        TLorentzVector lvJElPair;
        if(!areSame(fEtaJElPair,-1e3)) {
            lvJElPair.SetPtEtaPhiE(fPtJElPair,fEtaJElPair,fPhiJElPair,fEnJElPair);
            hJ2_ppeClPairRap_mtchRap->Fill(lvClPair.Rapidity(),lvJElPair.Rapidity());
            hJ2_ppeClPairEn_mtchEn->Fill(fEnClPair,fEnJElPair);
        } 
    }

    DrawHistoCOLZ(hJ2_ppeClEta_mtchEta,(sOut + "_").Data());
    DrawHistoCOLZ(hJ2_ppeClPairRap_mtchRap,(sOut + "_").Data());
    DrawHistoCOLZ(hJ2_ppeClPairEn_mtchEn,(sOut + "_").Data());

    return;
}