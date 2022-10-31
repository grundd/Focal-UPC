// FocalUpc_AnalysisSecondary.C
// David Grund, Oct 31, 2022

// root headers
// ...
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

    return;
}