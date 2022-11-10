// AnaMain.C
// David Grund, Nov 06, 2022

// root headers
#include "TSystem.h"
#include "TFileMerger.h"
#include "TFile.h"
#include "TList.h"
#include "TLorentzVector.h"
// my headers
#include "ConfigAnalysis.h"
#include "ConfigParameters.h"
#include "AnaMain.h"

TString outSubDir = "";

void MergeOutputFiles()
{
    TFileMerger m;
    m.OutputFile(Form("%smerged_%sanalysisResults.root",outDir.Data(),outSubDir.Data()));
    for(Int_t i = 0; i < nFiles; i++) 
    {
        TString sFile = Form("%s%03i/%sanalysisResults.root",outDir.Data(),i+1,outSubDir.Data());
        cout << "Adding file: " << sFile << endl;
        m.AddFile(Form("%s",sFile.Data()));
    }
    m.Merge();
    cout << endl << " FILES MERGED! " << endl;

    return;
}

void PrintHistograms(TString sFile)
{
    TFile* f = TFile::Open(sFile.Data(), "read");
    if(f) Printf("File %s loaded.", f->GetName());
    // get list with TH1F histograms
    TList *l1 = (TList*) f->Get("lTH1F");
    if(l1) Printf("List %s loaded.", l1->GetName()); 
    // get list with TH2F histograms
    TList *l2 = (TList*) f->Get("lTH2F");
    if(l2) Printf("List %s loaded.", l2->GetName()); 
    // get list with TProfile histograms
    TList *lP = (TList*) f->Get("lTPrf");
    if(lP) Printf("List %s loaded.", lP->GetName()); 
    // go over the lists and print all histograms
    TString sOut = Form("%smerged_%s",outDir.Data(),outSubDir.Data());
    // TH1F
    l1->ls();
	TIter next1(l1);
	TObject* object = NULL;
	while ((object = next1()))
	{
		cout << "Got an object " << object->GetName() << ":" << endl;
        DrawHisto(((TH1F*)object),sOut);
	}
    // TH2F
    l2->ls();
	TIter next2(l2);
	object = NULL;
	while ((object = next2()))
	{
		cout << "Got an object " << object->GetName() << ":" << endl;
        DrawHistoCOLZ(((TH2F*)object),sOut);
	}
    // TH1F
    lP->ls();
	TIter nextP(lP);
	object = NULL;
	while ((object = nextP()))
	{
		cout << "Got an object " << object->GetName() << ":" << endl;
        DrawHisto(((TProfile*)object),sOut);
	}

    return;
}

void AnalyzeClPairs(TString sDir, Bool_t debug = kFALSE)
{
    // access input file
    TFile* f = TFile::Open(sDir + "analysisResults.root", "read");
    if(f) Printf("File %s loaded.", f->GetName());
    TTree* tCls = dynamic_cast<TTree*>(f->Get("tCls"));
    if(tCls) Printf("Tree %s loaded.", tCls->GetName());
    SetBranchAddresses_tCls(tCls);

    // print first 100 tree entries
    if(debug)
    {
        for(Int_t i = 0; i < 100; i++)
        {
            tCls->GetEntry(i);
            cout << Form("Cl %i: ev %03i\n", i, fEvNumber);
        }
        cout << "\n";
    }

    // create output file containing a tree of cl pairs
    Float_t fEnClPair, fPtClPair, fEtaClPair, fPhiClPair;
    Float_t fEnJElPair, fPtJElPair, fEtaJElPair, fPhiJElPair;
    TFile* fClPairs = new TFile(sDir + "tClPairs.root","RECREATE");
    TTree* tClPairs = new TTree("tClPairs", "output tree containing cluster pairs");
    // pairs of summed clusters/superclusters
    tClPairs->Branch("fEnClPair", &fEnClPair, "fEnClPair/F");
    tClPairs->Branch("fPtClPair", &fPtClPair, "fPtClPair/F");
    tClPairs->Branch("fEtaClPair", &fEtaClPair, "fEtaClPair/F");
    tClPairs->Branch("fPhiClPair", &fPhiClPair, "fPhiClPair/F");
    // pairs of J/psi electrons (if cluster pairs was matched with it)
    tClPairs->Branch("fEnJElPair", &fEnJElPair, "fEnJElPair/F");
    tClPairs->Branch("fPtJElPair", &fPtJElPair, "fPtJElPair/F");
    tClPairs->Branch("fEtaJElPair", &fEtaJElPair, "fEtaJElPair/F");
    tClPairs->Branch("fPhiJElPair", &fPhiJElPair, "fPhiJElPair/F");
    gROOT->cd();

    // loop over entries in the tree
    Int_t iCl = 0;
    Int_t nCls = tCls->GetEntries();
    if(debug) cout << "Tree contains " << nCls << " cls" << endl;
    while(iCl < nCls-1)
    {
        // get the first cluster cluster
        tCls->GetEntry(iCl);
        Int_t EvNumberThisCl = fEvNumber;
        if(debug) cout << Form("   Cl %i: ev %03i\n", iCl, fEvNumber);
        // create its 4-vector
        TLorentzVector* thisCl = new TLorentzVector();
        thisCl->SetPtEtaPhiE(fPtCl,fEtaCl,fPhiCl,fEnCl);
        // create 4-vector of the matched physical primary electron
        TLorentzVector* thisEl = new TLorentzVector();
        thisEl->SetPtEtaPhiE(fPtJEl,fEtaJEl,fPhiJEl,fEnJEl);
        // create TLists where all clusters and pp electrons from the same event will be stored
        TList lCls;
        TList lEls;
        lCls.SetOwner(kTRUE);
        lEls.SetOwner(kTRUE);
        // using this command, the list is an owner of its content
        // content will be deleted whenever the list itself is deleted
        lCls.AddLast(thisCl);
        lEls.AddLast(thisEl);
        // create a vector to store indices of phys. prim. electrons matched with the clusters
        std::vector<Int_t> idxMtchJEl;
        idxMtchJEl.push_back(fIdxJEl);
        Int_t nClThisEvent(1.);
        // is next cluster from the same event?
        tCls->GetEntry(iCl+1);
        Int_t EvNumberNextCl = fEvNumber;
        while(EvNumberNextCl == EvNumberThisCl)
        {
            iCl++;
            nClThisEvent++;
            if(debug) cout << Form(" + Cl %i: ev %03i\n", iCl, fEvNumber);
            // create its 4-vector
            TLorentzVector* nextCl = new TLorentzVector();
            nextCl->SetPtEtaPhiE(fPtCl,fEtaCl,fPhiCl,fEnCl);
            // create 4-vector of the matched physical primary electron
            TLorentzVector* nextEl = new TLorentzVector();
            nextEl->SetPtEtaPhiE(fPtJEl,fEtaJEl,fPhiJEl,fEnJEl);
            // add them to the lists
            lCls.AddLast(nextCl);
            lEls.AddLast(nextEl);
            // store the electron index
            idxMtchJEl.push_back(fIdxJEl);
            // is next cluster from the same event?
            if(tCls->GetEntry(iCl+1) == 0) break;
            else EvNumberNextCl = fEvNumber;
        }
        // analyze the list of clusters from this event
        if(debug) cout << Form("-> Event %03i contains %i cls\n", EvNumberThisCl, nClThisEvent);
        // go over the cluster list 
        for(Int_t iCl1 = 0; iCl1 < nClThisEvent; iCl1++)
        {
            TLorentzVector* Cl1 = (TLorentzVector*)lCls.At(iCl1);
            if(!Cl1) continue;
            TLorentzVector* El1 = (TLorentzVector*)lEls.At(iCl1);
            if(!El1) continue;     
            // go over all possible pairs of clusters
            for(Int_t iCl2 = iCl1+1; iCl2 < nClThisEvent; iCl2++) 
            {
                TLorentzVector* Cl2 = (TLorentzVector*)lCls.At(iCl2);
                if(!Cl2) continue;
                TLorentzVector* El2 = (TLorentzVector*)lEls.At(iCl2);
                if(!El2) continue;
                // add the momentum vectors of the two clusters
                TLorentzVector Cl12 = *Cl1 + *Cl2;
                // cluster pair kinematics -> tree
                fEnClPair = Cl12.Energy();
                fPtClPair = Cl12.Pt();
                fEtaClPair = Cl12.Eta();
                fPhiClPair = Cl12.Phi();
                // add the momentum vector of the two electrons
                TLorentzVector El12 = *El1 + *El2;
                // electron pair kinematics -> tree
                fEnJElPair = El12.Energy();
                fPtJElPair = El12.Pt();
                fEtaJElPair = El12.Eta();
                fPhiJElPair = El12.Phi();
                // if both clusters are paired to a different physical primary electron
                if(idxMtchJEl[iCl1] != idxMtchJEl[iCl2]) {

                }
                tClPairs->Fill();
            }
        }
        // increase the iterator for the next cluster that will be loaded
        iCl++;
    }
    // save the file with the output tree
    fClPairs->Write("",TObject::kWriteDelete);
    delete fClPairs;

    /*
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
            tOutClPairs->Fill();
        } // end of for over iCl2
    } // end of for over iCl1
    */
    return;
}

void AnaMain(TString sim)
{
    ConfigLocalAnalysis(sim);
    outSubDir = CreateOutputSubDir();
    gSystem->Exec(Form("mkdir -p %smerged_%s",outDir.Data(),outSubDir.Data()));

    // merge all output files produced by FocalUpcGrid
    //MergeOutputFiles();

    // print all histograms from the merged output file
    //PrintHistograms(Form("%smerged_%sanalysisResults.root",outDir.Data(),outSubDir.Data()));

    // open merged output file, go over clusters and create cluster pairs
    // store cluster pairs in a tree
    //AnalyzeClPairs(Form("%s001/%s",outDir.Data(),outSubDir.Data()),kTRUE);
    AnalyzeClPairs(Form("%smerged_%s",outDir.Data(),outSubDir.Data()),kTRUE);
}