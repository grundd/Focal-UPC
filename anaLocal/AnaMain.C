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
#include "CreateHistograms.h"

TString outSubDir = "";

void MergeOutputFiles()
{
    TFileMerger m;
    m.OutputFile(Form("%smerged_%sanalysisResultsGrid.root",outDir.Data(),outSubDir.Data()));
    for(Int_t i = 0; i < nFiles; i++) 
    {
        TString sFile = Form("%s%03i/%sanalysisResultsGrid.root",outDir.Data(),i+1,outSubDir.Data());
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

void PrepareClPairs(TString sDir, Bool_t debug = kFALSE)
{
    // access input file
    TFile* f = TFile::Open(sDir + "analysisResultsGrid.root", "read");
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

    // array with output histograms
    TObjArray* arrHistos = NULL;
    arrHistos = new TObjArray(kMainJpsi_all);
    CreateHistos_MainJpsi(arrHistos);

    // create output file containing a tree of cl pairs
    Float_t fEnClPair, fPtClPair, fEtaClPair, fPhiClPair;
    Float_t fEnJElPair, fPtJElPair, fEtaJElPair, fPhiJElPair;
    TFile* fOut = new TFile(sDir + "analysisResultsMain.root","RECREATE");
    TTree* tOut = new TTree("tClPairs", "output tree containing cluster pairs");
    // pairs of summed clusters/superclusters
    tOut->Branch("fEnClPair", &fEnClPair, "fEnClPair/F");
    tOut->Branch("fPtClPair", &fPtClPair, "fPtClPair/F");
    tOut->Branch("fEtaClPair", &fEtaClPair, "fEtaClPair/F");
    tOut->Branch("fPhiClPair", &fPhiClPair, "fPhiClPair/F");
    // pairs of J/psi electrons (if cluster pairs was matched with it)
    tOut->Branch("fEnJElPair", &fEnJElPair, "fEnJElPair/F");
    tOut->Branch("fPtJElPair", &fPtJElPair, "fPtJElPair/F");
    tOut->Branch("fEtaJElPair", &fEtaJElPair, "fEtaJElPair/F");
    tOut->Branch("fPhiJElPair", &fPhiJElPair, "fPhiJElPair/F");
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
                Float_t clPairM = Cl12.M();
                // add the momentum vector of the two electrons
                TLorentzVector El12 = *El1 + *El2;
                // electron pair kinematics -> tree
                fEnJElPair = El12.Energy();
                fPtJElPair = El12.Pt();
                fEtaJElPair = El12.Eta();
                fPhiJElPair = El12.Phi();
                // fill some kinematic histograms
                ((TH1F*)arrHistos->At(kJ1_clPairEn))->Fill(fEnClPair);
                ((TH1F*)arrHistos->At(kJ1_clPairPt))->Fill(fPtClPair);
                if(clPairM > 2.8) ((TH1F*)arrHistos->At(kJ1_clPairPt_massCut))->Fill(fPtClPair);
                ((TH1F*)arrHistos->At(kJ1_clPairRap))->Fill(Cl12.Rapidity());
                ((TH1F*)arrHistos->At(kJ1_clPairM))->Fill(clPairM);
                // radial separation between the pairs of clusters
                /*
                Float_t sepCl = TMath::Sqrt(TMath::Power(xCl1-xCl2,2) + TMath::Power(yCl1-yCl2,2));
                ((TH1F*)arrTH1F->At(kJ1_clPairSep))->Fill(sepCl);
                // radial separation between cluster pair vs between the pair of ppe
                Float_t x1(0.), y1(0.);
                TrackCoordinatesAtZ(ppEl1,zCl1,x1,y1);
                Float_t x2(0.), y2(0.);
                TrackCoordinatesAtZ(ppEl2,zCl2,x2,y2);
                Float_t sepMC = TMath::Sqrt(TMath::Power(x1-x2,2) + TMath::Power(y1-y2,2));
                ((TH2F*)arrTH2F->At(kJ2_clPairSep_mcJElSep))->Fill(sepCl,sepMC);
                */
                // if both clusters are paired to a different physical primary electron
                if(idxMtchJEl[iCl1] != idxMtchJEl[iCl2]) {
                    ((TH1F*)arrHistos->At(kJ1_ppeClPairM))->Fill(clPairM);
                    //((TH1F*)arrHistos->At(kJ1_ppeClPairSep))->Fill(sepCl); 
                    ((TH2F*)arrHistos->At(kJ2_ppeClPairEn_mtchEn))->Fill(fEnClPair,fEnJElPair);
                    ((TH2F*)arrHistos->At(kJ2_ppeClPairRap_mtchRap))->Fill(Cl12.Rapidity(),El12.Rapidity());
                    ((TH2F*)arrHistos->At(kJ2_ppeClPairPt_mtchPt))->Fill(fPtClPair,fPtJElPair);
                    ((TH2F*)arrHistos->At(kJ2_ppeClPairM_mtchM))->Fill(clPairM,El12.M());
                } else {
                    //((TH1F*)arrHistos->At(kJ1_sameppeClPairSep))->Fill(sepCl); 
                }
                tOut->Fill();
            }
        }
        // increase the iterator for the next cluster that will be loaded
        iCl++;
    }

    fOut->cd();
    // output list with histograms
    TList* lOut1 = new TList(); // TH1F
    TList* lOut2 = new TList(); // TH2F
    TList* lOutP = new TList(); // TProfile
    // add all histograms to the lists
    for(Int_t i = kMainJpsi_firstTH1F+1; i < kMainJpsi_firstTH2F; i++) if((TH1F*)arrHistos->At(i)) lOut1->Add((TH1F*)arrHistos->At(i));
    for(Int_t i = kMainJpsi_firstTH2F+1; i < kMainJpsi_firstTH2F; i++) if((TH2F*)arrHistos->At(i)) lOut2->Add((TH2F*)arrHistos->At(i));
    for(Int_t i = kMainJpsi_firstTPrf+1; i < kMainJpsi_all; i++)   if((TProfile*)arrHistos->At(i)) lOutP->Add((TProfile*)arrHistos->At(i));
    lOut1->Write("lTH1F", TObject::kSingleKey);
    lOut2->Write("lTH2F", TObject::kSingleKey);
    lOutP->Write("lTPrf", TObject::kSingleKey);
    // print file content and close it
    fOut->ls();
    // save the file with the output tree
    fOut->Write("",TObject::kWriteDelete);
    delete fOut;

    return;
}

void AnalyzeClPairs()
{
    return;
}

void AnaMain(TString sim)
{
    ConfigLocalAnalysis(sim);
    outSubDir = CreateOutputSubDir();
    gSystem->Exec(Form("mkdir -p %smerged_%s",outDir.Data(),outSubDir.Data()));

    // merge all output files produced by FocalUpcGrid
    MergeOutputFiles();

    // print all histograms from the Grid analysis
    PrintHistograms(Form("%smerged_%sanalysisResultsGrid.root",outDir.Data(),outSubDir.Data()));

    // open merged output file, go over clusters and create cluster pairs
    // store cluster pairs in a tree
    //PrepareClPairs(Form("%s001/%s",outDir.Data(),outSubDir.Data()),kTRUE);
    PrepareClPairs(Form("%smerged_%s",outDir.Data(),outSubDir.Data()));

    // print all histograms from the main analysis
    PrintHistograms(Form("%smerged_%sanalysisResultsMain.root",outDir.Data(),outSubDir.Data()));
}