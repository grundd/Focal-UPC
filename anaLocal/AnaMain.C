// AnaMain.C
// David Grund, Nov 06, 2022

// root headers
#include "TSystem.h"
#include "TFileMerger.h"
#include "TFile.h"
#include "TList.h"
#include "TCollection.h"
#include "TObject.h"
#include "TLorentzVector.h"
#include "TLegend.h"
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

void DrawHistograms(TString sFile)
{
    TFile* f = TFile::Open(sFile.Data(), "read");
    if(f) Printf("File %s loaded.", f->GetName());
    // get list with TH1F histograms
    TList *lTH1F = (TList*) f->Get("lTH1F");
    if(lTH1F) Printf("List %s loaded.", lTH1F->GetName()); 
    // get list with TH2F histograms
    TList *lTH2F = (TList*) f->Get("lTH2F");
    if(lTH2F) Printf("List %s loaded.", lTH2F->GetName()); 
    // get list with TProfile histograms
    TList *lTPrf = (TList*) f->Get("lTPrf");
    if(lTPrf) Printf("List %s loaded.", lTPrf->GetName()); 
    // get list with TProfile2D histograms
    TList *lTP2D = (TList*) f->Get("lTP2D");
    if(lTP2D) Printf("List %s loaded.", lTP2D->GetName()); 
    // go over the lists and print all histograms
    TString sOut = Form("%smerged_%s",outDir.Data(),outSubDir.Data());

    // TH1F
    Printf("\n***");
    lTH1F->ls();
	TIter nextTH1F(lTH1F);
	TObject* object = NULL;
	while ((object = nextTH1F()))
	{
		cout << "Got an object " << object->GetName() << ":" << endl;
        DrawHisto1D(((TH1F*)object),sOut);
	}

    // TH2F
    Printf("\n***");
    lTH2F->ls();
	TIter nextTH2F(lTH2F);
	object = NULL;
	while ((object = nextTH2F()))
	{
		cout << "Got an object " << object->GetName() << ":" << endl;
        DrawHisto2D(((TH2F*)object),sOut);
	}

    // TProfile
    TProfile* hJP_ppeClX_mtchEn = NULL;
    TProfile* hJP_clX_clEn = NULL;
    TProfile* hJP_ppeClY_mtchEn = NULL;
    TProfile* hJP_clY_clEn = NULL;
    Printf("\n***");
    lTPrf->ls();
	TIter nextTPrf(lTPrf);
	object = NULL;
	while ((object = nextTPrf()))
	{
		cout << "Got an object " << object->GetName() << ":" << endl;
        DrawHisto1D(((TProfile*)object),sOut);
        if((TString)object->GetName() == "hJP_ppeClX_mtchEn") hJP_ppeClX_mtchEn = (TProfile*)object;
        if((TString)object->GetName() == "hJP_clX_clEn") hJP_clX_clEn = (TProfile*)object;
        if((TString)object->GetName() == "hJP_ppeClY_mtchEn") hJP_ppeClY_mtchEn = (TProfile*)object;
        if((TString)object->GetName() == "hJP_clY_clEn") hJP_clY_clEn = (TProfile*)object;
	}
    // draw XY energy distribution: clusters vs MC
    if(hJP_ppeClX_mtchEn) {
        hJP_ppeClX_mtchEn->SetName("h_XEnDist_ClvsMC");
        hJP_ppeClX_mtchEn->SetTitle("X-profile of energy distributions: matched MC tracks (blue) vs clusters (red);#it{x}_{cl};#it{E} [GeV]");
        DrawHisto1D(hJP_ppeClX_mtchEn,sOut,hJP_clX_clEn);
    }
    if(hJP_ppeClY_mtchEn) {
        hJP_ppeClY_mtchEn->SetName("h_YEnDist_ClvsMC");
        hJP_ppeClY_mtchEn->SetTitle("Y-profile of energy distributions: matched MC tracks (blue) vs clusters (red);#it{y}_{cl};#it{E} [GeV]");
        DrawHisto1D(hJP_ppeClY_mtchEn,sOut,hJP_clY_clEn);
    }

    // TProfile2D
    Printf("\n***");
    lTP2D->ls();
    TIter nextTP2D(lTP2D);
    object = NULL;
	while ((object = nextTP2D()))
	{
		cout << "Got an object " << object->GetName() << ":" << endl;
        DrawHisto3D(((TProfile2D*)object),sOut);
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
                // acceptance and efficiency
                ((TH1F*)arrHistos->At(kJ1_mcJElPairRap_rec))->Fill(El12.Rapidity());
                ((TH1F*)arrHistos->At(kJ1_clPairRap_rec))->Fill(Cl12.Rapidity());
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
    TList* lTH1F = new TList(); // TH1F
    TList* lTH2F = new TList(); // TH2F
    TList* lTPrf = new TList(); // TProfile
    TList* lTP2D = new TList(); // TProfile2D
    // add all histograms to the lists
    for(Int_t i = kMainJpsi_firstTH1F+1; i < kMainJpsi_firstTH2F; i++) if((TH1F*)arrHistos->At(i)) lTH1F->Add((TH1F*)arrHistos->At(i));
    for(Int_t i = kMainJpsi_firstTH2F+1; i < kMainJpsi_firstTPrf; i++) if((TH2F*)arrHistos->At(i)) lTH2F->Add((TH2F*)arrHistos->At(i));
    for(Int_t i = kMainJpsi_firstTPrf+1; i < kMainJpsi_firstTP2D; i++) if((TProfile*)arrHistos->At(i)) lTPrf->Add((TProfile*)arrHistos->At(i));
    for(Int_t i = kMainJpsi_firstTP2D+1; i < kMainJpsi_all; i++)     if((TProfile2D*)arrHistos->At(i)) lTP2D->Add((TProfile2D*)arrHistos->At(i));
    lTH1F->Write("lTH1F", TObject::kSingleKey);
    lTH2F->Write("lTH2F", TObject::kSingleKey);
    lTPrf->Write("lTPrf", TObject::kSingleKey);
    lTP2D->Write("lTP2D", TObject::kSingleKey);
    // print file content and close it
    fOut->ls();
    // save the file with the output tree
    fOut->Write("",TObject::kWriteDelete);
    delete fOut;

    return;
}

template <typename TH> // for TH2F and TProfile2D
TH* GetHistoFromList(TString filePath, TString listName, TString histoName, Bool_t debug = kFALSE)
{
    Printf("Accessing the histogram %s:", histoName.Data());
    // open the file
    TFile* f = TFile::Open(filePath.Data(), "read");
    if(f) Printf("File %s loaded.", f->GetName());
    // load the list
    TList *l = (TList*) f->Get(listName.Data());
    if(l) Printf("List %s loaded.", l->GetName());
    if(debug) l->ls();
    // go over the list and look for the histogram
    TIter next(l);
	TObject* object = NULL;
    TH* h = NULL;
	while ((object = next()))
	{
		if(debug) cout << "Found object " << object->GetName() << ". ";
        if((TString)object->GetName() == histoName) {
            h = (TH*)(((TH*)object)->Clone(histoName.Data()));
            cout << "Histogram found!" << endl;
            return h;
        } else if(debug) cout << "Going on..." << endl;
	}
    cout << "Requested histogram not found!" << endl;
    return NULL;
} 

void AxE()
{
    // rapidity distribution of AxE:
    // histogram with all JEl pairs generated in FOCAL rapidity coverage
    TH1F *hJElPairGen = GetHistoFromList<TH1F>(Form("%smerged_%sanalysisResultsGrid.root",outDir.Data(),outSubDir.Data()),
        "lTH1F","hJ1_mcJElPairRap_gen");
    hJElPairGen->SetName("hJElPairGen");
    // histogram with all JEl pairs with both electrons reaching FOCAL
    TH1F *hJElPairAcc = GetHistoFromList<TH1F>(Form("%smerged_%sanalysisResultsGrid.root",outDir.Data(),outSubDir.Data()),
        "lTH1F","hJ1_mcJElPairRap_acc");
    hJElPairAcc->SetName("hJElPairAcc");
    // histogram with rapidity dist of pp ele pairs matched with reconstructed cl pairs
    TH1F *hJElPairRec = GetHistoFromList<TH1F>(Form("%smerged_%sanalysisResultsMain.root",outDir.Data(),outSubDir.Data()),
        "lTH1F","hJ1_mcJElPairRap_rec");
    hJElPairRec->SetName("hJElPairRec");
    // histogram with rapidity dist of reconstructed cl pairs
    TH1F* hClPairRec = GetHistoFromList<TH1F>(Form("%smerged_%sanalysisResultsMain.root",outDir.Data(),outSubDir.Data()),
        "lTH1F","hJ1_clPairRap_rec");
    hClPairRec->SetName("hClPairRec");
    // calculate and draw AxE in given rapidity bins
    TH1F* hAxE = (TH1F*)(hJElPairRec->Clone("hAxE"));
    hAxE->SetTitle("Rapidity dependence of (Acc#times#varepsilon);#it{y} [-];#it{N}_{rec cl pairs} / #it{N}_{gen events}");
    hAxE->SetBit(TH1::kNoStats);
    hAxE->Sumw2();
    hAxE->Divide(hJElPairGen);
    // calculate integrated efficiency
    Float_t NGen(0.), NRec(0.);
    for(Int_t iBin = 1; iBin <= 30; iBin++) {
        NGen += hJElPairGen->GetBinContent(iBin);
        NRec += hClPairRec->GetBinContent(iBin);
    }
    Float_t totalAxE = NRec / NGen;
    // plot rapidity dist of all four histograms
    TCanvas c1("c1","c1",700,600);
    // canvas settings
    //c1.SetLogy();
    c1.SetLeftMargin(0.11);
    c1.SetRightMargin(0.03);
    // x-axis
    hJElPairGen->GetXaxis()->SetTitle("#it{y} [-]");
    hJElPairGen->GetXaxis()->SetDecimals(1);
    hJElPairGen->GetXaxis()->SetTitleOffset(1.2);
    // y-axis
    hJElPairGen->GetYaxis()->SetDecimals(1);
    hJElPairGen->GetYaxis()->SetMaxDigits(3);
    // ranges of axes
    hJElPairGen->GetYaxis()->SetRangeUser(1.,hJElPairGen->GetMaximum()*1.05);
    // print the histograms
    SetHistoLineFill(hJElPairGen,kBlue);
    hJElPairGen->SetBit(TH1::kNoStats);
    hJElPairGen->Draw("HIST");
    SetHistoLineFill(hJElPairAcc,kRed);
    hJElPairAcc->Draw("HIST SAME");
    SetHistoLineFill(hJElPairRec,kGreen);
    hJElPairRec->Draw("HIST SAME");
    SetHistoLineFill(hClPairRec,kViolet);
    hClPairRec->Draw("HIST SAME");
    // legend
    TLegend l1(0.50,0.70,0.95,0.90);
    l1.AddEntry(hJElPairGen,"generated e^{-}e^{+} pairs","L");
    l1.AddEntry(hJElPairAcc,"FOCAL acceptance: 3.4 < #eta^{e^{#pm}} < 5.8","L");
    l1.AddEntry(hJElPairRec,"e^{-}e^{+} pairs matched with rec cl pairs","L");
    l1.AddEntry(hClPairRec,"rec cl pairs","L");
    l1.SetTextSize(0.032);
    l1.SetBorderSize(0);
    l1.SetFillStyle(0);
    l1.SetMargin(0.15);
    l1.Draw();
    // print the canvas
    c1.Print(Form("%smerged_%shAxE_rapDists.pdf",outDir.Data(),outSubDir.Data()));
    // plot AxE
    TCanvas c2("c2","c2",700,600);
    // canvas settings
    c2.SetLeftMargin(0.11);
    c2.SetRightMargin(0.03);
    // x-axis
    hAxE->GetXaxis()->SetDecimals(1);
    hAxE->GetXaxis()->SetTitleOffset(1.2);
    // y-axis
    hAxE->GetYaxis()->SetDecimals(1);
    hAxE->GetYaxis()->SetMaxDigits(3);
    // ranges of axes
    hAxE->GetYaxis()->SetRangeUser(0.,hAxE->GetMaximum()*1.05);
    // print the histogram
    SetHistoLineFill(hAxE,kBlue);
    hAxE->SetBit(TH1::kNoStats);
    hAxE->Draw("");
    // legend
    TLegend l2(0.73,0.85,0.95,0.90);
    l2.AddEntry((TObject*)0,Form("total (Acc#times#varepsilon): %.2f%%",totalAxE*100.),"");
    l2.SetTextSize(0.032);
    l2.SetBorderSize(0);
    l2.SetFillStyle(0);
    l2.SetMargin(0.);
    l2.Draw();
    c2.Print(Form("%smerged_%shAxE.pdf",outDir.Data(),outSubDir.Data()));

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
    DrawHistograms(Form("%smerged_%sanalysisResultsGrid.root",outDir.Data(),outSubDir.Data()));

    // open merged output file, go over clusters and create cluster pairs
    // store cluster pairs in a tree
    PrepareClPairs(Form("%smerged_%s",outDir.Data(),outSubDir.Data()));

    // print all histograms from the main analysis
    DrawHistograms(Form("%smerged_%sanalysisResultsMain.root",outDir.Data(),outSubDir.Data()));

    // calculate rapidity dependence of AxE
    AxE();

    return;
}