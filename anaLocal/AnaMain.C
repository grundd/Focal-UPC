// AnaMain.C
// David Grund, Nov 06, 2022

// cpp headers
#include <fstream>
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
#include "CreateHistograms.h"
#include "FocalUpcGrid.h"
#include "AnaMain.h"

TString outSubDir = "";
Int_t iSim;

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
    TFile* fOut = new TFile(sDir + "analysisResultsMain.root","RECREATE");
    TTree* tOut = new TTree("tClPairs", "output tree containing cluster pairs");
    // pairs of summed clusters/superclusters
    tOut->Branch("fEnClPair", &fEnClPair, "fEnClPair/F");
    tOut->Branch("fPtClPair", &fPtClPair, "fPtClPair/F");
    tOut->Branch("fEtaClPair", &fEtaClPair, "fEtaClPair/F");
    tOut->Branch("fPhiClPair", &fPhiClPair, "fPhiClPair/F");
    // pairs of J/psi electrons (if cluster pairs were matched with it)
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
        TLorentzVector* thisCl = new TLorentzVector(ConvertXYZEtoLorVec(fXCl,fYCl,fZCl,fEnCl));
        // create 4-vector of the matched physical primary electron
        TParticle* thisEl = new TParticle(*fJEl);
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
        // vectors to store xCl, yCl and zCl
        std::vector<Float_t> vectX; std::vector<Float_t> vectY; std::vector<Float_t> vectZ;
        vectX.push_back(fXCl); vectY.push_back(fYCl); vectZ.push_back(fZCl);
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
            TLorentzVector* nextCl = new TLorentzVector(ConvertXYZEtoLorVec(fXCl,fYCl,fZCl,fEnCl));
            // create 4-vector of the matched physical primary electron
            TParticle* nextEl = new TParticle(*fJEl);
            // add them to the lists
            lCls.AddLast(nextCl);
            lEls.AddLast(nextEl);
            // store the electron index
            idxMtchJEl.push_back(fIdxJEl);
            // store xCl, yCl and zCl
            vectX.push_back(fXCl); vectY.push_back(fYCl); vectZ.push_back(fZCl);
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
            TParticle* El1 = (TParticle*)lEls.At(iCl1);
            if(!El1) continue;     
            // go over all possible pairs of clusters
            for(Int_t iCl2 = iCl1+1; iCl2 < nClThisEvent; iCl2++) 
            {
                TLorentzVector* Cl2 = (TLorentzVector*)lCls.At(iCl2);
                if(!Cl2) continue;
                TParticle* El2 = (TParticle*)lEls.At(iCl2);
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
                TLorentzVector lvEl1;
                lvEl1.SetPxPyPzE(El1->Px(),El1->Py(),El1->Pz(),El1->Energy());
                TLorentzVector lvEl2;
                lvEl2.SetPxPyPzE(El2->Px(),El2->Py(),El2->Pz(),El2->Energy());
                TLorentzVector lvEl12 = lvEl1 + lvEl2;
                // electron pair kinematics -> tree
                fEnJElPair = lvEl12.Energy();
                fPtJElPair = lvEl12.Pt();
                fEtaJElPair = lvEl12.Eta();
                fPhiJElPair = lvEl12.Phi();
                // fill some kinematic histograms
                ((TH1F*)arrHistos->At(kJ1_clPairEn))->Fill(fEnClPair);
                ((TH1F*)arrHistos->At(kJ1_clPairPt))->Fill(fPtClPair);
                if(clPairM > cutMLowPtDist && clPairM < cutMUppPtDist) ((TH1F*)arrHistos->At(kJ1_clPairPt_massCut))->Fill(fPtClPair);
                ((TH1F*)arrHistos->At(kJ1_clPairRap))->Fill(Cl12.Rapidity());
                ((TH1F*)arrHistos->At(kJ1_clPairM))->Fill(clPairM);      
                // acceptance and efficiency
                if(3.4 < El1->Eta() && El1->Eta() < 5.8 && 3.4 < El2->Eta() && El2->Eta() < 5.8) {
                    ((TH1F*)arrHistos->At(kJ1_rap_rec_genKine))->Fill(lvEl12.Rapidity());
                    ((TH1F*)arrHistos->At(kJ1_rap_rec))->Fill(Cl12.Rapidity());
                    // inv mass cut on (sup)cl pairs
                    if(clPairM > cutMLow && clPairM < cutMUpp) {
                        ((TH1F*)arrHistos->At(kJ1_rap_sig_genKine))->Fill(lvEl12.Rapidity());
                        ((TH1F*)arrHistos->At(kJ1_rap_sig))->Fill(Cl12.Rapidity());
                    }
                }
                // radial separation between the pairs of clusters
                Float_t sepCl = TMath::Sqrt(TMath::Power(vectX[iCl1]-vectX[iCl2],2) + TMath::Power(vectY[iCl1]-vectY[iCl2],2));
                ((TH1F*)arrHistos->At(kJ1_clPairSep))->Fill(sepCl);
                // radial separation between cluster pair vs between the pair of ppe
                Float_t x1(0.), y1(0.);
                TrackCoordinatesAtZ(El1,vectZ[iCl1],x1,y1);
                Float_t x2(0.), y2(0.);
                TrackCoordinatesAtZ(El2,vectZ[iCl2],x2,y2);
                Float_t sepMC = TMath::Sqrt(TMath::Power(x1-x2,2) + TMath::Power(y1-y2,2));
                ((TH2F*)arrHistos->At(kJ2_clPairSep_mcJElSep))->Fill(sepCl,sepMC);
                // if both clusters are paired to a different physical primary electron
                if(idxMtchJEl[iCl1] != idxMtchJEl[iCl2]) {
                    if(idxMtchJEl[iCl1] == -1 || idxMtchJEl[iCl2] == -1) continue;
                    ((TH1F*)arrHistos->At(kJ1_ppeClPairM))->Fill(clPairM);
                    ((TH1F*)arrHistos->At(kJ1_ppeClPairSep))->Fill(sepCl); 
                    ((TH2F*)arrHistos->At(kJ2_ppeClPairEn_mtchEn))->Fill(fEnClPair,lvEl12.Energy());
                    ((TH2F*)arrHistos->At(kJ2_ppeClPairRap_mtchRap))->Fill(Cl12.Rapidity(),lvEl12.Rapidity());
                    ((TH2F*)arrHistos->At(kJ2_ppeClPairPt_mtchPt))->Fill(fPtClPair,lvEl12.Pt());
                    ((TH2F*)arrHistos->At(kJ2_ppeClPairM_mtchM))->Fill(clPairM,lvEl12.M());
                } else {
                    ((TH1F*)arrHistos->At(kJ1_sameppeClPairSep))->Fill(sepCl); 
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
    // rapidity distributions contributing to AxE:
    // rapidity of all generated JEl pairs
    TH1F *hRap_gen = GetHistoFromList<TH1F>(Form("%smerged_%sanalysisResultsGrid.root",outDir.Data(),outSubDir.Data()),
        "lTH1F","hJ1_rap_gen");
    hRap_gen->SetName("hRap_gen");
    // rapidity of generated JEl pairs, where both electrons reach FOCAL acceptance (3.4 < eta < 5.8)
    TH1F *hRap_acc = GetHistoFromList<TH1F>(Form("%smerged_%sanalysisResultsGrid.root",outDir.Data(),outSubDir.Data()),
        "lTH1F","hJ1_rap_acc");
    hRap_acc->SetName("hRap_acc");
    // rapidity of generated JEl pairs in FOCAL acc matched to reconstructed cl pairs
    TH1F *hRap_rec_genKine = GetHistoFromList<TH1F>(Form("%smerged_%sanalysisResultsMain.root",outDir.Data(),outSubDir.Data()),
        "lTH1F","hJ1_rap_rec_genKine");
    hRap_rec_genKine->SetName("hRap_rec_genKine");
    // rapidity of reconstructed cl pairs matched to generated JEl pairs in FOCAL acc
    TH1F *hRap_rec_detKine = GetHistoFromList<TH1F>(Form("%smerged_%sanalysisResultsMain.root",outDir.Data(),outSubDir.Data()),
        "lTH1F","hJ1_rap_rec");
    hRap_rec_detKine->SetName("hRap_rec_detKine");
    // rapidity of generated JEl pairs in FOCAL acc matched to reconstructed cl pairs, with a cut on m_{cl pair}
    TH1F *hRap_sig_genKine = GetHistoFromList<TH1F>(Form("%smerged_%sanalysisResultsMain.root",outDir.Data(),outSubDir.Data()),
        "lTH1F","hJ1_rap_sig_genKine");
    hRap_sig_genKine->SetName("hRap_sig_genKine");
    // rapidity of reconstructed cl pairs matched to generated JEl pairs in FOCAL acc, with a cut on m_{cl pair}
    TH1F *hRap_sig_detKine = GetHistoFromList<TH1F>(Form("%smerged_%sanalysisResultsMain.root",outDir.Data(),outSubDir.Data()),
        "lTH1F","hJ1_rap_sig");
    hRap_sig_detKine->SetName("hRap_sig_detKine");

    // ********************************************
    // use distributions with generator kinematics?
    bool genKine = true;
    TH1F* hRap_rec = NULL;
    TH1F* hRap_sig = NULL;
    if(genKine) {
        hRap_rec = hRap_rec_genKine;
        hRap_sig = hRap_sig_genKine;
    } else {
        hRap_rec = hRap_rec_detKine;
        hRap_sig = hRap_sig_detKine;
    }
    // ********************************************

    // calculate and draw AxE in given rapidity bins
    // only acceptance
    TH1F* hA = (TH1F*)(hRap_acc->Clone("hA"));
    hA->SetTitle("Rapidity dependence of acceptance;#it{y} [-];Acc = #it{N}_{FOCAL acc} / #it{N}_{gen}");
    hA->Sumw2();
    hA->Divide(hRap_gen);
    // only efficiency
    TH1F* hERec = hERec = (TH1F*)(hRap_rec->Clone("hERec"));
    hERec->SetTitle("Rapidity dependence of rec. efficiency;#it{y} [-];#it{#varepsilon} = #it{N}_{rec} / #it{N}_{FOCAL acc}");
    hERec->Sumw2();
    hERec->Divide(hRap_acc);
    // only signal eff
    TH1F* hESig = (TH1F*)(hRap_sig->Clone("hESig"));
    hESig->SetTitle("Rapidity dependence of sig. efficiency;#it{y} [-];#it{#varepsilon} = #it{N}_{sig} / #it{N}_{rec}");
    hESig->Sumw2();
    hESig->Divide(hRap_rec);
    // full AxE
    TH1F* hAxE = (TH1F*)(hRap_sig->Clone("hAxE"));
    hAxE->SetTitle("Rapidity dependence of Acc#times#it{#varepsilon};#it{y} [-];Acc#times#it{#varepsilon} = #it{N}_{sig} / #it{N}_{gen}");
    hAxE->Sumw2();
    hAxE->Divide(hRap_gen);

    // calculate the integrated (total) efficiency
    Float_t NGen(0.), NSig(0.);
    for(Int_t iBin = 1; iBin <= 30; iBin++) {
        NGen += hRap_gen->GetBinContent(iBin);
        NSig += hRap_sig->GetBinContent(iBin);
    }
    Float_t totalAxE = NSig / NGen;
    // print the value
    ofstream of(Form("%smerged_%stotalAxE.txt",outDir.Data(),outSubDir.Data()));
    of << totalAxE;
    of.close();

    // plot rapidity dist of all four histograms
    TCanvas c1("c1","c1",700,600);
    SetCanvas(&c1,kFALSE);
    // x-axis
    hRap_gen->GetXaxis()->SetDecimals(1);
    hRap_gen->GetXaxis()->SetLabelOffset(0.01);
    hRap_gen->GetXaxis()->SetTitle("#it{y} (-)");
    hRap_gen->GetXaxis()->SetTitleOffset(1.2);
    // y-axis
    hRap_gen->GetYaxis()->SetDecimals(1);
    hRap_gen->GetYaxis()->SetMaxDigits(3);
    // ranges of axes
    hRap_gen->GetYaxis()->SetRangeUser(0.,hRap_gen->GetMaximum()*1.05);
    // print the histograms
    SetHistoLineFill(hRap_gen,kBlue,kTRUE);
    hRap_gen->SetBit(TH1::kNoStats);
    hRap_gen->SetBit(TH1::kNoTitle);
    hRap_gen->Draw("HIST");
    SetHistoLineFill(hRap_acc,kRed,kTRUE);
    hRap_acc->Draw("HIST SAME");
    SetHistoLineFill(hRap_rec,kGreen,kTRUE);
    hRap_rec->Draw("HIST SAME");
    SetHistoLineFill(hRap_sig,kMagenta,kTRUE);
    hRap_sig->Draw("HIST SAME");
    // title
    TLatex* ltx = new TLatex();
    ltx->SetTextSize(0.032);
    ltx->SetTextAlign(21);
    ltx->SetNDC();
    ltx->DrawLatex(0.55,0.96,"ALICE Run-4 Simulation: Pb#minusPb UPC at #sqrt{#it{s}_{NN}} = 5.02 TeV");
    // legend
    Int_t nRows1 = 5;
    TLegend l1(0.42,0.925-nRows1*0.04,0.95,0.925);
    l1.AddEntry((TObject*)0,Form("#bf{%s}",names[iSim].Data()),"");
    l1.AddEntry(hRap_gen,"generated e^{+}e^{-} pairs","L");
    l1.AddEntry(hRap_acc,"e^{+}e^{-} pairs with 3.4 < #eta^{e^{#pm}} < 5.8","L");
    if(genKine) {
        l1.AddEntry(hRap_rec,"rec e^{+}e^{-} pairs","L");
        l1.AddEntry(hRap_sig,Form("rec e^{+}e^{-} pairs with %.1f < #it{m} < %.1f GeV/#it{c}^{2}",cutMLow,cutMUpp),"L");
    } else {
        l1.AddEntry(hRap_rec,Form("rec %s pairs",sCl.Data()),"L");
        l1.AddEntry(hRap_sig,Form("rec %s pairs with %.1f < #it{m} < %.1f GeV/#it{c}^{2}",sCl.Data(),cutMLow,cutMUpp),"L");
    }
    l1.SetTextSize(0.032);
    l1.SetBorderSize(0);
    l1.SetFillStyle(0);
    l1.SetMargin(0.14);
    l1.Draw();
    // print the canvas
    if(genKine) c1.Print(Form("%smerged_%shAxE_genKine_rapDists.pdf",outDir.Data(),outSubDir.Data()));
    else        c1.Print(Form("%smerged_%shAxE_rapDists.pdf",outDir.Data(),outSubDir.Data()));
    
    // plot AxE
    TCanvas c2("c2","c2",700,600);
    SetCanvas(&c2,kFALSE);
    c2.SetTopMargin(0.265);
    // x-axis
    hA->GetXaxis()->SetDecimals(1);
    hA->GetXaxis()->SetLabelOffset(0.01);
    hA->GetXaxis()->SetTitleOffset(1.2);
    // y-axis
    hA->GetYaxis()->SetDecimals(1);
    hA->GetYaxis()->SetMaxDigits(3);
    hA->GetYaxis()->SetTitle("Acc#times#it{#varepsilon} = #it{N}_{sig} / #it{N}_{gen}");
    // style
    gStyle->SetEndErrorSize(1); 
    hA->SetBit(TH1::kNoStats);
    hA->SetBit(TH1::kNoTitle);
    hA->GetYaxis()->SetRangeUser(0.,1.05);
    SetMarkerProperties(hA,kRed+1);
    SetMarkerProperties(hERec,kGreen+1);
    SetMarkerProperties(hESig,kMagenta+1);
    SetMarkerProperties(hAxE,kBlue+1);
    // print the histogram
    hA->Draw("E1");
    hERec->Draw("E1 SAME");
    hESig->Draw("E1 SAME");
    hAxE->Draw("E1 SAME");
    // legends
    Int_t nRows2 = 5;
    TLegend l2(0.12,0.94-nRows2*0.04,0.40,0.94);
    l2.AddEntry((TObject*)0,Form("#bf{%s}",names[iSim].Data()),"");
    l2.AddEntry(hA,"acceptance","EPL");
    l2.AddEntry(hERec,"reconstruction eff","EPL");
    l2.AddEntry(hESig,"signal eff","EPL");
    l2.AddEntry(hAxE,Form("Acc#times#it{#varepsilon} (total = %.1f%%)",totalAxE*100.),"EPL");
    l2.SetTextSize(0.032);
    l2.SetBorderSize(0);
    l2.SetFillStyle(0);
    l2.SetMargin(0.16);
    l2.Draw();
    Int_t nRows3 = 2;
    if(doSupercls) nRows3++;
    if(cutE > 0.)  nRows3++;
    if(cutM > 0.)  nRows3++;
    TLegend l3(0.60,0.94-nRows3*0.04,0.90,0.94);
    l3.AddEntry((TObject*)0,"#bf{selections}:","");
    if(doSupercls) l3.AddEntry((TObject*)0,Form("superclusterizer (%.0f GeV, %.0f cm)",minSeedE,radius),"");
    if(cutE > 0.)  l3.AddEntry((TObject*)0,Form("#it{E}_{%s} > %.0f GeV",sCl.Data(),cutE),"");
    if(cutM > 0.)  l3.AddEntry((TObject*)0,Form("mass filtering (%.1f)",cutM),"");
    l3.AddEntry((TObject*)0,Form("%.1f < #it{m}_{%s pair} < %.1f GeV/#it{c}^{2}",cutMLow,sCl.Data(),cutMUpp),"");
    l3.SetTextSize(0.032);
    l3.SetBorderSize(0);
    l3.SetFillStyle(0);
    l3.SetMargin(0.0);
    l3.Draw();
    // title
    ltx->DrawLatex(0.5,0.96,"ALICE Run-4 Simulation: Pb#minusPb UPC at #sqrt{#it{s}_{NN}} = 5.02 TeV");
    if(genKine) c2.Print(Form("%smerged_%shAxE_genKine.pdf",outDir.Data(),outSubDir.Data()));
    else        c2.Print(Form("%smerged_%shAxE.pdf",outDir.Data(),outSubDir.Data()));
    
    return;
}

void AnaMain(TString sim)
{
    ConfigLocalAnalysis(sim);
    outSubDir = CreateOutputSubDir();
    iSim = SetSimIndex(sim);
    gSystem->Exec(Form("mkdir -p %smerged_%s",outDir.Data(),outSubDir.Data()));
    Bool_t isBoxSim = kFALSE;
    if(sim=="boxEle" || sim=="boxPho") isBoxSim = kTRUE;

    // merge all output files produced by FocalUpcGrid
    MergeOutputFiles();

    // print all histograms from the Grid analysis
    DrawHistograms(Form("%smerged_%sanalysisResultsGrid.root",outDir.Data(),outSubDir.Data()));

    if(!isBoxSim) 
    {
        // open merged output file, go over clusters and create cluster pairs
        // store cluster pairs in a tree
        PrepareClPairs(Form("%smerged_%s",outDir.Data(),outSubDir.Data()));

        // print all histograms from the main analysis
        DrawHistograms(Form("%smerged_%sanalysisResultsMain.root",outDir.Data(),outSubDir.Data()));

        // calculate rapidity dependence of AxE
        AxE();
    }

    return;
}