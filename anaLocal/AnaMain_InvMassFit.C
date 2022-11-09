// FocalUpc_InvMassFit.C
// David Grund, Nov 01, 2022

// roofit headers
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooGenericPdf.h"
#include "RooBinning.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
// my headers
#include "FocalUpc_Utilities.h"
#include "FocalUpc_ConfigAnalysis.h"

using namespace RooFit;

Float_t fM,fPt,fRap;
Bool_t fixNParameters = kFALSE;
Bool_t fitWithDSCB = kTRUE;

void PrepareTree(TString sIn)
{
    TFile* f = TFile::Open(sIn.Data(),"read");
    if(f)
    {
        Printf("File %s already created.", sIn.Data());
        return;
    } 
    else
    { 
        Printf("File %s will be created.", sIn.Data());

        TFile* fIn = TFile::Open((sOut + "merged/_mergedTrees.root").Data(), "read");
        if(fIn) Printf("Input file %s loaded.", fIn->GetName());
        TTree* tIn = dynamic_cast<TTree*>(fIn->Get("tMergedClPairs"));
        if(tIn) Printf("Input tree %s loaded.", tIn->GetName());
        SetBranchAddresses_tClPairs(tIn);

        // create new tree
        f = new TFile(sIn.Data(),"RECREATE");
        TTree *tOut = new TTree("tInvMassFit", "tInvMassFit");
        tOut->Branch("fM", &fM, "fM/F");
        tOut->Branch("fPt", &fPt, "fPt/F");
        tOut->Branch("fRap", &fRap, "fRap/F");

        Int_t nEntries = tIn->GetEntries();
        Printf("%i entries found in %s.", nEntries, tIn->GetName());
        Float_t progress = 0.; // perc

        for(Int_t iEntry = 0; iEntry < nEntries; iEntry++)
        {
            tIn->GetEntry(iEntry);
            // update the progress bar
            if((iEntry+1) % (Int_t)(nEntries/10.) == 0) {
                progress += 10.;
                cout << "[" << progress << "%] done." << endl;
            }
            TLorentzVector lvClPair;
            lvClPair.SetPtEtaPhiE(fPtClPair,fEtaClPair,fPhiClPair,fEnClPair);
            fM = lvClPair.M();
            fPt = lvClPair.Pt();
            fRap = lvClPair.Rapidity();
            tOut->Fill();
        }

        f->Write("",TObject::kWriteDelete);
        return;
    }
}

void DrawCM(TCanvas* cCM, RooFitResult* fResFit)
{
    cCM->SetTopMargin(0.05);
    cCM->SetRightMargin(0.12);
    cCM->SetLeftMargin(0.12);

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2f");

    TH2* hCorr = fResFit->correlationHist();

    if(!fixNParameters) hCorr->GetXaxis()->SetBinLabel(7,"#sigma");
    else                   hCorr->GetXaxis()->SetBinLabel(5,"#sigma");
    hCorr->GetYaxis()->SetBinLabel(1,"#sigma");
    hCorr->SetMarkerSize(2.0);
    hCorr->GetXaxis()->SetLabelSize(0.08); // 0.049
    hCorr->GetYaxis()->SetLabelSize(0.08);
    hCorr->Draw("colz,text");

    return;
}

void SetCanvas(TCanvas* c, Bool_t logScale)
{
    if(logScale) c->SetLogy();
    c->SetTopMargin(0.03);
    c->SetBottomMargin(0.12);
    c->SetRightMargin(0.03);
    c->SetLeftMargin(0.13);

    return;
}

void DoFit(TString sSim)
{
    // fit the invariant mass distribution of the signal using Double-sided CB function

    // cuts:
    char fStrReduce[120];
    Float_t fMCutLow  = 1.;
    Float_t fMCutUpp  = 4.;
    Float_t fPtCutLow = -1e3;
    Float_t fPtCutUpp = -1e3;
    Float_t fRapCutLow = fEtaFOClow;
    Float_t fRapCutUpp = fEtaFOCupp;

    if(sSim == "cohJpsi") {
        fPtCutLow = 0.0;
        fPtCutUpp = 1.0;
    }
    else if(sSim == "incJpsi") {
        fPtCutLow = 0.2;
        fPtCutUpp = 1.0;    
    }
    sprintf(fStrReduce,"fM>%f && fM<%f && fPt>%f && fPt<%f && fRap>%f && fRap<%f",fMCutLow,fMCutUpp,fPtCutLow,fPtCutUpp,fRapCutLow,fRapCutUpp);

    // binning:
    Int_t nBinsM = 30;
    RooBinning binningM(nBinsM,fMCutLow,fMCutUpp);
    Float_t binSize = (fMCutUpp - fMCutLow) * 1000 / nBinsM; // in MeV
    cout << "Bin size: " << binSize << " MeV" << endl;

    // roofit variables
    RooRealVar fM("fM","fM",fMCutLow,fMCutUpp);
    RooRealVar fPt("fPt","fPt",0.,10.);
    RooRealVar fRap("fRap","fRap",2.,7.);

    //fM.setBinning(binM);

    // get the tree
    TFile* fIn = TFile::Open((sOut + "invMassFit/tInvMassFit.root").Data(), "read");
    if(fIn) Printf("Input file %s loaded.", fIn->GetName());
    TTree* tIn = dynamic_cast<TTree*>(fIn->Get("tInvMassFit"));
    if(tIn) Printf("Input tree %s loaded.", tIn->GetName());

    RooDataSet* fDataIn = new RooDataSet("fDataIn","fDataIn",RooArgSet(fM,fPt,fRap),Import(*tIn));
    RooAbsData* fDataSet = fDataIn->reduce(fStrReduce);

    // print the number of entries in the dataset
    Int_t nEv = fDataSet->numEntries();
    cout << "Number of events in the dataset: " << nEv << endl;
    
    // roofit variables
    RooRealVar norm_L("norm_L","N_{L}(J/#psi)",nEv,0,1e04);
    RooRealVar norm_R("norm_R","N_{R}(J/#psi)",nEv,0,1e04);

    RooRealVar mean_L("m","m_{J/#psi}",3.097,2.8,3.3);
    RooRealVar sigma_L("sig","#sigma_{J/#psi}",0.1,0.01,0.2);
    RooRealVar alpha_L("#alpha_{L}","alpha_{L}",1.,0.,20.);
    RooRealVar n_L("n_{L}","n_{L}",1.,0.,40.);

    RooGenericPdf mean_R("mean_R","m_{J/#psi}","m",RooArgSet(mean_L));
    RooGenericPdf sigma_R("sigma_R","#sigma_{J/#psi}","sig",RooArgSet(sigma_L));
    RooRealVar alpha_R("#alpha_{R}","alpha_{R}",-1.,-20.,0.);
    RooRealVar n_R("n_{R}","n_{R}",1.,0.,40.);
    
    if(fixNParameters)
    {
        n_L.setVal(10.);
        n_R.setVal(10.);
        n_L.setConstant(kTRUE);
        n_R.setConstant(kTRUE);
    }

    RooCBShape CB_L("CB_L","CB_L",fM,mean_L,sigma_L,alpha_L,n_L);
    RooCBShape CB_R("CB_R","CB_R",fM,mean_R,sigma_R,alpha_R,n_R);
    RooRealVar frac("frac","fraction of CBs",0.5);
    RooAddPdf DoubleSidedCB("DoubleSidedCB","DoubleSidedCB",RooArgList(CB_L,CB_R),RooArgList(frac));

    // create the model
    RooRealVar norm("norm","N(J/#psi)",nEv,0,1e04);
    RooExtendPdf ExtendedDSCB("ExtendedDSCB","Extended DSCB function",DoubleSidedCB,norm);
    RooExtendPdf ExtendedCB("ExtendedCB","Extended CB function",CB_L,norm);
    // perform the fit
    RooFitResult* fResFit = NULL;
    if(fitWithDSCB) fResFit = ExtendedDSCB.fitTo(*fDataSet,Extended(kTRUE),Range(fMCutLow,fMCutUpp),Save());
    else            fResFit = ExtendedCB.fitTo(*fDataSet,Extended(kTRUE),Range(fMCutLow,fMCutUpp),Save());

    // ##########################################################
    // plot the results
    // draw correlation matrix
    TCanvas* cCM = new TCanvas("cCM","cCM",700,600);
    DrawCM(cCM,fResFit);

    // draw histogram and fit
    TCanvas *c = new TCanvas("c","c",700,600);
    SetCanvas(c,kFALSE);
    
    RooPlot* fr = fM.frame(Title("invariant mass fit")); 
    fDataSet->plotOn(fr,Name("fDataSet"),Binning(binningM),MarkerStyle(20),MarkerSize(1.));
    if(fitWithDSCB) ExtendedDSCB.plotOn(fr,Name("ExtendedDSCB"),LineColor(215),LineWidth(3),LineStyle(9));
    else            ExtendedCB.plotOn(fr,Name("ExtendedCB"),LineColor(215),LineWidth(3),LineStyle(9));
    // Y axis
    // title
    fr->GetYaxis()->SetTitle(Form("Counts per %.0f MeV/#it{c}^{2}", binSize));
    fr->GetYaxis()->SetTitleSize(0.045);
    fr->GetYaxis()->SetTitleOffset(1.35);
    // label
    fr->GetYaxis()->SetLabelSize(0.045);
    fr->GetYaxis()->SetLabelOffset(0.01);
    fr->GetYaxis()->SetMaxDigits(3);
    // X axis
    // title
    fr->GetXaxis()->SetTitle("#it{m}_{cl pair} [GeV/#it{c}^{2}]");
    fr->GetXaxis()->SetTitleSize(0.045);
    fr->GetXaxis()->SetTitleOffset(1.1);
    // label
    fr->GetXaxis()->SetLabelSize(0.045);
    fr->GetXaxis()->SetLabelOffset(0.01);
    fr->GetYaxis()->SetDecimals(0);
    fr->Draw();

    // get chi2 
    Float_t chi2;
    if(fitWithDSCB) chi2 = fr->chiSquare("ExtendedDSCB","fDataSet",fResFit->floatParsFinal().getSize());
    else            chi2 = fr->chiSquare("ExtendedCB","fDataSet",fResFit->floatParsFinal().getSize());
    Printf("********************");
    Printf("chi2/NDF = %.3f", chi2);
    Printf("NDF = %i", fResFit->floatParsFinal().getSize());
    Printf("chi2/NDF = %.3f/%i", chi2*fResFit->floatParsFinal().getSize(), fResFit->floatParsFinal().getSize());
    Printf("********************");   

    TLegend *l = new TLegend(0.10,0.38,0.40,0.96);
    l->AddEntry((TObject*)0,Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}",fPtCutLow,fPtCutUpp),"");
    l->AddEntry((TObject*)0,Form("%.1f < #it{y} < %.1f",fRapCutLow,fRapCutUpp),"");
    l->AddEntry((TObject*)0,Form("#it{E}_{cl} > %.0f GeV",cutE),"");
    l->AddEntry((TObject*)0,Form("#chi^{2}/NDF = %.3f",chi2),"");
    l->AddEntry((TObject*)0,Form("#it{N} = %.f #pm %.f", norm.getVal(), norm.getError()),"");
    l->AddEntry((TObject*)0,Form("#mu = %.2f GeV/#it{c}^{2}", mean_L.getVal()),""); // mean_L.getError()
    l->AddEntry((TObject*)0,Form("#sigma = %.2f GeV/#it{c}^{2}", sigma_L.getVal()),""); // sigma_L.getError()
    l->AddEntry((TObject*)0,Form("#alpha_{L} = %.2f #pm %.2f", alpha_L.getVal(), alpha_L.getError()),"");
    if(fitWithDSCB) l->AddEntry((TObject*)0,Form("#alpha_{R} = %.2f #pm %.2f", (-1)*(alpha_R.getVal()), alpha_R.getError()),"");
    if(!fixNParameters) {
        l->AddEntry((TObject*)0,Form("#it{n}_{L} = %.0f #pm %.0f", n_L.getVal(), n_L.getError()),"");
        if(fitWithDSCB) l->AddEntry((TObject*)0,Form("#it{n}_{R} = %.0f #pm %.0f", n_R.getVal(), n_R.getError()),"");
    } else {
        l->AddEntry((TObject*)0,Form("#it{n}_{L} = %.1f", n_L.getVal()),"");
        if(fitWithDSCB) l->AddEntry((TObject*)0,Form("#it{n}_{R} = %.1f", n_R.getVal()),"");
    }
    l->SetTextSize(0.042);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->Draw();

    // draw histogram with log scale
    TCanvas *cLog = new TCanvas("cLog","cLog",800,600);
    SetCanvas(cLog,kTRUE);
    fr->Draw();

    // save the plots
    c->Print(Form("%sinvMassFit/fit.pdf",sOut.Data()));
    cLog->Print(Form("%sinvMassFit/fit_log.pdf",sOut.Data()));
    cCM->Print(Form("%sinvMassFit/fit_CM.pdf",sOut.Data()));

    delete c;
    delete cLog;
    delete cCM;
    return;
}

void FocalUpc_InvMassFit(TString sSim)
{
    if(!ConfigAnalysis(sSim)) {
        cout << "Wrong configuration. Terminating..." << endl;
        return;
    }
    gSystem->Exec("mkdir -p " + sOut + "invMassFit/");

    // prepare the data tree that will be used for fitting
    TString sIn = sOut + "invMassFit/tInvMassFit.root";
    PrepareTree(sIn);

    // do the invariant mass fit
    DoFit(sSim);

    return;
}