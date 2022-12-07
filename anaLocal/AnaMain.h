// AnaMain.h
// David Grund, Nov 09, 2022

// root headers
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TTree.h"

// processes
TString processes[7] = {
    "cohJpsi",
    "incJpsi",
    "cohFD",
    "incFD",
    "cohPsi2s",
    "incPsi2s",
    "all"
};
// names
TString names[7] = {
    "coh J/#psi #rightarrow e^{+}e^{-}",
    "inc J/#psi #rightarrow e^{+}e^{-}",
    "coh #psi' #rightarrow J/#psi + #pi^{+}#pi^{-} #rightarrow e^{+}e^{-} + #pi^{+}#pi^{-}",
    "inc #psi' #rightarrow J/#psi + #pi^{+}#pi^{-} #rightarrow e^{+}e^{-} + #pi^{+}#pi^{-}",
    "coh #psi' #rightarrow e^{+}e^{-}",
    "inc #psi' #rightarrow e^{+}e^{-}",
    "ALICE Run-4 Simulation: Pb#minusPb UPC at #sqrt{#it{s}_{NN}} = 5.02 TeV, J/#psi and #psi' #rightarrow e^{+}e^{-}"
};

Int_t SetSimIndex(TString sim)
{
    Int_t iSim = -1;
    if(sim=="cohJpsi") iSim = 0;
    else if(sim=="incJpsi") iSim = 1;
    else if(sim=="cohFD") iSim = 2;
    else if(sim=="incFD") iSim = 3;
    else if(sim=="cohPsi2s") iSim = 4;
    else if(sim=="incPsi2s") iSim = 5;
    else cout << "Unknown process!\n"; 
    return iSim;
}

// ******************************************************************************************************************
// Functions to plot 1d and 2d histograms 
// ******************************************************************************************************************

template <typename TH> // for TH1 and TProfile
void SetHistoLineFill(TH *h, Color_t c, Bool_t filled = kFALSE)
{
    h->SetLineColor(c+1);
    h->SetLineWidth(2);
    if(filled) {
        h->SetFillColor(c);
        h->SetFillStyle(1001);
        h->SetFillColorAlpha(c,0.2);
    } 
}

void SetMarkerProperties(TH1F* h, Color_t c)
{
    gStyle->SetEndErrorSize(1); 
    h->SetMarkerStyle(kFullCircle);
    h->SetMarkerSize(0.7);
    h->SetMarkerColor(c);
    h->SetLineColor(c);
    h->SetLineWidth(2);
    return;
}

template <typename TH> // for TH1 and TProfile
void DrawHisto1D(TH* h, TString sDir, TH* h2 = NULL)
{
    TCanvas c("c","c",700,600);
    // canvas settings
    c.SetLeftMargin(0.11);
    c.SetRightMargin(0.03);
    // x-axis
    h->GetXaxis()->SetDecimals(1);
    h->GetXaxis()->SetTitleOffset(1.2);
    // y-axis
    h->GetYaxis()->SetDecimals(1);
    h->GetYaxis()->SetMaxDigits(3);
    // style
    SetHistoLineFill(h,kBlue,kTRUE);
    // ranges of axes
    h->GetYaxis()->SetRangeUser(0.,h->GetMaximum()*1.05);
    // print the histogram
    if(h2) h->SetBit(TH1::kNoStats);
    c.cd();
    h->Draw();
    if(h2) {
        SetHistoLineFill(h2,kRed,kTRUE);
        h2->SetBit(TH1::kNoStats);
        h2->Draw("SAME");
    }
    TString sName = sDir + h->GetName() + ".pdf";
    c.Print(sName.Data());
    return;
}

template <typename TH> // for TH2F and TProfile2D
void DrawHisto2D(TH* h, TString sDir) 
{
    TCanvas c("c","c",700,600);
    // canvas settings
    c.SetGrid();
    c.SetLogz();
    c.SetLeftMargin(0.11);
    c.SetRightMargin(0.12);
    // x-axis
    h->GetXaxis()->SetDecimals(1);
    h->GetXaxis()->SetTitleOffset(1.2);
    // y-axis
    h->GetYaxis()->SetDecimals(1);
    h->GetYaxis()->SetMaxDigits(3);
    // ranges of axes
    Float_t hMax = h->GetMaximum();
    h->GetZaxis()->SetRangeUser(1.,hMax);
    // print the histogram
    c.cd();
    h->Draw("COLZ");
    TString sName = sDir + h->GetName() + ".pdf";
    c.Print(sName.Data());
    return;
}

template <typename TH> // for TH3F and TProfile2D
void DrawHisto3D(TH* h, TString sDir)
{
    TCanvas c("c","c",700,600);
    // canvas settings
    c.SetGrid();
    c.SetLeftMargin(0.13);
    c.SetRightMargin(0.04);
    // x-axis
    h->GetXaxis()->SetTitleOffset(1.8);
    // y-axis
    h->GetYaxis()->SetTitleOffset(2.0);
    // z-axis
    h->GetZaxis()->SetMaxDigits(3);
    h->GetZaxis()->SetTitleOffset(1.6);
    // no info box
    h->SetBit(TH1::kNoStats);
    // ranges of axes
    Float_t hMax = h->GetMaximum();
    h->GetZaxis()->SetRangeUser(0.,hMax);
    // print the histogram
    c.cd();
    h->Draw("SURF1");
    TString sName = sDir + h->GetName() + ".pdf";
    c.Print(sName.Data());
    return;
}

void SetCanvas(TCanvas* c, Bool_t logScale)
{
    if(logScale) c->SetLogy();
    c->SetTopMargin(0.06);
    c->SetLeftMargin(0.11);
    c->SetRightMargin(0.03);
    return;
}

// ******************************************************************************************************************
// Functions to manage analysis trees
// ******************************************************************************************************************

Int_t fEvNumber, fIdxJEl;
Float_t fEnCl, fXCl, fYCl, fZCl;
TParticle* fJEl = new TParticle();

void SetBranchAddresses_tCls(TTree* t)
{
    t->SetBranchAddress("fEvNumber", &fEvNumber);
    t->SetBranchAddress("fEnCl", &fEnCl);
    t->SetBranchAddress("fXCl", &fXCl);
    t->SetBranchAddress("fYCl", &fYCl);
    t->SetBranchAddress("fZCl", &fZCl);
    t->SetBranchAddress("fIdxJEl", &fIdxJEl);
    t->SetBranchAddress("fJEl", &fJEl);
    Printf("Branch addresses of %s set.", t->GetName());
    return;
}

Float_t fEnClPair, fPtClPair, fEtaClPair, fPhiClPair;
Float_t fEnJElPair, fPtJElPair, fEtaJElPair, fPhiJElPair;

void SetBranchAddresses_tClPairs(TTree* t)
{
    t->SetBranchAddress("fEnClPair", &fEnClPair);
    t->SetBranchAddress("fPtClPair", &fPtClPair);
    t->SetBranchAddress("fEtaClPair", &fEtaClPair);
    t->SetBranchAddress("fPhiClPair", &fPhiClPair);
    t->SetBranchAddress("fEnJElPair", &fEnJElPair);
    t->SetBranchAddress("fPtJElPair", &fPtJElPair);
    t->SetBranchAddress("fEtaJElPair", &fEtaJElPair);
    t->SetBranchAddress("fPhiJElPair", &fPhiJElPair);
    Printf("Branch addresses of %s set.", t->GetName());
    return;
}