// AnaMain.h
// David Grund, Nov 09, 2022

// root headers
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"

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
        h->SetFillColorAlpha(c,0.3);
    } 
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

// ******************************************************************************************************************
// Functions to manage analysis trees
// ******************************************************************************************************************

Int_t fEvNumber, fIdxJEl;
Float_t fEnCl, fPtCl, fEtaCl, fPhiCl;
Float_t fEnJEl, fPtJEl, fEtaJEl, fPhiJEl;

void SetBranchAddresses_tCls(TTree* t)
{
    t->SetBranchAddress("fEvNumber", &fEvNumber);
    t->SetBranchAddress("fEnCl", &fEnCl);
    t->SetBranchAddress("fPtCl", &fPtCl);
    t->SetBranchAddress("fEtaCl", &fEtaCl);
    t->SetBranchAddress("fPhiCl", &fPhiCl);
    t->SetBranchAddress("fIdxJEl", &fIdxJEl);
    t->SetBranchAddress("fEnJEl", &fEnJEl);
    t->SetBranchAddress("fPtJEl", &fPtJEl);
    t->SetBranchAddress("fEtaJEl", &fEtaJEl);
    t->SetBranchAddress("fPhiJEl", &fPhiJEl);
    Printf("Branch addresses of %s set.", t->GetName());
    return;
}