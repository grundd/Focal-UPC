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
void DrawHisto1D(TH* h, TString sDir, TH* h2 = NULL)
{
    TCanvas c("c","c",700,600);
    // canvas settings
    c.SetRightMargin(0.03);
    // x-axis
    h->GetXaxis()->SetDecimals(1);
    h->GetXaxis()->SetTitleOffset(1.2);
    // y-axis
    h->GetYaxis()->SetDecimals(1);
    h->GetYaxis()->SetMaxDigits(3);
    // style
    h->SetLineColor(kBlue+1);
    h->SetLineWidth(2);
    h->SetFillColor(kBlue);
    h->SetFillStyle(1001);
    h->SetFillColorAlpha(kBlue,0.3);
    // ranges of axes
    Float_t hMax = h->GetMaximum();
    h->GetYaxis()->SetRangeUser(0.,hMax*1.05);
    // print the histogram
    TString sName = sDir + h->GetName() + ".pdf";
    c.cd();
    h->Draw();
    if(h2) {
        h->SetBit(TH1::kNoTitle);
        h->SetBit(TH1::kNoStats);
        h2->SetLineColor(kRed+1);
        h2->SetLineWidth(2);
        h2->SetFillColor(kRed);
        h2->SetFillStyle(1001);
        h2->SetFillColorAlpha(kRed,0.3);
        h2->Draw("SAME");
    }
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
    TString sName = sDir + h->GetName() + ".pdf";
    c.cd();
    h->Draw("COLZ");
    c.Print(sName.Data());
    return;
}

template <typename TH> // for TH3F and TProfile2D
void DrawHisto3D(TH* h, TString sDir)
{
    TCanvas c("c","c",700,600);
    // canvas settings
    c.SetGrid();
    c.SetRightMargin(0.04);
    c.SetLeftMargin(0.13);
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
    TString sName = sDir + h->GetName() + ".pdf";
    c.cd();
    h->Draw("SURF1");
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