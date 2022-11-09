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
void DrawHisto(TH* h, TString sDir)
{
    TCanvas c("c","c",700,600);
    h->GetYaxis()->SetMaxDigits(3);
    h->GetXaxis()->SetTitleOffset(1.2);
    h->SetLineColor(kBlue+1);
    h->SetFillColor(kBlue);
    h->SetFillStyle(3012);
    TString sName = sDir + h->GetName() + ".pdf";
    c.cd();
    h->Draw();
    c.Print(sName.Data());
    return;
}

void DrawHistoCOLZ(TH2F* h, TString sDir) // for TH2F
{
    TCanvas c("c","c",700,600);
    c.SetGrid();
    c.SetLogz();
    h->GetYaxis()->SetMaxDigits(3);
    h->GetXaxis()->SetTitleOffset(1.2);
    Float_t hMax = h->GetMaximum();
    h->GetZaxis()->SetRangeUser(1.,hMax);
    TString sName = sDir + h->GetName() + ".pdf";
    c.cd();
    h->Draw("COLZ");
    c.Print(sName.Data());
    return;
}

// ******************************************************************************************************************
// Functions to manage analysis trees
// ******************************************************************************************************************

Int_t fEvNumber;
Float_t fEnCl, fPtCl, fEtaCl, fPhiCl;
Float_t fEnJEl, fPtJEl, fEtaJEl, fPhiJEl;

void SetBranchAddresses_tCls(TTree* t)
{
    t->SetBranchAddress("fEvNumber", &fEvNumber);
    t->SetBranchAddress("fEnCl", &fEnCl);
    t->SetBranchAddress("fPtCl", &fPtCl);
    t->SetBranchAddress("fEtaCl", &fEtaCl);
    t->SetBranchAddress("fPhiCl", &fPhiCl);
    t->SetBranchAddress("fEnJEl", &fEnJEl);
    t->SetBranchAddress("fPtJEl", &fPtJEl);
    t->SetBranchAddress("fEtaJEl", &fEtaJEl);
    t->SetBranchAddress("fPhiJEl", &fPhiJEl);
    Printf("Branch addresses of %s set.", t->GetName());
    return;
}