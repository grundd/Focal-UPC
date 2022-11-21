// ClusterizeGridDebug_run.C
// David Grund, Nov 18, 2022

// root headers
#include "TSystem.h"
#include "TROOT.h"
#include "TFileMerger.h"

Int_t nFiles = 25;

void PlotHistogram1D(TH1F* h)
{
    TCanvas c("c","c",700,600);
    c.SetLeftMargin(0.12);
    c.SetRightMargin(0.03);
    //c.SetLogz();
    h->SetBit(TH1::kNoStats);
    h->SetLineColor(kBlue);
    h->SetLineWidth(2);
    h->Draw("");
    c.Print(Form("results/sim02_g02_p02/clusterizeGridDebug/%s.pdf",h->GetName()));
    return;
}

void PlotHistogram2D(TH2F* h)
{
    TCanvas c("c","c",700,600);
    c.SetLeftMargin(0.12);
    c.SetRightMargin(0.11);
    //c.SetLogz();
    h->SetBit(TH2::kNoStats);
    h->Draw("COLZ");
    c.Print(Form("results/sim02_g02_p02/clusterizeGridDebug/%s.pdf",h->GetName()));
    return;
}

void CompareBeforeAfterCalib(TString sBefore, TString sAfter, TString sEntries, Bool_t isX)
{
    TFile* f = TFile::Open("results/sim02_g02_p02/clusterizeGridDebug/mergedClusterizerHistograms.root", "read");
    if(f) Printf("File %s loaded.", f->GetName());
    TH2F* hBefore = (TH2F*) f->Get(sBefore.Data());
    TH2F* hAfter = (TH2F*) f->Get(sAfter.Data());
    TH1F* hNEntries = new TH1F(sEntries.Data(),"",100,-50.,+50.);
    if(isX) {
        hBefore->SetTitle("Energy distribution of cls per segment vs #it{x} before calibration;#it{x}_{cl per seg} [cm];#it{E}_{cl per seg} [GeV]");
        hAfter->SetTitle("Energy distribution of cls per segment vs #it{x} after calibration;#it{x}_{cl per seg} [cm];#it{E}_{cl per seg} [GeV]");
        hNEntries->SetTitle("Number of cls per segment vs #it{x};#it{x}_{cl per seg} [cm];#it{N}_{cl per seg}");
    } else {
        hBefore->SetTitle("Energy distribution of cls per segment vs #it{y} before calibration;#it{y}_{cl per seg} [cm];#it{E}_{cl per seg} [GeV]");
        hAfter->SetTitle("Energy distribution of cls per segment vs #it{y} after calibration;#it{y}_{cl per seg} [cm];#it{E}_{cl per seg} [GeV]");
        hNEntries->SetTitle("Number of cls per segment vs #it{y};#it{y}_{cl per seg} [cm];#it{N}_{cl per seg}");
    }   
    for(Int_t iRow = 1; iRow <= 100; iRow++) {
        Int_t nPerCol = 0;
        for(Int_t iCol = 1; iCol <= 100; iCol++) nPerCol += hBefore->GetBinContent(iRow,iCol);
        hNEntries->SetBinContent(iRow,nPerCol);
    }
    PlotHistogram2D(hBefore);
    PlotHistogram2D(hAfter);
    PlotHistogram1D(hNEntries);
    delete hNEntries;
    f->Close();
    return;
}

void ClusterizeGridDebug_run()
{
    gSystem->Load("libpythia6_4_28.so");

    gROOT->ProcessLine(".L ClusterizeGridDebug.C");

    for(Int_t i = 0; i < nFiles; i++) 
    {
        TString sIn = Form("inputData/aliDPG_v02/kCohJpsiToElRad/%03i/",i+1);
        TString sOut = Form("results/sim02_g02_p02/clusterizeGridDebug/%03i/",i+1);
        gSystem->Exec("mkdir -p " + sOut);
        cout << sOut << ": " << endl;
        if(gSystem->AccessPathName(Form("%sfocalClusters.root",sOut.Data()))) {
            cout << " Running clusterizer:" << endl;
            TString sCmd = Form("ClusterizeGridDebug(kTRUE,\"geometry_02.txt\",\"parameters_02.txt\",\"%s\",\"%s\")",sIn.Data(),sOut.Data());
            gROOT->ProcessLine(sCmd.Data());
        } else {
            cout << " Clusters already produced! " << endl;
        }
    }
    cout << endl << " FINISHED! " << endl;

    // merge all output files
    TFileMerger m;
    m.OutputFile("results/sim02_g02_p02/clusterizeGridDebug/mergedClusterizerHistograms.root");
    for(Int_t i = 0; i < nFiles; i++) 
    {
        TString sFile = Form("results/sim02_g02_p02/clusterizeGridDebug/%03i/clusterizerHistograms.root",i+1);
        cout << "Adding file: " << sFile << endl;
        m.AddFile(Form("%s",sFile.Data()));
    }
    m.Merge();
    cout << endl << " FILES MERGED! " << endl;

    // print all histograms
    TFile* f = TFile::Open("results/sim02_g02_p02/clusterizeGridDebug/mergedClusterizerHistograms.root", "read");
    if(f) Printf("File %s loaded.", f->GetName());
    TH2F *hCalibratedCls = (TH2F*) f->Get("hCalibratedCls");
    hCalibratedCls->SetTitle(";#it{x}_{cl per seg} [cm];#it{y}_{cl per seg} [cm]");
    PlotHistogram2D(hCalibratedCls);
    
    CompareBeforeAfterCalib("hXCl_EnCl_CoarseBefore","hXCl_EnCl_CoarseAfter","hXCl_NEntries_Coarse",kTRUE);
    CompareBeforeAfterCalib("hYCl_EnCl_CoarseBefore","hYCl_EnCl_CoarseAfter","hYCl_NEntries_Coarse",kFALSE);
    CompareBeforeAfterCalib("hXCl_EnCl_FineBefore","hXCl_EnCl_FineAfter","hXCl_NEntries_Fine",kTRUE);
    CompareBeforeAfterCalib("hYCl_EnCl_FineBefore","hYCl_EnCl_FineAfter","hYCl_NEntries_Fine",kFALSE);

    return;
}