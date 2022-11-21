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
    c.SetLeftMargin(0.13);
    c.SetRightMargin(0.05);
    //c.SetLogz();
    h->Draw("");
    c.Print(Form("results/sim02_g02_p02/clusterizeGridDebug/%s.pdf",h->GetName()));
    return;
}

void PlotHistogram2D(TH2F* h)
{
    TCanvas c("c","c",700,600);
    //c.SetLogz();
    h->Draw("COLZ");
    c.Print(Form("results/sim02_g02_p02/clusterizeGridDebug/%s.pdf",h->GetName()));
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
    TH2F *hXCl_EnCl_CoarseBefore = (TH2F*) f->Get("hXCl_EnCl_CoarseBefore");
    TH2F *hXCl_EnCl_CoarseAfter = (TH2F*) f->Get("hXCl_EnCl_CoarseAfter");
    TH1F *hXCl_NEntries_Coarse = new TH1F("hXCl_NEntries_Coarse","",100,-50.,+50.);
    hXCl_NEntries_Coarse->SetTitle("Number of cls per segment at given #it{x};#it{x}_{cl};#it{N}_{cl per segment}");
    for(Int_t iRow = 1; iRow <= 100; iRow++) {
        Int_t nPerCol = 0;
        for(Int_t iCol = 1; iCol <= 100; iCol++) nPerCol += hXCl_EnCl_CoarseBefore->GetBinContent(iRow,iCol);
        hXCl_NEntries_Coarse->SetBinContent(iRow,nPerCol);
    }
    TH2F *hXCl_EnCl_FineBefore = (TH2F*) f->Get("hXCl_EnCl_FineBefore");
    TH2F *hXCl_EnCl_FineAfter = (TH2F*) f->Get("hXCl_EnCl_FineAfter");
    TH1F *hXCl_NEntries_Fine = new TH1F("hXCl_NEntries_Fine","",100,-50.,+50.);
    hXCl_NEntries_Fine->SetTitle("Number of cls per segment at given #it{x};#it{x}_{cl};#it{N}_{cl per segment}");
    for(Int_t iRow = 1; iRow <= 100; iRow++) {
        Int_t nPerCol = 0;
        for(Int_t iCol = 1; iCol <= 100; iCol++) nPerCol += hXCl_EnCl_FineBefore->GetBinContent(iRow,iCol);
        hXCl_NEntries_Fine->SetBinContent(iRow,nPerCol);
    }
    TH2F *hYCl_EnCl_CoarseBefore = (TH2F*) f->Get("hYCl_EnCl_CoarseBefore");
    TH2F *hYCl_EnCl_CoarseAfter = (TH2F*) f->Get("hYCl_EnCl_CoarseAfter");
    TH1F *hYCl_NEntries_Coarse = new TH1F("hYCl_NEntries_Coarse","",100,-50.,+50.);
    hYCl_NEntries_Coarse->SetTitle("Number of cls per segment at given #it{y};#it{y}_{cl};#it{N}_{cl per segment}");
    for(Int_t iRow = 1; iRow <= 100; iRow++) {
        Int_t nPerCol = 0;
        for(Int_t iCol = 1; iCol <= 100; iCol++) nPerCol += hYCl_EnCl_CoarseBefore->GetBinContent(iRow,iCol);
        hYCl_NEntries_Coarse->SetBinContent(iRow,nPerCol);
    }
    TH2F *hYCl_EnCl_FineBefore = (TH2F*) f->Get("hYCl_EnCl_FineBefore");
    TH2F *hYCl_EnCl_FineAfter = (TH2F*) f->Get("hYCl_EnCl_FineAfter");
    TH1F *hYCl_NEntries_Fine = new TH1F("hYCl_NEntries_Fine","",100,-50.,+50.);
    hYCl_NEntries_Fine->SetTitle("Number of cls per segment at given #it{y};#it{y}_{cl};#it{N}_{cl per segment}");
    for(Int_t iRow = 1; iRow <= 100; iRow++) {
        Int_t nPerCol = 0;
        for(Int_t iCol = 1; iCol <= 100; iCol++) nPerCol += hYCl_EnCl_FineBefore->GetBinContent(iRow,iCol);
        hYCl_NEntries_Fine->SetBinContent(iRow,nPerCol);
    }

    // plot all histograms
    PlotHistogram2D(hCalibratedCls);
    PlotHistogram2D(hXCl_EnCl_CoarseBefore);
    PlotHistogram2D(hXCl_EnCl_CoarseAfter);
    PlotHistogram1D(hXCl_NEntries_Coarse);
    PlotHistogram2D(hXCl_EnCl_FineBefore);
    PlotHistogram2D(hXCl_EnCl_FineAfter);
    PlotHistogram1D(hXCl_NEntries_Fine);
    PlotHistogram2D(hYCl_EnCl_CoarseBefore);
    PlotHistogram2D(hYCl_EnCl_CoarseAfter);
    PlotHistogram1D(hYCl_NEntries_Coarse);
    PlotHistogram2D(hYCl_EnCl_FineBefore);
    PlotHistogram2D(hYCl_EnCl_FineAfter);
    PlotHistogram1D(hYCl_NEntries_Fine);

    return;
}