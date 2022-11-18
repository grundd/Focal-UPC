// ClusterizeGridDebug_run.C
// David Grund, Nov 18, 2022

// root headers
#include "TSystem.h"
#include "TROOT.h"
#include "TFileMerger.h"

Int_t nFiles = 25;

void PlotHistogram(TH2F* h)
{
    TCanvas c("c","c",700,600);
    c.SetLogz();
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
    TH2F *hYCl_EnCl_CoarseBefore = (TH2F*) f->Get("hYCl_EnCl_CoarseBefore");
    TH2F *hXCl_EnCl_FineBefore = (TH2F*) f->Get("hXCl_EnCl_FineBefore");
    TH2F *hYCl_EnCl_FineBefore = (TH2F*) f->Get("hYCl_EnCl_FineBefore");
    TH2F *hXCl_EnCl_CoarseAfter = (TH2F*) f->Get("hXCl_EnCl_CoarseAfter");
    TH2F *hYCl_EnCl_CoarseAfter = (TH2F*) f->Get("hYCl_EnCl_CoarseAfter");
    TH2F *hXCl_EnCl_FineAfter = (TH2F*) f->Get("hXCl_EnCl_FineAfter");
    TH2F *hYCl_EnCl_FineAfter = (TH2F*) f->Get("hYCl_EnCl_FineAfter");

    // plot all histograms
    PlotHistogram(hCalibratedCls);
    PlotHistogram(hXCl_EnCl_CoarseBefore);
    PlotHistogram(hYCl_EnCl_CoarseBefore);
    PlotHistogram(hXCl_EnCl_FineBefore);
    PlotHistogram(hYCl_EnCl_FineBefore);
    PlotHistogram(hXCl_EnCl_CoarseAfter);
    PlotHistogram(hYCl_EnCl_CoarseAfter);
    PlotHistogram(hXCl_EnCl_FineAfter);
    PlotHistogram(hYCl_EnCl_FineAfter);

    return;
}