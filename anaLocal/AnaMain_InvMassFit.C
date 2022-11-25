// FocalUpc_InvMassFit.C
// David Grund, Nov 01, 2022

// root headers
#include "TLatex.h"
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
#include "ConfigAnalysis.h"
#include "ConfigParameters.h"
#include "AnaMain.h"

using namespace RooFit;

// processes
TString processes[6] = {
    "cohJpsi",
    "incJpsi",
    "cohFD",
    "incFD",
    "cohPsi2s",
    "incPsi2s"
};
// names
TString names[6] = {
    "coh J/#psi #rightarrow e^{+}e^{-}",
    "inc J/#psi #rightarrow e^{+}e^{-}",
    "coh #psi' #rightarrow J/#psi + #pi^{+}#pi^{-} #rightarrow e^{+}e^{-} + #pi^{+}#pi^{-}",
    "inc #psi' #rightarrow J/#psi + #pi^{+}#pi^{-} #rightarrow e^{+}e^{-} + #pi^{+}#pi^{-}",
    "coh #psi' #rightarrow e^{+}e^{-}",
    "inc #psi' #rightarrow e^{+}e^{-}",
};
// expected number of photoproduced VMs in FOCAL rapidity converage:
// 3.4 < y < 5.8
// obtained from StarlightRapDep.C
Float_t expNFocRap[6] = {
    965000, // coh J/psi
    715000, // inc J/psi
    95000,  // coh FD
    74000,  // inc FD
    20000,  // coh psi'
    16000   // inc psi'
};
// paths
TString outSubDir = "";
TString sOut = "";
// mass binning
Int_t nBinsM;
Float_t fMLow = 1.4; // GeV
Float_t fMUpp = 4.2;
Float_t fMStep = 0.04;
// cuts
Float_t fPtLow = 0.0;
Float_t fPtUpp = 0.2;

Float_t roundFloat(Float_t N)
{
    Int_t divisor = 1;
    if(N > 1e6) divisor = 100000;
    else if(N > 1e5) divisor = 10000;
    else if(N > 1e4) divisor = 1000;
    else if(N > 1e3) divisor = 100;
    else if(N > 1e2) divisor = 10;
    return (Int_t)N - ((Int_t)N % divisor);
}

void PrintHistogram(TH1F *h)
{
    TCanvas c("c","c",700,600);
    // canvas settings
    c.SetLeftMargin(0.11);
    c.SetRightMargin(0.03);
    // x-axis
    h->GetXaxis()->SetTitle("#it{m}_{cl pair} [GeV/#it{c}^{2}]");
    h->GetXaxis()->SetDecimals(1);
    h->GetXaxis()->SetTitleOffset(1.2);
    // y-axis
    h->GetYaxis()->SetTitle(Form("counts per %.0f MeV [-]", fMStep*1e3));
    h->GetYaxis()->SetDecimals(1);
    h->GetYaxis()->SetMaxDigits(3);
    // ranges of axes
    h->GetYaxis()->SetRangeUser(0.,h->GetMaximum()*1.05);
    // print the histogram
    SetHistoLineFill(h,kBlue);
    h->SetBit(TH1::kNoStats);
    h->SetBit(TH1::kNoTitle);
    h->Draw("HIST");
    // title
    TLatex* latex = new TLatex();
    latex->SetTextSize(0.032);
    latex->SetTextAlign(21);
    latex->SetNDC();
    latex->DrawLatex(0.50,0.95,Form("%s",h->GetTitle()));
    // legend
    TLegend l(0.18,0.72,0.36,0.88);
    l.AddEntry((TObject*)0,"#it{L}_{int} = 6.9 nb^{-1}","");
    l.AddEntry((TObject*)0,Form("total: %.0f events",h->Integral()),"");
    l.AddEntry((TObject*)0,"3.4 < #it{y}_{VM} < 5.8","");
    l.AddEntry((TObject*)0,"3.4 < #eta_{e^{#pm}} < 5.8","");
    l.SetTextSize(0.032);
    l.SetBorderSize(0);
    l.SetFillStyle(0);
    l.SetMargin(0.);
    l.Draw();
    c.Print(Form("%s%s.pdf",sOut.Data(),h->GetName()));
    return;
}

void PrepareHist(Int_t iPcs, TH1F* h, TH1F* hRnd)
{
    ConfigLocalAnalysis(processes[iPcs]);
    // get the file with the tree
    TFile* f = TFile::Open(Form("%smerged_%sanalysisResultsMain.root",outDir.Data(),outSubDir.Data()),"read");
    if(f) Printf("File %s loaded.", f->GetName());
    TTree* t = dynamic_cast<TTree*>(f->Get("tClPairs"));
    if(t) Printf("Tree %s loaded.", t->GetName());
    SetBranchAddresses_tClPairs(t);
    // go over tree entries and fill a histogram with masses
    Int_t nClPairs = t->GetEntries();
    Float_t nNoRapCut = 0;
    Float_t nRapCut = 0;
    for(Int_t iClPair = 0; iClPair < nClPairs; iClPair++)
    {
        t->GetEntry(iClPair);
        TLorentzVector clPair;
        clPair.SetPtEtaPhiE(fPtClPair,fEtaClPair,fEtaClPair,fEnClPair);
        nNoRapCut++;
        if(clPair.Pt() > fPtLow && clPair.Pt() < fPtUpp) {
            nRapCut++;
            h->Fill(clPair.M());
        }
    }
    // get the total AxE
    Float_t totalAxE;
    ifstream ifs(Form("%smerged_%stotalAxE.txt",outDir.Data(),outSubDir.Data()));
    ifs >> totalAxE;
    ifs.close();
    // scale the histogram to the expected yields
    Float_t currIntegral = h->Integral();
    Float_t currNFocRap = currIntegral / (totalAxE * nRapCut / nNoRapCut);
    Float_t scale = expNFocRap[iPcs] / currNFocRap;
    // generate new histogram with the expected # of events
    Int_t nToGenerate = (Int_t)roundFloat(currIntegral * scale);
    for(Int_t i = 0; i < nToGenerate; i++) {
        hRnd->Fill(h->GetRandom());
    }
    // scale the histogram to the expected # of events
    h->Scale(roundFloat(currIntegral * scale) / h->Integral());
    Float_t newIntegral = h->Integral();
    Float_t newAxE = newIntegral / expNFocRap[iPcs];
    // print the info
    cout << "\n**********************************\n"
         << " process: " << processes[iPcs] << "\n"
         << " before scaling: \n"
         << "  - integral: " << currIntegral << "\n"
         << "  - total AxE: " << totalAxE << "\n"
         << "  - pT cut eff: " << nRapCut / nNoRapCut << "\n"
         << "  - N in FOC rap: " << currNFocRap << "\n"
         << " scale: " << scale << "\n"
         << " after scaling: \n"
         << "  - integral: " << newIntegral << "\n"
         << "  - total AxE * pT cut eff: " << newAxE << "\n"
         << "  - N in FOC rap: " << newIntegral / newAxE << "\n";
    // print the histograms
    PrintHistogram(h);
    PrintHistogram(hRnd);
    return;
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
    hCorr->SetMarkerSize(2.0);
    hCorr->GetXaxis()->SetLabelSize(0.08);
    hCorr->GetYaxis()->SetLabelSize(0.08);
    hCorr->Draw("colz,text");

    return;
}

void SetCanvas(TCanvas* c, Bool_t logScale)
{
    if(logScale) c->SetLogy();
    c->SetTopMargin(0.05);
    c->SetLeftMargin(0.11);
    c->SetRightMargin(0.03);
    return;
}

void DoFit(Int_t iPcs, TH1F* h)
{
    // fit the invariant mass distribution of the signal using Double sided CB function

    // roofit variables
    RooRealVar fM("fM","fM",fMLow,fMUpp);

    // histogram with data
    RooDataHist fDataHist("fDataHist","fDataHist",fM,Import(*h));

    // print the number of events in the histogram
    Float_t nEv = h->Integral();
    cout << "Number of events in the histogram: " << nEv << endl;

    // roofit variables
    Float_t meanVal, meanLow, meanUpp;
    if(iPcs == 0) { meanVal = 3.097; meanLow = 2.8; meanUpp = 3.3; }
    if(iPcs == 4) { meanVal = 3.686; meanLow = 3.4; meanUpp = 3.9; }
    RooRealVar mean("#mu","m",meanVal,meanLow,meanUpp);
    RooRealVar sigmaL("#sigma_{L}", "sigma_{L}",0.1, 0.01, 0.20);
    RooRealVar sigmaR("#sigma_{R}", "sigma_{R}",0.1, 0.01, 0.20);
    RooRealVar nL("#it{n}_{L}","n_{L}",1.,0.01,100.);
    RooRealVar nR("#it{n}_{R}","n_{R}",1.,0.01,100.);
    RooRealVar alphaL("#alpha_{L}","alpha_{L}",1.,0.01,40.);
    RooRealVar alphaR("#alpha_{R}","alpha_{R}",1.,0.01,40.);
    RooCrystalBall DSCB("DSCB", "Double sided CB function", fM, mean, sigmaL, sigmaR, alphaL, nL, alphaR, nR);
    // fix some parameters?
    Bool_t fix = kFALSE;
    if(fix) {
        // ...
    }
    // perform the fit
    // extended fit?
    Bool_t ext = kTRUE;
    RooRealVar N("#it{N}","N",nEv,nEv*0.8,nEv*1.2);
    RooExtendPdf ExtDSCB("ExtDSCB","Extended Double sided CB function",DSCB,N);
    RooFitResult* fResFit;
    if(ext) fResFit = ExtDSCB.fitTo(fDataHist,SumW2Error(kFALSE),Extended(kTRUE),Save()); 
    else    fResFit = DSCB.fitTo(fDataHist,SumW2Error(kFALSE),Save());

    // plot the results
    // draw correlation matrix
    TCanvas* cCM = new TCanvas("cCM","cCM",700,600);
    DrawCM(cCM,fResFit);

    // draw histogram and fit
    TCanvas* c = new TCanvas("c","c",700,600);
    SetCanvas(c,kFALSE);
    
    RooPlot* fr = fM.frame(Title("inv mass fit")); 
    fDataHist.plotOn(fr,Name("fDataHist"),MarkerStyle(20),MarkerSize(1.));
    if(ext) ExtDSCB.plotOn(fr,Name("ExtDSCB"),LineColor(215),LineWidth(3),LineStyle(9));
    else    DSCB.plotOn(fr,Name("DSCB"),LineColor(215),LineWidth(3),LineStyle(9));
    // Y axis
    Float_t hMax = h->GetMaximum();
    //fr->GetYaxis()->SetRangeUser(0.0,hMax*1.05);
    // title
    fr->GetYaxis()->SetTitle(Form("Counts per %.0f MeV/#it{c}^{2}", fMStep*1e3));
    fr->GetYaxis()->SetTitleSize(0.032);
    fr->GetYaxis()->SetTitleOffset(1.4);
    // label
    fr->GetYaxis()->SetLabelSize(0.032);
    fr->GetYaxis()->SetLabelOffset(0.01);
    fr->GetYaxis()->SetMaxDigits(3);
    // X axis
    // title
    fr->GetXaxis()->SetTitle("#it{m}_{cl pair} [GeV/#it{c}^{2}]");
    fr->GetXaxis()->SetTitleSize(0.032);
    fr->GetXaxis()->SetTitleOffset(1.25);
    // label
    fr->GetXaxis()->SetLabelSize(0.032);
    fr->GetXaxis()->SetLabelOffset(0.01);
    fr->GetYaxis()->SetDecimals(1);
    fr->Draw();

    // calculate chi2 
    Float_t chi2;
    if(ext) chi2 = fr->chiSquare("ExtDSCB","fDataHist",fResFit->floatParsFinal().getSize());
    else    chi2 = fr->chiSquare("DSCB","fDataHist",fResFit->floatParsFinal().getSize());
    Printf("********************");
    Printf("chi2/NDF = %.3f", chi2);
    Printf("NDF = %i", fResFit->floatParsFinal().getSize());
    Printf("chi2/NDF = %.3f/%i", chi2*fResFit->floatParsFinal().getSize(), fResFit->floatParsFinal().getSize());
    Printf("********************");
    
    TLegend *l = new TLegend(0.16,0.54,0.42,0.94);
    l->AddEntry((TObject*)0,Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}",fPtLow,fPtUpp),"");
    l->AddEntry((TObject*)0,Form("#chi^{2}/NDF = %.3f",chi2),"");
    l->AddEntry((TObject*)0,Form("#it{N} = %.f #pm %.f", N.getVal(), N.getError()),"");
    l->AddEntry((TObject*)0,Form("#mu = %.2f GeV/#it{c}^{2}", mean.getVal()),"");
    l->AddEntry((TObject*)0,Form("#sigma_{L} = %.2f GeV/#it{c}^{2}", sigmaL.getVal()),"");
    l->AddEntry((TObject*)0,Form("#sigma_{R} = %.2f GeV/#it{c}^{2}", sigmaR.getVal()),""); 
    l->AddEntry((TObject*)0,Form("#alpha_{L} = %.1f", alphaL.getVal()),"");
    l->AddEntry((TObject*)0,Form("#alpha_{R} = %.1f", alphaR.getVal()),"");
    l->AddEntry((TObject*)0,Form("#it{n}_{L} = %.1f", nL.getVal()),"");
    l->AddEntry((TObject*)0,Form("#it{n}_{R} = %.1f", nR.getVal()),"");
    l->SetTextSize(0.032);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetMargin(0.);
    l->Draw();
    // save the plots
    cCM->Print(Form("%sfit_%s_CM.pdf",sOut.Data(),processes[iPcs].Data()));
    c->Print(Form("%sfit_%s.pdf",sOut.Data(),processes[iPcs].Data()));

    // draw histogram with log scale
    TCanvas *cLog = new TCanvas("cLog","cLog",700,600);
    SetCanvas(cLog,kTRUE);
    fr->Draw();
    l->Draw();
    cLog->Print(Form("%sfit_%s_log.pdf",sOut.Data(),processes[iPcs].Data()));
    
    delete cCM;
    delete c;
    delete cLog;
    return;
}

void AnaMain_InvMassFit()
{
    // set the binning
    nBinsM = (fMUpp-fMLow) / fMStep;
    Printf("%i mass bins will be used.", nBinsM);

    // set the paths
    outSubDir = CreateOutputSubDir();
    sOut = Form("results/sim%02i_g%02i_p%02i/invMassFit/",simFiles,fileGeom,filePara);
    gSystem->Exec("mkdir -p " + sOut);

    // fill inv mass histograms
    TH1F* h[6] = { NULL };
    TH1F* hRnd[6] = { NULL };
    TH1F* hAll = new TH1F("h_all","ALICE Run-4 Simulation: Pb#minusPb UPC at #sqrt{#it{s}_{NN}} = 5.02 TeV, J/#psi and #psi' #rightarrow e^{+}e^{-}",nBinsM,fMLow,fMUpp);
    TH1F* hRndAll = new TH1F("hRnd_all","ALICE Run-4 Simulation: Pb#minusPb UPC at #sqrt{#it{s}_{NN}} = 5.02 TeV, J/#psi and #psi' #rightarrow e^{+}e^{-}",nBinsM,fMLow,fMUpp);
    for(Int_t i = 0; i < 6; i++) {
        // define the histograms
        h[i] = new TH1F(Form("h_%s",processes[i].Data()),names[i],nBinsM,fMLow,fMUpp);
        hRnd[i] = new TH1F(Form("hRnd_%s",processes[i].Data()),names[i],nBinsM,fMLow,fMUpp);
        // fill the histograms
        PrepareHist(i,h[i],hRnd[i]);
        hAll->Add(h[i]);
        hRndAll->Add(hRnd[i]);
    }
    PrintHistogram(hAll);
    PrintHistogram(hRndAll);

    // fit coherent J/psi
    DoFit(0,hRnd[0]);

    // fit coherent psi'
    DoFit(4,hRnd[4]);

    // fit the combined sample
    return;
}