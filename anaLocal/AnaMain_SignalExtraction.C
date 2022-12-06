// AnaMain_SignalExtraction.C
// David Grund, Nov 01, 2022

// cpp headers
#include <fstream>
#include <iomanip> // std::setprecision()
// root headers
#include "TSystem.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLorentzVector.h"
// roofit headers
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooCrystalBall.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
// my headers
#include "ConfigAnalysis.h"
#include "ConfigParameters.h"
#include "AnaMain.h"

using namespace RooFit;

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
// expected number of photoproduced VMs in FOCAL rapidity converage:
// 3.4 < y < 5.8
// obtained from StarlightRapDep.C
Float_t nFocRapExpected[6] = {
    965000, // coh J/psi
    715000, // inc J/psi
    95000,  // coh FD
    74000,  // inc FD
    20000,  // coh psi'
    16000   // inc psi'
};
// paths
TString outSubDir = "";
TString sOutMass = "";
TString sOutPt = "";
// mass binning
Int_t nBinsM;
Float_t fMEdgeLow = 1.4; // GeV
Float_t fMEdgeUpp = 4.2;
Float_t fMStep = 0.04;
// pT binning
Int_t nBinsPt;
Float_t fPtEdgeLow = 0.; // GeV
Float_t fPtEdgeUpp = 2.;
Float_t fPtStep = 0.05;
// cuts: inv mass fit
Float_t fPtLow = 0.0;
Float_t fPtUpp = 0.2;
// cuts: pT distribution
Float_t fMLow = 2.8;
Float_t fMUpp = 3.4;
// reduce the total scale by a factor of:
Float_t fReduce = 1.;

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

void PrintHistogram(TH1F *h, TString sOut)
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
    TLatex* ltx = new TLatex();
    ltx->SetTextSize(0.032);
    ltx->SetTextAlign(21);
    ltx->SetNDC();
    ltx->DrawLatex(0.50,0.95,Form("%s",h->GetTitle()));
    // legend
    TLegend l(0.35,0.72,0.55,0.88);
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

void ScaleHistos(TH1F* h, TH1F* hRnd, Int_t iPcs, Float_t fAxE, Float_t fCutEff, TString sOut)
{
    // scale the histogram to the expected yield
    Float_t intHistNow = h->Integral();
    Float_t nFocRapNow = intHistNow / fAxE / fCutEff;
    Float_t fScale = nFocRapExpected[iPcs] / nFocRapNow;
    fScale = fScale / fReduce;
    // generate new histogram with the expected # of events
    Int_t nToGene = (Int_t)roundFloat(intHistNow * fScale);
    for(Int_t i = 0; i < nToGene; i++) {
        hRnd->Fill(h->GetRandom());
    }
    // scale the histogram to the expected # of events
    h->Scale(roundFloat(intHistNow * fScale) / h->Integral());
    Float_t intHistNew = h->Integral();
    Float_t fAxENew = intHistNew / nFocRapExpected[iPcs];
    // print the info
    cout << "\n**********************************\n"
         << " process: " << processes[iPcs] << "\n"
         << " before scaling: \n"
         << "  - integral: " << intHistNow << "\n"
         << "  - total AxE: " << fAxE << "\n"
         << "  - pT/mass cut eff: " << fCutEff << "\n"
         << "  - N in FOC rap: " << nFocRapNow << "\n"
         << " scale: " << fScale << "\n"
         << " after scaling: \n"
         << "  - integral: " << intHistNew << "\n"
         << "  - total AxE * pT cut eff: " << fAxENew << "\n"
         << "  - N in FOC rap: " << intHistNew / fAxENew << "\n";
    // print the histograms
    PrintHistogram(h,sOut);
    PrintHistogram(hRnd,sOut);
    return;
}

void PrepareHistos(Int_t iPcs, TH1F* hMass, TH1F* hMassRnd, TH1F* hPt, TH1F* hPtRnd)
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
    Float_t nAll(0.), nPtCut(0.), nMassCut(0.);
    for(Int_t iClPair = 0; iClPair < nClPairs; iClPair++)
    {
        t->GetEntry(iClPair);
        TLorentzVector clPair;
        clPair.SetPtEtaPhiE(fPtClPair,fEtaClPair,fEtaClPair,fEnClPair);
        // inv mass dist:
        nAll++;
        if(clPair.Pt() > fPtLow && clPair.Pt() < fPtUpp) {
            nPtCut++;
            hMass->Fill(clPair.M());
        }
        // pT dist:
        if(clPair.M() > fMLow && clPair.M() < fMUpp) {
            nMassCut++;
            hPt->Fill(clPair.Pt());
        }
    }
    // get the total AxE
    Float_t fAxE;
    ifstream ifs(Form("%smerged_%stotalAxE.txt",outDir.Data(),outSubDir.Data()));
    ifs >> fAxE;
    ifs.close();
    // scale the histograms to the expected yields
    // inv mass fit hist
    Float_t fPtCutEff = nPtCut / nAll;
    ScaleHistos(hMass,hMassRnd,iPcs,fAxE,fPtCutEff,sOutMass);
    // pT hist
    Float_t fMassCutEff = nMassCut / nAll;
    ScaleHistos(hPt,hPtRnd,iPcs,fAxE,fMassCutEff,sOutPt);

    return;
}

void DrawCM(TCanvas* c, RooFitResult* fResFit)
{
    // canvas settings
    c->SetTopMargin(0.05);
    c->SetRightMargin(0.12);
    c->SetLeftMargin(0.12);
    // style
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2f");
    // corr matrix
    TH2* h = fResFit->correlationHist();
    h->SetMarkerSize(1.5);
    // x-axis
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetLabelOffset(0.01);
    // y-axis
    h->GetYaxis()->SetLabelSize(0.06);
    // z-axis
    h->GetZaxis()->SetDecimals(1);
    h->Draw("colz,text");
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

void LoadFitParams(TString s, Float_t& aL, Float_t& aR, Float_t& m, Float_t& nL, Float_t& nR, Float_t& sL, Float_t& sR)
{
    ifstream ifs(s.Data());
    Float_t aL_e, aR_e, m_e, nL_e, nR_e, sL_e, sR_e;
    ifs >> aL >> aL_e
        >> aR >> aR_e
        >> m >> m_e
        >> nL >> nL_e
        >> nR >> nR_e
        >> sL >> sL_e
        >> sR >> sR_e;
    ifs.close();
    cout << "Loaded values of the parameters: \n"
         << "aL: " << aL << "\n"
         << "aR: " << aR << "\n"
         << "m:  " << m << "\n"
         << "nL: " << nL << "\n"
         << "nR: " << nR << "\n"
         << "sL: " << sL << "\n"
         << "sR: " << sR << "\n";
    return;
}

void DoInvMassFit(Int_t iPcs, TH1F* h)
{
    // fit the invariant mass distribution of the signal using Double sided CB function
    // if iPcs == 6 => fit to the combined sample

    // roofit variables
    RooRealVar fM("fM","fM",fMEdgeLow,fMEdgeUpp);

    // histogram with data
    RooDataHist fDataHist("fDataHist","fDataHist",fM,Import(*h));

    // print the number of events in the histogram
    Float_t nEv = h->Integral();
    cout << "Number of events in the histogram: " << nEv << endl;

    // initial values of the parameters
    // J/psi CB
    Float_t aL(10.), aR(15.0), m, nL(10.), nR(0.5), sL(.2), sR(.2);
    // psi' CB
    Float_t aL2(1.), aR2(1.), m2, nL2(1.), nR2(1.), sL2(.1), sR2(.1);
    // fit to the combined sample => load the values from the previous fits
    if(iPcs==6) {
        LoadFitParams(Form("%sparam_cohJpsi.txt",sOutMass.Data()), aL, aR, m, nL, nR, sL, sR);
        LoadFitParams(Form("%sparam_cohPsi2s.txt",sOutMass.Data()), aL2, aR2, m2, nL2, nR2, sL2, sR2);
    }

    // roofit variables
    // first CB function
    TString VM;
    Float_t meanVal, meanLow, meanUpp;
    if(iPcs<4 || iPcs==6) { VM = "J/#psi", meanVal = 3.097; meanLow = 2.8; meanUpp = 3.3; }
    else                  { VM = "#psi'",  meanVal = 3.686; meanLow = 3.4; meanUpp = 3.9; }
    RooRealVar _aL(Form("#alpha_{L,%s}",VM.Data()),Form("#alpha_{L,%s}",VM.Data()),aL,0.001,40.);
    RooRealVar _aR(Form("#alpha_{R,%s}",VM.Data()),Form("#alpha_{R,%s}",VM.Data()),aR,0.001,40.);   
    RooRealVar _mu(Form("#it{m}_{%s}",VM.Data()),Form("#it{m}_{%s}",VM.Data()),meanVal,meanLow,meanUpp);
    RooRealVar _nL(Form("#it{n}_{L,%s}",VM.Data()),Form("#it{n}_{L,%s}",VM.Data()),nL,0.001,100.);
    RooRealVar _nR(Form("#it{n}_{R,%s}",VM.Data()),Form("#it{n}_{R,%s}",VM.Data()),nR,0.001,100.);
    RooRealVar _sL(Form("#sigma_{L,%s}",VM.Data()),Form("#sigma_{L,%s}",VM.Data()),sL, 0.01, 0.50);
    RooRealVar _sR(Form("#sigma_{R,%s}",VM.Data()),Form("#sigma_{R,%s}",VM.Data()),sR, 0.01, 0.50);
    RooCrystalBall DSCB("DSCB", "Double sided CB function", fM, _mu, _sL, _sR, _aL, _nL, _aR, _nR);
    
    // second CB function (for iPcs == 6)
    TString VM2 = "#psi'";
    RooRealVar _aL2(Form("#alpha_{L,%s}",VM2.Data()),Form("#alpha_{L,%s}",VM2.Data()),aL2,0.01,40.);
    RooRealVar _aR2(Form("#alpha_{R,%s}",VM2.Data()),Form("#alpha_{R,%s}",VM2.Data()),aR2,0.01,40.);   
    RooRealVar _mu2(Form("#it{m}_{%s}",VM2.Data()),Form("#it{m}_{%s}",VM2.Data()),3.686,3.4,3.9);
    RooRealVar _nL2(Form("#it{n}_{L,%s}",VM2.Data()),Form("#it{n}_{L,%s}",VM2.Data()),nL2,0.01,100.);
    RooRealVar _nR2(Form("#it{n}_{R,%s}",VM2.Data()),Form("#it{n}_{R,%s}",VM2.Data()),nR2,0.01,100.);
    RooRealVar _sL2(Form("#sigma_{L,%s}",VM2.Data()),Form("#sigma_{L,%s}",VM2.Data()),sL2, 0.01, 0.20);
    RooRealVar _sR2(Form("#sigma_{R,%s}",VM2.Data()),Form("#sigma_{R,%s}",VM2.Data()),sR2, 0.01, 0.20);
    RooCrystalBall DSCB2("DSCB2", "Double sided CB function", fM, _mu2, _sL2, _sR2, _aL2, _nL2, _aR2, _nR2);

    // extended fit?
    Bool_t ext = kTRUE;
    RooRealVar _N("#it{N}","#it{N}",nEv,nEv*0.8,nEv*1.2);
    RooExtendPdf ExtDSCB("ExtDSCB","Extended Double sided CB function",DSCB,_N);

    // fit to the combined sample - PDF:
    RooRealVar _N1(Form("#it{N}_{%s}",VM.Data()),Form("#it{N}_{%s}",VM.Data()),nEv,nEv*0.8,nEv*1.2);
    RooRealVar _N2(Form("#it{N}_{%s}",VM2.Data()),Form("#it{N}_{%s}",VM2.Data()),nEv*0.01,nEv*0.001,nEv*0.1);
    RooAddPdf CombinedPDF("CombinedPDF","J/psi DSCB and psi' DSCB", RooArgList(DSCB,DSCB2),RooArgList(_N1,_N2));

    // perform the fit
    RooFitResult* fResFit;
    if(iPcs<6) {
        // fix some parameters?
        Bool_t fix = kFALSE;
        if(fix) {
            _aR.setConstant(kTRUE);
            _nR.setConstant(kTRUE);
        }
        if(ext) fResFit = ExtDSCB.fitTo(fDataHist,SumW2Error(kFALSE),Extended(kTRUE),Save());
        else    fResFit = DSCB.fitTo(fDataHist,SumW2Error(kFALSE),Save());
    } else {
        // fix the parameters of the right tails
        //_aL.setConstant(kTRUE); _nL.setConstant(kTRUE); _sL.setConstant(kTRUE); _sR.setConstant(kTRUE);
        _aR.setConstant(kTRUE); 
        _nR.setConstant(kTRUE);
        //_aL2.setConstant(kTRUE); _nL2.setConstant(kTRUE); _sL2.setConstant(kTRUE);  _sR2.setConstant(kTRUE);
        _aR2.setConstant(kTRUE); 
        _nR2.setConstant(kTRUE); 
        fResFit = CombinedPDF.fitTo(fDataHist,SumW2Error(kFALSE),Extended(kTRUE),Save());
    }

    // plot the results
    // draw correlation matrix
    TCanvas* cCM = new TCanvas("cCM","cCM",700,600);
    DrawCM(cCM,fResFit);

    // draw histogram and fit and calculate chi2
    TCanvas* c = new TCanvas("c","c",700,600);
    SetCanvas(c,kFALSE);
    gStyle->SetEndErrorSize(1); 

    RooPlot* fr = fM.frame(Title("inv mass fit")); 
    fDataHist.plotOn(fr,Name("fDataHist"),MarkerStyle(kFullCircle),MarkerSize(0.8),MarkerColor(kBlack),LineColor(kBlack),LineWidth(2));
    Float_t chi2;
    if(iPcs<6) {
        if(ext) {
            ExtDSCB.plotOn(fr,Name("ExtDSCB"),LineColor(215),LineWidth(3),LineStyle(9));
            chi2 = fr->chiSquare("ExtDSCB","fDataHist",fResFit->floatParsFinal().getSize());
        } else {
            DSCB.plotOn(fr,Name("DSCB"),LineColor(215),LineWidth(3),LineStyle(9));
            chi2 = fr->chiSquare("DSCB","fDataHist",fResFit->floatParsFinal().getSize());
        }
    } else {
        CombinedPDF.plotOn(fr,Name("DSCB"),Components(DSCB),LineColor(kRed),LineWidth(3),LineStyle(2));
        CombinedPDF.plotOn(fr,Name("DSCB2"),Components(DSCB2),LineColor(kViolet),LineWidth(3),LineStyle(2));
        CombinedPDF.plotOn(fr,Name("CombinedPDF"),LineColor(215),LineWidth(3),LineStyle(9));
        chi2 = fr->chiSquare("CombinedPDF","fDataHist",fResFit->floatParsFinal().getSize());
    }
    // y-axis
    fr->GetYaxis()->SetRangeUser(0.0,h->GetMaximum()*1.05);
    // title
    fr->GetYaxis()->SetTitle(Form("Counts per %.0f MeV/#it{c}^{2}", fMStep*1e3));
    fr->GetYaxis()->SetTitleSize(0.032);
    fr->GetYaxis()->SetTitleOffset(1.4);
    // label
    fr->GetYaxis()->SetLabelSize(0.032);
    fr->GetYaxis()->SetLabelOffset(0.01);
    fr->GetYaxis()->SetDecimals(1);
    fr->GetYaxis()->SetMaxDigits(3);
    // x-axis
    // title
    fr->GetXaxis()->SetTitle("#it{m}_{cl pair} [GeV/#it{c}^{2}]");
    fr->GetXaxis()->SetTitleSize(0.032);
    fr->GetXaxis()->SetTitleOffset(1.25);
    // label
    fr->GetXaxis()->SetLabelSize(0.032);
    fr->GetXaxis()->SetLabelOffset(0.01);
    fr->GetXaxis()->SetDecimals(1);
    fr->Draw();

    // print calculated chi2 
    Printf("********************");
    Printf("chi2/NDF = %.3f", chi2);
    Printf("NDF = %i", fResFit->floatParsFinal().getSize());
    Printf("chi2/NDF = %.3f/%i", chi2*fResFit->floatParsFinal().getSize(), fResFit->floatParsFinal().getSize());
    Printf("********************");
    
    TLegend* l = NULL;
    TLegend* l1 = NULL;
    if(iPcs<6) {
        Int_t nRows = 10;
        l = new TLegend(0.14,0.925-nRows*0.04,0.56,0.925);
        l->AddEntry("fDataHist","simulation","EPL");
        if(ext) l->AddEntry("ExtDSCB",Form("model (#chi^{2}/NDF = %.3f)",chi2),"L");
        else    l->AddEntry("DSCB",Form("model (#chi^{2}/NDF = %.3f)",chi2),"L");
        l->AddEntry((TObject*)0,Form("%s = %.f #pm %.f",_N.getTitle().Data(),_N.getVal(),_N.getError()),"");
        l->AddEntry((TObject*)0,Form("%s = %.2f GeV/#it{c}^{2}",_mu.getTitle().Data(),_mu.getVal()),"");
        l->AddEntry((TObject*)0,Form("%s = %.2f GeV/#it{c}^{2}",_sL.getTitle().Data(),_sL.getVal()),"");
        l->AddEntry((TObject*)0,Form("%s = %.2f GeV/#it{c}^{2}",_sR.getTitle().Data(),_sR.getVal()),""); 
        l->AddEntry((TObject*)0,Form("%s = %.1f",_aL.getTitle().Data(),_aL.getVal()),"");
        l->AddEntry((TObject*)0,Form("%s = %.1f",_aR.getTitle().Data(),_aR.getVal()),"");
        l->AddEntry((TObject*)0,Form("%s = %.1f",_nL.getTitle().Data(),_nL.getVal()),"");
        l->AddEntry((TObject*)0,Form("%s = %.1f",_nR.getTitle().Data(),_nR.getVal()),"");
        Int_t nRows1 = 2;
        l1 = new TLegend(0.14,0.925-nRows1*0.04,0.56,0.925);
        l1->AddEntry("fDataHist","simulation","EPL");
        if(ext) l1->AddEntry("ExtDSCB",Form("model (#chi^{2}/NDF = %.3f)",chi2),"L");
        else    l1->AddEntry("DSCB",Form("model (#chi^{2}/NDF = %.3f)",chi2),"L");
    } else {
        Int_t nRows = 12;
        l = new TLegend(0.14,0.925-nRows*0.04,0.56,0.925);
        l->AddEntry("fDataHist","simulation","EPL");
        l->AddEntry("CombinedPDF",Form("model (#chi^{2}/NDF = %.3f)",chi2),"L");
        l->AddEntry("DSCB","double-sided CB: J/#psi","L");
        l->AddEntry((TObject*)0,Form("%s = %.f #pm %.f",_N1.getTitle().Data(),_N1.getVal(),_N1.getError()),"");
        l->AddEntry((TObject*)0,Form("%s = %.2f GeV/#it{c}^{2}",_mu.getTitle().Data(),_mu.getVal()),"");
        l->AddEntry((TObject*)0,Form("%s = %.2f GeV/#it{c}^{2}",_sL.getTitle().Data(),_sL.getVal()),"");
        l->AddEntry((TObject*)0,Form("%s = %.2f GeV/#it{c}^{2}",_sR.getTitle().Data(),_sR.getVal()),""); 
        l->AddEntry("DSCB2","double-sided CB: #psi'","L");
        l->AddEntry((TObject*)0,Form("%s = %.f #pm %.f",_N2.getTitle().Data(),_N2.getVal(),_N2.getError()),"");
        l->AddEntry((TObject*)0,Form("%s = %.2f GeV/#it{c}^{2}",_mu2.getTitle().Data(),_mu2.getVal()),"");
        l->AddEntry((TObject*)0,Form("%s = %.2f GeV/#it{c}^{2}",_sL2.getTitle().Data(),_sL2.getVal()),"");
        l->AddEntry((TObject*)0,Form("%s = %.2f GeV/#it{c}^{2}",_sR2.getTitle().Data(),_sR2.getVal()),"");
        Int_t nRows1 = 4;
        l1 = new TLegend(0.14,0.925-nRows1*0.04,0.56,0.925);
        l1->AddEntry("fDataHist","simulation","EPL");
        l1->AddEntry("CombinedPDF",Form("model (#chi^{2}/NDF = %.3f)",chi2),"L");
        l1->AddEntry("DSCB","double-sided CB: J/#psi","L");
        l1->AddEntry("DSCB2","double-sided CB: #psi'","L");
    }
    l->SetMargin(0.16); 
    l->SetTextSize(0.032);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->Draw();
    l1->SetMargin(0.16); 
    l1->SetTextSize(0.032);
    l1->SetBorderSize(0);
    l1->SetFillStyle(0);
    // second legend
    Int_t nRows2 = 2;
    TLegend l2(0.72,0.925-nRows2*0.04,0.92,0.925);
    l2.AddEntry((TObject*)0,"J/#psi and #psi' #rightarrow e^{+}e^{-}","");
    l2.AddEntry((TObject*)0,Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}",fPtLow,fPtUpp),"");
    l2.SetTextSize(0.032);
    l2.SetBorderSize(0);
    l2.SetFillStyle(0);
    l2.SetMargin(0.);
    l2.Draw();
    // title
    TLatex* ltx = new TLatex();
    ltx->SetTextSize(0.032);
    ltx->SetTextAlign(21);
    ltx->SetNDC();
    ltx->DrawLatex(0.55,0.96,"ALICE Run-4 Simulation: Pb#minusPb UPC at #sqrt{#it{s}_{NN}} = 5.02 TeV");
    // save the plots
    cCM->Print(Form("%sfit_%s_CM.pdf",sOutMass.Data(),processes[iPcs].Data()));
    c->Print(Form("%sfit_%s.pdf",sOutMass.Data(),processes[iPcs].Data()));

    // save the values of the parameters
    ofstream of(Form("%sparam_%s.txt",sOutMass.Data(),processes[iPcs].Data()));
    of << std::fixed << std::setprecision(3)
       << _aL.getVal() << "\t" << _aL.getError() << "\n"
       << _aR.getVal() << "\t" << _aR.getError() << "\n"
       << _mu.getVal() << "\t" << _mu.getError() << "\n"
       << _nL.getVal() << "\t" << _nL.getError() << "\n"
       << _nR.getVal() << "\t" << _nR.getError() << "\n"
       << _sL.getVal() << "\t" << _sL.getError() << "\n"
       << _sR.getVal() << "\t" << _sR.getError() << "\n";
    if(ext) of << _N.getVal() << "\t" << _N.getError() << "\n";
    of.close();

    // draw histogram with log scale
    TCanvas *cLog = new TCanvas("cLog","cLog",700,600);
    SetCanvas(cLog,kTRUE);
    fr->GetYaxis()->SetRangeUser(0.1,h->GetMaximum()*2.5);
    fr->Draw();
    l1->Draw();
    l2.Draw();
    ltx->DrawLatex(0.55,0.96,"ALICE Run-4 Simulation: Pb#minusPb UPC at #sqrt{#it{s}_{NN}} = 5.02 TeV");
    cLog->Print(Form("%sfit_%s_log.pdf",sOutMass.Data(),processes[iPcs].Data()));
    
    delete cCM;
    delete c;
    delete cLog;
    return;
}

void SetHistoProperties(TH1F* h, Color_t c, Int_t lineStyle)
{
    h->GetYaxis()->SetRangeUser(1.,h->GetMaximum()*1.05);
    h->SetLineWidth(3);
    h->SetLineColor(c);
    h->SetLineStyle(lineStyle);
    return;
}

void PlotPtDistribution(TH1F* hData, TH1F* hCohJ, TH1F* hIncJ, TH1F* hCohFD, TH1F* hIncFD)
{
    TCanvas* c = new TCanvas("c","c",700,600);
    SetCanvas(c,kFALSE);
    // plot data
    // y-axis
    hData->GetYaxis()->SetRangeUser(1.,hData->GetMaximum()*1.05);
    // title
    hData->GetYaxis()->SetTitle(Form("Counts per %.0f MeV/#it{c}", fPtStep*1e3));
    hData->GetYaxis()->SetTitleSize(0.032);
    hData->GetYaxis()->SetTitleOffset(1.4);
    // label
    hData->GetYaxis()->SetLabelSize(0.032);
    hData->GetYaxis()->SetLabelOffset(0.01);
    hData->GetYaxis()->SetDecimals(1);
    hData->GetYaxis()->SetMaxDigits(3);
    // x-axis
    // title
    hData->GetXaxis()->SetTitle("#it{p}_{T,cl pair} [GeV/#it{c}]");
    hData->GetXaxis()->SetTitleSize(0.032);
    hData->GetXaxis()->SetTitleOffset(1.25);
    // label
    hData->GetXaxis()->SetLabelSize(0.032);
    hData->GetXaxis()->SetLabelOffset(0.01);
    hData->GetXaxis()->SetDecimals(1);
    // marker style
    hData->SetMarkerStyle(kFullCircle);
    hData->SetMarkerSize(0.7);
    hData->SetMarkerColor(kBlack);
    hData->SetLineColor(kBlack);
    hData->SetLineWidth(2);
    hData->Draw("E1");
    // coh J/psi
    SetHistoProperties(hCohJ,215,1);
    hCohJ->Draw("HIST SAME");
    // inc J/psi
    SetHistoProperties(hIncJ,kRed,1);
    hIncJ->Draw("HIST SAME");
    // coh FD
    SetHistoProperties(hCohFD,215,2);
    hCohFD->Draw("HIST SAME");
    // inc FD
    SetHistoProperties(hIncFD,kRed,2);
    hIncFD->Draw("HIST SAME");
    // title
    TLatex* ltx = new TLatex();
    ltx->SetTextSize(0.032);
    ltx->SetTextAlign(21);
    ltx->SetNDC();
    ltx->DrawLatex(0.55,0.96,"ALICE Run-4 Simulation: Pb#minusPb UPC at #sqrt{#it{s}_{NN}} = 5.02 TeV");
    // legends
    Int_t nRows = 7;
    TLegend l(0.58,0.925-nRows*0.04,1.00,0.925);
    l.AddEntry((TObject*)0,"J/#psi and #psi' photoproduction","");
    l.AddEntry((TObject*)0,Form("%.1f < #it{m}_{ee} < %.1f GeV/#it{c}^{2}",fMLow,fMUpp),"");
    l.AddEntry(hData,"simulation","EPL");
    l.AddEntry(hCohJ,"coh J/#psi","L");
    l.AddEntry(hIncJ,"inc J/#psi","L");
    l.AddEntry(hCohFD,"feed-down from coh #psi'","L");
    l.AddEntry(hIncFD,"feed-down from inc #psi'","L");
    l.SetTextSize(0.032);
    l.SetBorderSize(0);
    l.SetFillStyle(0);
    l.SetMargin(0.16);
    l.Draw();
    c->Print(Form("%sptDist.pdf",sOutPt.Data()));
    // log scale
    TCanvas* cLog = new TCanvas("cLog","cLog",700,600);
    SetCanvas(cLog,kTRUE);
    hData->GetYaxis()->SetRangeUser(1.,hData->GetMaximum()*2.0);
    hData->Draw("E1");
    hCohJ->Draw("HIST SAME");
    hIncJ->Draw("HIST SAME");
    hCohFD->Draw("HIST SAME");
    hIncFD->Draw("HIST SAME");
    ltx->DrawLatex(0.55,0.96,"ALICE Run-4 Simulation: Pb#minusPb UPC at #sqrt{#it{s}_{NN}} = 5.02 TeV");
    l.Draw();
    cLog->Print(Form("%sptDist_log.pdf",sOutPt.Data()));
    return;
}

void AnaMain_SignalExtraction()
{
    // set the binning
    nBinsM = (fMEdgeUpp-fMEdgeLow) / fMStep;
    Printf("%i mass bins will be used.", nBinsM);
    nBinsPt = (fPtEdgeUpp-fPtEdgeLow) / fPtStep;
    Printf("%i pT bins will be used.", nBinsPt);

    // set the paths
    outSubDir = CreateOutputSubDir();
    sOutMass = Form("results/sim%02i_g%02i_p%02i/invMassFit/",simFiles,fileGeom,filePara);
    gSystem->Exec("mkdir -p " + sOutMass);
    sOutPt = Form("results/sim%02i_g%02i_p%02i/ptDist/",simFiles,fileGeom,filePara);
    gSystem->Exec("mkdir -p " + sOutPt);

    // fill inv mass and pT histograms
    TH1F* hMass[7] = { NULL };
    TH1F* hPt[7] = { NULL };
    TH1F* hMassRnd[7] = { NULL };
    TH1F* hPtRnd[7] = { NULL };
    // define the histograms
    for(Int_t i = 0; i < 7; i++) {
        hMass[i] = new TH1F(Form("hMass_%s",processes[i].Data()),names[i],nBinsM,fMEdgeLow,fMEdgeUpp);
        hPt[i] = new TH1F(Form("hPt_%s",processes[i].Data()),names[i],nBinsPt,fPtEdgeLow,fPtEdgeUpp);
        hMassRnd[i] = new TH1F(Form("hMassRnd_%s",processes[i].Data()),names[i],nBinsM,fMEdgeLow,fMEdgeUpp);
        hPtRnd[i] = new TH1F(Form("hPtRnd_%s",processes[i].Data()),names[i],nBinsPt,fPtEdgeLow,fPtEdgeUpp);       
    }
    // fill the histograms
    for(Int_t i = 0; i < 6; i++) {
        PrepareHistos(i,hMass[i],hMassRnd[i],hPt[i],hPtRnd[i]);
        hMass[6]->Add(hMass[i]);
        hPt[6]->Add(hPt[i]);
        hMassRnd[6]->Add(hMassRnd[i]);
        hPtRnd[6]->Add(hPtRnd[i]);
    }
    PrintHistogram(hMass[6],sOutMass);
    PrintHistogram(hPt[6],sOutPt);
    PrintHistogram(hMassRnd[6],sOutMass);
    PrintHistogram(hPtRnd[6],sOutPt);

    // fit coherent J/psi
    DoInvMassFit(0,hMass[0]);

    // fit coherent psi'
    DoInvMassFit(4,hMass[4]);

    // fit the combined sample
    DoInvMassFit(6,hMass[6]);

    // plot the pT distribution with highlighted components
    PlotPtDistribution(hPt[6],hPt[0],hPt[1],hPt[2],hPt[3]);

    return;
}