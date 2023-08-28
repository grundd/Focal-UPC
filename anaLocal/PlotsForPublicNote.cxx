// PlotsForPublicNote.cxx
// David Grund, Aug 16, 2023

int npixx = 900;
int npixy = 800;
float textsize = 0.05;
float textsizeleg = 0.038;
bool swappedAxes = true;

template <typename TH>
void setHistoAxes (TH* h)
{
    // histo axes settings
    h->SetBit(TH1::kNoStats);
    h->SetBit(TH1::kNoTitle);
    // x
    h->GetXaxis()->SetLabelSize(textsize);
    h->GetXaxis()->SetTitleSize(textsize);
    h->GetXaxis()->SetTitleOffset(1.05);
    // y
    h->GetYaxis()->SetLabelSize(textsize);
    h->GetYaxis()->SetTitleSize(textsize);
    h->GetYaxis()->SetTitleOffset(1.25);
    h->GetYaxis()->SetMaxDigits(3);
    return;
}

TCanvas* plotHisto1D (TH1F* h, TH1F* h1 = NULL)
{
    TCanvas* c = new TCanvas("c","",npixx,npixy);
    c->SetTopMargin(0.06);
    c->SetBottomMargin(0.12);
    c->SetRightMargin(0.04);
    c->SetLeftMargin(0.14);
    c->cd();
    // draw histo(s)
    setHistoAxes(h);
    h->SetLineColor(kBlue);
    h->SetLineWidth(2);
    h->Draw("HIST"); // or HIST E
    if(h1) {
        h1->SetLineColor(kRed);
        h1->SetLineWidth(2);
        h1->Draw("HIST SAME"); // or HIST E SAME
    }
    return c;
}

TCanvas* plotHisto2D (TH2F* h)
{
    TCanvas* c = new TCanvas("c","",npixx,npixy);
    c->SetTopMargin(0.06);
    c->SetBottomMargin(0.12);
    c->SetRightMargin(0.14);
    c->SetLeftMargin(0.14);
    c->cd();
    // draw histo:
    setHistoAxes(h);
    h->GetZaxis()->SetLabelSize(textsize);
    h->GetZaxis()->SetMaxDigits(3);
    h->Draw("COLZ");
    return c;
}

void swapHistos (TH2F* hOld, TH2F* hNew)
{
    for(int i = 1; i <= 100; i++) {
        for(int j = 1; j <= 100; j++) {
            float binContent = hOld->GetBinContent(i, j);
            hNew->SetBinContent(j, i, binContent);            
        }
    }
    return;
}

void PlotsForPublicNote () 
{
    gSystem->Exec("mkdir -p publicNotePlots/");
    // coh J/psi plots
    // grid analysis: no superclusterizer
    // grid analysis
    TFile *fGrd = TFile::Open("results/sim05_g05_p05/cohJpsi/merged_output/analysisResultsGrid.root","read");
    TList *lGrd_TH1F = (TList*) fGrd->Get("lTH1F");
    TList *lGrd_TH2F = (TList*) fGrd->Get("lTH2F");
    // grid analysis: superclusterizer
    TFile *fGrd_supcl= TFile::Open("results/sim05_g05_p05/cohJpsi/merged_output_supCl_cutE10.0/analysisResultsGrid.root","read");
    TList *lGrd_supcl_TH1F = (TList*) fGrd_supcl->Get("lTH1F");
    TList *lGrd_supcl_TH2F = (TList*) fGrd_supcl->Get("lTH2F");
    // main analysis: no superclusterizer
    TFile *fMn = TFile::Open("results/sim05_g05_p05/cohJpsi/merged_output/analysisResultsMain.root","read");
    TList *lMn_TH1F = (TList*) fMn->Get("lTH1F");
    TList *lMn_TH2F = (TList*) fMn->Get("lTH2F");
    // main analysis: superclusterizer
    TFile *fMn_supcl = TFile::Open("results/sim05_g05_p05/cohJpsi/merged_output_supCl_cutE10.0/analysisResultsMain.root","read");
    TList *lMn_supcl_TH1F = (TList*) fMn_supcl->Get("lTH1F");
    TList *lMn_supcl_TH2F = (TList*) fMn_supcl->Get("lTH2F");

    // pt dist of coh J/psi
    TH1F* hPt = (TH1F*)lGrd_supcl_TH1F->FindObject("hJ1_mcJElPairPt");
    hPt->SetTitle(";#it{p}_{T} (GeV/#it{c});counts");
    hPt->GetXaxis()->SetRangeUser(0.,0.6);
    TCanvas* cPt = plotHisto1D(hPt);
    TLegend* lGen = new TLegend(0.32,0.795,1.,0.93);
    lGen->AddEntry((TObject*)0,"ALICE Simulation, Pb#minusPb UPC #sqrt{#it{s}_{NN}} = 5.5 TeV","");
    lGen->AddEntry((TObject*)0,"STARlight coherent J/#psi #rightarrow e^{+}e^{-}","");
    lGen->AddEntry((TObject*)0,"3.4 < #it{y} < 6.0","");
    lGen->SetTextSize(textsizeleg);
    lGen->SetBorderSize(0);
    lGen->SetFillStyle(0);
    lGen->SetMargin(0.);
    lGen->Draw();
    cPt->Print("publicNotePlots/cohJpsi_genPt.pdf");
    cPt->Print("publicNotePlots/cohJpsi_genPt.C");
    delete cPt;

    // rap dist of coh J/psi
    TH1F* hRap = (TH1F*)lGrd_supcl_TH1F->FindObject("hJ1_mcJElPairRap");
    hRap->SetTitle(";#it{y} (-);counts");
    hRap->GetXaxis()->SetRangeUser(3.3,6.0);
    TCanvas* cRap = plotHisto1D(hRap);
    lGen->Draw();
    cRap->Print("publicNotePlots/cohJpsi_genRap.pdf");
    cRap->Print("publicNotePlots/cohJpsi_genRap.C");
    delete cRap;
    delete lGen;

    /*
    // rad separation of clusters
    TH1F* hSep = (TH1F*)lMn_TH1F->FindObject("hJ1_clPairSep");
    TH1F* hSep_supcl = (TH1F*)lMn_supcl_TH1F->FindObject("hJ1_clPairSep");
    TCanvas* cSep = plotHisto1D(hSep,hSep_supcl);
    //cSep->SetLogy();
    cSep->Print("publicNotePlots/cohJpsi_sep.pdf");
    cSep->Print("publicNotePlots/cohJpsi_sep.C");
    delete cSep;

    // correlation: radial sep
    TH2F* hSepCorr = (TH2F*)lMn_supcl_TH2F->FindObject("hJ2_clPairSep_mcJElSep");
    hSepCorr->SetTitle(";#Delta#it{R}_{sup-cl pairs} (cm);#Delta#it{R}_{electron pairs} (cm)");
    TCanvas* cSepCorr = plotHisto2D(hSepCorr);
    cSepCorr->Print("publicNotePlots/cohJpsi_sepCorr.pdf");
    cSepCorr->Print("publicNotePlots/cohJpsi_sepCorr.C");
    delete cSepCorr;
    */

    // correlation: energy of sup-cls
    TH2F* hEneCorr = (TH2F*)lGrd_supcl_TH2F->FindObject("hJ2_ppeClEn_mtchEn");
    hEneCorr->SetTitle(";#it{E}_{cl} (GeV);#it{E}_{matched prim el} (GeV)");
    hEneCorr->GetYaxis()->SetNdivisions(508,kFALSE);
    hEneCorr->GetXaxis()->SetNdivisions(508,kFALSE);
    // create a histogram with swapped axes
    TH2F* hEneCorrSwp = new TH2F("hEneCorrSwp","",100,0.,200.,100,0.,200.);
    swapHistos(hEneCorr,hEneCorrSwp);
    hEneCorrSwp->SetTitle(";#it{E}_{matched prim el} (GeV);#it{E}_{cl} (GeV)");
    hEneCorrSwp->GetYaxis()->SetNdivisions(508,kFALSE);
    hEneCorrSwp->GetXaxis()->SetNdivisions(508,kFALSE);
    // print it
    TCanvas* cEneCorr = NULL;
    if(swappedAxes) cEneCorr = plotHisto2D(hEneCorrSwp);
    else            cEneCorr = plotHisto2D(hEneCorr);
    cEneCorr->SetLogz();
    int nRows = 2;
    TLegend* lCorr = new TLegend(0.18,0.930-nRows*0.045,0.80,0.93);
    lCorr->AddEntry((TObject*)0,"ALICE Simulation, Pb#minusPb UPC #sqrt{#it{s}_{NN}} = 5.5 TeV","");
    lCorr->AddEntry((TObject*)0,"STARlight coherent J/#psi #rightarrow e^{+}e^{-}","");
    lCorr->SetTextSize(textsizeleg);
    lCorr->SetBorderSize(0);
    lCorr->SetFillStyle(0);
    lCorr->SetMargin(0.);
    lCorr->Draw();
    cEneCorr->Print("publicNotePlots/cohJpsi_energyCorr.pdf");
    cEneCorr->Print("publicNotePlots/cohJpsi_energyCorr.C");
    delete cEneCorr;
    delete lCorr;

    // correlation: energy of sup-cl pairs
    TH2F* hPairEneCorr = (TH2F*)lMn_supcl_TH2F->FindObject("hJ2_ppeClPairEn_mtchEn");
    hPairEneCorr->SetTitle(";#it{E}_{cl pair} (GeV);#it{E}_{matched prim el pair} (GeV)");
    hPairEneCorr->GetYaxis()->SetNdivisions(508,kFALSE);
    hPairEneCorr->GetXaxis()->SetNdivisions(508,kFALSE);
    // create a histogram with swapped axes
    TH2F* hPairEneCorrSwp = new TH2F("hPairEneCorrSwp","",100,0.,200.,100,0.,200.);
    swapHistos(hPairEneCorr,hPairEneCorrSwp);
    hPairEneCorrSwp->SetTitle(";#it{E}_{matched prim el pair} (GeV);#it{E}_{cl pair} (GeV)");
    hPairEneCorrSwp->GetYaxis()->SetNdivisions(508,kFALSE);
    hPairEneCorrSwp->GetXaxis()->SetNdivisions(508,kFALSE);
    // print it
    TCanvas* cPairEneCorr = NULL;
    if(swappedAxes) cPairEneCorr = plotHisto2D(hPairEneCorrSwp);
    else            cPairEneCorr = plotHisto2D(hPairEneCorr);
    cPairEneCorr->SetLogz();
    nRows = 3;
    TLegend* lCorrPair = new TLegend(0.18,0.93-nRows*0.045,0.56,0.93);
    lCorrPair->AddEntry((TObject*)0,"ALICE Simulation","");
    lCorrPair->AddEntry((TObject*)0,"Pb#minusPb UPC #sqrt{#it{s}_{NN}} = 5.5 TeV","");
    lCorrPair->AddEntry((TObject*)0,"STARlight coh. J/#psi #rightarrow e^{+}e^{-}","");
    lCorrPair->SetTextSize(textsizeleg);
    lCorrPair->SetBorderSize(0);
    //lCorrPair->SetFillStyle(0);
    //lCorrPair->SetFillColor(kRed);
    lCorrPair->SetMargin(0.);
    lCorrPair->Draw();
    cPairEneCorr->Print("publicNotePlots/cohJpsi_energyCorrPair.pdf");
    cPairEneCorr->Print("publicNotePlots/cohJpsi_energyCorrPair.C");
    delete cPairEneCorr;
    delete lCorrPair;

    return;
}