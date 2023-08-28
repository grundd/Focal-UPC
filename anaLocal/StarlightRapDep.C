// StarlightRapDep.C
// David Grund, Nov 02, 2022

// cpp headers
#include <fstream>
#include <vector>
// root headers
#include "TSystem.h"
#include "TList.h"
#include "TH1.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TLegend.h"
// my headers
#include "StarlightRapDep.h"

Int_t nRapBins;
Float_t fRapLow = -6.; // [-]
Float_t fRapUpp = +6.;
Float_t fRapStep = 0.2; 
Int_t binNegLow = 2;
Int_t binNegUpp = 13;
Int_t binPosLow = 48;
Int_t binPosUpp = 59;
Float_t fLumiRun4 = 7.0; // [nb^(-1)]
// direct starlight
Float_t fTotalCS[7] = {
    38.8, // [mb]
    17.8,
    7.5,
    3.1,
    0.094,
    0.041,
    83.1
}; 
TString SL_label[7] = {
    "coh J/#psi #rightarrow e^{+}e^{-}",
    "inc J/#psi #rightarrow e^{+}e^{-}",
    "coh #psi' #rightarrow e^{+}e^{-}",
    "inc #psi' #rightarrow e^{+}e^{-}",
    "coh #Upsilon(1S) #rightarrow e^{+}e^{-}",
    "inc #Upsilon(1S) #rightarrow e^{+}e^{-}",
    "medium-mass continuum (1.8 < #it{W} < 15 GeV): #gamma#gamma #rightarrow e^{+}e^{-}",
};
TString SL_folder[7] = {
    "CohJ",
    "IncJ",
    "CohP",
    "IncP",
    "CohU",
    "IncU",
    "Bkg"
};
Float_t fBR[7] = {
    0.0594, // J/psi to e^(+)e^(-)
    0.0594, 
    0.00772, // psi' to e^(+)e^(-)
    0.00772, 
    0.0238, // Upsilon to e^(+)e^(-)
    0.0238,
    1.0, // gamma-gamma -> no BR
};
// feed-down:
TString FD_label[2] = {
    "coh #psi' #rightarrow J/#psi + #pi^{+}#pi^{-} #rightarrow e^{+}e^{-} + #pi^{+}#pi^{-}",
    "inc #psi' #rightarrow J/#psi + #pi^{+}#pi^{-} #rightarrow e^{+}e^{-} + #pi^{+}#pi^{-}"
};
TString FD_folder[2] = {
    "CohFD",
    "IncFD"
};

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

TCanvas* CreateCanvas()
{
    TCanvas* c = new TCanvas("c","c",700,600);
    c->SetTopMargin(0.07);
    c->SetBottomMargin(0.1);
    c->SetRightMargin(0.02);
    c->SetLeftMargin(0.12);
    return c;
}

void CalculateRapDep(Int_t opt)
{
    Printf("\n*** Process: %s:", SL_label[opt].Data());

    TH1F* hRap_all = NULL;
    TH1F* hRap_twoEl = NULL;
    TH1F* hRap_accFo = NULL;
    TFile* fOut = TFile::Open("results/starlightRapDep/h" + SL_folder[opt] + ".root","read");
    if(fOut)
    {
        Printf("Histograms already exist and will be loaded.");
        TList* l = (TList*) fOut->Get("HistList");
        if(l) Printf("%s loaded.", l->GetName()); 
        hRap_all = (TH1F*) l->FindObject("hRap_all");
        if(hRap_all) Printf("%s loaded.", hRap_all->GetName());
        hRap_twoEl = (TH1F*) l->FindObject("hRap_twoEl");
        if(hRap_twoEl) Printf("%s loaded.", hRap_twoEl->GetName());
        hRap_accFo = (TH1F*) l->FindObject("hRap_accFo");
        if(hRap_accFo) Printf("%s loaded.", hRap_accFo->GetName());
    } 
    else 
    {
        Printf("Histograms will be created.");

        TFile* fSL = TFile::Open("inputData/starlight/" + SL_folder[opt] + "/tSTARlight.root", "read");
        if(fSL) Printf("File %s loaded.", fSL->GetName());

        TTree* tSL = dynamic_cast<TTree*> (fSL->Get("starlightTree"));
        if(tSL) Printf("Tree %s loaded.", tSL->GetName());
        SetBranchAddresses_tSL(tSL);

        Int_t nEntries = tSL->GetEntries();
        Printf("tSTARlight contains %i entries.", nEntries);

        gROOT->cd();
        TList* l = new TList();
        hRap_all = new TH1F("hRap_all","rap. dist.: gen. VM",nRapBins,fRapLow,fRapUpp);
        hRap_twoEl = new TH1F("hRap_twoEl","rap. dist.: gen. VM with two electrons",nRapBins,fRapLow,fRapUpp);
        hRap_accFo = new TH1F("hRap_accFo","rap. dist.: gen. VM with both electrons in FoCal acceptance",nRapBins,fRapLow,fRapUpp);
        l->Add(hRap_all);
        l->Add(hRap_twoEl);
        l->Add(hRap_accFo);

        Float_t progress = 0.; // perc
        for(Int_t iEntry = 0; iEntry < nEntries; iEntry++)
        {
            // update the progress bar
            if((iEntry+1) % (Int_t)(nEntries/10.) == 0) {
                progress += 10.;
                cout << "[" << progress << "%] done." << endl;
            }
            tSL->GetEntry(iEntry);
            if(daughters->GetEntries() == 2)
            {
                TLorentzVector* el1 = dynamic_cast<TLorentzVector*>(daughters->At(0));
                TLorentzVector* el2 = dynamic_cast<TLorentzVector*>(daughters->At(1));
                if(el1 && el2) {
                    hRap_twoEl->Fill(parent->Rapidity());
                    Float_t EtaEl1 = el1->Eta();
                    Float_t EtaEl2 = el2->Eta();
                    if(3.4 < EtaEl1 && EtaEl1 < 5.8 && 3.4 < EtaEl2 && EtaEl2 < 5.8) hRap_accFo->Fill(parent->Rapidity());
                }
            }
            hRap_all->Fill(parent->Rapidity());
        }
        fOut = new TFile("results/starlightRapDep/h" + SL_folder[opt] + ".root","RECREATE");
        l->Write("HistList", TObject::kSingleKey);
        l->ls();
        fOut->ls();
        fSL->Close();
    }

    // scale the histogram and create the histogram showing the cross section
    TH1F* hRap_CS = (TH1F*)hRap_all->Clone("hRap_CS");
    hRap_CS->Scale(fTotalCS[opt]/hRap_all->Integral(),"width");
    // create the graph
    TGraph* gr = new TGraph(nRapBins);
    gr->SetTitle(";#it{y} [-]; d#sigma/d#it{y} [mb]");
    // graph properties
    gr->SetLineStyle(9);
    gr->SetLineColor(kBlue);
    gr->SetLineWidth(3);
    // fill the graph
    for(Int_t i = 0; i < nRapBins; i++) {
        gr->SetPoint(i,hRap_CS->GetBinCenter(i+1),hRap_CS->GetBinContent(i+1));
    }
    // axis ranges
    TAxis *ax = gr->GetXaxis();
    ax->SetLimits(fRapLow-0.5,fRapUpp+0.5); 
    TH1F* h = gr->GetHistogram();
    h->SetMaximum(h->GetMaximum() * 1.25);      
    h->SetMinimum(0.); 
    h->GetXaxis()->SetTitleSize(0.042);
    h->GetXaxis()->SetLabelSize(0.042);
    h->GetYaxis()->SetTitleSize(0.042);
    h->GetYaxis()->SetLabelSize(0.042);
    h->GetYaxis()->SetDecimals(1);
    h->GetYaxis()->SetMaxDigits(3);
    // plot it
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TCanvas* c = CreateCanvas();
    gr->Draw("AC");

    TLatex* latex = new TLatex();
    latex->SetTextSize(0.036);
    latex->SetTextAlign(21);
    latex->SetNDC();
    latex->DrawLatex(0.55,0.95,Form("STARlight: %s",SL_label[opt].Data()));

    Printf("Low edge of bin %i: %.1f", binNegLow, hRap_CS->GetBinLowEdge(binNegLow));
    Printf("Upp edge of bin %i: %.1f", binNegUpp, hRap_CS->GetBinLowEdge(binNegUpp+1));
    Printf("Low edge of bin %i: %.1f", binPosLow, hRap_CS->GetBinLowEdge(binPosLow));
    Printf("Upp edge of bin %i: %.1f", binPosUpp, hRap_CS->GetBinLowEdge(binPosUpp+1));
    Float_t sigmaTotal = hRap_CS->Integral("width");
    Float_t sigmaForwd = hRap_CS->Integral(binNegLow,binNegUpp,"width") + hRap_CS->Integral(binPosLow,binPosUpp,"width");
    Float_t sigmaFocal = hRap_CS->Integral(binPosLow,binPosUpp,"width");

    // calculate expected yields of J/psi, psi' and Y in FoCal in Run 4
    // convert lumi from nb^(-1) to mb^(-1)
    Float_t N1 = hRap_all->Integral(); // # of generated events
    Float_t N2 = hRap_twoEl->Integral(); // # of events with two electrons in daughter particles
    Float_t N3 = hRap_all->Integral(binPosLow,binPosUpp); // # of events with VM rapidity from 3.4 to 5.8
    Float_t N4 = hRap_accFo->Integral(); // # of events with eta of both electrons from 3.4 to 5.8
    Float_t lumi = fLumiRun4 * 1e6; // mb^(-1)
    Float_t N_simulate = lumi * sigmaTotal * fBR[opt] * N3 / N2;
    Float_t N_yieldFoc = lumi * sigmaTotal * fBR[opt] * N4 / N2;
    ofstream of;
    of.open("results/starlightRapDep/log" + SL_folder[opt] + ".txt");
    of << Form("rapidity dependence of the %s cross section:", SL_label[opt].Data())
       << "*\n"
       << Form(" N1: # of gen ev: %.0f\n", N1)
       << Form(" N2: # of gen ev containing two daugter electrons: %.0f\n", N2)
       << Form(" N3: # of gen ev with J/psi rapidity from 3.4 to 5.8: %.0f\n", N3)
       << Form(" N4: # of gen ev with eta of both ele within 3.4 to 5.8: %.0f\n", N4)
       << "*\n"
       << "cross section values for respective rapidity intervals:\n"
       << Form(" |y| < 6.0: %.2e mb\n", sigmaTotal)
       << Form(" 3.4 < |y| < 5.8: %.2e mb\n", sigmaForwd)
       << Form(" 3.4 < y < 5.8: %.2e mb\n", sigmaFocal)
       << "*\n"
       << Form(" acceptance wrt |y| < 6.0 (N4 / N2): %.3f\n", N4 / N2)
       << Form(" acceptance wrt 3.4 < y < 5.8 (N4 / N3): %.3f\n", N4 / N3)
       << Form(" integrated lumi UPC Pb+Pb Run 4: %.2e mb^(-1)\n", lumi)
       << Form(" FOCAL cross section: %.2e mb\n", sigmaTotal * fBR[opt] * N4 / N2)
       << " | meaning: photopr. VM decays into ee both electrons reach FOCAL (eta of both within 3.4 to 5.8)\n"
       << " | equal to sigma(|y| < 6.0) * BR(VM -> ee) * Acc\n"
       << "*\n"
       << "expected yields:\n"
       << Form(" # of VM with 3.4 < y < 5.8: %.2e\n", N_simulate)
       << Form(" # of VM with both ele having 3.4 < eta < 5.8: %.2e\n", N_yieldFoc);
    of.close();

    // legends
    TLegend *l1 = new TLegend(0.15,0.73,0.55,0.93);
    l1->AddEntry((TObject*)0,Form("#sigma(|#it{y}| < 6.0) = %.3f mb", sigmaTotal),"");
    l1->AddEntry((TObject*)0,Form("#sigma(3.4 < |#it{y}| < 5.8) = %.3f mb", sigmaForwd),"");
    l1->AddEntry((TObject*)0,Form("#sigma(3.4 < #it{y} < 5.8) = %.3f mb", sigmaFocal),"");
    l1->AddEntry((TObject*)0,Form("Acc^{[1]}: %.2f %%", N4 / N3 * 1e2),"");
    l1->SetTextSize(0.036);
    l1->SetBorderSize(0);
    l1->SetFillStyle(0);
    l1->SetMargin(0.);
    l1->Draw();

    TLegend *l2 = new TLegend(0.58,0.73,0.97,0.93);
    l2->AddEntry((TObject*)0,Form("Run-4 expectations:"),"");
    l2->AddEntry((TObject*)0,Form("#it{L}_{int} (UPC Pb#minusPb): %.1f nb^{#minus1}", fLumiRun4),"");
    l2->AddEntry((TObject*)0,Form("N(3.4 < #it{y} < 5.8) #approx %.0f", roundFloat(N_simulate)),"");
    l2->AddEntry((TObject*)0,Form("N(3.4 < #it{#eta}^{e^{#pm}} < 5.8 #approx %.0f", roundFloat(N_yieldFoc)),"");
    l2->SetTextSize(0.036);
    l2->SetBorderSize(0);
    l2->SetFillStyle(0);
    l2->SetMargin(0.);
    l2->Draw();
    // save the canvas
    c->Print("results/starlightRapDep/gr" + SL_folder[opt] + ".pdf");
    delete c;
    delete gr;
    return;
}

void AnalyzeFeedDown(Int_t opt)
{
    Printf("\n*** Process: %s:", FD_label[opt].Data());

    TH1F* hRap_all = NULL;
    TH1F* hRap_twoEl = NULL;
    TH1F* hRap_accFo = NULL;
    TFile* fOut = TFile::Open("results/starlightRapDep/h" + FD_folder[opt] + ".root","read");
    if(fOut)
    {
        Printf("Histograms already exist and will be loaded.");
        TList* l = (TList*) fOut->Get("HistList");
        if(l) Printf("%s loaded.", l->GetName()); 
        hRap_all = (TH1F*) l->FindObject("hRap_all");
        if(hRap_all) Printf("%s loaded.", hRap_all->GetName());
        hRap_twoEl = (TH1F*) l->FindObject("hRap_twoEl");
        if(hRap_twoEl) Printf("%s loaded.", hRap_twoEl->GetName());
        hRap_accFo = (TH1F*) l->FindObject("hRap_accFo");
        if(hRap_accFo) Printf("%s loaded.", hRap_accFo->GetName());
    } 
    else 
    {
        Printf("Histograms will be created.");

        gROOT->cd();
        TList* l = new TList();
        hRap_all = new TH1F("hRap_all","rap. dist.: psi'",nRapBins,fRapLow,fRapUpp);
        hRap_twoEl = new TH1F("hRap_twoEl","rap. dist.: psi' with two electrons from J/psi decay",nRapBins,fRapLow,fRapUpp);
        hRap_accFo = new TH1F("hRap_accFo","rap. dist.: psi' with both electrons from J/psi decay in FoCal acceptance",nRapBins,fRapLow,fRapUpp);
        l->Add(hRap_all);
        l->Add(hRap_twoEl);
        l->Add(hRap_accFo);

        // go over galice files and fill the histograms
        for(Int_t iFile = 0; iFile < 100; iFile++)
        {
            TString galice = Form("inputData/starlight/%s/%03i/galice.root",FD_folder[opt].Data(),iFile+1);
            AliRunLoader* runLoader = AliRunLoader::Open(galice);
            runLoader->LoadKinematics();

            Int_t nEv = runLoader->GetNumberOfEvents();
            Printf("galice file %03i: %i entries.", iFile+1, nEv);

            for(Int_t iEv = 0; iEv < nEv; iEv++)
            {
                // get tracks in this event
                runLoader->GetEvent(iEv);	    
                AliStack* stack = runLoader->Stack();
                Int_t nTrks = stack->GetNtrack();
                // loop over tracks
                vector<Int_t> idxElTrk;
                Int_t idxPsi2s;
                for(Int_t iTrk = 0; iTrk < nTrks; iTrk++) 
                {
                    TParticle* part = stack->Particle(iTrk);
                    // if psi'
                    if(part->GetPdgCode() == 100443) idxPsi2s = iTrk;
                    // if an electron
                    if(TMath::Abs(part->GetPdgCode()) == 11) {
                        // if physical primary
                        if(stack->IsPhysicalPrimary(iTrk)) idxElTrk.push_back(iTrk);
                    } 
                }
                TParticle *psi2s = stack->Particle(idxPsi2s);
                if(idxElTrk.size() == 2)
                {
                    TParticle* el1 = stack->Particle(idxElTrk[0]);
                    TParticle* el2 = stack->Particle(idxElTrk[1]);
                    if(el1 && el2) {
                        hRap_twoEl->Fill(psi2s->Y());
                        Float_t EtaEl1 = el1->Eta();
                        Float_t EtaEl2 = el2->Eta();
                        if(3.4 < EtaEl1 && EtaEl1 < 5.8 && 3.4 < EtaEl2 && EtaEl2 < 5.8) hRap_accFo->Fill(psi2s->Y());
                    }
                }
                hRap_all->Fill(psi2s->Y());
            }
            runLoader->Delete();
        }
        fOut = new TFile("results/starlightRapDep/h" + FD_folder[opt] + ".root","RECREATE");
        l->Write("HistList", TObject::kSingleKey);
        l->ls();
        fOut->ls();
    }
    Printf("Low edge of bin %i: %.1f", binPosLow, hRap_all->GetBinLowEdge(binPosLow));
    Printf("Upp edge of bin %i: %.1f", binPosUpp, hRap_all->GetBinLowEdge(binPosUpp+1));
    // calculate the fraction of events with psi' with 3.4 < y < 5.8
    // where both electrons are within FoCal acceptance (3.4 < eta < 5.8)
    Float_t nTwoEl = hRap_twoEl->Integral(binPosLow,binPosUpp,"width");
    Float_t nAccFo = hRap_accFo->Integral(binPosLow,binPosUpp,"width");
    Float_t frac = nAccFo / nTwoEl;
        ofstream of;
    of.open("results/starlightRapDep/log" + FD_folder[opt] + ".txt");
    of << Form("%s:\n", FD_label[opt].Data())
       << "*\n"
       << Form(" N1: # of events with 3.4 < y(psi') < 5.8: %.0f\n", nTwoEl)
       << Form(" N2: # of (a & b) events where both J/psi electrons reached FoCal: %.0f\n", nAccFo)
       << " | a -> have 3.4 < y(psi') < 5.8\n"
       << " | b -> psi' decayed into J/psi + pi pi -> e e pi pi\n"
       << Form(" fraction (N2 / N1): %.3f\n", frac);
    of.close();

    // draw the histograms
    TCanvas* c = CreateCanvas();
    hRap_twoEl->SetTitle(";#it{y}(#psi') [-]; #it{N} [counts]");
    hRap_twoEl->SetLineColor(kBlue);
    hRap_twoEl->SetLineWidth(2);
    hRap_twoEl->GetXaxis()->SetRangeUser(3.0,6.0);
    hRap_twoEl->GetXaxis()->SetTitleSize(0.042);
    hRap_twoEl->GetXaxis()->SetLabelSize(0.042);
    hRap_twoEl->GetYaxis()->SetTitleSize(0.042);
    hRap_twoEl->GetYaxis()->SetLabelSize(0.042);
    hRap_twoEl->GetYaxis()->SetDecimals(1);
    hRap_twoEl->GetYaxis()->SetMaxDigits(3);
    hRap_twoEl->Draw("HIST");
    hRap_accFo->SetLineColor(kRed);
    hRap_accFo->SetLineWidth(2);
    hRap_accFo->Draw("HIST SAME");

    TLatex* latex = new TLatex();
    latex->SetTextSize(0.036);
    latex->SetTextAlign(21);
    latex->SetNDC();
    latex->DrawLatex(0.55,0.95,Form("STARlight: %s",FD_label[opt].Data()));

    // legends
    TLegend *l1 = new TLegend(0.50,0.75,0.90,0.90);
    l1->AddEntry(hRap_twoEl,"#it{y}(#psi') of generated events","L");
    l1->AddEntry(hRap_accFo,"#it{y}(#psi') of events with 3.4 < #it{#eta}^{e^{#pm}} < 5.8","L");
    l1->AddEntry((TObject*)0,Form("Acc^{[1]}: %.2f %%", frac * 1e2),"");
    l1->SetTextSize(0.036);
    l1->SetBorderSize(0);
    l1->SetFillStyle(0);
    l1->SetMargin(0.);
    l1->Draw();

    c->Print("results/starlightRapDep/h" + FD_folder[opt] + ".pdf");
    delete c;
    return;
}

void StarlightRapDep()
{
    gSystem->Exec("mkdir -p results/starlightRapDep/");

    nRapBins = (fRapUpp-fRapLow) / fRapStep;
    Printf("%i rapidity bins will be used.", nRapBins);

    // direct SL samples
    for(Int_t i = 0; i < 7; i++) {
        ConvertStarlightAsciiToTree(5e6,"inputData/starlight/" + SL_folder[i] + "/");
        CalculateRapDep(i);
    }
    // feed-down samples (from AliDPG)
    for(Int_t i = 0; i < 2; i++) {
        AnalyzeFeedDown(i);
    }
    return;
}