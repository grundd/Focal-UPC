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
#include "TCanvas.h"
#include "TLegend.h"
// my headers
#include "StarlightRapDep.h"

Float_t fRapLow = -6.;
Float_t fRapUpp = +6.;
Float_t fRapStep = 0.2; // [-]
Float_t fLumiRun4 = 6.9; // [nb^(-1)]
Float_t fTotalCS[6] = {38.8, // [mb]
                       17.8,
                       7.5,
                       3.1,
                       0.094,
                       0.041
}; 
TString sLabel[6] = {"coherent J/#psi",
                  "incoherent J/#psi",
                  "coherent #psi'",
                  "incoherent #psi'",
                  "coherent #Upsilon(1S)",
                  "incoherent #Upsilon(1S)"
};
TString sMC[6] = {"CohJ",
                  "IncJ",
                  "CohP",
                  "IncP",
                  "CohU",
                  "IncU"
};
Float_t BR[6] = {0.0594, // J/psi to e^(+)e^(-)
                 0.0594, 
                 0.00772, // psi' to e^(+)e^(-)
                 0.00772, 
                 0.0238, // Upsilon to e^(+)e^(-)
                 0.0238  
};

void CalculateRapDep(Int_t opt)
{
    gSystem->Exec("mkdir -p results/starlightRapDep/");
    Bool_t isUpsilon = kFALSE;
    if(opt >= 4) isUpsilon = kTRUE;

    Int_t nRapBins = (fRapUpp-fRapLow) / fRapStep;
    Printf("%i rapidity bins will be used.", nRapBins);

    TH1F* hRap_all = NULL;
    TH1F* hRap_twoEl = NULL;
    TH1F* hRap_accFo = NULL;
    TFile* fOut = TFile::Open("results/starlightRapDep/h" + sMC[opt] + ".root","read");
    if(fOut)
    {
        Printf("Histogram hRap_all already created.");
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
        Printf("Histogram hRap_all will be created.");

        TFile* fSL = TFile::Open("inputData/pureStarlightSim/" + sMC[opt] + "/tSTARlight.root", "read");
        if(fSL) Printf("File %s loaded.", fSL->GetName());

        TTree* tSL = dynamic_cast<TTree*> (fSL->Get("starlightTree"));
        if(tSL) Printf("Tree %s loaded.", tSL->GetName());
        SetBranchAddresses_tSL(tSL);

        Int_t nEntries = tSL->GetEntries();
        Printf("tSTARlight contains %i entries.", nEntries);

        gROOT->cd();
        TList* l = new TList();
        hRap_all = new TH1F("hRap_all","rap. distribution of generated events",nRapBins,fRapLow,fRapUpp);
        hRap_twoEl = new TH1F("hRap_twoEl","rap. distribution of events with two electrons",nRapBins,fRapLow,fRapUpp);
        hRap_accFo = new TH1F("hRap_accFo","rap. distribution of events with both electrons in FoCal acceptance",nRapBins,fRapLow,fRapUpp);
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
        fOut = new TFile("results/starlightRapDep/h" + sMC[opt] + ".root","RECREATE");
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
    gr->SetTitle(";|#it{y}| [-]; d#sigma/d#it{y} [mb]");
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
    TCanvas* c = new TCanvas("c","c",700,600);
    //c->SetLogy(); 
    c->SetTopMargin(0.05);
    c->SetBottomMargin(0.1);
    c->SetRightMargin(0.02);
    c->SetLeftMargin(0.12);
    gr->Draw("AC");

    // legend
    Int_t binNegLow = 2; Printf("Low edge of bin %i: %.1f", binNegLow, hRap_CS->GetBinLowEdge(binNegLow));
    Int_t binNegUpp = 13; Printf("Upp edge of bin %i: %.1f", binNegUpp, hRap_CS->GetBinLowEdge(binNegUpp+1));
    Int_t binPosLow = 48; Printf("Low edge of bin %i: %.1f", binPosLow, hRap_CS->GetBinLowEdge(binPosLow));
    Int_t binPosUpp = 59; Printf("Upp edge of bin %i: %.1f", binPosUpp, hRap_CS->GetBinLowEdge(binPosUpp+1));
    Float_t sigmaTotal = hRap_CS->Integral("width");
    Float_t sigmaForwd = hRap_CS->Integral(binNegLow,binNegUpp,"width") + hRap_CS->Integral(binPosLow,binPosUpp,"width");
    Float_t sigmaFocal = hRap_CS->Integral(binPosLow,binPosUpp,"width");

    TLegend *l1 = new TLegend(0.62,0.76,0.97,0.94);
    l1->AddEntry((TObject*)0,Form("#bf{STARlight: %s}",sLabel[opt].Data()),"");
    if(!isUpsilon) {
        l1->AddEntry((TObject*)0,Form("#sigma(|#it{y}| < 6.0) = %.2f mb",sigmaTotal),"");
        l1->AddEntry((TObject*)0,Form("#sigma(3.4 < |#it{y}| < 5.8) = %.2f mb",sigmaForwd),"");
        l1->AddEntry((TObject*)0,Form("#sigma(3.4 < #it{y} < 5.8) = %.2f mb",sigmaFocal),"");
    } else {
        l1->AddEntry((TObject*)0,Form("#sigma(|#it{y}| < 6.0) = %.2f #mub",sigmaTotal*1e3),"");
        l1->AddEntry((TObject*)0,Form("#sigma(3.4 < |#it{y}| < 5.8) = %.2f #mub",sigmaForwd*1e3),"");
        l1->AddEntry((TObject*)0,Form("#sigma(3.4 < #it{y} < 5.8) = %.2f #mub",sigmaFocal*1e3),"");
    }
    l1->SetTextSize(0.036);
    l1->SetBorderSize(0);
    l1->SetFillStyle(0);
    l1->SetMargin(0.);
    l1->Draw();
    // calculate expected yields of J/psi, psi' and Y in FoCal in Run 4
    // convert lumi from nb^(-1) to mb^(-1)
    Float_t nEvTot = hRap_all->Integral();
    Float_t nEv2El = hRap_twoEl->Integral();
    Float_t nEvAcc = hRap_accFo->Integral();
    Float_t acc = nEvAcc / nEv2El;
    Float_t lumi = fLumiRun4 * 1e6; // mb^(-1)
    Float_t yield = lumi * sigmaTotal * BR[opt] * acc;
    ofstream of;
    of.open("results/starlightRapDep/log" + sMC[opt] + ".txt");
    of << Form("generated events: %.0f\n", nEvTot);
    of << Form("gen ev with 2 ele: %.0f\n", nEv2El);
    of << Form("gen ev with 2 ele in Foc acc: %.0f\n", nEvAcc);
    of << Form("acceptance: %.3f\n", acc);
    of << Form("int lumi for UPC Pb+Pb in Run 4: %.2e mb^(-1)\n",lumi);
    of << Form("total sigma for %s: %.2e mb\n",sLabel[opt].Data(),sigmaTotal);
    of << Form("sigma for FoCal %s: %.2e mb\n",sLabel[opt].Data(),sigmaTotal * BR[opt] * acc * 1e3);
    of << Form("rough expected yield with VM rapidity within 3.4 and 5.8: %.2e\n",yield);
    of.close();

    TLegend *l2 = new TLegend(0.15,0.72,0.55,0.94);
    l2->AddEntry((TObject*)0,Form("gen. events: %.2e",nEv2El),"");
    l2->AddEntry((TObject*)0,Form("e^{+}e^{#minus} in 3.4 < |#it{#eta}| < 5.8: %.2e",nEvAcc),"");
    l2->AddEntry((TObject*)0,Form("Acc: %.2f %%",acc * 1e2),"");
    l2->AddEntry((TObject*)0,Form("Expected yield (#it{L}_{int} = 10 nb^{#minus1}):"),"");
    l2->AddEntry((TObject*)0,Form("#bf{%2.f}",yield),"");
    l2->SetTextSize(0.036);
    l2->SetBorderSize(0);
    l2->SetFillStyle(0);
    l2->SetMargin(0.);
    l2->Draw();
    // save the canvas
    c->Print("results/starlightRapDep/gr" + sMC[opt] + ".pdf");
    delete c;
    delete gr;
    return;
}

void StarlightRapDep()
{
    for(Int_t i = 0; i < 6; i ++) {
        ConvertStarlightAsciiToTree(5e6,"inputData/pureStarlightSim/" + sMC[i] + "/");
        CalculateRapDep(i);
    }
    return;
}