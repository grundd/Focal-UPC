// FocalUpcGrid_CreateHistograms.h
// David Grund, Nov 07, 2022

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TObjArray.h"

// binning
const Int_t nBinsEn = 100; // bin per 2 GeV
const Float_t lowEn = 0.;
const Float_t uppEn = 200.;
const Int_t nBinsPt = 80; // bin per 25 MeV
const Float_t lowPt = 0.;
const Float_t uppPt = 2.;
const Int_t nBinsYEta = 100; // bin per 0.05
const Float_t lowYEta = 2.;
const Float_t uppYEta = 7.;
const Int_t nBinsPhi = 80; // bin per 0.04
const Float_t lowPhi = 0.;
const Float_t uppPhi = 3.2;
const Int_t nBinsM = 50; // bin per 100 MeV
const Float_t lowM = 0.;
const Float_t uppM = 5.;
const Int_t nBinsSep = 50;
const Float_t lowSep = 0.;
const Float_t uppSep = 100.;

enum kGridBox {
    //* TH1F histograms:
    kGridBox_firstTH1F = 0,
    kB1_mcE,
    kB1_mcPt,
    kB1_all,
    //* TH2F histograms:
    kGridBox_firstTH2F,
    kB2_clMcDX_clMcDY,
    kB2_mcE_clMcSep,
    kB2_clMaxClDX_clMaxClDY,
    kB2_mcE_clMaxClSep,
    kB2_clX_clMaxClSep,
    kB2_clY_clMaxClSep,
    kB2_mcE_nCls,
    kB2_clE_mcE,
    kB2_totE_mcE,
    kB2_maxClE_mcE,
    //* TProfile histograms:
    kGridBox_firstTPrf,
    kBP_totE_mcE,
    kBP_maxClE_mcE,
    kGridBox_all
};

enum kGridJpsi {
    //* TH1F histograms:
    kGridJpsi_firstTH1F = 0,
    kJ1_clZ,
    // MC J/psi
    kJ1_mcJEn,
    kJ1_mcJPt,
    kJ1_mcJRap,
    kJ1_mcJM,
    // ppe pairs
    kJ1_mcJElPairEn,
    kJ1_mcJElPairPt,
    kJ1_mcJElPairRap,
    kJ1_mcJElPairM,
    //* TH2F histograms:
    kGridJpsi_firstTH2F,
    // MC J/psi kinematics
    kJ2_mcJRap_mcJPt,
    // ppe kinematics
    kJ2_mcJElEta_mcJElPt,
    kJ2_mcJEl1En_mcJEl2En,
    kJ2_mcJEl1Pt_mcJEl2Pt,
    // cluster kinematics corr
    kJ2_clEta_clPhi,
    kJ2_clEta_clPt,
    // cluster vs ppp
    kJ2_pppClEn_mtchEn,
    kJ2_pppClEta_mtchEta,
    kJ2_pppClPt_mtchPt,
    // cluster vs ppe
    kJ2_ppeClEn_mtchEn,
    kJ2_ppeClEta_mtchEta,
    kJ2_ppeClPt_mtchPt,
    //* TProfile histograms:
    kGridJpsi_firstTPrf,
    kGridJpsi_all
};

enum kMainJpsi {
    //* TH1F histograms
    kMainJpsi_firstTH1F = 0,
    // cluster pairs
    kJ1_clPairEn,
    kJ1_clPairPt,
    kJ1_clPairPt_massCut,
    kJ1_clPairRap,
    kJ1_clPairM,
    kJ1_clPairSep,
    // matched cl. pairs
    kJ1_ppeClPairM,
    kJ1_ppeClPairSep,
    kJ1_sameppeClPairSep,
    //* TH2F histograms:
    kMainJpsi_firstTH2F,
    // cluster pairs vs ppe pairs
    kJ2_ppeClPairEn_mtchEn,
    kJ2_ppeClPairRap_mtchRap,
    kJ2_ppeClPairPt_mtchPt,
    kJ2_ppeClPairM_mtchM,
    kJ2_clPairSep_mcJElSep,
    //* TProfile histograms:
    kMainJpsi_firstTPrf,
    kMainJpsi_all
};

// ******************************************************************************************************************
// GRID analysis -> box electron/photon events:
// ******************************************************************************************************************

void CreateHistos_GridBox(TObjArray* objArr)
{
    objArr->SetOwner();
    //* TH1F histograms:
    TH1F* hB1_mcE = new TH1F("hB1_mcE","",nBinsEn,lowEn,uppEn);
                    hB1_mcE->SetTitle("#it{E} of generated electron;#it{E}_{MC} [GeV];counts");
                    objArr->AddAt(hB1_mcE, kB1_mcE);
    TH1F* hB1_mcPt = new TH1F("hB1_mcPt","",40,lowPt,uppPt);
                    hB1_mcPt->SetTitle("#it{p}_{T} of generated electron;#it{p}_{T} [GeV/#it{c}];counts");
                    objArr->AddAt(hB1_mcPt, kB1_mcPt);
    //* TH2F histograms:
    // cluster XY position vs that of ppe / cluster with maximum energy
    TH2F* hB2_clMcDX_clMcDY = new TH2F("hB2_clMcDX_clMcDY","",100,-20.,20.,100,-20.,20.);
                    hB2_clMcDX_clMcDY->SetTitle("Distance (in the XY plane) between the cluster and the MC track;#Delta#it{x} = #it{x}_{cl} - #it{x}_{MC} [cm];#Delta#it{y} = #it{y}_{cl} - #it{y}_{MC} [cm]");
                    objArr->AddAt(hB2_clMcDX_clMcDY, kB2_clMcDX_clMcDY);
    TH2F* hB2_mcE_clMcSep = new TH2F("hB2_mcE_clMcSep","",nBinsEn,lowEn,uppEn,100,0.,20.);
                    hB2_mcE_clMcSep->SetTitle("Radial distance #Delta#it{R} = #sqrt{(#Delta#it{x})^{2} + (#Delta#it{y})^{2}} between the cluster and the MC track;#it{E}_{MC} [GeV];#Delta#it{R} [cm]");
                    objArr->AddAt(hB2_mcE_clMcSep, kB2_mcE_clMcSep);
    TH2F* hB2_clMaxClDX_clMaxClDY = new TH2F("hB2_clMaxClDX_clMaxClDY","",100,-20.,20.,100,-20.,20.);
                    hB2_clMaxClDX_clMaxClDY->SetTitle("Distance between the cluster and the cluster with maximum energy;#Delta#it{x} = #it{x}_{cl} - #it{x}_{cl,max} [cm];#Delta#it{y} = #it{y}_{cl} - #it{y}_{cl,max} [cm]");
                    objArr->AddAt(hB2_clMaxClDX_clMaxClDY, kB2_clMaxClDX_clMaxClDY);
    TH2F* hB2_mcE_clMaxClSep = new TH2F("hB2_mcE_clMaxClSep","",nBinsEn,lowEn,uppEn,100,0.,20.);
                    hB2_mcE_clMaxClSep->SetTitle("Radial distance #Delta#it{R} = #sqrt{(#Delta#it{x})^{2} + (#Delta#it{y})^{2}} between the cluster and the cluster with maximum energy;#it{E}_{MC} [GeV];#Delta#it{R} [cm]");
                    objArr->AddAt(hB2_mcE_clMaxClSep, kB2_mcE_clMaxClSep);
    TH2F* hB2_clX_clMaxClSep = new TH2F("hB2_clX_clMaxClSep","",100,-50.,50.,100,0.,20.);
                    hB2_clX_clMaxClSep->SetTitle("Radial dist. #Delta#it{R} = #sqrt{(#Delta#it{x})^{2} + (#Delta#it{y})^{2}} between the cluster and the cluster with maximum energy vs average X coordinate of the two;#it{X}_{cl} [GeV];#Delta#it{R} [cm]");
                    objArr->AddAt(hB2_clX_clMaxClSep, kB2_clX_clMaxClSep);
    TH2F* hB2_clY_clMaxClSep = new TH2F("hB2_clY_clMaxClSep","",100,-50.,50.,100,0.,20.);
                    hB2_clY_clMaxClSep->SetTitle("Radial dist. #Delta#it{R} = #sqrt{(#Delta#it{x})^{2} + (#Delta#it{y})^{2}} between the cluster and the cluster with maximum energy vs average Y coordinate of the two;#it{Y}_{cl} [GeV];#Delta#it{R} [cm]");
                    objArr->AddAt(hB2_clY_clMaxClSep, kB2_clY_clMaxClSep);
    // number of clusters, cluster energy, total energy and maximum cluster energy vs MC energy
    TH2F* hB2_mcE_nCls = new TH2F("hB2_mcE_nCls","",nBinsEn,lowEn,uppEn,8,-0.5,7.5);
                    hB2_mcE_nCls->SetTitle("(#it{E}_{MC}, #it{N}_{cls});#it{E}_{MC} [GeV];#it{N}_{cls} [-]");
                    objArr->AddAt(hB2_mcE_nCls, kB2_mcE_nCls);
    TH2F* hB2_clE_mcE = new TH2F("hB2_clE_mcE","",nBinsEn,lowEn,uppEn,nBinsEn,lowEn,uppEn);
                    hB2_clE_mcE->SetTitle("(#it{E}_{cl}, #it{E}_{MC});#it{E}_{cl} [GeV];#it{E}_{MC} [GeV]");
                    objArr->AddAt(hB2_clE_mcE, kB2_clE_mcE);
    TH2F* hB2_totE_mcE = new TH2F("hB2_totE_mcE","",nBinsEn,lowEn,uppEn,nBinsEn,lowEn,uppEn);
                    hB2_totE_mcE->SetTitle("(#it{E}_{total}, #it{E}_{MC});#it{E}_{total} [GeV];#it{E}_{MC} [GeV]");
                    objArr->AddAt(hB2_totE_mcE, kB2_totE_mcE);
    TH2F* hB2_maxClE_mcE = new TH2F("hB2_maxClE_mcE","",nBinsEn,lowEn,uppEn,nBinsEn,lowEn,uppEn);
                    hB2_maxClE_mcE->SetTitle("(#it{E}_{cl,max}, #it{E}_{MC});#it{E}_{cl,max} [GeV];#it{E}_{MC} [GeV]");
                    objArr->AddAt(hB2_maxClE_mcE, kB2_maxClE_mcE);
    //* TProfile histograms:
    TProfile* hBP_totE_mcE = new TProfile("hBP_totE_mcE","",nBinsEn,lowEn,uppEn,lowEn,uppEn);
                    hBP_totE_mcE->SetTitle("(#it{E}_{total}, #it{E}_{MC});#it{E}_{total} [GeV];#it{E}_{MC} [GeV]");
                    objArr->AddAt(hBP_totE_mcE, kBP_totE_mcE);
    TProfile* hBP_maxClE_mcE = new TProfile("hBP_maxClE_mcE","",nBinsEn,lowEn,uppEn,lowEn,uppEn);
                    hBP_maxClE_mcE->SetTitle("(#it{E}_{cl,max}, #it{E}_{MC});#it{E}_{cl,max} [GeV];#it{E}_{MC} [GeV]");
                    objArr->AddAt(hBP_maxClE_mcE, kBP_maxClE_mcE);
    return;
}

// ******************************************************************************************************************
// GRID analysis -> STARlight J/psi simulations:
// ******************************************************************************************************************

void CreateHistos_GridJpsi(TObjArray* objArr)
{
    objArr->SetOwner();
    //* TH1F histograms:
    TH1F* hJ1_clZ = new TH1F("hJ1_clZ","",60,700.,715.);
                    hJ1_clZ->SetTitle("#it{z} coordinate of FoCal clusters;#it{z}_{cl} [cm];counts");
                    objArr->AddAt(hJ1_clZ, kJ1_clZ);
    // MC kinematics
    // J/psi
    TH1F* hJ1_mcJEn = new TH1F("hJ1_mcJEn","",nBinsEn,lowEn,uppEn);
                    hJ1_mcJEn->SetTitle("#it{E} of generated J/#psi;#it{E}_{J/#psi} [GeV];counts");
                    objArr->AddAt(hJ1_mcJEn, kJ1_mcJEn);
    TH1F* hJ1_mcJPt = new TH1F("hJ1_mcJPt","",nBinsPt,lowPt,uppPt);
                    hJ1_mcJPt->SetTitle("#it{p}_{T} of generated J/#psi;#it{p}_{T,J/#psi} [GeV/#it{c}];counts");
                    objArr->AddAt(hJ1_mcJPt, kJ1_mcJPt);
    TH1F* hJ1_mcJRap = new TH1F("hJ1_mcJRap","",nBinsYEta,lowYEta,uppYEta);
                    hJ1_mcJRap->SetTitle("#it{y} of generated J/#psi;#it{y}_{J/#psi} [-];counts");
                    objArr->AddAt(hJ1_mcJRap, kJ1_mcJRap);
    TH1F* hJ1_mcJM = new TH1F("hJ1_mcJM","",nBinsM,lowM,uppM);
                    hJ1_mcJM->SetTitle("#it{m} of generated J/#psi;#it{m}_{J/#psi} [GeV/#it{c}^{2}];counts");
                    objArr->AddAt(hJ1_mcJM, kJ1_mcJM);
    // pairs of pp electrons
    TH1F* hJ1_mcJElPairEn = new TH1F("hJ1_mcJElPairEn","",nBinsEn,lowEn,uppEn);
                    hJ1_mcJElPairEn->SetTitle("#it{E} of pairs of pp electrons;#it{E}_{ppe pair} [GeV];counts");
                    objArr->AddAt(hJ1_mcJElPairEn, kJ1_mcJElPairEn);
    TH1F* hJ1_mcJElPairPt = new TH1F("hJ1_mcJElPairPt","",nBinsPt,lowPt,uppPt);
                    hJ1_mcJElPairPt->SetTitle("#it{p}_{T} of pairs of pp electrons;#it{p}_{T,ppe pair} [GeV/#it{c}];counts");
                    objArr->AddAt(hJ1_mcJElPairPt, kJ1_mcJElPairPt);
    TH1F* hJ1_mcJElPairRap = new TH1F("hJ1_mcJElPairRap","",nBinsYEta,lowYEta,uppYEta);
                    hJ1_mcJElPairRap->SetTitle("#it{y} of pairs of pp electrons;#it{y}_{ppe pair} [-];counts");
                    objArr->AddAt(hJ1_mcJElPairRap, kJ1_mcJElPairRap);
    TH1F* hJ1_mcJElPairM = new TH1F("hJ1_mcJElPairM","",nBinsM,lowM,uppM);
                    hJ1_mcJElPairM->SetTitle("#it{m} of pairs of pp electrons;#it{m}_{ppe pair} [GeV/#it{c}^{2}];counts");
                    objArr->AddAt(hJ1_mcJElPairM, kJ1_mcJElPairM);
    //* TH2F histograms:
    // correlation of MC kinematics
    // J/psi
    TH2F* hJ2_mcJRap_mcJPt = new TH2F("hJ2_mcJRap_mcJPt","",nBinsYEta,lowYEta,uppYEta,nBinsPt,lowPt,uppPt);
                    hJ2_mcJRap_mcJPt->SetTitle("(#it{y}, #it{p}_{T}) of generated J/#psi;#it{y}_{J/#psi} [-];#it{p}_{T,J/#psi} [GeV/#it{c}]");
                    objArr->AddAt(hJ2_mcJRap_mcJPt, kJ2_mcJRap_mcJPt);
    // pp electrons
    TH2F* hJ2_mcJElEta_mcJElPt = new TH2F("hJ2_mcJElEta_mcJElPt","",nBinsYEta,lowYEta,uppYEta,nBinsPt,lowPt,uppPt);
                    hJ2_mcJElEta_mcJElPt->SetTitle("(#eta, #it{p}_{T}) of a pp electron;#eta_{ppe^{#pm}} [-];#it{p}_{T,ppe^{#pm}} [GeV/#it{c}]");
                    objArr->AddAt(hJ2_mcJElEta_mcJElPt, kJ2_mcJElEta_mcJElPt);
    TH2F* hJ2_mcJEl1En_mcJEl2En = new TH2F("hJ2_mcJEl1En_mcJEl2En","",nBinsEn,lowEn,uppEn,nBinsEn,lowEn,uppEn);
                    hJ2_mcJEl1En_mcJEl2En->SetTitle("correlation between energy of two pp electrons;#it{E}_{ppe^{#pm},1} [GeV];#it{E}_{ppe^{#pm},2} [GeV]");
                    objArr->AddAt(hJ2_mcJEl1En_mcJEl2En, kJ2_mcJEl1En_mcJEl2En);
    TH2F* hJ2_mcJEl1Pt_mcJEl2Pt = new TH2F("hJ2_mcJEl1Pt_mcJEl2Pt","",nBinsPt,lowPt,uppPt,nBinsPt,lowPt,uppPt);
                    hJ2_mcJEl1Pt_mcJEl2Pt->SetTitle("correlation between transverse momenta of two pp electrons;#it{p}_{T,ppe^{#pm},1} [GeV/#it{c}];#it{p}_{T,ppe^{#pm},2} [GeV/#it{c}]");
                    objArr->AddAt(hJ2_mcJEl1Pt_mcJEl2Pt, kJ2_mcJEl1Pt_mcJEl2Pt);
    // in the following, "cluster" denotes a prefiltered cluster
    // correlation of cluster kinematics
    TH2F* hJ2_clEta_clPhi = new TH2F("hJ2_clEta_clPhi","",nBinsYEta,lowYEta,uppYEta,nBinsPhi,lowPhi,uppPhi);
                    hJ2_clEta_clPhi->SetTitle("(#eta, #phi) of cluster;#eta_{cl} [-];#phi_{cl} [-]");
                    objArr->AddAt(hJ2_clEta_clPhi, kJ2_clEta_clPhi);
    TH2F* hJ2_clEta_clPt = new TH2F("hJ2_clEta_clPt","",nBinsYEta,lowYEta,uppYEta,nBinsPt,lowPt,uppPt);
                    hJ2_clEta_clPt->SetTitle("(#eta, #it{p}_{T}) of cluster;#eta_{cl} [-];#it{p}_{T,cl} [GeV/#it{c}]");
                    objArr->AddAt(hJ2_clEta_clPt, kJ2_clEta_clPt);
    // cluster kinematics vs kinematics of matched ppp
    TH2F* hJ2_pppClEn_mtchEn = new TH2F("hJ2_pppClEn_mtchEn","",nBinsEn,lowEn,uppEn,nBinsEn,lowEn,uppEn);
                    hJ2_pppClEn_mtchEn->SetTitle("energy of a cluster matched with a pp particle vs energy ppp;#it{E}_{cl} [GeV];#it{E}_{matched ppp} [GeV]");
                    objArr->AddAt(hJ2_pppClEn_mtchEn, kJ2_pppClEn_mtchEn);
    TH2F* hJ2_pppClEta_mtchEta = new TH2F("hJ2_pppClEta_mtchEta","",nBinsYEta,lowYEta,uppYEta,nBinsYEta,lowYEta,uppYEta);
                    hJ2_pppClEta_mtchEta->SetTitle("#eta of a cluster matched with a pp particle vs #eta of ppp;#eta_{cl} [-];#eta_{matched ppp} [-]");
                    objArr->AddAt(hJ2_pppClEta_mtchEta, kJ2_pppClEta_mtchEta);
    TH2F* hJ2_pppClPt_mtchPt = new TH2F("hJ2_pppClPt_mtchPt","",nBinsPt,lowPt,uppPt,nBinsPt,lowPt,uppPt);
                    hJ2_pppClPt_mtchPt->SetTitle("#it{p}_{T} of a cluster matched with a pp particle vs #it{p}_{T} of ppp;#it{p}_{T,cl} [GeV/#it{c}];#it{p}_{T,matched ppp} [GeV/#it{c}]");
                    objArr->AddAt(hJ2_pppClPt_mtchPt, kJ2_pppClPt_mtchPt);
    // cluster kinematics vs kinematics of matched ppe
    TH2F* hJ2_ppeClEn_mtchEn = new TH2F("hJ2_ppeClEn_mtchEn","",nBinsEn,lowEn,uppEn,nBinsEn,lowEn,uppEn); 
                    hJ2_ppeClEn_mtchEn->SetTitle("energy of a cluster matched with a pp electron vs vs energy of ppe;#it{E}_{cl} [GeV];#it{E}_{matched ppe^{#pm}} [GeV]");
                    objArr->AddAt(hJ2_ppeClEn_mtchEn, kJ2_ppeClEn_mtchEn);
    TH2F* hJ2_ppeClEta_mtchEta = new TH2F("hJ2_ppeClEta_mtchEta","",nBinsYEta,lowYEta,uppYEta,nBinsYEta,lowYEta,uppYEta);
                    hJ2_ppeClEta_mtchEta->SetTitle("#eta of a cluster matched with a pp electron vs #eta of ppe;#eta_{cl} [-];#eta_{matched ppe^{#pm}} [-]");
                    objArr->AddAt(hJ2_ppeClEta_mtchEta, kJ2_ppeClEta_mtchEta);
    TH2F* hJ2_ppeClPt_mtchPt = new TH2F("hJ2_ppeClPt_mtchPt","",nBinsPt,lowPt,uppPt,nBinsPt,lowPt,uppPt);
                    hJ2_ppeClPt_mtchPt->SetTitle("#it{p}_{T} of a cluster matched with a pp electron vs #it{p}_{T} of ppe;#it{p}_{T,cl} [GeV/#it{c}];#it{p}_{T,matched ppe^{#pm}} [GeV/#it{c}]");
                    objArr->AddAt(hJ2_ppeClPt_mtchPt, kJ2_ppeClPt_mtchPt);
    //* TProfile histograms:
    return;
}

// ******************************************************************************************************************
// Main analysis -> STARlight J/psi simulations:
// ******************************************************************************************************************

void CreateHistos_MainJpsi(TObjArray* objArr)
{
    objArr->SetOwner();
    //* TH1F histograms:
    // cluster pair kinematics
    TH1F* hJ1_clPairEn = new TH1F("hJ1_clPairEn","",nBinsEn,lowEn,uppEn);
                    hJ1_clPairEn->SetTitle("#it{E} of cluster pairs;#it{E}_{cl pair} [GeV];counts");
                    objArr->AddAt(hJ1_clPairEn, kJ1_clPairEn);
    TH1F* hJ1_clPairPt = new TH1F("hJ1_clPairPt","",nBinsPt,lowPt,uppPt);
                    hJ1_clPairPt->SetTitle("#it{p}_{T} of cluster pairs;#it{p}_{T,cl pair} [GeV/#it{c}];counts");
                    objArr->AddAt(hJ1_clPairPt, kJ1_clPairPt);
    TH1F* hJ1_clPairPt_massCut = new TH1F("hJ1_clPairPt_massCut","",nBinsPt,lowPt,uppPt);
                    hJ1_clPairPt_massCut->SetTitle("#it{p}_{T} of cluster pairs with #it{m} > 2.8 GeV/#it{c}^{2};#it{p}_{T,cl pair} [GeV/#it{c}];counts");
                    objArr->AddAt(hJ1_clPairPt_massCut, kJ1_clPairPt_massCut);
    TH1F* hJ1_clPairM = new TH1F("hJ1_clPairM","",nBinsM,lowM,uppM);
                    hJ1_clPairM->SetTitle("#it{m} of cluster pairs;#it{m}_{cl pair} [GeV/#it{c}^{2}];counts");
                    objArr->AddAt(hJ1_clPairM, kJ1_clPairM);
    TH1F* hJ1_clPairRap = new TH1F("hJ1_clPairRap","",nBinsYEta,lowYEta,uppYEta);
                    hJ1_clPairRap->SetTitle("rapidity of cluster pairs;#it{y}_{cl pair} [-];counts");
                    objArr->AddAt(hJ1_clPairRap, kJ1_clPairRap);
    TH1F* hJ1_clPairSep = new TH1F("hJ1_clPairSep","",nBinsSep,lowSep,uppSep);
                    hJ1_clPairSep->SetTitle("Radial separation of cluster pairs;#Delta#it{R}_{cl pair} [cm];counts");
                    objArr->AddAt(hJ1_clPairSep, kJ1_clPairSep);
    // cluster pairs matched with ppe pairs
    TH1F* hJ1_ppeClPairM = new TH1F("hJ1_ppeClPairM","",nBinsM,lowM,uppM);
                    hJ1_ppeClPairM->SetTitle("#it{m} of cluster pairs matched with a pair of pp electrons;#it{m}_{matched cl pair} [GeV/#it{c}^{2}];counts");
                    objArr->AddAt(hJ1_ppeClPairM, kJ1_ppeClPairM);
    TH1F* hJ1_ppeClPairSep = new TH1F("hJ1_ppeClPairSep","",nBinsSep,lowSep,uppSep);
                    hJ1_ppeClPairSep->SetTitle("Radial separation of cluster pairs matched with a pair of pp electrons;#Delta#it{R}_{matched cl pair} [cm];counts");
                    objArr->AddAt(hJ1_ppeClPairSep, kJ1_ppeClPairSep);
    TH1F* hJ1_sameppeClPairSep = new TH1F("hJ1_sameppeClPairSep","",nBinsSep,lowSep,uppSep);
                    hJ1_sameppeClPairSep->SetTitle("Radial separation of cluster pairs matched with the same pp electron;#Delta#it{R}_{matched cl pair} [cm];counts");
                    objArr->AddAt(hJ1_sameppeClPairSep, kJ1_sameppeClPairSep);
    //* TH2F histograms:
    // in the following, "cluster" denotes a prefiltered cluster
    // cluster pair kinematics vs kinematics of matched ppe pair
    TH2F* hJ2_ppeClPairEn_mtchEn = new TH2F("hJ2_ppeClPairEn_mtchEn","",nBinsEn,lowEn,uppEn,nBinsEn,lowEn,uppEn);
                    hJ2_ppeClPairEn_mtchEn->SetTitle("energy of a cluster pair matched with a pair of pp electrons vs energy of ppe pair;#it{E}_{cl pair} [GeV];#it{E}_{matched ppe pair} [GeV]");
                    objArr->AddAt(hJ2_ppeClPairEn_mtchEn, kJ2_ppeClPairEn_mtchEn);
    TH2F* hJ2_ppeClPairRap_mtchRap = new TH2F("hJ2_ppeClPairRap_mtchRap","",nBinsYEta,lowYEta,uppYEta,nBinsYEta,lowYEta,uppYEta);
                    hJ2_ppeClPairRap_mtchRap->SetTitle("#it{y} of a cluster pair matched with a pair of pp electrons vs #it{y} of ppe pair;#it{y}_{cl pair} [-];#it{y}_{matched ppe pair} [-]");
                    objArr->AddAt(hJ2_ppeClPairRap_mtchRap, kJ2_ppeClPairRap_mtchRap);
    TH2F* hJ2_ppeClPairPt_mtchPt = new TH2F("hJ2_ppeClPairPt_mtchPt","",nBinsPt,lowPt,uppPt,nBinsPt,lowPt,uppPt);
                    hJ2_ppeClPairPt_mtchPt->SetTitle("#it{p}_{T} of a cluster pair matched with a pair of pp electrons vs #it{p}_{T} of ppe pair;#it{p}_{T,cl pair} [GeV/#it{c}];#it{p}_{T,matched ppe pair} [GeV/#it{c}]");
                    objArr->AddAt(hJ2_ppeClPairPt_mtchPt, kJ2_ppeClPairPt_mtchPt);
    TH2F* hJ2_ppeClPairM_mtchM = new TH2F("hJ2_ppeClPairM_mtchM","",nBinsM,lowM,uppM,nBinsM,lowM,uppM);
                    hJ2_ppeClPairM_mtchM->SetTitle("#it{m} of a cluster pair matched with a pair of pp electrons vs #it{m} of ppe pair;#it{m}_{cl pair} [GeV/#it{c}];#it{m}_{matched ppe pair} [GeV/#it{c}]");
                    objArr->AddAt(hJ2_ppeClPairM_mtchM, kJ2_ppeClPairM_mtchM);
    TH2F* hJ2_clPairSep_mcJElSep = new TH2F("hJ2_clPairSep_mcJElSep","",nBinsSep,lowSep,uppSep,nBinsSep,lowSep,uppSep);
                    hJ2_clPairSep_mcJElSep->SetTitle("Radial separation of cluster pairs vs radial separation of ppe pair;#Delta#it{R}_{cl pair} [cm];#Delta#it{R}_{ppe pair} [cm]");
                    objArr->AddAt(hJ2_clPairSep_mcJElSep, kJ2_clPairSep_mcJElSep);
    return;
}