// FocalUpcAnalysisJp_DefineOutput.h
// David Grund, Oct 16, 2022

// binning
const Int_t nBinsPt = 80; // bin per 25 MeV
const Float_t lowPt = 0.;
const Float_t uppPt = 2.;
const Int_t nBinsEta = 100; // bin per 0.05
const Float_t lowEta = 2.;
const Float_t uppEta = 7.;
const Int_t nBinsE = 100; // bin per 2 GeV
const Float_t lowE = 0.;
const Float_t uppE = 200.;
const Int_t nBinsM = 50; // bin per 100 MeV
const Float_t lowM = 0.;
const Float_t uppM = 5.;
const Int_t nBinsPhi = 80; // bin per 0.04
const Float_t lowPhi = 0.;
const Float_t uppPhi = 3.2;

enum kBx_TH1F {
    kBx_mcE = 0,
    kBx_mcPt,
    kBx_TH1F_all    
};

enum kBx_TH2F {
    kBx_mcE_nCls = 0,
    kBx_mcE_totE,
    kBx_totE_mcE,
    kBx_totEwHCalE_mcE,
    kBx_totEwIsoR2E_mcE,
    kBx_totEwIsoR4E_mcE,
    kBx_totE_totHCalE,
    kBx_totE_totIsoR2E,
    kBx_totE_totIsoR4E,
    kBx_totEFromSegCls_mcE,
    kBx_mcE_maxClE,
    kBx_maxClE_mcE,
    kBx_clE_mcE,
    kBx_TH2F_all
};

enum kBx_TPrf {
    kBx_mcE_totE_prof = 0,
    kBx_totE_mcE_prof,
    kBx_mcE_maxClE_prof,
    kBx_maxClE_mcE_prof,
    kBx_TPrf_all
};

enum kJp_TH1F {
    kJp_mcJPt = 0,
    kJp_clPairM,
    kJp_clPairEta,
    kJp_clPairPt,
    kJp_clPairPt_massCut,
    kJp_primElClPairM,
    kJp_primElClPairM_dir,
    //
    kJp_primElClPairEta,
    kJp_primElClPairPt,
    kJp_mcJEEall,
    kJp_mcJEEmtch,
    kJp_mcJEEratio,
    kJp_TH1F_all
};

enum kJp_TH2F {
    kJp_mcJRap_mcJPt = 0,
    kJp_mcJRap_mcJPt_acc,
    kJp_mcJEEta_mcJEPt,
    kJp_mcJE1E_mcJE2E,
    kJp_mcJE1Pt_mcJE2Pt,
    kJp_clEta_clPhi,
    kJp_clEta_clPt,
    // cluster vs ppp
    kJp_clE_mtchE,
    kJp_clE_mtchDirE,
    kJp_clEta_mtchEta,
    kJp_clEta_mtchDirEta,
    kJp_clPt_mtchPt,
    kJp_clPt_mtchDirPt,
    // cluster vs ppe
    kJp_primElClE_mtchE,
    kJp_primElClE_mtchDirE,
    kJp_primElClEta_mtchEta,
    kJp_primElClEta_mtchDirEta,
    kJp_primElClPt_mtchPt,
    kJp_primElClPt_mtchDirPt,
    //
    kJp_primElClPairEta_mcJEta,
    kJp_primElClPairPt_mcJPt,
    // matching pp (direct vs by mother)
    kJp_mtchE_mtchDirE,
    kJp_mtchE_mtchDirE_primEl,
    kJp_mtchEta_mtchDirEta_primEl,
    kJp_mtchPt_mtchDirPt_primEl,
    kJp_TH2F_all
};

enum kJp_TPrf {
    kJp_TPrf_all = 0
};

// ******************************************************************************************************************
// for analyses of box electron/photon: events:
// ******************************************************************************************************************

void DefineHisto_Bx_TH1F(TObjArray* objArr)
{
    objArr->SetOwner();

    TH1F* hBx_mcE = new TH1F("hBx_mcE","",nBinsE,lowE,uppE);
                    hBx_mcE->SetTitle("electron #it{E} generated;#it{E}_{MC} [GeV];counts");
                    objArr->AddAt(hBx_mcE, kBx_mcE);
    TH1F* hBx_mcPt = new TH1F("hBx_mcPt","",40,lowPt,uppPt);
                    hBx_mcPt->SetTitle("electron #it{p}_{T} generated;#it{p}_{T} [GeV/#it{c}];counts");
                    objArr->AddAt(hBx_mcPt, kBx_mcPt);
    return;
}

void DefineHisto_Bx_TH2F(TObjArray* objArr)
{
    objArr->SetOwner();

    TH2F* hBx_mcE_nCls = new TH2F("hBx_mcE_nCls","",nBinsE,lowE,uppE,8,-0.5,7.5);
                    hBx_mcE_nCls->SetTitle("(#it{E}_{MC}, #it{N}_{cls});#it{E}_{MC} [GeV];#it{N}_{cls} [-]");
                    objArr->AddAt(hBx_mcE_nCls, kBx_mcE_nCls);
    TH2F* hBx_mcE_totE = new TH2F("hBx_mcE_totE","",nBinsE,lowE,uppE,nBinsE,lowE,uppE);
                    hBx_mcE_totE->SetTitle("(#it{E}_{MC}, #it{E}_{total});#it{E}_{MC} [GeV];#it{E}_{total} [GeV]");
                    objArr->AddAt(hBx_mcE_totE, kBx_mcE_totE);
    TH2F* hBx_totE_mcE = new TH2F("hBx_totE_mcE","",nBinsE,lowE,uppE,nBinsE,lowE,uppE);
                    hBx_totE_mcE->SetTitle("(#it{E}_{total}, #it{E}_{MC});#it{E}_{total} [GeV];#it{E}_{MC} [GeV]");
                    objArr->AddAt(hBx_totE_mcE, kBx_totE_mcE);
    TH2F* hBx_totEwHCalE_mcE = new TH2F("hBx_totEwHCalE_mcE","",nBinsE,lowE,uppE,nBinsE,lowE,uppE);
                    hBx_totEwHCalE_mcE->SetTitle("(#it{E}_{total with HCal}, #it{E}_{MC});#it{E}_{total with HCal} [GeV];#it{E}_{MC} [GeV]");
                    objArr->AddAt(hBx_totEwHCalE_mcE, kBx_totEwHCalE_mcE);
    TH2F* hBx_totEwIsoR2E_mcE = new TH2F("hBx_totEwIsoR2E_mcE","",nBinsE,lowE,uppE,nBinsE,lowE,uppE);
                    hBx_totEwIsoR2E_mcE->SetTitle("(#it{E}_{total with HCal isoR2}, #it{E}_{MC});#it{E}_{total with HCal isoR2} [GeV];#it{E}_{MC} [GeV]");
                    objArr->AddAt(hBx_totEwIsoR2E_mcE, kBx_totEwIsoR2E_mcE);
    TH2F* hBx_totEwIsoR4E_mcE = new TH2F("hBx_totEwIsoR4E_mcE","",nBinsE,lowE,uppE,nBinsE,lowE,uppE);
                    hBx_totEwIsoR4E_mcE->SetTitle("(#it{E}_{total with HCal isoR4}, #it{E}_{MC});#it{E}_{total with HCal isoR4} [GeV];#it{E}_{MC} [GeV]");
                    objArr->AddAt(hBx_totEwIsoR4E_mcE, kBx_totEwIsoR4E_mcE);
    TH2F* hBx_totE_totHCalE = new TH2F("hBx_totE_totHCalE","",nBinsE,lowE,uppE,nBinsE,lowE,uppE);
                    hBx_totE_totHCalE->SetTitle("(#it{E}_{total}, #it{E}_{total HCal});#it{E}_{total} [GeV];#it{E}_{total HCal} [GeV]");
                    objArr->AddAt(hBx_totE_totHCalE, kBx_totE_totHCalE);
    TH2F* hBx_totE_totIsoR2E = new TH2F("hBx_totE_totIsoR2E","",nBinsE,lowE,uppE,nBinsE,lowE,uppE);
                    hBx_totE_totIsoR2E->SetTitle("(#it{E}_{total}, #it{E}_{total isoR2});#it{E}_{total} [GeV];#it{E}_{total isoR2} [GeV]");
                    objArr->AddAt(hBx_totE_totIsoR2E, kBx_totE_totIsoR2E);
    TH2F* hBx_totE_totIsoR4E = new TH2F("hBx_totE_totIsoR4E","",nBinsE,lowE,uppE,nBinsE,lowE,uppE);
                    hBx_totE_totIsoR4E->SetTitle("(#it{E}_{total}, #it{E}_{total isoR4});#it{E}_{total} [GeV];#it{E}_{total isoR4} [GeV]");
                    objArr->AddAt(hBx_totE_totIsoR4E, kBx_totE_totIsoR4E);
    TH2F* hBx_totEFromSegCls_mcE = new TH2F("hBx_totEFromSegCls_mcE","",nBinsE,lowE,2500,nBinsE,lowE,uppE);
                    hBx_totEFromSegCls_mcE->SetTitle("(#it{E}_{total from per seg cls}, #it{E}_{MC});#it{E}_{total from per seg cls} [GeV];#it{E}_{MC} [GeV]");
                    objArr->AddAt(hBx_totEFromSegCls_mcE, kBx_totEFromSegCls_mcE);
    TH2F* hBx_mcE_maxClE = new TH2F("hBx_mcE_maxClE","",nBinsE,lowE,uppE,nBinsE,lowE,uppE);
                    hBx_mcE_maxClE->SetTitle("(#it{E}_{MC}, #it{E}_{cl,max});#it{E}_{MC} [GeV];#it{E}_{cl,max} [GeV]");
                    objArr->AddAt(hBx_mcE_maxClE, kBx_mcE_maxClE);
    TH2F* hBx_maxClE_mcE = new TH2F("hBx_maxClE_mcE","",nBinsE,lowE,uppE,nBinsE,lowE,uppE);
                    hBx_maxClE_mcE->SetTitle("(#it{E}_{cl,max}, #it{E}_{MC});#it{E}_{cl,max} [GeV];#it{E}_{MC} [GeV]");
                    objArr->AddAt(hBx_maxClE_mcE, kBx_maxClE_mcE);
    TH2F* hBx_clE_mcE = new TH2F("hBx_clE_mcE","",nBinsE,lowE,uppE,nBinsE,lowE,uppE);
                    hBx_clE_mcE->SetTitle("(#it{E}_{cl}, #it{E}_{MC});#it{E}_{cl} [GeV];#it{E}_{MC} [GeV]");
                    objArr->AddAt(hBx_clE_mcE, kBx_clE_mcE);
    return;
}

void DefineHisto_Bx_TPrf(TObjArray* objArr)
{
    objArr->SetOwner();

    TProfile* hBx_mcE_totE_prof = new TProfile("hBx_mcE_totE_prof","",nBinsE,lowE,uppE,lowE,uppE);
                    hBx_mcE_totE_prof->SetTitle("(#it{E}_{MC}, #it{E}_{total});#it{E}_{MC} [GeV];#it{E}_{total} [GeV]");
                    objArr->AddAt(hBx_mcE_totE_prof, kBx_mcE_totE_prof); 
    TProfile* hBx_totE_mcE_prof = new TProfile("hBx_totE_mcE_prof","",nBinsE,lowE,uppE,lowE,uppE);
                    hBx_totE_mcE_prof->SetTitle("(#it{E}_{total}, #it{E}_{MC});#it{E}_{total} [GeV];#it{E}_{MC} [GeV]");
                    objArr->AddAt(hBx_totE_mcE_prof, kBx_totE_mcE_prof);   
    TProfile* hBx_mcE_maxClE_prof = new TProfile("hBx_mcE_maxClE_prof","",nBinsE,lowE,uppE,lowE,uppE);
                    hBx_mcE_maxClE_prof->SetTitle("(#it{E}_{MC}, #it{E}_{cl,max});#it{E}_{MC} [GeV];#it{E}_{cl,max} [GeV]");
                    objArr->AddAt(hBx_mcE_maxClE_prof, kBx_mcE_maxClE_prof);   
    TProfile* hBx_maxClE_mcE_prof = new TProfile("hBx_maxClE_mcE_prof","",nBinsE,lowE,uppE,lowE,uppE);
                    hBx_maxClE_mcE_prof->SetTitle("(#it{E}_{cl,max}, #it{E}_{MC});#it{E}_{cl,max} [GeV];#it{E}_{MC} [GeV]");
                    objArr->AddAt(hBx_maxClE_mcE_prof, kBx_maxClE_mcE_prof);
    return;
}

// ******************************************************************************************************************
// for analyses of J/psi events:
// ******************************************************************************************************************

/*
    TH2F* hMismatchXY = CreateTH2F("hMismatchXY", 
                    "XY distance between a cluster and a matched primary J/#psi electron vs matched energy;#it{E} (matched primary electron) [GeV];#sqrt{(#Delta#it{x})^{2} + (#Delta#it{y})^{2}} [cm]",
                    nBinsE,lowE,uppE,40,0.,8.);
*/

void DefineHisto_Jp_TH1F(TObjArray* objArr)
{
    objArr->SetOwner();
    TH1F* hJp_mcJPt = new TH1F("hJp_mcJPt","",nBinsPt,lowPt,uppPt);
                    hJp_mcJPt->SetTitle("J/#psi #it{p}_{T} generated;#it{p}_{T,J/#psi} [GeV/#it{c}];counts");
                    objArr->AddAt(hJp_mcJPt, kJp_mcJPt);
    TH1F* hJp_clPairM = new TH1F("hJp_clPairM","",nBinsM,lowM,uppM);
                    hJp_clPairM->SetTitle("inv. mass of cluster pairs;#it{m}_{cl pairs} [GeV/#it{c}^{2}];counts");
                    objArr->AddAt(hJp_clPairM, kJp_clPairM);
    TH1F* hJp_clPairEta = new TH1F("hJp_clPairEta","",nBinsEta,lowEta,uppEta);
                    hJp_clPairEta->SetTitle("#eta of cluster pairs;#eta_{cl pairs} [-];counts");
                    objArr->AddAt(hJp_clPairEta, kJp_clPairEta);
    TH1F* hJp_clPairPt = new TH1F("hJp_clPairPt","",nBinsPt,lowPt,uppPt);
                    hJp_clPairPt->SetTitle("#it{p}_{T} of cluster pairs;#it{p}_{T,cl pairs} [GeV/#it{c}];counts");
                    objArr->AddAt(hJp_clPairPt, kJp_clPairPt);
    TH1F* hJp_clPairPt_massCut = new TH1F("hJp_clPairPt_massCut","",nBinsPt,lowPt,uppPt);
                    hJp_clPairPt_massCut->SetTitle("#it{p}_{T} of cluster pairs having inv. mass above 2.5 GeV/#it{c}^{2};#it{p}_{T,cl pairs} [GeV/#it{c}];counts");
                    objArr->AddAt(hJp_clPairPt_massCut, kJp_clPairPt_massCut);
    TH1F* hJp_primElClPairM = new TH1F("hJp_primElClPairM","",nBinsM,lowM,uppM);
                    hJp_primElClPairM->SetTitle("inv. mass of cluster pairs matched with a pair of pp electrons;#it{m}_{cl pairs matched} [GeV/#it{c}^{2}];counts");
                    objArr->AddAt(hJp_primElClPairM, kJp_primElClPairM);
    TH1F* hJp_primElClPairM_dir = new TH1F("hJp_primElClPairM_dir","",nBinsM,lowM,uppM);
                    hJp_primElClPairM_dir->SetTitle("inv. mass of cluster pairs directly matched with a pair of pp electrons;#it{m}_{cl pairs dir. matched} [GeV/#it{c}^{2}];counts");
                    objArr->AddAt(hJp_primElClPairM_dir, kJp_primElClPairM_dir);
    /*
    TH1F* hJp_primElClPairEta = new TH1F("hJp_primElClPairEta","",nBinsEta,lowEta,uppEta);
                    hJp_primElClPairEta->SetTitle("#eta of cluster pairs matched with a pair of primary J/#psi electrons;#eta [-];counts");
                    objArr->AddAt(hJp_primElClPairEta, kJp_primElClPairEta); 
    TH1F* hJp_primElClPairPt = new TH1F("hJp_primElClPairPt","",nBinsPt,lowPt,uppPt);
                    hJp_primElClPairPt->SetTitle("#it{p}_{T} of cluster pairs matched with a pair of primary J/#psi electrons;#it{p}_{T} [GeV/#it{c}];counts");
                    objArr->AddAt(hJp_primElClPairPt, kJp_primElClPairPt); 
    TH1F* hJp_mcJEEall = new TH1F("hJp_mcJEEall","",nBinsE,lowE,uppE);
                    hJp_mcJEEall->SetTitle("energy of primary J/#psi electron tracks: all;#it{E} [GeV];counts");
                    objArr->AddAt(hJp_mcJEEall, kJp_mcJEEall); 
    TH1F* hJp_mcJEEmtch = new TH1F("hJp_mcJEEmtch","",nBinsE,lowE,uppE);
                    hJp_mcJEEmtch->SetTitle("energy of primary J/#psi electron tracks: only matched;#it{E} [GeV];counts");
                    objArr->AddAt(hJp_mcJEEmtch, kJp_mcJEEmtch); 
    TH1F* hJp_mcJEEratio = new TH1F("hJp_mcJEEratio","",nBinsE,lowE,uppE);
                    hJp_mcJEEratio->SetTitle("energy of primary J/#psi electron tracks: ratio matched/all;#it{E} [GeV];ratio matched/all [-]");
                    objArr->AddAt(hJp_mcJEEratio, kJp_mcJEEratio); 
    // matching to physical primary particles (direct vs finding physical primary mother of matched track of arbitrary type)
    //TH1F* hJp_mcJE_
    */
    return;
}

void DefineHisto_Jp_TH2F(TObjArray* objArr)
{
    objArr->SetOwner();
    TH2F* hJp_mcJRap_mcJPt = new TH2F("hJp_mcJRap_mcJPt","",nBinsEta,lowEta,uppEta,nBinsPt,lowPt,uppPt);
                    hJp_mcJRap_mcJPt->SetTitle("(#it{y}, #it{p}_{T}) of generated J/#psi;#it{y}_{J/#psi} [-];#it{p}_{T,J/#psi} [GeV/#it{c}]");
                    objArr->AddAt(hJp_mcJRap_mcJPt, kJp_mcJRap_mcJPt);
    /*
    TH2F* hJp_mcJRap_mcJPt_acc = new TH2F("hJp_mcJRap_mcJPt_acc","",nBinsEta,lowEta,uppEta,nBinsPt,lowPt,uppPt);
                    hJp_mcJRap_mcJPt_acc->SetTitle("(#it{y}, #it{p}_{T}) of generated J/#psi with electrons within FOCAL acceptance;#it{y} [-];#it{p}_{T} [GeV/#it{c}]");
                    objArr->AddAt(hJp_mcJRap_mcJPt_acc, kJp_mcJRap_mcJPt_acc);
    */
    TH2F* hJp_mcJEEta_mcJEPt = new TH2F("hJp_mcJEEta_mcJEPt","",nBinsEta,lowEta,uppEta,nBinsPt,lowPt,uppPt);
                    hJp_mcJEEta_mcJEPt->SetTitle("(#eta, #it{p}_{T}) of pp electron;#eta_{ppe^{#pm}} [-];#it{p}_{T,ppe^{#pm}} [GeV/#it{c}]");
                    objArr->AddAt(hJp_mcJEEta_mcJEPt, kJp_mcJEEta_mcJEPt);
    TH2F* hJp_mcJE1E_mcJE2E = new TH2F("hJp_mcJE1E_mcJE2E","",nBinsE,lowE,uppE,nBinsE,lowE,uppE);
                    hJp_mcJE1E_mcJE2E->SetTitle("correlation between energy of two pp electrons;#it{E}_{ppe^{#pm},1} [GeV];#it{E}_{ppe^{#pm},2} [GeV]");
                    objArr->AddAt(hJp_mcJE1E_mcJE2E, kJp_mcJE1E_mcJE2E);
    TH2F* hJp_mcJE1Pt_mcJE2Pt = new TH2F("hJp_mcJE1Pt_mcJE2Pt","",nBinsPt,lowPt,uppPt,nBinsPt,lowPt,uppPt);
                    hJp_mcJE1Pt_mcJE2Pt->SetTitle("correlation between transverse momentum of two pp electrons;#it{p}_{T,ppe^{#pm},1} [GeV/#it{c}];#it{p}_{T,ppe^{#pm},2} [GeV/#it{c}]");
                    objArr->AddAt(hJp_mcJE1Pt_mcJE2Pt, kJp_mcJE1Pt_mcJE2Pt);
    // in the following, a cluster means a prefiltered cluster
    // correlation of cluster kinematic variables
    TH2F* hJp_clEta_clPhi = new TH2F("hJp_clEta_clPhi","",nBinsEta,lowEta,uppEta,nBinsPhi,lowPhi,uppPhi);
                    hJp_clEta_clPhi->SetTitle("(#eta, #phi) of pref. clusters;#eta_{cl} [-];#phi_{cl} [-]");
                    objArr->AddAt(hJp_clEta_clPhi, kJp_clEta_clPhi);
    TH2F* hJp_clEta_clPt = new TH2F("hJp_clEta_clPt","",nBinsEta,lowEta,uppEta,nBinsPt,lowPt,uppPt);
                    hJp_clEta_clPt->SetTitle("(#eta, #it{p}_{T}) of pref. clusters;#eta_{cl} [-];#it{p}_{T,cl} [GeV/#it{c}]");
                    objArr->AddAt(hJp_clEta_clPt, kJp_clEta_clPt);
    // cluster kinematics vs matched kinematics (ppp)
    TH2F* hJp_clE_mtchE = new TH2F("hJp_clE_mtchE","",nBinsE,lowE,uppE,nBinsE,lowE,uppE); 
                    hJp_clE_mtchE->SetTitle("cl. energy vs energy of matched pp particle;#it{E}_{cl} [GeV];#it{E}_{matched ppp} [GeV]");
                    objArr->AddAt(hJp_clE_mtchE, kJp_clE_mtchE);
    TH2F* hJp_clE_mtchDirE = new TH2F("hJp_clE_mtchDirE","",nBinsE,lowE,uppE,nBinsE,lowE,uppE); 
                    hJp_clE_mtchDirE->SetTitle("cl. energy vs energy of directly matched pp particle;#it{E}_{cl} [GeV];#it{E}_{dir. matched ppp} [GeV]");
                    objArr->AddAt(hJp_clE_mtchDirE, kJp_clE_mtchDirE);
    TH2F* hJp_clEta_mtchEta = new TH2F("hJp_clEta_mtchEta","",nBinsEta,lowEta,uppEta,nBinsEta,lowEta,uppEta);
                    hJp_clEta_mtchEta->SetTitle("cl. #eta vs #eta of matched pp particle;#eta_{cl} [-];#eta_{matched ppp} [-]");
                    objArr->AddAt(hJp_clEta_mtchEta, kJp_clEta_mtchEta);
    TH2F* hJp_clEta_mtchDirEta = new TH2F("hJp_clEta_mtchDirEta","",nBinsEta,lowEta,uppEta,nBinsEta,lowEta,uppEta);
                    hJp_clEta_mtchDirEta->SetTitle("cl. #eta vs #eta of directly matched pp particle;#eta_{cl} [-];#eta_{dir. matched ppp} [-]");
                    objArr->AddAt(hJp_clEta_mtchDirEta, kJp_clEta_mtchDirEta);
    TH2F* hJp_clPt_mtchPt = new TH2F("hJp_clPt_mtchPt","",nBinsPt,lowPt,uppPt,nBinsPt,lowPt,uppPt);
                    hJp_clPt_mtchPt->SetTitle("cl. #it{p}_{T} vs #it{p}_{T} of matched pp particle;#it{p}_{T,cl} [GeV/#it{c}];#it{p}_{T,matched ppp} [GeV/#it{c}]");
                    objArr->AddAt(hJp_clPt_mtchPt, kJp_clPt_mtchPt);
    TH2F* hJp_clPt_mtchDirPt = new TH2F("hJp_clPt_mtchDirPt","",nBinsPt,lowPt,uppPt,nBinsPt,lowPt,uppPt);
                    hJp_clPt_mtchDirPt->SetTitle("cl. #it{p}_{T} vs #it{p}_{T} of directly matched pp particle;#it{p}_{T,cl} [GeV/#it{c}];#it{p}_{T,dir. matched ppp} [GeV/#it{c}]");
                    objArr->AddAt(hJp_clPt_mtchDirPt, kJp_clPt_mtchDirPt);
    // cluster kinematics vs matched kinematics (ppe only)
    TH2F* hJp_primElClE_mtchE = new TH2F("hJp_primElClE_mtchE","",nBinsE,lowE,uppE,nBinsE,lowE,uppE); 
                    hJp_primElClE_mtchE->SetTitle("cl. energy vs energy of matched pp electron;#it{E}_{cl} [GeV];#it{E}_{matched ppe^{#pm}} [GeV]");
                    objArr->AddAt(hJp_primElClE_mtchE, kJp_primElClE_mtchE);
    TH2F* hJp_primElClE_mtchDirE = new TH2F("hJp_primElClE_mtchDirE","",nBinsE,lowE,uppE,nBinsE,lowE,uppE); 
                    hJp_primElClE_mtchDirE->SetTitle("cl. energy vs energy of directly matched pp electron;#it{E}_{cl} [GeV];#it{E}_{dir. matched ppe^{#pm}} [GeV]");
                    objArr->AddAt(hJp_primElClE_mtchDirE, kJp_primElClE_mtchDirE);
    TH2F* hJp_primElClEta_mtchEta = new TH2F("hJp_primElClEta_mtchEta","",nBinsEta,lowEta,uppEta,nBinsEta,lowEta,uppEta);
                    hJp_primElClEta_mtchEta->SetTitle("cl. #eta vs #eta of matched pp electron;#eta_{cl} [-];#eta_{matched ppe^{#pm}} [-]");
                    objArr->AddAt(hJp_primElClEta_mtchEta, kJp_primElClEta_mtchEta);
    TH2F* hJp_primElClEta_mtchDirEta = new TH2F("hJp_primElClEta_mtchDirEta","",nBinsEta,lowEta,uppEta,nBinsEta,lowEta,uppEta);
                    hJp_primElClEta_mtchDirEta->SetTitle("cl. #eta vs #eta of directly matched pp electron;#eta_{cl} [-];#eta_{dir. matched ppe^{#pm}} [-]");
                    objArr->AddAt(hJp_primElClEta_mtchDirEta, kJp_primElClEta_mtchDirEta);
    TH2F* hJp_primElClPt_mtchPt = new TH2F("hJp_primElClPt_mtchPt","",nBinsPt,lowPt,uppPt,nBinsPt,lowPt,uppPt);
                    hJp_primElClPt_mtchPt->SetTitle("cl. #it{p}_{T} vs #it{p}_{T} of matched pp electron;#it{p}_{T,cl} [GeV/#it{c}];#it{p}_{T,matched ppe^{#pm}} [GeV/#it{c}]");
                    objArr->AddAt(hJp_primElClPt_mtchPt, kJp_primElClPt_mtchPt);
    TH2F* hJp_primElClPt_mtchDirPt = new TH2F("hJp_primElClPt_mtchDirPt","",nBinsPt,lowPt,uppPt,nBinsPt,lowPt,uppPt);
                    hJp_primElClPt_mtchDirPt->SetTitle("cl. #it{p}_{T} vs #it{p}_{T} of directly matched pp electron;#it{p}_{T,cl} [GeV/#it{c}];#it{p}_{T,dir. matched ppe^{#pm}} [GeV/#it{c}]");
                    objArr->AddAt(hJp_primElClPt_mtchDirPt, kJp_primElClPt_mtchDirPt);
    /*
    TH2F* hJp_primElClPairEta_mcJEta = new TH2F("hJp_primElClPairEta_mcJEta","",nBinsEta,lowEta,10.,nBinsEta,lowEta,10.);
                    hJp_primElClPairEta_mcJEta->SetTitle("#eta of a cluster pair matched with a pair of primary J/#psi electrons vs #eta of generated J/#psi;#eta (cluster pair) [-];#eta (matched primary electron pair) [-]");
                    objArr->AddAt(hJp_primElClPairEta_mcJEta, kJp_primElClPairEta_mcJEta);
    TH2F* hJp_primElClPairPt_mcJPt = new TH2F("hJp_primElClPairPt_mcJPt","",nBinsPt,lowPt,uppPt,nBinsPt,lowPt,uppPt);
                    hJp_primElClPairPt_mcJPt->SetTitle("#it{p}_{T} of a cluster pair matched with a pair of primary J/#psi electrons vs #it{p}_{T} of generated J/#psi;#it{p}_{T} (cluster pair) [GeV/#it{c}];#it{p}_{T} (matched primary electron pair) [GeV/#it{c}]");
                    objArr->AddAt(hJp_primElClPairPt_mcJPt, kJp_primElClPairPt_mcJPt);
    */

    // matching to physical primary particles (direct vs finding physical primary mother of matched track of arbitrary type)
    TH2F* hJp_mtchE_mtchDirE = new TH2F("hJp_mtchE_mtchDirE","",nBinsE,lowE,uppE,nBinsE,lowE,uppE); 
                    hJp_mtchE_mtchDirE->SetTitle("energy of matched pp particle vs energy of directly matched ppp;#it{E}_{matched ppp} [GeV];#it{E}_{dir. matched ppp} [GeV]");
                    objArr->AddAt(hJp_mtchE_mtchDirE, kJp_mtchE_mtchDirE);
    TH2F* hJp_mtchE_mtchDirE_primEl = new TH2F("hJp_mtchE_mtchDirE_primEl","",nBinsE,lowE,uppE,nBinsE,lowE,uppE); 
                    hJp_mtchE_mtchDirE_primEl->SetTitle("energy of matched pp electron vs energy of directly matched ppe;#it{E}_{matched ppe^{#pm}} [GeV];#it{E}_{dir. matched ppe^{#pm}} [GeV]");
                    objArr->AddAt(hJp_mtchE_mtchDirE_primEl, kJp_mtchE_mtchDirE_primEl);
    TH2F* hJp_mtchEta_mtchDirEta_primEl = new TH2F("hJp_mtchEta_mtchDirEta_primEl","",nBinsEta,lowEta,uppEta,nBinsEta,lowEta,uppEta); 
                    hJp_mtchEta_mtchDirEta_primEl->SetTitle("#eta of matched pp electron vs #eta of directly matched ppe;#it{eta}_{matched ppe^{#pm}} [GeV];#it{eta}_{dir. matched ppe^{#pm}} [GeV]");
                    objArr->AddAt(hJp_mtchEta_mtchDirEta_primEl, kJp_mtchEta_mtchDirEta_primEl);
    TH2F* hJp_mtchPt_mtchDirPt_primEl = new TH2F("hJp_mtchPt_mtchDirPt_primEl","",nBinsPt,lowPt,uppPt,nBinsPt,lowPt,uppPt); 
                    hJp_mtchPt_mtchDirPt_primEl->SetTitle("#it{p}_{T} of matched pp electron vs #it{p}_{T} of directly matched ppe;#it{p}_{T,matched ppe^{#pm}} [GeV];#it{p}_{T,dir. matched ppe^{#pm}} [GeV]");
                    objArr->AddAt(hJp_mtchPt_mtchDirPt_primEl, kJp_mtchPt_mtchDirPt_primEl);
    return;
}

void DefineHisto_Jp_TPrf(TObjArray* objArr)
{
    objArr->SetOwner();

    return;
}

// ******************************************************************************************************************
// Functions to plot 1d and 2d histograms 
// ******************************************************************************************************************

template <typename TH> // for TH1 and TProfile
void DrawHisto(TH* h, TString subfolder)
{
    TCanvas c("c","c",700,600);
    h->GetYaxis()->SetMaxDigits(3);
    h->GetXaxis()->SetTitleOffset(1.2);
    h->SetLineColor(kBlue+1);
    h->SetFillColor(kBlue);
    h->SetFillStyle(3012);
    TString path_out = subfolder + h->GetName() + ".pdf";
    c.cd();
    h->Draw();
    c.Print(path_out.Data());
    return;
}

void DrawHistoCOLZ(TH2F* h, TString subfolder) // for TH2F
{
    TCanvas c("c","c",700,600);
    c.SetGrid();
    //c.SetLogz();
    h->GetYaxis()->SetMaxDigits(3);
    h->GetXaxis()->SetTitleOffset(1.2);
    TString path_out = subfolder + h->GetName() + ".pdf";
    c.cd();
    h->Draw("COLZ");
    c.Print(path_out.Data());
    return;
}