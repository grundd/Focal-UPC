// FocalUpcGrid.C
// David Grund, Nov 05, 2022

#include "ConfigParameters.h"
#include "FocalUpcGrid.h"
#include "CreateHistograms.h"

void FocalUpcGrid(Bool_t isLocal, TString sim, Bool_t overwrite = kTRUE, TString sIn = "", TString sOut = "")
{
    // PDG code of expected mother of physical primary electrons in the dataset
    // (!) in AliDPG notation in feed-down files, the mother of electron pairs is psi', not J/psi
    // not needed for box simulations
    Int_t pdgMother(-1);
    // PDG code of the main vector meson in the dataset
    Int_t pdgMainVM(-1);
    if(sim == "cohJpsi" 
    || sim == "cohJpsiNoFIT"
    || sim == "incJpsi") { pdgMother = 443; pdgMainVM = 443; }
    else if(sim == "cohFD"  
         || sim == "incFD") { pdgMother = 100443; pdgMainVM = 443; }
    else if(sim == "cohPsi2s" 
         || sim == "incPsi2s") { pdgMother = 100443; pdgMainVM = 100443; }

    // is it box simulation?
    Bool_t isBoxSim(kFALSE);
    if(sim == "boxEle" || sim == "boxPho") isBoxSim = kTRUE;

    // mass filtering only for J/psi simulations
    if(cutM > 0 && isBoxSim) {
        cout << " ERROR: Cannot do mass cleaning for box simulations of electrons/photons. Terminating... " << endl;
        // terminate:
        return;
    }

    // prepare subdirectory in output directory
    TString outSubDir = "";
    if(isLocal) 
    {
        outSubDir = CreateOutputSubDir();
        gSystem->Exec("mkdir -p " + sOut + outSubDir);
    }

    // check if focalClusters.root were produced properly
    TString sClFile = Form("%sfocalClusters.root",sOut.Data());
    if(gSystem->AccessPathName(sClFile.Data()))
    {
        cout << " ERROR: cluster file not found! Terminating." << endl;
        // terminate:
        return;
    } 
    // open the file with clusters
    TFile* fCls = new TFile(sClFile.Data());
    cout << " MESSAGE: Loading clusters from: " << fCls->GetName() << endl;

    // define ALICE run loader: open galice.root
    AliRunLoader* runLoader = NULL;
    if(!isLocal) runLoader = AliRunLoader::Open(sIn + "root_archive.zip#galice.root");
    else         runLoader = AliRunLoader::Open(sIn + "galice.root");
    if(!runLoader) 
    {
        cout << " ERROR: AliRunLoader not good! Terminating." << endl;
        // terminate:
        return;   
    }
    if(!runLoader->GetAliRun()) runLoader->LoadgAlice();
    if(!runLoader->TreeE()) runLoader->LoadHeader();
    if(!runLoader->TreeK()) runLoader->LoadKinematics();

    // arrays with output histograms
    TObjArray* arrHistos = NULL;
    Int_t kFirstTH1F(-1), kFirstTH2F(-1), kFirstTPrf(-1), kFirstTP2D(-1), kAll(-1);
    // if box simulations
    if(isBoxSim) 
    {
        kFirstTH1F = kGridBox_firstTH1F; 
        kFirstTH2F = kGridBox_firstTH2F; 
        kFirstTPrf = kGridBox_firstTPrf; 
        kFirstTP2D = 
        kAll = kGridBox_all;
        arrHistos = new TObjArray(kAll);
        CreateHistos_GridBox(arrHistos);
    } 
    // if starlight J/psi simulations
    else
    {
        kFirstTH1F = kGridJpsi_firstTH1F; 
        kFirstTH2F = kGridJpsi_firstTH2F;
        kFirstTPrf = kGridJpsi_firstTPrf; 
        kFirstTP2D = kGridJpsi_firstTP2D;
        kAll = kGridJpsi_all;
        arrHistos = new TObjArray(kAll);
        CreateHistos_GridJpsi(arrHistos);
    }

    // output file
    TString sFile = Form("%s%sanalysisResultsGrid.root",sOut.Data(),outSubDir.Data());
    TFile* fOut = NULL;
    // local analysis: if the output .root file has already been produced,
    // this dataset won't be analyzed again if the option 'overwrite' is set to kFALSE
    if(isLocal == kTRUE && overwrite == kFALSE) {
        cout << sFile << ":" << endl;
        // output file not found:
        if(gSystem->AccessPathName(sFile.Data())) cout << " MESSAGE: output file not found! Runing the analysis now:" << endl;
        // output file found:
        else {
            cout << " MESSAGE: output file found. Using previous results..." << endl;
            // terminate:
            runLoader->Delete();
            return;
        }
    }
    
    fOut = new TFile(sFile.Data(),"RECREATE");
    // output tree
    TTree* tOut = new TTree("tCls", "output tree containing prefiltered clusters");
    // prefiltered clusters
    Int_t fEvNumber;
    Float_t fEnCl, fXCl, fYCl, fZCl;
    tOut->Branch("fEvNumber", &fEvNumber, "fEvNumber/I"); 
    tOut->Branch("fEnCl", &fEnCl, "fEnCl/F");
    tOut->Branch("fXCl", &fXCl, "fXCl/F");
    tOut->Branch("fYCl", &fYCl, "fYCl/F");
    tOut->Branch("fZCl", &fZCl, "fZCl/F");
    // J/psi electrons with which the clusters were matched
    Int_t fIdxJEl;
    tOut->Branch("fIdxJEl", &fIdxJEl, "fIdxJEl/I");
    TParticle* fJEl = new TParticle();
    tOut->Branch("fJEl","TParticle",&fJEl,32000,-1);
    gROOT->cd();    

    // output log file
    TString sLog = Form("%s%sanalysis.log",sOut.Data(),outSubDir.Data());
    ofstream of(sLog);
    of << "create supercls: " << doSupercls << "\n"
       << std::fixed << std::setprecision(1)
       << "  with min seed energy: " << minSeedE << " GeV\n"
       << "  radius parameter: " << radius << " cm\n"
       << "cut on min cl energy: " << cutE << " GeV\n"
       << "cut on min cl pair mass: " << cutM << " GeV/c^2\n"
       << "match directly: " << matchDirectly << "\n"
       << "**************************************\n\n";

    // loop over MC events:
    Int_t nEvents = runLoader->GetNumberOfEvents();
    of << "analyzing " << nEvents << " events:" << endl;
    // progress bar:
    Float_t progress = 0.; // perc
    for(Int_t iEv = 0; iEv < nEvents; iEv++) 
    {
        of << "Ev " << iEv+1 << ":" << endl;
        // update the progress bar
        if((iEv+1) % (Int_t)(nEvents/10.) == 0) {
            progress += 10.;
            cout << "[" << progress << "%] done." << endl;
        }

        // get current MC event
        Int_t isEventOk = runLoader->GetEvent(iEv);
        // the method GetEvent(i) returns zero if event is loaded succesfully, if not we continue
        if (isEventOk != 0) {
            of << "Ev " << iEv+1 << " not OK, skipping." << endl;
            continue;
        }

        // get the stack of MC tracks for this event
        AliStack* stack = runLoader->Stack();

        // get a tree with clusters for this event 
        // (separate tree for each event in subfolders "Event0, Event1, ...")
        TTree* tCls = NULL;
        if(fCls->GetDirectory(Form("Event%i",iEv))) fCls->GetDirectory(Form("Event%i",iEv))->GetObject("fTreeR",tCls);
        else {
            of << "  (!) cannot find Ev " << iEv+1 << " in a cluster file " << fCls->GetName() << ". Skipping..." << endl;
            fCls->ls();
            continue;
        }

        // get the branch with final summed clusters
        TBranch* bClsSum = NULL;
        bClsSum = tCls->GetBranch("AliFOCALCluster");
        TClonesArray* arrClsSum = NULL;
        bClsSum->SetAddress(&arrClsSum);
        bClsSum->GetEvent(0);
        Int_t nClsSum = arrClsSum->GetEntries();
        of << "  summed cls: " << nClsSum << endl;

        // prepare a list of prefiltered clusters 
        // for the moment, add all clusters (later we will apply the selections)
        TList* listClsPref = new TList();
        listClsPref->SetOwner(kFALSE);
        // find the total and maximum cluster energy in this event
        Float_t ETot(0.), EClMax(-1.);
        Int_t iClMaxE(-1); // index of a cluster with maximum energy
        for(Int_t iCl = 0; iCl < nClsSum; iCl++) 
        {
            AliFOCALCluster* clust = (AliFOCALCluster*) arrClsSum->At(iCl);
            listClsPref->Add(clust);
            // get the total energy summing over summed clusters
            ETot += clust->E();
            // get the maximum cluster energy
            if(clust->E() > EClMax) {
                EClMax = clust->E();
                iClMaxE = iCl;
            } 
        }
        Int_t nClsPref = listClsPref->GetEntries();

        // prepare superclusters
        if(doSupercls)
        {
            // create new list: listClsPrefSort
            TList listClsPrefSort;
            listClsPrefSort.SetOwner(kFALSE);
            // fill it with the clusters from listClsPref, sorting them according to their energy
            while(listClsPrefSort.GetEntries() < nClsPref)
            {
                Float_t EClMaxFound = -1.;
                AliFOCALCluster* clustMaxE = NULL;
                for(Int_t iCl = 0; iCl < listClsPref->GetEntries(); iCl++) 
                {
                    AliFOCALCluster* clust = (AliFOCALCluster*)listClsPref->At(iCl);
                    Float_t ECl = clust->E();
                    if(ECl > EClMaxFound) {
                        EClMaxFound = ECl;
                        clustMaxE = clust;
                    }
                }
                listClsPrefSort.Add(clustMaxE);
                listClsPref->Remove(clustMaxE);
            }
            // while the energy of the most energetic cluster in listClsPrefSort is above minSeedE, 
            // take this cluster as a seed for a new supercluster
            // then combine it with all other clusters lying within the predefined radius
            AliFOCALCluster* clustTop = NULL;
            Float_t EClTop(-1.);
            if(nClsPref > 0) {
                clustTop = (AliFOCALCluster*)listClsPrefSort.At(0);
                EClTop = clustTop->E();
            }
            while(EClTop > minSeedE)
            {
                Float_t xSeed = clustTop->X();
                Float_t ySeed = clustTop->Y();
                Float_t xSupCl = clustTop->X() * EClTop;
                Float_t ySupCl = clustTop->Y() * EClTop;
                Float_t zSupCl = clustTop->Z() * EClTop;
                Float_t ESupCl = EClTop;
                listClsPrefSort.Remove(clustTop);
                for(Int_t iCl = listClsPrefSort.GetEntries()-1; iCl >= 0; iCl--) 
                {
                    AliFOCALCluster* clust = (AliFOCALCluster*)listClsPrefSort.At(iCl);
                    Float_t xCl = clust->X();
                    Float_t yCl = clust->Y();
                    Float_t zCl = clust->Z();
                    Float_t ECl = clust->E();
                    Float_t distXY = TMath::Sqrt(TMath::Power(xSeed-xCl,2) + TMath::Power(ySeed-yCl,2));
                    // if this cluster is close enough to the supercluster seed, add it to the supercluster
                    if(distXY < radius) {
                        xSupCl += xCl * ECl;
                        ySupCl += yCl * ECl;
                        zSupCl += zCl * ECl;
                        ESupCl += ECl;
                        listClsPrefSort.Remove(clust);
                    }
                }
                if(listClsPrefSort.GetEntries() > 0) {
                    clustTop = (AliFOCALCluster*)listClsPrefSort.At(0);
                    EClTop = clustTop->E();
                } 
                else EClTop = 0;
                // add the new supercluster to listClsPref if abs(x_cl) > 5 cm and abs(y_cl) > 5 cm
                AliFOCALCluster* supClNew = new AliFOCALCluster(xSupCl/ESupCl,ySupCl/ESupCl,zSupCl/ESupCl,ESupCl,-1);
                if(TMath::Abs(xSupCl/ESupCl) > 5 || TMath::Abs(ySupCl/ESupCl) > 5) listClsPref->Add(supClNew);
            }
            nClsPref = listClsPref->GetEntries();
            of << "  identified supercls: " << nClsPref << endl;
        }

        // if J/psi analysis, do "mass cleaning"
        // i.e., remove all clusters that have a mass below cutM when pared with any other cluster
        if(cutM > 0)
        {
            std::vector<Bool_t> foundLowMass;
            for(Int_t i = 0; i < nClsPref; i++) foundLowMass.push_back(kFALSE);
            // go over the clusters
            for(Int_t iCl1 = 0; iCl1 < nClsPref; iCl1++) 
            {            
                AliFOCALCluster* clust1 = (AliFOCALCluster*)listClsPref->At(iCl1);
                if(!clust1) continue;
                // get energy and coordinates of this cluster
                Float_t xCl1 = clust1->X();
                Float_t yCl1 = clust1->Y();
                Float_t zCl1 = clust1->Z();
                Float_t ECl1 = clust1->E();
                TLorentzVector cl1 = ConvertXYZEtoLorVec(xCl1,yCl1,zCl1,ECl1);
                // go over all possible pairs and combine clusters
                for(Int_t iCl2 = iCl1+1; iCl2 < nClsPref; iCl2++) 
                {
                    AliFOCALCluster* clust2 = (AliFOCALCluster*)listClsPref->At(iCl2);
                    if(!clust2) continue;
                    // get energy and coordinates of this cluster
                    Float_t xCl2 = clust2->X();
                    Float_t yCl2 = clust2->Y();
                    Float_t zCl2 = clust2->Z();
                    Float_t ECl2 = clust2->E();
                    TLorentzVector cl2 = ConvertXYZEtoLorVec(xCl2,yCl2,zCl2,ECl2);
                    // calculate the mass of the two current clusters
                    TLorentzVector cl12 = cl1 + cl2;
                    Float_t mass = cl12.M();
                    if (mass < cutM) 
                    {
                        foundLowMass[iCl1] = kTRUE;
                        foundLowMass[iCl2] = kTRUE;
                        // no need to look for any more pairs for cluster at iCl1
                        break;
                    }
                }
            }
            // remove all clusters giving low mass when paired
            for(Int_t iCl = nClsPref-1; iCl >= 0; iCl--)
            {
                AliFOCALCluster* clust = (AliFOCALCluster*)listClsPref->At(iCl);
                if(foundLowMass[iCl]) listClsPref->Remove(clust);
            }
            // how many clusters there are after prefiltering
            nClsPref = listClsPref->GetEntries();
            of << "  after mass cleaning: " << nClsPref << endl;
        }

        // cut on minimum cluster energy
        // remove all remaining clusters (if any) with energy below threshold
        if(cutE > 0)
        {
            for(Int_t iCl = nClsPref-1; iCl >= 0; iCl--) 
            {
                AliFOCALCluster* clust = (AliFOCALCluster*)listClsPref->At(iCl);
                if(clust->E() < cutE) listClsPref->Remove(clust);
            }
            nClsPref = listClsPref->GetEntries();
            of << "  after energy cut: " << nClsPref << endl;
        }

        // ******************************************************************************************************************
        // analysis of box electron or photon simulations
        // ******************************************************************************************************************
        if(isBoxSim)
        {
            // general histograms
            ((TH1F*)arrHistos->At(kB1_mcE))->Fill(stack->Particle(0)->Energy());
            ((TH1F*)arrHistos->At(kB1_mcPt))->Fill(stack->Particle(0)->Pt());
            ((TH2F*)arrHistos->At(kB2_mcE_nCls))->Fill(stack->Particle(0)->Energy(), nClsPref);
            // total energy vs MC energy
            ((TH2F*)arrHistos->At(kB2_totE_mcE))->Fill(ETot, stack->Particle(0)->Energy());
            ((TProfile*)arrHistos->At(kBP_totE_mcE))->Fill(ETot, stack->Particle(0)->Energy());
            // maximum cluster energy vs MC energy
            ((TH2F*)arrHistos->At(kB2_maxClE_mcE))->Fill(EClMax, stack->Particle(0)->Energy());
            ((TProfile*)arrHistos->At(kBP_maxClE_mcE))->Fill(EClMax, stack->Particle(0)->Energy());
            // X and Y of the cluster with the highest energy
            Float_t xClMaxE(-1.), yClMaxE(-1.);
            Bool_t isClWithMaxE = kFALSE;
            if(iClMaxE != -1 && cutM == 0 && cutE == 0 && !doSupercls) {
                isClWithMaxE = kTRUE;
                AliFOCALCluster* clustMaxE = (AliFOCALCluster*) listClsPref->At(iClMaxE);
                xClMaxE = clustMaxE->X();
                yClMaxE = clustMaxE->Y();
            }
            // info per prefiltered (super)clusters
            for(Int_t iCl = 0; iCl < nClsPref; iCl++)
            {
                AliFOCALCluster* clust = (AliFOCALCluster*) listClsPref->At(iCl);
                Float_t xCl = clust->X();
                Float_t yCl = clust->Y();
                Float_t zCl = clust->Z();
                // physical primary electron:
                TParticle* part = stack->Particle(0);
                Float_t x(0.), y(0.);
                TrackCoordinatesAtZ(part,zCl,x,y);
                // correlation between MC energy and energies of prefiltered clusters
                ((TH2F*)arrHistos->At(kB2_clE_mcE))->Fill(clust->E(), part->Energy());
                // radial distance between the cluster and the MC track
                Float_t DeltaX = xCl - x;
                Float_t DeltaY = yCl - y;
                Float_t DeltaR = TMath::Sqrt(TMath::Power(DeltaX,2) + TMath::Power(DeltaY,2));
                ((TH2F*)arrHistos->At(kB2_clMcDX_clMcDY))->Fill(DeltaX, DeltaY);
                ((TH2F*)arrHistos->At(kB2_mcE_clMcSep))->Fill(part->Energy(), DeltaR);
                // radial distance between the cluster and the cluster with maximum energy
                if(isClWithMaxE && (iCl != iClMaxE)) {
                    Float_t DeltaX_maxE = xCl - xClMaxE;
                    Float_t DeltaY_maxE = yCl - yClMaxE;
                    Float_t DeltaR_maxE = TMath::Sqrt(TMath::Power(DeltaX_maxE,2) + TMath::Power(DeltaY_maxE,2));
                    ((TH2F*)arrHistos->At(kB2_clMaxClDX_clMaxClDY))->Fill(DeltaX_maxE, DeltaY_maxE);
                    ((TH2F*)arrHistos->At(kB2_mcE_clMaxClSep))->Fill(part->Energy(), DeltaR_maxE);
                    ((TH2F*)arrHistos->At(kB2_clX_clMaxClSep))->Fill((xCl + xClMaxE)/2., DeltaR_maxE);
                    ((TH2F*)arrHistos->At(kB2_clY_clMaxClSep))->Fill((yCl + yClMaxE)/2., DeltaR_maxE);
                }
            }
        }
        // ******************************************************************************************************************
        // analysis of starlight J/psi simulations
        // ******************************************************************************************************************
        else
        {
            // fill histograms with MC kinematics
            TObjArray ppElectrons; // array of physical primary electrons
            for(Int_t iTrk = 0; iTrk < stack->GetNtrack(); iTrk++)
            {
                TParticle* part = stack->Particle(iTrk);
                // if J/psi
                if(part->GetPdgCode() == pdgMainVM) {
                    ((TH1F*)arrHistos->At(kJ1_mcJEn))->Fill(part->Energy());
                    ((TH1F*)arrHistos->At(kJ1_mcJPt))->Fill(part->Pt());
                    ((TH1F*)arrHistos->At(kJ1_mcJRap))->Fill(part->Y());
                    ((TH1F*)arrHistos->At(kJ1_mcJM))->Fill(part->GetCalcMass());
                    ((TH2F*)arrHistos->At(kJ2_mcJRap_mcJPt))->Fill(part->Y(), part->Pt());
                }
                // if physical primary (pp) electron (i.e., J/psi direct decay product)
                if(isEleOrPos(part) && stack->IsPhysicalPrimary(iTrk)) {
                    // if a mother particle exists:
                    if(part->GetMother(0) >= 0) {
                        TParticle* mother = stack->Particle(part->GetMother(0));
                        // if it doesn't have the expected pdg code:
                        if(mother->GetPdgCode() != pdgMother) 
                            cout << " * MESSAGE: Unexpected mother of a pp electron." << endl;
                    }
                    ((TH2F*)arrHistos->At(kJ2_mcJElEta_mcJElPt))->Fill(part->Eta(), part->Pt());
                    ppElectrons.AddLast(part);
                }
            }
            // kinematics of the two physical primary electrons
            TParticle* ppEl1 = NULL;
            TParticle* ppEl2 = NULL;
            TLorentzVector* lvJElPair = new TLorentzVector();
            if(ppElectrons.GetEntries() != 2) {
                cout << "  (!) ERROR: Unexpected number of phys. prim. electrons: " << ppElectrons.GetEntries() << ", terminating..." << endl;
                // terminate:
                runLoader->Delete();
                fCls->Close();
                delete fCls;
                return;
            } else {
                // fill the histogram showing the correlation between energies and transverse momenta of the two pp electrons
                ppEl1 = (TParticle*) ppElectrons[0];
                ppEl2 = (TParticle*) ppElectrons[1];
                ((TH2F*)arrHistos->At(kJ2_mcJEl1En_mcJEl2En))->Fill(ppEl1->Energy(), ppEl2->Energy());
                ((TH2F*)arrHistos->At(kJ2_mcJEl1Pt_mcJEl2Pt))->Fill(ppEl1->Pt(), ppEl2->Pt());
                // fill the histograms showing kinematics of J/psi reconstructed from the two pp electrons
                TLorentzVector lvPPEl1;
                lvPPEl1.SetPxPyPzE(ppEl1->Px(),ppEl1->Py(),ppEl1->Pz(),ppEl1->Energy());
                TLorentzVector lvPPEl2;
                lvPPEl2.SetPxPyPzE(ppEl2->Px(),ppEl2->Py(),ppEl2->Pz(),ppEl2->Energy());
                TLorentzVector Jpsi = lvPPEl1 + lvPPEl2;
                lvJElPair->SetPxPyPzE(Jpsi.Px(),Jpsi.Py(),Jpsi.Pz(),Jpsi.Energy());
                ((TH1F*)arrHistos->At(kJ1_mcJElPairEn))->Fill(lvJElPair->Energy());
                ((TH1F*)arrHistos->At(kJ1_mcJElPairPt))->Fill(lvJElPair->Pt());
                ((TH1F*)arrHistos->At(kJ1_mcJElPairRap))->Fill(lvJElPair->Rapidity());
                ((TH1F*)arrHistos->At(kJ1_mcJElPairM))->Fill(lvJElPair->M());
                // acceptance
                ((TH1F*)arrHistos->At(kJ1_rap_gen))->Fill(lvJElPair->Rapidity());
                if(3.4 < ppEl1->Eta() && ppEl1->Eta() < 5.8 && 3.4 < ppEl2->Eta() && ppEl2->Eta() < 5.8)
                    ((TH1F*)arrHistos->At(kJ1_rap_acc))->Fill(lvJElPair->Rapidity());
                // number of (sup)cls
                ((TH2F*)arrHistos->At(kJ2_mcJElPairEn_nCls))->Fill(lvJElPair->Energy(),nClsPref);
            }

            // match clusters to MC tracks
            std::vector<Int_t> idxMtchPhysPrimParts; 
            TObjArray arrMtchPhysPrimParts;
            MatchClsToPhysPrimP(stack,listClsPref,idxMtchPhysPrimParts,&arrMtchPhysPrimParts);

            // fill the histograms
            for(Int_t iCl = 0; iCl < nClsPref; iCl++)
            {
                AliFOCALCluster* clust = (AliFOCALCluster*) listClsPref->At(iCl);
                if(!clust) continue;
                // get energy and coordinates of this cluster
                Float_t xCl = clust->X();
                Float_t yCl = clust->Y();
                Float_t zCl = clust->Z();
                Float_t ECl = clust->E();
                TLorentzVector cl = ConvertXYZEtoLorVec(xCl,yCl,zCl,ECl);
                // -> tree:
                fEvNumber = iEv;
                fEnCl = cl.Energy();
                fXCl = xCl;
                fYCl = yCl;
                fZCl = zCl;
                fIdxJEl = -1;
                // other cluster kinematics
                Float_t fPtCl = cl.Pt();
                Float_t fEtaCl = cl.Eta();
                Float_t fPhiCl = cl.Phi();

                ((TH1F*)arrHistos->At(kJ1_clZ))->Fill(zCl);
                // fill some histograms with cluster kinematics
                ((TH2F*)arrHistos->At(kJ2_clEta_clPhi))->Fill(fEtaCl, fPhiCl);
                ((TH2F*)arrHistos->At(kJ2_clEta_clPt))->Fill(fEtaCl, fPtCl);
                ((TH2F*)arrHistos->At(kJ2_clX_clY))->Fill(xCl, yCl);
                ((TH2F*)arrHistos->At(kJ2_clX_clEn))->Fill(xCl, ECl);
                ((TH2F*)arrHistos->At(kJ2_clY_clEn))->Fill(yCl, ECl);
                ((TProfile*)arrHistos->At(kJP_clX_clEn))->Fill(xCl, ECl);
                ((TProfile*)arrHistos->At(kJP_clY_clEn))->Fill(yCl, ECl);
                ((TProfile2D*)arrHistos->At(kJP2_clX_clY_clEn))->Fill(xCl, yCl, ECl);

                TParticle* mtchPhysPrimPart = (TParticle*)arrMtchPhysPrimParts[iCl];
                // if matched to a ppp
                if(mtchPhysPrimPart) 
                {
                    ((TH2F*)arrHistos->At(kJ2_pppClEn_mtchEn))->Fill(fEnCl, mtchPhysPrimPart->Energy());
                    ((TH2F*)arrHistos->At(kJ2_pppClEta_mtchEta))->Fill(fEtaCl, mtchPhysPrimPart->Eta());
                    ((TH2F*)arrHistos->At(kJ2_pppClPt_mtchPt))->Fill(fPtCl, mtchPhysPrimPart->Pt());
                    // if electron (ppe)
                    if(isEleOrPos(mtchPhysPrimPart)) {
                        fIdxJEl = idxMtchPhysPrimParts[iCl];
                        Float_t fEnJEl = mtchPhysPrimPart->Energy();
                        Float_t fPtJEl = mtchPhysPrimPart->Pt();
                        Float_t fEtaJEl = mtchPhysPrimPart->Eta();
                        *fJEl = TParticle(*mtchPhysPrimPart);
                        ((TH2F*)arrHistos->At(kJ2_ppeClEn_mtchEn))->Fill(fEnCl, fEnJEl);
                        ((TH2F*)arrHistos->At(kJ2_ppeClEta_mtchEta))->Fill(fEtaCl, fEtaJEl);
                        ((TH2F*)arrHistos->At(kJ2_ppeClPt_mtchPt))->Fill(fPtCl, fPtJEl);
                        ((TProfile*)arrHistos->At(kJP_ppeClX_mtchEn))->Fill(xCl, fEnJEl);
                        ((TProfile*)arrHistos->At(kJP_ppeClY_mtchEn))->Fill(yCl, fEnJEl);
                        ((TProfile2D*)arrHistos->At(kJP2_ppeClX_ppeClY_mtchEn))->Fill(xCl, yCl, fEnJEl);
                    } else {
                        *fJEl = TParticle();
                    }
                }
                tOut->Fill();
            }
            delete lvJElPair;
        }
        delete listClsPref;
    } // end of for over events in AliRunLoader

    fOut->cd();
    // output list with histograms
    TList* lTH1F = new TList(); // TH1F
    TList* lTH2F = new TList(); // TH2F
    TList* lTPrf = new TList(); // TProfile
    TList* lTP2D = new TList(); // TProfile2D
    // add all histograms to the lists
    for(Int_t i = kFirstTH1F+1; i < kFirstTH2F; i++) if((TH1F*)arrHistos->At(i)) lTH1F->Add((TH1F*)arrHistos->At(i));
    for(Int_t i = kFirstTH2F+1; i < kFirstTPrf; i++) if((TH2F*)arrHistos->At(i)) lTH2F->Add((TH2F*)arrHistos->At(i));
    for(Int_t i = kFirstTPrf+1; i < kFirstTP2D; i++) if((TProfile*)arrHistos->At(i)) lTPrf->Add((TProfile*)arrHistos->At(i));
    for(Int_t i = kFirstTP2D+1; i < kAll; i++)     if((TProfile2D*)arrHistos->At(i)) lTP2D->Add((TProfile2D*)arrHistos->At(i));
    lTH1F->Write("lTH1F", TObject::kSingleKey);
    lTH2F->Write("lTH2F", TObject::kSingleKey);
    lTPrf->Write("lTPrf", TObject::kSingleKey);
    lTP2D->Write("lTP2D", TObject::kSingleKey);
    // print file content and close it
    fOut->ls();
    fOut->Write("",TObject::kWriteDelete);
    delete fOut;

    // terminate:
    runLoader->Delete();
    fCls->Close();
    delete fCls;
    return;
}