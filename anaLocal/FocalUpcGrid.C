// FocalUpcGrid.C
// David Grund, Nov 05, 2022

#include "FocalUpcGrid_Config.h"
#include "FocalUpcGrid.h"
#include "FocalUpcGrid_CreateHistograms.h"

void FocalUpcGrid(Bool_t isLocal, Bool_t isBoxSim, TString inDir = "", TString outDir = "")
{
    TString outSubDir = "";
    // mass filtering only for J/psi simulations
    if(cutM > 0 && isBoxSim) {
        cout << " ERROR: Cannot do mass cleang for box simulations of electrons/photons. Terminating... " << endl;
        return;
    }
    if(isLocal) 
    {
        outSubDir = "output";
        if(doSupercls) outSubDir += Form("_supCl");
        if(!isBoxSim && matchDirectly) outSubDir += Form("_dirMtch");
        if(cutM > 0)   outSubDir += Form("_cutM%.1f",cutM);
        if(cutE > 0)   outSubDir += Form("_cutE%.1f",cutE);
        outSubDir += "/";
        gSystem->Exec("mkdir -p " + outDir + outSubDir);
    }

    // check if focalClusters.root were produced properly
    TString sClFile = Form("%sfocalClusters.root",outDir.Data());
    if(gSystem->AccessPathName(sClFile.Data()))
    {
        cout << " ERROR: cluster file not found! Terminating." << endl;
        return;
    } 
    // open the file with clusters
    TFile* fCls = new TFile(sClFile.Data());
    cout << " MESSAGE: Loading clusters from: " << fCls->GetName() << endl;

    // define ALICE run loader: open galice.root
    AliRunLoader* runLoader = NULL;
    if(!isLocal) runLoader = AliRunLoader::Open(inDir + "root_archive.zip#galice.root");
    else         runLoader = AliRunLoader::Open(inDir + "galice.root");
    if(!runLoader) 
    {
        cout << " ERROR: AliRunLoader not good! Terminating." << endl;
        return;   
    }
    if(!runLoader->GetAliRun()) runLoader->LoadgAlice();
    if(!runLoader->TreeE()) runLoader->LoadHeader();
    if(!runLoader->TreeK()) runLoader->LoadKinematics();

    // arrays with output histograms
    TObjArray* arrTH1F = NULL;
    TObjArray* arrTH2F = NULL;
    TObjArray* arrTPrf = NULL;
    Int_t nTH1F(-1), nTH2F(-1), nTPrf(-1);
    // if box simulations
    if(isBoxSim) 
    {
        arrTH1F = new TObjArray(kB1_all); Bx_DefineTH1F(arrTH1F); nTH1F = kB1_all;
        arrTH2F = new TObjArray(kB2_all); Bx_DefineTH2F(arrTH2F); nTH2F = kB2_all;
        arrTPrf = new TObjArray(kBP_all); Bx_DefineTPrf(arrTPrf); nTPrf = kBP_all;
    } 
    // if starlight J/psi simulations
    else
    {
        arrTH1F = new TObjArray(kJ1_all); Jp_DefineTH1F(arrTH1F); nTH1F = kJ1_all;
        arrTH2F = new TObjArray(kJ2_all); Jp_DefineTH2F(arrTH2F); nTH2F = kJ2_all;
        arrTPrf = new TObjArray(kJP_all); Jp_DefineTPrf(arrTPrf); nTPrf = kJP_all;
    }

    // output file
    TString sFile = Form("%s%sanalysisResults.root",outDir.Data(),outSubDir.Data());
    TFile* fOut = new TFile(sFile.Data(),"RECREATE");
    // output tree
    TTree* tOut = new TTree("tOut", "output tree containing prefiltered clusters");
    // prefiltered clusters
    Int_t fEvNumber;
    Float_t fEnCl, fPtCl, fEtaCl, fPhiCl;
    tOut->Branch("fEvNumber", &fEvNumber, "fEvNumber/I"); 
    tOut->Branch("fEnCl", &fEnCl, "fEnCl/F");
    tOut->Branch("fPtCl", &fPtCl, "fPtCl/F");
    tOut->Branch("fEtaCl", &fEtaCl, "fEtaCl/F");
    tOut->Branch("fPhiCl", &fPhiCl, "fPhiCl/F");
    // J/psi electrons with which the clusters were matched
    Float_t fEnJEl, fPtJEl, fEtaJEl, fPhiJEl;
    tOut->Branch("fEnJEl", &fEnJEl, "fEnJEl/F");
    tOut->Branch("fPtJEl", &fPtJEl, "fPtJEl/F");
    tOut->Branch("fEtaJEl", &fEtaJEl, "fEtaJEl/F");
    tOut->Branch("fPhiJEl", &fPhiJEl, "fPhiJEl/F");
    gROOT->cd();    

    // output log file
    TString sLog = Form("%s%sanalysis.log",outDir.Data(),outSubDir.Data());
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
                // add the new supercluster to listClsPref
                AliFOCALCluster* supClNew = new AliFOCALCluster(xSupCl/ESupCl,ySupCl/ESupCl,zSupCl/ESupCl,ESupCl,-1);
                listClsPref->Add(supClNew);
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
            ((TH1F*)arrTH1F->At(kB1_mcE))->Fill(stack->Particle(0)->Energy());
            ((TH1F*)arrTH1F->At(kB1_mcPt))->Fill(stack->Particle(0)->Pt());
            ((TH2F*)arrTH2F->At(kB2_mcE_nCls))->Fill(stack->Particle(0)->Energy(), nClsPref);
            // total energy vs MC energy
            ((TH2F*)arrTH2F->At(kB2_totE_mcE))->Fill(ETot, stack->Particle(0)->Energy());
            ((TProfile*)arrTPrf->At(kBP_totE_mcE))->Fill(ETot, stack->Particle(0)->Energy());
            // maximum cluster energy vs MC energy
            ((TH2F*)arrTH2F->At(kB2_maxClE_mcE))->Fill(EClMax, stack->Particle(0)->Energy());
            ((TProfile*)arrTPrf->At(kBP_maxClE_mcE))->Fill(EClMax, stack->Particle(0)->Energy());
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
                ((TH2F*)arrTH2F->At(kB2_clE_mcE))->Fill(clust->E(), part->Energy());
                // radial distance between the cluster and the MC track
                Float_t DeltaX = xCl - x;
                Float_t DeltaY = yCl - y;
                Float_t DeltaR = TMath::Sqrt(TMath::Power(DeltaX,2) + TMath::Power(DeltaY,2));
                ((TH2F*)arrTH2F->At(kB2_clMcDX_clMcDY))->Fill(DeltaX, DeltaY);
                ((TH2F*)arrTH2F->At(kB2_mcE_clMcSep))->Fill(part->Energy(), DeltaR);
                // radial distance between the cluster and the cluster with maximum energy
                if(isClWithMaxE && (iCl != iClMaxE)) {
                    Float_t DeltaX_maxE = xCl - xClMaxE;
                    Float_t DeltaY_maxE = yCl - yClMaxE;
                    Float_t DeltaR_maxE = TMath::Sqrt(TMath::Power(DeltaX_maxE,2) + TMath::Power(DeltaY_maxE,2));
                    ((TH2F*)arrTH2F->At(kB2_clMaxClDX_clMaxClDY))->Fill(DeltaX_maxE, DeltaY_maxE);
                    ((TH2F*)arrTH2F->At(kB2_mcE_clMaxClSep))->Fill(part->Energy(), DeltaR_maxE);
                    ((TH2F*)arrTH2F->At(kB2_clX_clMaxClSep))->Fill((xCl + xClMaxE)/2., DeltaR_maxE);
                    ((TH2F*)arrTH2F->At(kB2_clY_clMaxClSep))->Fill((yCl + yClMaxE)/2., DeltaR_maxE);
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
                if(part->GetPdgCode() == 443) {
                    ((TH1F*)arrTH1F->At(kJ1_mcJEn))->Fill(part->Energy());
                    ((TH1F*)arrTH1F->At(kJ1_mcJPt))->Fill(part->Pt());
                    ((TH1F*)arrTH1F->At(kJ1_mcJRap))->Fill(part->Y());
                    ((TH1F*)arrTH1F->At(kJ1_mcJM))->Fill(part->GetCalcMass());
                    ((TH2F*)arrTH2F->At(kJ2_mcJRap_mcJPt))->Fill(part->Y(), part->Pt());
                }
                // if physical primary (pp) electron (i.e., J/psi direct decay product)
                if(isEleOrPos(part) && stack->IsPhysicalPrimary(iTrk)) {
                    TParticle* mother = stack->Particle(part->GetMother(0));
                    if(mother->GetPdgCode() != 443) cout << " * MESSAGE: Unexpected mother of a pp electron." << endl;
                    ((TH2F*)arrTH2F->At(kJ2_mcJElEta_mcJElPt))->Fill(part->Eta(), part->Pt());
                    ppElectrons.AddLast(part);
                }
            }
            // kinematics of the two physical primary electrons
            TParticle* ppEl1 = NULL;
            TParticle* ppEl2 = NULL;
            TLorentzVector* lvJElPair = new TLorentzVector();
            if(ppElectrons.GetEntries() != 2) {
                cout << "  (!) ERROR: Unexpected number of phys. prim. electrons: " << ppElectrons.GetEntries() << ", terminating..." << endl;
                runLoader->Delete();
                fCls->Close();
                return;
            } else {
                // fill the histogram showing the correlation between energies and transverse momenta of the two pp electrons
                ppEl1 = (TParticle*) ppElectrons[0];
                ppEl2 = (TParticle*) ppElectrons[1];
                ((TH2F*)arrTH2F->At(kJ2_mcJEl1En_mcJEl2En))->Fill(ppEl1->Energy(), ppEl2->Energy());
                ((TH2F*)arrTH2F->At(kJ2_mcJEl1Pt_mcJEl2Pt))->Fill(ppEl1->Pt(), ppEl2->Pt());
                // fill the histograms showing kinematics of J/psi reconstructed from the two pp electrons
                TLorentzVector lvPPEl1;
                lvPPEl1.SetPxPyPzE(ppEl1->Px(),ppEl1->Py(),ppEl1->Pz(),ppEl1->Energy());
                TLorentzVector lvPPEl2;
                lvPPEl2.SetPxPyPzE(ppEl2->Px(),ppEl2->Py(),ppEl2->Pz(),ppEl2->Energy());
                TLorentzVector Jpsi = lvPPEl1 + lvPPEl2;
                lvJElPair->SetPxPyPzE(Jpsi.Px(),Jpsi.Py(),Jpsi.Pz(),Jpsi.Energy());
                ((TH1F*)arrTH1F->At(kJ1_mcJElPairEn))->Fill(lvJElPair->Energy());
                ((TH1F*)arrTH1F->At(kJ1_mcJElPairPt))->Fill(lvJElPair->Pt());
                ((TH1F*)arrTH1F->At(kJ1_mcJElPairRap))->Fill(lvJElPair->Rapidity());
                ((TH1F*)arrTH1F->At(kJ1_mcJElPairM))->Fill(lvJElPair->M());
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

                fEvNumber = iEv;
                // cluster kinematics -> tree
                fEnCl = cl.Energy();
                fPtCl = cl.Pt();
                fEtaCl = cl.Eta();
                fPhiCl = cl.Phi();
                // pp electron kinematics -> tree
                fEnJEl = -1e3;
                fPtJEl = -1e3;
                fEtaJEl = -1e3;
                fPhiJEl = -1e3;

                ((TH1F*)arrTH1F->At(kJ1_clZ))->Fill(zCl);
                // fill some histograms with cluster kinematics
                ((TH2F*)arrTH2F->At(kJ2_clEta_clPhi))->Fill(fEtaCl, fPhiCl);
                ((TH2F*)arrTH2F->At(kJ2_clEta_clPt ))->Fill(fEtaCl, fPtCl);

                TParticle* mtchPhysPrimPart = (TParticle*)arrMtchPhysPrimParts[iCl];
                // if matched to a ppp
                if(mtchPhysPrimPart) 
                {
                    ((TH2F*)arrTH2F->At(kJ2_pppClEn_mtchEn))->Fill(fEnCl, mtchPhysPrimPart->Energy());
                    ((TH2F*)arrTH2F->At(kJ2_pppClEta_mtchEta))->Fill(fEtaCl, mtchPhysPrimPart->Eta());
                    ((TH2F*)arrTH2F->At(kJ2_pppClPt_mtchPt))->Fill(fPtCl, mtchPhysPrimPart->Pt());
                    // if electron (ppe)
                    if(isEleOrPos(mtchPhysPrimPart)) {
                        fEnJEl = mtchPhysPrimPart->Energy();
                        fPtJEl = mtchPhysPrimPart->Pt();
                        fEtaJEl = mtchPhysPrimPart->Eta();
                        fPhiJEl = mtchPhysPrimPart->Phi();
                        ((TH2F*)arrTH2F->At(kJ2_ppeClEn_mtchEn))->Fill(fEnCl, fEnJEl);
                        ((TH2F*)arrTH2F->At(kJ2_ppeClEta_mtchEta))->Fill(fEtaCl, fEtaJEl);
                        ((TH2F*)arrTH2F->At(kJ2_ppeClPt_mtchPt))->Fill(fPtCl, fPtJEl);
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
    TList* lOut1 = new TList(); // TH1F
    TList* lOut2 = new TList(); // TH2F
    TList* lOutP = new TList(); // TProfile
    // add all histograms to this list
    for(Int_t i = 0; i < nTH1F; i++) if((TH1F*)arrTH1F->At(i)) lOut1->Add((TH1F*)arrTH1F->At(i));
    for(Int_t i = 0; i < nTH2F; i++) if((TH2F*)arrTH2F->At(i)) lOut2->Add((TH2F*)arrTH2F->At(i));
    for(Int_t i = 0; i < nTPrf; i++) if((TProfile*)arrTPrf->At(i)) lOutP->Add((TProfile*)arrTPrf->At(i));
    lOut1->Write("lTH1F", TObject::kSingleKey);
    lOut2->Write("lTH2F", TObject::kSingleKey);
    lOutP->Write("lTPrf", TObject::kSingleKey);
    // print file content and close it
    fOut->ls();
    fOut->Write("",TObject::kWriteDelete);
    delete fOut;

    runLoader->Delete();
    fCls->Close();

    return;
}