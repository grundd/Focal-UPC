// FocalUpcAnalysis_SimulateClusters.C

#include <iostream>
// root headers
#include "TSystem.h"
#include "TRandom3.h"
// aliroot headers
#include "AliRunLoader.h"
// my headers
#include "FocalUpcAnalysis_Utilities.h"
#include "FocalUpcAnalysis_DefineOutput.h"

void DoSimulation(TString dataset, Int_t opt, Int_t nEv = -1);
// dataset -> name of a folder where the input files are stored
// nEv -> number of events from each input file to be analyzed 
// (if -1, then analyze all available events)

void FocalUpcAnalysis_SimulateClusters(TString dataset, Int_t nEv = -1)
{
    DoSimulation(dataset,0,nEv);
    DoSimulation(dataset,1,nEv);
    DoSimulation(dataset,2,nEv);
    DoSimulation(dataset,3,nEv);
    return;
}

void DoSimulation(TString dataset, Int_t opt, Int_t nEv)
// opt = 0 -> kinematics taken from MC
//     = 1 -> positions random, energy from MC
//     = 2 -> positions from MC, energy random
//     = 3 -> both random
{
    gSystem->Load("libpythia6_4_28.so"); 
    TString outputFolder = "results/simulateClusters/" + dataset;
    gSystem->Exec(("mkdir -p " + outputFolder).Data());

    // define ALICE run loader: open galice.root
    AliRunLoader *runLoader = NULL;
    runLoader = AliRunLoader::Open(("inputData/" + dataset + "galice.root").Data());
    if(!runLoader) 
    {
        cout << " ERROR: AliRunLoader not good! Terminating." << endl;
        return;   
    }
    if(!runLoader->GetAliRun()) runLoader->LoadgAlice();
    if(!runLoader->TreeE()) runLoader->LoadHeader();
    if(!runLoader->TreeK()) runLoader->LoadKinematics();

    // define histograms
    TH1F* h_clPairM = new TH1F(Form("sim%i_h_clPairM",opt),"inv. mass of simulated cluster pairs;#it{m}_{cl}",50,0.,5.);
    TH2F* h_mcE_clE = new TH2F(Form("sim%i_h_mcE_clE",opt),"MC energy vs simulated cluster energy;#it{E}_{MC} [GeV];#it{E}_{cl} [GeV]",60,0.,120.,60,0.,120.);

    Int_t nEvents(nEv);
    if(nEv == -1) nEvents = runLoader->GetNumberOfEvents();
    cout << "Analyzing " << nEvents << " events:" << endl;
    // loop over MC events contained within Kinematics.root
    Double_t progress = 0.; // perc
    for(Int_t iEv = 0; iEv < nEvents; iEv++) 
    {
        // update the progress bar
        if((iEv+1) % (Int_t)(nEvents/10.) == 0){
            progress += 10.;
            cout << "[" << progress << "%] done." << endl;
        }

        // get current MC event
        Int_t isEventOk = runLoader->GetEvent(iEv);
        // the method GetEvent(i) returns zero if event is loaded succesfully, if not we continue
        if (isEventOk != 0)
        {
            cout << "Event " << iEv << " not OK, skipping." << endl;
            continue;
        }
        
        // get the stack of MC tracks for this event
        AliStack *stack = runLoader->Stack();

        // array to store kinematics of two physical primary (pp) electrons
        TObjArray ppElectrons;

        // go over MC tracks and select two primary electrons
        for(Int_t iTrk = 0; iTrk < stack->GetNtrack(); iTrk++) 
        {
            TParticle *part = stack->Particle(iTrk);
            if(!stack->IsPhysicalPrimary(iTrk)) continue;
            if(TMath::Abs(part->GetPdgCode()) == 11) ppElectrons.AddLast(part);
        }

        // check that we have exactly two primary electrons, otherwise return
        if(!(ppElectrons.GetEntries() == 2))
        {
            cout << ppElectrons.GetEntries() << "primary el. tracks found (!= 2), terminating..." << endl;
            return;
        } 

        // simulate two clusters, both to be near the position of MC tracks in the FOCAL assuming straight trajectories
        // z position will be always taken at 710 cm
        // x, y position will be sampled from gaussian distribution with mean equal to anticipated MC position and std dev of 2 cm
        // energy will be sampled from Landau distribution
        TLorentzVector vecCls[2];
        for(Int_t i = 0; i < 2; i++)
        {
            Float_t xyzE[4] = { 0 };
            Float_t xyzEMC[4] = { 0 };
            TParticle* part = (TParticle*) ppElectrons[i];
            Float_t dz = 710. - part->Vz(); // cm
            // fill MC kinematics
            xyzEMC[0] = part->Vx() + part->Px()/part->Pz() * dz; // cm
            xyzEMC[1] = part->Vy() + part->Py()/part->Pz() * dz; // cm
            xyzEMC[2] = 710.; // cm
            xyzEMC[3] = part->Energy();
            // simulate cluster position and energy
            xyzE[0] = xyzEMC[0];
            xyzE[1] = xyzEMC[1];
            xyzE[2] = xyzEMC[2];
            xyzE[3] = xyzEMC[3];
            if(opt == 1 || opt == 3) {
                xyzE[0] = gRandom->Gaus(xyzEMC[0], 2.);
                xyzE[1] = gRandom->Gaus(xyzEMC[1], 2.);
            }
            if(opt == 2 || opt == 3) {
                Float_t mean = 2. +  10./ xyzEMC[3];
                Float_t x = gRandom->Landau(mean,1.); // by default: mu = 0; sigma = 1 
                xyzE[3] = xyzEMC[3] * (15. - x) / 15.;
            }
            // calculate the cluster four vector
            vecCls[i] = ConvertXYZEtoLorVec(xyzE[0],xyzE[1],xyzE[2],xyzE[3]);
            h_mcE_clE->Fill(xyzEMC[3],xyzE[3]);
        }
        // calculate the invariant mass of the two simulated clusters
        TLorentzVector vecClsCombined = vecCls[0] + vecCls[1];
        Double_t mass = vecClsCombined.M();
        h_clPairM->Fill(mass);
    }
    DrawHisto<TH1F>(h_clPairM,outputFolder);
    DrawHistoCOLZ(h_mcE_clE,outputFolder);
    
    runLoader->Delete();

    return;
}