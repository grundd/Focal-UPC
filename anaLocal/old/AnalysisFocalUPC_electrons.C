// AnalysisFocalUPC_electrons.C

#include <iostream>
#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()
// root headers
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TMarker.h"
#include "TClonesArray.h"
#include "TObjArray.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TROOT.h"
// aliroot headers
#include "AliRunLoader.h"
#include "AliFOCALCluster.h"
#include "AliFOCALClusterizerv2.h"
// my headers
#include "Utilities.h"

// *************************
// options to set:
Bool_t makePlots = kTRUE;
const Float_t cut_E = 1.0; // [GeV] cut on cluster energy
const Float_t cut_m = 0.1; // [GeV] cut on mass of cluster pairs
const Float_t cut_dPhi = 0.4; // [-] cut on the difference in phi angles for primary MC particles and prefiltered clusters
const Int_t nBinsPt = 40;
const Float_t PtLow = 0.;
const Float_t PtUpp = 2.;
const Int_t nBinsEta = 60;
const Float_t EtaLow = 3.;
const Float_t EtaUpp = 6.;
const Int_t nBinsE = 60;
const Float_t ELow = 0.;
const Float_t EUpp = 120.;
// *************************

void AnalysisFocalUPC_electrons(Int_t iAnalysis) 
{
    if(iAnalysis == 0) subfolder = "10-12-2022_BoxElectrons_1000ev/";

    gSystem->Load("libpythia6_4_28.so"); 
    gSystem->Exec(("mkdir -p " + subfolder + "histograms/").Data());
    gSystem->Exec(("mkdir -p " + subfolder + "eventInfo/").Data());

    // if clusters not yet produced, run the clusterizer
    if(gSystem->AccessPathName((subfolder + "focalClusters.root").Data()))
    {
        cout << "AnalyzeJpsi_UPC() MESSAGE: focalClusters.root not found! Running the clusterizer now." << endl;
        gROOT->ProcessLine(".x ClusterizeGrid.C");
    }
    // check if focalClusters.root were produced properly
    if(gSystem->AccessPathName((subfolder + "focalClusters.root").Data()))
    {
        cout << "AnalyzeJpsi_UPC() ERROR: focalClusters.root not found! Terminating." << endl;
        return;
    } 
    // open the file with clusters
    TFile *clusterFile = new TFile((subfolder + "focalClusters.root").Data());
    cout << "AnalyzeJpsi_UPC() MESSAGE: Loading clusters from: " << clusterFile->GetName() << endl;
    // create the output lists
    l1h = new TList;
    l2h = new TList;
    // define the histograms to be stored in the output file
    /*
    TH1F* hMCJpsiPt = CreateHisto1D("hMCJpsiPt",
                    "J/#psi #it{p}_{T} generated;#it{p}_{T} [GeV/#it{c}];counts",
                    nBinsPt,PtLow,PtUpp);
    TH2F* hMCJpsiRapPt = CreateHisto2D("hMCJpsiRapPt",
                    "(#it{y}, #it{p}_{T}) of generated J/#psi ;#it{y} [-];#it{p}_{T} [GeV/#it{c}]",
                    nBinsEta,EtaLow,EtaUpp,nBinsPt,PtLow,PtUpp);
    TH2F* hMCJpsiRapPtAcc = CreateHisto2D("hMCJpsiRapPtAcc",
                    "(#it{y}, #it{p}_{T}) of generated J/#psi with electrons within FOCAL acceptance;#it{y} [-];#it{p}_{T} [GeV/#it{c}]",
                    nBinsEta,EtaLow,EtaUpp,nBinsPt,PtLow,PtUpp);
    TH2F* hMCJpsiEleEtaPt = CreateHisto2D("hMCJpsiEleEtaPt",
                    "(#eta, #it{p}_{T}) of electron from generated J/#psi;#eta [-];#it{p}_{T} [GeV/#it{c}]",
                    nBinsEta,EtaLow,EtaUpp,nBinsPt,PtLow,PtUpp);
    */
    // in the following histograms, a cluster means a prefiltered cluster
    TH2F* hClEtaPhi = CreateHisto2D("hClEtaPhi",
                    "(#eta, #phi) of prefiltered cluster;#eta [-];#phi [-]",
                    nBinsEta,EtaLow,EtaUpp,50,0.,3.2);
    TH2F* hClEtaPt = CreateHisto2D("hClEtaPt", 
                    "(#eta, #it{p}_{T}) of prefiltered cluster;#eta [-];#it{p}_{T} [GeV/#it{c}]",
                    nBinsEta,EtaLow,EtaUpp,nBinsPt,PtLow,PtUpp);
    /*
    TH2F* hClVsMtchEta = CreateHisto2D("hClVsMtchEta", 
                    "#eta of cluster vs. matched #eta;#eta (cluster) [-];#eta (matched primary MC track) [-]",
                    nBinsEta,EtaLow,EtaUpp,nBinsEta,EtaLow,EtaUpp);
    TH2F* hClVsMtchPt = CreateHisto2D("hClVsMtchPt", 
                    "#it{p}_{T} of cluster vs. matched #it{p}_{T};#it{p}_{T} (cluster) [GeV/#it{c}];#it{p}_{T} (matched primary MC track) [GeV/#it{c}]",
                    nBinsPt,PtLow,PtUpp,nBinsPt,PtLow,PtUpp);
    TH2F* hClVsMtchEn = CreateHisto2D("hClVsMtchEn", 
                    "energy of cluster vs. matched energy;#it{E} (cluster) [GeV];#it{E} (matched primary MC track) [GeV]",
                    nBinsE,ELow,EUpp,nBinsE,ELow,EUpp);
    TH2F* hClVsMtchJpsiEleEta = CreateHisto2D("hClVsMtchJpsiEleEta", 
                    "#eta of cluster matched with a primary J/#psi electron vs. matched #eta;#eta (cluster) [-];#eta (matched primary electron) [-]",
                    nBinsEta,EtaLow,EtaUpp,nBinsEta,EtaLow,EtaUpp);
    TH2F* hClVsMtchJpsiElePt = CreateHisto2D("hClVsMtchJpsiElePt", 
                    "#it{p}_{T} of cluster matched with a primary J/#psi electron vs. matched #it{p}_{T};#it{p}_{T} (cluster) [GeV/#it{c}];#it{p}_{T} (matched primary electron) [GeV/#it{c}]",
                    nBinsPt,PtLow,PtUpp,nBinsPt,PtLow,PtUpp);
    TH2F* hClVsMtchJpsiEleEn = CreateHisto2D("hClVsMtchJpsiEleEn", 
                    "energy of cluster matched with a primary J/#psi electron vs. matched energy;#it{E} (cluster) [GeV];#it{E} (matched primary electron) [GeV]",
                    nBinsE,ELow,EUpp,nBinsE,ELow,EUpp);
    TH1F* hClPairMass = CreateHisto1D("hClPairMass", 
                    "inv. mass of cluster pairs;#it{m} [GeV/#it{c}^{2}];counts", 
                    125,0.,5.);
    TH1F* hClPairEta = CreateHisto1D("hClPairEta", 
                    "#eta of cluster pairs;#eta [-];counts",
                    nBinsEta,EtaLow,EtaUpp);
    TH1F* hClPairPt = CreateHisto1D("hClPairPt", 
                    "#it{p}_{T} of cluster pairs;#it{p}_{T} [GeV/#it{c}];counts",
                    nBinsPt,PtLow,PtUpp);
    TH1F* hClPairJpsiEleMass = CreateHisto1D("hClPairJpsiEleMass", 
                    "inv. mass of cluster pairs matched with a pair of primary J/#psi electrons;#it{m} [GeV/#it{c}^{2}];counts", 
                    125,0.,5.);
    TH1F* hClPairJpsiEleEta = CreateHisto1D("hClPairJpsiEleEta", 
                    "#eta of cluster pairs matched with a pair of primary J/#psi electrons;#eta [-];counts",
                    nBinsEta,EtaLow,EtaUpp);
    TH1F* hClPairJpsiElePt = CreateHisto1D("hClPairJpsiElePt", 
                    "#it{p}_{T} of cluster pairs matched with a pair of primary J/#psi electrons;#it{p}_{T} [GeV/#it{c}];counts",
                    nBinsPt,PtLow,PtUpp);
    TH2F* hClPairJpsiEleVsMCJpsiEta = CreateHisto2D("hClPairJpsiEleVsMCJpsiEta", 
                    "#eta of a cluster pair matched with a pair of primary J/#psi electrons vs. #eta of generated J/#psi;#eta (cluster pair) [-];#eta (matched primary electron pair) [-]",
                    nBinsEta,EtaLow,10.,nBinsEta,EtaLow,10.);
    TH2F* hClPairJpsiEleVsMCJpsiPt = CreateHisto2D("hClPairJpsiEleVsMCJpsiPt", 
                    "#it{p}_{T} of a cluster pair matched with a pair of primary J/#psi electrons vs. #it{p}_{T} of generated J/#psi;#it{p}_{T} (cluster pair) [GeV/#it{c}];#it{p}_{T} (matched primary electron pair) [GeV/#it{c}]",
                    nBinsPt,PtLow,PtUpp,nBinsPt,PtLow,PtUpp);
    TH2F* hRelDiffVsMtchEn = CreateHisto2D("hRelDiffVsMtchEn", 
                    "rel. diff. in energy of a matched primary J/#psi electron and cluster energy vs. energy of primary J/#psi electron;(#it{E}_{match} - #it{E}_{cluster}) / #it{E}_{match} [-];counts",
                    nBinsE,ELow,EUpp,40,0.,1.);
    TH2F* hMismatchXY = CreateHisto2D("hMismatchXY", 
                    "XY distance between a cluster and a matched primary J/#psi electron vs. matched energy;#it{E} (matched primary electron) [GeV];#sqrt{(#Delta#it{x})^{2} + (#Delta#it{y})^{2}} [cm]",
                    nBinsE,ELow,EUpp,40,0.,8.);
    */
    TH1F* hMCJpsiEnAll = CreateHisto1D("hMCJpsiEnAll", 
                    "energy of primary J/#psi electron tracks: all;#it{E} [GeV];counts",
                    nBinsE,ELow,EUpp);
    TH1F* hMCJpsiEnMtch = CreateHisto1D("hMCJpsiEnMtch", 
                    "energy of primary J/#psi electron tracks: only matched;#it{E} [GeV];counts",
                    nBinsE,ELow,EUpp);
    TH1F* hMCJpsiEnRatio = CreateHisto1D("hMCJpsiEnRatio", 
                    "energy of primary J/#psi electron tracks: ratio matched/all;#it{E} [GeV];ratio matched/all [-]",
                    nBinsE,ELow,EUpp);

    // load the clusterizer
    AliFOCALClusterizerv2 *clusterizer = new AliFOCALClusterizerv2();
    clusterizer->InitGeometry("geometry.txt");
    AliFOCALGeometry *geometry = clusterizer->GetGeometry();

    // define ALICE run loader: open galice.root
    AliRunLoader *runLoader = NULL;
    runLoader = AliRunLoader::Open((simFolder + subfolder + "galice.root").Data());
    if(!runLoader) 
    {
        cout << "AnalyzeJpsi_UPC() ERROR: AliRunLoader not good! Terminating." << endl;
        return;   
    }
    if(!runLoader->GetAliRun()) runLoader->LoadgAlice();
    if(!runLoader->TreeE()) runLoader->LoadHeader();
    if(!runLoader->TreeK()) runLoader->LoadKinematics();
    // create an array of objects where we will store the primary tracks in each event
    TObjArray primTracks;
    // loop over MC events contained within Kinematics.root
    //for(Int_t iEv = 0; iEv < runLoader->GetNumberOfEvents(); iEv++) 
    for(Int_t iEv = 0; iEv < 10; iEv++) 
    {
        // output text file
        ofstream of(Form("%seventInfo/Ev%i.txt",subfolder.Data(),iEv));
        of << std::fixed << std::setprecision(0);

        // get current MC event
        Int_t isEventOk = runLoader->GetEvent(iEv);
        cout << "Event " << iEv << ":" << endl;
        // the method GetEvent(i) returns zero if event is loaded succesfully, if not we continue
        if (isEventOk != 0)
        {
            cout << "Event " << iEv << " not OK, skipping." << endl;
            continue;
        }
        // get the stack of MC tracks for this event
        AliStack *stack = runLoader->Stack();
        // clean the TObjArray from tracks from the previous event
        primTracks.Clear();
        // vectors that will be filled with eta and energy of tracks being electrons and immediate J/psi daughters
        std::vector<Float_t> eleEta;
        std::vector<Float_t> eleEne;

        // loop over the tracks
        for (Int_t iTrk = 0; iTrk < stack->GetNtrack(); iTrk++) 
        {
            TParticle *part = stack->Particle(iTrk);
            // if not a physical primary, continue
            if(!stack->IsPhysicalPrimary(iTrk)) continue;
            // if primary electron, add it to the array
            if(part->GetPdgCode() == 11)
            {
                primTracks.AddLast(part);
            }
            // here: fill the histograms            
        }
        // get the number of primary tracks: should be only a single electron
        if(primTracks.GetEntries() != 1) 
        {
            cout << "AnalyzeJpsi_UPC() ERROR: More/less than one primary electron: " << primTracks.GetEntries() << endl;
        }

        // canvas showing a 2d view of an event (z-r space)
        TCanvas* cMain = NULL;
        if(makePlots)
        {
            cMain = EventDisplay_PrepareCanvas();
            PlotFOCAL(geometry,cMain);
            // plot the tracks of particles anticipating straight tracks in the direction of the initial momentum vector
            // plot primary and their direct daughters       
            EventDisplay_PlotMCTracks(stack,cMain);
        }

        // get a tree with clusters for this event 
        // (separate tree for each event in subfolders "Event0, Event1, ...")
        TTree *tClusters = NULL;
        if(clusterFile->GetDirectory(Form("Event%i",iEv))) 
        {
            clusterFile->GetDirectory(Form("Event%i",iEv))->GetObject("fTreeR",tClusters);
        }
        else 
        {
            cout << " * AnalyzeJpsi_UPC() MESSAGE: Cannot find event " << iEv << " in a cluster file " << clusterFile->GetName() << ". Skipping this event." << endl;
            clusterFile->ls();
            continue;
        }
        // get the branch with segment-by-segment clusters
        TBranch *bClustersSeg = NULL;
        bClustersSeg = tClusters->GetBranch("AliFOCALClusterItr");
        TClonesArray *clustersSegArray = NULL;
        bClustersSeg->SetAddress(&clustersSegArray);
        bClustersSeg->GetEvent(0);
        Int_t nClustersSeg = clustersSegArray->GetEntries();
        of << "segment clusters: " << nClustersSeg << endl;
        // plot the positions of segment-by-segment clusters
        if(makePlots)
        {
            for(Int_t iClust = 0; iClust < nClustersSeg; iClust++)
            {
                AliFOCALCluster *clust = (AliFOCALCluster*) clustersSegArray->At(iClust);
                Float_t z(clust->Z()), x(clust->X()), y(clust->Y());
                // z-x plane
                TMarker *mX = new TMarker(z,x,0);
                mX->SetMarkerStyle(kOpenSquare);
                mX->SetMarkerColor(kRed);
                mX->SetMarkerSize(.5);
                cMain->cd(1);
                mX->Draw("SAME");
                // z-y plane
                TMarker *mY = new TMarker(z,y,0);
                mY->SetMarkerStyle(kOpenSquare);
                mY->SetMarkerColor(kRed);
                mY->SetMarkerSize(.5);
                cMain->cd(2);
                mY->Draw("SAME");
            }
        }

        // get the branch with final clusters
        TBranch *bClusters = NULL;
        bClusters = tClusters->GetBranch("AliFOCALCluster");
        TClonesArray *clustersArray = NULL;
        bClusters->SetAddress(&clustersArray);
        bClusters->GetEvent(0);
        Int_t nClusters = clustersArray->GetEntries();
        of << "summed clusters: " << nClusters << endl;

        // create a list of prefiltered clusters (with energy higher than the cut)
        TList prefilteredClusters;
        prefilteredClusters.SetOwner(kFALSE);
        for(Int_t iClust = 0; iClust < nClusters; iClust++) 
        {
            AliFOCALCluster *clust = (AliFOCALCluster*) clustersArray->At(iClust);
            if(clust->E() >= cut_E)
            prefilteredClusters.Add(clust);
        }
        Int_t nClustersPref = prefilteredClusters.GetEntries();
        of << "prefil. clusters (with E > " << cut_E << " GeV): " << nClustersPref << endl;

        // plot the positions of prefiltered clusters
        if(makePlots)
        {
            TLegend* lCl = new TLegend(0.35,0.4-(prefilteredClusters.GetEntries()+1)*0.04,0.78,0.4);
            lCl->SetMargin(0.10);
            lCl->SetTextSize(0.025);
            lCl->SetBorderSize(0);
            lCl->AddEntry((TObject*)0,Form("%i summed clusters with E > %.0f GeV and m > %.2f GeV when paired:", prefilteredClusters.GetEntries(), cut_E, cut_m),"");
            for(Int_t iClust = 0; iClust < prefilteredClusters.GetEntries(); iClust++) 
            {
                AliFOCALCluster *clust = (AliFOCALCluster*)prefilteredClusters.At(iClust);
                Float_t z(clust->Z()), x(clust->X()), y(clust->Y());
                Float_t r = TMath::Sqrt(TMath::Power(x,2) + TMath::Power(y,2));
                // z-x plane
                TMarker *mX = new TMarker(z,x,0);
                mX->SetMarkerStyle(kOpenCircle);
                mX->SetMarkerColor(kGreen+1);
                mX->SetMarkerSize(1);
                cMain->cd(1);
                mX->Draw("SAME");
                // z-y plane
                TMarker *mY = new TMarker(z,y,0);
                mY->SetMarkerStyle(kOpenCircle);
                mY->SetMarkerColor(kGreen+1);
                mY->SetMarkerSize(1);
                cMain->cd(2);
                mY->Draw("SAME");
                lCl->AddEntry(mX,Form("x = %.2f cm, y = %.2f cm, E = %.1f GeV", x, y, clust->E()),"P");
            }
            cMain->cd(1);
            lCl->Draw();
        }

        /*
        // go over prefiltered clusters and match them with MC particles
        for(Int_t iCluster = 0; iCluster < prefilteredClusters.GetEntries(); iCluster++) 
        {
            AliFOCALCluster *clust = (AliFOCALCluster*)prefilteredClusters.At(iCluster);
            // get energy and coordinates of this cluster
            Float_t e = clust->E();
            Float_t x = clust->X();
            Float_t y = clust->Y();
            Float_t z = clust->Z();
            // create its 4-vector
            TLorentzVector cl = ConvertXYZEtoLorVec(x,y,z,e);
            
            hClEtaPhi->Fill(cl.Eta(),cl.Phi());
            hClEtaPt->Fill(cl.Eta(),cl.Pt());
        
            // now find a matching MC particle
            // matchedXY [cm] ... cut on the difference in the position in xy plane between a particle with 
            // a straight trajectory in direction of its initial momentum and a cluster 
            // its value will improve based on best match so far
            Float_t matchedXY(1e3), matchedPt(-1), matchedE(-1), matchedEta(-1);
            // go over tracks of primary MC particles
            for(Int_t iPrim = 0; iPrim < 1; iPrim++) 
            {
                TParticle *part = (TParticle*) primTracks[iPrim];
                if((TMath::Abs(part->GetPdgCode()) == 11) && ((stack->Particle(part->GetMother(0)))->GetPdgCode() == 443))
                {
                    hMCJpsiEnAll->Fill(part->Energy());
                } 
                // first do eta-phi matching for isolation cut
                Float_t dEta = TMath::Abs(part->Eta() - cl.Eta());
                Float_t dPhi = TMath::Abs(part->Phi() - cl.Phi());
                if(dPhi > TMath::Pi()) dPhi = 2*TMath::Pi() - dPhi;
                // make the cut on the difference in azimuthal angle
                if(TMath::Abs(dPhi) > cut_dPhi) continue;

                // Vx, Vy and Vz: coordinates of the production vertex of TParticle
                // z-distance between the particle production vertex and the cluster
                Float_t dz = z - part->Vz();
                // x- and y-coordinates where the particle should have created the cluster assuming a straight
                // trajectory in the direction of its initial momentum
                Float_t xPart = part->Vx() + part->Px()/part->Pz() * dz;
                Float_t yPart = part->Vy() + part->Py()/part->Pz() * dz;
                // distance between the cluster coordinates and xPart, yPart at dz
                Float_t dist = TMath::Sqrt(TMath::Power(xPart-x,2) + TMath::Power(yPart-y,2));
                // do the matching
                if(dist < matchedXY) 
                {
                    matchedXY = dist;
                    matchedPt = part->Pt();
                    matchedE = part->Energy();
                    matchedEta = part->Eta();
                }
            }
            // here: fill the histograms
        }  // end of for over iclust
        */
        // save the event info canvas and text file
        if(makePlots)
        {
             cout << " * ";
            cMain->SaveAs(Form("%seventInfo/Ev%i.pdf",subfolder.Data(),iEv));
        } 
        delete cMain;
        of.close();
    }  // end of for over events in AliRunLoader

    //hMCJpsiEnRatio->Sumw2();
    //hMCJpsiEnRatio->Divide(hMCJpsiEnAll);
    // print the histograms
    for(Int_t i = 0; i < l1h->GetSize(); i++)
    {
        TH1F *h = (TH1F*)l1h->At(i);
        PrintHisto1D(h);
    }
    for(Int_t i = 0; i < l2h->GetSize(); i++)
    {
        TH2F *h = (TH2F*)l2h->At(i);
        PrintHisto2D(h);
    }

    // create the output file
    TFile *f = new TFile((subfolder + "analysisElectrons.root").Data(),"RECREATE");
    f->cd();
    l1h->Write("histograms1d", TObject::kSingleKey);
    l2h->Write("histograms2d", TObject::kSingleKey);
    //l1h->ls();
    //l2h->ls();
    //f->ls();
    f->Close();
    clusterFile->Close();
    runLoader->Delete();

    delete clusterizer;
}