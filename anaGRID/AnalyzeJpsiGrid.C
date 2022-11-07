#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TProfile.h"
#include "TH2.h"
#include "TH3.h"
#include "TParticle.h"
#include "TSystem.h"
#include "Riostream.h"

#include "AliRunLoader.h"
#include "AliStack.h"
#include "AliFOCALCluster.h"

using std::cout;
using std::endl;

Float_t eta(TParticle* ); // declaration, implementation at end of file
Float_t phi(TParticle* );

void AnalyzeJpsiGrid() {
        
  gSystem->Load("libpythia6_4_28.so");
  // Run the clusterizer ------------------------------------------------------------------      
  gROOT->ProcessLine(".x ClusterizeGrid.C");

  const Float_t Ecut = 50.0; //[GeV] cut on cluster energy for pi0 candidate pairs
  const Float_t kPi0Mass = 0.135;
  
  //NOTE: try asymmetry cuts on candidate cluster pairs

  if(gSystem->AccessPathName("focalClusters.root")) {
    cout << "AnalyzeJpsiGrid() ERROR: focalClusters.root not found!" << endl;      
    return;        
  }      
  TFile *clusterFile = new TFile("focalClusters.root");
  cout << "Clusters from: " << clusterFile->GetName() << endl;
    
  TFile *fout = new TFile("analysisJpsi.root","RECREATE");
  fout->cd();
    
  TH2F* histEtaPtJpsiMC = new TH2F("etaPtJpsi", "J/#psi (#eta,p_{T}) generated", 30, 3.0, 6.0, 20, 0.0, 10.0);
  TH2F* histEtaPtJpsiMC_acc = new TH2F("etaPtJpsi_acc", "J/#psi (#eta,p_{T}) generated x acceptance", 30, 3.0, 6.0, 20, 0.0, 10.0);
  TH2F* histEtaPtJpsiEleMC = new TH2F("etaPtJpsiEle", "Electron from J/#psi (#eta,p_{T}) generated", 30, 3.0, 6.0, 20, 0.0, 10.0);
  TH2F* histEtaPhiCluster = new TH2F("etaPhiCluster", "Cluster (#eta,#phi)", 30, 3.0, 6.0, 50, 0.0, 3.2);
  TH2F* histEtaPtCluster = new TH2F("etaPtCluster", "Cluster (#eta,p_{T})", 30, 3.0, 6.0, 50, 0.0, 10.0);
  TH2F* histEtaRecoVsMatch = new TH2F("etaRecoVsMatch", "Cluster #eta reco vs matched", 30, 3.0, 6.0, 30, 3.0, 6.0);
  TH2F* histPtRecoVsMatch = new TH2F("ptRecoVsMatch", "Cluster p_{T} reco vs matched", 50, 0.0, 10.0, 50, 0.0, 10.0);
  TH2F* histERecoVsMatch = new TH2F("eRecoVsMatch", "Cluster energy reco vs matched", 100, 0.0, 500.0, 100, 0.0, 500.0);
  TH2F* histEtaRecoVsMatch_jpsiEle = new TH2F("etaRecoVsMatch_jpsiEle", "Cluster #eta reco vs matched (jpsi ele)", 30, 3.0, 6.0, 30, 3.0, 6.0);
  TH2F* histPtRecoVsMatch_jpsiEle = new TH2F("ptRecoVsMatch_jpsiEle", "Cluster p_{T} reco vs matched (jpsi ele)", 50, 0.0, 10.0, 50, 0.0, 10.0);
  TH2F* histERecoVsMatch_jpsiEle = new TH2F("eRecoVsMatch_jpsiEle", "Cluster energy reco vs matched (jpsi ele)", 100, 0.0, 500.0, 100, 0.0, 500.0);
  TH1F* histMass = new TH1F("mass", "inv mass", 125, 0.0, 5.0);
  TH1F* histEtaClsPair = new TH1F("etaClsPair", "#eta for cluster pair", 30, 3.0, 6.0);
  TH1F* histPtClsPair = new TH1F("ptClsPair", "p_{T} for cluster pair", 50, 0.0, 10.0);
  TH1F* histMass_jpsiEle = new TH1F("mass_jpsiEle", "inv mass (jpsi electrons)", 125, 0.0, 5.0);
  TH1F* histEtaClsPair_jpsiEle = new TH1F("etaClsPair_jpsiEle", "#eta for cluster pair (jpsi electrons)", 30, 3.0, 6.0);
  TH1F* histPtClsPair_jpsiEle = new TH1F("ptClsPair_jpsiEle", "p_{T} for cluster pair (jpsi electrons)", 50, 0.0, 10.0);
    
  AliRunLoader *runLoader = 0;
    
  //Alice run loader
  runLoader = AliRunLoader::Open("root_archive.zip#galice.root");
        
  if (!runLoader) {
    cout << "AnalyzeJpsiGrid() ERROR: runLoader not good!" << endl;
    return;   
  }
  if (!runLoader->GetAliRun()) runLoader->LoadgAlice();
  if (!runLoader->TreeE()) runLoader->LoadHeader();
  if (!runLoader->TreeK()) runLoader->LoadKinematics();

  TObjArray primTracks;
  for (Int_t ievt = 0; ievt < runLoader->GetNumberOfEvents(); ievt++) {

    // Get MC Stack
    Int_t ie = runLoader->GetEvent(ievt);
    cout << "Event: " << ievt << endl;

    if (ie != 0)
       continue;

    AliStack *stack = runLoader->Stack();

    // Build TObjArray of tracks for matching
    primTracks.Clear();
    std::vector<float> eleEta;
    std::vector<float> elePt;
    std::vector<float> eleE;
      
    Int_t jpsiMotherLabel = -1;
    for (Int_t iTrk = 0; iTrk < stack->GetNtrack(); iTrk++) {
      TParticle *part = stack->Particle(iTrk);
        
      if(part->GetPdgCode()==443) {
        histEtaPtJpsiMC->Fill(part->Eta(), part->Pt());
        jpsiMotherLabel = iTrk;           // NOTE: this only works if there is just one jpsi per event
      }
      if(TMath::Abs(part->GetPdgCode())==11) {
        TParticle* mother = stack->Particle(part->GetMother(0));
        if (mother->GetPdgCode()==443) {
          histEtaPtJpsiEleMC->Fill(part->Eta(), part->Pt());
          eleEta.push_back(part->Eta());
          elePt.push_back(part->Pt());
          eleE.push_back(part->Energy());
        }
      }
        
      if (!stack->IsPhysicalPrimary(iTrk))
        continue;
      if (part->GetFirstDaughter() >= 0 && stack->IsPhysicalPrimary(part->GetFirstDaughter()))  // if decay daughters are phys prim, only use those
        continue;

      primTracks.AddLast(part);
      // could also store extrapolated positions here...
    }
    if (eleEta.size()>=2) {
      int switches = 1;
      while (switches==0) {
        switches = 0;      
        for (int i=0; i<eleEta.size()-1; i++) {
          if (eleE[i]<eleE[i+1]) {
            float e = eleE[i];
            eleE[i] = eleE[i+1];
            eleE[i+1] = e;
            float eta = eleEta[i];
            eleEta[i] = eleEta[i+1];
            eleEta[i+1] = eta;
            switches++;
          }        
        }
      }
      // apply kinematic cut to the jpsi electrons (use the most energetic 2 of them)
      if (eleEta[0]>3.5 && eleEta[0]<5.8 && eleEta[1]>3.5 && eleEta[1]<5.8) {
        TParticle* jpsiMother = stack->Particle(jpsiMotherLabel);      
        histEtaPtJpsiMC_acc->Fill(jpsiMother->Eta(), jpsiMother->Pt());
      }
    }
      
    Int_t nPrim = primTracks.GetEntries();
    // Get Clusters
    TTree *tClusters = 0;
    if (clusterFile->GetDirectory(Form("Event%i",ievt))) {
      clusterFile->GetDirectory(Form("Event%i",ievt))->GetObject("fTreeR",tClusters);
    }
    else {
      cout << "Cannot find event " << ievt << " in cluster file " << clusterFile->GetName() << endl;
      clusterFile->ls();
    }

    TBranch * bClusters;
    //if (segmentToAnalyze == -1)
    bClusters = tClusters->GetBranch("AliFOCALCluster");  // Branch for final clusters
    // bClusters = tClusters->GetBranch("AliFOCALClusterItr"); // Branch for segment-by-segment clusters
    TClonesArray * clustersArray = 0;
    bClusters->SetAddress(&clustersArray);
    bClusters->GetEvent(0);
      
    TList prefilteredClusters;       // list of prefiltered clusters
    prefilteredClusters.SetOwner(kFALSE);
    for (Int_t iClust = 0; iClust < clustersArray->GetEntries(); iClust++) {
      AliFOCALCluster *clust = (AliFOCALCluster*) clustersArray->At(iClust);
      if (clust->E() >= Ecut)
        prefilteredClusters.Add(clust);
    }
    TIter nextCluster(&prefilteredClusters);
    TIter nextCluster2(&prefilteredClusters);
    Int_t nclusters = clustersArray->GetEntries();
    cout << "nclusters:: " << nclusters << endl;
    for(Int_t iClust = 0; iClust < nclusters; iClust++) {
      AliFOCALCluster* clust = (AliFOCALCluster*)nextCluster();
      if(!clust)
        continue;
      TVector3 p1(clust->X(),clust->Y(),clust->Z());
      p1.SetMag(clust->E());
      nextCluster2.Reset();
      Bool_t foundLowMass = kFALSE;
      for(Int_t iClust2 = 0; iClust2 < nclusters; iClust2++) {
        AliFOCALCluster* clust2 = (AliFOCALCluster*)nextCluster2();
        if(!clust2)
          continue;
        TVector3 p2(clust2->X(),clust2->Y(),clust2->Z());
        p2.SetMag(clust2->E());
        Double_t mass = 2*clust->E()*clust2->E()-2*p1*p2;
        if(mass<0.0)
          mass = 0.0;
        else
          mass = TMath::Sqrt(mass);
        if (mass<0.000) {
          prefilteredClusters.Remove(clust);
          prefilteredClusters.Remove(clust2);
          foundLowMass = kTRUE;
          break;
        }
      }
      if(foundLowMass)
        continue;
    }
      
    for (Int_t iClust = 0; iClust < prefilteredClusters.GetEntries(); iClust++) {
      AliFOCALCluster *clust = (AliFOCALCluster*) prefilteredClusters.At(iClust);
      Float_t e = clust->E();
      Float_t x = clust->X();
      Float_t y = clust->Y();
      Float_t z = clust->Z();
      Float_t pt = clust->E();

      Float_t trans = TMath::Sqrt(x*x + y*y);      
      Float_t theta = TMath::Pi()/2;
      if (z != 0) {
        pt = trans/z*e;
        theta = TMath::ATan(trans/z);
      }
      Float_t etaCls = 1.0e+6;
      if (theta > 1e-6)
        etaCls = -TMath::Log(TMath::Tan(theta/2.));
      Float_t phiCls = TMath::ATan2(y,-x);
      histEtaPhiCluster->Fill(etaCls,phiCls);
      histEtaPtCluster->Fill(etaCls,pt);
        
      // Now find matching particle
      Float_t matchDist = 1000;
      Float_t ptMatch;
      Float_t eMatch;
      Float_t etaMatch;
      Bool_t isJpsiEle = kFALSE;

      for (Int_t iPrim = 0; iPrim < nPrim; iPrim++) {
          
        TParticle *part = (TParticle*) primTracks[iPrim];
        // first do eta-phi matching for isolation cut
        Float_t deta = eta(part)-etaCls;
        Float_t dphi = phi(part)-phiCls;
        if (dphi < -TMath::Pi())
          dphi += 2*TMath::Pi();
        if (dphi > TMath::Pi())
          dphi -= 2*TMath::Pi();
        if (TMath::Abs(dphi) > 0.4)
          continue;

        // Only works if Pz > 0
        Float_t dz = z - part->Vz();
        Float_t xpart =  part->Vx() + part->Px()/part->Pz() * dz;
        Float_t ypart =  part->Vy() + part->Py()/part->Pz() * dz;
        Float_t dist = TMath::Sqrt((xpart-x)*(xpart-x)+(ypart-y)*(ypart-y));
        //cout << "x " << x << ", y " << y << ", dist " << dist << endl;
        if (dist < matchDist) {
          matchDist = dist;
          ptMatch = TMath::Sqrt(part->Px()*part->Px()+part->Py()*part->Py());
          eMatch = part->Energy();
          etaMatch = eta(part);
          if(TMath::Abs(part->GetPdgCode()) == 11) {
            TParticle* mother = stack->Particle(part->GetMother(0));
            if (mother->GetPdgCode()==443) {
              isJpsiEle = kTRUE;
            }
          }
          else
            isJpsiEle = kFALSE;
        }
      }
      histERecoVsMatch->Fill(e,eMatch);
      histPtRecoVsMatch->Fill(pt,ptMatch);
      histEtaRecoVsMatch->Fill(etaCls,etaMatch);
      if(isJpsiEle) {
        histERecoVsMatch_jpsiEle->Fill(e,eMatch);
        histPtRecoVsMatch_jpsiEle->Fill(pt,ptMatch);
        histEtaRecoVsMatch_jpsiEle->Fill(etaCls,etaMatch);        
      }

      if (clust->E() >= Ecut) { 
        TVector3 p1(x,y,z);
        p1.SetMag(e);
        for (Int_t iClust2 = iClust+1; iClust2 < prefilteredClusters.GetEntries(); iClust2++) {
          AliFOCALCluster *clust2 = (AliFOCALCluster*) prefilteredClusters.At(iClust2);
          if (clust2->E() < Ecut) {
            continue;
          }
          Float_t e2 = clust2->E();
          Float_t x2 = clust2->X();
          Float_t y2 = clust2->Y();
          Float_t z2 = clust2->Z();
          Float_t pt2 = clust2->E();

          Float_t trans2 = TMath::Sqrt(x2*x2 + y2*y2);      
          Float_t theta2 = TMath::Pi()/2;
          if (z2 != 0) {
            pt2 = trans2/z2*e2;
            theta2 = TMath::ATan(trans2/z2);
          }
          Float_t etaCls2 = 1.0e+6;
          if (theta2 > 1e-6)
            etaCls2 = -TMath::Log(TMath::Tan(theta2/2.));
          Float_t phiCls2 = TMath::ATan2(y2,-x2);
            
          TVector3 p2(clust2->X(),clust2->Y(),clust2->Z());
          p2.SetMag(clust2->E());
          Double_t mass = 2*clust->E()*clust2->E()-2*p1*p2;
          if (mass < 0) {
            mass = 0;
          }
          else {
            mass = TMath::Sqrt(mass);
          }
          TVector3 ptot(p1+p2);
          histMass->Fill(mass);
          histEtaClsPair->Fill(ptot.Eta());
          histPtClsPair->Fill(ptot.Perp());
         
          // Now find matching particle for the second cluster
          Float_t matchDist2 = 1000;
          Float_t ptMatch2;
          Float_t eMatch2;
          Float_t etaMatch2;
          Bool_t isJpsiEle2 = kFALSE;

          for (Int_t iPrim = 0; iPrim < nPrim; iPrim++) {
            TParticle *part = (TParticle*) primTracks[iPrim];
            // first do eta-phi matching for isolation cut
            Float_t deta = eta(part)-etaCls2;
            Float_t dphi = phi(part)-phiCls2;
            if (dphi < -TMath::Pi())
              dphi += 2*TMath::Pi();
            if (dphi > TMath::Pi())
              dphi -= 2*TMath::Pi();
            if (TMath::Abs(dphi) > 0.4)
              continue;

            // Only works if Pz > 0
            Float_t dz = z2 - part->Vz();
            Float_t xpart =  part->Vx() + part->Px()/part->Pz() * dz;
            Float_t ypart =  part->Vy() + part->Py()/part->Pz() * dz;
            Float_t dist = TMath::Sqrt((xpart-x2)*(xpart-x2)+(ypart-y2)*(ypart-y2));
            //cout << "x " << x << ", y " << y << ", dist " << dist << endl;
            if (dist < matchDist2) {
              matchDist2 = dist;
              ptMatch2 = TMath::Sqrt(part->Px()*part->Px()+part->Py()*part->Py());
              eMatch2 = part->Energy();
              etaMatch2 = eta(part);
              if(TMath::Abs(part->GetPdgCode()) == 11) {
                TParticle* mother = stack->Particle(part->GetMother(0));
                if (mother->GetPdgCode()==443) {
                  isJpsiEle2 = kTRUE;
                }
              }
              else
                isJpsiEle2 = kFALSE;
            }
          }
          if (isJpsiEle && isJpsiEle2) {
            histMass_jpsiEle->Fill(mass);
            histEtaClsPair_jpsiEle->Fill(ptot.Eta());
            histPtClsPair_jpsiEle->Fill(ptot.Perp());      
          }
            
        }  // end for over iClust2
      }
    }  // end for over iClust
  }  // end for over events

  fout->Write();
  fout->Close();
  clusterFile->Close();
  runLoader->Delete();
}

//____________________________________________________________________________
Float_t eta(TParticle *part) {
  Double_t pt = sqrt (part->Px()*part->Px()+part->Py()*part->Py());
  if (pt == 0)
    return 9999;
  return -log(tan(TMath::ATan2(pt,part->Pz())/2));
}

//____________________________________________________________________________
Float_t phi(TParticle *part) {
  return TMath::ATan2(part->Py(),-part->Px());
}


