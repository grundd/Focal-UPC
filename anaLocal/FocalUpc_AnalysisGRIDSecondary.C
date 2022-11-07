// FocalUpc_AnalysisGRIDSecondary.C
// David Grund, Nov 06, 2022

void FocalUpc_AnalysisGRIDSecondary()
{
    // combine clusters into pairs
    for(Int_t iCl1 = 0; iCl1 < nClsPref; iCl1++) 
    {
        AliFOCALCluster* clust1 = (AliFOCALCluster*) listClsPref->At(iCl1);
        if(!clust1) continue;
        // get energy and coordinates of this cluster
        Float_t xCl1 = clust1->X();
        Float_t yCl1 = clust1->Y();
        Float_t zCl1 = clust1->Z();
        Float_t ECl1 = clust1->E();
        TLorentzVector cl1 = ConvertXYZEtoLorVec(xCl1,yCl1,zCl1,ECl1);
        // go over all possible pairs of clusters
        for(Int_t iCl2 = iCl1+1; iCl2 < nClsPref; iCl2++) 
        {
            AliFOCALCluster* clust2 = (AliFOCALCluster*) listClsPref->At(iCl2);
            if(!clust2) continue;
            // get energy and coordinates of this cluster
            Float_t xCl2 = clust2->X();
            Float_t yCl2 = clust2->Y();
            Float_t zCl2 = clust2->Z();
            Float_t ECl2 = clust2->E();
            TLorentzVector cl2 = ConvertXYZEtoLorVec(xCl2,yCl2,zCl2,ECl2);

            // add the momentum vectors of the two clusters to get the total momentum
            TLorentzVector cl12 = cl1 + cl2;
            // cluster pair kinematics -> tree
            fEnClPair = cl12.Energy();
            fPtClPair = cl12.Pt();
            fEtaClPair = cl12.Eta();
            fPhiClPair = cl12.Phi();
            Float_t mass = cl12.M();
            // pp electron pair kinematics -> tree
            fEnJElPair = -1e3;
            fPtJElPair = -1e3;
            fEtaJElPair = -1e3;
            fPhiJElPair = -1e3;
            // fill some kinematic histograms
            ((TH1F*)arrTH1F->At(kJ1_clPairEn))->Fill(fEnClPair);
            ((TH1F*)arrTH1F->At(kJ1_clPairPt))->Fill(fPtClPair);
            if(mass > 2.8) ((TH1F*)arrTH1F->At(kJ1_clPairPt_massCut))->Fill(fPtClPair);
            ((TH1F*)arrTH1F->At(kJ1_clPairRap))->Fill(cl12.Rapidity());
            ((TH1F*)arrTH1F->At(kJ1_clPairM))->Fill(mass);
            // radial separation between the pairs of clusters
            Float_t sepCl = TMath::Sqrt(TMath::Power(xCl1-xCl2,2) + TMath::Power(yCl1-yCl2,2));
            ((TH1F*)arrTH1F->At(kJ1_clPairSep))->Fill(sepCl);
            // radial separation between cluster pair vs between the pair of ppe
            Float_t x1(0.), y1(0.);
            TrackCoordinatesAtZ(ppEl1,zCl1,x1,y1);
            Float_t x2(0.), y2(0.);
            TrackCoordinatesAtZ(ppEl2,zCl2,x2,y2);
            Float_t sepMC = TMath::Sqrt(TMath::Power(x1-x2,2) + TMath::Power(y1-y2,2));
            ((TH2F*)arrTH2F->At(kJ2_clPairSep_mcJElSep))->Fill(sepCl,sepMC);
            // if both paired to a different pp electron:
            if((idxMtchPhysPrimParts[iCl1] != idxMtchPhysPrimParts[iCl2]) && idxMtchPhysPrimParts[iCl1] != -1 && idxMtchPhysPrimParts[iCl2] != -1)
            {
                fEnJElPair = lvJElPair->Energy();
                fPtJElPair = lvJElPair->Pt();
                fEtaJElPair = lvJElPair->Eta();
                fPhiJElPair = lvJElPair->Phi();
                ((TH1F*)arrTH1F->At(kJ1_ppeClPairM))->Fill(mass);
                ((TH1F*)arrTH1F->At(kJ1_ppeClPairSep))->Fill(sepCl); 
                ((TH2F*)arrTH2F->At(kJ2_ppeClPairEn_mtchEn))->Fill(fEnClPair,fEnJElPair);
                ((TH2F*)arrTH2F->At(kJ2_ppeClPairRap_mtchRap))->Fill(cl12.Rapidity(),lvJElPair->Rapidity());
                ((TH2F*)arrTH2F->At(kJ2_ppeClPairPt_mtchPt))->Fill(fPtClPair,fPtJElPair);
                ((TH2F*)arrTH2F->At(kJ2_ppeClPairM_mtchM))->Fill(mass,lvJElPair->M());
            } else {
                ((TH1F*)arrTH1F->At(kJ1_sameppeClPairSep))->Fill(sepCl); 
            }
            tOutClPairs->Fill();
        } // end of for over iCl2
    } // end of for over iCl1
    return;
}