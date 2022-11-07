// FocalUpcGrid.h
// David Grund, Oct 15, 2022

// cpp headers
#include <iostream>
#include <fstream>
// root headers
#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TParticle.h"
// aliroot headers
#include "AliStack.h"
#include "AliRunLoader.h"
// focal headers
#include "AliFOCALCluster.h"

using std::cout;
using std::endl;

Bool_t areSame(Float_t a, Float_t b)
{
    // desired precision: ~0.00001
    Float_t epsilon = 1e-5;
    return TMath::Abs(a - b) < epsilon;
}

Bool_t isEleOrPos(TParticle* part)
{
    // is particle an electron or a positron
    return (TMath::Abs(part->GetPdgCode()) == 11);
}

TLorentzVector ConvertXYZEtoLorVec(Float_t x, Float_t y, Float_t z, Float_t e, Bool_t debug = kFALSE)
{
    TLorentzVector vec4;
    // calculate kinematic variables
    Float_t r = TMath::Sqrt(x*x + y*y + z*z);
    Float_t r_xy = TMath::Sqrt(x*x + y*y);
    // polar angle: [0,pi]
    Float_t theta = TMath::ACos(z/r);
    // assume that mass is negligible compared to energy
    Float_t p = e;
    // z-projection of the momentum and transverse momentum
    Float_t pz = p * TMath::Cos(theta);
    Float_t pt = p * TMath::Sin(theta);
    // pseudorapidity
    Float_t eta = 1.0e+6;
    if(theta > 1e-6) eta = -TMath::Log(TMath::Tan(theta/2.));
    // azimuthal angle: [0,2pi)
    Float_t phi = TMath::ATan2(y,x); // returns values between -pi < phi <= pi
    // x- and y-projections of the momentum
    Float_t px = pt * TMath::Cos(phi);
    Float_t py = pt * TMath::Sin(phi);
    // create the 4-vector
    vec4.SetPxPyPzE(px,py,pz,e);
    // check that we calculated everything properly
    if(debug)
    {
        TVector3 vec3_1(x,y,z);
        vec3_1.SetMag(1.);
        TVector3 vec3_2 = vec4.Vect();
        vec3_2.SetMag(1.);
        for(Int_t i = 0; i < 3; i++) if(!areSame(vec3_1[i],vec3_2[i])) Printf("Problem with a coordinate %i: %.10f, %.10f",i+1,vec3_1[i],vec3_2[i]);

        Printf("-------------------------------------------------------");
        Printf("(x,y,z) [cm]: (%.3f,%.3f,%.3f)",vec4.X(),vec4.Y(),vec4.Z());
        Printf("(px,py,pz,E) [GeV]: (%.3f,%.3f,%.3f,%.3f)",vec4.Px(),vec4.Py(),vec4.Pz(),vec4.E());
        Printf("theta [-]: %.4f",vec4.Theta());
        Printf("phi [-]: %.4f",vec4.Phi());
        Printf("eta [-]: %.4f",vec4.Eta());
        Printf("");
    }
    return vec4;
}

// ******************************************************************************************************************
// Function to match clusters with physical primary MC particles
// ******************************************************************************************************************

void TrackCoordinatesAtZ(TParticle* p, Float_t z, Float_t &x, Float_t &y)
{
    Float_t x0(p->Vx()), y0(p->Vy()), z0(p->Vz());
    x = x0 + p->Px() / p->Pz() * (z - z0);
    y = y0 + p->Py() / p->Pz() * (z - z0);
    return;
}

void MatchClsToPhysPrimP(AliStack* stack, TList* listClsPref, vector<Int_t>& idxMtchPhysPrimP, TObjArray* arrMtchPhysPrimP)
{
    arrMtchPhysPrimP->Clear();
    // initialize the arrays
    for(Int_t iCl = 0; iCl < listClsPref->GetEntries(); iCl++)
    {
        idxMtchPhysPrimP.push_back(-1);
        arrMtchPhysPrimP->AddLast(NULL);
    }
    // now do the actual matching
    for(Int_t iCl = 0; iCl < listClsPref->GetEntries(); iCl++)
    {
        AliFOCALCluster *clust = (AliFOCALCluster*) listClsPref->At(iCl);
        if(!clust) continue;
        // get energy and coordinates of this cluster
        Float_t xCl = clust->X();
        Float_t yCl = clust->Y();
        Float_t zCl = clust->Z();
        Float_t ECl = clust->E();
        TLorentzVector cl = ConvertXYZEtoLorVec(xCl,yCl,zCl,ECl);

        // do the matching
        Float_t mtchXY(1e3);
        Float_t mtchXY_dir(1e3);
        TParticle* mtchP_any = NULL;
        TParticle* mtchP_prim = NULL;
        TParticle* mtchP_primDir = NULL;
        Int_t iMtchP_any = -1;
        Int_t iMtchP_prim = -1;
        Int_t iMtchP_primDir = -1;
        for(Int_t iTrk = 0; iTrk < stack->GetNtrack(); iTrk++)
        {
            TParticle *part = stack->Particle(iTrk);
            // if J/psi, continue (is also a physical primary and could cause problems during matching)
            if(part->GetPdgCode() == 443) continue;
            // during matching, skip photons with very low energy (< 0.2 MeV)
            // trajectories of these often overlap with those of pp electron (when projected as straight lines)
            // and clusters could very likely be matched with them instead of pp electrons
            if(part->GetPdgCode() == 22 && part->Energy() < 0.2) continue;
            // eta cut
            Float_t dEta = TMath::Abs(part->Eta() - cl.Eta());
            if(dEta > cutdEta) continue;
            // phi cut
            Float_t dPhi = TMath::Abs(part->Phi() - cl.Phi());
            if(dPhi > TMath::Pi()) dPhi = 2*TMath::Pi() - dPhi;
            if(dPhi > cutdPhi) continue;
            // do matching
            Float_t x(0.), y(0.);
            TrackCoordinatesAtZ(part,zCl,x,y);
            Float_t distXY = TMath::Sqrt(TMath::Power(x-xCl,2) + TMath::Power(y-yCl,2));
            if(distXY < mtchXY) {
                // a new closest MC particle found
                mtchP_any = part;
                iMtchP_any = iTrk;
                mtchXY = distXY;
            }
            if((distXY < mtchXY_dir) && stack->IsPhysicalPrimary(iTrk)) {
                // a new closest physical primary MC particle found
                mtchP_primDir = part;
                iMtchP_primDir = iTrk;
                mtchXY_dir = distXY;
            }
        }
        // if a matching particle found
        if(mtchP_any) {
            // find the physical primary particle from which the matched particle originates
            mtchP_prim = mtchP_any;
            Bool_t physPrimFound = stack->IsPhysicalPrimary(iMtchP_any);
            while(!physPrimFound) {
                Int_t iNewMother = mtchP_prim->GetMother(0);
                mtchP_prim = stack->Particle(iNewMother);
                physPrimFound = stack->IsPhysicalPrimary(iNewMother);
                iMtchP_any = iNewMother;
            }
            iMtchP_prim = iMtchP_any;
        }
        // if going by mothers
        if(!matchDirectly) {
            idxMtchPhysPrimP[iCl] = iMtchP_prim;
            arrMtchPhysPrimP->AddAt(mtchP_prim, iCl);
        // if matching directly to the closest ppp
        } else {
            idxMtchPhysPrimP[iCl] = iMtchP_primDir;
            arrMtchPhysPrimP->AddAt(mtchP_primDir, iCl);
        }
    }    
    return;
}