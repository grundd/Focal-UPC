// author: Ionut Arsene, i.c.arsene@gsi.de
// date:   11 august 2010

#include <iostream>
using namespace std;

#include <TF1.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TH1.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TLorentzVector.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TParticle.h>
#include <TMath.h>
#include <TVector3.h>
#include <TGraph.h>

#include "AliDecayerPolarized.h"
//#include "AliDecayerPythia.h"

enum gkHistos {
  kMotherPtThrown = 0,
  kMotherPtAcc,
  kMotherYThrown,
  kMotherYAcc,
  kMotherPtYThrown,
  kMotherPtYAcc,
  kPosLegPtThrown,
  kPosLegPtAcc,
  kPosLegPThrown,
  kPosLegPAcc,
  kPosLegPYThrown,
  kPosLegPYAcc,
  kPosLegYThrown,
  kPosLegYAcc,
  kPosLegPtYThrown,
  kPosLegPtYAcc,
  kPosLegPhiThrown,
  kPosLegPhiAcc,
  kNegLegPtThrown,
  kNegLegPtAcc,
  kNegLegPThrown,
  kNegLegPAcc,
  kNegLegPYThrown,
  kNegLegPYAcc,
  kNegLegYThrown,
  kNegLegYAcc,
  kNegLegPtYThrown,
  kNegLegPtYAcc,
  kNegLegPhiThrown,
  kNegLegPhiAcc,
  kLegCorrelationPxThrown,
  kLegCorrelationPxAcc,
  kLegCorrelationPyThrown,
  kLegCorrelationPyAcc,
  kLegCorrelationPzThrown,
  kLegCorrelationPzAcc,
  kLegCorrelationPtThrown,
  kLegCorrelationPtAcc,
  kOpAngleThrown,
  kOpAngleAcc,
  kOpAnglevsptThrown,
  kOpAnglevsptAcc,
  kOpAnglevsmotherptThrown,
  kOpAnglevsmotherptAcc,
  kCosThetaCSThrown,
  kCosThetaCSAcc,
  kCosThetaHEThrown,
  kCosThetaHEAcc,
  kNhistos = kCosThetaHEAcc,
  kEffPt,
  kEffPosLegP,
  kEffNegLegP,
  kEffPosLegPY,
  kEffNegLegPY,
  kEffY,
  kEffPtY,
  kEffCosThetaCS,
  kEffCosThetaHE,
  kEffPtProf,
  kEffYProf,
  kEffPtYProf,
  kEffCosThetaCSProf,
  kEffCosThetaHEProf,
  kEffOpAngle,
  kEffIntegratedProf,
  kNall
};

// TRD supermodules_____________________________________________________________________
Int_t gSmStatus[18] = {1, 1, 1, 1, 1, 1, 
		       1, 1, 1, 1, 1, 1, 
		       1, 1, 1, 1, 1, 1};
Double_t gPhiLimits[19] = {0.000000, 0.349066, 0.698132, 1.047197, 1.396263, 1.745329, 
			   2.094395, 2.443461, 2.792527, 3.141593, 3.490659, 3.839724,
			   4.188790, 4.537856, 4.886922, 5.235988, 5.585054, 5.934119,
			   6.283185};
// _____________________________________________________________________________________

TF1 *gSpectrumPt;
TF1 *gSpectrumY=new TF1("jpsiSpectrumY","gaus",-8.,8.);
TH1F* gReferencePt;       // reference Pt from the LHC10e2 simulation
TH1D* gEffCosThetaStar;   // full J/Psi efficiency vs cos(theta*) CS/HE
TH1F* gPIDEffVsP;            // electron from jpsi decay PID eff as a function of p 

const Double_t gkMassJpsi    = 3.096916;
const Double_t gkMassPsiPrime = 3.686097;
const Double_t gkMassUpsilon = 9.4603;
const Double_t kElectronMass = 0.0005109989;
const Double_t kMuonMass     = 0.105658;
const Double_t kProtonMass   = 0.93827231;
const Double_t kBeamEnergy   = 3500.;

//_______________________________________________________________________________________
// prototypes
Bool_t InTRDAcceptance(TParticle* part);
void ShootParticle(TLorentzVector* mother, Int_t pdgCode=443, 
		   Double_t lowRap=-0.8, Double_t highRap=0.8, 
		   Double_t lowPt = 0., Double_t highPt = 50.);
Double_t CosThetaCS(TParticle* negPart, TParticle* posPart, Int_t daughterPdg);
Double_t CosThetaHE(TParticle* negPart, TParticle* posPart, Int_t daughterPdg);
void DefineHistograms(TObjArray* array, Int_t parMotherCode, Double_t parHighPt);
void FillHistograms(TObjArray* histos, TClonesArray* particles, Bool_t accepted, Int_t parDaughterCode);
Bool_t CheckTRDtriggerHQU(TParticle* particle);
Bool_t CheckTRDtriggerHSE(TParticle* particle);
void BuildLHCb13TeVPtSpectrum();
void BuildPP5TeVInputSpectrum(const Char_t* filename);
void BuildPP5TeVSingleLegPIDeff(const Char_t* filename, const Char_t* histName);
void BuildPP5TeVItsTpcMatchingEffSyst();
Bool_t CheckTPCpidEff(TParticle* particle);
Double_t DecayPolarizedQuarkonia(const Char_t* parOutFile = "jpsiKinematics.root", 
				 Double_t parPrecision = 0.01,      // relative precision required
				 Int_t    parMotherCode = 443,      // only J/Psi(443) and Upsilon(553) to electrons implemented
				 Int_t    parDaughterCode = 11,     // 11-electrons; 13-muons
				 Double_t parAlpha=1.0,             // alpha = (L-2T)/(L+2T)
				 Int_t    parPolarizType = 1,       // 1- Collins-Soper; 2-Helicity
				 Double_t parLowY=3.20,             // low y for generated J/Psi's       
				 Double_t parHighY=6.00,            // high y  --"--
				 Double_t parLowPt=0.0,             // low pt
				 Double_t parHighPt=0.4,            // high pt
				 Double_t parLegCutLowEta=3.40,     // low eta cut for the legs
				 Double_t parLegCutHighEta=5.80,    // high eta cut for the legs
				 Double_t parLeg1CutLowPt=0.0,      // low pt cut on legs
				 Double_t parLeg1CutHighPt=10.0,    // high pt cut on legs
				 Double_t parLeg2CutLowPt=0.0,      // low pt cut on legs
				 Double_t parLeg2CutHighPt=10.0,    // high pt cut on legs
				 Bool_t   parRequireTRD = kFALSE,   // switch on/off the TRD acceptance requirement
				 Bool_t   parFillHistos = kTRUE, 
                                 const Char_t* effFile="", const Char_t* effName="");
     

//_________________________________________________________________________________________
void CheckStatistics(Int_t parEntries,            // how many J/Psi's
		     Float_t parAlpha,            // alpha (L-2T)/(L+2T)
		     Int_t parType = 1) {         // 1:Collins-Soper, 2:Helicity
  // Small routine to check the response of the detector for a given input alpha parameter.
  //  Input (well known alpha) --> Output (statistics and efficiency affected alpha)
  //    What is done here:
  //  a - Generate a number of cos(theta star) values according to the parAlpha
  //  b - Use the detector eficiency (vs cos theta star) to decide wheter to accept or not the J/Psi
  //  c - The result is a uncorrected cos theta star distr.
  //  d - The result is multiplied by the efficiency histogram to get a "real" detector response
  //      which is then fitted by the 1+alpha*cos2(theta)
  
  // initialize the random generator ++++++++++++++++++++++++++++++++++++++++++++++
  TDatime time;
  gRandom->SetSeed(time.GetTime());

  // dN/d(cos theta) function +++++++++++++++++++++++++++++++++++++++++++++++++++++
  TF1* dNdCosTheta = new TF1("dNdCosTheta", "[0]*(1+[1]*x*x)", -1., 1.);
  dNdCosTheta->SetParameter(0, 1.0);
  dNdCosTheta->SetParameter(1, parAlpha);
  dNdCosTheta->SetParName(0, "norm");
  dNdCosTheta->SetParName(1, "alpha");

  // get the efficiency histograms ++++++++++++++++++++++++++++++++++++++++++++++++
  TFile* effFile=TFile::Open("Efficiencies_iter2_bcd.root");
  gEffCosThetaStar = (parType == 1 ? (TH1D*)effFile->Get("fullCorrectionSPDfirst_ThetaCS") : 
		      (TH1D*)effFile->Get("fullCorrectionSPDfirst_ThetaHE"));
  // histo with the "real" cos theta distribution +++++++++++++++++++++++++++++++++
  TH1F *histo = new TH1F("histo", "cos(#theta^{*}_{CS/HE}), statistics study",20,-1,1);	
  histo->Sumw2();

  Int_t entries = 0;
  while(entries<parEntries) {
    Double_t cosThetaStar = dNdCosTheta->GetRandom();
    Double_t probability = gRandom->Rndm();
    Double_t eff = gEffCosThetaStar->GetBinContent(gEffCosThetaStar->FindBin(cosThetaStar));

    if(probability<eff) {
      histo->Fill(cosThetaStar, 1.0/eff);
      entries++;
    }
  }

  histo->Draw();
  TF1* fit = new TF1("fit", "[0]*(1+[1]*x*x)", -1., 1.);
  histo->Fit(fit, "ME", "", -1., 1.);
}

//_______________________________________________________________________________________
void MakeAlphaDependence(Double_t parMinAlpha = -1.0,       // minimum alpha
			 Double_t parMaxAlpha = +1.0,       // maximum alpha
			 Double_t parDeltaAlpha = 0.1,      // step for alpha
			 const Char_t* parOutFile = "jpsiKineEff.root",    // output file
			 Double_t parPrecision = 0.01,      // relative precision required
			 Int_t    parMotherCode = 443,      // only J/Psi(443) and Upsilon(553) to electrons implemented
			 Int_t    parDaughterCode = 11,     // 11-electrons; 13-muons
			 Int_t    parPolarizType = 1,       // 1- Collins-Soper; 2-Helicity
			 Double_t parLowY=-0.90,            // low y for generated J/Psi's       
			 Double_t parHighY=+0.90,           // high y  --"--
			 Double_t parLowPt=0.0,             // low pt
			 Double_t parHighPt=10.0,           // high pt
			 Double_t parLegCutLowEta=-0.90,    // low eta cut for the legs
			 Double_t parLegCutHighEta=+0.90,   // high eta cut for the legs
			 Double_t parLeg1CutLowPt=1.0,       // low pt cut on legs
			 Double_t parLeg1CutHighPt=100.0,    // high pt cut on legs
			 Double_t parLeg2CutLowPt=1.0,       // low pt cut on legs
			 Double_t parLeg2CutHighPt=100.0,    // high pt cut on legs
			 Bool_t   parRequireTRD = kFALSE    // switch on/off the TRD acceptance requirement
      ) {
  //
  // Calculate the alpha dependence of the quarkonia kinematic acceptance. 
  // The interval and step for the alpha parameter are required, as well as the parameters for the 
  // DecayPolarizedQuarkonia() function
  //

  TGraph* effVsAlpha = new TGraph();
  effVsAlpha->SetName("effVsAlpha");
  Int_t iPoint = 0;
  for(Double_t iAlpha=parMinAlpha+0.5*parDeltaAlpha; iAlpha<parMaxAlpha; iAlpha+=parDeltaAlpha) {
    Double_t eff = DecayPolarizedQuarkonia("dummy.root", parPrecision, parMotherCode, parDaughterCode, 
					   iAlpha, parPolarizType, 
					   parLowY, parHighY, parLowPt, parHighPt, parLegCutLowEta, parLegCutHighEta, 
					   parLeg1CutLowPt, parLeg1CutHighPt, parLeg2CutLowPt, parLeg2CutHighPt, 
					   parRequireTRD, kFALSE);
    effVsAlpha->SetPoint(iPoint++, iAlpha, eff);
  }
  //  effVsAlpha->Draw("A*");
  TFile* save=new TFile(parOutFile, "RECREATE");
  effVsAlpha->Write();
  save->Close();
}

//_______________________________________________________________________________________
void MakeLegPtDependence(Double_t parMinLegPt = 1.0,       // minimum one leg pt
			 Double_t parMaxLegPt = 3.0,       // maximum one leg pt
			 Double_t parDeltaLegPt = 0.1,     // step for the one leg pt
			 const Char_t* parOutFile = "jpsiKineEff.root",    // output file
			 Double_t parPrecision = 0.01,      // relative precision required
			 Int_t    parMotherCode = 443,      // only J/Psi(443) and Upsilon(553) to electrons implemented
			 Int_t    parDaughterCode = 11,     // 11-electrons; 13-muons
			 Double_t parPolarizAlpha = 0.0,    // alpha parameter for the polarization
			 Int_t    parPolarizType = 1,       // 1- Collins-Soper; 2-Helicity
			 Double_t parLowY=-0.90,            // low y for generated J/Psi's       
			 Double_t parHighY=+0.90,           // high y  --"--
			 Double_t parLowPt=0.0,             // low pt
			 Double_t parHighPt=10.0,           // high pt
			 Double_t parLegCutLowEta=-0.90,    // low eta cut for the legs
			 Double_t parLegCutHighEta=+0.90,   // high eta cut for the legs
			 Double_t parLegCutLowPt=1.0,       // low pt cut on legs
			 Double_t parLegCutHighPt=100.0,    // high pt cut on legs
			 Bool_t   parRequireTRD = kFALSE    // switch on/off the TRD acceptance requirement
  ) {
    
  //
  // Calculate the one leg minimum pt cut dependence of the quarkonia kinematic acceptance. 
  // The interval and step for the cut are required, as well as the other parameters for the 
  // DecayPolarizedQuarkonia() function
  //

  TGraph* effVsLegPt = new TGraph();
  effVsLegPt->SetName("effVsLegPt");
  Int_t iPoint = 0;
  for(Double_t iLegPt=parMinLegPt+0.5*parDeltaLegPt; iLegPt<parMaxLegPt; iLegPt+=parDeltaLegPt) {
    Double_t eff = DecayPolarizedQuarkonia("dummy.root", parPrecision, parMotherCode, parDaughterCode, 
					   parPolarizAlpha, parPolarizType, 
					   parLowY, parHighY, parLowPt, parHighPt, parLegCutLowEta, parLegCutHighEta, 
					   iLegPt, parLegCutHighPt, parLegCutLowPt, parLegCutHighPt, 
					   parRequireTRD, kFALSE);
    effVsLegPt->SetPoint(iPoint++, iLegPt, eff);
  }
  //  effVsAlpha->Draw("A*");
  TFile* save=new TFile(parOutFile, "RECREATE");
  effVsLegPt->Write();
  save->Close();

}


//_______________________________________________________________________________________
Double_t DecayPolarizedQuarkonia(const Char_t* parOutFile/* = "jpsiKinematics.root"*/, 
				 Double_t parPrecision/* = 0.01*/,      // relative precision required
				 Int_t    parMotherCode/* = 443*/,      // only J/Psi(443) and Upsilon(553) to electrons implemented
				 Int_t    parDaughterCode/* = 11*/,     // 11-electrons; 13-muons
				 Double_t parAlpha/*=0.0*/,             // alpha = (L-2T)/(L+2T)
				 Int_t    parPolarizType/* = 1*/,       // 1- Collins-Soper; 2-Helicity
				 Double_t parLowY/*=-0.90*/,            // low y for generated J/Psi's       
				 Double_t parHighY/*=+0.90*/,           // high y  --"--
				 Double_t parLowPt/*=0.0*/,             // low pt
				 Double_t parHighPt/*=10.0*/,           // high pt
				 Double_t parLegCutLowEta/*=-0.90*/,    // low eta cut for the legs
				 Double_t parLegCutHighEta/*=+0.90*/,   // high eta cut for the legs
				 Double_t parLeg1CutLowPt/*=1.0*/,       // low pt cut on leg
				 Double_t parLeg1CutHighPt/*=100.0*/,    // high pt cut on leg
				 Double_t parLeg2CutLowPt/*=1.0*/,       // low pt cut on leg
				 Double_t parLeg2CutHighPt/*=100.0*/,    // high pt cut on leg
				 Bool_t   parRequireTRD/* = kFALSE*/,   // switch on/off the TRD acceptance requirement
				 Bool_t   parFillHistos/* = kTRUE*/,  // do not fill/save histograms
                                 const Char_t* effFile, const Char_t* effName) {
  //
  // Shoot quarkonia in defined rapidity and pt windows and compute kinematic efficiencies
  //

  // Function to smear the pseudo-rapidity of the quarkonia legs
  TF1* etaLegSmearFunc = 0x0;
  //  etaLegSmearFunc = new TF1("etaLegSmearFunc", "gaus", -3., 3.);
  //  etaLegSmearFunc->SetParameters(1.0, 0.0, 0.1);   // width of 0.1 rapidity units

  // initialize the random number generator +++++++++++++++++++++++++++++++++++++++
  TDatime time;
  gRandom->SetSeed(time.GetTime());

  // set the rapidity probability distribution function +++++++++++++++++++++++++++
  // TODO: Use a realistic y distribution
  gSpectrumY->SetParameter(0, 1.);
  gSpectrumY->SetParameter(1, 0.);
  gSpectrumY->SetParameter(2, 2.);

  // protect against unsuported mother code parameters
  if(parMotherCode!=443 && parMotherCode!=553) {
    cout << "Unsuported mother code!" << endl;
    return -1.;
  }
  
  // initialize the polarized particle decayer ++++++++++++++++++++++++++++++++++++
  if(parPolarizType<1 || parPolarizType>2) {
    cout << "Bad polarization frame!" << endl;
    cout << "  1 - Collins-Soper; 2 - Helicity" << endl;
    return -1.;
  }
  AliDecayerPolarized::Polar_t polarizFlag = AliDecayerPolarized::kColSop;
  if(parPolarizType==2) polarizFlag = AliDecayerPolarized::kHelicity;
  if(TMath::Abs(parDaughterCode)!=11 && TMath::Abs(parDaughterCode)!=13) {
    cout << "Bad polarization frame!" << endl;
    cout << "  1 - Collins-Soper; 2 - Helicity" << endl;
    return -1.;
  }
  AliDecayerPolarized::FinState_t finState = AliDecayerPolarized::kElectron;
  if(parDaughterCode==13) finState = AliDecayerPolarized::kMuon;
  AliDecayerPolarized *decayer = new AliDecayerPolarized(parAlpha, polarizFlag, finState);
  
  
  /*AliDecayerPythia* decayer = new AliDecayerPythia();
  decayer->Init();
  decayer->SetForceDecay(kPsiPrimeJpsiDiElectron);*/
  
  
  // get the reference Pt spectrum (hardcoded for J/Psi) ++++++++++++++++++++++++++
  // now the spectrum is the CDF scaled to 7 TeV
  // TODO: remove the hardcoded file and histogram name
  BuildPP5TeVInputSpectrum(effFile);
  //BuildPP5TeVSingleLegPIDeff(effFile, effName);
  //BuildPP5TeVItsTpcMatchingEffSyst();
  
  // fit the spectrum with a double exponential in mt
  //  Double_t mass = (parMotherCode==443 ? gkMassJpsi : gkMassUpsilon);
  // parameterization inspired from F.Bossu et al., arxiv 1103.2394
  
  // open output file and define histograms
  TFile* outputFile;
  TObjArray* histos;
  if(parFillHistos) {
    outputFile = new TFile(parOutFile, "RECREATE");
    histos = new TObjArray(kNall);
    DefineHistograms(histos, parMotherCode, parHighPt);
  }
  
  // prepare for the loop over generated quarkonia
  TLorentzVector mother;
  TClonesArray *array=new TClonesArray(TParticle::Class());
  array->SetOwner();
  Long_t nPassed = 0;
  Long_t nGenerated = 0;
  Bool_t passed = kFALSE;
  Double_t currPrecision = 1;
  Int_t nPredefinedEvents = -1;
  if(parPrecision>1.0) {
    nPredefinedEvents = parPrecision;
    parPrecision = 0.00001;
  }
  while(currPrecision>parPrecision) {
    if(nGenerated%10000==0)
      cout << "event/precision :  " << nGenerated << " / " << currPrecision << "\r" << flush;

    ShootParticle(&mother, parMotherCode, parLowY, parHighY, parLowPt, parHighPt);
    decayer->Decay(parMotherCode, &mother);
    decayer->ImportParticles(array);
    //cout << "multiplicity :: " << array->GetEntries() << endl;
    TParticle* part1=(TParticle*)array->At(1);
    TParticle* part2=(TParticle*)array->At(2);
    TParticle* partPos=part1;
    TParticle* partNeg=part2;
    if(part1->GetPdgCode()>0 && part2->GetPdgCode()<0) {
      partPos = part2;
      partNeg = part1;
    }

    // smear the psudo-rapidity of the legs to imitate a detector effect
    Double_t partNegEta = partNeg->Eta();
    Double_t partPosEta = partPos->Eta();
    if(etaLegSmearFunc) {
      partNegEta += etaLegSmearFunc->GetRandom();
      partPosEta += etaLegSmearFunc->GetRandom();
    }    
    // here apply the kinematic cuts +++++++++++++++++++++++++++++++++++++++++
    passed = (parLegCutLowEta<partNegEta && partNegEta<parLegCutHighEta &&   // negative leg eta
	      parLegCutLowEta<partPosEta && partPosEta<parLegCutHighEta &&   // positive leg eta
              ((partNeg->Pt()>parLeg1CutLowPt && partNeg->Pt()<parLeg1CutHighPt &&     // pt cut on legs  
	       partPos->Pt()>parLeg2CutLowPt && partPos->Pt()<parLeg2CutHighPt) ||
	      (partNeg->Pt()>parLeg2CutLowPt && partNeg->Pt()<parLeg2CutHighPt &&       
	       partPos->Pt()>parLeg1CutLowPt && partPos->Pt()<parLeg1CutHighPt)) ?      
	      kTRUE : kFALSE);
    if(parRequireTRD) {
       passed = passed && CheckTPCpidEff(partNeg)  && CheckTPCpidEff(partPos);
       //passed = passed && InTRDAcceptance(partNeg) && InTRDAcceptance(partPos);   // TRD
       //passed = passed && (CheckTRDtriggerHQU(partNeg) || CheckTRDtriggerHQU(partPos));    // TRD HQU trigger
       //passed = passed && (CheckTRDtriggerHSE(partNeg) || CheckTRDtriggerHSE(partPos));    // TRD HSE trigger
    }

    //    passed = IsSelected(partPos, partNeg, parRequireTRD);
    
    if(parFillHistos)
      FillHistograms(histos, array, passed, parDaughterCode);

    if(passed) ++nPassed;
    ++nGenerated;
    array->Clear();

    currPrecision = (nPassed>0 && nGenerated>0 ? TMath::Sqrt((1.0/Double_t(nPassed)) + (1.0/Double_t(nGenerated))) : 1);
    if(nPredefinedEvents>0 && nGenerated>=nPredefinedEvents) break; 
  }   // end while

  cout << "Integrated efficiency :  " << nPassed << " / " << nGenerated 
       << " = " << Double_t(nPassed)/Double_t(nGenerated)
       << " +/- " << (nPassed>0 && nGenerated>0 ? (Double_t(nPassed)/Double_t(nGenerated))*TMath::Sqrt((1.0/Double_t(nPassed)) + 
												       (1.0/Double_t(nGenerated))) : 0.0) << endl;
  if(parFillHistos)
    cout << "Mother <pt> : " << ((TH1F*)histos->At(kMotherPtThrown))->GetMean() << endl;
  
  // normalize the histograms to the bins size
  if(parFillHistos) {
    for(Int_t iHist=kMotherPtThrown; iHist<=kNhistos; iHist++) {
      TH1* hist = (TH1*)histos->At(iHist);
      if(hist->GetDimension()==1) {
	for(Int_t iX=1; iX<=((TH1F*)hist)->GetXaxis()->GetNbins(); iX++) {
	  hist->SetBinContent(iX, hist->GetBinContent(iX)/hist->GetBinWidth(iX));
	  hist->SetBinError(iX, hist->GetBinError(iX)/hist->GetBinWidth(iX));
	}
      }
      if(hist->GetDimension()==2) {
	for(Int_t iX=1; iX<=hist->GetXaxis()->GetNbins(); iX++) {
	  for(Int_t iY=1; iY<=hist->GetYaxis()->GetNbins(); iY++) {
	    hist->SetBinContent(iX, iY, hist->GetBinContent(iX,iY)/hist->GetXaxis()->GetBinWidth(iX)/hist->GetYaxis()->GetBinWidth(iY));
	    hist->SetBinError(iX, iY, hist->GetBinError(iX,iY)/hist->GetXaxis()->GetBinWidth(iX)/hist->GetYaxis()->GetBinWidth(iY));
	  }
	}
      }
      //    hist->Scale(Double_t(nGenerated)/hist->Integral());
    }     // end for-loop over histograms 

    ((TH1F*)histos->At(kEffPt))->Divide((TH1F*)histos->At(kMotherPtAcc), (TH1F*)histos->At(kMotherPtThrown));
    ((TH1F*)histos->At(kEffY))->Divide((TH1F*)histos->At(kMotherYAcc), (TH1F*)histos->At(kMotherYThrown));
    ((TH2F*)histos->At(kEffPtY))->Divide((TH2F*)histos->At(kMotherPtYAcc), (TH2F*)histos->At(kMotherPtYThrown));
    ((TH1F*)histos->At(kEffCosThetaCS))->Divide((TH1F*)histos->At(kCosThetaCSAcc), (TH1F*)histos->At(kCosThetaCSThrown));
    ((TH1F*)histos->At(kEffCosThetaHE))->Divide((TH1F*)histos->At(kCosThetaHEAcc), (TH1F*)histos->At(kCosThetaHEThrown));
    outputFile->Write();
    outputFile->Close();
  }    // end if(parFillHistograms)  

  return Double_t(nPassed)/Double_t(nGenerated);
}


//_________________________________________________________________________
Bool_t InTRDAcceptance(TParticle* particle) {
  //
  // Check if a particle falls in the TRD acceptance
  //
  
  Double_t phi = particle->Phi();
  Bool_t in=kFALSE;
  for(Int_t i=0; i<18; i++) {
    if(gSmStatus[i] && phi>gPhiLimits[i] && phi<gPhiLimits[i+1])
      in = kTRUE;
  }
  return in;
}


//_________________________________________________________________________
Bool_t CheckTRDtriggerHQU(TParticle* particle) {
   //
   // Check if a particle fires the TRD HQU trigger
   //
   if(particle->Pt()<2.0) return kFALSE;
   
   const Double_t kPidEff = 0.46;         // online PID > 135
   const Double_t kSagittaEff = 0.6;  // sagitta cut of < 0.2
   const Double_t kTrackingEffPos = 0.5;  // for positrons, 5 or 6 tracklets
   const Double_t kTrackingEffNeg = 0.61;  // for electrons, 5 or 6 tracklets
      
   Double_t trackingEffPosVsPt = 0.5;
   if(particle->Pt()<3.0) trackingEffPosVsPt = 0.2 + 0.1*particle->Pt();
      
   Double_t totalEff = kPidEff * kSagittaEff;
   //if(particle->GetPdgCode()<0) totalEff *= kTrackingEffPos;     // this is the positron
   if(particle->GetPdgCode()<0) totalEff *= trackingEffPosVsPt;     // this is the positron
   if(particle->GetPdgCode()>0) totalEff *= kTrackingEffNeg;     // this is the electron
   
   if(gRandom->Rndm()<totalEff) return kTRUE;
   else return kFALSE;
}


//_________________________________________________________________________
Bool_t CheckTPCpidEff(TParticle* particle) {
   //
   // Check if a particle passes the single leg PID cut
   //
   Double_t p = particle->P();
   //cout << p << endl;
   if(particle->P()<1.0) return kTRUE;
   //if(p>=6.0) p=5.99;    // assume the eff to be saturated above this threshold to avoid fluctuations due to limited stats in the conversion sample
  
   Int_t bin = gPIDEffVsP->GetXaxis()->FindBin(p);
  //cout << "  " << bin << endl;
   Double_t eff = gPIDEffVsP->GetBinContent(bin);
   //cout << "eff1 = " << eff << endl;
   if(eff>1.0) eff = 1.0/eff;
   //cout << "eff2 = " << eff << endl;
   Double_t effErr = gPIDEffVsP->GetBinError(bin);
   //cout << "effErr = " << effErr << endl;
   eff = gRandom->Gaus(eff,effErr);
   //cout << "eff3 = " << eff << endl;
   if(eff>1.0) eff=1.0;
   //cout << "eff4 = " << eff << endl;
   if(gRandom->Rndm()>eff) {
      //cout << "failed" << endl;
      return kFALSE;
   }
   else {
      //cout << "passed" << endl;
      return kTRUE;
   }
}


//_________________________________________________________________________
Bool_t CheckTRDtriggerHSE(TParticle* particle) {
   //
   // Check if a particle fires the TRD HSE trigger
   //
   if(particle->Pt()<3.0) return kFALSE;
   
   const Double_t kPidEff = 0.68;         // online PID > 120
   const Double_t kSagittaEff = 0.6;  // sagitta cut of < 0.2
   const Double_t kTrackingEffPos = 0.5;  // for positrons, 5 or 6 tracklets
   const Double_t kTrackingEffNeg = 0.61;  // for electrons, 5 or 6 tracklets
   
   Double_t totalEff = kPidEff * kSagittaEff;
   if(particle->GetPdgCode()<0) totalEff *= kTrackingEffPos;     // this is the positron
   if(particle->GetPdgCode()>0) totalEff *= kTrackingEffNeg;     // this is the electron
   
   if(gRandom->Rndm()<totalEff) return kTRUE;
   else return kFALSE;
}


//__________________________________________________________________________
void ShootParticle(TLorentzVector* mother, Int_t pdgCode, 
		   Double_t lowRap, Double_t highRap, 
		   Double_t lowPt, Double_t highPt) {
  // rapidity -----------------------
  Double_t y = gSpectrumY->GetRandom(lowRap,highRap);
  // pt -----------------------------
  //  Double_t pt = WeightedRandom(gSpectrumPt,lowPt,highPt);
  Double_t pt=gSpectrumPt->GetRandom(lowPt, highPt);
  Double_t mass = 0;
  if(pdgCode==443) mass = gkMassJpsi;
  if(pdgCode==100443) mass = gkMassPsiPrime;
  if(pdgCode==553) mass = gkMassUpsilon;
  // scale pt to the mother mass (in lack of Upsilon pt spectrum)
  // TODO: reference Upsilon pt spectrum
  Double_t mt = TMath::Sqrt(pt*pt+gkMassJpsi*gkMassJpsi);
  mt -= gkMassJpsi;     // mt - m0
  mt += mass;
  pt = TMath::Sqrt(mt*mt-mass*mass);
  // phi (uniform) -------------------
  Double_t phi = 2.0*TMath::Pi()*gRandom->Rndm();
  mother->SetPx(pt*TMath::Cos(phi));
  mother->SetPy(pt*TMath::Sin(phi));
  mother->SetPz(TMath::Sqrt(pt*pt+mass*mass)*TMath::SinH(y));
  mother->SetE(TMath::Sqrt(pt*pt+mass*mass)*TMath::CosH(y));
  return;
}

//____________________________________________________________________________
Double_t WeightedRandom(TF1* distribution, Double_t min, Double_t max) {
  // von Neumann procedure to generate weighted random numbers according to some distribution

  // normalize the spectrum function
  // NOTE: Make sure this is a positive definite function !!!!!!!!
  if(TMath::Abs(distribution->GetMaximum(min,max))>1e-10)
    distribution->SetParameter(0, 1.0/distribution->GetMaximum(min,max));

  Double_t number = -999.;
  while(1) {
    // generate a random number
    number = min + (max-min)*gRandom->Rndm();
    if(gRandom->Rndm() < distribution->Eval(number)) {
      return number;
    }
  }
  return -999.;   // we should not get here!!
}

//_____________________________________________________________________________
Double_t CosThetaCS(TParticle* negPart, TParticle* posPart, 
		    Int_t daughterPdg) {
  // Calculate cos(theta) in the J/Psi rest frame in the CS picture

  Double_t pxNeg = negPart->Px();  Double_t pyNeg = negPart->Py(); Double_t pzNeg = negPart->Pz();
  Double_t pxPos = posPart->Px();  Double_t pyPos = posPart->Py(); Double_t pzPos = posPart->Pz();
  Double_t mass = kElectronMass;
  if(TMath::Abs(daughterPdg)==13)
    mass = kMuonMass;
  
  // laboratory frame
  TLorentzVector projectile4Mom(0.,0.,-kBeamEnergy,TMath::Sqrt(kBeamEnergy*kBeamEnergy+kProtonMass*kProtonMass)); 
  TLorentzVector target4Mom(0.,0.,kBeamEnergy,TMath::Sqrt(kBeamEnergy*kBeamEnergy+kProtonMass*kProtonMass)); 
  TLorentzVector negPart4Mom(pxNeg,pyNeg,pzNeg,TMath::Sqrt(pxNeg*pxNeg+pyNeg*pyNeg+pzNeg*pzNeg+mass*mass));
  TLorentzVector posPart4Mom(pxPos,pyPos,pzPos,TMath::Sqrt(pxPos*pxPos+pyPos*pyPos+pzPos*pzPos+mass*mass));
  TLorentzVector mother4Mom=negPart4Mom+posPart4Mom;

  // boost back to the J/Psi - Upsilon rest frame 
  TVector3 beta;
  beta = (-1.0/mother4Mom.E())*mother4Mom.Vect();
  posPart4Mom.Boost(beta);
  negPart4Mom.Boost(beta);
  projectile4Mom.Boost(beta);
  target4Mom.Boost(beta);

  // z axis for the CS angle (difference of the 2 beams unit vectors in the J/Psi rest frame) 
  TVector3 zAxisCS=(((projectile4Mom.Vect()).Unit())-((target4Mom.Vect()).Unit())).Unit();

  // calculate the cos(theta*) for the positive particle
  Double_t cosThetaStar = zAxisCS.Dot((posPart4Mom.Vect()).Unit());
  return cosThetaStar;
  //  cout << "cosThetaStar = " << cosThetaStar << endl;
}

//_____________________________________________________________________________
Double_t CosThetaHE(TParticle* negPart, TParticle* posPart,
		    Int_t daughterPdg) {
  // Calculate cos(theta) in the J/Psi rest frame in the HE picture

  Double_t pxNeg = negPart->Px();  Double_t pyNeg = negPart->Py(); Double_t pzNeg = negPart->Pz();
  Double_t pxPos = posPart->Px();  Double_t pyPos = posPart->Py(); Double_t pzPos = posPart->Pz();
  Double_t mass = kElectronMass;
  if(TMath::Abs(daughterPdg)==13)
    mass = kMuonMass;

  // laboratory frame
  TLorentzVector projectile4Mom(0.,0.,-kBeamEnergy,TMath::Sqrt(kBeamEnergy*kBeamEnergy+kProtonMass*kProtonMass)); 
  TLorentzVector target4Mom(0.,0.,kBeamEnergy,TMath::Sqrt(kBeamEnergy*kBeamEnergy+kProtonMass*kProtonMass)); 
  TLorentzVector negPart4Mom(pxNeg,pyNeg,pzNeg,TMath::Sqrt(pxNeg*pxNeg+pyNeg*pyNeg+pzNeg*pzNeg+mass*mass));
  TLorentzVector posPart4Mom(pxPos,pyPos,pzPos,TMath::Sqrt(pxPos*pxPos+pyPos*pyPos+pzPos*pzPos+mass*mass));
  TLorentzVector mother4Mom=negPart4Mom+posPart4Mom;  

  // boost back to the J/Psi rest frame 
  TVector3 beta;
  beta = (-1.0/mother4Mom.E())*mother4Mom.Vect();
  posPart4Mom.Boost(beta);
  negPart4Mom.Boost(beta);
  projectile4Mom.Boost(beta);
  target4Mom.Boost(beta);

  // z axis for the CS angle (difference of the 2 beams unit vectors in the J/Psi rest frame)
  TVector3 zAxisHE = (mother4Mom.Vect()).Unit(); 
  
  // calculate the cos(theta*) for the positive particle
  Double_t cosThetaStar = zAxisHE.Dot((posPart4Mom.Vect()).Unit());
  return cosThetaStar;
  //  cout << "cosThetaStar = " << cosThetaStar << endl;
}

//_______________________________________________________________________________________
Double_t OpeningAngle(TParticle* negPart, TParticle* posPart,
		  Int_t daughterPdg) {

    Double_t pxNeg = negPart->Px();  Double_t pyNeg = negPart->Py(); Double_t pzNeg = negPart->Pz();
    Double_t pxPos = posPart->Px();  Double_t pyPos = posPart->Py(); Double_t pzPos = posPart->Pz();
    Double_t mass = kElectronMass;
    if(TMath::Abs(daughterPdg)==13)
	mass = kMuonMass;

    // laboratory frame
    TLorentzVector negPart4Mom(pxNeg,pyNeg,pzNeg,TMath::Sqrt(pxNeg*pxNeg+pyNeg*pyNeg+pzNeg*pzNeg+mass*mass));
    TLorentzVector posPart4Mom(pxPos,pyPos,pzPos,TMath::Sqrt(pxPos*pxPos+pyPos*pyPos+pzPos*pzPos+mass*mass));

    Double_t opangle = TMath::RadToDeg() * negPart4Mom.Angle(posPart4Mom.Vect());
    return opangle;
}

//_______________________________________________________________________________________
void BuildLHCb13TeVPtSpectrum() {
   //
   // Use LHCb 13 TeV J/psi prompt spectrum published in arXiv 1509.00771
   //
   gReferencePt = new TH1F("LHCb13TeVInclusive", "LHCb inclusive J/psi spectrum in pp collisions at 13 TeV", 14, 0., 14.);
   gReferencePt->SetBinContent(1,   1059+113.2);
   gReferencePt->SetBinContent(2,   2079+276.3);
   gReferencePt->SetBinContent(3,   1863+300.6);
   gReferencePt->SetBinContent(4,   1218+238.6);
   gReferencePt->SetBinContent(5,   755+159.9);
   gReferencePt->SetBinContent(6,   419+103.3);
   gReferencePt->SetBinContent(7,   236.7+67.9);
   gReferencePt->SetBinContent(8,   130.2+42.9);
   gReferencePt->SetBinContent(9,   76.0+25.3);
   gReferencePt->SetBinContent(10, 47.0+18.9);
   gReferencePt->SetBinContent(11, 28.7+14.2);
   gReferencePt->SetBinContent(12, 16.9+9.0);
   gReferencePt->SetBinContent(13, 11.5+6.5);
   gReferencePt->SetBinContent(14, 7.2+4.9);
   
   gReferencePt->SetBinError(1, 35);
   gReferencePt->SetBinError(2, 30);
   gReferencePt->SetBinError(3, 27);
   gReferencePt->SetBinError(4, 20);
   gReferencePt->SetBinError(5, 14);
   gReferencePt->SetBinError(6, 8);
   gReferencePt->SetBinError(7, 8);
   gReferencePt->SetBinError(8, 5);
   gReferencePt->SetBinError(9, 2);
   gReferencePt->SetBinError(10, 1.6);
   gReferencePt->SetBinError(11, 1.2);
   gReferencePt->SetBinError(12, 1.0);
   gReferencePt->SetBinError(13, 0.7);
   gReferencePt->SetBinError(14, 0.6);
}

//_______________________________________________________________________________________
void BuildPP5TeVInputSpectrum(const Char_t* filename) {
   //
   //
   //
   gSpectrumPt = new TF1("fpl","[0]*x/TMath::Power((1+TMath::Power(x/[1],[3])),[2])",0.,30.);
   gSpectrumPt->SetParameters(1.0,3.67987,3.00625,2.0);
   
   gReferencePt = new TH1F("gReferencePt","", 300, 0.0, 30.0);
   gReferencePt->Add(gSpectrumPt);
}


//_______________________________________________________________________________________
void BuildPP5TeVSingleLegPIDeff(const Char_t* filename, const Char_t* histName) {
   //
   //
   //
   TFile* effFile = TFile::Open(filename);
   TH1F* eff = (TH1F*)effFile->Get(histName);
   
   gPIDEffVsP = (TH1F*)eff->Clone("gPIDEffVsP");
}

void BuildPP5TeVItsTpcMatchingEffSyst() {
   //
   //
   //
   Double_t ptLims[10] = {0.0, 0.8, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 15.0, 50.0};
   Double_t syst[9] = {0.99, 0.99, 2.03, 2.60, 2.62, 2.44, 2.66, 2.72, 2.72};

   gPIDEffVsP = new TH1F("gPIDEffVsP", "", 9, ptLims);
   for(Int_t i=1; i<=9; ++i) gPIDEffVsP->SetBinContent(i, 1.0-0.01*syst[i-1]);  
}


//________________________________________________________________________________________
void FillHistograms(TObjArray* histos, TClonesArray* particles, Bool_t accepted,
		    Int_t parDaughterCode) {

  // fill histograms
  TParticle* mother=(TParticle*)particles->At(0);
  TParticle* part1=(TParticle*)particles->At(1);
  TParticle* part2=(TParticle*)particles->At(2);
  TParticle* partPos=part1;
  TParticle* partNeg=part2;
  if(part1->GetPdgCode()>0 && part2->GetPdgCode()<0) {
    partPos = part2;
    partNeg = part1;
  }
  
  Double_t cosThetaCS = CosThetaCS(partNeg, partPos, parDaughterCode);
  Double_t cosThetaHE = CosThetaHE(partNeg, partPos, parDaughterCode);
  Double_t opangle = OpeningAngle(partNeg, partPos, parDaughterCode);
  ((TH1F*)histos->At(kMotherPtThrown))->Fill(mother->Pt());
  ((TH1F*)histos->At(kMotherYThrown))->Fill(mother->Y());
  ((TH2F*)histos->At(kMotherPtYThrown))->Fill(mother->Y(), mother->Pt());
  ((TH1F*)histos->At(kPosLegPtThrown))->Fill(partPos->Pt());
  ((TH1F*)histos->At(kPosLegPThrown))->Fill(partPos->P());
  ((TH2F*)histos->At(kPosLegPYThrown))->Fill(partPos->Y(),partPos->P());
  ((TH1F*)histos->At(kPosLegYThrown))->Fill(partPos->Y());
  ((TH2F*)histos->At(kPosLegPtYThrown))->Fill(partPos->Y(), partPos->Pt());
  ((TH1F*)histos->At(kNegLegPtThrown))->Fill(partNeg->Pt());
  ((TH1F*)histos->At(kNegLegPThrown))->Fill(partNeg->P());
  ((TH2F*)histos->At(kNegLegPYThrown))->Fill(partNeg->Y(),partNeg->P());
  ((TH1F*)histos->At(kNegLegYThrown))->Fill(partNeg->Y());
  ((TH2F*)histos->At(kNegLegPtYThrown))->Fill(partNeg->Y(), partNeg->Pt());
  ((TH1F*)histos->At(kPosLegPhiThrown))->Fill(partPos->Phi());
  ((TH1F*)histos->At(kNegLegPhiThrown))->Fill(partNeg->Phi());
  ((TH2F*)histos->At(kLegCorrelationPxThrown))->Fill(partPos->Px(), partNeg->Px());
  ((TH2F*)histos->At(kLegCorrelationPyThrown))->Fill(partPos->Py(), partNeg->Py());
  ((TH2F*)histos->At(kLegCorrelationPzThrown))->Fill(partPos->Pz(), partNeg->Pz());
  ((TH2F*)histos->At(kLegCorrelationPtThrown))->Fill(partPos->Pt(), partNeg->Pt());
  ((TH1F*)histos->At(kCosThetaCSThrown))->Fill(cosThetaCS);
  ((TH1F*)histos->At(kCosThetaHEThrown))->Fill(cosThetaHE);
  ((TH1F*)histos->At(kOpAngleThrown))->Fill(opangle);
  ((TH2F*)histos->At(kOpAnglevsptThrown))->Fill(opangle, partNeg->Pt());
  ((TH2F*)histos->At(kOpAnglevsmotherptThrown))->Fill(opangle, mother->Pt());
  if(accepted) {
    ((TH1F*)histos->At(kMotherPtThrown+1))->Fill(mother->Pt());
    ((TH1F*)histos->At(kMotherYThrown+1))->Fill(mother->Y());
    ((TH2F*)histos->At(kMotherPtYThrown+1))->Fill(mother->Y(), mother->Pt());
    ((TH1F*)histos->At(kPosLegPtThrown+1))->Fill(partPos->Pt());
    ((TH1F*)histos->At(kPosLegPThrown+1))->Fill(partPos->P());
    ((TH2F*)histos->At(kPosLegPYThrown+1))->Fill(partPos->Y(), partPos->P());
    ((TH1F*)histos->At(kPosLegYThrown+1))->Fill(partPos->Y());
    ((TH2F*)histos->At(kPosLegPtYThrown+1))->Fill(partPos->Y(), partPos->Pt());
    ((TH1F*)histos->At(kNegLegPtThrown+1))->Fill(partNeg->Pt());
    ((TH1F*)histos->At(kNegLegPThrown+1))->Fill(partNeg->P());
    ((TH2F*)histos->At(kNegLegPYThrown+1))->Fill(partNeg->Y(),partNeg->P());
    ((TH1F*)histos->At(kNegLegYThrown+1))->Fill(partNeg->Y());
    ((TH2F*)histos->At(kNegLegPtYThrown+1))->Fill(partNeg->Y(), partNeg->Pt());
    ((TH1F*)histos->At(kPosLegPhiThrown+1))->Fill(partPos->Phi());
    ((TH1F*)histos->At(kNegLegPhiThrown+1))->Fill(partNeg->Phi());
    ((TH2F*)histos->At(kLegCorrelationPxThrown+1))->Fill(partPos->Px(), partNeg->Px());
    ((TH2F*)histos->At(kLegCorrelationPyThrown+1))->Fill(partPos->Py(), partNeg->Py());
    ((TH2F*)histos->At(kLegCorrelationPzThrown+1))->Fill(partPos->Pz(), partNeg->Pz());
    ((TH2F*)histos->At(kLegCorrelationPtThrown+1))->Fill(partPos->Pt(), partNeg->Pt());
    ((TH1F*)histos->At(kCosThetaCSThrown+1))->Fill(cosThetaCS);
    ((TH1F*)histos->At(kCosThetaHEThrown+1))->Fill(cosThetaHE);
    ((TH1F*)histos->At(kOpAngleThrown+1))->Fill(opangle);
    ((TH2F*)histos->At(kOpAnglevsptThrown+1))->Fill(opangle, partNeg->Pt());
    ((TH2F*)histos->At(kOpAnglevsmotherptThrown+1))->Fill(opangle, mother->Pt());
  }
  ((TProfile*)histos->At(kEffPtProf))->Fill(mother->Pt(), accepted);
  ((TProfile*)histos->At(kEffPosLegP))->Fill(partPos->P(), accepted);
  ((TProfile*)histos->At(kEffNegLegP))->Fill(partNeg->P(), accepted);
  ((TProfile2D*)histos->At(kEffPosLegPY))->Fill(partPos->Y(), partPos->P(), accepted);
  ((TProfile2D*)histos->At(kEffNegLegPY))->Fill(partNeg->Y(), partNeg->P(), accepted);
  ((TProfile*)histos->At(kEffYProf))->Fill(mother->Y(), accepted);
  ((TProfile2D*)histos->At(kEffPtYProf))->Fill(mother->Y(), mother->Pt(), accepted);
  ((TProfile*)histos->At(kEffCosThetaCSProf))->Fill(cosThetaCS, accepted);
  ((TProfile*)histos->At(kEffCosThetaHEProf))->Fill(cosThetaHE, accepted);
  ((TProfile*)histos->At(kEffOpAngle))->Fill(opangle, accepted);
  ((TProfile*)histos->At(kEffIntegratedProf))->Fill(0.5, accepted);
}

//_____________________________________________________________________________________________
void DefineHistograms(TObjArray* array, Int_t parMotherCode, Double_t parHighPt) {
  // Define histograms here

  array->SetOwner();

  // change as needed
  /*Int_t nPtBins=0;
  for(Double_t xBin=0; xBin<3.999; xBin+=0.05) nPtBins++;
  for(Double_t xBin=4.0; xBin<5.999; xBin+=0.1) nPtBins++;
  for(Double_t xBin=6.0; xBin<7.999; xBin+=0.2) nPtBins++;
  for(Double_t xBin=8.0; xBin<10.001; xBin+=0.4) nPtBins++;
  Double_t *ptBinLimits = new Double_t[nPtBins];
  nPtBins=0;
  for(Double_t xBin=0; xBin<3.999; xBin+=0.05) ptBinLimits[nPtBins++] = xBin;
  for(Double_t xBin=4.0; xBin<5.999; xBin+=0.1) ptBinLimits[nPtBins++] = xBin;
  for(Double_t xBin=6.0; xBin<7.999; xBin+=0.2) ptBinLimits[nPtBins++] = xBin;
  for(Double_t xBin=8.0; xBin<10.001; xBin+=0.4) ptBinLimits[nPtBins++] = xBin;
  */
  
  Int_t nPtBins = 8;
  Double_t ptBinLimits[8] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0};
  
  // barrel rapidity bins
  Int_t nYBins=0;
  for(Double_t yBin=-5.0; yBin<5.001; yBin+=0.01) nYBins++;
  Double_t *yBinLimits = new Double_t[nYBins];
  nYBins=0;
  for(Double_t yBin=-5.0; yBin<5.001; yBin+=0.01) yBinLimits[nYBins++] = yBin;
  
  TH1F *hPtMotherThrown = new TH1F("hPtMotherThrown",
				   Form("P_{t} distribution of the mother (%d), thrown", parMotherCode), 
				   nPtBins-1, ptBinLimits);
  array->AddAt(hPtMotherThrown, kMotherPtThrown);
  TH1F *hYMotherThrown = new TH1F("hYMotherThrown", 
				  Form("Rapidity distribution of the mother (%d), thrown", parMotherCode), 
				  nYBins-1, yBinLimits);
  array->AddAt(hYMotherThrown, kMotherYThrown);				  
  TH2F *hPtYMotherThrown = new TH2F("hPtYMotherThrown", 
				    Form("(P_{t},y) distribution of the mother (%d), thrown", parMotherCode),
				    nYBins-1, yBinLimits, nPtBins-1, ptBinLimits);
  array->AddAt(hPtYMotherThrown, kMotherPtYThrown);
  TH1F *hPtPosLegThrown = new TH1F("hPtPosLegThrown","P_{t} positive leg, thrown", 1000, 0., parHighPt+2.0);
  array->AddAt(hPtPosLegThrown, kPosLegPtThrown);
  TH1F *hPPosLegThrown = new TH1F("hPPosLegThrown","P positive leg, thrown", 1000, 0., parHighPt*2.0+2.0);
  array->AddAt(hPPosLegThrown, kPosLegPThrown);
  TH2F *hPYPosLegThrown = new TH2F("hPYPosLegThrown","P,y positive leg, thrown", nYBins-1, yBinLimits, 1000, 0., parHighPt*2.0+2.0);
  array->AddAt(hPYPosLegThrown, kPosLegPYThrown);
  TH1F *hPtNegLegThrown = new TH1F("hPtNegLegThrown","P_{t} negative leg, thrwon", 1000, 0., parHighPt+2.0);
  array->AddAt(hPtNegLegThrown, kNegLegPtThrown);
  TH1F *hPNegLegThrown = new TH1F("hPNegLegThrown","P negative leg, thrown", 1000, 0., parHighPt*2.0+2.0);
  array->AddAt(hPNegLegThrown, kNegLegPThrown);
  TH2F *hPYNegLegThrown = new TH2F("hPYNegLegThrown","P,y negative leg, thrown", nYBins-1, yBinLimits, 1000, 0., parHighPt*2.0+2.0);
  array->AddAt(hPYNegLegThrown, kNegLegPYThrown);
  TH2F *hLegCorrelationPxThrown = new TH2F("hLegCorrelationPxThrown", "positive vs. negative leg correlation for P_{x}, thrown",
				     400, -10., 10., 400, -10., 10.);
  array->AddAt(hLegCorrelationPxThrown, kLegCorrelationPxThrown);
  TH2F *hLegCorrelationPyThrown = new TH2F("hLegCorrelationPyThrown", "positive vs. negative leg correlation for P_{y}, thrown",
				     400, -10., 10., 400, -10., 10.);
  array->AddAt(hLegCorrelationPyThrown, kLegCorrelationPyThrown);
  TH2F *hLegCorrelationPzThrown = new TH2F("hLegCorrelationPzThrown", "positive vs. negative leg correlation for P_{z}, thrown",
				     400, -100., 100., 400, -100., 100.);
  array->AddAt(hLegCorrelationPzThrown, kLegCorrelationPzThrown);
  TH2F *hLegCorrelationPtThrown = new TH2F("hLegCorrelationPtThrown", "positive vs. negative leg correlation for P_{t}, thrown",
                                           200, 0., 10., 200, 0., 10.);
  array->AddAt(hLegCorrelationPtThrown, kLegCorrelationPtThrown);
  TH1F *hYPosLegThrown = new TH1F("hYPosLegThrown","y leg, thrown", 1600, -10.,+10.);
  array->AddAt(hYPosLegThrown, kPosLegYThrown);
  TH1F *hYNegLegThrown = new TH1F("hYNegLegThrown","y negative leg, thrown", 1600, -10.,+10.);
  array->AddAt(hYNegLegThrown, kNegLegYThrown);
  TH2F *hPtYPosLegThrown = new TH2F("hPtYPosLegThrown", "(P_{t},#eta) positive leg, thrown", 
				    1600, -8.,+8.,1000, 0., parHighPt+2.0);
  array->AddAt(hPtYPosLegThrown, kPosLegPtYThrown);				    
  TH2F *hPtYNegLegThrown = new TH2F("hPtYNegLegThrown", "(P_{t},#eta) negative leg, thrown", 
				    1600, -8.,+8.,1000, 0., parHighPt+2.0);
  array->AddAt(hPtYNegLegThrown, kNegLegPtYThrown);
  TH1F *hPhiPosLegThrown = new TH1F("hPhiPosLegThrown", "#phi positive leg, thrown", 600, 0., 6.3);
  array->AddAt(hPhiPosLegThrown, kPosLegPhiThrown);
  TH1F *hPhiNegLegThrown = new TH1F("hPhiNegLegThrown", "#phi negative leg, thrwon", 600, 0., 6.3);
  array->AddAt(hPhiNegLegThrown, kNegLegPhiThrown);
  TH1F *hCosThetaStarCSThrown = new TH1F("hCosThetaStarCSThrown", "cos(#theta^{*}_{CS}), thrown",100,-1,1);	
  array->AddAt(hCosThetaStarCSThrown, kCosThetaCSThrown);
  TH1F *hCosThetaStarHEThrown = new TH1F("hCosThetaStarHEThrown", "cos(#theta^{*}_{HE}), thrown",100,-1,1);
  array->AddAt(hCosThetaStarHEThrown, kCosThetaHEThrown);
  TH1F *hOpAngleThrown = new TH1F("hOpAngleThrown", "hopangle",360,-360,360);
  array->AddAt(hOpAngleThrown, kOpAngleThrown);
  TH1F *hOpAngleAcc = new TH1F("hOpAngleAcc", "hopangle",360,-360,360);
  array->AddAt(hOpAngleAcc, kOpAngleAcc);
  TH2F *hOpAnglevsptThrown = new TH2F("hOpAnglevsptThrown", "hopanglevspt",360,0,360,200,0,20);
  array->AddAt(hOpAnglevsptThrown, kOpAnglevsptThrown);
  TH2F *hOpAnglevsptAcc = new TH2F("hOpAnglevsptAcc", "hopanglevspt",360,0,360,200,0,20);
  array->AddAt(hOpAnglevsptAcc, kOpAnglevsptAcc);
  TH2F *hOpAnglevsmotherptThrown = new TH2F("hOpAnglevsmotherptThrown", "hopanglevsmotherpt",360,0,360,200,0,20);
  array->AddAt(hOpAnglevsmotherptThrown, kOpAnglevsmotherptThrown);
  TH2F *hOpAnglevsmotherptAcc = new TH2F("hOpAnglevsmotherptAcc", "hopanglevsmotherpt",360,0,360,200,0,20);
  array->AddAt(hOpAnglevsmotherptAcc, kOpAnglevsmotherptAcc);
  TH1F *hPtMotherAcc = new TH1F("hPtMotherAcc",
				Form("P_{t} distribution of the mother (%d), accepted", parMotherCode), 
				nPtBins-1, ptBinLimits);
  array->AddAt(hPtMotherAcc, kMotherPtAcc);
  TH1F *hYMotherAcc = new TH1F("hYMotherAcc", 
				Form("Rapidity distribution of the mother (%d), accepted", parMotherCode), 
				nYBins-1, yBinLimits);
  array->AddAt(hYMotherAcc, kMotherYAcc);
  TH2F *hPtYMotherAcc = new TH2F("hPtYMotherAcc", 
				 Form("(P_{t},y) distribution of the mother (%d), accepted", parMotherCode),
				 nYBins-1, yBinLimits, nPtBins-1, ptBinLimits);
  array->AddAt(hPtYMotherAcc, kMotherPtYAcc);
  TH1F *hPtPosLegAcc = new TH1F("hPtPosLegAcc","P_{t} positive leg, accepted", 1000, 0., parHighPt+2.0);
  array->AddAt(hPtPosLegAcc, kPosLegPtAcc);
  TH1F *hPPosLegAcc = new TH1F("hPPosLegAcc","P positive leg, accepted", 1000, 0., parHighPt*2.0+2.0);
  array->AddAt(hPPosLegAcc, kPosLegPAcc);
  TH2F *hPYPosLegAcc = new TH2F("hPYPosLegAcc","P,y positive leg, accepted", nYBins-1, yBinLimits, 1000, 0., parHighPt*2.0+2.0);
  array->AddAt(hPYPosLegAcc, kPosLegPYAcc);
  TH1F *hPtNegLegAcc = new TH1F("hPtNegLegAcc","P_{t} negative leg, accepted", 1000, 0., parHighPt+2.0);
  array->AddAt(hPtNegLegAcc, kNegLegPtAcc);
  TH1F *hPNegLegAcc = new TH1F("hPNegLegAcc","P negative leg, accepted", 1000, 0., parHighPt*2.0+2.0);
  array->AddAt(hPNegLegAcc, kNegLegPAcc);
  TH2F *hPYNegLegAcc = new TH2F("hPYNegLegAcc","P,y negative leg, accepted", nYBins-1, yBinLimits, 1000, 0., parHighPt*2.0+2.0);
  array->AddAt(hPYNegLegAcc, kNegLegPYAcc);
  TH2F *hLegCorrelationPxAcc = new TH2F("hLegCorrelationPxAcc", "positive vs. negative leg correlation for P_{x}, accepted",
				     400, -10., 10., 400, -10., 10.);
  array->AddAt(hLegCorrelationPxAcc, kLegCorrelationPxAcc);
  TH2F *hLegCorrelationPyAcc = new TH2F("hLegCorrelationPyAcc", "positive vs. negative leg correlation for P_{y}, accepted",
				     400, -10., 10., 400, -10., 10.);				     
  array->AddAt(hLegCorrelationPyAcc, kLegCorrelationPyAcc);
  TH2F *hLegCorrelationPzAcc = new TH2F("hLegCorrelationPzAcc", "positive vs. negative leg correlation for P_{z}, accepted",
				     400, -10., 10., 400, -10., 10.);
  array->AddAt(hLegCorrelationPzAcc, kLegCorrelationPzAcc);
  TH2F *hLegCorrelationPtAcc = new TH2F("hLegCorrelationPtAcc", "positive vs. negative leg correlation for P_{t}, accepted",
                                        200, 0., 10., 200, 0., 10.);
  array->AddAt(hLegCorrelationPtAcc, kLegCorrelationPtAcc);
  TH1F *hYPosLegAcc = new TH1F("hYPosLegAcc","#eta positive leg, accepted", 1600, -8.,+8.);
  array->AddAt(hYPosLegAcc, kPosLegYAcc);
  TH1F *hYNegLegAcc = new TH1F("hYNegLegAcc","#eta negative leg, accepted", 1600, -8.,+8.);
  array->AddAt(hYNegLegAcc, kNegLegYAcc);
  TH2F *hPtYPosLegAcc = new TH2F("hPtYPosLegAcc", "(P_{t},#eta) positive leg, accepted", 
				    1600, -8.,+8., 1000, 0., parHighPt+2.0);
  array->AddAt(hPtYPosLegAcc, kPosLegPtYAcc);
  TH2F *hPtYNegLegAcc = new TH2F("hPtYNegLegAcc", "(P_{t},#eta) negative leg, accepted", 
				    1600, -8.,+8., 1000, 0., parHighPt+2.0);
  array->AddAt(hPtYNegLegAcc, kNegLegPtYAcc);
  TH1F *hPhiPosLegAcc = new TH1F("hPhiPosLegAcc", "#phi positive leg, accepted", 600, 0., 6.3);
  array->AddAt(hPhiPosLegAcc, kPosLegPhiAcc);
  TH1F *hPhiNegLegAcc = new TH1F("hPhiNegLegAcc", "#phi negative leg, accepted", 600, 0., 6.3);
  array->AddAt(hPhiNegLegAcc, kNegLegPhiAcc);				    
  TH1F *hCosThetaStarCSAcc = new TH1F("hCosThetaStarCSAcc", "cos(#theta^{*}_{CS}), accepted",100,-1,1);
  array->AddAt(hCosThetaStarCSAcc, kCosThetaCSAcc);
  TH1F *hCosThetaStarHEAcc = new TH1F("hCosThetaStarHEAcc", "cos(#theta^{*}_{HE}), accepted",100,-1,1);
  array->AddAt(hCosThetaStarHEAcc, kCosThetaHEAcc);
  
  TH1F* hPtEff = new TH1F("hPtEff", "pt efficiency", nPtBins-1, ptBinLimits);
  array->AddAt(hPtEff, kEffPt);
  TH1F* hYEff = new TH1F("hYEff", "y efficiency", nYBins-1, yBinLimits);
  array->AddAt(hYEff, kEffY);
  TH2F* hPtYEff = new TH2F("hPtYEff", "(pt,Y) efficiency", 
		           nYBins-1, yBinLimits, nPtBins-1, ptBinLimits);
  array->AddAt(hPtYEff, kEffPtY);
  TH1F *hCosThetaStarEffCS = new TH1F("hCosThetaStarEffCS", "cos(#theta^{*}_{CS}) efficiency",100,-1,1);
  array->AddAt(hCosThetaStarEffCS, kEffCosThetaCS);
  TH1F *hCosThetaStarEffHE = new TH1F("hCosThetaStarEffHE", "cos(#theta^{*}_{HE}) efficiency",100,-1,1);
  array->AddAt(hCosThetaStarEffHE, kEffCosThetaHE);
  
  //TProfile* hPtEffProf = new TProfile("hPtEffProf","pt efficiency", nPtBins-1, ptBinLimits);
  TProfile* hPtEffProf = new TProfile("hPtEffProf","pt efficiency", 1000, 0., 10.);
  array->AddAt(hPtEffProf, kEffPtProf);
  TProfile* hYEffProf = new TProfile("hYEffProf","y efficiency", nYBins-1, yBinLimits);
  array->AddAt(hYEffProf, kEffYProf);
  TProfile* hPEffPosLegProf = new TProfile("hPEffPosLegProf","P efficiency, positive leg", nPtBins-1, ptBinLimits);
  array->AddAt(hPEffPosLegProf, kEffPosLegP);
  TProfile* hPEffNegLegProf = new TProfile("hPEffNegLegProf","P efficiency, negative leg", nPtBins-1, ptBinLimits);
  array->AddAt(hPEffNegLegProf, kEffNegLegP);
  TProfile *hCosThetaStarEffCSProf = new TProfile("hCosThetaStarEffCSProf", "cos(#theta^{*}_{CS}) efficiency",100,-1,1);
  array->AddAt(hCosThetaStarEffCSProf, kEffCosThetaCSProf);
  TProfile *hCosThetaStarEffHEProf = new TProfile("hCosThetaStarEffHEProf", "cos(#theta^{*}_{HE}) efficiency",100,-1,1);
  array->AddAt(hCosThetaStarEffHEProf, kEffCosThetaHEProf);
  TProfile2D* hPtYEffProf = new TProfile2D("hPtYEffProf", "Mother acceptance, (pt,Y)", 
				           nYBins-1, yBinLimits, nPtBins-1, ptBinLimits);
  array->AddAt(hPtYEffProf, kEffPtYProf);
  TProfile2D* hPYPosLegEffProf = new TProfile2D("hPYPosLegEffProf", "Positive leg acceptance, (P,Y)", 
						nYBins-1, yBinLimits, nPtBins-1, ptBinLimits);
  array->AddAt(hPYPosLegEffProf, kEffPosLegPY);
  TProfile2D* hPYNegLegEffProf = new TProfile2D("hPYNegLegEffProf", "Negative leg acceptance, (P,Y)", 
						nYBins-1, yBinLimits, nPtBins-1, ptBinLimits);
  array->AddAt(hPYNegLegEffProf, kEffNegLegPY);

  TProfile* hEffIntegrated = new TProfile("hEffIntegrated", "Integrated efficiency", 1, 0., 1.);
  array->AddAt(hEffIntegrated, kEffIntegratedProf);
  
  TProfile *hEffOpAngle = new TProfile("hEffOpAngle", "hopangle",360,-360,360);
  array->AddAt(hEffOpAngle, kEffOpAngle);
  // ----------------------------------------------------------------------------------------
}
