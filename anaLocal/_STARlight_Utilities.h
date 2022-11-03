// _STARlight_Utilities.h
// David Grund, Nov 02, 2022

// cpp headers
#include <iostream>
#include <fstream> 
#include <sstream> 
// root headers
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

TLorentzVector* parent;
TClonesArray* daughters;

// ******************************************************************************************************************
// To connect branch addresses of STARlight trees:
// ******************************************************************************************************************

void SetBranchAddresses_tSL(TTree *t)
{
    t->SetBranchAddress("parent", &parent);
    t->SetBranchAddress("daughters", &daughters);
    Printf("Branch addresses of %s set.", t->GetName());
    return;
}

// ******************************************************************************************************************
// To create tree from STARlight output:
// ******************************************************************************************************************

double IDtoMass(int particleCode)
{
    double mass;
    if (particleCode == 2 || particleCode==3) {mass = 0.00051099907;} // electron
    else if (particleCode == 5 || particleCode==6) {mass = 0.105658389;} // muon
    else if (particleCode == 8 || particleCode==9)  {mass = 0.13956995;} // charged pion
    else if (particleCode == 7) {mass = 0.1345766;} // neutral pion
    else if (particleCode == 11|| particleCode==12) {mass = 0.493677;} // charged kaon
    else if (particleCode == 10 || particleCode == 16)  {mass = 0.497614;} // neutral kaon
    else if (particleCode == 14)	{mass = 0.93827231;} // proton
    else {
        std::cout << "unknown daughter particle (ID = " << particleCode << "), please modify code to accomodate" << std::endl;
        mass = -1.0;
        //exit(0); 
    } 
    return mass;
}

void ConvertStarlightAsciiToTree(Int_t nGenEv, TString sInput)
{
	TString sOut = sInput + "tSTARlight.root";
    TFile* fOut = TFile::Open(sOut.Data(),"read");
    if(fOut)
    {
        Printf("STARlight tree %s already created.", sOut.Data());
        return;
    } 
    else 
    {   
        Printf("STARlight tree %s will be created.", sOut.Data());

		// create the output file and tree
		fOut = new TFile(sOut.Data(), "RECREATE");
		if(!fOut){
			Printf("Could not create output file %s. Terminating...", sOut.Data());
			return;
		}

		TTree*          outTree           = new TTree("starlightTree", "starlightTree");
		TLorentzVector* parentParticle    = new TLorentzVector();
		TClonesArray*   daughterParticles = new TClonesArray("TLorentzVector");
		outTree->Branch("parent",    "TLorentzVector", &parentParticle,    32000, -1);
		outTree->Branch("daughters", "TClonesArray",   &daughterParticles, 32000, -1);

		Int_t nEnAnalyzed = 0;
		Int_t nEnProgress = (Double_t)nGenEv / 100.;
		Int_t nPerc = 0;
		Int_t i = 0;

		ifstream ifs;
		ifs.open((sInput + "slight.out").Data());
		unsigned int countLines = 0;
		while (ifs.good()) {
			string       line;
			stringstream lineStream;
			
			// read EVENT
			string label;
			int    eventNmb, nmbTracks;
			// no more lines => end
			if (!getline(ifs, line))
				break;
			++countLines;
			lineStream.str(line);
			lineStream >> label >> eventNmb >> nmbTracks;
			//Printf("%s", line.data()); // DGRUND
			//cout << countLines << "\t" << eventNmb << "\t" << nmbTracks << endl; // DGRUND
			if (!(label == "EVENT:"))
				continue;
			
			// read VERTEX
			// no more lines => end
			if (!getline(ifs, line))
				break;
			++countLines;
			lineStream.str(line);
			lineStream >> label;
			//Printf("%s", line.data()); // DGRUND
			assert(label == "VERTEX:");
				
			*parentParticle = TLorentzVector(0, 0, 0, 0);
			for (int i = 0; i < nmbTracks; ++i) {
				// read tracks
				int    particleCode;
				double momentum[3];
				// no more lines => end
				if (!getline(ifs, line))
					break;
				++countLines;
				lineStream.str(line);
				lineStream >> label >> particleCode >> momentum[0] >> momentum[1] >> momentum[2];
				//Printf("%s", line.data()); // DGRUND
				assert(label == "TRACK:");
				Double_t daughterMass = IDtoMass(particleCode);
				if (daughterMass < 0) {break;}
				const double E = sqrt(  momentum[0] * momentum[0] + momentum[1] * momentum[1]
									+ momentum[2] * momentum[2] + daughterMass * daughterMass);
				new ( (*daughterParticles)[i] ) TLorentzVector(momentum[0], momentum[1], momentum[2], E);
				*parentParticle += *(static_cast<TLorentzVector*>(daughterParticles->At(i)));
			}
			daughterParticles->Compress();
			outTree->Fill();

			if((i+1) % nEnProgress == 0){
			nPerc += 1;
			nEnAnalyzed += nEnProgress;
			Printf("[%i%%] %i entries analysed.", nPerc, nEnAnalyzed);
			}
			i++;
		}

		outTree->Write("",TObject::kWriteDelete);
		if(fOut){
			fOut->Close();
			delete fOut;
		}

		Printf("\n*****");
		Printf("Done.");
		Printf("*****\n");
		return;
	} 
}
//###############################################################################