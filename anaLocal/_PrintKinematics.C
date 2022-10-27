// _PrintKinematics.C
// David Grund, Oct 20, 2022

// cpp headers
#include "fstream"
#include "stdio.h"
#include "iomanip" // std::setprecision()
// aliroot headers
#include "AliRunLoader.h"
#include "AliStack.h"
#include "AliFOCALCluster.h"
// root headers
#include "TString.h"
#include "TParticle.h"

void _PrintKinematics() 
{
    TString dir = "inputData/sim02/kIncohJpsiToElRad_001_1000ev/"; // add subfolder, if needed
    //dir = "inputData/sim01/BoxElectrons_001_1000ev/";
    // get the directory
    TString galiceFile = dir.Data();
    galiceFile += "galice.root";

    // open the working directory
    AliRunLoader* runLoader = AliRunLoader::Open(galiceFile);
    runLoader->LoadKinematics();

    // create output text file
    ofstream out_txt;
    out_txt.open(dir + "printKinematics.txt");

    // loop over events
    Int_t nEvents = runLoader->GetNumberOfEvents();
    out_txt << " ++++ TOTAL NO. OF EVENTS: " << nEvents << endl;

    for(Int_t iEvent = 0; iEvent < nEvents; iEvent++)
    {
        if((iEvent+1 % 1000) == 0) Printf("%i events printed.", iEvent+1);

        out_txt << " ++++ EVENT no. " << iEvent+1 << endl;

        runLoader->GetEvent(iEvent);	    
        AliStack* stack = runLoader->Stack();
        Int_t nTracks = stack->GetNtrack();
        out_txt << " ++++ NO. OF TRACKS:" << nTracks << endl;
        out_txt << "iTrk\tmother\tpdg\tname\t1stDgh\tVz[cm]\tE[GeV]\tm[MeV]\n";

        // loop over tracks in the event
        for(Int_t iTrk = 0; iTrk < nTracks; iTrk++)
        {
            // get particle from the current track
            TParticle *part = stack->Particle(iTrk);

            // print info to the text file
            out_txt << iTrk << "\t"
                    << part->GetMother(0) << "\t"
                    << part->GetPdgCode() << "\t"
                    << part->GetName() << "\t"
                    << part->GetDaughter(0) << "\t"
                    << std::fixed << std::setprecision(1)
                    << part->Vz() << "\t"
                    << part->Energy() << "\t"
                    << std::fixed << std::setprecision(3)
                    << part->GetCalcMass()*1e3 << endl;
        } // track loop
    } // event loop

    runLoader->UnloadAll();
    out_txt.close();

    return;
}