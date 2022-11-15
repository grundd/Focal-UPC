// PrintKinematics.C
// David Grund, Oct 20, 2022

// cpp headers
#include "fstream"
#include "stdio.h"
#include "iomanip" // std::setprecision()
// root headers
#include "TString.h"
#include "TParticle.h"
// aliroot headers
#include "AliRunLoader.h"
#include "AliStack.h"
#include "AliFOCALCluster.h"

void PrintKinematics() 
{
    TString dir = "inputData/starlight/CohFD/"; // add subfolder, if needed
    // get the directory
    TString galiceFile = dir.Data();
    galiceFile += "galice.root";

    // open the working directory
    AliRunLoader* runLoader = AliRunLoader::Open(galiceFile);
    runLoader->LoadKinematics();

    // create output text file
    ofstream of;
    of.open(dir + "printKinematics.txt");

    // loop over events
    Int_t nEvents = runLoader->GetNumberOfEvents();
    of << " ++++ TOTAL NO. OF EVENTS: " << nEvents << endl;

    for(Int_t iEvent = 0; iEvent < 10; iEvent++)
    {
        if((iEvent+1 % 1000) == 0) Printf("%i events printed.", iEvent+1);

        of << " ++++ EVENT no. " << iEvent+1 << endl;

        runLoader->GetEvent(iEvent);	    
        AliStack* stack = runLoader->Stack();
        Int_t nTracks = stack->GetNtrack();
        of << " ++++ NO. OF TRACKS:" << nTracks << endl;
        of << "iTrk\tmother\tpdg\tname\t1stDgh\tVz[cm]\tE[GeV]\tm[MeV]\n";

        // loop over tracks in the event
        for(Int_t iTrk = 0; iTrk < nTracks; iTrk++)
        {
            // get particle from the current track
            TParticle *part = stack->Particle(iTrk);

            // print info to the text file
            of << iTrk << "\t"
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
    of.close();

    return;
}