// _ReadKinematics.C
// David Grund, Sep 23, 2022

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

void _ReadKinematics() 
{
    TString dir = ""; // add subfolder, if needed
    // get the directory
    TString galiceFile = dir.Data();
    galiceFile += "galice.root";

    // open the working directory
    AliRunLoader* runLoader = AliRunLoader::Open(galiceFile);
    runLoader->LoadKinematics();

    // create output text file
    ofstream out_txt;
    out_txt.open(dir + "ParticlesInKinematics.txt");

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
        out_txt << "iTrack\tmother\tpdg\tname\t1st dg\teta \tpx \tpy \t pz \tmass \n";

        // loop over tracks in the event
        for(Int_t iTrk = 0; iTrk < nTracks; iTrk++)
        {
            // get particle from the current track
            TParticle *part = stack->Particle(iTrk);

            // print info to the text file
            out_txt << iTrk
                    << "\t" << part->GetMother(0)
                    << "\t" << part->GetPdgCode()
                    << "\t" << part->GetName() 
                    << "\t" << part->GetDaughter(0)
                    << std::fixed << std::setprecision(2)
                    << "\t" << part->Eta()
                    << std::fixed << std::setprecision(3) 
                    << "\t" << part->Px()
                    << "\t" << part->Py() 
                    << "\t" << part->Pz()
                    << std::fixed << std::setprecision(5)
                    << "\t" << part->GetCalcMass() << endl;
        } // track loop
    } // event loop

    runLoader->UnloadAll();
    out_txt.close();

    return;
}