// DisplayTracksMC.C

#include <iostream>
// root headers
#include "TSystem.h"
// aliroot headers
#include "AliRunLoader.h"
// my headers
#include "Utilities.h"

void DisplayTracksMC(Int_t iAnalysis)
{
    z_min = -10.;
    z_max = 890.;
    if(iAnalysis == 0) subfolder = "09-26-2022_kCohJpsiToElRad_1000ev/";
    if(iAnalysis == 1) subfolder = "10-05-2022_kCohJpsiToElRad_3000ev/";

    gSystem->Load("libpythia6_4_28.so"); 
    gSystem->Exec(("mkdir -p " + subfolder + "tracksMC/").Data());

    // define ALICE run loader: open galice.root
    AliRunLoader *runLoader = NULL;
    runLoader = AliRunLoader::Open((simFolder + subfolder + "galice.root").Data());
    if(!runLoader) 
    {
        cout << "DisplayTracksMC() ERROR: AliRunLoader not good! Terminating." << endl;
        return;   
    }
    if(!runLoader->GetAliRun()) runLoader->LoadgAlice();
    if(!runLoader->TreeE()) runLoader->LoadHeader();
    if(!runLoader->TreeK()) runLoader->LoadKinematics();

    // loop over MC events contained within Kinematics.root
    for(Int_t iEv = 0; iEv < runLoader->GetNumberOfEvents(); iEv++) 
    //for(Int_t iEv = 0; iEv < 10; iEv++) 
    {
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

        // plot the tracks
        TCanvas* c = EventDisplay_PrepareCanvas();
        EventDisplay_PlotMCTracks(stack,c,kFALSE);
        c->SaveAs(Form("%stracksMC/Ev%i.pdf",subfolder.Data(),iEv));
        delete c;
    }

    return;
}