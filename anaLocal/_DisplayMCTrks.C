// _DisplayMCTrks.C

#include <iostream>
// root headers
#include "TSystem.h"
// aliroot headers
#include "AliRunLoader.h"
// my headers
#include "FocalUpcAnalysis_Utilities.h"

void _DisplayMCTrks()
{
    z_min = -10.;
    z_max = 890.;
    TString sSubfolder = "sim01/kCohJpsiToElRad_001_1000ev/";

    gSystem->Load("libpythia6_4_28.so"); 
    gSystem->Exec(("mkdir -p results/displayMCTrks/" + sSubfolder).Data());

    // define ALICE run loader: open galice.root
    AliRunLoader *runLoader = NULL;
    runLoader = AliRunLoader::Open(("inputData/" + sSubfolder + "galice.root").Data());
    if(!runLoader) 
    {
        cout << "_DisplayMCTrks() ERROR: AliRunLoader not good! Terminating." << endl;
        return;   
    }
    if(!runLoader->GetAliRun()) runLoader->LoadgAlice();
    if(!runLoader->TreeE()) runLoader->LoadHeader();
    if(!runLoader->TreeK()) runLoader->LoadKinematics();

    // loop over MC events contained within Kinematics.root
    //for(Int_t iEv = 0; iEv < runLoader->GetNumberOfEvents(); iEv++) 
    for(Int_t iEv = 0; iEv < 10; iEv++) 
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
        TCanvas* c = PrepareCanvas();
        DrawTracksMC(stack,c,kFALSE);
        c->SaveAs(Form("results/displayMCTrks/%sEv%i.pdf",sSubfolder.Data(),iEv));
        delete c;
    }

    return;
}