// FocalUpcGrid_Config.C 
// David Grund, Nov 07, 2022

// configuration of the FoCal UPC analysis:

// superclusterizer
Bool_t doSupercls = kTRUE;
// if true, create superclusters using:
const Float_t minSeedE = 5; // [GeV]
const Float_t radius = 8; // [cm]
TString sCl = "";

// selections:
// applied on (sup)clusters:
const Float_t cutE = 10.0; // [GeV]; if > 0, filter out all (super)cls having energy below cutE
const Float_t cutM = 0.0; // [GeV]; if > 0, filter out all (super)cls having mass below cutM when paired to any other (sup)cl
// applied on (sup)cl pairs:
const Float_t cutMLow = 1.4; // [GeV]
const Float_t cutMUpp = 4.2; // [GeV]
// pT distribution of (sup)cl pairs:
const Float_t cutMLowPtDist = 2.8; // [GeV]
const Float_t cutMUppPtDist = 3.4; // [GeV]

// matching to MC particles:
Bool_t matchDirectly = kFALSE;
// if true, match every cl to the closest physical primary MC track
// otherwise, match it to the closest track of any type and then find the physical primary track (...)
// which is a mother of the matched MC track
// selections used when matching:
const Float_t cutdEta = 0.4; // [-] difference in eta of a MC particle and a (super)cluster
const Float_t cutdPhi = 0.4; // [-] difference in phi angles of a MC particle and a (super)cluster

TString CreateOutputSubDir()
{
    // set cl name
    if(doSupercls) sCl = "supcl";
    else           sCl = "cl";
    // output dir
    TString outSubDir = "";
    outSubDir = "output";
    if(doSupercls)    outSubDir += Form("_supCl");
    if(matchDirectly) outSubDir += Form("_dirMtch");
    if(cutM > 0)      outSubDir += Form("_cutM%.1f",cutM);
    if(cutE > 0)      outSubDir += Form("_cutE%.1f",cutE);
    outSubDir += "/";
    return outSubDir;
}