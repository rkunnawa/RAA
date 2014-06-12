{
TString cmsswbase = getenv("CMSSW_BASE");
    if(cmsswbase.Length() > 0){
        cout << "Loading FW Lite Setup" << endl;
        gSystem->Load("libFWCoreFWLite.so");
        AutoLibraryLoader::enable();
        gSystem->Load("libDataFormatsFWLite.so");
        gSystem->Load("libDataFormatsPatCandidates.so");
	gSystem->Load("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_8_HI_patch2/src/Headers/RooUnfold-1.1.1/libRooUnfold.so");
    }

cout<<"HI Raghav, Try to not get frustrated with me"<<endl;
cout<<"afterall you are only human"<<endl;
cout<<"have a nice day!! :)"<<endl;
}
