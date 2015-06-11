{
TString cmsswbase = getenv("CMSSW_BASE");
    if(cmsswbase.Length() > 0){
        cout << "Loading FW Lite Setup" << endl;
        gSystem->Load("libFWCoreFWLite.so");
        AutoLibraryLoader::enable();
        gSystem->Load("libDataFormatsFWLite.so");
        gSystem->Load("libDataFormatsPatCandidates.so");
	gSystem->Load("../../../Headers/RooUnfold-1.1.1/libRooUnfold.so");
    }

    cout<<"Hi Raghav"<<endl;
    cout<<"Please dont get angry/annoyed with me, after all you are only human"<<endl;
    cout<<"good luck :) "<<endl;

}
