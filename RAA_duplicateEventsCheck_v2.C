#include "HiForestAnalysis/hiForest.h"
#include <TFile.h>
#include <TH1D.h>
#include <TNtuple.h>
#include <iostream>
#include <unordered_set>


// Example of forest skim

void removeDuplicates(char *infname = "/mnt/hadoop/cms/store/user/dgulhan/HIMC/MB/Track8_Jet26_STARTHI53_LV1/merged2/HiForest_HYDJET_Track8_Jet26_STARTHI53_LV1_merged_forest_0.root", char *outfname = "skim_Jet.root")
{
  // Define the input file and HiForest
  HiForest *c = new HiForest(infname);
  c->hasHltTree=0;
  c->hasPFTree=0;
  c->hasPhotonTree=0;
  c->hasTowerTree=0;
  c->hasHbheTree=0;
  c->hasEbTree=0;
  c->hasGenpTree=0;
  c->hasGenParticleTree=0;   
  c->hasAk5CaloJetTree=0;
  c->hasAkPu2CaloJetTree=0;
  c->hasAkPu3CaloJetTree=0;
  c->hasAkPu4CaloJetTree=0;
  c->hasAkPu5CaloJetTree=0;
  c->hasAkPu2JetTree=0;
  c->hasAkPu3JetTree=0;
  c->hasAkPu4JetTree=0;
  c->hasAkPu5JetTree=0;
  c->hasAkVs2PFJetTree=0;
  c->hasAkVs3PFJetTree=0;
  c->hasAkVs4PFJetTree=0;
  c->hasAkVs5PFJetTree=0;
  c->SetOutputFile(outfname);

  int filtered=0;
  std::unordered_set<long long int> visitedevents;
  // Main loop
  for (int i=0;i<c->GetEntries();i++)
  {
    c->GetEntry(i);
    if (i%1000==0) cout <<filtered<<" "<<i<<" / "<<c->GetEntries()<<endl;
    long long big = 1000000000;
    long long int thisid = (big*c->evt.run) + c->evt.evt; // format is run|0|evt 
    auto search = visitedevents.find(thisid);
    if(search != visitedevents.end()) {
      cout<<"this data sample has duplicate events :( , but we're not analyzing them :) "<<endl;
      filtered++;
      continue; // found duplicate
    }
    else // no duplicate found, add this to visited events
    {
      visitedevents.insert(thisid);
    }




    // int flag=0;    // # of jets with |eta|<2 and pt > 120
    // int flag2=0;   // # of jets with |eta|<2 and pt > 50
    // for (int j=0;j<c->akVs3Calo.nref;j++) {
    //   if (fabs(c->akVs3Calo.jteta[j])>2) continue;
    //   if (c->akVs3Calo.jtpt[j]>120) flag=1;
    //   if (c->akVs3Calo.jtpt[j]>50) flag2++;
    //   if (flag>=1&&flag2>=2) break;
    // }
    // if (flag>=1&&flag2>=2) {
    //   filtered++;
    // }  
    c->FillOutput(); // Write output forest
  }

  delete c;
}



int main(int argc, char *argv[])
{
  if(argc != 3)
  {
    std::cout << "Usage: ./removeDuplicates.exe <input> <output>" << std::endl;
    return 1;
  }
  removeDuplicates(argv[1], argv[2]);

  return 0;
}
