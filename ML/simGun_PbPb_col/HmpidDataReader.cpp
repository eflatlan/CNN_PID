#ifndef HMPID_DATA_READER_H
#define HMPID_DATA_READER_H

#include "HmpidDataReader.cpp"

#include "DataFormatsHMP/Cluster.h"
#include "DataFormatsHMP/Digit.h"
#include "DataFormatsHMP/Trigger.h"
#include "HMPIDReconstruction/Clusterer.h"
#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLine.h>
#include <TList.h>
#include <TMath.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TROOT.h> // gRoot

#include <TRandom.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <fairlogger/Logger.h>
#include <fstream>
#include <vector>

#include "CommonDataFormat/InteractionRecord.h"

// C++ header files and libraries
#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <gsl/gsl>
#include <iostream>
#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <sys/stat.h>
#include <thread>
#include <unistd.h>

using o2::hmpid::Cluster, o2::hmpid::Digit, o2::hmpid::Trigger,
    o2::hmpid::Clusterer;
using std::vector, std::cout, std::cin, std::endl;
using std::this_thread::sleep_for;

using Clusters = o2::hmpid::Cluster;
using Cluster = o2::hmpid::Cluster; //, o2::hmpid::Digit, o2::hmpid::Trigger,
                                    //o2::hmpid::Clusterer;

class HmpidDataReader {
private:
    TFile fileKine;
    TFile fileClu;
    TFile fileMatch;
    
    TTree* treeMatch;
    TTree* treeClu;
    TTree* treeKine;
    
    std::vector<o2::MCTrack> mcArr, *mcArrPtr = &mcArr;
    std::vector<o2::dataformats::MatchInfoHMP> mMatches, *matchArrPtr = &mMatches;
  	std::vector<Cluster> cluArr, *cluArrPtr = &cluArr;
  	/*std::vector<o2::hmpid::Topology> mTopologyFromFile,
      *mTopologyFromFilePtr = &mTopologyFromFile;*/ 

  	std::vector<Trigger> *trigArr = nullptr;

    std::vector<o2::MCCompLabel> mLabelHMP, *mLabelHMPPtr = &mLabelHMP;



public:



    HmpidDataReader(const char* matchFileName, const char* cluFileName, const char* mcFileName)
    : fileKine(mcFileName, "READ"), fileClu(cluFileName, "READ"), fileMatch(matchFileName, "READ") {
        // Check if files are opened correctly
        if (fileKine.IsZombie() || fileClu.IsZombie() || fileMatch.IsZombie()) {
            throw std::runtime_error("Failed to open one or more files.");
        }

        // Initialize TTree objects

        	                

        initializeClusterTree();
				initializeMatchTree(0, 0, 0);


				if (!treeClu) {
					Printf("Error accessing TTree");

				}

        // Check if TTree objects are initialized correctly
        if (!treeMatch || !treeClu || !treeKine) {
            throw std::runtime_error("Failed to initialize one or more TTree objects.");
        }

        /*        void HmpidDataReader::initializeClusterTree(/ * std::vector<Cluster> *&cluArrPtr, std::vector<Trigger> *&trigArr* / ) {
          // TTree *treeClu = (TTree *)fileClu.Get("o2sim");
          
          if (!treeClu)
            treeClu = (TTree *)fileClu.Get("o2hmp");
          if (!treeClu) {
            Printf("Error accessing TTree");
            fileClu->Close();
            delete fileClu;
            return nullptr;
          }
          treeClu->Print("toponly");

          treeClu->SetBranchAddress("HMPIDclusters", &cluArrPtr);
          treeClu->SetBranchAddress("InteractionRecords", &trigArr);

          treeClu->GetEntry(0);
          // return treeClu;
        } */

				// int eventID, int trackID, int pdg
    }


  ~HmpidDataReader() {

  }

	std::vector<Cluster> getClusInEvent(int event) const {
    auto pTgr = &trigArr->at(event);
    std::vector<Cluster>  oneEventClusters;
		for (int j = pTgr->getFirstEntry(); j <= pTgr->getLastEntry(); j++) {
			const auto &clu = static_cast<o2::hmpid::Cluster>(cluArrPtr->at(j));
			oneEventClusters.push_back(clu);      
		}    

	 	return oneEventClusters; 
	}
	
	std::vector<Cluster> * getcluArrPtr() const { return cluArrPtr ;}

	std::vector<Trigger> * getTrigArr() const {return trigArr ;}


  void initializeClusterTree(/*std::vector<Cluster> *&cluArrPtr, std::vector<Trigger> *&trigArr*/ );


  void initializeMatchTree(int eventID, int trackID, int pdg);
  void readMatch(int eventID, int &startIndex, std::vector<o2::dataformats::MatchInfoHMP>& filteredMatches, std::vector<o2::MCCompLabel>& filteredLblMatches);

  // std::vector<o2::dataformats::MatchInfoHMP> readMatch(int eventID, int &startIndex);
  /*static TTree *initializeClusterTree(
      std::vector<Cluster> *&cluArrPtr, std::vector<Trigger> *&trigArr);*/ 



  TTree *initializeMCTree(std::vector<o2::MCTrack> *&mcArrPtr);
  std::vector<o2::MCTrack> *readMC(std::vector<o2::MCTrack> *&mcArrPtr,
                                          TTree *treeKine, int eventId);


  static const o2::MCTrack *getMCEntry(std::vector<o2::MCTrack> *mcArrPtr,
                                       int trackID);
  static void readTreeEntries();
};

void HmpidDataReader::initializeMatchTree(int eventID,
    int trackID, int pdg) {
    
   
  // std::unique_ptr<TFile> fileMatch(TFile::Open("o2match_hmp.root", "READ"));
  
  treeMatch = dynamic_cast<TTree*>(fileMatch.Get("matchHMP"));

  if (!treeMatch)
    treeMatch = dynamic_cast<TTree*>(fileMatch.Get("o2hmp"));
    
  if (!treeMatch)
    treeMatch = dynamic_cast<TTree*>(fileMatch.Get("o2sim"));

  if (!treeMatch)
		throw std::runtime_error("treeMatch nullptr");



  // std::vector<o2::dataformats::MatchInfoHMP>* matchArrPtr = nullptr;
  treeMatch->SetBranchAddress("HMPMatchInfo", &matchArrPtr);

  if (mUseMC) {
    treeMatch->SetBranchAddress("MatchHMPMCTruth", &mLabelHMPPtr);
  }

  treeMatch->GetEntry(0);
  treeMatch->Print("toponly");
  if (matchArrPtr == nullptr) {
    Printf("HmpidDataReader::initializeMatchTree matchArrPtr== nullptr");
		// fileMatch.Close();
		throw std::runtime_error("matchArrPtr nullptr");
  }

  if(mUseMC && !mLabelHMPPtr) {
    Printf("HmpidDataReader::initializeMatchTree mLabelHMPPtr== nullptr");
		// fileMatch.Close();
		throw std::runtime_error("mLabelHMPPtr nullptr");
  }
}

// eventId = eventID to be searched for;
// startIndex : index of where matchArrPtr is to be searched
// newStartIndex startIndex for next event

void HmpidDataReader::readMatch(int eventID, int &startIndex, std::vector<o2::dataformats::MatchInfoHMP>& filteredMatches, std::vector<o2::MCCompLabel>& filteredLblMatches) {
               
               
  Printf("Call HmpidDataReader::readMatch");                           
	//std::vector<o2::dataformats::MatchInfoHMP> filteredMatches;// = new std::vector<o2::dataformats::MatchInfoHMP>;                           
                           
  if (!treeMatch) {
    Printf("TTree not initialized");
    return;// filteredMatches;
  } else {
      Printf("HmpidDataReader::readMatch : TTree  initialized");
  }

  // Prepare to store filtered matches
  // std::vector<o2::dataformats::MatchInfoHMP> filteredMatches;// = new std::vector<o2::dataformats::MatchInfoHMP>;
  // tracks should be stored in "time" --> when we find our event we can then
  // switch this condition "of" when the event changes:

  bool found = false;

  if (matchArrPtr == nullptr) {
    Printf("HmpidDataReader::readMatch :: matchArrPtr== nullptr");
    return;// filteredMatches;
  } else {
      Printf("HmpidDataReader::readMatch : matchArrPtr ok");
  } 
  
  Printf("readMatch : (*matchArrPtr) size : %zu ", (*matchArrPtr).size());
  Printf("readMatch : startIndex %d", startIndex);

         
  if((*matchArrPtr).size() < 1) {
  	LOGP(info, "matchArrPtr was 0");
    return;// filteredMatches;
  }       
         
         
  Printf("readMatch : (*matchArrPtr)[startIndex].getEvent() %d eventID %d",
         (*matchArrPtr)[startIndex].getEvent(), eventID);
         
  if ((*matchArrPtr)[startIndex].getEvent() != eventID) {
    Printf("This shouldnt happen");
  } else
    found = true;

  Printf("readMatch : event %d startIndex %d", eventID, startIndex);
	LOGP(info, "matchArrPtr->size() {}", matchArrPtr->size());

	LOGP(info, "mLabelHMPPtr->size() {}", mLabelHMPPtr->size());

  for (int i = startIndex; i < matchArrPtr->size(); i++) {
    const auto &track = (*matchArrPtr)[1];
    const auto &lblTrk = (*mLabelHMPPtr)[1];

    // filteredMatches->push_back(track);
    // startIndex = i;
    //  ef : we cant use this atm, since the clusters from the same trigger
    //  sometimes have different eventNumber!
    
    
    
    LOGP(info, "trackEvent {}, i {} | getMipClusEvent {}", track.getEvent(), i, track.getMipClusEvent());
    
    
    // ef: TODO is this now needed when checking MCtruth from MClabel?
    if (track.getEvent() != eventID) {
      startIndex = i; // new startIndex for next event
      Printf("readMatch : eventID changed - %d; end of loop ",
             track.getEvent());

      break;
    } else {

      filteredMatches.push_back(track);
      filteredLblMatches.push_back(lblTrk);
    }
  }

  Printf("readMatch : new startIndex %d", startIndex);

}

void HmpidDataReader::initializeClusterTree(/*std::vector<Cluster> *&cluArrPtr, std::vector<Trigger> *&trigArr*/ ) {
  // TTree *treeClu = (TTree *)fileClu.Get("o2sim");
  treeClu = dynamic_cast<TTree*>(fileClu.Get("o2hmp"));
  if (!treeClu)
    treeClu = dynamic_cast<TTree*>(fileClu.Get("o2sim"));

  if (!treeClu) {
    Printf("Error accessing TTree");
  }
  
  treeClu->Print("toponly");

  treeClu->SetBranchAddress("HMPIDclusters", &cluArrPtr);
  treeClu->SetBranchAddress("InteractionRecords", &trigArr);

  // ef : add MC clus information
  if(treeClu->GetBranchStatus("") {
    treeClu->SetBranchAddress("HMPIDclusters", &cluArrPtr);

  }


  treeClu->GetEntry(0);
  // return treeClu;
} 

TTree *HmpidDataReader::initializeMCTree(std::vector<o2::MCTrack> *&mcArrPtr) {

  treeKine = dynamic_cast<TTree*>(fileKine.Get("o2sim"));
  if (!treeKine)
    treeKine = dynamic_cast<TTree*>(fileKine.Get("o2sim"));

  if (!treeKine) {
    Printf("Error accessing TTree");
    /*fileKine->Close();
    delete fileKine;*/
    return nullptr;
  }


  treeKine->SetBranchAddress("MCTrack", &mcArrPtr);
  //treeKine->SetBranchAddress("TrackRefs", &mcArrPtr);

  treeKine->GetEntry(0);
  treeKine->Print("toponly");
  return treeKine;
}

std::vector<o2::MCTrack> *
HmpidDataReader::readMC(std::vector<o2::MCTrack> *&mcArrPtr, TTree *treeKine,
                        int eventID) {
  if (treeKine == nullptr) {
    Printf("Error : treeKine == nullptr");
    return nullptr;
  }

  if (eventID < treeKine->GetEntries()) {
    treeKine->GetEntry(eventID);
    Printf("readMC at entry %d", eventID);
    return mcArrPtr;
  } else {
    Printf("eventId > treeKine->GetEntries()");
    return nullptr;
  }
}

// for the given eventID; read trackID
const o2::MCTrack *HmpidDataReader::getMCEntry(std::vector<o2::MCTrack> *mcArrPtr,
                                               int trackID) {

  if (trackID < 0 || trackID >= mcArrPtr->size()) {
    return nullptr;
  } else {
    // const auto& track = mcArrPtr->at(trackID);
    for (int i = 0; i < mcArrPtr->size(); ++i) {
      const auto &mcTrack = (*mcArrPtr)[i];
      if (i == trackID) {

        TParticlePDG* pPDG = TDatabasePDG::Instance()->GetParticle(mcTrack->GetPdgCode());
        auto pT = TMath::Abs(mcTrack.GetStartVertexMomentumX() *
                              mcTrack.GetStartVertexMomentumX() +
                          mcTrack.GetStartVertexMomentumY() *
                              mcTrack.GetStartVertexMomentumY());
        auto pdgCode = mcTrack.GetPdgCode();
        LOGP(info, "Particle {}: pdg = {}, pT = {}", i, pdgCode, pT);

        /*
        Printf("Particle %d: pdg = %d, pT = %f, px = %f, py = %f, pz = %f, vx "
               "= %f, vy = %f, vz = %f",
               i, mcTrack.GetPdgCode(),
               TMath::Abs(mcTrack.GetStartVertexMomentumX() *
                              mcTrack.GetStartVertexMomentumX() +
                          mcTrack.GetStartVertexMomentumY() *
                              mcTrack.GetStartVertexMomentumY()),
               mcTrack.GetStartVertexMomentumX(),
               mcTrack.GetStartVertexMomentumY(),
               mcTrack.GetStartVertexMomentumZ(),
               mcTrack.GetStartVertexCoordinatesX(),
               mcTrack.GetStartVertexCoordinatesY(),
               mcTrack.GetStartVertexCoordinatesZ());*/
        return &mcTrack;
      }
    }
  }
  return nullptr;
}

void HmpidDataReader::readTreeEntries() {
  // Open the ROOT file
  /*
  int i;
  auto trackVector = readTrack(0,0,i);
  auto clusterVector = readClu(0,0,i);
  std::vector<o2::MCTrack>* mcVector = readMC(0,0,i);
  for() {

  }

        Printf("numTracks %d numClusters %d numMC %d", trackVector->size(),
  clusterVector->size(), mcVector->size()); */
}

#endif // HMPID_DATA_READER_H
