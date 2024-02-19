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
    
    std::vector<o2::MCTrack> *mcArr = nullptr;
    std::vector<o2::dataformats::MatchInfoHMP> *matchArr = nullptr;
  	std::vector<Cluster> *cluArr = nullptr;
  	/*std::vector<o2::hmpid::Topology> mTopologyFromFile,
      *mTopologyFromFilePtr = &mTopologyFromFile;*/ 

  	std::vector<Trigger> *trigArr = nullptr;



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

        /*        void HmpidDataReader::initializeClusterTree(/ * std::vector<Cluster> *&cluArr, std::vector<Trigger> *&trigArr* / ) {
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

          treeClu->SetBranchAddress("HMPIDclusters", &cluArr);
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
			const auto &clu = static_cast<o2::hmpid::Cluster>(cluArr->at(j));
			oneEventClusters.push_back(clu);      
		}    

	 	return oneEventClusters; 
	}
	
	std::vector<Cluster> * getCluArr() const { return cluArr ;}

	std::vector<Trigger> * getTrigArr() const {return trigArr ;}


  void initializeClusterTree(/*std::vector<Cluster> *&cluArr, std::vector<Trigger> *&trigArr*/ );


  void initializeMatchTree(int eventID, int trackID, int pdg);

  std::vector<o2::dataformats::MatchInfoHMP> readMatch(int eventID, int &startIndex);
  /*static TTree *initializeClusterTree(
      std::vector<Cluster> *&cluArr, std::vector<Trigger> *&trigArr);*/ 



  TTree *initializeMCTree(std::vector<o2::MCTrack> *&mcArr);
  std::vector<o2::MCTrack> *readMC(std::vector<o2::MCTrack> *&mcArr,
                                          TTree *treeKine, int eventId);
  static const o2::MCTrack *getMCEntry(std::vector<o2::MCTrack> *mcArr,
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



  // std::vector<o2::dataformats::MatchInfoHMP>* matchArr = nullptr;
  treeMatch->SetBranchAddress("HMPMatchInfo", &matchArr);
  treeMatch->GetEntry(0);
  treeMatch->Print("toponly");
  if (matchArr == nullptr) {
    Printf("HmpidDataReader::initializeMatchTree matchArr== nullptr");
		// fileMatch.Close();
		throw std::runtime_error("matchArr nullptr");

  }
}

// eventId = eventID to be searched for;
// startIndex : index of where matchArr is to be searched
// newStartIndex startIndex for next event
std::vector<o2::dataformats::MatchInfoHMP> 
HmpidDataReader::readMatch(int eventID, int &startIndex) {
               
               
  Printf("Call HmpidDataReader::readMatch");                           
	std::vector<o2::dataformats::MatchInfoHMP> filteredMatches;// = new std::vector<o2::dataformats::MatchInfoHMP>;                           
                           
  if (!treeMatch) {
    Printf("TTree not initialized");
    return filteredMatches;
  } else {
      Printf("HmpidDataReader::readMatch : TTree  initialized");
  }

  // Prepare to store filtered matches
  // std::vector<o2::dataformats::MatchInfoHMP> filteredMatches;// = new std::vector<o2::dataformats::MatchInfoHMP>;
  // tracks should be stored in "time" --> when we find our event we can then
  // switch this condition "of" when the event changes:

  bool found = false;

  if (matchArr == nullptr) {
    Printf("HmpidDataReader::readMatch :: matchArr== nullptr");
    return filteredMatches;
  } else {
      Printf("HmpidDataReader::readMatch : matchArr ok");
  } 
  
  Printf("readMatch : (*matchArr) size : %zu ", (*matchArr).size());
  Printf("readMatch : startIndex %d", startIndex);

         
  if((*matchArr).size() < 1) {
  	LOGP(info, "matchArr was 0");
    return filteredMatches;
  }       
         
         
  Printf("readMatch : (*matchArr)[startIndex].getEvent() %d eventID %d",
         (*matchArr)[startIndex].getEvent(), eventID);
         
  if ((*matchArr)[startIndex].getEvent() != eventID) {
    Printf("This shouldnt happen");
  } else
    found = true;

  Printf("readMatch : event %d startIndex %d", eventID, startIndex);

  for (int i = startIndex; i < matchArr->size(); i++) {
    const auto &track = (*matchArr)[i];
    // filteredMatches->push_back(track);
    // startIndex = i;
    //  ef : we cant use this atm, since the clusters from the same trigger
    //  sometimes have different eventNumber!
    if (track.getEvent() != eventID) {
      startIndex = i; // new startIndex for next event
      Printf("readMatch : eventID changed - %d; end of loop ",
             track.getEvent());

      break;
    } else {

      filteredMatches.push_back(track);
    }
  }

  Printf("readMatch : new startIndex %d", startIndex);

  return filteredMatches;
}

void HmpidDataReader::initializeClusterTree(/*std::vector<Cluster> *&cluArr, std::vector<Trigger> *&trigArr*/ ) {
  // TTree *treeClu = (TTree *)fileClu.Get("o2sim");
  treeClu = dynamic_cast<TTree*>(fileClu.Get("o2hmp"));
  if (!treeClu)
    treeClu = dynamic_cast<TTree*>(fileClu.Get("o2sim"));

  if (!treeClu) {
    Printf("Error accessing TTree");
  }
  
  treeClu->Print("toponly");

  treeClu->SetBranchAddress("HMPIDclusters", &cluArr);
  treeClu->SetBranchAddress("InteractionRecords", &trigArr);

  treeClu->GetEntry(0);
  // return treeClu;
} 

TTree *HmpidDataReader::initializeMCTree(std::vector<o2::MCTrack> *&mcArr) {

  treeKine = dynamic_cast<TTree*>(fileKine.Get("o2hmp"));
  if (!treeKine)
    treeKine = dynamic_cast<TTree*>(fileKine.Get("o2sim"));

  if (!treeKine) {
    Printf("Error accessing TTree");
    /*fileKine->Close();
    delete fileKine;*/
    return nullptr;
  }

  treeKine->SetBranchAddress("MCTrack", &mcArr);
  treeKine->GetEntry(0);
  treeKine->Print("toponly");
  return treeKine;
}

std::vector<o2::MCTrack> *
HmpidDataReader::readMC(std::vector<o2::MCTrack> *&mcArr, TTree *treeKine,
                        int eventID) {
  if (treeKine == nullptr) {
    Printf("Error : treeKine == nullptr");
    return nullptr;
  }

  if (eventID < treeKine->GetEntries()) {
    treeKine->GetEntry(eventID);
    Printf("readMC at entry %d", eventID);
    return mcArr;
  } else {
    Printf("eventId > treeKine->GetEntries()");
    return nullptr;
  }
}

// for the given eventID; read trackID
const o2::MCTrack *HmpidDataReader::getMCEntry(std::vector<o2::MCTrack> *mcArr,
                                               int trackID) {

  if (trackID < 0 || trackID >= mcArr->size()) {
    return nullptr;
  } else {
    // const auto& track = mcArr->at(trackID);
    for (int i = 0; i < mcArr->size(); ++i) {
      const auto &mcTrack = (*mcArr)[i];
      if (i == trackID) {
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
               mcTrack.GetStartVertexCoordinatesZ());
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
