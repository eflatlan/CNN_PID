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
public:
  HmpidDataReader() {
    // Constructor to initialize any member variables or resources
  }

  ~HmpidDataReader() {
    // Destructor to clean up any allocated resources
  }

  static TTree *
  initializeMatchTree(std::vector<o2::dataformats::MatchInfoHMP> *&matchArr,
                      int eventID, int trackID, int pdg);

  static std::vector<o2::dataformats::MatchInfoHMP> *
  readMatch(TTree *tMatch, std::vector<o2::dataformats::MatchInfoHMP> *matchArr,
            int eventID, int &startIndex);
  static TTree *initializeClusterTree(
      std::vector<Cluster> *&cluArr, std::vector<Trigger> *&trigArr,
      std::vector<o2::hmpid::Topology> *&mTopologyFromFilePtr);

  static TTree *initializeMCTree(std::vector<o2::MCTrack> *&mcArr);
  static std::vector<o2::MCTrack> *readMC(std::vector<o2::MCTrack> *&mcArr,
                                          TTree *treeKine, int eventId);
  static const o2::MCTrack *getMCEntry(std::vector<o2::MCTrack> *mcArr,
                                       int trackID);
  static void readTreeEntries();
};

// should they be static here ? :
TTree *HmpidDataReader::initializeMatchTree(
    std::vector<o2::dataformats::MatchInfoHMP> *&matchArr, int eventID,
    int trackID, int pdg) {
  TFile *fMatch = new TFile("o2match_hmp.root");
  TTree *tMatch = (TTree *)fMatch->Get("matchHMP");

  if (!tMatch)
    tMatch = (TTree *)fMatch->Get("o2hmp");
  if (!tMatch)
    tMatch = (TTree *)fMatch->Get("o2sim");

  // std::vector<o2::dataformats::MatchInfoHMP>* matchArr = nullptr;
  tMatch->SetBranchAddress("HMPMatchInfo", &matchArr);
  tMatch->GetEntry(0);
  tMatch->Print("toponly");
  if (matchArr == nullptr) {
    Printf("HmpidDataReader::initializeClusterTree matchArr== nullptr");
    return nullptr;
  }

  return tMatch;
}

// eventId = eventID to be searched for;
// startIndex : index of where matchArr is to be searched
// newStartIndex startIndex for next event
std::vector<o2::dataformats::MatchInfoHMP> *
HmpidDataReader::readMatch(TTree *tMatch,
                           std::vector<o2::dataformats::MatchInfoHMP> *matchArr,
                           int eventID, int &startIndex) {
  if (!tMatch) {
    Printf("TTree not initialized");
    return nullptr;
  }

  // Prepare to store filtered matches
  std::vector<o2::dataformats::MatchInfoHMP> *filteredMatches =
      new std::vector<o2::dataformats::MatchInfoHMP>;
  // tracks should be stored in "time" --> when we find our event we can then
  // switch this condition "of" when the event changes:

  bool found = false;

  if (matchArr == nullptr) {
    Printf("HmpidDataReader::readMatch :: matchArr== nullptr");
    return nullptr;
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

      filteredMatches->push_back(track);
    }
  }

  Printf("readMatch : new startIndex %d", startIndex);

  return filteredMatches;
}

TTree *HmpidDataReader::initializeClusterTree(
    std::vector<Cluster> *&cluArr, std::vector<Trigger> *&trigArr,
    std::vector<o2::hmpid::Topology> *&mTopologyFromFilePtr) {
  TFile *fileClu = TFile::Open("hmpidclusters.root");
  if (!fileClu || fileClu->IsZombie()) {
    Printf("Error opening file");
    return nullptr;
  }

  TTree *treeClu = (TTree *)fileClu->Get("o2sim");
  if (!treeClu)
    treeClu = (TTree *)fileClu->Get("o2hmp");
  if (!treeClu) {
    Printf("Error accessing TTree");
    fileClu->Close();
    delete fileClu;
    return nullptr;
  }
  treeClu->Print("toponly");

  treeClu->SetBranchAddress("HMPIDDigitTopology", &mTopologyFromFilePtr);
  treeClu->SetBranchAddress("HMPIDclusters", &cluArr);
  treeClu->SetBranchAddress("InteractionRecords", &trigArr);

  treeClu->GetEntry(0);
  return treeClu;
}

TTree *HmpidDataReader::initializeMCTree(std::vector<o2::MCTrack> *&mcArr) {
  TFile *fileKine = TFile::Open("o2sim_Kine.root");
  if (!fileKine || fileKine->IsZombie()) {
    Printf("Error opening file");
    return nullptr;
  }

  TTree *treeKine = (TTree *)fileKine->Get("o2sim");
  if (!treeKine) {
    Printf("Error accessing TTree");
    fileKine->Close();
    delete fileKine;
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
