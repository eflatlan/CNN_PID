
#include "DataFormatsHMP/Cluster.h"
#include "DataFormatsHMP/Digit.h"
#include "DataFormatsHMP/Trigger.h"
#include "HMPIDReconstruction/Clusterer.h"
#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <iostream>
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "DataFormatsHMP/Trigger.h"
#include "ReconstructionDataFormats/MatchInfoHMP.h"

// Function to process "hmpidclusters.root" with "o2hmp" tree
void ProcessHmpidClustersFile(const char* fileName, const char* treeName) {
    TFile* file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Could not open " << fileName << std::endl;
        return;
    }

    TTree* tree = (TTree*)file->Get(treeName);
    if (!tree) {
        std::cerr << "Could not find the tree " << treeName << " in " << fileName << std::endl;
        file->Close();
        delete file;
        return;
    }

    std::vector<o2::hmpid::Cluster>* hmpidClusters = nullptr;
    std::vector<o2::hmpid::Trigger>* interactionRecords = nullptr;

    tree->SetBranchAddress("HMPIDclusters", &hmpidClusters);
    tree->SetBranchAddress("InteractionRecords", &interactionRecords);

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        if (hmpidClusters) {
            std::cout << fileName << " [" << treeName << "] Entry " << i << ", HMPIDclusters size: " << hmpidClusters->size() << std::endl;
        }
        if (interactionRecords) {
            std::cout << fileName << " [" << treeName << "] Entry " << i << ", InteractionRecords size: " << interactionRecords->size() << std::endl;
        }
    }

    file->Close();
    delete file;
}

// Function to process "o2match_HMP.root" with "matchHMP" tree
void ProcessMatchHMPFile(const char* fileName, const char* treeName) {
    TFile* file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Could not open " << fileName << std::endl;
        return;
    }

    TTree* tree = (TTree*)file->Get(treeName);
    if (!tree) {
        std::cerr << "Could not find the tree " << treeName << " in " << fileName << std::endl;
        file->Close();
        delete file;
        return;
    }

    std::vector<o2::dataformats::MatchInfoHMP>* hmpMatchInfo = nullptr;
    std::vector<o2::MCCompLabel>* matchHMPMCTruth = nullptr;

    tree->SetBranchAddress("HMPMatchInfo", &hmpMatchInfo);
    tree->SetBranchAddress("MatchHMPMCTruth", &matchHMPMCTruth);

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        if (hmpMatchInfo) {
            std::cout << fileName << " [" << treeName << "] Entry " << i << ", HMPMatchInfo size: " << hmpMatchInfo->size() << std::endl;
        }
        if (matchHMPMCTruth) {
            std::cout << fileName << " [" << treeName << "] Entry " << i << ", MatchHMPMCTruth size: " << matchHMPMCTruth->size() << std::endl;
        }
    }

    file->Close();
    delete file;
}

// Master function to process all files and trees
void PrintTreeInfo() {
    ProcessHmpidClustersFile("hmpidclusters.root", "o2hmp");
    ProcessMatchHMPFile("o2match_hmp.root", "matchHMP");
}

