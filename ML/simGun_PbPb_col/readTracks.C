#if !defined(__CLING__) || defined(__ROOTCLING__)

//#include "Framework/DataProcessorSpec.h"
#include <TTree.h>
#include <TFile.h>
#include "DPLUtils/MakeRootTreeWriterSpec.h"

#include "Framework/InputSpec.h"

#include "DataFormatsHMP/Digit.h"

#include "DataFormatsHMP/Trigger.h"

#include "SimulationDataFormat/MCTruthContainer.h"

#include "SimulationDataFormat/MCCompLabel.h"
#include "ReconstructionDataFormats/MatchInfoHMP.h"
#include "ReconstructionDataFormats/TrackHMP.h"


#include "HMPIDSimulation/HMPIDDigitizer.h"
#include "DataFormatsHMP/Digit.h"
#include "DataFormatsHMP/Trigger.h"

#include "Framework/Logger.h"
#endif

void testDigit(string mOutRootFileName = "hmpiddigits.root")

{

  // ExecutionTimer mExTimer;

  std::vector<o2::hmpid::Trigger>* mTriggers;

  std::vector<o2::hmpid::Digit>* mDigits;

  // std::vector<o2::MCCompLabel>* mDigitLabels; // ef : added...

  

  //using LabelsType =

  std::vector<o2::dataformats::MCTruthContainer<o2::MCCompLabel>>* mDigitLabelsVec;

  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* mDigitLabels;

  // or o2::dataformats::MCTruthContainer<o2::MCCompLabel>??

  bool mUseMC = true;


  std::cout<<"Try to set" << std::endl;

  std::string mDigitMCTruthBranchName = "HMPIDDigitMCTruth";




  TString filename = TString::Format("%s", mOutRootFileName.c_str());

  TString tit = TString::Format("HMPID Digits File Decoding");
    std::cout<<"okset" << std::endl;
    //std::unique_ptr<TFile> mfileOut;

  std::unique_ptr<TTree> mDigitTree;
  TFile *mfileOut = TFile::Open(filename);
std::cout<<"okset" << std::endl;
  if(!mfileOut)
    std::cout<<"mfileOut nptr";
  //TTree *mDigitTree = (TTree*)mfileOut->Get("o2sim");
  
  
  mDigitTree.reset((TTree*)mfileOut->Get("o2sim"));
  if(!mDigitTree)
    std::cout<<"mDigitTree nptr";
  //tTree.reset(new TTree("o2hmp", 




// BranchDefinition<o2::dataformats::MCTruthContainer<o2::emcal::MCLabel>>{InputSpec{"emcaldigitlabels", "EMC", "DIGITSMCTR"}, "EMCALDigitMCTruth", mctruth ? 1 : 0})();





  mDigitTree->Print("toponly");

//TBranch* branch = mDigitTree->GetBranch(mDigitMCTruthBranchName.c_str());
// HMPDigit

TBranch* branch = mDigitTree->GetBranch("HMPDigit");
if (branch) {
    const char* className = branch->GetClassName();
    std::cout << "Class type of the objects in this branch is: " << className << std::endl;
    
mDigitTree->SetBranchAddress("HMPDigit", &mDigits);
} else {
    std::cout << "Branch not found!" << std::endl;
}


  if ( mDigitTree->GetBranch(mDigitMCTruthBranchName.c_str())) {
    //mDigitTree->SetBranchAddress(mDigitMCTruthBranchName.c_str(), &mDigitLabels);  
    std::cout<<"got branch-adress";
  }
  else {    std::cout<<"didnt get branch-adress";}
    

  if(!mDigitLabels)
    std::cout<<"mDigitLabels nptr";
  else
    std::cout<<"dmDigitLabels ok";
  
  std::cout<<mDigitTree->GetEntries();
  mDigitTree->GetEntry(0);
  
  
  if(!mDigitLabels)
    std::cout<<"mDigitLabels nptr";
  
  
  
  
  mDigitTree->SetBranchAddress("InteractionRecords", &mTriggers);

  mDigitTree->SetBranchAddress("HMPDigit", &mDigits);


/******************************************************************************
branch: HMPDigit                  2682

branch: InteractionRecords         582

branch: HMPIDDigitMCTruth          783
*/
 
    
  if(!mTriggers)
    std::cout<<"mTriggers nptr";


  if(!mDigits)
    std::cout<<"mDigits nptr";




  
  mTriggers->clear();

  mDigits->clear();


}

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <memory>
#include <iostream>

void readTracks() {

    

    std::unique_ptr<TFile> file(TFile::Open("o2match_hmp.root", "READ"));
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file" << std::endl;
        return;
    } else {
    
    	std::cout<<"got file ok" << std::endl;
    }

    // Get the TTree from the file
    TTree* tree = nullptr;
    file->GetObject("matchHMP", tree);
    if (!tree) {
        std::cerr << "Error reading the TTree 'matchHMP'" << std::endl;
        return;
    } else {
    	std::cout<<"got tree ok" << std::endl;
    }

    // Print basic information about the tree
    tree->Print();

    // Loop over all branches in the tree
    TObjArray* branchList = tree->GetListOfBranches();
    
    std::cout << "Number of branches : "<< branchList->GetEntries() << std::endl;
    

    for (int i = 0; i < branchList->GetEntries(); ++i) {
   
		
      TBranch* branch = static_cast<TBranch*>(branchList->At(i));
				    
			if (!branch) {
					std::cerr << "Branch at index "<< i <<"does not exist! " << std::endl;
					return;
			}
		      
		      
      if (branch) {
				const char* className = branch->GetClassName();
				std::cout << "Class type of the objects in this branch is: " << className << std::endl;
				// mDigitTree->SetBranchAddress("HMPDigit", &mDigits);
			}
		      
      std::cout << "Branch name: " << branch->GetName() << ", Entries: " << branch->GetEntries() << std::endl;
    }



    std::vector<o2::dataformats::MatchInfoHMP> *matchInfo = nullptr;
    
    if(tree->GetBranchStatus("HMPMatchInfo")) {
      tree->SetBranchAddress("HMPMatchInfo", &matchInfo);
    	std::cerr << "Branch HMPMatchInfo  found !" << std::endl;
    } else {
    	std::cerr << "Branch HMPMatchInfo not found !" << std::endl;
			return;
    }

		if(!matchInfo) {
	   	std::cerr << "matchInfo was nullptr !" << std::endl;
			return;
    } else {
	   	std::cout << "matchInfo was not nullptr " << std::endl;
  	  auto ent = tree->GetReadEntry() + 1;
  		
  		if(ent > tree->GetEntries()) {  		
		   	std::cerr << "ent " << ent <<" tree->GetReadEntry() "<<tree->GetReadEntry() <<" tree->GetEntries() "<<tree->GetEntries()<< std::endl;
  		} else {
  			std::cout << "ent " << ent <<" tree->GetReadEntry() "<<tree->GetReadEntry() <<" tree->GetEntries() "<<tree->GetEntries()<< std::endl;
			
				std::cout<< " reading at entry " <<ent<< std::endl;
				
				try {
					tree->GetEntry(ent);
					std::cout<< " read from entry " <<ent <<" ok"<< std::endl;
					
					
					if(matchInfo == nullptr) {
					 	std::cerr << "matchInfo was nullptr !" << std::endl;
						return;
					} else {
						std::cout<< " read from entry " <<ent <<" ok"<< std::endl;					
						// auto a = *matchInfo; // check if can dereference
					} 
					
					
				}  catch (const std::exception& e) {
				   // Catch and handle standard exceptions
					 std::cerr << "Standard exception caught: " << e.what() << std::endl;
				} catch (...) {
						// Catch all other exceptions
						std::cerr << "Unknown exception caught" << std::endl;
				}
			}
  		//assert(ent < tree->GetEntries()); // this should not happen
  		//tree->GetEntry(ent);    	
    	
    	
    
    }
    
    
    std::vector<o2::MCCompLabel> *mcTruthTrack = nullptr; 
    
    if(tree->GetBranchStatus("MatchHMPMCTruth")) {
      tree->SetBranchAddress("MatchHMPMCTruth", &mcTruthTrack);
    } else {
    	std::cerr << "Branch MatchHMPMCTruth not found !" << std::endl;
			return;
    }
		
		if(!matchInfo) {
	   	std::cerr << "matchInfo was nullptr !" << std::endl;
			return;
    } else {			  
    	std::cout << "length of mcTruthTrack-vector " << mcTruthTrack->size() << std::endl;    
    }	
		


		/* 
    // Example of accessing data from one branch
    // Assuming 'HMPMatchInfo.mxRa' is a branch you are interested in
    Float_t* mxRa = nullptr;
    TBranch* branchMxRa = tree->GetBranch("HMPMatchInfo.mxRa");
    if (branchMxRa) {
        branchMxRa->SetAddress(&mxRa); // Set the address of the local variable to the branch
        Int_t nEntries = tree->GetEntries(); // Get the number of entries in the tree
        for (Int_t j = 0; j < nEntries; ++j) {
            branchMxRa->GetEntry(j); // Load the data for entry j
            // Now, mxRa points to the array of data for this entry
            // Process the data as needed, for example, print the first element
            std::cout << "First element of mxRa for entry " << j << ": " << mxRa[0] << std::endl;
        }
    } */

    // The file will be automatically closed and memory will be cleaned up by the unique_ptr
}

