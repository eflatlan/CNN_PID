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
	
	TBranch* branch2 = mDigitTree->GetBranch(mDigitMCTruthBranchName.c_str());
	if (branch2) {
	  const char* className = branch2->GetClassName();
	  std::cout << "Class type of the objects in this branch is: " << className << std::endl;
		  
		mDigitTree->SetBranchAddress(mDigitMCTruthBranchName.c_str(), &mDigitLabels);
	} else {
		  std::cout << "Branch not found!" << std::endl;
	}



	
	/*
  if (mDigitTree->GetBranch(mDigitMCTruthBranchName.c_str())) {
    mDigitTree->SetBranchAddress(mDigitMCTruthBranchName.c_str(), &mDigitLabels);  
    std::cout<<"got branch-adress" << std::endl;
  }
  else {    
  	std::cout<<"didnt get branch-adress"<< std::endl;
	} */
    

  if(!mDigitLabels)
    std::cout<<"mDigitLabels nptr"<< std::endl;
  else { 
   	std::cout << "matchInfo was not nullptr " << std::endl;
	  auto ent = mDigitTree->GetReadEntry() + 1;
				
				
		std::cout<< "entries " << mDigitTree->GetEntries();
		
		if(ent > mDigitTree->GetEntries()) {  		
		 	std::cerr << "ent " << ent <<" tree->GetReadEntry() "<<mDigitTree->GetReadEntry() <<" tree->GetEntries() "<<mDigitTree->GetEntries()<< std::endl;
		} else {
			std::cout << "ent " << ent <<" tree->GetReadEntry() "<<mDigitTree->GetReadEntry() <<" tree->GetEntries() "<<mDigitTree->GetEntries()<< std::endl;
		
			std::cout<< " reading at entry " <<ent<< std::endl;
			
			try {
				// mDigitTree->GetEntry(ent);
				std::cout<< " read from entry " << ent <<" ok"<< std::endl;				
				
				if(mDigitLabels == nullptr) {
				 	std::cerr << "mDigitLabels was nullptr !" << std::endl;
					return;
				} else {
				
				 	std::cout << "mDigitLabels was not nullptr " << std::endl;
					//o2::dataformats::MCTruthContainer<o2::MCCompLabel>* a = mDigitLabels; // check if can dereference
				} 
				
				
			}  catch (const std::exception& e) {
				 // Catch and handle standard exceptions
				 std::cerr << "Standard exception caught: " << e.what() << std::endl;
			} catch (...) {
					// Catch all other exceptions
					std::cerr << "Unknown exception caught" << std::endl;
			}
		}
	}

  //mDigitTree->GetEntry(0);
  
    /*  
  if(!mDigitLabels)
    std::cout<<"mDigitLabels nptr";
  
  
  
  
  mDigitTree->SetBranchAddress("InteractionRecords", &mTriggers);

  mDigitTree->SetBranchAddress("HMPDigit", &mDigits);


/******************************************************************************
branch: HMPDigit                  2682

branch: InteractionRecords         582

branch: HMPIDDigitMCTruth          783
* /
 
    
  if(!mTriggers)
    std::cout<<"mTriggers nptr";


  if(!mDigits)
    std::cout<<"mDigits nptr";




  
  mTriggers->clear();

  mDigits->clear(); */ 


}
