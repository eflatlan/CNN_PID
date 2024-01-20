#pragma once
#include "CommonDataFormat/InteractionRecord.h"
#include "CommonDataFormat/RangeReference.h"

#include "HMPIDBase/Param.h"
#include "DataFormatsHMP/Cluster.h"
#include "DataFormatsHMP/Digit.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"

#include "HmpidDataReader.cpp"
#include "CkovToolsSingle.cpp"
//#include "CkovToolsErlend.cpp"
#include "ParticleUtils2.cpp"

#include <utility>
#include <vector>

#include <TMath.h>
#include <cmath>

float calcCkovFromMass(float p, float n, int pdg);


/*
struct ShallowDigit {
  int mMotherTrackId;
  int mSourceId;
  int mEventNumber;
  int mTrackId;
  int mParticlePdg;
  uint16_t mQ = 0;
  uint8_t mCh = 0; // 0xFF indicates invalid digit
  uint8_t mPh = 0;
  uint8_t mX = 0;
  uint8_t mY = 0;
  
  uint16_t getQ() const { return mQ; }
  uint8_t getCh() const { return mCh; }
  uint8_t getPh() const { return mPh; }
  uint8_t getX() const { return mX; }
  uint8_t getY() const { return mY; }


  void setEventNumber (int eventNumber) {mEventNumber = eventNumber;}
  int getEventNumber () const  {return mEventNumber;}


  void setMotherId (int motherTrackId) {mMotherTrackId = motherTrackId;}
  int getMotherId () const  {return mMotherTrackId;}


  int getSourceId () const  {return mSourceId;}
    ShallowDigit(uint16_t q, uint8_t x, uint8_t y, Int_t trackId, Int_t particlePdg) : mQ(q), mX(x), mY(y), trackId(mTrackId),particlePdg(mParticlePdg) {}
}; */ 


void evaluateClusterTrack(std::vector<o2::hmpid::ClusterCandidate>& clusterPerChamber, const o2::dataformats::MatchInfoHMP& track, const std::vector<float>& mipCharges, int mcTrackPdg, int trackNumber, int& plotNumber);


std::array<float, 3> calcCherenkovHyp(float p, float n);

/*
struct ClusterCandidate {
   
    int mCh = 0;
    double mX = 0., mY = 0.;
    double mQ = 0;
    double mChi2 = 0;
    double mXe = 0., mYe = 0.;
    int mPDG = -1;

    // vector e.l. som holder truth? // i.e., for hver track, set MIP og trackIndex fra track
    int trackId = -1;
    bool isMip = false;

    std::vector<std::pair<int,int>>* mCandidateStatusVector = nullptr;
    // std::vector<o2::hmpid::Cluster::Topology> mTopologyVector = nullptr;

    // Constructor based on the order and types you provided
    ClusterCandidate(int ch, double x, double y, double q, double chi2, 
                     double xe, double ye, /*std::vector<Topology>* topologyVector,* /, int pdg, 
                     std::vector<std::pair<int,int>>* candidateStatusVector) 
        : mCh(ch), mX(x), mY(y), mQ(q), mChi2(chi2), mXe(xe), mYe(ye), 
          /*mTopologyVector(topologyVector),* / mPDG(pdg), mCandidateStatusVector(candidateStatusVector) {}


    //obj.ch, obj.x, obj.y, obj.q, shallowDigits, obj.chi2, obj.xE, obj.yE, candStatus


    /*
    void setDigits(const std::vector<Topology>*& topologyVector) 
    {
        if(!mTopologyVector) {
            mTopologyVector = new std::vector<Topology>;
        }
        *mTopologyVector = topologyVector;
    } * /

    void addCandidateStatus(int iTrack, int hadronCandidateBit) 
    {
        if(!mCandidateStatusVector) {
            mCandidateStatusVector = new std::vector<std::pair<int,int>>;
        }
        mCandidateStatusVector->emplace_back(iTrack, hadronCandidateBit);
    }

    std::vector<std::pair<int,int>>* getCandidateStatus()
    {    
        if(!mCandidateStatusVector) {
            mCandidateStatusVector = new std::vector<std::pair<int,int>>;
        }
        return mCandidateStatusVector;
    }
};*/



void Segmentation()
{


  int plotNumber = 0;
  ParticleUtils2 mParticleEvents; 

   auto myTree = std::make_unique<TTree>("myTree", "Tree to store clusters, trackInfo, and mcPDG");

	std::vector<o2::hmpid::ClusterCandidate>* clusterBranch = nullptr;
	std::vector<o2::hmpid::Cluster>* clusterBranch2 = nullptr;
	o2::dataformats::MatchInfoHMP* trackInfoBranch = nullptr;
	int mcPDGBranch = 0;
	
	float reconCkovBranch = 0;




	myTree->Branch("clusters2", &clusterBranch2);
	myTree->Branch("clusters", &clusterBranch);
	myTree->Branch("trackInfo", &trackInfoBranch);
	myTree->Branch("mcPDG", &mcPDGBranch);



  int cluChargeBranch, cluSizeBranch;
  		  
  		  
  		  
  float refIndexBranch, xRadBranch, yRadBranch, xMipBranch, yMipBranch;
	float thBranch, phBranch, pBranch;



	myTree->Branch("reconCkov", &reconCkovBranch);
	myTree->Branch("cluCharge", &cluChargeBranch);
	myTree->Branch("cluSize", &cluSizeBranch);
	myTree->Branch("refIndex", &refIndexBranch);
	myTree->Branch("xRad", &xRadBranch);
	myTree->Branch("yRad", &yRadBranch);
	myTree->Branch("xMip", &xMipBranch);
	myTree->Branch("yMip", &yMipBranch);
	myTree->Branch("th", &thBranch);
	myTree->Branch("ph", &phBranch);
	myTree->Branch("p", &pBranch);
		
    // clusters and triggers 
    std::vector<Cluster>* clusterArr = nullptr;
    std::vector<o2::hmpid::Topology> mTopologyFromFile, *mTopologyFromFilePtr = &mTopologyFromFile;


    std::vector<Trigger>* trigArr = nullptr;


    TTree* tCluster = HmpidDataReader::initializeClusterTree(clusterArr, trigArr, mTopologyFromFilePtr);
    // clusterArr now initialized correctly


    // MathcInfoHMP : holding trackinfo
    std::vector<o2::dataformats::MatchInfoHMP>* matchArr = nullptr;
    TTree* tMatch = HmpidDataReader::initializeMatchTree(matchArr, 0 ,0 ,0);


    // McTrack : holding PDG code of track
    std::vector<o2::MCTrack>* mcArr = nullptr;
    TTree* tMcTrack = HmpidDataReader::initializeMCTree(mcArr);


    int startIndexTrack = 0;
    if(trigArr== nullptr) {
    	Printf("HmpidDataReader::initializeClusterTree trigArr== nullptr"); 
    	return ;
  	}	
    
    
    for(int i = 0; i < trigArr->size(); i++) //for(const auto& clusters : clustersVector) // "events loop"
    { 
		
        auto pTgr = &trigArr->at(i);
        if(pTgr== nullptr) {Printf("pTgr== nullptr"); continue;}	

        const int firstEntry = pTgr->getFirstEntry();
        const int lastEntry = pTgr->getLastEntry();
        
        //Printf("Checking trigger number %d Range Clu = %d :: %d", i, firstEntry, lastEntry);

        std::vector<Cluster> oneEventClusters;

        auto fClu = static_cast<o2::hmpid::Cluster>(clusterArr->at(firstEntry));
        auto s1Clu = static_cast<o2::hmpid::Cluster>(clusterArr->at(lastEntry-1));
        auto sClu = static_cast<o2::hmpid::Cluster>(clusterArr->at(lastEntry));
        int eventNumber1 = fClu.getEventNumber();
        int eventNumberLast = sClu.getEventNumber();
        int eventNumberLast1 = s1Clu.getEventNumber();




        /* ef : TODO some clusters in the same trigger have different 
                event-number!?
        if(eventNumberLast != eventNumber1) {
            Printf("eventNumberLast%d != eventNumber1%d", eventNumberLast, eventNumber1);
            Printf("eventNumberLast1%d",eventNumberLast1);
        } // TODO: throw error? ef: 
        */


        for(int j = pTgr->getFirstEntry(); j <= pTgr->getLastEntry(); j++) {      
            const auto& clu = static_cast<o2::hmpid::Cluster>(clusterArr->at(j));
            //std::cout << j << " evNum " <<  clu.getEventNumber() << " |";
        }


        for(int j = pTgr->getFirstEntry(); j <= pTgr->getLastEntry(); j++) {      
		      const auto& clu = static_cast<o2::hmpid::Cluster>(clusterArr->at(j));


					oneEventClusters.push_back(clu);
					
					
					// ef: dont check this atm, resolve error wrt different eNum for 
					// same trigger later
					/*
          Printf("============\n Cluster Loop \n ===============");
		      if(clu.getEventNumber() != eventNumber1) {
		      	Printf("Eventnumber changed??");
		      	Printf("clu.getEventNumber()%d", clu.getEventNumber());
		      } else {
			      oneEventClusters.push_back(clu);
		           std::cout << " clu " << j << " eventNum " << clu.getEventNumber();
		      } */
        }  
          Printf("============\n Cluster Loop \n ===============");
        // find entries in tracksOneEvent which corresponds to correct eventNumber
        Printf("Reading match vector<MatchInfoHMP> for startIndexTrack %d", startIndexTrack);


				// read instead MatchInfoHMP field mEvent? (i)
        std::vector<o2::dataformats::MatchInfoHMP>* tracksOneEvent = HmpidDataReader::readMatch(tMatch, matchArr, i, startIndexTrack);

        // get MC tracks for given event from mc;
        Printf("Reading vector<o2::MCTrack>* mcTracks for  eventNumber %d", eventNumber1);
        std::vector<o2::MCTrack>* mcTracks = HmpidDataReader::readMC(mcArr, tMcTrack, eventNumber1);


        Printf("tracksOneEvent size %d", tracksOneEvent->size());

        // for this event 





        Printf("Sorting events by chamber");

        std::sort((*tracksOneEvent).begin(), (*tracksOneEvent).end(), [](const o2::dataformats::MatchInfoHMP &a, const o2::dataformats::MatchInfoHMP &b) {
            return a.getChamber() < b.getChamber();
        });

        Printf("Sorting clusteres by chamber");

        std::sort((oneEventClusters).begin(), (oneEventClusters).end(), [](const Cluster &a, const Cluster &b) {
            return a.ch() < b.ch();
        });


        
        std::vector<o2::dataformats::MatchInfoHMP> sortedTracks[7];
        // Assign MLinfoHMP objects to corresponding vectors based on iCh value
        for (const auto &obj : *tracksOneEvent) {

	    			const auto& iCh = obj.getChamber();
            if (iCh < 0 || iCh > 6) {
                std::cerr << "Warning: iCh value out of expected range: " << iCh << std::endl;
            }            
            else {
							// check if there was a match:
							if (obj.getMatchStatus()) {
                sortedTracks[iCh].push_back(obj);
              } else {
                //Printf("track didnt match MIP skipping");
              }
                //sstd::cerr << "sortedTracks[iCh] " << iCh << " pushback" << std::endl;
            }
        }
				std::cout << "Length of sortedTracks vector ";
        for (int i = 0; i < 7; i++) {
           if (sortedTracks[i].size() != 0)
             std::cout <<  i << " : " << sortedTracks[i].size() << " |";
        }
        std::cout << std::endl;

        // Assuming the range of iCh values is from 0 to 6 (inclusive)
        std::vector<o2::hmpid::ClusterCandidate> sortedClusters[7];
        // Assign MLinfoHMP objects to corresponding vectors based on iCh value
        for (const auto &obj : oneEventClusters) {


					// from each track; we assign a label to each of the clusters in 
					// corresponding to the type of hadron it can be 
	    		const auto& iCh = obj.ch();
          if (iCh >= 0 && iCh <= 6) {
						if(sortedTracks[iCh].size() > 0) {
						  				    
				      //const std::vector<o2::hmpid::Cluster::Topology>& topology = obj.getClusterTopology();  // some info about digits associated w cluster*/


				      //std::vector<std::pair<int,int>> candStatus;
				      //candStatus.resize(sortedTracks[iCh].size()); // should now be initialized to (0,0) x numTracks
							// Printf("ClusterCandidate Ch %d", iCh);
							o2::hmpid::ClusterCandidate temp(obj.ch(), obj.x(), obj.y(), obj.q(), obj.chi2(), obj.xe(), obj.ye(), obj.getPDG(), obj.size());
							sortedClusters[iCh].emplace_back(temp);
	          }  else {
            	//std::cerr << "sortedTracks[iCh] " << iCh << " empty " << sortedTracks[iCh].size() << std::endl;
            }
             
          } else {
            std::cerr << "Warning: iCh value out of expected range: " << iCh << std::endl;
          }
          
        } // end for
        

        for(int i = 0; i < 7; i++) {
            
            // check if has more than one track --> this means there is no candidates
            if(sortedTracks[i].size() < 1) {
            	  //Printf("sortedTracks[iCh%d].size() %d", i, sortedTracks[i].size());
                continue;
            }

            auto& clusterPerChamber = sortedClusters[i];


            std::vector<float> mipCharges;
            // fill charges of MIPs
            for(const auto& track : sortedTracks[i]) {

                auto q = track.getMipClusQ();
                mipCharges.emplace_back(q);
            }


            int tNum = 0;
            for(const auto& track : sortedTracks[i]) {
            	  //Printf("TrackNumber%d track[iCh%d].size() %d", tNum++, i, sortedTracks[i].size());
            	  
            	  //Printf("TrackNumberMom %f", track.getHmpMom());
            	  
                // pass clusters (and track) by reference, and add its status per track (adding to candStatus vector )

                // for each clusterPerChamber we will have a "candidate-status" for each of the photons, this is a vector of length of sortedTracks[i].size();
                // and holds the fields 

                

								// checked earleir also but..
								if (track.getMatchStatus()) { 
								//if (true) { 
		              // get MCTrack correspondign to trackId
		              const auto mcTrackIndex = track.getTrackIndex();

		              // find the PDG code in teh o2Kine_sim.root file by matching the mcTrackIndex for the current event ; 
		              //--((const o2::MCTrack* mcTrack = HmpidDataReader::getMCEntry(mcTracks, mcTrackIndex);

		              //const int mcTrackPdg = mcTrack->GetPdgCode();
		              

									//mcPDGBranch = mcTrack->GetPdgCode(); 
	


								  const double winThick = 0.5, radThick = 1.5; const int gapThick = 8;
  								const double gapIdx = 1.0005, winIdx = 1.5787;  
  								auto dz = 9.25;
  								
  
		              float xRad,  yRad,  xPc,  yPc,  th,  ph; // make for these 
		              float xMip = track.getMipX(), yMip = track.getMipY(); // and thse 
		              track.getHMPIDtrk(xRad,  yRad,  xPc,  yPc,  th,  ph);
		              
		              auto ypc = yRad + dz * TMath::Tan(th)*  TMath::Sin(ph);		              
		              auto xpc = xRad +  dz * TMath::Tan(th)*  TMath::Cos(ph);		              
		              
		              TVector2 mip(xMip, yMip);
		              
		              TVector2 pcC(xpc, ypc);		              
		              
									if((mip-pcC).Mod() > 1) {
								  //	Printf("(mip-pcC).Mod() > 1 %.2f", (mip-pcC).Mod());
										continue;
									}

		              refIndexBranch = track.getRefIndex(); 
		              cluChargeBranch = track.getMipClusQ();
		              cluSizeBranch = track.getMipClusSize();
		              		              
									xRadBranch = xRad; // also for xPc yPc 
									yRadBranch = yRad; // also for xPc yPc 
									
									xMipBranch = xMip; // also for xPc yPc 
									yMipBranch = yMip; // also for xPc yPc


									thBranch = th; // also for xPc yPc 
									phBranch = ph; // also for xPc yPc
					
					
									// get the calculated recon track ckov signal
									reconCkovBranch = track.getHMPsignal();
									
									
									// clusterPerChamber : vector of ClusterCandidate-objects
									// track : object with 10 scalar values
									// mcTradckPDG : MC truth
									// 
									int pdg = track.getMipClusEventPDG();
                  evaluateClusterTrack(clusterPerChamber, track, mipCharges, pdg, tNum, plotNumber);
                  
                  
                  
                  // saving...
                  // branch : mcTrackPdg
                  // branch : clusterPerChamber --> vector of clusterCandidates
                  // branch : track --> trackInformation (Ra, MIP, refIndex, momentum, theta, phi)





									clusterBranch = const_cast<std::vector<o2::hmpid::ClusterCandidate>*>(&clusterPerChamber); // Make sure your ClusterCandidate class is compatible with ROOT I/O
									trackInfoBranch = const_cast<o2::dataformats::MatchInfoHMP*>(&track); // Make sure your MatchInfoHMP class is compatible with ROOT I/O
									
									//mcPDGBranch = mcTrack->GetPdgCode(); 
									
									

									pBranch = track.getHmpMom();  									
									// Fill the tree
									myTree->Fill();
									
									

									mParticleEvents.fillCandidate(clusterPerChamber, track, 0);
                  
                  // for(auto& clusterPerChamber)
                }
								else  
								{
									//Printf("Track didnt match!");
  							}
            }

            // save ClusterPerChamber

            // save trackInfo : 
            //                  trackIndex : mIdxTrack
            //                  scalar : p, refIndex,
            //                  trackInfo theta, phi, xRad, yRad; {xPc yPc is redundant}
            //                  MIP            : x, y, q
            //                  MCTruth : pdg code of track || pdg code of MIP?
            



            // now assigned all types of candidates for all clusters in teh given chamber
            // match sortedClusters[i] with sortedTracks[i] --> 
            // 

        }
        
    }
    
    auto fOut = std::make_unique<TFile>("MLOUTPUT.root", "RECREATE");
		myTree->Write();
		fOut->Close();
		
		
		Long64_t nEntries = myTree->GetEntries();
		for (Long64_t i = 0; i < nEntries; ++i) {
				myTree->GetEntry(i);

				// The TBranches you are interested in (clusterBranch, trackInfoBranch, mcPDGBranch)
				// should now be filled with the data from the i-th entry.

				// Perform some checks or analysis using these variables.
				// For example, to print the size of clusterBranch (assuming it is a std::vector)
				//std::cout << "Size of clusterBranch: " << clusterBranch->size() << std::endl;
				

		}
		

	   mParticleEvents.writeH5();
}







void read_tree() {
    TCanvas *canvas1 = new TCanvas("canvas1", "Particle Distributions", 1200, 400);

    // Initialize histograms
    TH1F *momHistPion = new TH1F("momHistPion", ";Cherenkov Angle;Frequency", 50, 0.5, 0.728);
    TH1F *momHistKaon = new TH1F("momHistKaon", ";Cherenkov Angle;Frequency", 50, 0.5, 0.728);
    TH1F *momHistProton = new TH1F("momHistProton", ";Cherenkov Angle;Frequency", 50, 0.5, 0.728);

    // Reading for pions
    auto fPion = std::make_unique<TFile>(Form("MLOUTPUTPion.root"));
    std::unique_ptr<TTree> treePion((TTree*)fPion->Get("myTree"));

    // Reading for kaons
    auto fKaon = std::make_unique<TFile>(Form("MLOUTPUTKaon.root"));
    std::unique_ptr<TTree> treeKaon((TTree*)fKaon->Get("myTree"));

    // Reading for protons
    auto fProton = std::make_unique<TFile>(Form("MLOUTPUTProton.root"));
    std::unique_ptr<TTree> treeProton((TTree*)fProton->Get("myTree"));

    // Variables for reading tree branches
    float pPion, pKaon, pProton;
    float ckovPion, ckovKaon, ckovProton;

    // Setup branches
    treePion->SetBranchAddress("p", &pPion);
    treePion->SetBranchAddress("reconCkov", &ckovPion);

    treeKaon->SetBranchAddress("p", &pKaon);
    treeKaon->SetBranchAddress("reconCkov", &ckovKaon);

    treeProton->SetBranchAddress("p", &pProton);
    treeProton->SetBranchAddress("reconCkov", &ckovProton);

    // Assume that all trees have the same number of entries
    Long64_t nEntries = treePion->GetEntries();

    for (Long64_t i = 0; i < nEntries; ++i) {
        // Fill histograms based on momentum range
        treePion->GetEntry(i);
        if (pPion > 2.8 && pPion < 2.85 && ckovPion > 0.66 && 0.66 + 0.048 > ckovPion) {

            momHistPion->Fill(ckovPion);
        }

        treeKaon->GetEntry(i);
        if (pKaon > 2.8 && pKaon < 2.85 && ckovKaon > 0.644 && 0.69 > ckovKaon) {
            momHistKaon->Fill(ckovKaon);
        }

        treeProton->GetEntry(i);
        if (pProton > 2.8 && pProton < 2.85 && ckovProton > 0.59 && 0.645 > ckovProton) {
            momHistProton->Fill(ckovProton);
        }
    }

    // Draw and fit histograms
    canvas1->Divide(3);

    // Pion
    canvas1->cd(1);
    momHistPion->Draw();
    TF1 *fitPion = new TF1("fitPion", "gaus", 0.4, 0.8);
    momHistPion->Fit("fitPion", "R");
    TLatex lPion;
    //lPion.DrawLatex(0.45, momHistPion->GetMaximum() * 0.9, Form("#entries = %lld, #mu = %.3f, #sigma = %.3f", momHistPion->GetEntries(), fitPion->GetParameter(1), fitPion->GetParameter(2)));

    // Kaon
    canvas1->cd(2);
    momHistKaon->Draw();
    TF1 *fitKaon = new TF1("fitKaon", "gaus", 0.4, 0.8);
    momHistKaon->Fit("fitKaon", "R");
    TLatex lKaon;
    //lKaon.DrawLatex(0.45, momHistKaon->GetMaximum() * 0.9, Form("#entries = %lld, #mu = %.3f, #sigma = %.3f", momHistKaon->GetEntries(), fitKaon->GetParameter(1), fitKaon->GetParameter(2)));

    // Proton
    canvas1->cd(3);
    momHistProton->Draw();
    TF1 *fitProton = new TF1("fitProton", "gaus", 0.4, 0.8);
    momHistProton->Fit("fitProton", "R");
    TLatex lProton;
    //lProton.DrawLatex(0.45, momHistProton->GetMaximum() * 0.9, Form("#entries = %lld, #mu = %.3f, #sigma = %.3f", momHistProton->GetEntries(), fitProton->GetParameter(1), fitProton->GetParameter(2)));

    // Save plots
    canvas1->Update();
    canvas1->SaveAs("ParticlePlots.png");
    canvas1->SaveAs("ParticlePlots.jpg");
    canvas1->SaveAs("ParticlePlots.pdf");
    canvas1->SaveAs("ParticlePlots.eps");
}



float calcCkovFromMass(float p, float n, int pdg)
{
    // Constants for particle masses (in GeV/c^2)
    const float mass_Pion = 0.1396, mass_Kaon = 0.4937, mass_Proton = 0.938;

    float m; // variable to hold the mass

    // Switch based on absolute value of PDG code
    switch (std::abs(pdg))
    {
        case 211:
            m = mass_Pion;
            break;
        case 321:
            m = mass_Kaon;
            break;
        case 2212:
            m = mass_Proton;
            break;
        default:
            return 0;  // return 0 if PDG code doesn't match any known codes
    }

    const float p_sq = p * p;
    const float refIndexFreon = n;  // Assuming n is the refractive index
    const float cos_ckov_denom = p * refIndexFreon;

    // Sanity check
    if (p_sq + m * m < 0) {
        return 0;
    }

    const auto cos_ckov = static_cast<float>(TMath::Sqrt(p_sq + m * m) / cos_ckov_denom);

    // Sanity check
    if (cos_ckov > 1 || cos_ckov < -1) {
        return 0;
    }

    const auto ckovAngle = static_cast<float>(TMath::ACos(cos_ckov));

    return ckovAngle;
}




void evaluateClusterTrack(std::vector<o2::hmpid::ClusterCandidate>& clusterPerChamber, const o2::dataformats::MatchInfoHMP& track, const std::vector<float>& mipCharges, int mcTrackPdg, int trackNumber, int& plotNumber)
{

        const auto eventCnt = track.getEvent(); // check it corresponds to entry in loop of events?
        const auto momentum = track.getHmpMom();

        const auto nF = 1.2929 - 0.0025;//track.getRefIndex(); // ef: aon it only gets the mean value ; TODO: in MatchHMP.cxx get calibration value
        // https://github.com/eflatlan/AliceO2/blob/811fcc6d00b363b1e96e0aa8269d46eed95d879b/Detectors/GlobalTracking/src/MatchHMP.cxx#L432C28-L432C38
        // https://github.com/AliceO2Group/AliceO2/blob/e6b603e4c92f98733ff9f7954100140e72bd99f6/Detectors/HMPID/base/include/HMPIDBase/Param.h#L159C37-L159C64
          



        const auto nQ =1.583;
        const auto nG = 1.0005;

        // make static?


				// make account for varying nF ? +- 2 std ? 
        const auto& ckovHypsMin = calcCherenkovHyp(momentum, nF); 
						
		
        const auto& ckovHypsMax = calcCherenkovHyp(momentum, nF); 

        float xRad, yRad, xPc, yPc, thetaP, phiP;
        track.getHMPIDtrk(xRad, yRad, xPc, yPc, thetaP, phiP);

        float xMip, yMip;
        int nph, q;
        track.getHMPIDmip(xMip, yMip, q, nph);

        const auto& L = 0.5; // 
        double radParams[7] = {xRad, yRad, L, thetaP, phiP, momentum, static_cast<double>(mcTrackPdg*1.0)}; // ef : TODO add PID to MLinfoHMP?
        
        double refIndexes[3] = {nF, nQ, nG};


        double MIP[3] = {xMip, yMip, static_cast<double>(q*1.0)};
				// (double radParams[7], double refIndexes[3], double MIP[3],
        //  std::array<float, 3> ckovHyps, float trackCkov, int eventCnt)
        // ef: TODO use MIP to get radius and phi in CkovTools: 
				auto ckovAngle = 0.;
	
	
				
        CkovTools ckovTools(radParams, refIndexes, MIP, ckovHypsMin, ckovHypsMax, ckovAngle, eventCnt, mcTrackPdg);

        Printf(" Event%d Track%d  : ckovHyps = <%.3f, %.3f> | <%.3f, %.3f> | <%.3f, %.3f>", eventCnt, ckovTools.getMinCkovPion(),ckovTools.getMaxCkovPion(),ckovTools.getMinCkovKaon(), ckovTools.getMaxCkovKaon(),ckovTools.getMinCkovProton(), ckovTools.getMaxCkovProton(), eventCnt); 

       
        // only consider if adequate momentum? 
        //LOGP(info, "Momentum {} ckovHyps {} {} {}", ckov[0], ckov[1], ckov[2]);
        
        if(TMath::IsNaN(ckovHypsMax[0])) { 
        	Printf("was isnan!!!");
        	return;
        }
        //numBackgroundPhotons, numFoundActualCkov, numActualCkov, numBackgroundLabeledCkov
        std::array<int, 4> arrayInfo;



        // add charge of MIPS here? to skip them in candidates ? 

        // clusterPerChamber by reference, add element to vector {trackNumber, bitHadronStatus}


        // mcTrackPdg check that it matches with clusterPDG?
        // 
        ckovTools.segment(clusterPerChamber, arrayInfo, track.getTrackIndex(), mipCharges, xMip, yMip, q /*MIP-charge*/, mcTrackPdg, track, trackNumber, plotNumber); // temp --> mapBins
}


const float mass_Pion = 0.1396, mass_Kaon = 0.4937, mass_Proton = 0.938; // masses in
std::array<float, 3> masses = {mass_Pion, mass_Kaon, mass_Proton};
const float mass_Pion_sq = mass_Pion*mass_Pion, mass_Kaon_sq = mass_Kaon*mass_Kaon,  mass_Proton_sq = mass_Proton*mass_Proton;
std::array<float, 3> calcCherenkovHyp(float p, float n)
{
  const float p_sq = p*p;
  const float cos_ckov_denom = p*n;
  const auto cos_ckov_Pion = static_cast<float>(TMath::Sqrt(p_sq + mass_Pion_sq)/(cos_ckov_denom)); // n = refIndexFreon 1.289 later make it random?

  const float cos_ckov_Kaon = static_cast<float>(TMath::Sqrt(p_sq + mass_Kaon_sq)/(cos_ckov_denom)); 
  const float cos_ckov_Proton = static_cast<float>(TMath::Sqrt(p_sq + mass_Proton_sq)/(cos_ckov_denom));
  

  const float ckovAnglePion = static_cast<float>( TMath::ACos(cos_ckov_Pion)); 
  const float ckovAngleKaon = static_cast<float>(TMath::ACos(cos_ckov_Kaon)); 
  const float ckovAngleProton = static_cast<float>(TMath::ACos(cos_ckov_Proton)); 

  Printf("Pion %.3f Kaon %.3f Proton %.3f", ckovAnglePion, ckovAngleKaon, ckovAngleProton);

  return {ckovAnglePion, ckovAngleKaon, ckovAngleProton};
}
