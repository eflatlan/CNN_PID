#pragma once
#include "CommonDataFormat/InteractionRecord.h"
#include "CommonDataFormat/RangeReference.h"

#include "HMPIDBase/Param.h"
#include "DataFormatsHMP/Cluster.h"
#include "DataFormatsHMP/Digit.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"

#include "HmpidDataReader.cpp"
#include "CkovTools.cpp"

#include <utility>
#include <vector>
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


void evaluateClusterTrack(std::vector<ClusterCandidate>& clusterPerChamber, const o2::dataformats::MatchInfoHMP& track, const std::vector<float>& mipCharges, int mcTrackPdg, int trackNumber);


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



void process()
{

   auto myTree = std::make_unique<TTree>("myTree", "Tree to store clusters, trackInfo, and mcPDG");

	std::vector<ClusterCandidate>* clusterBranch = nullptr;
	o2::dataformats::MatchInfoHMP* trackInfoBranch = nullptr;
	int mcPDGBranch = 0;

	myTree->Branch("clusters", &clusterBranch);
	myTree->Branch("trackInfo", &trackInfoBranch);
	myTree->Branch("mcPDG", &mcPDGBranch);



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
    if(trigArr== nullptr) {Printf("HmpidDataReader::initializeClusterTree trigArr== nullptr"); return ;}	
    
    
    for(int i = 0; i < trigArr->size(); i++) //for(const auto& clusters : clustersVector) // "events loop"
    { 
		
        auto pTgr = &trigArr->at(i);
        if(pTgr== nullptr) {Printf("pTgr== nullptr"); continue;}	

        const int firstEntry = pTgr->getFirstEntry();
        const int lastEntry = pTgr->getLastEntry();
        
        Printf("Checking trigger number %d Range Clu = %d :: %d", i, firstEntry, lastEntry);

        std::vector<Cluster> oneEventClusters;

        auto fClu = static_cast<o2::hmpid::Cluster>(clusterArr->at(firstEntry));
        auto s1Clu = static_cast<o2::hmpid::Cluster>(clusterArr->at(lastEntry-1));
        auto sClu = static_cast<o2::hmpid::Cluster>(clusterArr->at(lastEntry));
        int eventNumber1 = fClu.getEventNumber();
        int eventNumberLast = sClu.getEventNumber();
        int eventNumberLast1 = s1Clu.getEventNumber();


        if(eventNumberLast != eventNumber1) {
            Printf("eventNumberLast%d != eventNumber1%d", eventNumberLast, eventNumber1);
            Printf("eventNumberLast1%d",eventNumberLast1);
        } // TODO: throw error? ef:


        for(int j = pTgr->getFirstEntry(); j <= pTgr->getLastEntry(); j++) {      
            const auto& clu = static_cast<o2::hmpid::Cluster>(clusterArr->at(j));
            std::cout << j << " evNum " <<  clu.getEventNumber() << " |";
        }


        for(int j = pTgr->getFirstEntry(); j <= pTgr->getLastEntry(); j++) {      
		      const auto& clu = static_cast<o2::hmpid::Cluster>(clusterArr->at(j));

          Printf("============\n Cluster Loop \n ===============");
		      if(clu.getEventNumber() != eventNumber1) {
		      	Printf("Eventnumber changed??");
		      	Printf("clu.getEventNumber()%d", clu.getEventNumber());
		      } else {
		          oneEventClusters.push_back(clu); std::cout << " clu " << j << " eventNum " << clu.getEventNumber();
		      }
        }  
          Printf("============\n Cluster Loop \n ===============");
        // find entries in tracksOneEvent which corresponds to correct eventNumber
        Printf("Reading match vector<MatchInfoHMP> for startIndexTrack %d", startIndexTrack);

        std::vector<o2::dataformats::MatchInfoHMP>* tracksOneEvent = HmpidDataReader::readMatch(tMatch, matchArr, eventNumber1, startIndexTrack);

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
                Printf("track didnt match MIP skipping");
              }
                //sstd::cerr << "sortedTracks[iCh] " << iCh << " pushback" << std::endl;
            }
        }

        for (int i = 0; i < 7; i++) {
            std::cout << "Length of sortedTracks vector " << i << ": " << sortedTracks[i].size() << 			std::endl;
        }
        

        // Assuming the range of iCh values is from 0 to 6 (inclusive)
        std::vector<ClusterCandidate> sortedClusters[7];
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
							ClusterCandidate temp(obj.ch(), obj.x(), obj.y(), obj.q(), obj.chi2(), obj.xe(), obj.ye(), obj.getPDG());
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
            	  Printf("sortedTracks[iCh%d].size() %d", i, sortedTracks[i].size());
                continue;
            }

            auto& clusterPerChamber = sortedClusters[i];


            std::vector<float> mipCharges;
            // fill charges of MIPs
            for(const auto& track : sortedTracks[i]) {
                float xMip, yMip; int q, nph;
                track.getHMPIDmip(xMip, yMip, q, nph);
                mipCharges.emplace_back(q);
            }


            int tNum = 0;
            for(const auto& track : sortedTracks[i]) {
            	  Printf("TrackNumber%d track[iCh%d].size() %d", tNum++, i, sortedTracks[i].size());
                // pass clusters (and track) by reference, and add its status per track (adding to candStatus vector )

                // for each clusterPerChamber we will have a "candidate-status" for each of the photons, this is a vector of length of sortedTracks[i].size();
                // and holds the fields 

                

								// checked earleir also but..
								if (track.getMatchStatus()) { 
									
		              // get MCTrack correspondign to trackId
		              const auto mcTrackIndex = track.getTrackIndex();

		              // find the PDG code in teh o2Kine_sim.root file by matching the mcTrackIndex for the current event ; 
		              const o2::MCTrack* mcTrack = HmpidDataReader::getMCEntry(mcTracks, mcTrackIndex);

		              const int mcTrackPdg = mcTrack->GetPdgCode();
		              
		              const int momentum = mcTrack->GetPdgCode();


		              float xRad,  yRad,  xPc,  yPc,  th,  ph;
		              float xMip = track.getMipX(), yMip = track.getMipY();
		              track.getHMPIDtrk(xRad,  yRad,  xPc,  yPc,  th,  ph);


									// clusterPerChamber : vector of ClusterCandidate-objects
									// track : object with 10 scalar values
									// mcTradckPDG : MC truth
									// 
                  evaluateClusterTrack(clusterPerChamber, track, mipCharges, mcTrackPdg, tNum);
                  
                  
                  
                  // saving...
                  // branch : mcTrackPdg
                  // branch : clusterPerChamber --> vector of clusterCandidates
                  // branch : track --> trackInformation (Ra, MIP, refIndex, momentum, theta, phi)





									clusterBranch = const_cast<std::vector<ClusterCandidate>*>(&clusterPerChamber); // Make sure your ClusterCandidate class is compatible with ROOT I/O
									trackInfoBranch = const_cast<o2::dataformats::MatchInfoHMP*>(&track); // Make sure your MatchInfoHMP class is compatible with ROOT I/O
									mcPDGBranch = mcTrack->GetPdgCode(); // This assumes mcTrack is properly initialized
									
									// Fill the tree
									myTree->Fill();
                  
                  // for(auto& clusterPerChamber)
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
				std::cout << "Size of clusterBranch: " << clusterBranch->size() << std::endl;
				
				std::cout << "Size of clusterBranch: " << clusterBranch->size() << std::endl;
		}
}

void read_tree() {
    TFile *f = new TFile("MLOUTPUT.root");
    TTree *tree = (TTree*) f->Get("myTree");

    std::vector<ClusterCandidate>* clusterBranch = nullptr;
    tree->SetBranchAddress("clusters", &clusterBranch);
    
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        std::cout << "Size of clusterBranch: " << clusterBranch->size() << std::endl;
    }
}



void evaluateClusterTrack(std::vector<ClusterCandidate>& clusterPerChamber, const o2::dataformats::MatchInfoHMP& track, const std::vector<float>& mipCharges, int mcTrackPdg, int trackNumber)
{

        const auto eventCnt = track.getEvent(); // check it corresponds to entry in loop of events?
        const auto momentum = track.getHmpMom();

        const auto nF = track.getRefIndex(); // ef: aon it only gets the mean value ; TODO: in MatchHMP.cxx get calibration value
        // https://github.com/eflatlan/AliceO2/blob/811fcc6d00b363b1e96e0aa8269d46eed95d879b/Detectors/GlobalTracking/src/MatchHMP.cxx#L432C28-L432C38
        // https://github.com/AliceO2Group/AliceO2/blob/e6b603e4c92f98733ff9f7954100140e72bd99f6/Detectors/HMPID/base/include/HMPIDBase/Param.h#L159C37-L159C64
        

        const auto nQ = 1.5787;
        const auto nG = 1.0005;

        // make static?
        const auto& ckovHyps = calcCherenkovHyp(momentum, nF); 

        // TODO: this must be changed??


        //double radParams[7] = {xRad,yRad,L,thetaP, phiP, randomValue.momentum,    randomValue.mass};


			

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
	
	
				
        CkovTools ckovTools(radParams, refIndexes, MIP, ckovHyps, ckovAngle, eventCnt, mcTrackPdg);

        Printf(" Event%d Track%d  : ckovHyps = <%.3f, %.3f> | <%.3f, %.3f> | <%.3f, %.3f>", eventCnt, ckovTools.getMinCkovPion(),ckovTools.getMaxCkovPion(),ckovTools.getMinCkovKaon(), ckovTools.getMaxCkovKaon(),ckovTools.getMinCkovProton(), ckovTools.getMaxCkovProton(), eventCnt); 

       
        // only consider if adequate momentum? 
        //LOGP(info, "Momentum {} ckovHyps {} {} {}", ckov[0], ckov[1], ckov[2]);
        
        if(TMath::IsNaN(ckovHyps[0])) { 
        	Printf("was isnan!!!");
        	return;
        }
        //numBackgroundPhotons, numFoundActualCkov, numActualCkov, numBackgroundLabeledCkov
        std::array<int, 4> arrayInfo;



        // add charge of MIPS here? to skip them in candidates ? 

        // clusterPerChamber by reference, add element to vector {trackNumber, bitHadronStatus}


        // mcTrackPdg check that it matches with clusterPDG?
        // 
        ckovTools.segment(clusterPerChamber, arrayInfo, track.getTrackIndex(), mipCharges, xMip, yMip, q /*MIP-charge*/, mcTrackPdg, track, trackNumber); // temp --> mapBins
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
