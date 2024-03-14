#pragma once
#include "CommonDataFormat/InteractionRecord.h"
#include "CommonDataFormat/RangeReference.h"

#include "DataFormatsHMP/Cluster.h"
#include "DataFormatsHMP/Digit.h"
#include "HMPIDBase/Param.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

// #include "CkovToolsSingle.cpp"
#include "HmpidDataReader.cpp"
#include "ParticleUtils2.cpp"

#include <utility>
#include <vector>

#include <TMath.h>
#include <cmath>

float calcCkovFromMass(float p, float n, int pdg);


void evaluateClusterTrack(
    std::vector<o2::hmpid::ClusterCandidate> &clusterPerChamber,
    const o2::dataformats::MatchInfoHMP &track,
    const std::vector<float> &mipCharges, int mcTrackPdg, int trackNumber,
    int &plotNumber);

std::array<float, 3> calcCherenkovHyp(float p, float n);

float sigmaSep = 1.5;

void SegmentationCkov(double _sigmaSep = 1.5) {
  sigmaSep = _sigmaSep;

  int plotNumber = 0;
  ParticleUtils2 mParticleEvents;
  
  
  auto fileOut = std::make_unique<TFile>("MLOUTPUT.root", "RECREATE");
  
  auto treeOut = std::make_unique<TTree>(
      "treeOut", "Tree to store clusters, trackInfo, and mcPDG");

  std::vector<o2::hmpid::ClusterCandidate> *clusterBranch = nullptr;
  o2::dataformats::MatchInfoHMP *trackInfoBranch = nullptr;
  int mcPDGBranch = 0;

  float reconCkovBranch = 0;

  treeOut->Branch("clusters", &clusterBranch);
  treeOut->Branch("trackInfo", &trackInfoBranch);
  treeOut->Branch("mcPDG", &mcPDGBranch);

  int cluChargeBranch, cluSizeBranch;

  float refIndexBranch, xRadBranch, yRadBranch, xMipBranch, yMipBranch;
  float thBranch, phBranch, pBranch;

  treeOut->Branch("reconCkov", &reconCkovBranch);
  treeOut->Branch("cluCharge", &cluChargeBranch);
  treeOut->Branch("cluSize", &cluSizeBranch);
  treeOut->Branch("refIndex", &refIndexBranch);
  treeOut->Branch("xRad", &xRadBranch);
  treeOut->Branch("yRad", &yRadBranch);
  treeOut->Branch("xMip", &xMipBranch);
  treeOut->Branch("yMip", &yMipBranch);
  treeOut->Branch("th", &thBranch);
  treeOut->Branch("ph", &phBranch);
  treeOut->Branch("p", &pBranch);

  // clusters and triggers
  std::vector<Cluster> *clusterArr = nullptr;
  /*std::vector<o2::hmpid::Topology> mTopologyFromFile,
      *mTopologyFromFilePtr = &mTopologyFromFile;*/ 




	auto matchFileName = "o2match_hmp.root";
	auto cluFileName = "hmpidclusters.root";
	auto mcFileName = "o2sim_Kine.root";
	
	
	/*
	TFile *fileKine = TFile::Open("o2sim_Kine.root");
	 	  
  TFile *fileClu = TFile::Open("hmpidclusters.root");
  */
  
  
	HmpidDataReader hmpidDataReader(matchFileName, cluFileName, mcFileName);


	/*
  TTree *treeCluster = HmpidDataReader::initializeClusterTree(
      clusterArr, trigArr);
  // clusterArr now initialized correctly

  // MathcInfoHMP : holding trackinfo
  std::vector<o2::dataformats::MatchInfoHMP> *matchArr = nullptr;
  TTree *treeMatch = HmpidDataReader::initializeMatchTree(matchArr, 0, 0, 0);

  // McTrack : holding PDG code of track
  //std::vector<o2::MCTrack> *mcArr = nullptr;
  //TTree *tMcTrack = HmpidDataReader::initializeMCTree(mcArr);
	*/ 
	
	
	// (int eventID, int trackID, int pdg)
	
	// must be called in loop ? 
	// hmpidDataReader.initializeMatchTree(0, 0, 0);
	
	
	std::vector<Trigger>* trigArr = hmpidDataReader.getTrigArr();
	
	
  int startIndexTrack = 0;
  if (trigArr == nullptr) {
    Printf(" trigArr== nullptr");
    return;
  }


  // Arrays of triggers file 
  for (int i = 0; i < trigArr->size();
       i++) // for(const auto& clusters : clustersVector) // "events loop"
  {
    auto pTgr = &trigArr->at(i);
    if (pTgr == nullptr) {
      Printf("pTgr== nullptr");
      continue;
    }
    
  	/* OLD HmpidDataReader usage :  
    std::vector<Cluster>* clusterArr = hmpidDataReader.getCluArr();
    const int firstCluInTrig = pTgr->getFirstEntry();
    const int lastCluInTrig = pTgr->getLastEntry();
    */ 

    // Printf("Checking trigger number %d Range Clu = %d :: %d", i, firstCluInTrig,
    // lastCluInTrig);


    // loop over clusters in the trigger
    /*for (int j = pTgr->getFirstEntry(); j <= pTgr->getLastEntry(); j++) {
      const auto &clu = static_cast<o2::hmpid::Cluster>(clusterArr->at(j));
      oneEventClusters.push_back(clu);      
    } */ 

	  auto oneEventClusters = hmpidDataReader.getClusInEvent(/*int event*/ i);
	  
	  if(oneEventClusters.size() < 1) {
	  	continue;
	  }
	  
	  auto firstCluInEvent = static_cast<o2::hmpid::Cluster>(oneEventClusters[0]);
	  
    auto lastClu = static_cast<o2::hmpid::Cluster>(oneEventClusters[-1]);
    
    int eventNumberFirstClu = firstCluInEvent.getEventNumber();
    int eventNumberLast = lastClu.getEventNumber();
	  
    Printf("============\n Cluster Loop \n ===============");
    // find entries in tracksOneEvent which corresponds to correct eventNumber
    Printf("Reading match vector<MatchInfoHMP> for startIndexTrack %d",
           startIndexTrack);


    // Read Information from  ROOT files:
    // 1
    
    LOGP(info, "now reading ReadMatch for trigger {} of {}", i+1, trigArr->size());
    
    
    std::vector<o2::dataformats::MatchInfoHMP> tracksOneEvent =
        hmpidDataReader.readMatch(i, startIndexTrack);
		// readMatch(int eventID, int &startIndex)


    // get MC tracks for given event from mc;
    /*
    Printf("Reading vector<o2::MCTrack>* mcTracks for  eventNumber %d",
           eventNumberFirstClu);
    std::vector<o2::MCTrack> *mcTracks =
        hmpidDataReader.readMC(eventNumberFirstClu); */ 
		/*
		if(!tracksOneEvent) {
		    Printf("tracksOneEvent nullptr");
    	continue;
		}	*/ 
		
		

    LOGP(info, "tracksOneEvent size {}", tracksOneEvent.size());
	
	
    if(tracksOneEvent.size() < 1) {
    	continue;
    }


    // sort the tracksOneEvent vector based on chamber number
    std::sort((tracksOneEvent).begin(), (tracksOneEvent).end(),
              [](const o2::dataformats::MatchInfoHMP &a,
                 const o2::dataformats::MatchInfoHMP &b) {
                return a.getChamber() < b.getChamber();
              });

    // sort the oneEventClusters vector based on chamber number
    std::sort(
        (oneEventClusters).begin(), (oneEventClusters).end(),
        [](const Cluster &a, const Cluster &b) { return a.ch() < b.ch(); });


    std::vector<o2::dataformats::MatchInfoHMP> sortedTracks[7];
    // Assign MatchInfoHMP objects to corresponding vectors based on iCh value
    for (const auto &track : tracksOneEvent) {
      const auto &iCh = track.getChamber();
      if (iCh < 0 || iCh > 6) {
        std::cerr << "Warning: iCh value out of expected range: " << iCh
                  << std::endl;
      } else {
        if (track.getMatchStatus()) {
          sortedTracks[iCh].push_back(track);
        } ;
      }
    }

    std::cout << "Length of sortedTracks vector ";
    for (int iCh = 0; iCh < 7; iCh++) {
      if (sortedTracks[iCh].size() != 0)
        std::cout << iCh << " : " << sortedTracks[iCh].size() << " |";
    }
    
    
    std::cout << std::endl;


    std::vector<o2::hmpid::ClusterCandidate> sortedClusters[7];

    // Assign ClusterCandidate objects to corresponding vectors based on iCh value
    
    LOGP(info, "now looping over clusters ");
    for (const auto &cluster : oneEventClusters) {

      // from each track; we assign a label to each of the clusters in
      // corresponding to the type of hadron it can be
      const auto &iCh = cluster.ch();
      auto& clusterInChamber =  sortedClusters[iCh];
      if (iCh >= 0 && iCh <= 6) {

        // if there is matched tracks for the chamber fill clusters for chamber 
        if (clusterInChamber.size() > 0) {

          o2::hmpid::ClusterCandidate cluCandidate(cluster.ch(), cluster.x(), cluster.y(), cluster.q(),
                                           cluster.chi2(), cluster.xe(), cluster.ye(),
                                           cluster.getPDG(), cluster.size());

          clusterInChamber.emplace_back(cluCandidate);
        } 
      } 
    } // end of loop over clusters in theOP(()			
    // loop over the sortedTracks and sortedClusters
    LOGP(info, "now looping over sortedTracks and sortedClusters ");
    for (int iCh = 0; iCh < 7; iCh++) {

      // if no tracks in the chamber, skip
      if (sortedTracks[iCh].size() < 1) {
        // Printf("sortedTracks[iCh%d].size() %d", i, sortedTracks[i].size());
        continue;
      }


      // std::vector<o2::hmpid::ClusterCandidate> clusterPerChamber
      auto &clusterPerChamber = sortedClusters[iCh];

      // fill charges of MIPs

      /*
      std::vector<float> mipCharges;

      for (const auto &track : sortedTracks[i]) {

        auto q = track.getMipClusQ();
        mipCharges.emplace_back(q);
      }*/

      int tNum = 0;


    	LOGP(info, "now looping over sortedTracks[iCh]");
      //  loop over the tracks in the event
      int trackCount = 1;
      for (const auto &track : sortedTracks[iCh]) {
      
      	LOGP(info, "Track {} of {}", trackCount++, sortedTracks[iCh].size());
        // Printf("TrackNumber%d track[iCh%d].size() %d", tNum++, i,
        // sortedTracks[i].size());

        // Printf("TrackNumberMom %f", track.getHmpMom());

        // pass clusters (and track) by reference, and add its status per track
        // (adding to candStatus vector )

        // for each clusterPerChamber we will have a "candidate-status" for each
        // of the photons, this is a vector of length of sortedTracks[i].size();
        // and holds the fields

        if (track.getMatchStatus() ) {
        
        	
          const auto mcTrackIndex = track.getTrackIndex();
					// LOGP(info, "Track mathced, mcTrackIndex: {} ", mcTrackIndex);
          // find the PDG code in teh o2Kine_sim.root file by matching the
          // mcTrackIndex for the current event ;
          //--((const o2::MCTrack* mcTrack =
          //HmpidDataReader::getMCEntry(mcTracks, mcTrackIndex);

          // const int mcTrackPdg = mcTrack->GetPdgCode();

          // mcPDGBranch = mcTrack->GetPdgCode();

          float xRad, yRad, xPc, yPc, th, ph; // make for these
          float xMip = track.getMipX(), yMip = track.getMipY(); // and thse
          track.getHMPIDtrk(xRad, yRad, xPc, yPc, th, ph);

					LOGP(info, "Track mathced, xRad: {}", xRad);
					LOGP(info, "Track mathced, xPc: {}", xPc);
					LOGP(info, "Track mathced, xMip: {}", xMip);
					LOGP(info, "Track mathced, mipSize: {}", track.getMipClusSize());

          TVector2 mip(xMip, yMip);

          refIndexBranch = track.getRefIndex();
          cluChargeBranch = track.getMipClusQ();
          cluSizeBranch = track.getMipClusSize();

          xRadBranch = xRad; 
          yRadBranch = yRad; 

          xMipBranch = xMip;
          yMipBranch = yMip; 

          thBranch = th; 
          phBranch = ph; 

          // get the calculated recon track ckov signal
          reconCkovBranch = track.getHMPsignal();

          // clusterPerChamber : vector of ClusterCandidate-objects
          // track : object with 10 scalar values
          // mcTradckPDG : MC truth
          //
          int pdg = track.getMipClusEventPDG();
          
          // ikke call denne, istedet run updateH5.C 
          //evaluateClusterTrack(clusterPerChamber, track, mipCharges, pdg, tNum,
          //                     plotNumber);

          // saving...
          // branch : mcTrackPdg
          // branch : clusterPerChamber --> vector of clusterCandidates
          // branch : track --> trackInformation (Ra, MIP, refIndex, momentum,
          // theta, phi)

          clusterBranch =
              const_cast<std::vector<o2::hmpid::ClusterCandidate> *>(
                  &clusterPerChamber); 
          trackInfoBranch = const_cast<o2::dataformats::MatchInfoHMP *>(
              &track); 


					LOGP(info, "Track mathced, clusterBranch");

          mParticleEvents.fillCandidate(clusterPerChamber, track, 0);
					LOGP(info, "Track mathced mParticleEvents fileld ! ");
					
	        // mcPDGBranch = mcTrack->GetPdgCode();

          pBranch = track.getHmpMom();
          // Fill the tree
          
          
          
          treeOut->Fill();

					
          // for(auto& clusterPerChamber)
        } else {
          LOGP(info, "Track didnt match!");
        } // end if matched 
        LOGP(info, "end if matched");
        
        // end if matched
      } // end for ( track : sortedTracks[iCh] )
        LOGP(info, "end end for ( track : sortedTracks[iCh] )");
      // save ClusterPerChamber

      // save trackInfo :
      //                  trackIndex : mIdxTrack
      //                  scalar : p, refIndex,
      //                  trackInfo theta, phi, xRad, yRad; {xPc yPc is
      //                  redundant} MIP            : x, y, q MCTruth : pdg code
      //                  of track || pdg code of MIP?

      // now assigned all types of candidates for all clusters in teh given
      // chamber match sortedClusters[i] with sortedTracks[i] -->
      //
    }  // end for (int iCh = 0; iCh < 7; iCh++)
    
    LOGP(info, "end end for (int iCh = 0; iCh < 7; iCh++) )");
  } // end for (int i = 0; i < trigArr->size();
  LOGP(info, "end for (int i = 0; i < trigArr->size()");

	/*
  treeOut->Write();
  
  
  LOGP(info, " treeOut->Write()");

  
  Long64_t nEntries = treeOut->GetEntries();
  LOGP(info, " nEntries {}", nEntries);
  
  for (Long64_t i = 0; i < nEntries; ++i) {
    treeOut->GetEntry(i);


    // std::cout << "Size of clusterBranch: " << clusterBranch->size() <<
    // std::endl;
  }
  fileOut->Close();
  
  
  LOGP(info, "fileOut->Close();"); */ 
  mParticleEvents.writeH5();
  LOGP(info, " mParticleEvents.writeH5();");
  
}

float calcCkovFromMass(float p, float n, int pdg) {
  // Constants for particle masses (in GeV/c^2)
  const float mass_Pion = 0.1396, mass_Kaon = 0.4937, mass_Proton = 0.938;

  float m; // variable to hold the mass

  // Switch based on absolute value of PDG code
  switch (std::abs(pdg)) {
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
    return 0; // return 0 if PDG code doesn't match any known codes
  }

  const float p_sq = p * p;
  const float refIndexFreon = n; // Assuming n is the refractive index
  const float cos_ckov_denom = p * refIndexFreon;

  // Sanity check
  if (p_sq + m * m < 0) {
    return 0;
  }

  const auto cos_ckov =
      static_cast<float>(TMath::Sqrt(p_sq + m * m) / cos_ckov_denom);

  // Sanity check
  if (cos_ckov > 1 || cos_ckov < -1) {
    return 0;
  }

  const auto ckovAngle = static_cast<float>(TMath::ACos(cos_ckov));

  return ckovAngle;
}

void evaluateClusterTrack(
    std::vector<o2::hmpid::ClusterCandidate> &clusterPerChamber,
    const o2::dataformats::MatchInfoHMP &track,
    const std::vector<float> &mipCharges, int mcTrackPdg, int trackNumber,
    int &plotNumber) {

  const auto eventCnt =
      track.getEvent(); // check it corresponds to entry in loop of events?
  const auto momentum = track.getHmpMom();

  const auto nF =
      1.2929 - 0.0025; // track.getRefIndex(); // ef: aon it only gets the mean
                       // value ; TODO: in MatchHMP.cxx get calibration value
  // https://github.com/eflatlan/AliceO2/blob/811fcc6d00b363b1e96e0aa8269d46eed95d879b/Detectors/GlobalTracking/src/MatchHMP.cxx#L432C28-L432C38
  // https://github.com/AliceO2Group/AliceO2/blob/e6b603e4c92f98733ff9f7954100140e72bd99f6/Detectors/HMPID/base/include/HMPIDBase/Param.h#L159C37-L159C64

  const auto nQ = 1.583;
  const auto nG = 1.0005;

  // make static?

  // make account for varying nF ? +- 2 std ?
  const auto &ckovHypsMin = calcCherenkovHyp(momentum, nF);

  const auto &ckovHypsMax = calcCherenkovHyp(momentum, nF);

  float xRad, yRad, xPc, yPc, thetaP, phiP;
  track.getHMPIDtrk(xRad, yRad, xPc, yPc, thetaP, phiP);

  float xMip, yMip;
  int nph, q;
  track.getHMPIDmip(xMip, yMip, q, nph);

  const auto &L = 0.5; //
  double radParams[7] = {
      xRad,
      yRad,
      L,
      thetaP,
      phiP,
      momentum,
      static_cast<double>(mcTrackPdg * 1.0)}; // ef : TODO add PID to MLinfoHMP?

  double refIndexes[3] = {nF, nQ, nG};

  double MIP[3] = {xMip, yMip, static_cast<double>(q * 1.0)};
  // (double radParams[7], double refIndexes[3], double MIP[3],
  //  std::array<float, 3> ckovHyps, float trackCkov, int eventCnt)
  // ef: TODO use MIP to get radius and phi in CkovTools:
  auto ckovAngle = 0.;
  /*
  CkovTools ckovTools(radParams, refIndexes, MIP, ckovHypsMin, ckovHypsMax,
                      ckovAngle, eventCnt, mcTrackPdg, sigmaSep);

  Printf(" Event%d Track%d  : ckovHyps = <%.3f, %.3f> | <%.3f, %.3f> | <%.3f, "
         "%.3f>",
         eventCnt, ckovTools.getMinCkovPion(), ckovTools.getMaxCkovPion(),
         ckovTools.getMinCkovKaon(), ckovTools.getMaxCkovKaon(),
         ckovTools.getMinCkovProton(), ckovTools.getMaxCkovProton(), eventCnt);
  */

  // only consider if adequate momentum?
  // LOGP(info, "Momentum {} ckovHyps {} {} {}", ckov[0], ckov[1], ckov[2]);


  // numBackgroundPhotons, numFoundActualCkov, numActualCkov,
  // numBackgroundLabeledCkov
  std::array<int, 4> arrayInfo;


  // mcTrackPdg check that it matches with clusterPDG?
  //
  /*  ckovTools.segment(clusterPerChamber, arrayInfo, track.getTrackIndex(),
                    xMip, yMip, q / * MIP-charge* /, mcTrackPdg, track,
                    trackNumber, plotNumber); // temp --> mapBins*/
}

const float mass_Pion = 0.1396, mass_Kaon = 0.4937,
            mass_Proton = 0.938; // masses in
std::array<float, 3> masses = {mass_Pion, mass_Kaon, mass_Proton};
const float mass_Pion_sq = mass_Pion * mass_Pion,
            mass_Kaon_sq = mass_Kaon * mass_Kaon,
            mass_Proton_sq = mass_Proton * mass_Proton;


std::array<float, 3> calcCherenkovHyp(float p, float n) {
  const float p_sq = p * p;
  const float cos_ckov_denom = p * n;
  const auto cos_ckov_Pion = static_cast<float>(
      TMath::Sqrt(p_sq + mass_Pion_sq) /
      (cos_ckov_denom)); // n = refIndexFreon 1.289 later make it random?

  const float cos_ckov_Kaon =
      static_cast<float>(TMath::Sqrt(p_sq + mass_Kaon_sq) / (cos_ckov_denom));
  const float cos_ckov_Proton =
      static_cast<float>(TMath::Sqrt(p_sq + mass_Proton_sq) / (cos_ckov_denom));

  const float ckovAnglePion = static_cast<float>(TMath::ACos(cos_ckov_Pion));
  const float ckovAngleKaon = static_cast<float>(TMath::ACos(cos_ckov_Kaon));
  const float ckovAngleProton =
      static_cast<float>(TMath::ACos(cos_ckov_Proton));

  Printf("Pion %.3f Kaon %.3f Proton %.3f", ckovAnglePion, ckovAngleKaon,
         ckovAngleProton);

  return {ckovAnglePion, ckovAngleKaon, ckovAngleProton};
}
