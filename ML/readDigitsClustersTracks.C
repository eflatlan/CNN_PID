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



void process()
{
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
            } else {
                sortedTracks[iCh].push_back(obj);
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

	    		const auto& iCh = obj.ch();
          if (iCh >= 0 && iCh <= 6) {
						if(sortedTracks[iCh].size() > 0) {
				      // make a light copy of digits, just holding the fields charge, x, y
				      /*std::vector<ShallowDigit> shallowDigits;


				      //const std::vector<o2::hmpid::Cluster::Topology>& topology = obj.getClusterTopology();  // some info about digits associated w cluster*/


				      std::vector<std::pair<int,int>> candStatus = {{0,0}};
							// Printf("ClusterCandidate Ch %d", iCh);
                        ClusterCandidate temp(obj.ch(), obj.x(), obj.y(), obj.q(), obj.chi2(), obj.xe(), obj.ye(), obj.getPDG(), candStatus);
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

                
                // get MCTrack correspondign to trackId
                const auto mcTrackIndex = track.getTrackIndex();

                // find the PDG code in teh o2Kine_sim.root file by matching the mcTrackIndex for the current event ; 
                const o2::MCTrack* mcTrack = HmpidDataReader::getMCEntry(mcTracks, mcTrackIndex);

                const int mcTrackPdg = mcTrack->GetPdgCode();
                
                const int momentum = mcTrack->GetPdgCode();
                
                
                float xRad,  yRad,  xPc,  yPc,  th,  ph;
                float xMip = track.getMipX(), yMip = track.getMipY();
                track.getHMPIDtrk(xRad,  yRad,  xPc,  yPc,  th,  ph);

                thTrackMip[i]->Fill(xMip, yMip);
                thTrackPc[i]->Fill(xPc, yPc);
                thTrackRad[i]->Fill(xRad, yRad);

                
                // add mcTrackPdg
                evaluateClusterTrack(clusterPerChamber, track, mipCharges, mcTrackPdg);
            }



        }
    }
}
