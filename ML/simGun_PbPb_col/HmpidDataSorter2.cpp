#ifndef HMPID_DATA_SORTER2_H
#define HMPID_DATA_SORTER2_H
#include <cxxabi.h>  // For abi::__cxa_demangle
#include <cstdlib>   // For free
#include <typeinfo>  // For typeid
#include <iomanip>
#include "Alisigma2_.cpp"
#include "SimpleCluster.cpp"


#include "DataFormatsHMP/Trigger.h"
#include "DataFormatsHMP/Cluster.h"
#include "ReconstructionDataFormats/MatchInfoHMP.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include <vector>
#include <map>
#include <algorithm>
// #include <omp.h>
#include <iostream>
#include "Steer/MCKinematicsReader.h"


class HmpidDataSorter2 {

    // o2::steer::MCKinematicsReader mcReader; // reader of MC information
    std::unique_ptr<o2::steer::MCKinematicsReader> mcReader;
    using MatchInfoHMP = o2::dataformats::MatchInfoHMP;
    using MCLabel = o2::MCCompLabel;
    using Cluster = o2::hmpid::Cluster;

    //"3D" structure for Clusters
    using ChamberClustersMap = std::map<int, std::vector<o2::hmpid::Cluster>>;
    using EventChamberClustersMap = std::map<int, ChamberClustersMap>;

    //"3D" structure for mathcinfo
    using ChamberMatchInfoMap = std::map<int, std::vector<o2::dataformats::MatchInfoHMP>>;
    using EventChamberMatchInfoMap = std::map<int, ChamberMatchInfoMap>;

    //"3D" structure for mathcinfo-MC-truth
    using ChamberMCLabelMap = std::map<int, std::vector<o2::MCCompLabel>>;
    using EventChamberMCLabelMap = std::map<int, ChamberMCLabelMap>;

    std::vector<o2::MCTrack> mcTrackArr, *mcTrackArrPtr = &mcTrackArr;


  public:
    //HMPIDDataSorter2(/*const std::vector<o2::hmpid::Cluster>& allClusters, const std::vector<o2::dataformats::MatchInfoHMP>& allMatchInfo, const std::vector<o2::MCCompLabel>& allMcLabels*/) {
    HmpidDataSorter2() { 
      mcReader = std::make_unique<o2::steer::MCKinematicsReader>("collisioncontext.root");        
    }

    void iterateOverMatchedTracks() {
        LOGP(info, "\n=======================================");
        LOGP(info, "iterateOverMatchedTracks");

        // to get sigmaRing
        const double refIndexTmp = 1.2904;
        Alisigma2_ alisSigma2(refIndexTmp);
        int trigNum = 0;
        for (const auto& eventEntry : matchInfoByEventChamber) {
            LOGP(info, "\n=======================================");
            LOGP(info, "  Event {}", eventEntry.first);
            LOGP(info, "=======================================");
        
            // const auto& trig = triggers[trigNum++]; // eventEntry.first?

            const auto& trig = triggers[eventEntry.first]; // eventEntry.first?


            std::vector<o2::hmpid::Cluster> clustersInEvent;

            // std::ordered_map<int, std::vector<o2::DataFormatsHMP::cluster>> clusterMaps;
            std::array<std::vector<o2::hmpid::Cluster>, 7> clusterArray;
            int tnum = 0;
            int cluNumPre = 0;
            for(int cluNum = trig.getFirstEntry(); cluNum < trig.getLastEntry(); cluNum++)
            {

                if(trig.getNumberOfObjects()<1) {
                  continue;
                }


                if (cluNum < mClusters.size() ) {
                  const auto& clu = mClusters[cluNum]; 

                  const int evNum = clu.getEventNumber();

                  const int chNum = clu.ch();
                  if (chNum >= o2::hmpid::Param::EChamberData::kMinCh && chNum <= o2::hmpid::Param::EChamberData::kMaxCh)
                  {

                    clustersInEvent.emplace_back(clu);
                    clusterArray[chNum].emplace_back(clu);
                  }
                  cluNumPre = evNum;
                }
                tnum++;
            }


            for (const auto& chamberEntry : eventEntry.second) {

                // here all the matched and unmatched tracks for a given event and chamber
            
                // get clusters for the given event and chamber
                // const auto& clusters = clustersByEventChamber[eventEntry.first][chamberEntry.first];
                // const auto& clusters = getClusters(eventEntry.first, chamberEntry.first);
                const auto& clustersInChamber = clusterArray[chamberEntry.first];
                const auto& mcMatchInfoArr = mcMatchInfoByEventChamber[eventEntry.first][chamberEntry.first];



                std::vector<SimpleCluster> simpleClusters;
                for (const auto& cluster : clustersInChamber) {
                    simpleClusters.emplace_back(cluster.x(), cluster.y(), cluster.q(), cluster.ch());                    
                } 

                LOGP(info, "clustersInEvent size {} : in chamber-- simpleClusters size {}", clustersInEvent.size(), simpleClusters.size());


                int trackNum = 0;
                LOGP(info, "=======================================");
                LOGP(info, "    Event {}, chamber {} : num Clusters = {}", eventEntry.first, chamberEntry.first, clustersInChamber.size());
 
                const o2::MCTrack* mcTrack = nullptr;
                for (const auto& matchInfo : chamberEntry.second) {

                    float xMip, yMip;
                    int qMip, nph;

                    float xPcUnc, yPcUnc;
                    float xPcCon, yPcCon, xRa, yRa, theta, phi;
                    matchInfo.getHMPIDtrk(xRa, yRa, xPcCon, yPcCon, theta, phi);
                    matchInfo.getUnconstrainedPc(xPcUnc, yPcUnc);
                    matchInfo.getHMPIDmip(xMip, yMip, qMip, nph);

                    const auto dist = TMath::Sqrt((xPcCon - xMip)*(xPcCon - xMip) + (yPcCon - yMip)*(yPcCon - yMip));
                    
                    
                   
                    if(dist > 6.) {
                        // Printf("               Too large distance %.0f : Skip", dist);
                        trackNum++;
                        continue;                    
                    }

                    LOGP(info, "        =========================================================");
                    LOGP(info, "        Track number {} : matched Status {}", trackNum, matchInfo.getMatchStatus());
                      


                    // LOGP(info, "        Dist PC CONC-UNC: deltaX {} deltaY {}", xPcCon - xPcUnc, yPcCon - yPcUnc);
                    Printf("               PC_CONC-MIP : deltaX %.1f deltaY %.1f", xPcCon - xMip, yPcCon - yMip);
                    Printf("               dist %.1f", dist);
                    Printf("               xPc %.1f, mip %.1f, ra %.1f | yPc %.1f mip %.1f ra %.1f", xPcCon, xMip, xRa, yPcCon, yMip, yRa);


                    // get MC truth 
                    int trackIdKine = -1;
                    int eventIdKine = -1;
                    int sourceIdKine = -1;

                    const auto& mcMatchInfo = mcMatchInfoArr[trackNum];

                    if (!mcMatchInfo.isValid() ) {
                        LOGP(info, "        mcMatchInfo not valid");
                    }

                    trackIdKine = mcMatchInfo.getTrackID();
                    eventIdKine = mcMatchInfo.getEventID();
                    sourceIdKine = mcMatchInfo.getSourceID();   


                    // readMcTrack(eventIdKine, trackIdKine, sourceIdKine);
                    // treeKine->GetEntry(eventIdKine);
                    LOGP(info, "        From mcMatchInfo | Event: {}, track: {}, source: {}", eventIdKine, trackIdKine, sourceIdKine);
                    // LOGP(info, "        Class name of mcMatchInfo: {}", typeid(mcMatchInfo).name());


                    LOGP(info, "        DataSorter2 : try to read MC-track");

                    // mcTrack = mcReader->getTrack(mcMatchInfo);        // mcTrack = mcReader->getTrack(lbl)
                    int status;
                    char* demangled = abi::__cxa_demangle(typeid(mcMatchInfo).name(), 0, 0, &status);
                    //LOGP(info, "        Type of mcMatchInfo: {}", (status == 0 ? demangled : typeid(mcMatchInfo).name()));
                    free(demangled);  

                    if(mcReader == nullptr) {
                      LOGP(info, "      mcReader == nullptr");
                    }

                    try {
                        LOGP(info, "        try mcReader w IP types: trackIdKine: {} {}, eventIdKine: {} {}, sourceIdKine: {} {}", trackIdKine, typeid(trackIdKine).name(), eventIdKine, typeid(eventIdKine).name(), sourceIdKine, typeid(sourceIdKine).name());
                         
                        mcTrack = mcReader->getTrack(mcMatchInfo);                        
                        // mcTrack = mcReader->getTrack(sourceIdKine, eventIdKine, trackIdKine);
                    } catch (const std::exception& e) {
                        LOGP(error, "       Exception caught while trying to read MC track: %s", e.what());
                    } catch (...) {
                        LOGP(error, "       Unknown exception caught while trying to read MC track");
                    }
                    
                    LOGP(info, "        From mcTrack | Event: {}, track: {}, source: {}", eventIdKine, trackIdKine, sourceIdKine);

                    if(mcTrack == nullptr) {
                        LOGP(info, "        MC track not found for event {} and track {}", eventIdKine, trackIdKine);
                        continue;
                    }  

                    TParticlePDG* pPDG = TDatabasePDG::Instance()->GetParticle(mcTrack->GetPdgCode());
                    
                    auto pdgCode = mcTrack->GetPdgCode();
                    
                    auto pT = TMath::Abs(mcTrack->GetStartVertexMomentumX() *
                                        mcTrack->GetStartVertexMomentumX() +
                                    mcTrack->GetStartVertexMomentumY() *
                                        mcTrack->GetStartVertexMomentumY());
               
                    int pdgFromMip = matchInfo.getMipClusEventPDG();
                    int pdgFromTrack = matchInfo.getParticlePdg();

                    LOGP(info, "        Particle {}: mcTrack pdgCode = {} | pdgFromMip {} pdgFromTrack {}", trackNum, pdgCode, pdgFromMip, pdgFromTrack);

                    if(pdgCode != pdgFromMip) {
                      LOGP(info, "      PDG code ulik! : num clustersInChamber = {}", clustersInChamber.size());

                      LOGP(info, "      **********************************");
                      LOGP(info, "      Check clusters \n");

                      const int mipcluSize = matchInfo.getMipClusSize();

                      //void setIdxHMPClus(int ch, int idx) { mIdxHMPClus = ch * 1000000 + idx; }
                      //matching.setIdxHMPClus(iCh, index + 1000 * cluSize); // set chamber, index of cluster + cluster size

                      const int mipIndex = matchInfo.getIdxHMPClus();
                      
                      
                      int chTrack = matchInfo.getChamber();

                      int mipch = mipIndex / 1000000; // Extract chamber
                      int remainder = mipIndex % 1000000; // Remainder after removing chamber information
                      int mipSz = remainder / 1000; // Extract cluster size
                      int index = remainder % 1000; // Extract index of cluster                        
                      
                      

                      const int cluTrigStartIndex = trig.getFirstEntry();

                      const int numCluTotal = mClusters.size();

                      const int numCluInTrig = trig.getNumberOfObjects();
                      
                      LOGP(info, "clustersInEvent size {} simpleClusters size {} numCluInTrig", clustersInEvent.size(), simpleClusters.size(), numCluInTrig);

                      LOGP(info,"       cluTrigStartIndex {} numCluTotal {} : numCluInTrig {} indexMIP {}",cluTrigStartIndex,numCluTotal,numCluInTrig,index);

                      const int mipcluCharge = matchInfo.getMipClusCharge();
                      Printf("              MIP PDG %d; x %.1f y %.1f q %d size %d", matchInfo.getMipClusEventPDG(), xMip, yMip, mipcluCharge, mipcluSize);

                      // for aa sjekke index :
                      const int indexTotal = cluTrigStartIndex  + index;
                      Printf("              mipIndex %d mipch %d mipSz %d index (%d/%d)", mipIndex, mipch, mipSz, indexTotal, numCluTotal);

                      const auto& mipFromMatch = mClusters[indexTotal];
                      Printf("              mipFromMatch PDG %d; Chamber %d x %.1f y %.1f q %d size %d", mipFromMatch.getPDG(), mipFromMatch.ch(), mipFromMatch.x(), mipFromMatch.y(), mipFromMatch.q(), mipFromMatch.size());

                      
                      const auto& mipFromMatch2 = clustersInEvent[index];
                      Printf("              mipFromMatch2 PDG %d; Chamber %d x %.1f y %.1f q %d size %d", mipFromMatch2.getPDG(), mipFromMatch2.ch(), mipFromMatch2.x(), mipFromMatch2.y(), mipFromMatch2.q(), mipFromMatch2.size());




                      //for(const auto& clu : clustersInChamber){
                      for(int ind = trig.getFirstEntry(), ind < trig.getLastEntry(); ind++) {

                       const auto& clu = mClusters[ind];
                       const auto cluDist2Mip = TMath::Sqrt((clu.x() - xMip)*(clu.x() - xMip) + (clu.y() - yMip)*(clu.y() - yMip));

                       if(cluDist2Mip < 20) {
                         
                        if(clu.ch() != mipch)
                          continue;

                         const auto tid = clu.getTrackId();
                         const auto mid = clu.getMotherId();
                         const auto sid = clu.getSourceId();
                          
                         const auto mcFromClu = mcReader->getTrack(eventIdKine, tid);
                         const auto mcFromCluMother = mcReader->getTrack(eventIdKine, mid);
                        
                        if(clu.q() < 100. || clu.size() < 2)
                        {
                            continue;
                        }

                        Printf("\n                cluDist2Mip : %.1f", cluDist2Mip);

                        //Printf("              Cluster PDG %d, q %.2f size %d", clu.getPDG(), clu.q(), clu.size());


                        Printf("              Cluster PDG %d; ch %d x %.1f y %.1f q %.2f size %d", clu.getPDG(), clu.ch(),clu.x(), clu.y(), clu.q(), clu.size());


                        //void setIdxHMPClus(int ch, int idx) { mIdxHMPClus = ch * 1000000 + idx; }
                        //matching.setIdxHMPClus(iCh, index + 1000 * cluSize); // set chamber, index of cluster + cluster size
                        Printf("              mipIndex %d mipch %d mipSz %d index %d", mipIndex, mipch, mipSz, index);
                    
                        Printf("               matchInfo ; chTrack %d", chTrack);
                        Printf("              MIP PDG %d; x %.1f y %.1f q %d size %d", matchInfo.getMipClusEventPDG(), xMip, yMip, mipcluCharge, mipcluSize);

                        // LOGP
                        // xMip, yMip, qMip, nph

                        LOGP(info, "        Cluster MC INFO : tid {} mid {} sid {}", tid, mid, sid);

                        if(mcFromClu == nullptr) {
                            LOGP(info, "        mcFromClu nullptr");
                            continue;
                        }  



                        if(mcFromCluMother == nullptr) {
                            LOGP(info, "        mcFromCluMother nullptr");
                            continue;
                        }  

                         if(mcFromClu!=nullptr && mcFromCluMother!=nullptr) {
                            auto pdgTid = mcFromClu->GetPdgCode();
                            auto pdgMid = mcFromCluMother->GetPdgCode();
                            LOGP(info, "        Cluster MC INFO PDG: From tid {} mid {}", pdgTid, pdgMid);
                         }


                       }

                      }
                      LOGP(info, "      No valid clusters");
                    }


                    // ef: TODO: add fields and stream to NUMPY/PANDAS
                    //  
                    /*                 
                    // to get thetaCer, phiCer...
                    Recon reconObj(thetaP, phiP, xMip, yMip, xRad, yRad, simpleClusters);

                    phi_cer_values[i] = -13.;
                    theta_cer_values[i] = -13.;
                    sigma_ring_values[i] = -13.;

                    // bool findPhotCkov(double cluX, double cluY, double& thetaCer, double&
                    // phiCer)

                    for(const auto& cluster : simpleClusters) {
                        double x = cluster.x();
                        double y = cluster.y();
                        double q = cluster.q();
                        double ch = cluster.ch();
                        double thetaCer, phiCer, sigmaRing;
                        if (reconObj.findPhotCkov(x, y, thetaCer, phiCer)) {
                            sigmaRing = alisSigma2.sigma2(thetaP, phiP, thetaCer, phiCer);
                        } else {
                            sigmaRing = -13.;
                            thetaCer = -13.;
                            phiCer = -13.;
                        }

                    }
                    */   



                
                    trackNum++;
                    

                }
            }
        }
    }
    void organizeAndSortClusters(const std::vector<o2::hmpid::Cluster>& allClusters) {
        

        for (const auto& cluster : allClusters) {
            mClusters.emplace_back(cluster);
            int eventNum = cluster.getEventNumberFromTrack();
            int chamberNum = cluster.ch();

            if(eventNum < 0) {
                LOGP(info, "evenNumber {}", eventNum);
            }

            if(chamberNum < 0 || chamberNum > 6) {
                LOGP(info, "chamberNum {}", chamberNum);
            }       

            if(cluster.q() < 0) {
                LOGP(info, "cluster.q() {}", cluster.q());
            }  

            if(cluster.size() < 0) {
                LOGP(info, "cluster.size() {}", cluster.size());
            }  

            clustersByEventChamber[eventNum][chamberNum].push_back(cluster);
            //LOGP(info, "organizeAndSortClusters : eventNum {} chamberNum {} : size {}", eventNum, chamberNum, clustersByEventChamber[eventNum][chamberNum].size());
            //LOGP(info, "organizeAndSortClusters : eventNum {} chamberNum {} : size {}", eventNum, chamberNum, getClusters(eventNum,chamberNum).size());

        }



        int tnum = 0;

        LOGP(info, "organizeAndSortClusters : nTriggers {}", triggers.size());

        for(const auto& trig : triggers)
        {


            /*
            if(trig.getNumberOfObjects() < 1) {
                LOGP(info, "no clus in trig");
                continue;
            }*/ 
            LOGP(info, "numClus in trig {}", trig.getNumberOfObjects());

            std::vector<o2::hmpid::Cluster> clustersInEvent;

            // std::ordered_map<int, std::vector<o2::DataFormatsHMP::cluster>> clusterMaps;
            int cluNumPre = 0;
            LOGP(info, " trigRange {}--{}", trig.getFirstEntry(), trig.getLastEntry());

            for(int cluNum = trig.getFirstEntry(); cluNum < trig.getLastEntry(); cluNum++)
            {

                if (cluNum < mClusters.size() ) {
                    const auto& clu = mClusters[cluNum]; 

                    const int evNum = clu.getEventNumber();

                    //if(evNum!=cluNumPre && cluNumPre!=0) 
                    LOGP(info, "tnum {}, evNum {}, cluNumPre {}", tnum, evNum, cluNumPre);
                    
                    /*
                    const int chNum = clu.ch();
                    if (chNum >= o2::hmpid::Param::EChamberData::kMinCh && chNum <= o2::hmpid::Param::EChamberData::kMaxCh)
                    {

                    }*/
                    cluNumPre = evNum;
                } else {
                    LOGP(info, "clusters out of range");
                }
                
            }
            tnum++;
        }




        for (const auto& cluster : mClusters) {
            int eventNum = cluster.getEventNumber();
            int chamberNum = cluster.ch();
            //if(cluster.q()>100 && cluster.size()>2) 
            //  Printf("evenNumber %d chamberNum %d x %.1f y %.1f q %.1f size %d", eventNum, chamberNum,  cluster.x(), cluster.y(), cluster.q(), cluster.size());

            if(eventNum < 0) {
                LOGP(info, "evenNumber {}", eventNum);
            }

            if(chamberNum < 0 || chamberNum > 6) {
                LOGP(info, "chamberNum {}", chamberNum);
            }       
                 
            if(cluster.q() < 0) {
                LOGP(info, "cluster.q() {}", cluster.q());
            }  

            if(cluster.size() < 0) {
                LOGP(info, "cluster.size() {}", cluster.size());
            }  

        }        
        LOGP(info, "organizeAndSortClusters : numClusters in = {}; out = {} ", allClusters.size(), mClusters.size());

    }

    void organizeAndSortMatchInfo(const std::vector<o2::dataformats::MatchInfoHMP>& allMatchInfo, const std::vector<o2::MCCompLabel>& allMcLabels) {
        
        for (size_t i = 0; i < allMatchInfo.size(); ++i) {
            const auto& matchInfo = allMatchInfo[i];
            const auto& mcMatchInfo = allMcLabels[i];

            int eventNum = matchInfo.getEventNumberFromTrack();
            int chamberNum = matchInfo.getChamber();

            matchInfoByEventChamber[eventNum][chamberNum].push_back(matchInfo);
            mcMatchInfoByEventChamber[eventNum][chamberNum].push_back(mcMatchInfo);
        }


        for (size_t i = 0; i < allMatchInfo.size(); ++i) {
            const auto& matchInfo = allMatchInfo[i];
            const auto& mcMatchInfo = allMcLabels[i];

            int eventNum = matchInfo.getEventNumberFromTrack();
            int chamberNum = matchInfo.getChamber();

            matchInfoByEventChamber[eventNum][chamberNum].push_back(matchInfo);
            mcMatchInfoByEventChamber[eventNum][chamberNum].push_back(mcMatchInfo);
        }        

    }


    std::vector<o2::hmpid::Cluster> getClusters(int eventNum, int chamberNum) {
        return clustersByEventChamber[eventNum][chamberNum];
    }

    std::vector<o2::dataformats::MatchInfoHMP> getMatchInfo(int eventNum, int chamberNum) {
        return matchInfoByEventChamber[eventNum][chamberNum];
    }

    std::vector<o2::MCCompLabel> getMcMatchInfo(int eventNum, int chamberNum) {
        return mcMatchInfoByEventChamber[eventNum][chamberNum];
    }


    void setTriggers(std::vector<o2::hmpid::Trigger>* trigArrPtr) {
      if(!trigArrPtr) {
        throw std::runtime_error("trigArrPtr nullptr");
        return;

      }
      for(const auto trig : *trigArrPtr) {
        triggers.emplace_back(trig);
      } 
    }

  private:
    std::vector<o2::hmpid::Cluster> mClusters;
    std::vector<o2::hmpid::Trigger> triggers;
    EventChamberMatchInfoMap matchInfoByEventChamber;
    EventChamberMCLabelMap mcMatchInfoByEventChamber;    
    EventChamberClustersMap clustersByEventChamber;


};
#endif