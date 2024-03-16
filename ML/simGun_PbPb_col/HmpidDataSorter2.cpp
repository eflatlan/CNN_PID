#ifndef HMPID_DATA_SORTER2_H
#define HMPID_DATA_SORTER2_H
#include <cxxabi.h>  // For abi::__cxa_demangle
#include <cstdlib>   // For free
#include <typeinfo>  // For typeid
#include <iomanip>
#include "Alisigma2_.cpp"
#include "SimpleCluster.cpp"

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

        for (const auto& eventEntry : matchInfoByEventChamber) {
            LOGP(info, "\n=======================================");
            LOGP(info, "  Event {}", eventEntry.first);
            LOGP(info, "=======================================");

            for (const auto& chamberEntry : eventEntry.second) {
                LOGP(info, "=======================================");

                // here all the matched and unmatched tracks for a given event and chamber
            
                // get clusters for the given event and chamber
                // const auto& clusters = clustersByEventChamber[eventEntry.first][chamberEntry.first];
                const auto& clusters = getClusters(eventEntry.first, chamberEntry.first);

                const auto& mcMatchInfoArr = mcMatchInfoByEventChamber[eventEntry.first][chamberEntry.first];

                std::vector<SimpleCluster> simpleClusters;
                for (const auto& cluster : clusters) {
                    simpleClusters.emplace_back(cluster.x(), cluster.y(), cluster.q(), cluster.ch());                    
                } 

                LOGP(info, "    Event {}, chamber {} : num Clusters = {}", eventEntry.first, chamberEntry.first, clusters.size());

                int trackNum = 0;
                
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
                    
                    
                    LOGP(info, "        =========================================================");
                    LOGP(info, "        Track number {} : matched Status {}", trackNum, matchInfo.getMatchStatus());
                                        
                    if(dist > 10) {
                        Printf("               Too large distance %.0f : Skip", dist);
                        trackNum++;
                        continue;                    
                    }



                    // LOGP(info, "        Dist PC CONC-UNC: deltaX {} deltaY {}", xPcCon - xPcUnc, yPcCon - yPcUnc);
                    Printf("               PC_CONC-MIP : deltaX %.1f deltaY %.1f", xPcCon - xMip, yPcCon - yMip);
                    Printf("               dist %.1f", dist);
                    Printf("               xPc %.1f, mip %.1f, ra %.1f | yPc %.1f mip %.1f ra %.1f", xPcCon, xMip, xRa, yPcCon, yMip, yRa);


                    // get MC truth 
                    int trackIdKine = -1;
                    int eventIdKine = -1;
                    int sourceIdKine = -1;

                    const auto& mcMatchInfo = mcMatchInfoArr[trackNum];

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
                    LOGP(info, "        Type of mcMatchInfo: {}", (status == 0 ? demangled : typeid(mcMatchInfo).name()));
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
                      LOGP(info, "      PDG code ulik! : num Clusters = {}", clusters.size());
                      
                      for(const auto& clu : clusters){
                       const auto cluDist2Mip = TMath::Sqrt((clu.x() - xMip)*(clu.x() - xMip) + (clu.y() - yMip)*(clu.y() - yMip));
                       Printf("      cluDist2Mip : %.1f", cluDist2Mip);

                       if(cluDist2Mip < 10) {
                         const auto tid = clu.getTrackId();
                         const auto mid = clu.getMotherId();
                         const auto sid = clu.getSourceId();
                          
                         const auto mcFromClu = mcReader->getTrack(eventIdKine, tid);
                         const auto mcFromCluMother = mcReader->getTrack(eventIdKine, mid);

                         Printf("              Cluster PDG %d, q %d size %d", clu.getPDG(), clu.q(), clu.size());


                         LOGP(info, "        Cluster MC INFO : tid {} mid {} sid {}", tid, mid, sid);

                         if(mcFromClu && mcFromCluMother) {
                            auto pdgTid = mcFromClu->GetPdgCode();
                            auto pdgMid = mcFromCluMother->GetPdgCode();

                            LOGP(info, "        Cluster MC INFO : From tid {} mid {}", pdgTid, pdgMid);
                         }


                       }

                      }
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
        
        LOGP(info, "organizeAndSortClusters : numClusters = {}", allClusters.size());

        for (const auto& cluster : allClusters) {
            int eventNum = cluster.getEventNumber();
            int chamberNum = cluster.ch();
            clustersByEventChamber[eventNum][chamberNum].push_back(cluster);
            LOGP(info, "organizeAndSortClusters : eventNum {} chamberNum {} : size {}", eventNum, chamberNum, clustersByEventChamber[eventNum][chamberNum].size());
            LOGP(info, "organizeAndSortClusters : eventNum {} chamberNum {} : size {}", eventNum, chamberNum, getClusters(eventNum,chamberNum).size());

        }

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


  private:
    EventChamberMatchInfoMap matchInfoByEventChamber;
    EventChamberMCLabelMap mcMatchInfoByEventChamber;    
    EventChamberClustersMap clustersByEventChamber;


};
#endif