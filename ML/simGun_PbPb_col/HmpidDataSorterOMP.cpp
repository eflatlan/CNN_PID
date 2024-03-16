#ifndef HMPID_DATA_SORTER_H
#define HMPID_DATA_SORTER_H

#include "DataFormatsHMP/Cluster.h"
#include "ReconstructionDataFormats/MatchInfoHMP.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include <vector>
#include <map>
#include <algorithm>
#include <omp.h>
#include <iostream>
#include "Steer/MCKinematicsReader.h"


class HMPIDDataSorter2 {

    o2::steer::MCKinematicsReader mcReader; // reader of MC information

    using MatchInfoHMP = o2::dataformats::MatchInfoHMP;
    using MCLabel = o2::MCCompLabel;
    using Cluster = o2::hmpid::Cluster;

    //"3D" structure for Clusters
    using ChamberClustersMap = std::map<int, std::vector<Cluster>>;
    using EventChamberClustersMap = std::map<int, ChamberClustersMap>;

    //"3D" structure for mathcinfo
    using ChamberMatchInfoMap = std::map<int, std::vector<MatchInfoHMP>>;
    using EventChamberMatchInfoMap = std::map<int, ChamberMatchInfoMap>;

    //"3D" structure for mathcinfo-MC-truth
    using ChamberMCLabelMap = std::map<int, std::vector<MCLabel>>;
    using EventChamberMCLabelMap = std::map<int, ChamberMCLabelMap>;

    std::vector<o2::MCTrack> mcTrackArr, *mcTrackArrPtr = &mcTrackArr;


  public:
    HMPIDDataSorter2(const std::vector<Cluster>& allClusters, const std::vector<MatchInfoHMP>& allMatchInfo, const std::vector<MCLabel>& allMcLabels) {
                
        organizeAndSortMatchInfo(allMatchInfo, allMcLabels, matchInfoByEventChamber, mcMatchInfoByEventChamber);
        organizeAndSortClusters(allClusters, clustersByEventChamber);
        
    }

    void iterateOverMatchedTracks() {



        // to get sigmaRing
        const double refIndexTmp = 1.2904;
        Alisigma2_ alisSigma2(refIndexTmp);

        for (const auto& eventEntry : matchInfoByEventChamber) {
            for (const auto& chamberEntry : eventEntry.second) {
                
                // here all the matched and unmatched tracks for a given event and chamber
            
                // get clusters for the given event and chamber
                const auto& clusters = clustersByEventChamber[eventEntry.first][chamberEntry.first];
                const auto& mcMatchInfoArr = mcMatchInfoByEventChamber[eventEntry.first][chamberEntry.first];

                std::vector<SimpleCluster> simpleClusters;
                for (const auto& cluster : clusters) {
                    simpleClusters.emplace_back(cluster.x(), cluster.y(), cluster.q(), cluster.ch());                    
                } 

                int trackcNum = 0;
                
                const o2::MCTrack* mcTrack = nullptr;
                for (const auto& matchInfo : chamberEntry.second) {

                    LOGP(info, "Event {}, chamber | Track number {} : matched Status {}", eventEntry.first, chamberEntry.first, trackcNum, matchInfo.getMatchedStatus());

                    // get MC truth 
                    int trackIdKine = -1;
                    int eventIdKine = -1;
                    int sourceIdKine = -1;
                    trackIdKine = matchInfo.getTrackID();
                    eventIdKine = matchInfo.getEventID();
                    sourceIdKine = matchInfo.getSourceID();   


                    // readMcTrack(eventIdKine, trackIdKine, sourceIdKine);
                    // treeKine->GetEntry(eventIdKine);
                    LOGP(info, "From matchInfo | Event: {}, track: {}, source: {}", eventIdKine, trackIdKine, sourceIdKine);

                    const auto& mcMatchInfo = mcMatchInfoArr[trackNum];
                    mcTrack = mcReader.getTrack(mcMatchInfo);        // mcTrack = mcReader.getTrack(lbl)

                    if(mcTrack == nullptr) {
                        LOGP(info, "MC track not found for event {} and track {}", eventIdKine, trackIdKine);
                        continue;
                    }  

                    TParticlePDG* pPDG = TDatabasePDG::Instance()->GetParticle(mcTrack->GetPdgCode());
                    
                    auto pdgCode = mcTrack.GetPdgCode();
                    
                    auto pT = TMath::Abs(mcTrack.GetStartVertexMomentumX() *
                                        mcTrack.GetStartVertexMomentumX() +
                                    mcTrack.GetStartVertexMomentumY() *
                                        mcTrack.GetStartVertexMomentumY());
               
                    int pdgFromMip = track.getMipClusEventPDG();
                    int pdgFromTrack = track.getParticlePdg();

                    

                    LOGP(info, "Particle {}: mcTrack pPDG = {}, pdgCode = {} | pdgFromMip {} pdgFromTrack {}", i, pPDG, pdgCode, pdgFromMip, pdgFromTrack);


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

    // Define the "3D" structures for MatchInfoHMP and MC labels


    // Function to organize and sort MatchInfoHMP objects
    EventChamberMatchInfoMap organizeAndSortMatchInfo(const std::vector<MatchInfoHMP>& allMatchInfo, const std::vector<MCLabel>& allMcLabels, EventChamberMatchInfoMap& matchInfoByEventChamber, EventChamberMCLabelMap& mcMatchInfoByEventChamber) {

        #pragma omp parallel
        {
            EventChamberMatchInfoMap localMatchInfoByEventChamber;
            EventChamberMCLabelMap localMcMatchInfoByEventChamber;

            #pragma omp for nowait
            for (size_t i = 0; i < allMatchInfo.size(); ++i) {
                const auto& matchInfo = allMatchInfo[i];
                const auto& mcMatchInfo = allMcLabels[i];

                int eventNum = matchInfo.getEventNumberFromTrack();
                int chamberNum = matchInfo.getChamber();

                localMatchInfoByEventChamber[eventNum][chamberNum].push_back(matchInfo);
                localMcMatchInfoByEventChamber[eventNum][chamberNum].push_back(mcMatchInfo);


            }



            // insert to global nested MC map

            #pragma omp critical
            for (const auto& eventEntry : localMatchInfoByEventChamber) {
                for (const auto& chamberEntry : eventEntry.second) {
                    matchInfoByEventChamber[eventEntry.first][chamberEntry.first].insert(
                        matchInfoByEventChamber[eventEntry.first][chamberEntry.first].end(),
                        chamberEntry.second.begin(),
                        chamberEntry.second.end()
                    );
                }
            }

            

            // insert to global nested MC map
            #pragma omp critical
            for (const auto& eventEntry : localMcMatchInfoByEventChamber) {
                for (const auto& chamberEntry : eventEntry.second) {
                    mcMatchInfoByEventChamber[eventEntry.first][chamberEntry.first].insert(
                        mcMatchInfoByEventChamber[eventEntry.first][chamberEntry.first].end(),
                        chamberEntry.second.begin(),
                        chamberEntry.second.end()
                    );
                }
            }


        }

        // Optionally, sort MatchInfoHMP objects within each chamber of each event
        // Implement sorting logic as per your requirements

        return matchInfoByEventChamber, mcMatchInfoByEventChamber;
    }


    std::vector<Cluster> getClusters(int eventNum, int chamberNum) {
        return clustersByEventChamber[eventNum][chamberNum];
    }

    std::vector<MatchInfoHMP> getMatchInfo(int eventNum, int chamberNum) {
        return matchInfoByEventChamber[eventNum][chamberNum];
    }

    std::vector<MCLabel> getMcMatchInfo(int eventNum, int chamberNum) {
        return mcMatchInfoByEventChamber[eventNum][chamberNum];
    }


    // Function to organize and sort Clusters
    EventChamberClustersMap organizeAndSortClusters(const std::vector<Cluster>& allClusters) {
        EventChamberClustersMap clustersByEventChamber;

        #pragma omp parallel
        {
            EventChamberClustersMap localClustersByEventChamber;

            #pragma omp for nowait
            for (size_t i = 0; i < allClusters.size(); ++i) {
                const auto& cluster = allClusters[i];
                int eventNum = cluster.getEventNumber(); // Assuming getEventNumber() exists in Cluster
                int chamberNum = cluster.getChamber(); // Assuming getChamber() exists in Cluster

                localClustersByEventChamber[eventNum][chamberNum].push_back(cluster);
            }

            #pragma omp critical
            for (const auto& eventEntry : localClustersByEventChamber) {
                for (const auto& chamberEntry : eventEntry.second) {
                    clustersByEventChamber[eventEntry.first][chamberEntry.first].insert(
                        clustersByEventChamber[eventEntry.first][chamberEntry.first].end(),
                        chamberEntry.second.begin(),
                        chamberEntry.second.end()
                    );
                }
            }
        }

        // Optionally, sort Clusters within each chamber of each event
        for (auto& eventEntry : clustersByEventChamber) {
            for (auto& chamberEntry : eventEntry.second) {
                std::sort(chamberEntry.second.begin(), chamberEntry.second.end(), [](const Cluster& a, const Cluster& b) {
                    // Define your sorting criteria for Clusters within the same chamber
                    return a.getSomeProperty() < b.getSomeProperty(); // Replace 'getSomeProperty' with your actual sorting criterion
                });
            }
        }

        return clustersByEventChamber;
    }


  private:
    EventChamberMatchInfoMap matchInfoByEventChamber;
    EventChamberMCLabelMap mcMatchInfoByEventChamber;    
    EventChamberClustersMap clustersByEventChamber;


};
#endif