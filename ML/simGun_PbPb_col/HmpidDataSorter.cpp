#include "HMPIDBase/Cluster.h" 
#include <vector>
#include <map>
#include <algorithm>
#include <omp.h>
#include <iostream> // For demo purposes
#include "ReconstructionDataFormats/MatchInfoHMP.h"

class HMPIDDataSorter {

using MatchInfoHMP = o2::dataformats::MatchInfoHMP;

public:
    using ClusterMap = std::map<int, std::vector<Cluster>>;
    // o2::dataformats::MCTruthContainer<o2::MCCompLabel> cluLblArr, *cluLblArrPtr = &cluLblArr;

    // using mcCluMap = std::map<int, std::vector<Cluster>>;


    ClusterMap organizeAndSortClusters(const std::vector<Cluster>& allClusters) {
        ClusterMap clustersByEvent;

        #pragma omp parallel
        {
            ClusterMap localClustersByEvent;

            #pragma omp for nowait
            for (size_t i = 0; i < allClusters.size(); ++i) {
                const auto& cluster = allClusters[i];
                const auto& eventNum = cluster.getEventNumberFromTrack();
                localClustersByEvent[eventNum].push_back(cluster);
            }

            #pragma omp critical
            for (const auto& entry : localClustersByEvent) {                
                clustersByEvent[entry.first].insert(clustersByEvent[entry.first].end(), entry.second.begin(), entry.second.end());
            }
        }

        // Sequentially sort clusters within each event
        for (auto& entry : clustersByEvent) {
            std::sort(entry.second.begin(), entry.second.end(), [](const Cluster& a, const Cluster& b) {
                return a.x() < b.x(); // Replace `x()` with the actual property representing the chamber number
            });
        }

        return clustersByEvent;
    }


    using MatchInfoMap = std::map<int, std::vector<MatchInfoHMP>>;
    using McMatchInfoMap = std::map<int, std::vector<o2::MCCompLabelP>>;

    MatchInfoMap organizeAndSortMatchInfo(const std::vector<MatchInfoHMP>& allMatchInfo, const std::vector<o2::MCCompLabelP>& allMcMatchInfo) {

        MatchInfoMap matchInfoByEvent;
        McMatchInfoMap mcMatchInfoByEvent;

        #pragma omp parallel
        {
            MatchInfoMap localMatchInfoByEvent;
            McMatchInfoMap localMcMatchInfoByEvent;

            #pragma omp for nowait
            for (size_t i = 0; i < allMatchInfo.size(); ++i) {
                const auto& matchInfo = allMatchInfo[i];
                const auto& mcMatchInfo = allMcMatchInfo[i];

                const auto& eventNum = matchInfo.getEventNumberFromTrack();
                localMatchInfoByEvent[eventNum].push_back(matchInfo);
                localMcMatchInfoByEvent[eventNum].push_back(mcMatchInfo);

            }

            #pragma omp critical
            for (const auto& entry : localMatchInfoByEvent) {
                matchInfoByEvent[entry.first].insert(matchInfoByEvent[entry.first].end(), entry.second.begin(), entry.second.end());
            }
            
            #pragma omp critical
            for (const auto& entry : localMcMatchInfoByEvent) {
                mcMatchInfoByEvent[entry.first].insert(mcMatchInfoByEvent[entry.first].end(), entry.second.begin(), entry.second.end());
            }

        }

        // Sequentially sort MatchInfo within each event
        for (auto& entry : matchInfoByEvent) {
            std::sort(entry.second.begin(), entry.second.end(), [](const dataformats::MatchInfoHMP& a, const dataformats::MatchInfoHMP& b) {
                return a.getIdxHMPClus() < b.getIdxHMPClus(); // Assuming mIdxHMPClus can serve as a proxy for chamber sorting
            });
        }

        return matchInfoByEvent;
    }
};