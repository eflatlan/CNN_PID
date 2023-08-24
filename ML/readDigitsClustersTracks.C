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


Int_t pos[]= {9,8,6,5,4,2,1};


TH2F* thTrackMip[7], *thTrackPc[7], *thTrackRad[7];

void readDigits(char* filename, int nEvent)
{
  TFile *fileDigits = TFile::Open(filename);
    std::map<int, int> pdgCounts;  // Map to store PDG values and their counts

  TTree *treeDigits = (TTree*)fileDigits->Get("o2sim");
  
  TH1F *hCharge[7];
  TH2F *hMap[7];
      
    TH2F *hParticlePdgMapAbove[7];
    TH2F *hParticlePdgMapBelow[7];


TH2F* hParticlePdgElectron[7], *hParticlePdgPion[7], *hParticlePdgKaon[7],  *hParticlePdgProton[7], *hParticlePdgPhoton[7], *hParticlePdgPhoton2[7];


//11 111 321 2212 22/50000050
for (int i = 0; i < 7; i++) {
    // Electron histograms
    hParticlePdgElectron[i] = new TH2F(Form("Electron_ParticlePdg_Chamber%i", i),
                                       Form("Electron PDG Map Chamber%i", i),
                                       160, 0, 159, 144, 0, 143);
    hParticlePdgElectron[i]->SetXTitle("pad X");
    hParticlePdgElectron[i]->SetYTitle("pad Y");
    hParticlePdgElectron[i]->SetMarkerColor(kRed); // Example color

    // Pion histograms
    hParticlePdgPion[i] = new TH2F(Form("Pion_ParticlePdg_Chamber%i", i),
                                   Form("Pion PDG Map Chamber%i", i),
                                   160, 0, 159, 144, 0, 143);
    hParticlePdgPion[i]->SetXTitle("pad X");
    hParticlePdgPion[i]->SetYTitle("pad Y");
    hParticlePdgPion[i]->SetMarkerColor(kBlue-2); // Example color
    hParticlePdgPion[i]->SetMarkerStyle(2); // Example color
    // Kaon histograms
    hParticlePdgKaon[i] = new TH2F(Form("Kaon_ParticlePdg_Chamber%i", i),
                                   Form("Kaon PDG Map Chamber%i", i),
                                   160, 0, 159, 144, 0, 143);
    hParticlePdgKaon[i]->SetXTitle("pad X");
    hParticlePdgKaon[i]->SetYTitle("pad Y");
    hParticlePdgKaon[i]->SetMarkerColor(kBlue); // Example color
    hParticlePdgKaon[i]->SetMarkerStyle(2); // Example color
    // Proton histograms
    hParticlePdgProton[i] = new TH2F(Form("Proton_ParticlePdg_Chamber%i", i),
                                     Form("Proton PDG Map Chamber%i", i),
                                     160, 0, 159, 144, 0, 143);
    hParticlePdgProton[i]->SetXTitle("pad X");
    hParticlePdgProton[i]->SetYTitle("pad Y");
    hParticlePdgProton[i]->SetMarkerColor(kBlue+2); // Example color
    hParticlePdgProton[i]->SetMarkerStyle(2); // Example color
    // Photon histograms
    hParticlePdgPhoton[i] = new TH2F(Form("Photon_ParticlePdg_Chamber%i", i),
                                     Form("Photon PDG Map Chamber%i", i),
                                     160, 0, 159, 144, 0, 143);
    hParticlePdgPhoton[i]->SetXTitle("pad X");
    hParticlePdgPhoton[i]->SetYTitle("pad Y");
    hParticlePdgPhoton[i]->SetMarkerColor(kBlack); // Example color

    hParticlePdgPhoton2[i] = new TH2F(Form("Photon2_ParticlePdg_Chamber%i", i),
                                     Form("Photon2 PDG Map Chamber%i", i),
                                     160, 0, 159, 144, 0, 143);
    hParticlePdgPhoton2[i]->SetXTitle("pad X");
    hParticlePdgPhoton2[i]->SetYTitle("pad Y");
    hParticlePdgPhoton2[i]->SetMarkerColor(kCyan); // Example color

}
// Add a new histogram to collect PDG values
TH1F* hPDG = new TH1F("hPDG", "Particle PDG values", 5000, -2500, 2500);
hPDG->SetXTitle("PDG Value");
hPDG->SetYTitle("Entries");


    for (int i = 0; i < 7; i++) {
        hParticlePdgMapAbove[i] = new TH2F(Form("Digits ParticlePdg Map Above 10000 chamber%i", i), 
                                          Form("Digits ParticlePdg Map Above 10000 chamber%i", i), 
                                          160, 0, 159, 144, 0, 143);
        hParticlePdgMapAbove[i]->SetXTitle("pad X");
        hParticlePdgMapAbove[i]->SetYTitle("pad Y");
        hParticlePdgMapAbove[i]->SetMarkerColor(kRed);

        hParticlePdgMapBelow[i] = new TH2F(Form("Digits ParticlePdg Map Below 10000 chamber%i", i), 
                                          Form("Digits ParticlePdg Map Below 10000 chamber%i", i), 
                                          160, 0, 159, 144, 0, 143);
        hParticlePdgMapBelow[i]->SetXTitle("pad X");
        hParticlePdgMapBelow[i]->SetYTitle("pad Y");
        hParticlePdgMapBelow[i]->SetMarkerColor(kGreen);
    }

  for(int i=0; i<7; i++) {
     hMap[i] = new TH2F(Form("Digits Map chmaber%i",i),Form("Digits Map chmaber%i",i), 160, 0, 159, 144, 0, 143);
     hMap[i]->SetXTitle("pad X");
     hMap[i]->SetYTitle("pad Y"); 
     
     hCharge[i] = new TH1F(Form("Digits Charge chamber%i",i),Form("Digits Charge chamber%i",i),2000, 100., 2100.);
     hCharge[i]->SetXTitle("Charge (ADC channel)");               
     hCharge[i]->SetYTitle("Entries");               
   }
   
  TCanvas *c1 = new TCanvas("c1","c1",1000,800); 
  TCanvas *c2 = new TCanvas("c2","c2",1000,800); 
   
  Int_t pos[]= {9,8,6,5,4,2,1};
  
  o2::hmpid::Trigger *pTgr;
  o2::hmpid::Digit *pDig;
  o2::hmpid::Digit *pDigEvt;
  o2::hmpid::Digit digit;
  
  std::vector<o2::hmpid::Digit> *digits = nullptr;
  std::vector<o2::hmpid::Digit> oneEventDigits;
  std::vector<o2::hmpid::Trigger> *trigger = nullptr;
  
  treeDigits->SetBranchAddress("HMPDigit",&digits);
  treeDigits->SetBranchAddress("InteractionRecords",&trigger);
  
  Printf("tree entries = %i", treeDigits->GetEntries());
  
  treeDigits->GetEntry(0);

  Printf("digit size = %i", digits->size());
  
  int padChX = 0, padChY = 0, module = 0;
    



//11 111 321 2212 22/50000050
  std::vector<int> particlePDGs;
  for(unsigned int j = 0; j < digits->size(); j++) {
    pDig = (o2::hmpid::Digit*)&digits->at(j);
	  if(j<15)
	    std::cout << j << " evnum" << pDig->getEventNumber()<< std::endl; 




    hPDG->Fill(pDig->mParticlePdg);
  pdgCounts[pDig->mParticlePdg]++;
    o2::hmpid::Digit::pad2Absolute(pDig->getPadID(), &module, &padChX, &padChY);    
                         
    hCharge[module]->Fill(pDig->getQ());

        if (pDig->mParticlePdg > 10000) {
            hParticlePdgMapAbove[module]->Fill(padChX, padChY);
        } else {
            hParticlePdgMapBelow[module]->Fill(padChX, padChY);
        }


    switch (TMath::Abs(pDig->mParticlePdg)) {
      case 11 : 
	hParticlePdgElectron[module]->Fill(padChX, padChY); 
	break;
      case 211 : 
	hParticlePdgPion[module]->Fill(padChX, padChY);
	break;
      case 321 : 
	hParticlePdgKaon[module]->Fill(padChX, padChY);
	break;
      case 2212 : 
	hParticlePdgProton[module]->Fill(padChX, padChY);
	break;
      case 50000050 : 
	hParticlePdgPhoton[module]->Fill(padChX, padChY);
	break;
      case 22 : 
	hParticlePdgPhoton2[module]->Fill(padChX, padChY);
	break;
    }
 
    bool isInPdg = false;
    for(const auto& pdg: particlePDGs) { if(pdg == pDig->mParticlePdg) {isInPdg = true;} }
    if(!isInPdg) {particlePDGs.push_back(pDig->mParticlePdg);}    

   }

  Printf("PDG code counts:");
  for (const auto& pair : pdgCounts) {
    Printf("PDG: %d, Count: %d", pair.first, pair.second);
  }
for(const auto& pdg: particlePDGs) { 
  Printf("pdg = %d", pdg); 
 }
       
  Printf("trigger size from digits = %i", trigger->size()); 
   
  for(int i = 0; i < trigger->size(); i++) {
    
    oneEventDigits.clear();
    pTgr = (o2::hmpid::Trigger*)&trigger->at(i);
		Printf("\n === Event %i", i); 
    for(int j = pTgr->getFirstEntry(); j <= pTgr->getLastEntry(); j++) {
      
      digit = (o2::hmpid::Digit)digits->at(j);
      oneEventDigits.push_back(digit);

	    std::cout << j << " evnum" << pDig->getEventNumber()<< std::endl; 
      
     }
      
    if(i==nEvent) {
       
       for(unsigned int k = 0; k < oneEventDigits.size(); k++) {
      
       pDigEvt = (o2::hmpid::Digit*)&oneEventDigits.at(k);
                     
       o2::hmpid::Digit::pad2Absolute(pDigEvt->getPadID(), &module, &padChX, &padChY);
         
       hMap[module]->Fill(padChX, padChY, pDigEvt->getQ());
                    
      }              
     }
   }
       
   c1->Divide(3,3);
   c2->Divide(3,3);
   TCanvas* c3 = new TCanvas("c3", "Particle Pdg Map", 800, 600);
   c3->Divide(3, 3); 

   TCanvas* c4 = new TCanvas("c4", "Particle Pdg Map Below", 800, 600);
   c4->Divide(3, 3);

   TCanvas* c5 = new TCanvas("c5", "Particle Pdg Map", 800, 600);
   c5->Divide(3, 3);



TCanvas* c6 = new TCanvas("c6", "PDG Histogram", 800, 600);
	c6->cd();
	hPDG->Draw();

   for(int iCh = 0; iCh<7; iCh++){
     c1->cd(pos[iCh]); 
     hMap[iCh]->Draw("Colz");
     
     c2->cd(pos[iCh]); 
     hCharge[iCh]->Draw();

     c3->cd(pos[iCh]); 
     hParticlePdgMapBelow[iCh]->Draw();
     //hParticlePdgMapAbove[iCh]->Draw("Colz");

     c4->cd(pos[iCh]); 
     hParticlePdgMapAbove[iCh]->Draw();
     //hParticlePdgMapAbove[iCh]->Draw("Colz");


     c5->cd(pos[iCh]); 
     hParticlePdgElectron[iCh]->Draw();
     hParticlePdgPion[iCh]->Draw("same");
     hParticlePdgKaon[iCh]->Draw("same");	
     hParticlePdgProton[iCh]->Draw("same");
     hParticlePdgPhoton[iCh]->Draw("same");
     hParticlePdgPhoton2[iCh]->Draw("same");
   }      
}
 TH2F *hCluChargeMap[7]; // add this
void process()
{

    // Loop over each array index to initialize the TH2F objects
    for (int i = 0; i < 7; ++i) {
        // You will want to give each TH2F object a unique name and title
        // possibly by incorporating the index `i` into the strings.

        // For thTrackMip
        thTrackMip[i] = new TH2F(Form("thTrackMip_%d", i), Form("MIP Track %d; x-axis; y-axis", i), 100, 0, 10, 100, 0, 10);
        thTrackMip[i]->SetMarkerColor(kBlue-2); // Example color
        thTrackMip[i]->SetMarkerStyle(2); // Example color
        // For thTrackPc
        thTrackPc[i] = new TH2F(Form("thTrackPc_%d", i), Form("PC Track %d; x-axis; y-axis", i), 100, 0, 10, 100, 0, 10);
        thTrackPc[i]->SetMarkerColor(kBlue); // Example color
        thTrackPc[i]->SetMarkerStyle(2); // Example color
        // For thTrackRad
        thTrackRad[i] = new TH2F(Form("thTrackRad_%d", i), Form("Rad Track %d; x-axis; y-axis", i), 100, 0, 10, 100, 0, 10);
        thTrackRad[i]->SetMarkerColor(kBlue+3); // Example color
        thTrackRad[i]->SetMarkerStyle(2); // Example color


        // For hCluChargeMap
        hCluChargeMap[i] = new TH2F(Form("hCluChargeMap_%d", i), Form("Cluster Charge Map %d; x-axis; y-axis", i), 100, 0, 10, 100, 0, 10);
    }



    // clusters and triggers
    std::vector<Cluster> *clusterArr = nullptr;
    std::vector<o2::hmpid::Topology> mTopologyFromFile, *mTopologyFromFilePtr = &mTopologyFromFile;

    std::vector<Trigger> *trigArr = nullptr;

    TTree *tCluster = HmpidDataReader::initializeClusterTree(clusterArr, trigArr, mTopologyFromFilePtr);
    // clusterArr now initialized correctly

    // MathcInfoHMP : holding trackinfo
    std::vector<o2::dataformats::MatchInfoHMP> *matchArr = nullptr;
    TTree *tMatch = HmpidDataReader::initializeMatchTree(matchArr, 0, 0, 0);

    // McTrack : holding PDG code of track
    std::vector<o2::MCTrack> *mcArr = nullptr;
    TTree *tMcTrack = HmpidDataReader::initializeMCTree(mcArr);

    int startIndexTrack = 0;
    if (trigArr == nullptr)
    {
        Printf("HmpidDataReader::initializeClusterTree trigArr== nullptr");
        return;
    }

    for (int i = 0; i < trigArr->size(); i++) // for(const auto& clusters : clustersVector) // "events loop"
    {

        auto pTgr = &trigArr->at(i);
        if (pTgr == nullptr)
        {
            Printf("pTgr== nullptr");
            continue;
        }

        const int firstEntry = pTgr->getFirstEntry();
        const int lastEntry = pTgr->getLastEntry();

        Printf("Checking trigger number %d Range Clu = %d :: %d", i, firstEntry, lastEntry);

        std::vector<Cluster> oneEventClusters;

        auto fClu = static_cast<o2::hmpid::Cluster>(clusterArr->at(firstEntry));
        auto s1Clu = static_cast<o2::hmpid::Cluster>(clusterArr->at(lastEntry - 1));
        auto sClu = static_cast<o2::hmpid::Cluster>(clusterArr->at(lastEntry));
        int eventNumber1 = fClu.getEventNumber();
        int eventNumberLast = sClu.getEventNumber();
        int eventNumberLast1 = s1Clu.getEventNumber();
        if (eventNumberLast != eventNumber1)
        {
            Printf("eventNumberLast%d != eventNumber1%d", eventNumberLast, eventNumber1);
            Printf("eventNumberLast1%d", eventNumberLast1);
        } // TODO: throw error? ef:

        for (int j = pTgr->getFirstEntry(); j <= pTgr->getLastEntry(); j++)
        {
            const auto &clu = static_cast<o2::hmpid::Cluster>(clusterArr->at(j));
            std::cout << j << " evNum " << clu.getEventNumber() << " |";
        }

        for (int j = pTgr->getFirstEntry(); j <= pTgr->getLastEntry(); j++)
        {
            const auto &clu = static_cast<o2::hmpid::Cluster>(clusterArr->at(j));

            Printf("============\n Cluster Loop \n ===============");
            if (clu.getEventNumber() != eventNumber1)
            {
                Printf("Eventnumber changed??");
                Printf("clu.getEventNumber()%d", clu.getEventNumber());
            }
            else
            {
                oneEventClusters.push_back(clu);
                std::cout << " clu " << j << " eventNum " << clu.getEventNumber();
            }
        }

        Printf("============\n Cluster Loop \n ===============");
        // find entries in tracksOneEvent which corresponds to correct eventNumber
        Printf("Reading match vector<MatchInfoHMP> for startIndexTrack %d", startIndexTrack);

        std::vector<o2::dataformats::MatchInfoHMP> *tracksOneEvent = HmpidDataReader::readMatch(tMatch, matchArr, eventNumber1, startIndexTrack);

        // get MC tracks for given event from mc;
        Printf("Reading vector<o2::MCTrack>* mcTracks for  eventNumber %d", eventNumber1);
        std::vector<o2::MCTrack> *mcTracks = HmpidDataReader::readMC(mcArr, tMcTrack, eventNumber1);

        Printf("tracksOneEvent size %d", tracksOneEvent->size());

        // for this event

        Printf("Sorting events by chamber");

        std::sort((*tracksOneEvent).begin(), (*tracksOneEvent).end(), [](const o2::dataformats::MatchInfoHMP &a, const o2::dataformats::MatchInfoHMP &b)
                  { return a.getChamber() < b.getChamber(); });

        Printf("Sorting clusteres by chamber");

        std::sort((oneEventClusters).begin(), (oneEventClusters).end(), [](const Cluster &a, const Cluster &b)
                  { return a.ch() < b.ch(); });

        std::vector<o2::dataformats::MatchInfoHMP> sortedTracks[7];
        // Assign MLinfoHMP objects to corresponding vectors based on iCh value
        for (const auto &obj : *tracksOneEvent)
        {

            const auto &iCh = obj.getChamber();
            if (iCh < 0 || iCh > 6)
            {
                std::cerr << "Warning: iCh value out of expected range: " << iCh << std::endl;
            }
            else
            {
                sortedTracks[iCh].push_back(obj);
                // sstd::cerr << "sortedTracks[iCh] " << iCh << " pushback" << std::endl;
            }
        }

        for (int i = 0; i < 7; i++)
        {
            std::cout << "Length of sortedTracks vector " << i << ": " << sortedTracks[i].size() << std::endl;
        }

        // Assuming the range of iCh values is from 0 to 6 (inclusive)
        std::vector<ClusterCandidate> sortedClusters[7];
        // Assign MLinfoHMP objects to corresponding vectors based on iCh value
        for (const auto &obj : oneEventClusters)
        {

            const auto &iCh = obj.ch();
            if (iCh >= 0 && iCh <= 6)
            {
                if (sortedTracks[iCh].size() > 0)
                {
                    // make a light copy of digits, just holding the fields charge, x, y
                    /*std::vector<ShallowDigit> shallowDigits;

                
                    //const std::vector<o2::hmpid::Cluster::Topology>& topology = obj.getClusterTopology();  // some info about digits associated w cluster*/

                    hCluChargeMap[bj.ch()]->Fill(obj.x(), obj.y(), obj.q());
                    
                    
                    
                    std::vector<std::pair<int, int>> candStatus = {{0, 0}};
                    // Printf("ClusterCandidate Ch %d", iCh);
                    ClusterCandidate temp(obj.ch(), obj.x(), obj.y(), obj.q(), obj.chi2(), obj.xe(), obj.ye(), obj.getPDG(), candStatus);
                    sortedClusters[iCh].emplace_back(temp);
                }
                else
                {
                    // std::cerr << "sortedTracks[iCh] " << iCh << " empty " << sortedTracks[iCh].size() << std::endl;
                }
            }
            else
            {
                std::cerr << "Warning: iCh value out of expected range: " << iCh << std::endl;
            }

        } // end for

        for (int i = 0; i < 7; i++)
        {

            // check if has more than one track --> this means there is no candidates
            if (sortedTracks[i].size() < 1)
            {
                Printf("sortedTracks[iCh%d].size() %d", i, sortedTracks[i].size());
                continue;
            }

            auto &clusterPerChamber = sortedClusters[i];

            std::vector<float> mipCharges;
            // fill charges of MIPs
            for (const auto &track : sortedTracks[i])
            {
                float xMip, yMip;
                int q, nph;
                track.getHMPIDmip(xMip, yMip, q, nph);
                mipCharges.emplace_back(q);
            }

            int tNum = 0;
            for (const auto &track : sortedTracks[i])
            {
                Printf("TrackNumber%d track[iCh%d].size() %d", tNum++, i, sortedTracks[i].size());
                // pass clusters (and track) by reference, and add its status per track (adding to candStatus vector )

                // for each clusterPerChamber we will have a "candidate-status" for each of the photons, this is a vector of length of sortedTracks[i].size();
                // and holds the fields

                // get MCTrack correspondign to trackId
                const auto mcTrackIndex = track.getTrackIndex();

                // find the PDG code in teh o2Kine_sim.root file by matching the mcTrackIndex for the current event ;
                const o2::MCTrack *mcTrack = HmpidDataReader::getMCEntry(mcTracks, mcTrackIndex);

                const int mcTrackPdg = mcTrack->GetPdgCode();

                const int momentum = mcTrack->GetPdgCode();

                float xRad, yRad, xPc, yPc, th, ph;
                float xMip = track.getMipX(), yMip = track.getMipY();
                track.getHMPIDtrk(xRad, yRad, xPc, yPc, th, ph);

                thTrackMip[i]->Fill(xMip, yMip);
                thTrackPc[i]->Fill(xPc, yPc);
                thTrackRad[i]->Fill(xRad, yRad);



                if(mcTrackPdg)
                thCluPdg[i]->Fill(xRad, yRad, );

                // add mcTrackPdg
                evaluateClusterTrack(clusterPerChamber, track, mipCharges, mcTrackPdg);
            }
        }
    }


    TCanvas* cCluCharge = new TCanvas("cCluCharge", "Cluster Charge MAP", 800, 600);
    TCanvas* cTrackInfo = new TCanvas("cTrackInfo", "Track impact and mip", 800, 600);

    cCluCharge->Divide(3,3);
    cTrackInfo->Divide(3,3);

    for (int i = 0; i < 7; i++) {
        int iCh = pos[i];
        cCluCharge->cd(pos[iCh]); 

        hCluChargeMap->Draw("Colz");

        cTrackInfo->cd(pos[iCh]);        
        thTrackMip[i]->Draw();
        thTrackPc[i]->Draw("same");
        thTrackRad[i]->Draw("same");


        /*

        draw maps of pdg codes for clu and track also??
        c5->cd(pos[iCh]); 
        hParticlePdgElectron[iCh]->Draw();
        hParticlePdgPion[iCh]->Draw("same");
        hParticlePdgKaon[iCh]->Draw("same");	
        hParticlePdgProton[iCh]->Draw("same");
        hParticlePdgPhoton[iCh]->Draw("same");
        hParticlePdgPhoton2[iCh]->Draw("same");
        */

    }

}
