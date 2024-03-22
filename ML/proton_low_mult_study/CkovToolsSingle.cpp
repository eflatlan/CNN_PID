

#ifndef TEST_POPULATE
#define TEST_POPULATE

//#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <random>
#include "populate.cpp"
#include "Sigma.cpp"
#include "Sigma2.cpp"
#include "populate2.cpp" // TODO: change name of class and file here 

#include <iostream>
#include <thread>
#include <chrono>

#include <vector>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <math.h>
//#include "ReconE.cpp"
//#include "ReconG.cpp"

//namespace ParticleUtils
using namespace o2;

struct Bin {
    float x;
    float y;
};




class ArrAndMap
{


  
  public :
  

  	
  	
  using vecArray2 = std::vector<std::array<double,2>>;
  
  
  vecArray2 arrMaxPionPos, arrMaxKaonPos, arrMaxProtonPos;		
  vecArray2 arrMinPionPos, arrMinKaonPos, arrMinProtonPos;		
  TVector2* errPos = new TVector2(0,0);
  	
  void setMaxArrays(const vecArray2& _arrMaxPionPos, const vecArray2& _arrMaxKaonPos, const vecArray2& _arrMaxProtonPos) {
  	arrMaxPionPos = _arrMaxPionPos;
  	arrMaxKaonPos = _arrMaxKaonPos;
  	arrMaxProtonPos = _arrMaxProtonPos;
  }	
		
 	std::unique_ptr<TH2F> limMin; 	
  std::unique_ptr<TH2F> limMax;
   void setDrawLimits(const TH2F* limMinIn, const TH2F* limMaxIn) 
   {

      
			//// limMin.reset(new TH2F(*// limMinIn));

			//// limMax.reset(new TH2F(*// limMaxIn));
   }

		

		void setErrorPos(const TVector2& posPhoton)
		{
			errPos->Set(posPhoton);
		}	
  	
  void setMinArrays(const vecArray2& _arrMinPionPos, const vecArray2& _arrMinKaonPos, const vecArray2& _arrMinProtonPos) {
    arrMinPionPos = _arrMinPionPos;
  	arrMinKaonPos = _arrMinKaonPos;
  	arrMinProtonPos = _arrMinProtonPos;
  }	
		
  
	void fillMapFromVec(TH2F* map, const vecArray2& vecArr)
	{
		for(const auto vE : vecArr) {

			//Printf("fillMapFromVec El %.2f %.2f", vE[0], vE[1]);
			map->Fill(vE[0], vE[1]);
		}
	}

  int scale = 10;
  /*
	TH2F* ckovCandMapRange = nullptr;

	TH2F* mipSizeFilter  = nullptr;

	TH2F* ckovCandMapOutRange = nullptr;

	TH2F* mipChargeFilter = nullptr;


 TH2F* hSignalAndNoiseMap =  nullptr;

	TH2F*	hSignalMIP =  nullptr;
	TH2F*	hSignalMIPpc =  nullptr;

		TH2F* arrayMaxMin[8];
	TH2F*	hMaxProton =  nullptr;
	TH2F*	hMaxPion =  nullptr;
	TH2F*	hMaxPionMinL =  nullptr;
	TH2F*	hMinPionMaxL =  nullptr;
		TH2F* hMaxKaon =  nullptr;


	TH2F*	hMinProton =  nullptr;
	TH2F*	hMinPion =  nullptr;
	TH2F*	hMinKaon = nullptr; 

  Populate2* populatePtr = nullptr; */ 
		

  void fillCkovCandMapRange(double x, double y) {
   	ckovCandMapRange.emplace_back(std::array<double,2>{x,y});
  }
  
  void fillmipSizeFilter(double x, double y) {
   	mipSizeFilter.emplace_back(std::array<double,2>{x,y});
  }

  void fillckovCandMapOutRange(double x, double y) {
   	ckovCandMapOutRange.emplace_back(std::array<double,2>{x,y});
  }
  
  void fillMipCharge(double x, double y) {
   	mipChargeFilter.emplace_back(std::array<double,2>{x,y});
  }
  

  
   


	vecArray2 mipChargeFilter, ckovCandMapRange, mipSizeFilter, ckovCandMapOutRange;
		
		/*
		~ArrAndMap() {
	hSignalMIP =  nullptr;
		hSignalMIPpc =  nullptr;
		hMaxProton =  nullptr;
		hMaxPion =  nullptr;
		hMaxPionMinL =  nullptr;
		hMinPionMaxL =  nullptr;
		hMaxKaon =  nullptr;


		hMinProton =  nullptr;
		hMinPion =  nullptr;
		hMinKaon = nullptr;

			delete ckovCandMapRange;

			delete mipSizeFilter;

		delete 	ckovCandMapOutRange;

		delete	mipChargeFilter;


delete 			hSignalAndNoiseMap;

		delete 	hSignalMIP;
		delete 	hSignalMIPpc;


			delete hMaxProton;
			delete hMaxPion;
			delete hMaxPionMinL;
		delete 	hMinPionMaxL;
			delete hMaxKaon;


			delete hMinProton;
			delete hMinPion;
			delete hMinKaon;

			delete populatePtr;
		
			delete ckovCandMapRange;

			delete mipSizeFilter;

			delete ckovCandMapOutRange;

			delete mipChargeFilter;


			delete hSignalAndNoiseMap;

			delete hSignalMIP;
			delete hSignalMIPpc;


			delete hMaxProton;
			delete hMaxPion;
			delete hMaxPionMinL;
			delete hMinPionMaxL;
			delete hMaxKaon;


			delete hMinProton;
			delete hMinPion;
			delete hMinKaon;

			delete populatePtr;
			
		}*/

		int eventCnt;
		
  std::unique_ptr<Populate2> populatePtr;


  void setPopulatePtr(std::unique_ptr<Populate2> &&_populatePtr) { 
      populatePtr = std::move(_populatePtr); 
  }
	
	int eventCount = 0;


	void setEventCount(int _eventCount) { 
		eventCount = _eventCount; 
	}	
	



  
  
	void drawTotalMap(std::vector<o2::hmpid::ClusterCandidate>& clusterTrack, int& plotNumber, float xMip, float yMip, vecArray2 pionCandidates, vecArray2 kaonCandidates, vecArray2 protonCandidates, vecArray2 canCombined, int trackPdg, vecArray2 allCand, float sigSep)
	{







		const auto trkRad = populatePtr->getTrackPos();
		const auto trkPC = populatePtr->getPcImp();

		auto xr = trkRad.X();   auto yr = trkRad.Y();
		TVector2 mip(xMip, yMip);

		auto  mipPhi = (mip-trkRad).Phi();
		auto pcPhi = (trkPC-trkRad).Phi();


		//Printf(" mip { %.2f,  %.2f}  trkPC { %.2f,  %.2f} trkRad { %.2f,  %.2f}", mip.X(), mip.Y(), trkPC.X(), trkPC.Y(), trkRad.X(), trkRad.Y());


		//Printf(" (mip-trkRad) %.2f(trkPC-trkRad)  %.2f", (mip-trkRad).Phi(), (trkPC-trkRad).Phi());
		//Printf(" (mip-trkRad) { %.2f,  %.2f}  (trkPC-trkRad) { %.2f,  %.2f}", (mip-trkRad).Phi(), (mip-trkRad).Theta(), (trkPC-trkRad).Phi(), (trkPC-trkRad).Theta());
		 auto len1 = 15.; auto len2 = 10.;


		std::unique_ptr<TLine> tLineMIP(new TLine(xr -len2 * TMath::Cos(mipPhi),yr - len2 * TMath::Sin(mipPhi),xr + len1 * TMath::Cos(mipPhi),yr + len1 * TMath::Sin(mipPhi)));

		std::unique_ptr<TLine> tLineTRK(new TLine(xr -len2 * TMath::Cos(pcPhi),yr - len2 * TMath::Sin(pcPhi),xr + len1 * TMath::Cos(pcPhi),yr + len1 * TMath::Sin(pcPhi)));



		auto distPC2MIP = (mip-trkPC).Mod();

	        int numPion = pionCandidates.size();
		int numKaon = kaonCandidates.size();
		int numProton = protonCandidates.size();
		int numTotal = canCombined.size();


		const auto trkdir = populatePtr->getTrkDir();
		auto theta = trkdir.Theta();
		  auto st = Form("\#Pi %d K %d Pr %d T %d | theta %.3f pPhi %.3f | MIP %.2f %.2f",  numPion, numKaon, numProton, numTotal, theta, pcPhi, mip.X(), mip.Y());
		  
		  std::unique_ptr<TH2F> hCkovCandMapRange(new TH2F(st, st, 1600, 0, 300, 1440, 0, 143));
		  std::unique_ptr<TH2F> hmipSizeFilter(new TH2F("mipSizeFilter", "mipSizeFilter", 1600, 0, 300, 1440, 0, 143));
		  std::unique_ptr<TH2F> hCkovCandMapOutRange(new TH2F("ckovCandMapOutRange", "ckovCandMapOutRange", 1600, 0, 300, 1440, 0, 143));
		  std::unique_ptr<TH2F> hmipChargeFilter(new TH2F("mipChargeFilter", "mipChargeFilter", 1600, 0, 300, 1440, 0, 143));
		  std::unique_ptr<TH2F> hSignalAndNoiseMap(new TH2F("Signal and Noise ", "Signal and Noise ; x [cm]; y [cm]",1600, 0, 300, 1440, 0, 1433));
		  std::unique_ptr<TH2F> hSignalMIP(new TH2F("hmip ", "hmip ; x [cm]; y [cm]", 1600, 0., 300., 1440, 0, 143));
		  std::unique_ptr<TH2F> hSignalMIPpc(new TH2F("hmip pc", "hmip pc; x [cm]; y [cm]", 1600, 0., 300., 1440, 0, 143));
		  std::unique_ptr<TH2F> hMaxProton(new TH2F("maxPoss Ckov Proton", "maxPoss Ckov Proton; x [cm]; y [cm]",1600, 0, 300, 1440, 0, 143));
		  std::unique_ptr<TH2F> hMaxPion(new TH2F("maxPoss Ckov Pion", "maxPoss Ckov Pion; x [cm]; y [cm]",1600, 0, 300, 1440, 0, 143));
		  std::unique_ptr<TH2F> hMaxPionMinL(new TH2F("maxPoss Ckov Pion min L", "maxPoss Ckov Pion min L; x [cm]; y [cm]",1600, 0, 300, 1440, 0, 143));
		  std::unique_ptr<TH2F> hMinPionMaxL(new TH2F("minPoss Ckov Pion max L", "minPoss Ckov Pion max L; x [cm]; y [cm]",1600, 0, 300, 1440, 0, 143));
		  std::unique_ptr<TH2F> hMaxKaon(new TH2F("maxPoss Ckov Kaon", "maxPoss Ckov Kaon; x [cm]; y [cm]",1600, 0, 300, 1440, 0, 143));
		  std::unique_ptr<TH2F> hMinProton(new TH2F("minPoss Ckov Proton", "minPoss Ckov Proton; x [cm]; y [cm]",1600, 0, 300, 1440, 0, 143));
		  std::unique_ptr<TH2F> hMinPion(new TH2F("minPoss Ckov Pion", "minPoss Ckov Pion; x [cm]; y [cm]",1600, 0, 300, 1440, 0, 143));
		  std::unique_ptr<TH2F> hMinKaon(new TH2F("minPoss Ckov Kaon", "minPoss Ckov Kaon; x [cm]; y [cm]",1600, 0, 300, 1440, 0, 143));
		  // ... similarly for other TH2F objects you might have


		  std::unique_ptr<TH2F> herrPos(new TH2F("errPos", "errPos; x [cm]; y [cm]",1600, 0, 300, 1440, 0, 143));




			//printf("		fillMapFromVec(mArrAndMap->hMaxPion, arrMaxPion);// map, array");
			fillMapFromVec(hMaxPion.get(), arrMaxPionPos);// map, array
			fillMapFromVec(hMaxKaon.get(), arrMaxKaonPos);// map, array
			fillMapFromVec(hMaxProton.get(), arrMaxProtonPos);// map, array

			fillMapFromVec(hMinPion.get(), arrMinPionPos);// map, array
			fillMapFromVec(hMinKaon.get(), arrMinKaonPos);// map, array
			fillMapFromVec(hMinProton.get(), arrMinProtonPos);
			



		//// limMin->SetMarkerStyle(2);
		//// limMax->SetMarkerStyle(2);

		// limMin->SetMarkerColor(kOrange-3);
		// limMax->SetMarkerColor(kOrange+2);

		// limMin->SetMarkerColor(kCyan);    


		hCkovCandMapRange->SetMarkerColor(kGreen );
		hCkovCandMapOutRange->SetMarkerColor(kGreen + 2);

		//herrPos->SetMarkerColor(kRed+3);
		//herrPos->SetMarkerStyle(3);
		
		
		//herrPos->Fill(errPos->X(), errPos->Y());
	

		hmipSizeFilter->SetMarkerColor(kOrange-1);
		hmipChargeFilter->SetMarkerColor(kOrange+1);

		hCkovCandMapRange->SetMarkerStyle(3);
		hCkovCandMapOutRange->SetMarkerStyle(2);

		hmipSizeFilter->SetMarkerStyle(2);
		hmipChargeFilter->SetMarkerStyle(2);

		int xrange = 100;
		int yrange = 100;

		int xMin = xMip - xrange;
		int xMax = xMip + xrange;
		int yMin = yMip - yrange;
		int yMax = yMip + yrange;

		//xMin = 0;
		//xMax = 160*0.8;
		//yMin = 0;
		//yMax = 144*0.84;


auto st2 = Form("All clusters pdg %d ; x [cm]; y [cm]",trackPdg);
		std::unique_ptr<TH2F> totalCluMap(new TH2F("All clusters", st2,320, 0, 159, 288, 0, 143));



		totalCluMap->SetAxisRange(xMin, xMax, "X");
		totalCluMap->SetAxisRange(yMin, yMax, "Y");
	  totalCluMap->SetMarkerStyle(2);


		for(const auto& c : clusterTrack) {
			totalCluMap->Fill(c.mX, c.mY, c.mQ);
		}


		hCkovCandMapRange->SetAxisRange(xMin, xMax, "X");
		hCkovCandMapRange->SetAxisRange(yMin, yMax, "Y");

		hmipSizeFilter->SetAxisRange(xMin, xMax, "X");
		hmipSizeFilter->SetAxisRange(yMin, yMax, "Y");

		hCkovCandMapOutRange->SetAxisRange(xMin, xMax, "X");
		hCkovCandMapOutRange->SetAxisRange(yMin, yMax, "Y");

		hmipChargeFilter->SetAxisRange(xMin, xMax, "X");
		hmipChargeFilter->SetAxisRange(yMin, yMax, "Y");

		hSignalAndNoiseMap->SetAxisRange(xMin, xMax, "X");
		hSignalAndNoiseMap->SetAxisRange(yMin, yMax, "Y");

		hSignalMIP->SetAxisRange(xMin, xMax, "X");
		hSignalMIP->SetAxisRange(yMin, yMax, "Y");

		hSignalMIPpc->SetAxisRange(xMin, xMax, "X");
		hSignalMIPpc->SetAxisRange(yMin, yMax, "Y");

		hMaxProton->SetAxisRange(xMin, xMax, "X");
		hMaxProton->SetAxisRange(yMin, yMax, "Y");

		hMaxPion->SetAxisRange(xMin, xMax, "X");
		hMaxPion->SetAxisRange(yMin, yMax, "Y");

		hMaxPionMinL->SetAxisRange(xMin, xMax, "X");
		hMaxPionMinL->SetAxisRange(yMin, yMax, "Y");

		hMinPionMaxL->SetAxisRange(xMin, xMax, "X");
		hMinPionMaxL->SetAxisRange(yMin, yMax, "Y");

		hMaxKaon->SetAxisRange(xMin, xMax, "X");
		hMaxKaon->SetAxisRange(yMin, yMax, "Y");

		hMinProton->SetAxisRange(xMin, xMax, "X");
		hMinProton->SetAxisRange(yMin, yMax, "Y");

		hMinPion->SetAxisRange(xMin, xMax, "X");
		hMinPion->SetAxisRange(yMin, yMax, "Y");

		hMinKaon->SetAxisRange(xMin, xMax, "X");
		hMinKaon->SetAxisRange(yMin, yMax, "Y");

		
		
		fillMapFromVec(hCkovCandMapRange.get(), ckovCandMapRange);// map, array
		fillMapFromVec(hCkovCandMapOutRange.get(), ckovCandMapOutRange);// map, array
		fillMapFromVec(hmipSizeFilter.get(), mipSizeFilter);// map, array
		fillMapFromVec(hmipChargeFilter.get(), mipChargeFilter);


		hMaxProton->SetMarkerColor(kGreen+4);
		hMaxKaon->SetMarkerColor(kRed +1);
		hMinProton->SetMarkerColor(kGreen+3);
		hMinKaon->SetMarkerColor(kRed);
		hMaxPion->SetMarkerColor(kBlue + 4);    // max ckov max L
		hMinPion->SetMarkerColor(kBlue +3 ); 				// min ckov min L/**/
		hMinPionMaxL->SetMarkerColor(kBlue); 		// min ckov max L
		hMaxPionMinL->SetMarkerColor(kBlue + 4);// max ckov min L
		/*
		hMaxPion->SetMarkerStyle(3);
		hMinPionMaxL->SetMarkerStyle(3);*/
		hSignalMIPpc->SetMarkerStyle(3);
		hSignalMIP->SetMarkerStyle(3);
		hSignalMIPpc->SetMarkerStyle(3);
		hSignalMIPpc->SetMarkerColor(kRed);
		hSignalAndNoiseMap->SetMarkerStyle(2);


		auto trkPCMap = std::make_unique<TH2F>("trkPCMap ", "trkPCMap; x [cm]; y [cm]",160*10,0.,300.,144*10,0,143);
		auto trkRadMap = std::make_unique<TH2F>("trkRadMap ", "trkRadMap; x [cm]; y [cm]",160*10,0.,300.,144*10,0,143);


    auto MIP = std::make_unique<TH2F>("MIP ", "MIP; x [cm]; y [cm]",160*10,0.,300.,144*10,0,143);


		trkRadMap->Fill(trkRad.X(), trkRad.Y());
		trkPCMap->Fill(trkPC.X(), trkPC.Y());

		MIP->Fill(xMip, yMip);
		MIP->SetMarkerColor(kCyan);    // max ckov max L
		MIP->SetMarkerStyle(3);    // max ckov max L



		trkRadMap->SetMarkerStyle(3); 
		trkPCMap->SetMarkerStyle(2); 
		tLineMIP->SetLineColor(kCyan); 
		tLineTRK->SetLineColor(kRed); 



		auto tcnvRane = std::make_unique<TCanvas>(Form("tcnvRane%d", plotNumber), Form("tcnvRane%d", plotNumber), 1600, 800);
	tcnvRane->Divide(2);


		tcnvRane->cd(1);


		hMaxProton->SetMarkerColor(kGreen+4);
		hMaxKaon->SetMarkerColor(kRed +1);
		hMinProton->SetMarkerColor(kGreen+3);
		hMinKaon->SetMarkerColor(kRed);
		hMaxPion->SetMarkerColor(kBlue + 4);    // max ckov max L
		hMinPion->SetMarkerColor(kBlue +3 ); 				// min ckov min L/**/
		hMinPionMaxL->SetMarkerColor(kBlue); 		// min ckov max L
		hMaxPionMinL->SetMarkerColor(kBlue + 4);// max ckov min L


	TLegend *legend = new TLegend(0.758,0.758, 0.975,0.925); // 

	// Add entries to the legend
	legend->AddEntry(trkRadMap.get(), "rad", "p"); //
	legend->AddEntry(trkPCMap.get(), "MIP", "p");
	legend->AddEntry(hCkovCandMapOutRange.get(), "Out of region", "p");
	legend->AddEntry(hmipSizeFilter.get(), "Size > 2 && Charge > 200", "p");
	legend->AddEntry(hCkovCandMapRange.get(), "Cluster Candidates", "p");


	legend->AddEntry(hMaxProton.get(), "hMaxProton", "p"); //
	legend->AddEntry(hMaxKaon.get(), "hMaxKaon", "p");
	legend->AddEntry(hMinProton.get(), "OhMinProton", "p");



		hCkovCandMapRange->Draw();
	  hCkovCandMapRange->SetStats(kFALSE);
		hCkovCandMapOutRange->Draw("same");
		hmipSizeFilter->Draw("same");
		hmipChargeFilter->Draw("same");
		//MIP->Draw("same");
legend->Draw("same");

		hMaxPion->Draw("same");     //  Printf("makeEvent()  hMaxPion->Draw");
		hMinPion->Draw("same");
		hMaxProton->Draw("same");   //Printf("makeEvent()  hMaxProton->Draw");
		hMinProton->Draw("same");
		hMinKaon->Draw("same");
		hMaxKaon->Draw("same");
		tLineTRK->Draw("same");
		// tLineMIP->Draw("same");



		trkRadMap->SetMarkerStyle(20);  // E.g., a filled circle
		trkRadMap->SetMarkerColor(kRed); // E.g., a red color
 		trkRadMap->SetContour(99); // Use a high number for detailed contours
     trkRadMap->SetLineWidth(4); ///
		trkRadMap->SetLineColor(kRed); // Setting a distinct color for the contour lines
		trkRadMap->Draw("CONT3 same"); 
		trkPCMap->Draw("same");
		//herrPos->Draw("same");
						//printf("opulatePtr->// limMin->GetEntries() num BinEntries = %f ", populatePtr->// limMin->GetEntries());
			//printf("end// limMax num BinEntries = %f ", // limMax->GetEntries());
			//printf("emd// limMin num BinEntries = %f ", // limMin->GetEntries());



    auto textNumPion = new TLatex(xMin, yMin + 5, ("Pions: " + std::to_string(numPion)).c_str()); 
    textNumPion->SetTextSize(0.04);
   //  textNumPion->Draw("same");

    auto textNumKaon = new TLatex(xMin + 10, yMin +5, (" K: " + std::to_string(numKaon)).c_str()); 
    textNumKaon->SetTextSize(0.04);
   //  textNumKaon->Draw("same");

    auto textNumProton = new TLatex(xMin + 20, yMin +5, (" Pr: " + std::to_string(numProton)).c_str()); 
    textNumProton->SetTextSize(0.04);
   //  textNumProton->Draw("same");

    auto textNumTotal = new TLatex(xMin + 30, yMin +5, ("To : " + std::to_string(numTotal)).c_str()); 
    textNumTotal->SetTextSize(0.04);
    // textNumTotal->Draw("same");





    TPad* pad2 = static_cast<TPad*>(tcnvRane->cd(2));



    pad2->SetRightMargin(.055+pad2->GetRightMargin());
    pad2->SetLeftMargin(-.055+pad2->GetLeftMargin());
    pad2->SetLogz(1);
	  totalCluMap->SetStats(kFALSE);
		totalCluMap->Draw("Colz");
  //	pad2->Clear();
     trkRadMap->SetLineWidth(4); ///
		trkRadMap->SetMarkerStyle(20);  // E.g., a filled circle
		trkRadMap->SetMarkerColor(kRed); // E.g., a red color
		trkRadMap->SetContour(99); // Use a high number for detailed contours
		trkRadMap->SetLineColor(kRed); // Setting a distinct color for the contour lines
		trkRadMap->Draw("CONT3 same"); // Drawing the contours on the same pad



    tcnvRane->SaveAs(Form("Segmented_%dCkov%.2f.png", plotNumber, sigSep));
    // limMin->Draw("same");
    // limMax->Draw("same");

    //tcnvRane->SaveAs(Form("Segmented%d.png", plotNumber));
    plotNumber++;






	}
  /*
	void drawMaxRegions()
	{


		const auto trkPC = populatePtr->getPcImp();
		const auto trkRad = populatePtr->getTrackPos();
		TH2F* trkPCMap = new TH2F("trkPCMap ", "trkPCMap; x [cm]; y [cm]",160*20,0.,300.,144*20,0,143);
		TH2F* trkRadMap = new TH2F("trkRadMap ", "trkRadMap; x [cm]; y [cm]",160*20,0.,300.,144*20,0,143);
	  trkRadMap->Fill(trkRad.X(), trkRad.Y());
	  trkPCMap->Fill(trkPC.X(), trkPC.Y());


		trkRadMap->SetMarkerStyle(3);  trkRadMap->SetMarkerColor(kGreen+4);
		trkPCMap->SetMarkerStyle(3);  trkPCMap->SetMarkerColor(kGreen+2);

		TCanvas *thSignalAndNoiseMap = new TCanvas(Form("hSignalAndNoiseMap%d", eventCnt),Form("hSignalAndNoiseMap%d", eventCnt),800,800);  
		thSignalAndNoiseMap->cd();
		hSignalAndNoiseMap->Draw();


		hMaxPion->Draw("same");     //  Printf("makeEvent()  hMaxPion->Draw");
		hMinPion->Draw("same");
		hMaxProton->Draw("same");   //Printf("makeEvent()  hMaxProton->Draw");
		hMinProton->Draw("same");
		hMinKaon->Draw("same");
		hMaxKaon->Draw("same");
		
		hSignalMIP->Draw("same");
		hSignalMIPpc->Draw("same");
		
		trkRadMap->Draw("same");
		trkPCMap->Draw("same");
	}



	void drawTotalMapAndMaxRegions()
	{
		TCanvas* tcnvRane = new TCanvas(Form("TotalMap%d", eventCnt), Form("TotalMap%d", eventCnt), 1600, 800);

				const auto trkPC = populatePtr->getPcImp();
		const auto trkRad = populatePtr->getTrackPos();
		TH2F* trkPCMap = new TH2F("trkPCMap ", "trkPCMap; x [cm]; y [cm]",160*20,0.,159.,144*20,0,143);
		TH2F* trkRadMap = new TH2F("trkRadMap ", "trkRadMap; x [cm]; y [cm]",160*20,0.,159.,144*20,0,143);
	  trkRadMap->Fill(trkRad.X(), trkRad.Y());
	  trkPCMap->Fill(trkPC.X(), trkPC.Y());


		trkRadMap->SetMarkerStyle(3);  trkRadMap->SetMarkerColor(kGreen+4);
		trkPCMap->SetMarkerStyle(3);  trkPCMap->SetMarkerColor(kGreen+2);

		tcnvRane->cd();
		ckovCandMapRange->Draw();
		ckovCandMapOutRange->Draw("same");
		mipSizeFilter->Draw("same");
		mipChargeFilter->Draw("same");

		hSignalAndNoiseMap->Draw();


		hMaxPion->Draw("same");     //  Printf("makeEvent()  hMaxPion->Draw");
		hMinPion->Draw("same");
		hMaxProton->Draw("same");   //Printf("makeEvent()  hMaxProton->Draw");
		hMinProton->Draw("same");
		hMinKaon->Draw("same");
		hMaxKaon->Draw("same");
		
		hSignalMIP->Draw("same");
		hSignalMIPpc->Draw("same");
				trkRadMap->Draw("same");
		trkPCMap->Draw("same");
	} */ 

};

class CkovTools {


private: 

  
	// for storing candidates...
	struct Candidate {
		double x, y = 0.;
		double R = 0.;
		double phiL = 0., phi = 0.;
		bool isCandidate = false;
	};


  bool print = false; // ef: TODO : later pass this in ctor 


  std::unique_ptr<Populate2> populatePtr, populatePtrCp, populate2Ptr, populatePtrOuter, populatePtrInner;


	// using array = std::array;
	using vecArray4 = std::vector<std::array<double,4>>;
	using vecArray3 = std::vector<std::array<double,3>>;
	using vecArray2 = std::vector<std::array<double,2>>;

	using vecArrayPair = std::vector<std::pair<double,double>>;


	// using arrArray3 = std::vector<std::array<double,3>>;


	using segType = std::vector<std::pair<double, double>>;
	segType segPionLocal = {{}};

	bool kaonStatus = true, pionStatus = true, protonStatus = true;
	//TLine* tlinePion;


        float mSigmaSep = 1.5;
	// used in SigCrom
	//  double f = 0.00928*(7.75-5.635)/TMath::Sqrt(12.);
  // static constexpr double f = 0.0172*(7.75-5.635)/TMath::Sqrt(24.);
  static constexpr double sq6 = 2.44948974278;
  static constexpr double f = 0.0172*(7.75-5.635)/(2 * sq6);

	static constexpr double PI = M_PI;
	static constexpr double halfPI = M_PI/2;
	static constexpr double twoPI = M_PI*2;


	static constexpr double stdDevPion = 0.001; 
	static constexpr double stdDevKaon = 0.001; 
	static constexpr double stdDevProton = 0.001;


	static constexpr float tGap = 8;
	static constexpr float  rW = 1.5; // was 1?, 

	static constexpr float  qW = 0.5;
	static constexpr float lMax = 1.5;


	static constexpr float CH4GapWidth = 8;
	static constexpr float  RadiatorWidth = 1.5; // was 1?
	static constexpr float  QuartzWindowWidth = 0.5;

	static constexpr float  L_CONST = rW/2;



	static constexpr double dnDE = 0.0172; // σE = (dn/dE)σdet dn/dE is 0.0172 eV−1. 
																				 // σdet represents the standard deviation of the detected Cherenkov photon spectrum


  // for simulation :  mean photonergny assumed to be cosntant
	static constexpr double photEnergyMean = 6.82/1000.; // jsut based on some runs  
	
	

  // 1: chromacity error
	// std-dev of single-photon angular resolution, contribution from fluctuation in photon-energy
	static constexpr double sigmaE = dnDE * photEnergyMean; // given in Rad
	


  // 4 : 
  // The track incidence angle error, related to the particle angle θp and to the precision of the tracking
  // devices. In the following, the θp error, assumed to be of the order of 2 mrad at the considered incidence angles, will not be quoted in tables and plots but simply included in the calculation of the
  // total angular resolution.
	static constexpr double sigmaThetaP = 0.002; // given in Rad



	// combined angular resolution 
	static constexpr float sigmaAnglesSq = sigmaE * sigmaE + sigmaThetaP * sigmaThetaP;
	const float sigmaAngles = TMath::Sqrt(sigmaE * sigmaE + sigmaThetaP * sigmaThetaP);



  // (3) The localization error, related to the precision with which the photon and particle impact coordinates can be measured. It is determined by the detector characteristics (pad size, sense wire pitch)
  // and by the photon feedback.

	// error of position of clusters 
	static constexpr float sigmaR = .2;    // TDR : sizeY / sqrt(12)
	static constexpr float sigmaRsq = .2 * .2;    // TDR : sizeY / sqrt(12)



  // (2) The geometric error, related to the spread of the emission point along the particle path in the Cherenkov radiator.
  // It depends on the ratio RW/GAP between the radiator thickness, RW, and he proximity gap width, GAP; 
  // it can be minimized by increasing GAP and reducing RW, provided the number of photoelectrons per ring is sufficient for pattern recognition (Chapter 4).
 const float  sigmaLfactor = rW/TMath::Sqrt(12);    // sigmaL = sigmaLfactor/ cosThetaP
 float sigmaL = 0.0;
 float sigmaLsq = 0.0;


	// ef : set this constexpr ins
	float L = rW/2;


  // L value for reconstruction
  static constexpr float  EmissionLenght = RadiatorWidth/2;

  float thetaP, phiP, xPC, yPC, xRad, yRad; 
 float nF, nQ, nG;  
 std::array<float, 3> ckovHypsMin, ckovHypsMax;
 std::vector<std::pair<double, double>> photons;



 double ckovProton = 999, ckovKaon = 999, ckovPion = 999;


 // instantiate ckovhyp boundaries to "invalid values" 
 // --> later set them if they exceed momentum threshold
 double ckovPionMin = 999, ckovPionMax = 0, ckovKaonMin = 999, ckovKaonMax = 0, ckovProtonMin = 999,ckovProtonMax = 0;

 // double mRMax, mL1Max, mL2Max, mRMin, mL1Min, mL2Min;

 double momentum, mass;

 float cosThetaP, sinThetaP, tanThetaP;
 float cosPhiP, sinPhiP, tanPhiP;
 float trackCkov;
 double xMipLocal, yMipLocal;
 float phiLCurrent, etaCCurrent;


 float xMip, yMip, qMip;

TVector2 trkPos; // trk at RAD
TVector3 trkDir; // trk mag theta phi
TVector2 mipPos; // MIP PC


 //std::unique_ptr<ArrAndMap> mArrAndMap;// = nullptr;// new ArrAndMap(eventCnt);
 int trackPdg; string trackPdgString;
 int eventCnt;
 TVector2 trkPC;

public:
typedef std::vector<std::pair<double, double>> MapType;

void addPhoton(double phiL, double etaC)
{
  photons.push_back(std::make_pair(phiL, etaC));
}


// set new values of phiL and etaC for new photon
void setPhoton(double phiL, double etaC)
{
  phiLCurrent = phiL;
  etaCCurrent = etaC;
}


//   double pc[5] = {xRad,yRad,L,thetaP, phiP};
double refIndexes[3] = {nF, nQ, nG};


~CkovTools() {
}


// ef: let trackCkov be empty aon; TODO: add this based on pdg!


CkovTools (double radParams[7], double refIndexes[3], double MIP[3],
           std::array<float, 3> ckovHypsMin,  std::array<float, 3> ckovHypsMax, float trackCkov, int eventCnt, int _trackPdg, float sigmaSep)
  : 
    ckovHypsMin(ckovHypsMin),  ckovHypsMax(ckovHypsMax), trackCkov(trackCkov), eventCnt(eventCnt) { 
    
      
  mSigmaSep = sigmaSep;
  trackPdg = _trackPdg;
  
  trackPdgString = getPDG(trackPdg);
  // double radParams[6] = {xRad,yRad,L,thetaP, phiP, randomValue.momentum};
  
  xMip = MIP[0], yMip = MIP[1], qMip = MIP[2]; 


  
  xRad= radParams[0];
  yRad= radParams[1];
  L = radParams[2]; 
  thetaP = radParams[3];
  phiP = radParams[4];
  momentum = radParams[5];	 
  mass = radParams[6]; // ef; let this be empty aon! TODO: get this from pdg

  nF = refIndexes[0];
	trkPos.Set(xRad, yRad); 								 // track positon in LORS at RAD   // XY mag

  const double winThick = 0.5, radThick = 1.5; const int gapThick = 8;
  const double gapIdx = 1.0005, winIdx = 1.583; // inIdx = 1.5787;

  double zRad = -0.5 * radThick - 0.5 * winThick;     // z position of middle of RAD
  TVector3 rad(trkPos.X(), trkPos.Y(), zRad);                           // impact point at middle of RAD
  TVector3 pc(xMip, yMip, 0.5 * winThick + gapThick); // mip at PC
  
  

 // Printf("Phi  %.2f Thta %.2f of Track", phiP, thetaP);
  phiP = (pc-rad).Phi();
  thetaP = (pc-rad).Theta();

  Printf("Phi  %.2f Thta %.2f of rad--MIP", phiP, thetaP);
  
  trkPC.Set(xMip, yMip); // MIP pos at PC
        

  mipPos.Set(xMip, yMip); // MIP pos at PC

	 

	trkDir; 
	trkDir.SetMagThetaPhi(1, thetaP, phiP);  // track direction in LORS at RAD

  nF = 1.2928  - 0.0025; // ef got this from 1 run .. assuming T = 20 for sim
  auto nFstd2 = 0.01; // 2 times std-dev for tejh run


  populatePtrInner.reset(new Populate2(trkPos, trkDir, nF));
  populatePtrOuter.reset(new Populate2(trkPos, trkDir, nF ));
  populatePtr.reset(new Populate2(trkPos, trkDir, nF));
  populatePtrCp.reset(new Populate2(trkPos, trkDir, nF));
	populate2Ptr.reset(new Populate2(trkPos, trkDir, nF));

	populate2Ptr->setPcImp(trkPC);
	populatePtrCp->setPcImp(trkPC);	
	populatePtr->setPcImp(trkPC);	
	populatePtrOuter->setPcImp(trkPC);
	populatePtrInner->setPcImp(trkPC);  
  
	// set tehse to be constant?
  nQ = 1.583; // inIdx = 1.5787;
  nG = 1.005;
 
 
   trackPdgString = getPDG(trackPdg);
  //printf(" Track PDG %d %s | CkovTools momentum = %.2f, refFreon = %.2f; ckovHyps : %.2f %.2f %.2f", trackPdg, trackPdgString.c_str(),momentum, nF, ckovHyps[0], ckovHyps[1], ckovHyps[2]);




  // TODO ; later check how to deal w this: set it to highest value in range that makes it valid ? 
	if(TMath::IsNaN(ckovHypsMax[0])){
 	  //printf("Pion CkovHyps is Nan!");
	  setPionStatus(false);
	  // setIsNan()?
	} else {
		setPionStatus(true);
 	  //printf("Pion CkovHyps %.2f", ckovHyps[0]);
  }

  if(TMath::IsNaN(ckovHypsMax[1])){
	  //printf("Kaon CkovHyps is Nan!");
   	setKaonStatus(false);
	  // setIsNan()?
	} 
	else { 
 	  //printf("Kaon CkovHyps %.2f", ckovHyps[1]);
  }

  if(TMath::IsNaN(ckovHypsMax[2])){
  	  //printf("Proton CkovHyps is Nan!");
	  	setProtonStatus(false);
	} else {
 	  //printf("Proton CkovHyps %.2f", ckovHyps[2]);
	}

	if(getPionStatus()){
	  ckovPionMin = ckovHypsMin[0] - sigmaAngles;
	  ckovPionMax = ckovHypsMax[0] + sigmaAngles;
    ckovPion =  ckovHypsMax[0];
 	  //printf("init CkovTools constructor : getPionStatus() true ! minPion %.2f, maxPion %.2f ", ckovPionMin, ckovPionMax);
  }	else {
 	  //printf("init CkovTools constructor : getPionStatus() was false !");
  }
  
  if(getKaonStatus()){
	  ckovKaonMin = ckovHypsMin[1] - sigmaAngles;
  	ckovKaonMax = ckovHypsMax[1] + sigmaAngles;
    ckovKaon =  ckovHypsMax[1];
  } if(getProtonStatus()){
		ckovProtonMin = ckovHypsMin[2] - sigmaAngles;
		ckovProtonMax = ckovHypsMax[2] + sigmaAngles;
    ckovProton =  ckovHypsMax[2];
 }

	cosThetaP = TMath::Cos(thetaP);
	sinThetaP = TMath::Sin(thetaP);
	tanThetaP = TMath::Tan(thetaP);

	cosPhiP = TMath::Cos(phiP);
	sinPhiP = TMath::Sin(phiP);




  // set the geometric error std-dev 
  sigmaL = sigmaLfactor/ cosThetaP;
  sigmaLsq = sigmaL * sigmaL;

	TRotation rotZ; rotZ.RotateZ(phiP);
	TRotation rotY; rotY.RotateY(thetaP);

	TVector3 ip(0,0,rW-L+tGap+qW);
	TVector3 op; op = rotZ*rotY*ip;

	xMipLocal =  tanThetaP*cosPhiP*(rW-L + tGap + qW);
	yMipLocal =  tanThetaP*sinPhiP*(rW-L + tGap + qW);

	auto dX = tanThetaP*cosPhiP*(rW-L + tGap + qW);
	auto dY = tanThetaP*sinPhiP*(rW-L + tGap + qW);


	xPC = dX + xRad; // bruke PC eller MIP?
	yPC = dY + yRad;


}



void setPionStatus(bool status)
{ 
  pionStatus = status;
}

bool getPionStatus() const
{
  return pionStatus;
}

void setKaonStatus(bool status)
{ 
  kaonStatus = status;
}

bool getKaonStatus() const
{
  return kaonStatus;
}

void setProtonStatus(bool status)
{ 
  protonStatus = status;
}

bool getProtonStatus() const
{
  return protonStatus;
}




double getCkovPion()
{
  return ckovPion;
}

double getCkovKaon()
{
  return ckovKaon;
}

double getCkovProton()
{
  return ckovProton;
}


double getMinCkovKaon()
{
  return ckovKaonMin;
}

double getMaxCkovKaon()
{
  return ckovKaonMax;
}

double getMinCkovPion()
{
  return ckovPionMin;
}

double getMaxCkovPion()
{
  return ckovPionMax;
}

double getMinCkovProton()
{
  return ckovProtonMin;
}

double getMaxCkovProton()
{
  return ckovProtonMax;
}


double getPhiP()
{
  return phiP;
}

double getThetaP()
{
  return thetaP;
}

// TODO: add nQ
// only consider photons in the correct range:


// TODO: add nQ
// only consider photons in the correct range:


// x, y, etaC


// std::array<int, 4> arrayInfo;

	// ;

//ckovTools.segment(clusterPerChamber, arrayInfo, track.getTrackIndex(), mipCharges, xMip, yMip, q /*MIP-charge*/, mcTrackPdg, track);


std::vector<std::pair<double, double>> segment(std::vector<o2::hmpid::ClusterCandidate>& clusterTrack, std::array<int, 4>& arrayInfo, int trackIndex, const std::vector<float>& mipCharges, float mipX, float mipY, float mipCharge, const int mcTrackPdg, const o2::dataformats::MatchInfoHMP& track, int trackNumber,  int& plotNumber)
{ 

  

  vecArray2 pionCandidates, kaonCandidates, protonCandidates, canCombined, allCand;
  
  
  ArrAndMap mArrAndMap;// new ArrAndMap(eventCnt);
  
  //printf("ckovTools enter  ckovTools.segment"); 	 
  //const auto infString = Form("localRef #Theta_{p}  = %.4f #Phi_{p} = %.4f L = %.2f \n #Theta_{C} = %.4f maxCkov = %.4f ; x [cm]; y [cm]", thetaP,phiP, L,trackCkov,ckovPionMax); 

  //TH2F *localRefMIP = new TH2F("localRefMIP ", infString,800,-40.,-40.,900,-40.,40.);

  
  /*
  TVector2 trkPos(xRad, yRad);
  TVector3 trkDir; 
  trkDir.SetMagThetaPhi(1, thetaP, phiP);
  Populate* populatePtr = new Populate(trkPos, trkDir, nF);
  Populate populate(trkPos, trkDir, nF);*/


  // void setArrayMin(Populate* populate, double etaTRS, vecArray3 inPutVector) 


	const size_t kN = 200;

	//array<array<double, 3> ,kN> arrMaxPion; 

  

	// disse kan brukes istedet for maxKaonVec...?
	vecArray4 arrMaxPion, arrMinPion, arrMinProton, arrMaxProton, arrMinKaon, arrMaxKaon;

	vecArray2 arrMaxPionPos, arrMinPionPos, arrMinProtonPos, arrMaxProtonPos, arrMinKaonPos, arrMaxKaonPos;

	// do not reserve size? NO prob just dont resize :) 
	//arrMaxPion.reserve(kN);
	//arrMinPion.reserve(kN);
	//arrMaxPion.resize(kN);
	//arrMinPion.resize(kN);

	arrMaxPionPos.reserve(kN);
	arrMinPionPos.reserve(kN);


	// arrMaxPion.resize(kN);
		
	// check if candidate can be proton (i.e., that it exceeds momentum threshold)
	if(getProtonStatus()) {
		//arrMaxProton.reserve(kN);
		//arrMinProton.reserve(kN);
		//arrMaxProton.resize(kN);
		//arrMinProton.resize(kN);

		arrMinProtonPos.reserve(kN);
		arrMaxProtonPos.reserve(kN);

		Sigma2 	sigmaProton(getCkovProton(), phiP, thetaP, nF);
	  //Printf(" getCkovProton() %.4f", getCkovProton());
    calculateDifference(arrMaxProtonPos, arrMinProtonPos, getCkovProton(), sigmaProton); 
	}
	

	// check if candidate can be Kaon (i.e., that it exceeds momentum threshold)
	if(getKaonStatus()) {
		//arrMaxKaon.reserve(kN);
		//arrMinKaon.reserve(kN);				

		//arrMaxKaon.resize(kN);
		//arrMinKaon.resize(kN);	

		arrMaxKaonPos.reserve(kN);
		arrMinKaonPos.reserve(kN);				


	Printf(" getCkovKaon() %.4f", getCkovKaon());
		Sigma2 	sigmaKaon(getCkovKaon(), phiP, thetaP, nF); //(double thetaC, double phiP, double thetaP, double _nF)  
		calculateDifference(arrMaxKaonPos, arrMinKaonPos, getCkovKaon(), sigmaKaon); 
	}

	


	//printf("BF : Length of elem vectors : arrMaxPion %zu", arrMaxPion.size());	
  //printf("calling setArrayMax w getMaxCkovPion() = %.2f", getMaxCkovPion());

	Printf("BF : Length of elem vectors : arrMinPionPos %zu", arrMinPionPos.size());
  // void calculateDifference(vecArray2& maxVec, vecArray2& minVec, double etaTrsMax, double etaTrsMin) {

  Sigma2 	sigmaPion(getCkovPion(), phiP, thetaP, nF); //(double thetaC, double phiP, double thetaP, double _nF)  


	Printf(" getCkovPion() %.4f", getCkovPion());
  calculateDifference(arrMaxPionPos, arrMinPionPos, getCkovPion(), sigmaPion); 

	Printf("AFTER : Length of elem vectors : arrMinPionPos %zu", arrMinPionPos.size());
	// fill all the values in the maps 
	
	//Printf("Length of elem vectors : arrMinPion %zu", arrMinPion.size());
	//Printf("Length of elem vectors : arrMaxPion %zu", arrMaxPion.size());
	
  // get track impacrt point at PC


	const int scale = 2;



  Printf("dX %f dY %f ", xMipLocal, yMipLocal);

  int numPhotons= 0;


  const auto area = 144*156; 
 
 
  MapType filledBins;
 
  std::vector<std::pair<double, double>> photonCandidates;


 
  double xML = xMipLocal, yML = yMipLocal;

  //int iPhotCount = 0;

  int iPhotCount = 0;
  //for(const auto& photons : ckovAndBackground) {


  //LOGP(info, "CkovTools : number of photons = {} ", clusterTrack.size());

  int rOverMax = 0;
  
  int photNum = 0;

  int cStatusCkov = 0, cStatus = 0;
  
  TVector2 errorPos;
  bool drawMap = false;
  for(auto& photons : clusterTrack) 
  {

    cStatusCkov = 0;

		auto pdgString = getPDG(photons.mPDG);


		trackPdgString = getPDG(trackPdg);

  
    //printf("\n\n =========== Photon num %d of %zu PDG : %d %s ========================== \n", photNum++, clusterTrack.size(), photons.mPDG, pdgString.c_str());	
    
    
    
    photons.setCandidateStatus(0);
        photons.setCandidateStatusCkov(0);
    const auto& x = photons.mX, y = photons.mY;
    allCand.push_back(std::array<double,2>{x,y});

    const auto& dist = (photons.mX - mipX)*(photons.mX - mipX) + (photons.mY - mipY)*(photons.mY - mipY);
    

    // match track PDG w MIP cluster : 

    // ef : MIP charge from where ? is it given by index?
   
    if(photons.mX == mipX && photons.mY == mipY /*&& photons.mQ ==  mipCharge*/) { // current cluster is mip
    
      const auto photonPDG = photons.mPDG; // check also other indexes?
      photons.setCandidateStatus(1);
      photons.setCandidateStatusCkov(1);
      cStatusCkov = 1;
	cStatus = 1;
      //LOGP(info, "Cluster PDG {} Track {}", photonPDG, mcTrackPdg);
 


 			mArrAndMap.fillmipSizeFilter(x,y);
      // check if PDG code matches track's PDG-code
      // : photons has field digits w pdg-code; get this and match w track: 

      //std::vector<Topology> topologyVector = photons.getClusterTopology();

      if(true) {

				// this just checks the first digit in the pDigs of the cluster


        if(photonPDG != mcTrackPdg) {
          //Printf("photonPDG != mcTrackPdg");		
					//continue; // not really important atm, see below
					
          
        }
        else {
          //printf("photonPDG matched mcTrackPdg!"); 
          
        } // trackIndex, const std::vector<float>& mipCharges, float mipX, float mipY, const int mcmcTrackPdg
      }
           
    } //<> if()




		auto mipCut = 200.;
		int mipSizeCut = 2;



    
    if(photons.mQ > mipCut && photons.mSize > mipSizeCut && cStatusCkov != 1) {
    
      //LOGP(info, "skipping photon beacause its the mip of another track");

      photons.setCandidateStatus(1);
      photons.setCandidateStatusCkov(1);
      cStatusCkov = 1;
      cStatus= 1;
 			mArrAndMap.fillmipSizeFilter(x,y);
      //continue; // this photon charge is a MIP--> dont consider as candidate 

    }


 
    iPhotCount++;



    //const auto& x = photons[0], y = photons[1];

    


    const auto& etaC = 0;// photons[2];


    double xL = x, yL = y;
    double xG = x, yG = y;


    double thetaCer, phiCer;
    auto xAbs = TMath::Abs(x);
    auto yAbs = TMath::Abs(y);

    bool withinRange = true; 

    const TVector2 posPhoton(x, y);
    
    const auto rPhoton = (posPhoton - mipPos).Mod();
    
    
    // ef : TODO here i put a radius of 40, should this be done?
	bool isPhotonProtonCand = false, isPhotonKaonCand = false, isPhotonPionCand = false;
	bool isPhotonProtonCandCkov = false, isPhotonKaonCandCkov = false, isPhotonPionCandCkov = false;

    const auto rMax = 60.;
    const auto rMin = 2.;
    if((rPhoton < rMax && rPhoton > rMin) && cStatusCkov == 0){
        double thetaCer, phiCer;
	if (findPhotCkov(photons.mX, photons.mY, thetaCer, phiCer)) { // find ckov angle for this  photon candidate

                                  // increment counter of photon candidates
            auto sigmaRing = aliSigma2(thetaP, phiP, thetaCer,
 phiCer);	 // rms of all contributing errors 




	    if(sigmaRing < 0.025) {
	     	Printf("thetaCer %.4f phiCer %.2f, sigmaRing %.5f" , thetaCer, phiCer, sigmaRing);
	     
	     
	        //Printf("getCkovPion %.4f getCkovKaon %.4f, getCkovProton %.4f" , getCkovPion(), getCkovKaon(), getCkovProton());	     
	     }
	     
	    /*if(TMath::Cos(thetaCer) > 1/nF)
    	    { 
	      Printf("cosThetaCer %.2f > 1/nF", TMath::Cos(thetaCer), 1/nF);
	    } */

            // 
	    if (thetaCer > TMath::ASin(1. / (nF - 3*0.0055))) {
	      Printf("jævla stor thetaCer !");
	    } 

	    else {


		    photons.mThetaCer = thetaCer;
		    photons.mPhiCer = phiCer;
		    photons.mSigmaRing = sigmaRing;
		    if(sigmaRing > 0.2) {
		      //Printf("jævla stor sigmaRing !");
		    } else {

		    if(sigmaRing > 0.02) { sigmaRing = 0.02;
		      //Printf("sigmaRing = 0.02 !");
		    }

		      if(TMath::Abs(thetaCer-getCkovPion()) < mSigmaSep*sigmaRing)  {isPhotonPionCandCkov = true; //Printf("ckovPionOk");
pionCandidates.push_back(std::array<double,2>{x,y});
}
		      if(TMath::Abs(thetaCer-getCkovKaon()) < mSigmaSep*sigmaRing) {isPhotonKaonCandCkov = true; //Printf("ckovKaonOk");
kaonCandidates.push_back(std::array<double,2>{x,y});
}
		      if(TMath::Abs(thetaCer-getCkovProton()) < mSigmaSep*sigmaRing) { isPhotonProtonCandCkov = true; //Printf("ckovProtonOk");
protonCandidates.push_back(std::array<double,2>{x,y});} 
		    }
	  }

        }



        const auto phiPhoton = (posPhoton - trkPos).Phi(); 
        																									

        const auto pc = populatePtrOuter->getPcImp();

        

        bool isMaxProtonOk = false, isMinProtonOk = false, isMaxKaonOk = false, isMinKaonOk = false, isMaxPionOk = false, isMinPionOk = false;
  

        
        
        //Printf(" Track PDG %d %s", trackPdg, trackPdgString.c_str());	
        if(getProtonStatus() and getKaonStatus() and getPionStatus()) {
            
            //printf("Pion%d can be Pion Kaon and Proton from p-hyp \n", iPhotCount);	
            // verify thatrPhoton > rMinProton(@ phiEstimated = phiPhoton)

            //printf("\n\npopulate2Ptr->checkUnder(posPhoton, rPhoton, phiPhoton, above, arrMinProton, getMinCkovProton());");
            
            bool shutDownOnOpen = false;
            //isMinProtonOk = populate2Ptr->checkUnder(posPhoton, rPhoton, phiPhoton, arrMinProton, getMinCkovProton(), "Proton", shutDownOnOpen);


								//printf("isMinProtonOk = evaluatePointContour(posPhoton, arrMinProtonPos, false);"); 
						// to be true; point needs to be outside of contour of minProton
            isMinProtonOk = evaluatePointContour(posPhoton, arrMinProtonPos, false); // outside = intersection!
						//printf("isMinProtonOk %d", isMinProtonOk);





            if (shutDownOnOpen) {
        			//printf("shutDownOnOpen!!");	            
            	drawMap = true;
            	errorPos.Set(posPhoton);
            	//break;
            }      
            

						/* checkUnder(const TVector2& posPhoton, const double& rPhoton, const double& phiPhoton, vecArray4& vec, const double& etaCkov, const char* hadronType)
						*/

            // check if rPhoton > rMax(@ phiEstimated = phiPhoton)
            if(isMinProtonOk) {
                //printf("isMinProtonOk = true");
                //printf("\n\npopulate2Ptr->checkOver(posPhoton, rPhoton, phiPhoton, arrMaxPion, getMaxCkovPion());");






                //isMaxPionOk = populate2Ptr->checkOver(posPhoton, rPhoton, phiPhoton, arrMaxPion, getMaxCkovPion(), "Pion");
								//printf("isMaxPionOk = !evaluatePointContour(posPhoton,  arrMaxPionPos, true);"); 

								// to be true; point needs to be inside of contour of isMaxPionOk
    						isMaxPionOk = !evaluatePointContour(posPhoton,  arrMaxPionPos, true); // inside = no intersection!
								//printf("isMaxPionOk %d", isMaxPionOk);
                // this means rProtonMin < rPhoton < rMaxPion
                if(isMaxPionOk) {
                    
                    
                    // this means rProtonMax > rPhoton; then also the rPh < rPionMax and rPh < rKaonMax

                    //isMaxProtonOk = populate2Ptr->checkOver(posPhoton, rPhoton, phiPhoton, arrMaxProton, getMaxCkovProton(), "Proton");


										// to be true; point needs to be inside of contour of isMaxProtonOk
    								//printf("		isMaxProtonOk = !evaluatePointContour(posPhoton,  arrMaxProtonPos, true)"); 
    								isMaxProtonOk = !evaluatePointContour(posPhoton,  arrMaxProtonPos, true);// inside = no intersection!
										//printf("isMaxProtonOk %d", isMaxProtonOk);
                    //printf("===================================="); 
                    if(isMaxProtonOk) {
                        // we have proton-candiate
                        isPhotonProtonCand = true;

                        //protonCandidates.push_back(std::array<double,2>{x,y});
                        //printf("Photon%d is a Proton Candidate", iPhotCount); 
                        // isMaxPionOk = true;
                        isMaxKaonOk = true;
                    }	else {
                        //printf("Photon%d not a Proton Candidate", iPhotCount); 
                    }
                    //printf("====================================\n"); 
                    //printf("===================================="); 
                    // this means rPionMin < rPhoton --> and also rKaonMin < rPhoton



								    bool shutDownOnOpen = false;
								    

      								//printf("	isMinPionOk = evaluatePointContour(posPhoton,  arrMinPionPos, false);"); 
                    //isMinPionOk = populate2Ptr->checkUnder(posPhoton, rPhoton, phiPhoton, arrMinPion, getMinCkovPion(), "Pion", shutDownOnOpen);
    								isMinPionOk = evaluatePointContour(posPhoton,  arrMinPionPos, false); 										//printf("isMinPionOk %d", isMinPionOk);
								    if (shutDownOnOpen) {
            					errorPos.Set(posPhoton);								    
											//printf("shutDownOnOpen!!");	            
								    	drawMap = true;
								    	//break;
								    }                    
                    
                    if(isMinPionOk) {
                        // we have pion-candiate
                        isPhotonPionCand = true;
                        //printf("Photon%d is a Pion Candidate", iPhotCount); 
                                            
                        //pionCandidates.push_back(std::array<double,2>{x,y});
                        
                        isMinKaonOk = true;
                    }	else {
                        //printf("Photon%d not a Pion Candidate", iPhotCount); 
                    }			
                        //printf("===================================="); 				
                    
                    // if rPhoton < rPhotonMin, check if rPhoton > rKaonMin
                    if(!isMinKaonOk) {
                    
		                  
										  bool shutDownOnOpen = false;
										  
											//isMinKaonOk = populate2Ptr->checkUnder(posPhoton, rPhoton, phiPhoton, arrMinKaon, getMinCkovKaon(), "Kaon", shutDownOnOpen);
      								//printf("isMinKaonOk = evaluatePointContour(posPhoton,  arrMinKaonPos, false);"); 
    									isMinKaonOk = evaluatePointContour(posPhoton,  arrMinKaonPos, false); //printf("isMinKaonOk %d", isMinKaonOk);

										  if (shutDownOnOpen) {
					            	errorPos.Set(posPhoton);										  
												//printf("shutDownOnOpen!!");	            
										  	drawMap = true;
										  	//break;
										  }                          

                    }
                    
                    //printf("\n===================================="); 
                    // only check isMaxKaon if isMinKaonOk
                    if(isMinKaonOk) {
                        if(!isMaxKaonOk) {
		        								//printf("isMaxKaonOk = !evaluatePointContour(posPhoton,  arrMaxKaonPos, true)"); 
                            //isMaxKaonOk = populate2Ptr->checkOver(posPhoton, rPhoton, phiPhoton, arrMaxKaon, getMaxCkovKaon(), "Kaon");
    												isMaxKaonOk = !evaluatePointContour(posPhoton,  arrMaxKaonPos, true); // inside = no intersection!
                            // isMaxKaonOk = populate2Ptr->checkCond()
	 														//printf("isMaxKaonOk %d", isMaxKaonOk);
                        } 

                        if(isMaxKaonOk and isMinKaonOk) {
                            //printf("Photon%d is a Kaon Candidate", iPhotCount); 
                            // we have kaon cand
                            isPhotonKaonCand = true;
                            //kaonCandidates.push_back(std::array<double,2>{x,y});
                        } else {
                            //printf("Photon%d not a Kaon Candidate", iPhotCount); 
                        }
                    } else {
                        //printf("Photon%d not a Kaon Candidate", iPhotCount); 
                    }
                    //printf("====================================\n"); 
                                
                }  // end if isMaxPionOk == true
                            
                // else isMaxPionOk == false
                else {
                    //printf("===============================================================");	
                    //printf("\n Photon%d not a Hadron Candidate :\n rPhoton %.2f > rPionMax", iPhotCount,  rPhoton);
                    //printf("===============================================================");	
                }
                
            }	// end if isMinProtonOk == true
            
            // else isMinProtonOk == false
            // this means  rPhoton < rProtonMin : 
            else {
                //printf("===============================================================");	
                //printf("\n Photon%d not a Hadron Candidate :\n rPhoton %.2f < rProtonMin",iPhotCount, rPhoton);
                //printf("===============================================================");	
            }

        
        } // end if getProtonStatus() and getKaonStatus and getPionStatus

        // if Proton not is possible because momentum-threshold is not exceeded:
        else if(!getProtonStatus() and getKaonStatus() and getPionStatus()) {
            //printf("\n====================================================="); 
            //printf("Photon%d can be Pion and Kaon from p-hyp", iPhotCount);	
            // verify that rPhoton > rMinKaon(@ phiEstimated = phiPhoton)




					  bool shutDownOnOpen = false;
          	//printf("isMinKaonOk = evaluatePointContour(posPhoton,  arrMinKaonPos, false)"); 
						//isMinKaonOk = populate2Ptr->checkUnder(posPhoton, rPhoton, phiPhoton, arrMinKaon, getMinCkovKaon(), "Kaon", shutDownOnOpen);
						isMinKaonOk = evaluatePointContour(posPhoton,  arrMinKaonPos, false); // outside =  intersection!
					  if (shutDownOnOpen) {            	errorPos.Set(posPhoton);
							//printf("shutDownOnOpen!!");	            
					  	drawMap = true;
					  	//break;
					  }     

            // this means rProtonMin < rPhoton < rMaxPion
            if(isMinKaonOk) {
            
            
                // denne var feil!!! her var det satt arrRMinPion og gtckovMin
                //isMaxPionOk = populate2Ptr->checkOver(posPhoton, rPhoton, phiPhoton, arrMaxPion, getMaxCkovPion(), "Pion");
                  //printf("isMaxPionOk = !evaluatePointContour(posPhoton,  arrMaxPionPos, true)"); 
								isMaxPionOk = !evaluatePointContour(posPhoton,  arrMaxPionPos, true); // inside =  no intersection!
                if(isMaxPionOk) {
										bool shutDownOnOpen = false;
                    //isMinPionOk = populate2Ptr->checkUnder(posPhoton, rPhoton, phiPhoton, arrMinPion, getMinCkovPion(), "Pion", shutDownOnOpen);
                        //printf("isMinPionOk = evaluatePointContour(posPhoton,  arrMinPionPos, false);"); 
										isMinPionOk = evaluatePointContour(posPhoton,  arrMinPionPos, false); // outside =  intersection!


										if (shutDownOnOpen) {            	errorPos.Set(posPhoton);
											//printf("shutDownOnOpen!!");	            
											drawMap = true;
											//break;
										}                      
                    
                    
                //isMinPionOk = ...
                    // this means rPionMin < rPhoton
                    if(isMinPionOk) {
                        isPhotonPionCand = true;
                        //printf("Photon%d is Pion Candiate", iPhotCount); 
                        // we have pion-candiate
                    }	

                    //isMaxKaonOk = populate2Ptr->checkOver(posPhoton, rPhoton, phiPhoton, arrMaxKaon, getMaxCkovKaon(), "Kaon");
                        //printf("isMaxKaonOk = !evaluatePointContour(posPhoton,  arrMaxKaonPos, true)"); 
										isMaxKaonOk = !evaluatePointContour(posPhoton,  arrMaxKaonPos, true); // inside =  no intersection!


                    // isMaxKaonOk = populate2Ptr->
                    if(isMaxKaonOk) {
                        isPhotonKaonCand = true;
                        //printf("Photon%d is Kaon Candiate", iPhotCount); 
                        // we have kaon-candiate
                    }							
                }	
                else {
                    //printf("Photon%d could be Pion/Kaon (from p-hyp) but was out of radius-range",iPhotCount);
                }	
            } // end if isMinKaonOk
            else {
                //printf("Photon%d could be Pion/Kaon (from p-hyp) but was out of radius-range",iPhotCount);
            }
        } // end else if (!getProtonStatus() and getKaonStatus() and getPionStatus())



        // denne ser litt rar ut?

        // if neither Proton or Kaon is possible because momentum-threshold is not exceeded:
        else if(!getProtonStatus() and !getKaonStatus() and getPionStatus()) {
            //printf("Photon%d can be Pion from p-hyp", iPhotCount);	

            //printf("populate2Ptr->checkUnder(posPhoton, rPhoton, phiPhotob, arrMaxPion, getMaxCkovPion());");

            //isMaxPionOk = populate2Ptr->checkOver(posPhoton, rPhoton, phiPhoton, arrMaxPion, getMaxCkovPion(), "Pion");
						isMaxPionOk = !evaluatePointContour(posPhoton,  arrMaxPionPos, true);      // inside =  no intersection!      
            
            if(isMaxPionOk) {

                //printf("populate2Ptr->checkOver(posPhoton, rPhoton, phiPhoton, arrMinPion, getMaxCkovPion());",iPhotCount);
                
								bool shutDownOnOpen = false;                
                //bool isMinPionOk = populate2Ptr->checkUnder(posPhoton, rPhoton, phiPhoton, arrMinPion, getMinCkovPion(), "Pion", shutDownOnOpen);
								bool isMinPionOk = evaluatePointContour(posPhoton,  arrMinPionPos, false); // outisde = intersection!      

								if (shutDownOnOpen) {            	errorPos.Set(posPhoton);
									//Printf("shutDownOnOpen!!");	            
									drawMap = true;
									//break;
								}       
                if(isMinPionOk) {
                    // we have pion candidate
                    //printf("Photon%d is Pion Candidate", iPhotCount);	
                    isPhotonPionCand = true; 
                    
                } else {
                    //printf("Photon%d could only have been Pion, but didnt fall within radius-range",iPhotCount);
                }
            } 
            else {
                //printf("Photon%d could only have been Pion, but didnt fall within radius-range",iPhotCount);
            }			
            
        } // end else if(!getProtonStatus() and !getKaonStatus() and getPionStatus())

        // this should not really be "possible", but can be if we have another candidate that is not pion/kaon/proton
        else {
            //printf("Photon%d was not Pion Kaon or Proton from p-hyp",iPhotCount);
            // we got particle that should not be able to be pion, kaon or proton
        }


        
       /*if (isPhotonPionCandCkov || isPhotonKaonCandCkov || isPhotonProtonCandCkov) { printf("\n======================");	
        printf("	Photon%d  Candidate Ckov: Pion = %d Kaon = %d Proton %d", iPhotCount, isPhotonPionCandCkov, isPhotonKaonCandCkov, isPhotonProtonCandCkov);
		printf("	Photon%d  Candidate segm: Pion = %d Kaon = %d Proton %d", iPhotCount, isPhotonPionCand, isPhotonKaonCand, isPhotonProtonCand);

        printf("============================\n");	
}*/


	      if(cStatusCkov == 0) {

auto pionBit = 4*static_cast<int>(isPhotonPionCandCkov); 
auto kaonBit =2*static_cast<int>(isPhotonKaonCandCkov); 
auto protonBit = 1*static_cast<int>(isPhotonProtonCandCkov);

          	cStatusCkov = 1  + pionBit + kaonBit +protonBit;

//printf("	Photon%d  Candidate Ckov: pionBit = %d kaonBit = %d protonBit %d", iPhotCount, pionBit, kaonBit, protonBit);


	      }


  	      if(cStatus == 0){
          cStatus = 1  + 4*static_cast<int>(isPhotonPionCand) + 2*static_cast<int>(isPhotonKaonCand) +


 1*static_cast<int>(isPhotonProtonCand);}



	      if(cStatusCkov == 1)
	        mArrAndMap.fillckovCandMapOutRange(x,y);
				//LOGP(info, "cStatusCkov {}", cStatusCkov);

            //o2::hmpid::ClusterCandidate :setCandidateStatus(int iTrack, int hadronCandidateBit)
            photons.setCandidateStatus(cStatus);

            photons.setCandidateStatusCkov(cStatusCkov);

	

	//LOGP(info, "Set Track {} : cStatusCkov {} | cStatus {} ", trackIndex, cStatusCkov, cStatus);
       


        }// end if radius ok
        
        
        // radius was to big to consider : 

        
        else if((rPhoton > rMax || rPhoton < rMin ) && cStatusCkov == 0){
            rOverMax++;
            cStatusCkov = 0; // set other value to indicate out of region?

	    cStatus = 0; // set other value to indicate out of region?
            photons.setCandidateStatus(cStatusCkov);
            photons.setCandidateStatusCkov(cStatusCkov);
						//LOGP(info, "Radius {} too high ! cStatusCkov {}", rPhoton, cStatusCkov);

		      	//o2::hmpid::ClusterCandidate :setCandidateStatus(int iTrack, int hadronCandidateBit)
		      					
        }

        // this means it falls within range
        if(cStatusCkov > 1) {
							canCombined.push_back(std::array<double,2>{x,y});
              // this means  candidate is a ckov-photon
              if(true) {
                  //(mArrAndMap->ckovCandMapRange)->Fill(x,y);
                  mArrAndMap.fillCkovCandMapRange(x,y); 
              }	else { // falls within range, but is bg
              //(mArrAndMap->mipSizeFilter)->Fill(x,y);
                  mArrAndMap.fillmipSizeFilter(x,y);
              }
          }
          // falls out of range
          else if (cStatusCkov == 0){
              // this means  candidate is a ckov-photon, but out of range
              if(true) {
                  mArrAndMap.fillckovCandMapOutRange(x,y);
                  //(mArrAndMap->ckovCandMapOutRange)->Fill(x,y);
              }	else { // falls out of range, and is bg
                  mArrAndMap.fillMipCharge(x,y);
                  //(mArrAndMap->mipChargeFilter)->Fill(x,y);			
              }
          }

        
    }   // end for ckovPhotons


	// iterate through photons to check : 


    /*
    Printf("	Identified the following  : \n");
    Printf("	numBackgroundPhotons %d", numBackgroundPhotons);
    Printf("	numFoundActualCkov %d", numFoundActualCkov);
    Printf("	numActualCkov %d", numActualCkov);
    Printf("	numBackgroundLabeledCkov %d", numBackgroundLabeledCkov);        
    */ 
        
        /*
    if(numFoundActualCkov != numActualCkov) {
                throw std::invalid_argument("wtf value does numActualCkov have?");
    }*/ 

	//printf("=========================================================================");
	//printf("CkovHyps, possible candidates from Momentum | Pion = %d, Kaon = %d, Proton = %d", getPionStatus(), getKaonStatus(), getProtonStatus());


    /*int cntTemp = 0;
  for(const auto& c : candidatesCombined){
		////printf("Phot%d : x = %.2f, y = %.2f || statusCand = %d", cntTemp++, c.x, c.y, c.candStatus); 
	}*/
	//printf("=========================================================================");


    //printf("number of candidates : proton %zu, kaon %zu, pion %zu | total %zu, radius over Limit %d ", protonCandidates.size(), kaonCandidates.size(), pionCandidates.size(), clusterTrack.size(), rOverMax);
 
      //printf("MIP(%.1f %.1f) RAD (%.1f %.1f) Track (%.1f %.1f)", mipX, mipY, xRad, yRad, thetaP, phiP);
      
      
      
    //for(const auto& pair: filledBins)
    //	//printf("CkovTools segment candidates: x%f y%f", pair.first, pair.second);    


    // hNoiseMap->SetMarkerColor(kRed);
    //printf("CkovTools segment filledBins Size %zu", filledBins.size());


   /*
    // get impact points of track RAD and PC
    const auto trkPC = populate.getPcImp();
    const auto trkRad = populate.getTrackPos();
    TH2F* trkPCMap = new TH2F("trkPCMap ", "trkPCMap; x [cm]; y [cm]",160*20,0.,159.,144*20,0,143);
    TH2F* trkRadMap = new TH2F("trkRadMap ", "trkRadMap; x [cm]; y [cm]",160*20,0.,159.,144*20,0,143);
    trkRadMap->Fill(trkRad.X(), trkRad.Y());
    trkPCMap->Fill(trkPC.X(), trkPC.Y());





   
    c->SetMarkerColor(kGreen + 4);
    localRefMIP2->Draw("same");  

    TCanvas *segm = new TCanvas("semg","semg",800,800);  
    segm->Divide(2,2);

    segm->cd(1);
    mapPhotons->SetMarkerStyle(2); 

    mapPhotons->Draw();
    localRefMIP2->Draw("same");  
    trkPCMap->Draw("same"); trkRadMap->Draw("same"); 


    segm->cd(2);
    mapPion->SetTitle(Form("mapPion : min%d max%d photons%d", mapPionMin->GetEntries(), mapPionMax->GetEntries(), mapPhotons->GetEntries()));
    mapPion->Draw();
    mapPhotons->Draw("same");
    localRefMIP2->Draw("same");
    mapPionMin->Draw("same");
    mapPionMax->Draw("same");
    trkPCMap->Draw("same"); trkRadMap->Draw("same"); 

    segm->cd(3);

    mapKaon->SetTitle(Form("mapKaon: min%d max%d photons%d", mapKaonMin->GetEntries(), mapKaonMax->GetEntries(), mapPhotons->GetEntries()));
    mapKaon->Draw();
    mapPhotons->Draw("same");
    mapKaonMin->Draw("same");
    mapKaonMax->Draw("same");
    trkPCMap->Draw("same"); trkRadMap->Draw("same"); 

    segm->cd(4);
    mapProton->SetTitle(Form("mapProton : min%d max%d photons%d", mapProtonMin->GetEntries(), mapProtonMax->GetEntries(), mapPhotons->GetEntries()));
    mapProton->Draw();

    mapPhotons->Draw("same");
    mapProtonMin->Draw("same");
    mapProtonMax->Draw("same");
    trkPCMap->Draw("same"); trkRadMap->Draw("same"); 

    gPad->Update();
    segm->Show(); / */





  // to draw maps: false
	if(/*drawMap*/ false) {

		
 		mArrAndMap.setEventCount(eventCnt);
		
		mArrAndMap.setErrorPos(errorPos);

		// set all clusters to see plot .
		//mArrAndMap.setAllClusters(allCand);



mArrAndMap.setPopulatePtr(std::move(populatePtrCp));

		
		if(false) { //printf("populate2Ptr was nullptr!");
    }
	  else {	



    int numPion = pionCandidates.size();
		int numKaon = kaonCandidates.size();
		int numProton = protonCandidates.size();
		int numTotal = canCombined.size();
		 
		Printf("SigmaSep %.2f  | nPi %d nK %d nPr %d == nT %d | PDG %d | thetaP %.2f phiP %.2f", mSigmaSep, numPion, numKaon,numProton,numTotal, trackPdg, thetaP, phiP);

		


		mArrAndMap.setMinArrays(arrMinPionPos, arrMinKaonPos, arrMinProtonPos);
		mArrAndMap.setMaxArrays(arrMaxPionPos, arrMaxKaonPos, arrMaxProtonPos);		

			//printf("populatePtr->GetHistogram()->GetEntries() num BinEntries = %f ", populate2Ptr->getLimMin()->GetEntries());
			//printf("mArrAndMap->// limMin->GetEntries() num BinEntries = %f ", populate2Ptr->// limMin->GetEntries());
    // to draw the "split-phi"...
		//mArrAndMap.setDrawLimits(populate2Ptr->// limMin, populate2Ptr->// limMax);

	  			//printf("mArrAndMap->// limMin->GetEntries() num BinEntries = %f ", populate2Ptr->// limMin->GetEntries());
	  
		auto d = (rW- L + tGap + qW);


		auto tanP = TMath::Tan(thetaP);
		
		auto cosP = TMath::Cos(phiP);
		auto sinP = TMath::Sin(phiP);
	
		//printf("delta %.2f | x %.2f y %.2f ", d*tanP, d*tanP*cosP, d*tanP*sinP);
		//mArrAndMap->drawTotalMapAndMaxRegions();
				




    mArrAndMap.drawTotalMap(clusterTrack, plotNumber, xMip, yMip, pionCandidates, kaonCandidates, protonCandidates, canCombined, trackPdg, allCand, mSigmaSep);
    // to drqw the maps :: 

		

   //mArrAndMap->drawMaxRegions();		  
  		//printf("==========================================================================");
  		//printf("==========================================================================");
  		//printf("==========================================================================");

		//printf("\n\n\nEvent Number%d, momentum %.2f | mass %.2f | thetaP %.2f | phiP %.2f | L %.2f", eventCnt,momentum, mass, thetaP, phiP, L);
					
		// const double& ckovThe, const double& ckovPhi, const double & L
		const auto l1 = populatePtr->tracePhot(getMaxCkovPion(), 0, L);
		const auto l2 = populatePtr->tracePhot(getMinCkovPion(), 3.14159285, L);

		//printf("phiP %.5f: angle rad --> pc = %.5f | Acos %.5f Asin %.5f"  , phiP, (trkPos-trkPC).Phi(), TMath::ACos(d*tanP*cosP/(trkPos-trkPC).Mod()), TMath::ASin(d*tanP*sinP/(trkPos-trkPC).Mod()));
		
		//printf("PHI : l1 %.5f | l2 %.5f ", (l1-trkPos).Phi(),(l2-trkPos).Phi());
		//printf("RADIUS : l1 %.5f | l2 %.5f ", (l1-trkPC).Mod(),(l2-trkPC).Mod());

		//printf("CkovHyps %.2f %.2f %.2f", ckovHyps[0], ckovHyps[1], ckovHyps[2]);
		//Printf("CkovHyps %.2f %.2f %.2f", ckovHyps[0], ckovHyps[1], ckovHyps[2]);
    //std::this_thread::sleep_for(std::chrono::seconds(0.1));

    //return;
		//throw std::invalid_argument("print invoked"); 
    }
  // drawTotalMap / drawMaxRegions
  } else {
    // NB! lagre denne istedet : s

    /*
    candCombined = candidatesCombined; // x, y, cStatusCkov

    protonCands = protonCandidates;//.emplace_back()
    kaonCands = kaonCandidates;//.emplace_back()
    pionCands = pionCandidates;//.emplace_back()
   */ 


    // returner disse ogsaa: kun for debugging :

    //arrayInfo = {numBackgroundPhotons, numFoundActualCkov, numActualCkov, numBackgroundLabeledCkov}; 

   	//delete mArrAndMap;	
   }
	return filledBins;
} // end segment


bool lineIntersects(const TVector2& B, const std::array<double, 2>& C, const std::array<double, 2>& D, double rUncertainty, TVector2& intersection, bool& pointIsInside  ) { 

    // Existing code for calculating line intersection
    double a1 = B.Y() - mipPos.Y();
    double b1 = mipPos.X() - B.X();
    double c1 = a1 * mipPos.X() + b1 * mipPos.Y();

    double a2 = D[1] - C[1];
    double b2 = C[0] - D[0];
    double c2 = a2 * C[0] + b2 * C[1];

    double det = a1 * b2 - a2 * b1;
    if (det == 0) {

				// Printf("dlineIntersects :: parellel lines, returnign false");
        return false;  // Lines are parallel
    }

    double x = (b2 * c1 - b1 * c2) / det;
    double y = (a1 * c2 - a2 * c1) / det;




    intersection.Set(x, y);

		bool isOutOfBounds = (x < std::min(mipPos.X(), B.X()) || x > std::max(mipPos.X(), B.X()) ||
    	                  x < std::min(C[0], D[0]) || x > std::max(C[0], D[0]) ||
    	                  y < std::min(mipPos.Y(), B.Y()) || y > std::max(mipPos.Y(), B.Y()) ||
    	                  y < std::min(C[1], D[1]) || y > std::max(C[1], D[1]));



	if (isOutOfBounds) {
			// Printf("Out of bounds, adj nec");
		  TVector2 unitDirMIPtoIntersection = (intersection - mipPos).Unit();
		  TVector2 adjustedIntersection = intersection + unitDirMIPtoIntersection * rUncertainty;
		  x = adjustedIntersection.X();
		  y = adjustedIntersection.Y();
		  isOutOfBounds = (x < std::min(mipPos.X(), B.X()) || x > std::max(mipPos.X(), B.X()) ||
		                   x < std::min(C[0], D[0]) || x > std::max(C[0], D[0]) ||
		                   y < std::min(mipPos.Y(), B.Y()) || y > std::max(mipPos.Y(), B.Y()) ||
		                   y < std::min(C[1], D[1]) || y > std::max(C[1], D[1]));

			auto dist1 = (B - mipPos).Mod();
			auto dist2 = (adjustedIntersection - mipPos).Mod();


		  if (isOutOfBounds) {
					// Printf("dlineIntersects :: adjusted point is still out of bounds, returning false");
		      return false;  // Adjusted intersection is still out of bounds
		  } else {
		      intersection = adjustedIntersection;
    			pointIsInside = (adjustedIntersection - mipPos).Mod() < (B - mipPos).Mod();
					//printf("OK :: dist2 Point %.2f | dist2 adj %.2f |  adj %.2f", dist1, dist2, rUncertainty);
		  }
	} else { 

			//printf("Not out of bounds, adj not nec");
			pointIsInside = (intersection - mipPos).Mod() < (B - mipPos).Mod();
			auto dist1 = (B - mipPos).Mod();
			auto dist2 = (intersection - mipPos).Mod();
			//printf("dist2 Point %.2f | dist2 adj %.2f", dist1, dist2, rUncertainty);
  }


    return true;
}

bool evaluatePointContour(const TVector2& posPhoton, const vecArray2& contour, bool checkInside) {

    double rUncertainty = 2;
    if (!checkInside) {
        rUncertainty = -2;
    }

    bool intersects = false;
    bool pointIsInside = false; // To hold the result from lineIntersects
    const auto& shortSidePoint = contour[0];
    TVector2 intersection; // To store the intersection point

    double d = (posPhoton.X() - shortSidePoint[0]) * (mipPos.Y() - shortSidePoint[1]) -
               (posPhoton.Y() - shortSidePoint[1]) * (mipPos.X() - shortSidePoint[0]);

    if (d > 0) {
        for (size_t i = 0; i < contour.size() - 1; ++i) {
            size_t next = i + 1;
            if (lineIntersects(posPhoton, contour[i], contour[next], rUncertainty, intersection, pointIsInside)) {
                intersects = true;
                if (pointIsInside) return true; // Point is inside the contour
                else return false; // Point is outside the contour
            }
        }
    } else {
        for (size_t i = contour.size() - 1; i > 0; --i) {
            size_t prev = i - 1;
            if (lineIntersects(posPhoton, contour[i], contour[prev], rUncertainty, intersection, pointIsInside)) {
                intersects = true;
                if (pointIsInside) return true; // Point is inside the contour
                else return false; // Point is outside the contour
            }
        }
    }
		//printf("evaluatePointContour checkInside %d  : pointIsInside == %d", checkInside, pointIsInside);
    return intersects;
}





// Trace a single Ckov photon from emission point somewhere in radiator up to photocathode taking into account ref indexes of materials it travereses
// Arguments: ckovThe,ckovPhi- photon ckov angles in TRS, [rad]
//   Returns: distance between photon point on PC and track projection
// placeholder...
void setArrayMax(double etaTRS, vecArray4& inPutVectorAngle, vecArray2& inPutVectorPos, const size_t kN)
{
  // const size_t kN = inPutVector.size();
  const float lMin = 0.;
  const auto trkPC2 = populatePtrOuter->getPcImp();      // track at PC
  const auto trkRad2 = populatePtrOuter->getTrackPos();  // track at RAD


  //Printf("\n\n setArrayMax() enter --> kN %zu, etaTRS %.2f", kN, etaTRS);
	for(int i = 0; i < kN; i++){

		const auto phiL = double(TMath::TwoPi()*(i+1)/kN);
		TVector3 dirTrs, dirLORS;
		
		// set TRS values :
		dirTrs.SetMagThetaPhi(1, etaTRS, phiL);

		//Printf("\n setArrayMax() dirTrs.SetMagThetaPhi(1, etaTRS = x %.2f, phiL = x %.2f); ", etaTRS, phiL);
		

		double thetaR, phiR; // phiR is value of phi @ estimated R in LORS


		populatePtrOuter->trs2Lors(dirTrs, thetaR, phiR);
		
		//Printf("setArrayMax() called trs2Lors on dirTRS : return { thetaR = %.2f, phiR = %.2f}", thetaR, phiR);
		

		dirLORS.SetMagThetaPhi(1, thetaR, phiR);

		//Printf("setArrayMax() dirLORS {x %.2f y %.2f z %.2f}", dirLORS.X(), dirLORS.Y(), dirLORS.Z());


		// temp
		// ckovThe, const double& ckovPhi, const double & L
		// this should return the same as max
		//const auto t = populatePtrOuter->tracePhot(etaTRS, phiL, lMin);

		// temp
		
		//Printf("setArrayMax() called  populatePtrOuter->traceForward(dirLORS (thetaR %.2f, phiR %.2f) lMin =  %.2f", thetaR, phiR, lMin);
		const auto& max = populatePtrOuter->traceForward(dirLORS, lMin); 
		// max = pos of tracked photon at PC
		
		const auto r = (max - trkPC2).Mod(); // trkPC2 : track impact @ PC
				// from alinot_paattrec : this is MIP2photon distance?
				// const auto r = (max - trkPC2).Mod(); // trkPC2 : track impact @ PC
		
		
		
		auto phiPC = (max - trkPC2).Phi();   // MIP phi value/ PC track value?
		// vector fra photon til MIP/PC (--> verdi som kommer fra phiL)
		
		
		//Printf("setArrayMax() traceForward returned TVector2 {x %.2f y %.2f} == > R  = %.2f", max.X(), max.Y(), r);
		//Printf("setArrayMax() tracePhot returned TVector2 {x %.2f y %.2f} == > R  = %.2f", t.X(), t.Y(), (t-trkPC2).Mod());
		
		// add protection if traceForward returns 0 or -999?


 		//Printf("setArrayMax() --> i %d, {x %.2f y %.2f} - MIP {x %.2f y %.2f}", i, max.X(), max.Y(), trkPC.X(), trkPC.Y());
 		

		// check at phiR  er det samme som (t-rad).phi og (max-rad).phi?
		



		//Printf("setArrayMax2() --> i %d, {x %.2f y %.2f} - MIP {x %.2f y %.2f}", i, max.X(), max.Y(), trkPC2.X(), trkPC2.Y());


		// phiR in [-pi, pi]? set to 0..2pi?
		// inPutVector.emplace_back(std::array<double, 3>{phiL, phiR, r});
		if(phiR < 0) {
			phiR = TMath::TwoPi() + phiR;

    }
     
		if(phiPC < 0) {
			phiPC = TMath::TwoPi() + phiPC;

    }   

    
		/*if( (max-trkRad2).Phi() != (t-trkRad2).Phi()) {
			throw std::invalid_argument("(max-trkRad2).Phi() != (t-trkRad2).Phi())");
		}

		if(TMath::Abs(phiR - (t-trkRad2).Phi()) > 0.0001 && t.X() != -999 && t.Y() != -999)   
		{
		
		  Printf("phiR at %.6f | (t-rad) %.6f | (max-rad) %.6f ", phiR, (t-trkRad2).Phi(), (max-trkRad2).Phi());
			throw std::invalid_argument("phiR != (t-trkRad2).Phi())");
		}
		// protections if r > value?
		//if(r > )*/


		// TODO: if it goes out of map, find intersection with chamber-edges??
		// really jsut check wether TraceForward returned x, y = -999.
 		// to set points out of the map here is fine	


		if((max.Y() == -999) or (max.X() == -999)) {
    	//inPutVector[i] = {0,0,0, 0, 0};
    	//inPutVector.emplace_back(std::array<double, 3>{0, 0, 0});
			// placeholder, find better solution? 
			if(max.Y() == -999) {
				//Printf("setArrayMax() max.Y() %.1f == -999", max.Y());
			}
			if(max.X() == -999) {
				//Printf("setArrayMax() max.X() %.1f == -999", max.X());
			}
    } else {
    	
    	//inPutVector.emplace_back(std::array<double, 3>{phiL, phiR, r});
    	
    	
    	// phiL : photon_phi i TRS system 
    	// phiR : photon_phi i LORS system
    	// phiPC : (photon - MIP).Phi();
    	 
    	inPutVectorPos.emplace_back(std::array<double, 2>{max.X(), max.Y()}); 
    	inPutVectorAngle.emplace_back(std::array<double, 4>{phiL, phiR, phiPC, r});   	
			const auto r = (max - trkPC2).Mod(); // trkPC2 : track impact @ PC
			// from alinot_paattrec : this is MIP2photon distance?
			// const auto r = (max - trkPC2).Mod(); // trkPC2 : track impact @ PC

    	
    	
    	//inPutVectorAngle.emplace_back(std::array<double, 3>{phiL, phiR, r});
    	//Printf("setArrayMax() emplacing element %d : phiL %.2f, phiR %.2f, r %.2f", i, phiL, phiR, r);
    }
    	//inPutVector[i] = {phiL, phiR, r};
    // Printf("setArrayMax() emplacing element %d : phiL %.2f, phiR %.2f, r %.2f", i, phiL, phiR, r);

		/*if(maxPion.X() > 0 && maxPion.X() < 156.0 && maxPion.Y() > 0 && maxPion.Y() < 144) {
			// hMaxPion->Fill(maxPion.X(), maxPion.Y());
			// maxPionVec[i] = std::make_pair(maxPion.X(), maxPion.Y());
			Printf("maxPion loop i = %d, maxSize = %zu", i, maxPionVec.size()); 
		}*/
	}/*
  	Printf("\n");
	for(const auto& ip : inPutVectorAngle) {
		const auto& phiL_ = ip[0];
		const auto& phiR_ = ip[1]; 
		const auto& r_ = ip[2];  
		//Printf("setArrayMax() --> checking inputVector | : phiL %.2f, phiR %.2f, r %.2f", phiL_, phiR_, r_);
		if(r_ == 0 ) {throw std::invalid_argument("r====??????;");}
  } */
}



double calculatePhi(const std::array<double, 2>& a, const std::array<double, 2>& b) {
    double dx = a[0] - b[0];
    double dy = a[1] - b[1];
    return std::atan2(dy, dx);  // atan2 gives the angle in [-π, π] range
}



// chane function name 
void calculateDifference(vecArray2& maxVec, vecArray2& minVec, double ckovHypVal, const Sigma2& sigma) {
//vecArray2 calculateDifference(const vecArray2& maxVec, const vecArray2& minVec) {
    assert(maxVec.size() == minVec.size());

		const int kN = 200;
    for (size_t i = 0; i < kN; ++i) {
			  const auto phiL = double(TMath::TwoPi()*(i+1)/ kN);
				const double sigma2 = sigma.sigma2(phiL); // calculate the combined 	 
				TVector2 min, max;

				

        auto sigma2Ali = 2*aliSigma2(thetaP, phiP, ckovHypVal, phiL);

//Printf(" phiP TRS_deg %.2f | thetaP_deg %.2f ckovHypVal %.2f | sigma %.5f sigma2Ali %.5f ", phiL*180/3.14, thetaP*180/3.14,  ckovHypVal, sigma2, sigma2Ali);

				setMeanSingle(min, ckovHypVal - sigma2Ali, phiL);
				setMeanSingle(max, ckovHypVal + sigma2Ali, phiL);

				if((max.Y() == -999) or (max.X() == -999)) {
					Printf("max.Y() == -999) or (max.X() == -999");
				} else {
					Printf("maxVec %.2f  %.2f" , max.X(), max.Y());							 
					maxVec.emplace_back(std::array<double, 2>{max.X(), max.Y()});
				}
				if((min.Y() == -999) or (min.X() == -999)) {
					Printf("min.Y() == -999) or (min.X() == -999");
				} else {										 
					//Printf("maminVecxVec %.2f  %.2f" , min.X(), min.Y());							 
					minVec.emplace_back(std::array<double, 2>{min.X(), min.Y()});
				}
				
    } 

}


double getTotalSigmaPos(double sigmaPositions)
{

	// sigmaAngles sigmaE sigmaThetaP --- sigmaPositions
	Printf(" getTotalSigmaPos : sigmaPositionssq %.3f  sigmaLsq %.23f sigmaRsq %.3f", sigmaPositions * sigmaPositions, sigmaLsq, sigmaRsq);
	return TMath::Sqrt(sigmaPositions * sigmaPositions + sigmaLsq + sigmaRsq);

}


void setMeanSingle(TVector2& vectorOut, double etaTRS, double phiL)
{
  const float lMean =  0.75;

	TVector3 dirTrs, dirLORS;
		
	// set TRS values :
	dirTrs.SetMagThetaPhi(1, etaTRS, phiL);
	double thetaR, phiR; // phiR is value of phi @ estimated R in LORS
	populatePtrOuter->trs2Lors(dirTrs, thetaR, phiR);
	dirLORS.SetMagThetaPhi(1, thetaR, phiR);		
	vectorOut = populatePtrInner->traceForward(dirLORS, lMean); 
}


void setMinSingle(TVector2& vectorOut, double etaTRS, double phiL)
{
  const float lMax = 1.5;
  const float lMean =  0.75;
	TVector3 dirTrs, dirLORS;
		

	dirTrs.SetMagThetaPhi(1, etaTRS, phiL);
	double thetaR, phiR; // phiR is value of phi @ estimated R in LORS
	populatePtrOuter->trs2Lors(dirTrs, thetaR, phiR);
	dirLORS.SetMagThetaPhi(1, thetaR, phiR);		
	vectorOut = populatePtrOuter->traceForward(dirLORS, lMean); 
}

void setMaxSingle(TVector2& vectorOut, double etaTRS, double phiL)
{
  const float lMin = 0.;
  const float lMean =  0.75;
	TVector3 dirTrs, dirLORS;
		
	// set TRS values :
	dirTrs.SetMagThetaPhi(1, etaTRS, phiL);
	double thetaR, phiR; // phiR is value of phi @ estimated R in LORS
	populatePtrOuter->trs2Lors(dirTrs, thetaR, phiR);
	dirLORS.SetMagThetaPhi(1, thetaR, phiR);		
	vectorOut = populatePtrOuter->traceForward(dirLORS, lMean); 
}


// Trace a single Ckov photon from emission point somewhere in radiator up to photocathode taking into account ref indexes of materials it travereses
// Arguments: ckovThe,ckovPhi- photon ckov angles in TRS, [rad]
//   Returns: distance between photon point on PC and track projection
// placeholder...
void setArrayMax(double etaTRS, vecArray2& inPutVectorPos, const size_t kN)
{
  // const size_t kN = inPutVector.size();
  const float lMin = 0.;
  const auto trkPC2 = populatePtrOuter->getPcImp();      // track at PC
  const auto trkRad2 = populatePtrOuter->getTrackPos();  // track at RAD


  //Printf("\n\n setArrayMax() enter --> kN %zu, etaTRS %.2f", kN, etaTRS);
	for(int i = 0; i < kN; i++){

		const auto phiL = double(TMath::TwoPi()*(i+1)/kN);
		TVector3 dirTrs, dirLORS;
		
		// set TRS values :
		dirTrs.SetMagThetaPhi(1, etaTRS, phiL);

		//Printf("\n setArrayMax() dirTrs.SetMagThetaPhi(1, etaTRS = x %.2f, phiL = x %.2f); ", etaTRS, phiL);
		

		double thetaR, phiR; // phiR is value of phi @ estimated R in LORS


		populatePtrOuter->trs2Lors(dirTrs, thetaR, phiR);

		dirLORS.SetMagThetaPhi(1, thetaR, phiR);
		
		const auto& max = populatePtrOuter->traceForward(dirLORS, lMin); 

  	inPutVectorPos.emplace_back(std::array<double, 2>{max.X(), max.Y()}); 


	}
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t aliSigma2(Double_t trkTheta,Double_t trkPhi,Double_t ckovTh, Double_t ckovPh)
{
// Analithical calculation of total error (as a sum of localization, geometrical and chromatic errors) on Cerenkov angle for a given Cerenkov photon 
// created by a given MIP. Fromulae according to CERN-EP-2000-058 
// Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]
//            dip and azimuthal angles for MIP taken at the entrance to radiator, [radians]        
//            MIP beta
//   Returns: absolute error on Cerenkov angle, [radians]    
  
  TVector3 v(-999,-999,-999);
  Double_t trkBeta = 1./(TMath::Cos(ckovTh)*GetRefIdx());
  
	if(TMath::Cos(ckovTh) )
  if (ckovTh > TMath::ASin(1. / GetRefIdx())) {  
    ckovTh = TMath::ASin(1. / GetRefIdx()) -0.01;
	}


  if(trkBeta > 1) trkBeta = 1;                 //protection against bad measured thetaCer  
  if(trkBeta < 0) trkBeta = 0.0001;            //

	double x = SigLoc(trkTheta, trkPhi, ckovTh, ckovPh, trkBeta);
	double y = SigGeom(trkTheta, trkPhi, ckovTh, ckovPh, trkBeta);
	double z = SigCrom(trkTheta, trkPhi, ckovTh, ckovPh, trkBeta);

	v.SetX(x);
	v.SetY(y);
	v.SetZ(z);


	Double_t sigRMS =  TMath::Sqrt(v.Mag2() + 0.002 * 0.002);


	if(sigRMS  < 0.00000002)
 		Printf("SigLoc %.5f, SigGeom %.5f, z (SigCrom): %.5f  == %.5f (sigRMS): | Inputs: trkTheta: %.5f, trkPhi: %.5f, ckovTh: %.5f, ckovPh: %.5f, trkBeta: %.5f", 
       x, y, z, sigRMS, trkTheta, trkPhi, ckovTh, ckovPh, trkBeta);


	return sigRMS;
  // adding 0.002 ckov from thetaP unc

}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t SigLoc(Double_t trkTheta,Double_t trkPhi,Double_t thetaC, Double_t phiC,Double_t betaM)
{
// Analitical calculation of localization error (due to finite segmentation of PC) on Cerenkov angle for a given Cerenkov photon 
// created by a given MIP. Fromulae according to CERN-EP-2000-058 
// Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]
//            dip and azimuthal angles for MIP taken at the entrance to radiator, [radians]        
//            MIP beta
//   Returns: absolute error on Cerenkov angle, [radians]    
  
  Double_t phiDelta = phiC - trkPhi;

  Double_t sint     = TMath::Sin(trkTheta);
  Double_t cost     = TMath::Cos(trkTheta);
  Double_t sinf     = TMath::Sin(trkPhi);
  Double_t cosf     = TMath::Cos(trkPhi);
  Double_t sinfd    = TMath::Sin(phiDelta);
  Double_t cosfd    = TMath::Cos(phiDelta);
  Double_t tantheta = TMath::Tan(thetaC);
  
  Double_t alpha =cost-tantheta*cosfd*sint;                                                 // formula (11)
  Double_t k = 1.-GetRefIdx()*GetRefIdx()+alpha*alpha/(betaM*betaM);        // formula (after 8 in the text)
  if (k<0) return 1e10;
  Double_t mu =sint*sinf+tantheta*(cost*cosfd*sinf+sinfd*cosf);                             // formula (10)
  Double_t e  =sint*cosf+tantheta*(cost*cosfd*cosf-sinfd*sinf);                             // formula (9)

  Double_t kk = betaM*TMath::Sqrt(k)/(GapThick()*alpha);                            // formula (6) and (7)
  Double_t dtdxc = kk*(k*(cosfd*cosf-cost*sinfd*sinf)-(alpha*mu/(betaM*betaM))*sint*sinfd); // formula (6)           
  Double_t dtdyc = kk*(k*(cosfd*sinf+cost*sinfd*cosf)+(alpha* e/(betaM*betaM))*sint*sinfd); // formula (7)            pag.4

  Double_t errX = 0.2,errY=0.25;                                                            //end of page 7
  return  TMath::Sqrt(errX*errX*dtdxc*dtdxc + errY*errY*dtdyc*dtdyc);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t SigCrom(Double_t trkTheta,Double_t trkPhi,Double_t thetaC, Double_t phiC,Double_t betaM)
{
// Analitical calculation of chromatic error (due to lack of knowledge of Cerenkov photon energy) on Cerenkov angle for a given Cerenkov photon 
// created by a given MIP. Fromulae according to CERN-EP-2000-058 
// Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]
//            dip and azimuthal angles for MIP taken at the entrance to radiator, [radians]        
//            MIP beta
//   Returns: absolute error on Cerenkov angle, [radians]    
  
  Double_t phiDelta = phiC - trkPhi;

  Double_t sint     = TMath::Sin(trkTheta);
  Double_t cost     = TMath::Cos(trkTheta);
  Double_t cosfd    = TMath::Cos(phiDelta);
  Double_t tantheta = TMath::Tan(thetaC);
  
  Double_t alpha =cost-tantheta*cosfd*sint;                                                 // formula (11)
  Double_t dtdn = cost*GetRefIdx()*betaM*betaM/(alpha*tantheta);                    // formula (12)
            
//  Double_t f = 0.00928*(7.75-5.635)/TMath::Sqrt(12.);
  Double_t f = 0.0172*(7.75-5.635)/TMath::Sqrt(24.);

  return f*dtdn;
}//SigCrom()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t SigGeom(Double_t trkTheta,Double_t trkPhi,Double_t thetaC, Double_t phiC,Double_t betaM)
{
// Analitical calculation of geometric error (due to lack of knowledge of creation point in radiator) on Cerenkov angle for a given Cerenkov photon 
// created by a given MIP. Formulae according to CERN-EP-2000-058 
// Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]
//            dip and azimuthal angles for MIP taken at the entrance to radiator, [radians]        
//            MIP beta
//   Returns: absolute error on Cerenkov angle, [radians]    

  Double_t phiDelta = phiC - trkPhi;

  Double_t sint     = TMath::Sin(trkTheta);
  Double_t cost     = TMath::Cos(trkTheta);
  Double_t sinf     = TMath::Sin(trkPhi);
  Double_t cosfd    = TMath::Cos(phiDelta);
  Double_t costheta = TMath::Cos(thetaC);
  Double_t tantheta = TMath::Tan(thetaC);
  
  Double_t alpha =cost-tantheta*cosfd*sint;                                                  // formula (11)
  
  Double_t k = 1.-GetRefIdx()*GetRefIdx()+alpha*alpha/(betaM*betaM);         // formula (after 8 in the text)
  if (k<0) return 1e10;

  Double_t eTr = 0.5*RadThick()*betaM*TMath::Sqrt(k)/(GapThick()*alpha);     // formula (14)
  Double_t lambda = (1.-sint*sinf)*(1.+sint*sinf);                                                  // formula (15)

  Double_t c1 = 1./(1.+ eTr*k/(alpha*alpha*costheta*costheta));                              // formula (13.a)
  Double_t c2 = betaM*TMath::Power(k,1.5)*tantheta*lambda/(GapThick()*alpha*alpha);  // formula (13.b)
  Double_t c3 = (1.+eTr*k*betaM*betaM)/((1+eTr)*alpha*alpha);                                // formula (13.c)
  Double_t c4 = TMath::Sqrt(k)*tantheta*(1-lambda)/(GapThick()*betaM);               // formula (13.d)
  Double_t dtdT = c1 * (c2+c3*c4);
  Double_t trErr = RadThick()/(TMath::Sqrt(12.)*cost);

  return trErr*dtdT;
}//SigGeom()



// Trace a single Ckov photon from emission point somewhere in radiator up to photocathode taking into account ref indexes of materials it travereses
// Arguments: ckovThe,ckovPhi- photon ckov angles in TRS, [rad]
//   Returns: distance between photon point on PC and track projection
// placeholder...
void setArrayMin(double etaTRS, vecArray2& inPutVectorPos, const size_t kN)
{
  // const size_t kN = inPutVector.size();
  const float lMax = 1.5;
  const auto trkPC2 = populatePtrInner->getPcImp();      // track at PC
  const auto trkRad2 = populatePtrInner->getTrackPos();  // track at RAD


  //Printf("\n\n setArrayMax() enter --> kN %zu, etaTRS %.2f", kN, etaTRS);
	for(int i = 0; i < kN; i++){

		const auto phiL = double(TMath::TwoPi()*(i+1)/kN);
		TVector3 dirTrs, dirLORS;
		
		// set TRS values :
		dirTrs.SetMagThetaPhi(1, etaTRS, phiL);

		//Printf("\n setArrayMax() dirTrs.SetMagThetaPhi(1, etaTRS = x %.2f, phiL = x %.2f); ", etaTRS, phiL);
		

		double thetaR, phiR; // phiR is value of phi @ estimated R in LORS


		populatePtrInner->trs2Lors(dirTrs, thetaR, phiR);
		
		//Printf("setArrayMax() called trs2Lors on dirTRS : return { thetaR = %.2f, phiR = %.2f}", thetaR, phiR);
		

		dirLORS.SetMagThetaPhi(1, thetaR, phiR);



		const auto& max = populatePtrInner->traceForward(dirLORS, lMax); 


    inPutVectorPos.emplace_back(std::array<double, 2>{max.X(), max.Y()});
		if((max.Y() == -999) or (max.X() == -999)) {
    	//inPutVector[i] = {0,0,0, 0, 0};
    	//inPutVector.emplace_back(std::array<double, 3>{0, 0, 0});
			// placeholder, find better solution? 
			if(max.Y() == -999) {
				//Printf("setArrayMax() max.Y() %.1f == -999", max.Y());
			}
			if(max.X() == -999) {
				//Printf("setArrayMax() max.X() %.1f == -999", max.X());
			}
    } else {
    	
    	//inPutVector.emplace_back(std::array<double, 3>{phiL, phiR, r});
    	
    	
    	// phiL : photon_phi i TRS system 
    	// phiR : photon_phi i LORS system
    	// phiPC : (photon - MIP).Phi();
    	
    	
    }

	}

}






// Trace a single Ckov photon from emission point somewhere in radiator up to photocathode taking into account ref indexes of materials it travereses
// Arguments: ckovThe,ckovPhi- photon ckov angles in TRS, [rad]
//   Returns: distance between photon point on PC and track projection
// placeholder...
void setArrayMin(double etaTRS, vecArray4& inPutVectorAngle, vecArray2& inPutVectorPos, const size_t kN)
{
  // const size_t kN = inPutVector.size();
  const float lMax = 1.5;
  const auto trkPC2 = populatePtrInner->getPcImp();      // track at PC
  const auto trkRad2 = populatePtrInner->getTrackPos();  // track at RAD


  //Printf("\n\n setArrayMax() enter --> kN %zu, etaTRS %.2f", kN, etaTRS);
	for(int i = 0; i < kN; i++){

		const auto phiL = double(TMath::TwoPi()*(i+1)/kN);
		TVector3 dirTrs, dirLORS;
		
		// set TRS values :
		dirTrs.SetMagThetaPhi(1, etaTRS, phiL);

		//Printf("\n setArrayMax() dirTrs.SetMagThetaPhi(1, etaTRS = x %.2f, phiL = x %.2f); ", etaTRS, phiL);
		

		double thetaR, phiR; // phiR is value of phi @ estimated R in LORS


		populatePtrInner->trs2Lors(dirTrs, thetaR, phiR);
		
		//Printf("setArrayMax() called trs2Lors on dirTRS : return { thetaR = %.2f, phiR = %.2f}", thetaR, phiR);
		

		dirLORS.SetMagThetaPhi(1, thetaR, phiR);

		//Printf("setArrayMax() dirLORS {x %.2f y %.2f z %.2f}", dirLORS.X(), dirLORS.Y(), dirLORS.Z());


		// temp
		// ckovThe, const double& ckovPhi, const double & L
		// this should return the same as max
		//const auto t = populatePtr->tracePhot(etaTRS, phiL, lMin);

		// temp
		
		//Printf("setArrayMax() called  populatePtrInner->traceForward(dirLORS (thetaR %.2f, phiR %.2f) lMin =  %.2f", thetaR, phiR, lMin);
		const auto& max = populatePtrInner->traceForward(dirLORS, lMax); 
		// max = pos of tracked photon at PC
		
		const auto r = (max - trkPC2).Mod(); // trkPC2 : track impact @ PC
				// from alinot_paattrec : this is MIP2photon distance?
				// const auto r = (max - trkPC2).Mod(); // trkPC2 : track impact @ PC
		
		
		
		auto phiPC = (max - trkPC2).Phi();   // MIP phi value/ PC track value?
		// vector fra photon til MIP/PC (--> verdi som kommer fra phiL)
		
		
		//Printf("setArrayMax() traceForward returned TVector2 {x %.2f y %.2f} == > R  = %.2f", max.X(), max.Y(), r);
		//Printf("setArrayMax() tracePhot returned TVector2 {x %.2f y %.2f} == > R  = %.2f", t.X(), t.Y(), (t-trkPC2).Mod());
		
		// add protection if traceForward returns 0 or -999?


 		//Printf("setArrayMax() --> i %d, {x %.2f y %.2f} - MIP {x %.2f y %.2f}", i, max.X(), max.Y(), trkPC.X(), trkPC.Y());
 		

		// check at phiR  er det samme som (t-rad).phi og (max-rad).phi?
		



		//Printf("setArrayMax2() --> i %d, {x %.2f y %.2f} - MIP {x %.2f y %.2f}", i, max.X(), max.Y(), trkPC2.X(), trkPC2.Y());


		// phiR in [-pi, pi]? set to 0..2pi?
		// inPutVector.emplace_back(std::array<double, 3>{phiL, phiR, r});
		if(phiR < 0) {
			phiR = TMath::TwoPi() + phiR;

    }
     
		if(phiPC < 0) {
			phiPC = TMath::TwoPi() + phiPC;

    }   

    
		/*if( (max-trkRad2).Phi() != (t-trkRad2).Phi()) {
			throw std::invalid_argument("(max-trkRad2).Phi() != (t-trkRad2).Phi())");
		}

		if(TMath::Abs(phiR - (t-trkRad2).Phi()) > 0.0001 && t.X() != -999 && t.Y() != -999)   
		{
		
		  Printf("phiR at %.6f | (t-rad) %.6f | (max-rad) %.6f ", phiR, (t-trkRad2).Phi(), (max-trkRad2).Phi());
			throw std::invalid_argument("phiR != (t-trkRad2).Phi())");
		}
		// protections if r > value?
		//if(r > )*/


		// TODO: if it goes out of map, find intersection with chamber-edges??
		// really jsut check wether TraceForward returned x, y = -999.
 		// to set points out of the map here is fine	


		if((max.Y() == -999) or (max.X() == -999)) {
    	//inPutVector[i] = {0,0,0, 0, 0};
    	//inPutVector.emplace_back(std::array<double, 3>{0, 0, 0});
			// placeholder, find better solution? 
			if(max.Y() == -999) {
				//Printf("setArrayMax() max.Y() %.1f == -999", max.Y());
			}
			if(max.X() == -999) {
				//Printf("setArrayMax() max.X() %.1f == -999", max.X());
			}
    } else {
    	
    	//inPutVector.emplace_back(std::array<double, 3>{phiL, phiR, r});
    	
    	
    	// phiL : photon_phi i TRS system 
    	// phiR : photon_phi i LORS system
    	// phiPC : (photon - MIP).Phi();
    	 
    	inPutVectorPos.emplace_back(std::array<double, 2>{max.X(), max.Y()}); 
    	inPutVectorAngle.emplace_back(std::array<double, 4>{phiL, phiR, phiPC, r});   	
			const auto r = (max - trkPC2).Mod(); // trkPC2 : track impact @ PC
			// from alinot_paattrec : this is MIP2photon distance?
			// const auto r = (max - trkPC2).Mod(); // trkPC2 : track impact @ PC

    	
    	
    	//inPutVectorAngle.emplace_back(std::array<double, 3>{phiL, phiR, r});
    	//Printf("setArrayMax() emplacing element %d : phiL %.2f, phiR %.2f, r %.2f", i, phiL, phiR, r);
    }
    	//inPutVector[i] = {phiL, phiR, r};
    // Printf("setArrayMax() emplacing element %d : phiL %.2f, phiR %.2f, r %.2f", i, phiL, phiR, r);

		/*if(maxPion.X() > 0 && maxPion.X() < 156.0 && maxPion.Y() > 0 && maxPion.Y() < 144) {
			// hMaxPion->Fill(maxPion.X(), maxPion.Y());
			// maxPionVec[i] = std::make_pair(maxPion.X(), maxPion.Y());
			Printf("maxPion loop i = %d, maxSize = %zu", i, maxPionVec.size()); 
		}*/
	}/*
  	Printf("\n");
	for(const auto& ip : inPutVectorAngle) {
		const auto& phiL_ = ip[0];
		const auto& phiR_ = ip[1]; 
		const auto& r_ = ip[2];  
		//Printf("setArrayMax() --> checking inputVector | : phiL %.2f, phiR %.2f, r %.2f", phiL_, phiR_, r_);
		if(r_ == 0 ) {throw std::invalid_argument("r====??????;");}
  } */
}



// map i.e., hMaxPion | vecArray i.e. arrMaxPion
 

// placeholder...
/*
void setArrayMin2(Populate* populate, double etaTRS, vecArray3& inPutVector)
{
  const size_t kN = inPutVector.size();
  const auto lMax = 1.5;
		
	for(int i = 0; i < kN; i++){

		const auto phiL = double(TMath::TwoPi()*(i+1)/kN);
		TVector3 dirTrs, dirLORS;
		dirTrs.SetMagThetaPhi(1, etaTRS, phiL);
		double thetaR, phiR; // phiR is value of phi @ estimated R in LORS
		populatePtr->trs2Lors(dirTrs, thetaR, phiR);
		dirLORS.SetMagThetaPhi(1, thetaR, phiR);

		const auto& max = populatePtr->traceForward(dirLORS, lMax);

		const auto r = (max - trkPC).Mod();
		inPutVector[i] = {phiL, phiR, r};


	}
} */ 


void populateRegions(std::vector<std::pair<double, double>>& vecArr, TH2F* map, const double& eta, const double& l) {	 
	 const int kN = vecArr.size();
	 for(int i = 0; i < vecArr.size(); i++){
		  const auto& value = populatePtr->tracePhot(eta, double(TMath::TwoPi()*(i+1)/kN), l);
		  if(/*value.X() > 0 && value.X() < 156.0 && value.Y() > 0 && value.Y() < 144*/true) {
		    map->Fill(value.X(), value.Y());
		    vecArr[i] = std::make_pair(value.X(), value.Y());
		  }
	 }
 		
 }


string getPDG(int pdg)
{
  	std::string pdgString;
    switch (TMath::Abs(pdg)) {
			case 11 : 
				pdgString = "Electron"; 
				break;
			case 211: 
				pdgString = "Pion"; 
				break;
			case 321: 
				pdgString = "Kaon"; 
				break;
			case 2212 : 
				pdgString = "Proton"; 
				break;
			case 50000050 : 
				pdgString = "Photon"; 
			case 50000051 : 
				pdgString = "Photon"; 
			case 22 : 
				pdgString = "Photon"; 

				break;
	  }
   	return pdgString;
 	}
    Double_t RadThick           (                                                                    ) const {return 1.5;}                                                        //Radiator thickness
    Double_t WinThick           (                                                                    ) const {return 0.5;}                                                        //Window thickness
    Double_t GapThick           (                                                                    ) const {return 8.0;}                                                        //Proximity gap thicknes
    Double_t GetRefIdx          (                                                                    ) const {return 1.2905;}                                                    //running refractive index
    Double_t WinIdx             (                                                                    ) const {return 1.583;}                                                     //Mean refractive index of WIN material (SiO2) 
    Double_t GapIdx             (                                                                    ) const {return 1.0005;}        







bool findPhotCkov(double cluX, double cluY, double& thetaCer, double& phiCer)
{
  // Finds Cerenkov angle  for this photon candidate
  // Arguments: cluX,cluY - position of cadidate's cluster
  // Returns: Cerenkov angle


 /*
TVector2 trkPos; // trk at RAD
TVector3 trkDir; // trk mag theta phi
TVector2 mipPos; // MIP PC
 */

  TVector3 dirCkov;

  double zRad = -0.5 * RadThick() - 0.5 * WinThick();     // z position of middle of RAD
  TVector3 rad(trkPos.X(), trkPos.Y(), zRad);                           // impact point at middle of RAD
  TVector3 pc(cluX, cluY, 0.5 * WinThick() +GapThick()); // mip at PC
  double cluR = TMath::Sqrt((cluX - mipPos.X()) * (cluX - mipPos.X()) +
                            (cluY - mipPos.Y()) * (cluY - mipPos.Y())); // ref. distance impact RAD-CLUSTER
  double phi = (pc - rad).Phi();                                  // phi of photon

  double ckov1 = 0;
  double ckov2 = 0.75 + thetaP; // start to find theta cerenkov in DRS
  const double kTol = 0.01;
  Int_t iIterCnt = 0;
  while (1) {
    if (iIterCnt >= 50) {
      return kFALSE;
    }
    double ckov = 0.5 * (ckov1 + ckov2);
    dirCkov.SetMagThetaPhi(1, ckov, phi);
    TVector2 posC = populatePtrInner->traceForward(dirCkov, 0.75);   // trace photon with actual angles
    double dist = cluR - (posC - mipPos).Mod(); // get distance between trial point and cluster position
    if (posC.X() == -999) {
      dist = -999;
    }           // total reflection problem
    iIterCnt++; // counter step
    if (dist > kTol) {
      ckov1 = ckov;
    } // cluster @ larger ckov
    else if (dist < -kTol) {
      ckov2 = ckov;
    }                                       // cluster @ smaller ckov
    else {                                  // precision achived: ckov in DRS found
      dirCkov.SetMagThetaPhi(1, ckov, phi); //
      populatePtrInner->lors2Trs(dirCkov, thetaCer, phiCer);  // find ckov (in TRS:the effective Cherenkov angle!)
      return kTRUE;
    }
  }
} // FindPhotTheta()

}; // end class CkovTools

#endif
