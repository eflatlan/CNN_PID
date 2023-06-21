#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <random>

#include <math.h>
//#include "ReconE.cpp"
#include "ReconG.cpp"

//namespace ParticleUtils
using namespace o2;
class CkovTools {


 private: 


  static constexpr double PI = M_PI;
  static constexpr double halfPI = M_PI/2;
  static constexpr double twoPI = M_PI*2;


  static constexpr double stdDevPion = 0.008; 
  static constexpr double stdDevKaon = 0.008; 
  static constexpr double stdDevProton = 0.008;
  static constexpr float tGap = 8;
  static constexpr float  rW = 1.; 
  static constexpr float  qW = 0.5;
  static constexpr float  L = rW/2;
  static constexpr float CH4GapWidth = 8;
  static constexpr float  RadiatorWidth = 1.;
  static constexpr float  QuartzWindowWidth = 0.5;
  static constexpr float  EmissionLenght = RadiatorWidth/2;
  float thetaP, phiP, xP, yP;
   float nF, nQ, nG;  
   double occupancy;
   std::array<float, 3> ckovHyps;
   std::vector<std::pair<double, double>> photons;

	 double ckovPionMin, ckovPionMax, ckovKaonMin, ckovKaonMax, ckovProtonMin,ckovProtonMax, mRMax, mL1Max, mL2Max, mRMin, mL1Min, mL2Min;

   float cosThetaP, sinThetaP, tanThetaP;
   float cosPhiP, sinPhiP, tanPhiP;
   float trackCkov;
   double xMipLocal, yMipLocal;
   float phiLCurrent, etaCCurrent;

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

  CkovTools (double xP, double yP, double thetaP, double phiP, 
             std::array<float, 3> ckovHyps, double nF, double nQ, double nG, double occupancy, float trackCkov)
    : xP(xP), yP(yP), thetaP(thetaP), phiP(phiP), 
      ckovHyps(ckovHyps), nF(nF), nQ(nQ), nG(nG), occupancy(occupancy) , trackCkov(trackCkov){
			  
		

			  ckovPionMin = ckovHyps[0] - 3 * stdDevPion;
			  ckovPionMax = ckovHyps[0] + 3 * stdDevPion;



			  ckovKaonMin = ckovHyps[1] - 3 * stdDevKaon;
			  ckovKaonMax = ckovHyps[1] + 3 * stdDevKaon;

			  ckovProtonMin = ckovHyps[2] - 3 * stdDevProton;
			  ckovProtonMax = ckovHyps[2] + 3 * stdDevProton;

      	cosThetaP = TMath::Cos(thetaP);
				sinThetaP = TMath::Sin(thetaP);
				tanThetaP = TMath::Tan(thetaP);

				cosPhiP = TMath::Cos(phiP);
				sinPhiP = TMath::Sin(phiP);


	TRotation rotZ; rotZ.RotateZ(phiP);
	TRotation rotY; rotY.RotateY(thetaP);

        TVector3 ip(0,0,rW-L+tGap+qW);
	TVector3 op; op = rotZ*rotY*ip;

				xMipLocal = tanThetaP*cosPhiP*(rW-L + tGap + qW);
				yMipLocal = tanThetaP*sinPhiP*(rW-L + tGap + qW);
	Printf("init Ckovtools \n MIP Root : %f %f %f \n MIP local %f %f",op.Px(),op.Py(),op.Pz(),xMipLocal,yMipLocal);
        // constructor body goes here, if needed


				mRMax = getR(ckovPionMax, halfPI);
				mL2Max = getR(ckovPionMax, 0);
				mL1Max = getR(ckovPionMax, PI);

				mRMin = getR(ckovProtonMin, halfPI);
        mL1Min = getR(ckovProtonMin, PI);
        mL2Min = getR(ckovProtonMin, 0);
			
				for(const auto& c : ckovHyps) {
					auto R = getR(c, halfPI);
					auto l2 = getR(c, 0);
					auto l1 = getR(c, PI);
					Printf("init CkovTools : CkovHyp %f, R %f, l1 %f l2 %f", c, R, l2, l1);
        }   
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
  std::vector<std::pair<double, double>> segment(std::vector<std::pair<double, double>>& cherenkovPhotons, MapType& bins) { 
     

    const auto infString = Form("localRef #Theta_{p}  = %.4f #Phi_{p} = %.4f \n #Theta_{C} = %.4f maxCkov = %.4f \n maxR = %.1f | l1 = %.1f | l2 = %.1f; x [cm]; y [cm]", thetaP,phiP,trackCkov,ckovPionMax, mRMax, mL1Max, mL2Max); 

    TH2F *localRefMIP = new TH2F("localRefMIP ", infString,400,-20.,20.,400,-20,20);
    TH2F *localRefMIPUnrot = new TH2F("localRefMIPUnrot ", infString,400,-20.,20.,400,-20,20);
    TH2F *localRefUnrot = new TH2F("localRef ", infString,400,-20.,20.,400,-20,20);
    TH2F *localRef = new TH2F("localRef ", infString,400,-20.,20.,400,-20,20);
  
    // TODO: ckovHyps : get std-dev for Theta_ckov of pion kaon and proton from the values theta_i

   // initialize recon with track input params 
   

   // TODO: change this, xMipLocal just placeholder
   // not sure if xPC simply is obtained like this
   double xPC = xMipLocal, yPC = yMipLocal;
   local2GlobalRef(xPC, yPC);


   ReconG reconG(thetaP, phiP, xP, yP, xPC, yPC, nF);



   TH2F *hNoiseMap = new TH2F("  Noise ", "  Noise ; x [cm]; y [cm]",160,0.,159.,144,0,143);
   int numPhotons= 0;
    for(const auto& p :cherenkovPhotons){

			const auto& xDif = p.first - xMipLocal; 
			const auto& yDif = p.second - yMipLocal; 
      auto R = TMath::Sqrt(xDif*xDif+yDif*yDif);
      Printf(" Ckovtools segments cherenkovPhotons : x %f y %f R = %f", p.first, p.second, R);
      numPhotons++;
    }

    TH2F *hSignalMap = new TH2F("Signal", Form("Cherenkov Photons = %d  #Theta_{p}  = %f #Phi_{p} = %f ; x [cm]; y [cm]", numPhotons, thetaP,phiP),160,0.,159.,144,0,143);   

    Printf("ckovtools segmetns cherenkovPhotons size  = %zu", cherenkovPhotons.size()); 
   
    MapType filledBins;
    // placeholders for the above : 

    // these are in local coordinate system
    // create bbox of minimal possible ckov
    // Proton with eta_c = theta_c_proton - 3*std_dev_proton

    double xLU = -mL1Max;
    double yLU = mRMax;

    double xRU = mL1Max;
    double yRU = mRMax;

    double xLD = -mL1Max;
    double yLD = -mL1Max;

    double xRD = mL1Max;
    double yRD = -mRMax;


		phiRing2Local(xLU, yLU, xMipLocal,yMipLocal);
		phiRing2Local(xRU, yRU, xMipLocal,yMipLocal);
		phiRing2Local(xRD, yRD, xMipLocal,yMipLocal);
		phiRing2Local(xLD, yLD, xMipLocal,yMipLocal);

    const auto coordsMinPhi0 = makeCkovPhoton(0., ckovProtonMin);
    const auto xMinPhi0 = coordsMinPhi0.first;
    const auto yMinPhi0 = coordsMinPhi0.second;

    // these are in local coordinate system
    const auto coordsMinPhiPi = makeCkovPhoton(static_cast<double>(3.1415), ckovProtonMin);
    
    const auto xMinPhiPi = coordsMinPhiPi.first;
    const auto yMinPhiPi = coordsMinPhiPi.second;

    const auto l1Min = mL1Min;//TMath::Sqrt((xMinPhiPi-  xMipLocal)*(xMinPhiPi-xMipLocal) + (yMinPhiPi-yMipLocal)*(yMinPhiPi-yMipLocal));
    const auto l2Min = mL2Min;//TMath::Sqrt((xMinPhi0-xMipLocal)*(xMinPhi0-xMipLocal) + (yMinPhi0-yMipLocal)*(yMinPhi0-yMipLocal));
    //const auto rMin = getR(ckovProtonMin);

    // these are in local coordinate system
    // create bbox of maximal possible ckov
    // Proton with eta_c = theta_c_pion + 3*std_dev_pion
    auto coordsMaxPhi0 = makeCkovPhoton(0, ckovPionMax);
    
    auto xMaxPhi0 = coordsMaxPhi0.first;
    auto yMaxPhi0 = coordsMaxPhi0.second;

    // these are in local coordinate system
    auto coordsMaxPhiPi = makeCkovPhoton(3.1415, ckovPionMax);
    auto xMaxPhiPi = coordsMaxPhiPi.first;
    auto yMaxPhiPi = coordsMaxPhiPi.second;
    
    local2PhiRing(xMaxPhiPi, yMaxPhiPi, xMipLocal,yMipLocal);
    local2PhiRing(xMaxPhi0, yMaxPhi0, xMipLocal,yMipLocal);
    double xMipRing = xMipLocal, yMipRing = yMipLocal;
    local2PhiRing(xMipRing, yMipRing, xMipRing, yMipRing); 
    //Printf("CkovTools segment : phiP %f thetaP %f xP %f yP %f ", phiP, thetaP, xP, yP);

    //Printf("CkovTools segment : xMaxPhiPi %f xMipLocal %f yMaxPhiPi %f yMipLocal %f ", xMaxPhiPi, xMipLocal, yMaxPhiPi, yMipLocal);

    // all retusn -nan?
    
    const auto l1Max = mL1Max;//TMath::Sqrt((xMaxPhiPi-xMipRing)*(xMaxPhiPi-xMipRing) + (yMaxPhiPi-yMipRing)*(yMaxPhiPi-yMipRing));
    const auto l2Max = mL2Max;//TMath::Sqrt((xMaxPhi0-xMipRing)*(xMaxPhi0-xMipRing) + (yMaxPhi0-yMipRing)*(yMaxPhi0-yMipRing)); 
    
    
      TBox* localBox = new TBox(xMaxPhi0, yMaxPhiPi - mRMax, xMaxPhiPi, yMaxPhiPi+mRMax);
      localBox->SetLineColor(kGreen);
      
      const auto& coords1 = local2Global(xMaxPhi0, yMaxPhi0 - mRMax);
      const auto& coords2 = local2Global(xMaxPhiPi, yMaxPhiPi - mRMax);
			const auto& coords3 = local2Global(xMaxPhi0, yMaxPhi0+mRMax);
			const auto& coords4 = local2Global(xMaxPhiPi, yMaxPhiPi+mRMax);
	    
	    TBox* globalBoxSignal = new TBox(coords1.first, coords1.second, coords2.first, coords2.second);
	    

			local2GlobalRef(xLU,yLU);
			local2GlobalRef(xRU,yRU);
			local2GlobalRef(xLD,yLD);
			local2GlobalRef(xRD,yRD);

	    TLine *tlineUpGlobal = new TLine(xLU, yLU, xRU, yRU);

	    
      TLine *tlineDownGlobal = new TLine(xLD, yLD, xRD, yRD);
      
		  tlineUpGlobal->SetLineColor(kGreen);
		  tlineDownGlobal->SetLineColor(kBlue);
	       
      globalBoxSignal->SetLineColor(kGreen);
      /*
      TLine *tlineUpLocal = new TLine(xMaxPhi0,yMaxPhi0 - mRMax, xMaxPhiPi, yMaxPhiPi - mRMax);
      TLine *tlineDownLocal = new TLine(xMaxPhi0,yMaxPhi0 + mRMax, xMaxPhiPi, yMaxPhiPi + mRMax);
			*/
      TLine *tlineUpLocal = new TLine(-mL1Max,mRMax, mL2Max, mRMax);

      TLine *tlineDownLocal = new TLine(-mL1Max,-mRMax, mL2Max, -mRMax);


			auto l1 = mL1Max; auto r1 = mRMax; auto r2 = mRMax; auto l2 = mL2Max;
			phiRing2Local(l1,r1,xMipLocal,yMipLocal);
			phiRing2Local(l2,r2,xMipLocal,yMipLocal);
      TLine *tlineUpLocalR = new TLine(-l1, r1,l2,r1);
      TLine *tlineDownLocalR = new TLine(-l1,-r2,l2,-r2);

		  tlineUpLocal->SetLineColor(kGreen);
		  tlineDownLocal->SetLineColor(kBlue);
		
		  tlineUpLocalR->SetLineColor(kGreen);
		  tlineDownLocalR->SetLineColor(kBlue);
      
       const auto rMax2 = getR3(ckovPionMax);
    //Printf("rMax2 = %f w ckov = %f", rMax2, ckovPionMax);
    // populate with background:
    const auto area = mRMax*2*(mL1Max+mL2Max);
    const auto numBackgroundPhotons = static_cast<int>(area*occupancy); 

    //Printf("CkovTools segment : rMax %f l1Max %f l2Max %f Area %f ", rMax, l1Max, l2Max, area);


   // Printf("CkovTools segment backGroundPhotons vector numBackgroundPhotons = %d", numBackgroundPhotons);
    std::vector<std::pair<double, double>> backGroundPhotons(numBackgroundPhotons);
     //   Printf("CkovTools segment backGroundPhotons vector created");
    std::random_device rd;
    std::mt19937 gen(rd());

    // in  phiRing ref sys
    std::uniform_real_distribution<> dis1(-mL1Max, mL2Max);
    std::uniform_real_distribution<> dis2(-mRMax, mRMax);

    for (auto& pair : backGroundPhotons) {
      pair.first = dis1(gen);
      pair.second = dis2(gen);
       // Printf("CkovTools segment created bg x%f y%f", pair.first, pair.second);
    }

    std::vector<std::pair<double, double>> photonCandidates;


    const auto infString2 = Form("globalREf #Theta_{p}  = %f #Phi_{p} = %f ; x [cm]; y [cm]", thetaP,phiP); 
    TH2F *globalREfMIP = new TH2F("globalREfMIP ", infString2,160,0.,159.,144,0,143);

		 const auto& coordsMIP = local2Global(xMipLocal, yMipLocal);
		 globalREfMIP->Fill(coordsMIP.first, coordsMIP.second);
   
      double xML = xMipLocal, yML = yMipLocal;
      localRefMIPUnrot->Fill(xMipLocal, yMipLocal);
      local2PhiRing(xML,yML,xML,yML);  
      localRefMIP->Fill(xML,yML);
    // here : loop through all backGroundPhotons and cherenkovPhotons
    for(const auto& photons : cherenkovPhotons) {  
      auto x = photons.first;
      auto y = photons.second;

      double xG = x, yG = y;
      localRefUnrot->Fill(x,y);

      // transform to phiRing ref-system
      local2PhiRing(x,y,xMipLocal,yMipLocal);
      localRef->Fill(x,y);

			/*
      Printf("ckovtools cherenkov photons x > xMaxPhiPi && x < xMaxPhi0 && y > yMaxPhiPi && y < yMaxPhi0");
      Printf("ckovtools cherenkov photons x  %f > xMaxPhiPi %f && x %f < xMaxPhi0 %f && y %f > yMaxPhiPi  %f && y %f < yMaxPhi0 %f",  x, xMaxPhiPi, x , xMaxPhi0 , y , yMaxPhiPi , y , yMaxPhi0);*/
      
		
	Printf("\nckovtools cherenkov photons x  %f > -mL1Max %f && x %f < mL2Max %f && y %f > -mRMax  %f && y %f < mRMax %f\n",  x, -mL1Max, x , mL2Max , y , -mRMax , y , mRMax);

      double thetaCer, phiCer;
      local2GlobalRef(xG, yG);
      // double cluX, double cluY, double& thetaCer, double& phiCer
      reconG.findPhotCkov(xG, yG, thetaCer, phiCer);	
      auto ckov = thetaCer;
      Printf("CkovTools segment thetaCer %f phiCer %f", thetaCer, phiCer);
	

      if(x > -mL1Max && x < mL2Max && y*y > mRMax*mRMax  && y < mRMax){
      //if(true){


	double thetaCer, phiCer;
	reconG.findPhotCkov(xG, yG, thetaCer, phiCer);	
	auto ckov = thetaCer;
        Printf("CkovTools segment thetaCer %f phiCer %f", thetaCer, phiCer);
	
        bool withinRange = true; 
        // check if the coordinates also corresponds to one of the possible cherenkov candidates


        // TODO : check if this method is wrong??
        // use here instead method from Recon.cxx
        //const auto& ckov = getCkovFromCoords(xP, yP, x, y, phiP, thetaP, nF, nQ, nG);      

       // Printf("CkovTools segment ckov %f", ckov);
        //Printf("CkovTools segment ckovPionMin %f ckovPionMax %f", ckovPionMin, ckovPionMax);
   
        // TODO: later, also add to candidates (i.e., pionCandidates, kaonCandidates...)
        if( ckovPionMin <  ckov & ckov < ckovPionMax ){
          withinRange = true;
        }
        
        if( ckovKaonMin <  ckov & ckov < ckovKaonMax ){
          withinRange = true;
        }
        
        if( ckovProtonMin <  ckov & ckov < ckovProtonMax ){
          withinRange = true;
        }
        
        if(withinRange){
          // transform to global coordinates:
          filledBins.push_back(std::make_pair(xG, yG));
          hSignalMap->Fill(xG, yG);

    	  //Printf("CkovTools segment : x%f y%f --> xG %f yG %f ", x,y, coords.first, coords.second); 
        } // end if withinRange
          
        // Fill map
        //hSignalAndNoiseMap->Fill(Xcen[n1], Ycen[n1]);
      } // end bbox if
    } // end for

    for(const auto& photons : backGroundPhotons) {  
      const auto& x = photons.first;
      const auto& y = photons.second;
        //Printf("CkovTools segment : backGroundPhotons %f x", x);
      
      //if(x > mL1Max && x < mL2Max && y > -mRMax && y < mRMax){
      if(x > -mL1Max && x < mL2Max && y*y > mRMax*mRMax  && y < mRMax){
        bool withinRange = true; 
        // check if the coordinates also corresponds to one of the possible cherenkov candidates

        // change to method from Recon.cxx:
        const auto& ckov = getCkovFromCoords(xP, yP, x, y, phiP, thetaP, nF, nQ, nG); 
        Printf("CkovTools segment :ckov%f ", ckov);
        // TODO: later, also add to candidates (i.e., pionCandidates, kaonCandidates...)
        if( ckovPionMin <  ckov & ckov < ckovPionMax ){
          withinRange = true;
        }
        
        if( ckovKaonMin <  ckov & ckov < ckovKaonMax ){
          withinRange = true;
        }
        
        if( ckovProtonMin <  ckov & ckov < ckovProtonMax ){
          withinRange = true;
        }
        
        if(withinRange){
          // transform to global coordinates:
          const auto coords = local2Global(x, y);
    	    //Printf("CkovTools segment : x%f y%f --> xG %f yG %f ", x,y, coords.first, coords.second);    
	        hNoiseMap->Fill(coords.first,coords.second);
          filledBins.push_back(coords);
          //Printf("CkovTools segment backGroundPhotons  x %f y %f", x,y);
        }      
      } // end if    
    } // end for
    
    //for(const auto& pair: filledBins)
    //	Printf("CkovTools segment candidates: x%f y%f", pair.first, pair.second);    

  hNoiseMap->SetMarkerColor(kRed);
    Printf("CkovTools segment filledBins Size %zu", filledBins.size());

  TCanvas *thSignalNoiseMap = new TCanvas("hSignalNoiseMap","hSignalNoiseMap",800,800);  // thSignalNoiseMap->Divide(2,1);
  thSignalNoiseMap->cd(1);
  //globalBoxSignal->Draw();

  globalREfMIP->SetMarkerStyle(3);
    globalREfMIP->SetMarkerColor(kRed);


  hSignalMap->SetMarkerStyle(2);
  hSignalMap->Draw();
  globalREfMIP->Draw("same");

  hNoiseMap->SetMarkerStyle(2);
  hNoiseMap->Draw("same");
	tlineUpGlobal->Draw();
	tlineDownGlobal->Draw();

  //hSignalAndNoiseMap->Show();
  /*thSignalNoiseMap->cd(2);
  //globalBoxSignal->Draw();
  hNoiseMap->Draw("same");*/
  thSignalNoiseMap->SaveAs("hSignalNoiseMap.png");
  thSignalNoiseMap->Show();
	tlineUpGlobal->Draw();
	tlineDownGlobal->Draw();

  TCanvas *thLocal = new TCanvas("thLocal","thLocal",800,800);  
  thLocal->cd();
  //localBox->Draw();
  //localRef->SetMarkerColor(kBlue);

  localRefMIP->SetMarkerStyle(3);
  localRef->SetMarkerStyle(2);
  localRef->Draw();
  localRefUnrot->SetMarkerStyle(2);
  localRefUnrot->SetMarkerColor(kRed);
  localRefUnrot->Draw("same");
  
  
  localRefMIPUnrot->SetMarkerColor(kRed);
  localRefMIPUnrot->SetMarkerStyle(3);
  localRefMIPUnrot->Draw("same");
    
  localRefMIP->Draw("same");
  tlineUpLocal->Draw();
	tlineDownLocal->Draw();
	  tlineUpLocalR->Draw();
	tlineDownLocalR->Draw();
	gPad->Update();

  return filledBins;
  } // end segment


  //r0G = rGL + rLP' ; rLP' = Rz(phi)[rLP]
  void phiRing2Local(double &xL, double &yL, const double& xMipL, const double& yMipL)
	{	  
	  TRotation mThetaRot;
	  mThetaRot.RotateZ(phiP);
	  
    TVector3 mip(xMipL, yMipL, 0);
	  TVector3 phiRingPos(xL, yL, 0);
	  TVector3 op = mip + mThetaRot*phiRingPos;
	  xL = op.Px();
	  yL = op.Py();
	}
    
  //r0G = rGL + rLP' ; rLP' = Rz(phi)[rLP]
  // rLP = Rz(-phi) * [r0G-rGL]
  void local2PhiRing(double &xL, double &yL, const double& xMipL, const double& yMipL)
	{	  
	  TRotation mThetaRot;
	  mThetaRot.RotateZ(-phiP);
	  
	  TVector3 pos(xL-xMipL, yL-yMipL, 0);
	  TVector3 op = mThetaRot*pos;
	  xL = op.Px();
	  yL = op.Py();
	}

	void local2GlobalRef(double& xL, double& yL)
	{	  
	  xL = xL  + xP;
	  yL = -yL  + yP;	  
	}


	std::pair<double, double> local2Global(double xL, double yL)
	{
	  
	  const auto x = xL  + xP;
	  const auto y = -yL  + yP;	  
	  return {x, y};
	}

	std::pair<double, double> global2Local(double xG, double yG)
	{
	
	  //mTheta.RotateY
	  
	  const auto x = xG - xP;
	  const auto y = - yG + yP;
	  return std::make_pair(x, y);
	}


 // get R at phiLocal = pi/2 V = 3pi/2
	double getR3(double etaC)
	{
		const auto cosEtaC = TMath::Cos(etaC);
		const auto sinEtaC = TMath::Sin(etaC);

		const auto rwlDeltaR = (rW - L)/(TMath::Sqrt(1-sinEtaC*sinEtaC));

		const auto qwDeltaR = (qW*nF)/(TMath::Sqrt(nQ*nQ-sinEtaC*sinEtaC*nF*nF));

		const auto tGapDeltaR = (tGap*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));

		const auto R = sinEtaC*(rwlDeltaR+qwDeltaR+tGapDeltaR)/cosThetaP;
    /*Printf("getR2 rwlDeltaR %f qwDeltaR %f tGapDeltaR %f", rwlDeltaR, qwDeltaR, tGapDeltaR);
    Printf("getR2 rW%f L %f qW %f tGap %f", rW, L, qW, tGap);
    Printf("getR2 nF %f nG %f nQ %f", nF, nG, nQ);*/
		return R;
	}

// get Radiu
	double getRadius(double etaC, double phiL)
	{
         
		const auto cosPhiL = TMath::Cos(phiL); // phiLocal = pi/2
		const auto sinPhiL = TMath::Sin(phiL); // --||-- 
		
		const auto cosEtaC = TMath::Cos(etaC);
		const auto sinEtaC = TMath::Sin(etaC);

		const auto rwlDeltaR = (rW - L)/(TMath::Sqrt(1-sinEtaC*sinEtaC));

		const auto qwDeltaR = (qW*nF)/(TMath::Sqrt(nQ*nQ-sinEtaC*sinEtaC*nF*nF));

		const auto num = (tGap + tanThetaP*cosPhiL*sinEtaC * (rwlDeltaR + qwDeltaR));


   // ef :error was on this line :
   // 		const auto denum = 1- (tanThetaP*cosPhiL*sinPhiP*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));

		const auto denum = 1 - (tanThetaP*cosPhiL*sinEtaC*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));

		const auto tZ = num/denum;

		const auto Lz = (rW-L) + qW + tZ;

		const auto tGapDeltaR = (tZ*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));



		const auto R = sinEtaC*(rwlDeltaR+qwDeltaR+tGapDeltaR)/cosThetaP;
		return R;
	}

// get R at phiLocal = pi/2 V = 3pi/2
	double getR2(double etaC)
	{
         
		const auto cosPhiL = 1; // phiLocal = pi/2
		const auto sinPhiL = 0; // --||-- 
		
		const auto cosEtaC = TMath::Cos(etaC);
		const auto sinEtaC = TMath::Sin(etaC);

		const auto rwlDeltaR = (rW - L)/(TMath::Sqrt(1-sinEtaC*sinEtaC));

		const auto qwDeltaR = (qW*nF)/(TMath::Sqrt(nQ*nQ-sinEtaC*sinEtaC*nF*nF));

		const auto num = (tGap + tanThetaP*cosPhiL*sinEtaC * (rwlDeltaR + qwDeltaR));


   // ef :error was on this line :
   // 		const auto denum = 1- (tanThetaP*cosPhiL*sinPhiP*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));

		const auto denum = 1 - (tanThetaP*cosPhiL*sinEtaC*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));

		const auto tZ = num/denum;

		const auto Lz = (rW-L) + qW + tZ;

		const auto tGapDeltaR = (tZ*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));



		const auto R = sinEtaC*(rwlDeltaR+qwDeltaR+tGapDeltaR)/cosThetaP;
    //Printf("getR : R %f |  wlGap %f qwDeltaR %f tGapDeltaR %f", R, rwlDeltaR,qwDeltaR,tGapDeltaR);
		return R;
	}



  // get radius from MIP to a specific phiL, etaC pair
	double getR(double etaC, double phiL)
	{
         
		const auto cosPhiL = TMath::Cos(phiL); // 
		const auto sinPhiL = TMath::Cos(phiL); // --||-- 
		
		const auto cosEtaC = TMath::Cos(etaC);
		const auto sinEtaC = TMath::Sin(etaC);

		const auto rwlDeltaR = (rW - L)/(TMath::Sqrt(1-sinEtaC*sinEtaC));
		const auto qwDeltaR = (qW*nF)/(TMath::Sqrt(nQ*nQ-sinEtaC*sinEtaC*nF*nF));
		const auto num = (tGap + tanThetaP*cosPhiL*sinEtaC * (rwlDeltaR + qwDeltaR));


   // ef :error was on this line :
   // 		const auto denum = 1- (tanThetaP*cosPhiL*sinPhiP*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));

		const auto denum = 1 - (tanThetaP*cosPhiL*sinEtaC*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));

		const auto tZ = num/denum;

		const auto Lz = (rW-L) + qW + tZ;

		const auto tGapDeltaR = (tZ*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));

		const auto R = sinEtaC*(rwlDeltaR+qwDeltaR+tGapDeltaR)/cosThetaP;
    //Printf("getR : R %f |  wlGap %f qwDeltaR %f tGapDeltaR %f", R, rwlDeltaR,qwDeltaR,tGapDeltaR);
		return R;
	}


  // populate map in local reference system
	// based on one etaC value (i.e., one photon at a time)
	const std::pair<double, double> makeCkovPhoton2(double phiL, double etaC)
	{

	auto thetaL=phiL;
	  //Printf("\nmakeCkovPhoton : phiL %f etaC %f", phiL, etaC);
		const auto cosPhiL = TMath::Cos(phiL);
		const auto sinPhiL = TMath::Sin(phiL);
		
    // if isNan; find etaC that corresponds to rwlDeltaR ^ qwDeltaR NOT nan
    //if(isNan){}
    
		const auto cosEtaC = TMath::Cos(etaC);
		const auto sinEtaC = TMath::Sin(etaC);

    auto thetaC = etaC;
    double firstTerm = ((rW - L) / std::cos(thetaP)) * std::tan(thetaC / 2.0);
    double secondTerm = (qW / std::cos(thetaP)) * std::tan(std::asin(std::sin(thetaC) * nF * nG / (nQ * nG)) / 2.0);
    double thirdTerm = ((rW - L + qW + tGap - rW + L - qW) / std::cos(thetaP)) * std::tan(std::asin((std::sin(thetaC) * nF) / nG) / 2.0);
		
		double R = firstTerm + secondTerm + thirdTerm;
		
		double x = R * std::cos(thetaL) * std::cos(phiP) * std::cos(thetaP) - R * std::sin(thetaL) * std::sin(thetaP);
		double y = R * std::cos(thetaL) * std::cos(phiP) * std::sin(thetaP) + R * std::sin(thetaL) * std::cos(thetaP);
		double z = R*std::cos(thetaL)*std::cos(phiP)*std::sin(thetaP) + R*std::sin(thetaL)*std::cos(thetaP);
			
		/*
		z = R*cos(thetaL)*sin(phiP) + l
		x = R*cos(thetaL)*cos(phiP)cos(thetaP) - Rsin(thetaL)sin(thetaP)
		y = Rcos(thetaL)*cos(phiP)sin(thetaP) + Rsin(thetaL)*cos(thetaP)*/
															
		//Printf("makeCkovPhoton : Lz %f tGapDeltaR %f T %f", Lz, tGapDeltaR, T);
		//Printf("makeCkovPhoton : x %f, y %f %z || z != %f \n", x, y, z, rW - L + tGap + qW);
				
		//		return std::make_pair(x, y);
		return std::make_pair(x,y);
	}


	// populate map in local reference system
	// based on one etaC value (i.e., one photon at a time)
	const std::pair<double, double> makeCkovPhoton4(double phiL, double etaC)
	{

	  //Printf("\nmakeCkovPhoton : phiL %f etaC %f", phiL, etaC);
		const auto cosPhiL = TMath::Cos(phiL);
		const auto sinPhiL = TMath::Sin(phiL);
		
    // if isNan; find etaC that corresponds to rwlDeltaR ^ qwDeltaR NOT nan
    //if(isNan){}
    
		const auto cosEtaC = TMath::Cos(etaC);
		const auto sinEtaC = TMath::Sin(etaC);

		const auto Lz = (rW-L) + qW + tGap;

		const auto rwlDeltaR = (rW - L)/(TMath::Sqrt(1-sinEtaC*sinEtaC));
		const auto qwDeltaR = (qW*nF)/(TMath::Sqrt(nQ*nQ-sinEtaC*sinEtaC*nF*nF));
		const auto tGapDeltaR = (tGap*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));

		const auto T = sinEtaC*(rwlDeltaR+qwDeltaR+tGapDeltaR);


    // stemte når denne var + ??
		auto z = - T*tanThetaP * cosPhiL + Lz;
    TRotation mZ1;
    mZ1.RotateZ(phiP);

    TRotation my;
    my.RotateY(thetaP);

    TRotation mZ2;
    mZ2.RotateZ(phiL);

    TVector3 vz(0,0,Lz); TVector3 vx(T,0,0);
    
    TVector3 pos; pos =  vz + mZ2 * vx;
    TVector3 posOut; posOut = mZ1*my*pos*(1/cosThetaP);
    
		const auto x = T*(cosPhiP*cosPhiL - sinPhiP*sinPhiL/cosThetaP) + Lz*tanThetaP*cosPhiP;

		const auto y = T*(sinPhiP*cosPhiL + cosPhiP*sinPhiL/cosThetaP) + Lz*tanThetaP*sinPhiP;

		/*
    Printf("makeCkovPhoton : Lz %f tGapDeltaR %f T %f", Lz, tGapDeltaR, T);

		Printf("makeCkovPhoton Rot x %f, y %f z %f || z != %f", posOut.Px(),posOut.Py(),posOut.Pz(), rW-L+tGap+qW);
		Printf("makeCkovPhoton Harcoded x %f, y %f z %f || z != %f \n", x,y,z, rW-L+tGap+qW);*/
		/*
		const auto xR = T*(cosPhiP*cosPhiL - sinPhiP*sinPhiL/cosThetaP) + Lz*tanThetaP*cosPhiP;

		const auto y = T*(sinPhiP*cosPhiL - cosPhiP*sinPhiL/cosThetaP) + Lz*tanThetaP*sinPhiP;*/
				                            

    //Printf("makeCkovPhoton : x %f, y %f || rot : x %f y %f \n", x, y, posOut.Px(), posOut.Py());
        
    //return std::make_pair(x, y);
		return std::make_pair(posOut.Px(), posOut.Py());
	}
  

	// populate map in local reference system
	// based on one etaC value (i.e., one photon at a time)
	const std::pair<double, double> makeCkovPhoton(double phiL, double etaC)
	{

	  Printf("\nmakeCkovPhoton : phiL %f etaC %f", phiL, etaC);
		const auto cosPhiL = TMath::Cos(phiL);
		const auto sinPhiL = TMath::Sin(phiL);
		
    // if isNan; find etaC that corresponds to rwlDeltaR ^ qwDeltaR NOT nan
    //if(isNan){}
    
		const auto cosEtaC = TMath::Cos(etaC);
		const auto sinEtaC = TMath::Sin(etaC);


		const auto rwlDeltaR = (rW - L)/(TMath::Sqrt(1-sinEtaC*sinEtaC));

		const auto qwDeltaR = (qW*nF)/(TMath::Sqrt(nQ*nQ-sinEtaC*sinEtaC*nF*nF));


		const auto num = (tGap + tanThetaP*cosPhiL*sinEtaC*(rwlDeltaR + qwDeltaR));

		const auto denum = 1 - (tanThetaP*cosPhiL*sinEtaC*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));


    Printf("makeCkovPhoton : rwlDeltaR %f qwDeltaR %f num %f", rwlDeltaR, qwDeltaR, num);

		const auto tZ = num/denum;


		const auto Lz = (rW-L) + qW + tZ;


		const auto tGapDeltaR = (tZ*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));

		const auto T = sinEtaC*(rwlDeltaR+qwDeltaR+tGapDeltaR);


    // stemte når denne var + ??
		auto z = - T*tanThetaP * cosPhiL + Lz;
    TRotation mZ1;
    mZ1.RotateZ(phiP);

    TRotation my;
    my.RotateY(thetaP);

    TRotation mZ2;
    mZ2.RotateZ(phiL);

    TVector3 vz(0,0,Lz); TVector3 vx(T,0,0);
    
    TVector3 pos; pos =  vz + mZ2 * vx;
    TVector3 posOut; posOut = mZ1*my*pos*(1/cosThetaP);
    
		const auto x = T*(cosPhiP*cosPhiL - sinPhiP*sinPhiL/cosThetaP) + Lz*tanThetaP*cosPhiP;

		const auto y = T*(sinPhiP*cosPhiL + cosPhiP*sinPhiL/cosThetaP) + Lz*tanThetaP*sinPhiP;


    Printf("makeCkovPhoton : Lz %f tGapDeltaR %f T %f", Lz, tGapDeltaR, T);

		Printf("makeCkovPhoton Rot x %f, y %f z %f || z != %f", posOut.Px(),posOut.Py(),posOut.Pz(), rW-L+tGap+qW);
		Printf("makeCkovPhoton Harcoded x %f, y %f z %f || z != %f \n", x,y,z, rW-L+tGap+qW);
		/*
		const auto xR = T*(cosPhiP*cosPhiL - sinPhiP*sinPhiL/cosThetaP) + Lz*tanThetaP*cosPhiP;

		const auto y = T*(sinPhiP*cosPhiL - cosPhiP*sinPhiL/cosThetaP) + Lz*tanThetaP*sinPhiP;*/
				                            

    //Printf("makeCkovPhoton : x %f, y %f || rot : x %f y %f \n", x, y, posOut.Px(), posOut.Py());
        
    //return std::make_pair(x, y);
		return std::make_pair(posOut.Px(), posOut.Py());
	}


	double getCkovFromCoords(double xP, double yP, double xL, double yL, double phiP, double thetaP, float nF, float nQ, float nG)
	{       
    const auto& coords = local2Global(xL, yL);
    const auto x = coords.first;
		const auto y = coords.second;
		double phiF = 0;
		double thetaF1,thetaF2,thetaF=0,thetaLimite;
		float xPi=0,yPi=0,xF=0,yF=0,xF1=0,yF1=0,xF2=0,yF2=0; 
		float ThetaCherenkov, PhiCherenkov, DegPhiCherenkov;

		
		float deltaX = (rW-L+qW+tGap)*tanThetaP*cosPhiP;
		float deltaY = (rW-L+qW+tGap)*tanThetaP*sinPhiP;
			
		xPi = xP - deltaX;
		yPi = yP - deltaY;

		TVector3 v2(x-xPi-L*tanThetaP*cosPhiP, y-yPi-L*tanThetaP*sinPhiP,rW+qW+tGap-L); 

		phiF = v2.Phi();      

		thetaLimite = TMath::ASin(nG/nQ);

		double thetaF0 = TMath::ASin(nQ/nF*TMath::Sin(thetaLimite))-0.00001;

		//  Printf("thetaF0 = %f",thetaF0*TMath::RadToDeg());

		double thetaF01 = TMath::ASin((nF/nQ)*(TMath::Sin(thetaF0)));      

		double thetaF02 = TMath::ASin((nQ/nG)*(TMath::Sin(thetaF01)));

		float x01 = L*tanThetaP*cosPhiP;

		float y01 =  L*tanThetaP*sinPhiP;

		float x02 = (rW - L)*TMath::Tan(thetaF0)*TMath::Cos(phiF)+qW*TMath::Tan(thetaF01)*TMath::Cos(phiF)+tGap*TMath::Tan(thetaF02)*TMath::Cos(phiF);

		float y02 = (rW - L)*TMath::Tan(thetaF0)*TMath::Sin(phiF) + qW*TMath::Tan(thetaF01)*TMath::Sin(phiF) + tGap*TMath::Tan(thetaF02)*TMath::Sin(phiF);  

		float x0 = x01 + x02;
		float y0 = y01 + y02;

		double ThetaMin = 0;
		//  Double_t ThetaMax = 0.75+thetaP;
		double ThetaMax = thetaF0;

		xF = 999;
		yF = 999;

		Int_t nWhile = 0;

		while(TMath::Sqrt((xF-x+xPi)*(xF-x+xPi)+(yF-y+yPi)*(yF-y+yPi))>0.0001)
		{ 
		  nWhile++;

		  thetaF = static_cast<double>( (0.5*(ThetaMax - ThetaMin) + ThetaMin));

		  thetaF1 = TMath::ASin((nF/nQ)*(TMath::Sin(thetaF)));     
		  thetaF2 = TMath::ASin((nQ/nG)*(TMath::Sin(thetaF1)));

		  xF1 = EmissionLenght*tanThetaP*cosPhiP;
		  yF1 =  EmissionLenght*tanThetaP*sinPhiP;

		  xF2 = (rW-L)*TMath::Tan(thetaF)*TMath::Cos(phiF)+qW*TMath::Tan(thetaF1)*TMath::Cos(phiF)+tGap*TMath::Tan(thetaF2)*TMath::Cos(phiF);
		  yF2 = (rW-L)*TMath::Tan(thetaF)*TMath::Sin(phiF) + qW*TMath::Tan(thetaF1)*TMath::Sin(phiF) + tGap*TMath::Tan(thetaF2)*TMath::Sin(phiF);  

		  xF = xF1 + xF2;
		  yF = yF1 + yF2;
		  
		  if(TMath::Sqrt((xF-x0)*(xF-x0)+(yF-y0)*(yF-y0))>TMath::Sqrt((x-xPi-x0)*(x-xPi-x0)+(y-yPi-y0)*(y-yPi-y0)))
		  {
		    ThetaMin = thetaF;
		  }

		  else 
		  {
		    ThetaMax = thetaF;   
		  }  

		  if(nWhile>30) break;

		} // while 
		      
			    
		TVector3 vP((TMath::Sin(thetaP))*(TMath::Cos(phiP)),(TMath::Sin(thetaP))*(TMath::Sin(phiP)),(TMath::Cos(thetaP)));
		TVector3 vz(0.,0.,1.);

		TVector3 v1 = vP.Cross(vz);

		TVector3 vF((TMath::Sin(thetaF))*(TMath::Cos(phiF)),(TMath::Sin(thetaF))*(TMath::Sin(phiF)),(TMath::Cos(thetaF)));

		
		
		if(thetaP==0)
		{	      
		  ThetaCherenkov = thetaF;		  
		  PhiCherenkov = phiF;	      
		}	  
		else		
		{
		  vF.Rotate(thetaP,v1);
		  ThetaCherenkov = vF.Theta();
		  PhiCherenkov = vF.Phi();
		}

		DegPhiCherenkov = 180*PhiCherenkov/(TMath::Pi());

		if(DegPhiCherenkov<0) DegPhiCherenkov+=360;
    return static_cast<double>(ThetaCherenkov);	
	}

};

