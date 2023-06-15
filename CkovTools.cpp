#include <iostream>
#include <cmath>
#include <random>




//namespace ParticleUtils

class CkovTools {


 private: 
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

   float cosThetaP, sinThetaP, tanThetaP;
   float cosPhiP, sinPhiP, tanPhiP;
   float xMipLocal, yMipLocal;
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
             std::array<float, 3> ckovHyps, double nF, double nQ, double nG, double occupancy)
    : xP(xP), yP(yP), thetaP(thetaP), phiP(phiP), 
      ckovHyps(ckovHyps), nF(nF), nQ(nQ), nG(nG), occupancy(occupancy) {
      	cosThetaP = TMath::Cos(thetaP);
				sinThetaP = TMath::Sin(thetaP);
				tanThetaP = TMath::Tan(thetaP);

				cosPhiP = TMath::Cos(phiP);
				sinPhiP = TMath::Sin(phiP);

				xMipLocal = tanThetaP*cosPhiP;
				yMipLocal = tanThetaP*sinPhiP;
        // constructor body goes here, if needed
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
    // TODO: ckovHyps : get std-dev for Theta_ckov of pion kaon and proton from the values theta_i
    

    MapType filledBins;
    // placeholders for the above : 
    double stdDevPion = 0.08; 
    double stdDevKaon = 0.08; 
    double stdDevProton = 0.08;
    double ckovPionMin = ckovHyps[0] - 3 * stdDevPion;
    double ckovPionMax = ckovHyps[0] + 3 * stdDevPion;

    double ckovKaonMin = ckovHyps[1] - 3 * stdDevKaon;
    double ckovKaonMax = ckovHyps[1] + 3 * stdDevKaon;

    double ckovProtonMin = ckovHyps[2] - 3 * stdDevProton;
    double ckovProtonMax = ckovHyps[2] + 3 * stdDevProton;

    // these are in local coordinate system
    // create bbox of minimal possible ckov
    // Proton with eta_c = theta_c_proton - 3*std_dev_proton
    const auto coordsMinPhi0 = makeCkovPhoton(0., ckovProtonMin);
    const auto xMinPhi0 = coordsMinPhi0.first;
    const auto yMinPhi0 = coordsMinPhi0.first;

    // these are in local coordinate system
    const auto coordsMinPhiPi = makeCkovPhoton(static_cast<double>(3.1415), ckovProtonMin);
    const auto xMinPhiPi = coordsMinPhiPi.first;
    const auto yMinPhiPi = coordsMinPhiPi.first;

    const auto l1Min = TMath::Sqrt((xMinPhiPi-  xMipLocal)*(xMinPhiPi-xMipLocal) + (yMinPhiPi-yMipLocal)*(yMinPhiPi-yMipLocal));
    const auto l2Min = TMath::Sqrt((xMinPhi0-xMipLocal)*(xMinPhi0-xMipLocal) + (yMinPhi0-yMipLocal)*(yMinPhi0-yMipLocal));
    const auto rMin = getR(ckovProtonMin);

    // these are in local coordinate system
    // create bbox of maximal possible ckov
    // Proton with eta_c = theta_c_pion + 3*std_dev_pion
    const auto coordsMaxPhi0 = makeCkovPhoton(0, ckovPionMax);
    const auto xMaxPhi0 = coordsMaxPhi0.first;
    const auto yMaxPhi0 = coordsMaxPhi0.first;

    // these are in local coordinate system
    const auto coordsMaxPhiPi = makeCkovPhoton(3.1415, ckovPionMax);
    const auto xMaxPhiPi = coordsMaxPhiPi.first;
    const auto yMaxPhiPi = coordsMaxPhiPi.first;
    Printf("CkovTools segment : phiP %f thetaP %f xP %f yP %f ", phiP, thetaP, xP, yP);

    Printf("CkovTools segment : xMaxPhiPi %f xMipLocal %f yMaxPhiPi %f yMipLocal %f ", xMaxPhiPi, xMipLocal, yMaxPhiPi, yMipLocal);

    // all retusn -nan?
    const auto l1Max = TMath::Sqrt((xMaxPhiPi-xMipLocal)*(xMaxPhiPi-xMipLocal) + (yMaxPhiPi-yMipLocal)*(yMaxPhiPi-yMipLocal));
    const auto l2Max = TMath::Sqrt((xMaxPhi0-xMipLocal)*(xMaxPhi0-xMipLocal) + (yMaxPhi0-yMipLocal)*(yMaxPhi0-yMipLocal));
    const auto rMax = getR(ckovPionMax);

    // populate with background:
    const auto area = rMax*2*(l1Max+l2Max);
    const auto numBackgroundPhotons = area*occupancy; 

    Printf("CkovTools segment : rMax %f l1Max %f l2Max %fArea %f ", rMax, l1Max, l2Max, area);


    Printf("CkovTools segment backGroundPhotons vector numBackgroundPhotons = %f", numBackgroundPhotons);
    std::vector<std::pair<double, double>> backGroundPhotons(numBackgroundPhotons);
        Printf("CkovTools segment backGroundPhotons vector created");
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis1(0, l1Max + l2Max);
    std::uniform_real_distribution<> dis2(-rMax, rMax);

    for (auto& pair : backGroundPhotons) {
      pair.first = dis1(gen);
      pair.second = dis2(gen);
    }

    std::vector<std::pair<double, double>> photonCandidates;

    // here : loop through all backGroundPhotons and cherenkovPhotons
    for(const auto& photons : cherenkovPhotons) {  
      const auto& x = photons.first;
      const auto& y = photons.second;
      if(x > xMaxPhiPi && x < xMaxPhi0 && y > yMaxPhiPi && y < yMaxPhi0){
        bool withinRange = false; 
        // check if the coordinates also corresponds to one of the possible cherenkov candidates
        const auto& ckov = getCkovFromCoords(xP, yP, x, y, phiP, thetaP, nF, nQ, nG);      

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
          const auto& coords = local2Global(x, y);
          filledBins.push_back(coords);
        }
          
        // Fill map
        //hSignalAndNoiseMap->Fill(Xcen[n1], Ycen[n1]);
      } // end bbox if
    } // end for

    for(const auto& photons : backGroundPhotons) {  
      const auto& x = photons.first;
      const auto& y = photons.second;
      if(x > xMaxPhiPi && x < xMaxPhi0 && y > yMaxPhiPi && y < yMaxPhi0){
        bool withinRange = false; 
        // check if the coordinates also corresponds to one of the possible cherenkov candidates
        const auto& ckov = getCkovFromCoords(xP, yP, x, y, phiP, thetaP, nF, nQ, nG); 
             
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
          filledBins.push_back(coords);
        }      
      } // end if    
    } // end for
    return filledBins;
  } // end segment

    
	std::pair<double, double> local2Global(double xL, double yL)
	{
	  const auto x = cosPhiP*xL + sinPhiP*yL - xP;
	  const auto y = sinPhiP*xL - cosPhiP*yL - yP;
	  return {x, y};
	}

	std::pair<double, double> global2Local(double xG, double yG)
	{
	  const auto x = cosPhiP*xG - sinPhiP*yG + xP;
	  const auto y = sinPhiP*xG + cosPhiP*yG + yP;
	  return std::make_pair(x, y);
	}




  // get R at phiLocal = pi/2 V = 3pi/2
	double getR(double etaC)
	{

		const auto cosPhiL = 0; // phiLocal = pi/2
		const auto sinPhiL = 1; // --||-- 
		
		const auto cosEtaC = TMath::Cos(etaC);
		const auto sinEtaC = TMath::Sin(etaC);

		const auto rwlGap = (rW - L)/(TMath::Sqrt(1-sinEtaC*sinEtaC));

		const auto qwGap = (qW*nF)/(TMath::Sqrt(nQ*nQ-sinEtaC*sinEtaC*nF*nF));

		const auto numerator = (tGap + tanThetaP*cosPhiL*sinPhiP * (rwlGap + qwGap));

		const auto denominator = 1- (tanThetaP*cosPhiL*sinPhiP*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));

		const auto tGapZ = numerator/denominator;


		const auto Lz = (rW-L) + qW + tGapZ;


		const auto tGapGap = (tGapZ*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));

		const auto R = sinEtaC*(rwlGap+qwGap+tGapGap)/cosThetaP;

		return R;
	}


  

	// populate map in local reference system
	// based on one etaC value (i.e., one photon at a time)
	const std::pair<double, double> makeCkovPhoton(double phiL, double etaC)
	{
		const auto cosPhiL = TMath::Cos(phiL);
		const auto sinPhiL = TMath::Sin(phiL);
		
		const auto cosEtaC = TMath::Cos(etaC);
		const auto sinEtaC = TMath::Sin(etaC);


		const auto rwlGap = (rW - L)/(TMath::Sqrt(1-sinEtaC*sinEtaC));

		const auto qwGap = (qW*nF)/(TMath::Sqrt(nQ*nQ-sinEtaC*sinEtaC*nF*nF));

		const auto numerator = (tGap + tanThetaP*cosPhiL*sinPhiP * (rwlGap + qwGap));

		const auto denominator = 1- (tanThetaP*cosPhiL*sinPhiP*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));

		const auto tGapZ = numerator/denominator;


		const auto Lz = (rW-L) + qW + tGapZ;


		const auto tGapGap = (tGapZ*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));

		const auto T = sinEtaC*(rwlGap+qwGap+tGapGap);

		const auto x = T*(cosPhiP*cosPhiL - sinPhiP*sinPhiL/cosThetaP) + Lz*tanThetaP*cosPhiP;

		const auto y = T*(sinPhiP*cosPhiL - cosPhiP*sinPhiL/cosThetaP) + Lz*tanThetaP*sinPhiP;

		return std::make_pair(x, y);
	}


	double getCkovFromCoords(double xP, double yP, double x, double y, double phiP, double thetaP, float nF, float nQ, float nG)
	{    
		double phiF = 0;
		double thetaF1,thetaF2,thetaF=0,thetaLimite;
		float xPi=0,yPi=0,xF=0,yF=0,xF1=0,yF1=0,xF2=0,yF2=0; 
		float ThetaCherenkov, PhiCherenkov, DegPhiCherenkov;

		
		float deltaX = (rW+qW+tGap-L)*TMath::Tan(thetaP)*TMath::Cos(phiP);
		float deltaY = (rW+qW+tGap-L)*TMath::Tan(thetaP)*TMath::Sin(phiP);
			
		xPi = xP - deltaX;
		yPi = yP - deltaY;

		TVector3 v2(x-xPi-L*TMath::Tan(thetaP)*TMath::Cos(phiP), y-yPi-L*TMath::Tan(thetaP)*TMath::Sin(phiP),rW+qW+tGap-L); 

		phiF = v2.Phi();      

		thetaLimite = TMath::ASin(nG/nQ);

		double thetaF0 = TMath::ASin(nQ/nF*TMath::Sin(thetaLimite))-0.00001;

		//  Printf("thetaF0 = %f",thetaF0*TMath::RadToDeg());

		double thetaF01 = TMath::ASin((nF/nQ)*(TMath::Sin(thetaF0)));      

		double thetaF02 = TMath::ASin((nQ/nG)*(TMath::Sin(thetaF01)));

		float x01 = EmissionLenght*TMath::Tan(thetaP)*TMath::Cos(phiP);

		float y01 =  EmissionLenght*TMath::Tan(thetaP)*TMath::Sin(phiP);

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

		  xF1 = EmissionLenght*TMath::Tan(thetaP)*TMath::Cos(phiP);
		  yF1 =  EmissionLenght*TMath::Tan(thetaP)*TMath::Sin(phiP);

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

