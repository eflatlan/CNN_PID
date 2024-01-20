#include <iostream>
#include <cmath>
#include <random>
#include "populate.cpp"
#include <math.h>

// populate.cpp
#include <TVector2.h>
#include <TVector3.h>
#include <TRotation.h>
#include <TRandom.h>
#include <vector>

#include "Math/GenVector/Rotation3D.h"
#include "Math/GenVector/Cartesian3D.h"
#include "Math/GenVector/DisplacementVector3D.h"

#include "Math/GenVector/RotationYfwd.h"
#include "Math/GenVector/RotationZfwd.h"


class getR {


   using Polar3D = ROOT::Math::Polar3D<double>;
   using Rotation3D = ROOT::Math::Rotation3D;
   using RotationY = ROOT::Math::RotationY;
   using RotationZ = ROOT::Math::RotationZ;


private:

   
    TVector2 fTrkPos; // track pos in LORS at RAD // xRad, yRad
    TVector3 fTrkDir; // track position in LORS at RAD // setMagThetaPhi(1, thetaP, phiP)
    TVector2 fPc; // track pos at PC

    TVector2 fMipPos; // MIP-pos at PC

public:

	void checkCond(const TVector2& posPhoton, bool above, const std::vector<std::array<double>, 3>& vec) {
		// vec : contains phiL, phi, R w etaMin/etaMax for hadron species

		
		double phi1L = vec[index][0];
		double phi2L = vec[index+1][0];
		
		double phi1 = vec[index][1];
		double phi2 = vec[index+1][1];

		// get radiuses
		double r1 = vec[index][2];
		double r2 = vec[index+1][2];
		
		auto phiC = phiPhoton - phiP;
		switch (phiC) {
			case phiC < TMath::Pi()/2:
				initValue = 0.;
			case phiC < TMath::Pi():
				initValue = kN/4;
			case phiC < (3./2.)*TMath::Pi():
				initValue = kN/2.;
			case phiC < TMath::TwoPi():
				initValue = (3./2.)*kN;
			default:
				initValue = 0;
				Printf("phiC %.2f", phiC);
				throw std::invalid_argument("wtf value does phiC have?");
		}
		int iCnt = 0;
		phi1 = maxPionVec[initValue]; // 


		// set increment opposite way if phi1 > phiPhoton
		int inc = 1;
		if phi1 > phiPhoton
			inc = -1;		

		// accesing the correct phi
		while(phi1 < phiPhoton) {
			phi1 = vec[initValue + iCnt][1];
			iCnt += inc;
		}

		// TODO: this has to be changed if inc = -1?
		int index = iCnt -1;

		// correct indexes are found 
		Printf("phiPhoton %.2f| phiPhoton -  1 %.2f, phiPhoton + 1 %.2f", phiPhoton, vec[initValue + iCnt][1], vec[initValue + iCnt+ 1][1]);
		
				

		// eta is min/max ckovHyps +- 3 std-dev

		if(r2 > r1) {
			rMin = r1;
			rMax = r2;
			phiMax = phi2;
			phiMin = phi1;
			phiLmax = phiL2;
			phiLmin = phiL1;
		}	

		else {
			rMin = r2;
			rMax = r1;
			phiMax = phi1;
			phiMin = phi2;
			phiLmax = phiL1;
			phiLmin = phiL2;
		}


		// check if under a radius 
		if(getAbove) {
			while(noDecisionTaken) {
				if((rPhoton > rMax)) {
					// stop iterating, condition is false
					noDecisionTaken = true;
					condition = false;
				} else if((rPhoton < rMin)) {
					// stop iterating, condition is ok!
					noDecisionTaken = true;
					condition = true;
				} else if((rPhoton > rMin && rPhoton < rMax)) {
					// iterate by splitting r1, r2 -> phi1 phi2
				
					// find new value for eta1, eta2; (passed by ref)
					// phiL1, phiL2 also passed by ref

					splitPhi(phiMax, phiLmax, phiMin, phiLmin, etaCkov, phiPhoton); // 
					Printf(=
					rMax = getR(etaCkov, phiLmax, L);
					rMin = getR(etaCkov, phiLmin, L);
				} else {
					throw std::invalid_argument("noDecisionTaken???");
				}
			}
		} 

		// check if outside of radius
		else {
			while(noDecisionTaken) {
				if((rPhoton > rMax)) {
					// stop iterating, condition is true
					noDecisionTaken = true;
					condition = true;
				} else if((rPhoton < rMin)) {
					// stop iterating, condition is false!
					noDecisionTaken = true;
					condition = false;
				} else if((rPhoton > rMin && rPhoton < rMax)) {
					// iterate by splitting r1, r2 -> phi1 phi2

					splitPhi(phiMax, phiLmax, phiMin, phiLmin, etaCkov, phiPhoton); // 
					rMax = getR(etaCkov, phiLmax, L);
					rMin = getR(etaCkov, phiLmin, L);
				} else {
					throw std::invalid_argument("noDecisionTaken???");
				}
			}
		} 
	}

	double getR(const double& etaTRS, const double& phiTRS, const double& L)
	{
		//TVector3 dirCkov;
		//dirCkov.SetMagThetaPhi(1, etaTRS, phiTRS);
		const Tvector2 rPosLORS = tracePhot(etaTRS, phiTRS, L); // pos of estimated rPos at LORS PC
		const auto R = (rPosLORS-fPC).Mod(); 
		return R;
	}
		
	void splitPhi(double& phiMax, double& phiLmax, double& phiMin, double& phiLmin, const double& etaCkov, const double& L, const double& phiPhoton)
	{
		const auto phiNewL = (phiMaxL + phiMinL)/2.0;	// phi in TRS


		// Get phi in LORS 
		Tvector3 dirNewTrs;
		dirNewTrs->SetMagThetaPhi(1, etaCkov, phiNewL); 
		
		
		phiNew = trs2Lors(dirNewTrs, thetaP, phiP); 
		
		
		if((phiPhoton) < TMath::Pi() + phiP){
			if((phiNew-phiP) > (phiPhoton-phiP)) {
				phiMin = phiNew; 
				phiLmin = phiNewL; 
			} else {
				phiMax = phiNew;
				phiLmax = phiNewL; 
			}
	} 
		else {
			if((phiNew-phiP) < (phiPhoton-phiP)) {
				phiMin = phiNew; 
				phiLmin = phiNewL; 
			} else {
				phiMax = phiNew;
				phiLmax = phiNewL; 
			}
		}
	} // end splitPhi()


};  // end class 

/*

// trkPC is track impact point at PC;
// TODO: later, better to use MIP pos?

if(getPionStatus()) {



  // check if rMax at phiL > R of photon
  std::vector<double> phiLvector, phiVector;
  
  auto phiL = switchLogic(thetaP, theta, phi, getMaxCkovPion(), phiLvector, phiVector);
  TVector2 posEstimatedPhi = tracePhoton(getMaxCkovPion(), phiL, lMin); // 

  auto r = (posEstimatedPhi - trkPC).Mod(); // find R of rMax (dist mip 2 {x,y} of estimated phi);
  
  if(r < rPhoton) { // ok, is in range
  }
} 

if(getProtonStatus()) {
  // check if rMax at phiL > R of photon
  std::vector<double> phiLvector, phiVector;
  
  auto phiL = switchLogic(thetaP, theta, phi, getMinCkovProton(), phiLvector, phiVector);
  TVector2 posEstimatedPhi = tracePhoton(getMinCkovProton(), phiL, lMin); // 

  auto r = (posEstimatedPhi - trkPC).Mod(); // find R of rMax (dist mip 2 {x,y} of estimated phi);
  
  if(r > rPhoton) { // ok, is in range
  }
} 
*/
