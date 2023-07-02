
#ifndef TEST_POPULATE2
#define TEST_POPULATE2
// populate.cpp
#include <TVector2.h>
#include <TVector3.h>
#include <TRotation.h>
#include <TRandom.h>
#include <vector>


//#include "Math/GenVector/Rotation3D.h"
#include "Math/GenVector/Rotation3D.h"
#include "Math/GenVector/Cartesian3D.h"
#include "Math/GenVector/DisplacementVector3D.h"

#include "Math/GenVector/RotationYfwd.h"
#include "Math/GenVector/RotationZfwd.h"


class Populate2 {

	 // using array = std::array;
	 using vecArray3 = std::vector<std::array<double,3>>;
   using Polar3D = ROOT::Math::Polar3D<double>;
   using Rotation3D = ROOT::Math::Rotation3D;
   using RotationY = ROOT::Math::RotationY;
   using RotationZ = ROOT::Math::RotationZ;


private:

   
    TVector2 fTrkPos; // track pos in LORS at RAD // xRad, yRad
    TVector3 fTrkDir; // track position in LORS at RAD // setMagThetaPhi(1, thetaP, phiP)
    TVector2 fPc; // track pos at PC

    TVector2 fMipPos; // MIP-pos at PC

    double nF;	     // refIdnex of freon

    const double winThick = 0.5, radThick = 1.5; const int gapThick = 8;
    const double getRefIdx = nF,  gapIdx = 1.0005, winIdx = 1.5787;

    const double nQ = winIdx, nG = gapIdx;

    double thetaP, cosPhiP, sinPhiP, tanThetaP;
    double deltaX, deltaY;

    double phiP, sinThetaP, cosThetaP;


    TVector2 fTrkPos2D;
public:

    Populate2(TVector2 trkPos, TVector3 trkDir, double _nF) : fTrkPos(trkPos),  fTrkDir(trkDir), nF(_nF) 
    {


      fTrkPos2D.Set(trkPos.X(), trkPos.Y()); 


      //fPc.setX()
      Printf("init Populate class");

      Printf("Track pos at RAD : x %.3f y %.3f ", trkPos.X(), trkPos.Y());

      
    thetaP = trkDir.Theta();
    tanThetaP = TMath::Tan(thetaP);
    cosThetaP = TMath::Cos(thetaP);		
    sinThetaP = TMath::Sin(thetaP);



    phiP = trkDir.Phi();
		cosPhiP = TMath::Cos(phiP);		
		sinPhiP = TMath::Sin(phiP);	
		




			// xRa = xPC - deltaX



      // should be radThick/2 here if assuming half em-length 
      deltaX = (radThick/2 + winThick + gapThick) * tanThetaP * cosPhiP;
      deltaY = (radThick/2 + winThick + gapThick) * tanThetaP * sinPhiP;		
			

      // NB! TODO: here PC impact point based on L = rW/2!!
      setPcImp(fTrkPos.X() + deltaX, fTrkPos.Y() + deltaY);
      Printf("Track pos at PC : x %.3f y %.3f ", fPc.X(), fPc.Y());

      Printf("Track dir at RAD : theta %.3f phi %.3f ", trkDir.Theta(), trkDir.Phi());
				
    }

    // check if   rMin < r_photon  for a given photon {x,y} -->Phi in LORS
    bool checkRangeAbove(const TVector2& photonPos, const double& etaMin,  TVector2& rPos)
    {

      const auto rPhoton = (photonPos - fPc).Mod();
      // TODO: better to use MIP-pos than fPC? 
      auto lMax = 1.5; 
      Printf("checkRangeAbove : etaMin %.3f ---> enter getRatPhi() ", etaMin);
      auto rMin = getRatPhi(photonPos,  etaMin, lMax, rPos);
      Printf("CheckRangeAbove : rMin %.3f < rPhoton %.3f \n", rMin, rPhoton);
      return (rMin < rPhoton);
    }

    // check if   r_photon  < rMax  for a given photon {x,y} -->Phi in LORS
    bool checkRangeBelow(const TVector2& photonPos, const double& etaMax, TVector2& rPos)
    {
      const auto rPhoton = (photonPos - fPc).Mod();
      // TODO: better to use MIP-pos than fPC? 
      auto lMin = 0; 
      Printf("checkRangeBelow : etaMax %.3f ---> enter getRatPhi()", etaMax);
      auto rMax = getRatPhi(photonPos, etaMax, lMin, rPos);
      Printf("checkRangeBelow : rPhoton %.3f, rMax %.3f \n", rPhoton, rMax);
      return (rPhoton < rMax);
    }

    // check if   rMin < r_photon  < rMax for a given photon {x,y} -->Phi in LORS
    bool checkRange2(const TVector2& photonPos, const double& etaMax, const double& etaMin, TVector2& rPos)
    {

      const auto rPhoton = (photonPos - fPc).Mod();
      // TODO: better to use MIP-pos than fPC? 
      auto lMin = 0.0, lMax = 1.5; 
      auto rMax = getRatPhi(photonPos, etaMax, lMin, rPos);
      auto rMin = getRatPhi(photonPos, etaMin, lMax, rPos);
      Printf("CheckRange : etaMax %.3f , etaMin %.3f", etaMax, etaMin);
      Printf("CheckRange : rMin %.3f < rPhoton %.3f, rMax %.3f", rMin, rPhoton, rMax);
      return (rMin < rPhoton && rPhoton < rMax);
    }


    // find R at specific value of Phi to apply mass-hyp
    double getRatPhi(const TVector2& photonPos, const double& eta, const double& L, TVector2& rPos)
    {
      

      // TODO: better to use MIP-pos than fPC? 

      // denne skal nok vaere saann:
      const auto phi = (photonPos - fTrkPos2D).Phi();


      const auto sinPhi = TMath::Sin(phi);
      const auto cosPhi = TMath::Cos(phi);

      TVector3 dirPhotonR;

      
      // denne skal nok ogsaa vaere saann:
      const auto cosTheta = (TMath::Cos(eta) - sinThetaP * TMath::Cos(phi-phiP))/cosThetaP;
      const auto theta = TMath::ACos(cosTheta);

      // must scale a, b here !(?:)) :

      const auto tanTheta = TMath::Tan(theta);
      const auto sinTheta = TMath::Sin(theta);

      const auto a = (radThick - L + winThick + gapThick)*tanThetaP;


      const auto nFSinSq = nF*nF*sinTheta*sinTheta;
      
      const auto b = (radThick - L)*tanTheta + winThick * nF * sinTheta/TMath::Sqrt(nQ*nQ - nFSinSq) + gapThick * nF * sinTheta/TMath::Sqrt(nG*nG - nFSinSq);

      const auto xDiff = (a*cosPhiP - b*cosPhi);
      const auto yDiff = (a*sinPhiP - b*sinPhi);
      const auto R = TMath::Sqrt(xDiff*xDiff + yDiff*yDiff);


     Printf("eta %.2f , Math::Cos(eta) %.2f | sinThetaP * TMath::Cos(phi-phiP) %.2f | cosThetaP %.2f", eta, TMath::Cos(eta), sinThetaP * TMath::Cos(phi-phiP), cosThetaP);


     Printf("cosTheta %.2f | theta %.2f | L %.2f",cosTheta, theta, L);


      // set max/min etaC value
      // dirPhotonR.SetMagThetaPhi(1, eta, phi);
      dirPhotonR.SetMagThetaPhi(1, theta, phi);

      // set max/min L value 	
      rPos = traceForward(dirPhotonR, L);  
      // for the given phi value in local-ref-system, L{min/max}, eta{min/max}, 
      // create the point for the mass-hyp
      
      Printf("getRatPhi : fPC: x %.2f y %.2f | Photon  x %.2f y %.2f | rPos x %.2f y %.2f", fPc.X(), fPc.Y(), photonPos.X(), photonPos.Y(), rPos.X(), rPos.Y());
     	
      // as for findphotckov : cluR = sqrt([cluX - fPc.X()]^2 [y..])
      auto dist = (rPos - fPc).Mod();    

      Printf("getRatPhi : dist = %.2f | R(a,b) = %.2f", dist, R);

      return R; // TODO : change back to return dist?	
    }


    // L here is "simulated" within 0..1.5 range
    TVector2 tracePhot(const double& ckovThe, const double& ckovPhi, const double & L) const {
        double theta, phi;
        TVector3 dirTRS, dirLORS;
        dirTRS.SetMagThetaPhi(1, ckovThe, ckovPhi); // photon in TRS
        trs2Lors(dirTRS, theta, phi);
        dirLORS.SetMagThetaPhi(1, theta, phi); // photon in LORS
        return traceForward(dirLORS, L);          // now foward tracing
    }

    void propagate(const TVector3& dir, TVector3& pos, double z) const {
        static TVector3 nrm(0, 0, 1);
        TVector3 pnt(0, 0, z);

        TVector3 diff = pnt - pos;
        double sint = (nrm * diff) / (nrm * dir);
        pos += sint * dir;
    }

    void refract(TVector3& dir, const double& n1, const double& n2) const {
        double sinref = (n1 / n2) * TMath::Sin(dir.Theta());
        if (TMath::Abs(sinref) > 1.) {
            dir.SetXYZ(-999, -999, -999);
        } else {
            dir.SetTheta(TMath::ASin(sinref));
        }
    }




    TVector2 traceForward(TVector3& dirCkov, const double& L) const {

	auto getRefIdx = nF;

        TVector2 pos(-999, -999);
        double thetaCer = dirCkov.Theta();
        if (thetaCer > TMath::ASin(1. / getRefIdx)) {
            return pos;
        }

	// change radThick to other value to change L
 	// auto radThick' = (radThick - L);
        // double zRad = - radThick' - 0.5 * winThick;
        double zRad = - (radThick - L) - 0.5 * winThick; 
  
	// TODO: which value should be changed??

        TVector3 posCkov(fTrkPos.X(), fTrkPos.Y(), zRad);
        propagate(dirCkov, posCkov, -0.5 * winThick); // TODO giacomo spm :er ikke dette også feil!
        refract(dirCkov, getRefIdx, winIdx);
        propagate(dirCkov, posCkov, 0.5 * winThick);
        refract(dirCkov, winIdx, gapIdx);
        propagate(dirCkov, posCkov, 0.5 * winThick + gapThick);
        pos.Set(posCkov.X(), posCkov.Y());
        return pos;
    }

    void lors2Trs(const TVector3& dirCkov, double& thetaCer, double& phiCer) const {
        TRotation mtheta;
        mtheta.RotateY(-fTrkDir.Theta());

        TRotation mphi;
        mphi.RotateZ(-fTrkDir.Phi());

        TRotation mrot = mtheta * mphi;

        TVector3 dirCkovTRS;
        dirCkovTRS = mrot * dirCkov;
        phiCer = dirCkovTRS.Phi();     // actual value of the phi of the photon
        thetaCer = dirCkovTRS.Theta(); // actual value of thetaCerenkov of the photon
    }

    void trs2Lors(const TVector3& dirCkov, double& thetaCer, double& phiCer) const {

        TRotation mtheta;
        mtheta.RotateY(fTrkDir.Theta());

        TRotation mphi;
        mphi.RotateZ(fTrkDir.Phi());

        TRotation mrot = mphi * mtheta;


        TVector3 dirCkovLORS;
        dirCkovLORS = mrot * dirCkov;
	

        //Polar3D dirCkovLORS2;
        //dirCkovLORS2 = mrot * dirCkov;

        phiCer = dirCkovLORS.Phi();     // actual value of the phi of the photon
        thetaCer = dirCkovLORS.Theta(); // actual value of thetaCerenkov of the photon

	/*Printf("trs2Lors");
	Printf("	old : phi %.3f, theta %.3f", phiCer, thetaCer);
	Printf("	new : phi %.3f, theta %.3f", dirCkovLORS2.Phi(), dirCkovLORS2.Theta());*/
    }


    // getter and setter functions

    void setTrackPos(double x, double y)
    {
			fTrkPos.SetX(x); 
 			fTrkPos.SetY(y); 
    }

    void setTrackPos(const TVector2& fTrkIn)
    {
			fTrkPos = fTrkIn;
    }
    
    TVector2 getTrackPos() const
    {
			return fTrkPos;
    }


    void setPcImp(double x, double y)
    {
			fPc.SetX(x); 
 			fPc.SetY(y); 
    }

    void setPcImp(const TVector2& fPcImpIn)
    {
			fPc = fPcImpIn;
    }
    
    TVector2 getPcImp() const
    {
			return fPc;
    }



	void checkCond(const TVector2& posPhoton, const double& rPhoton, const double& phiPhoton, bool getAbove, vecArray3& vec, const double& etaCkov) {
		// vec : contains phiL, phi, R w etaMin/etaMax for hadron species

		const double lMin = 0., lMax = 1.5;
    bool noDecisionTaken = false, condition = false;

    //const std::array<double, 3>	fst = vec[index];	

    const size_t kN = vec.size();
    int initValue;
		
		auto phiC = phiPhoton - phiP;

		if (phiC < TMath::Pi()/2)
			initValue = 0;
		else if (phiC < TMath::Pi())
			initValue = static_cast<int>(kN/4);
		else if (phiC < (3/2)*TMath::Pi())
			initValue = static_cast<int>(kN/2);
		else if (phiC < TMath::TwoPi())
			initValue = static_cast<int>((3/2)*kN);
		else {
			initValue = 0;
			Printf("phiC %.2f", phiC);
			throw std::invalid_argument("wtf value does phiC have?");
		} 

		int iCnt = 0;
		auto phi1 = vec[initValue][1]; // 


		// set increment opposite way if phi1 > phiPhoton
		int inc = 1;
		if (phi1 > phiPhoton)
			inc = -1;		

		// accesing the correct phi
		while(phi1 < phiPhoton) {
			phi1 = vec[initValue + iCnt][1];
			iCnt += inc;
		}

		// TODO: this has to be changed if inc = -1?
		int index = iCnt -1;


		//double phiL1 = vec.at(index).at(0);
		double phiL1 = vec[index][0];
		double phiL2 = vec[index+1][0];
		
		phi1 = vec[index][1];
		double phi2 = vec[index+1][1];

		// get radiuses
		double r1 = vec[index][2];
		double r2 = vec[index+1][2];

		// correct indexes are found 
		Printf("phiPhoton %.2f| phiPhoton -  1 %.2f, phiPhoton + 1 %.2f", phiPhoton, vec[initValue + iCnt][1], vec[initValue + iCnt+ 1][1]);
		
				

		// eta is min/max ckovHyps +- 3 std-dev
		
	  double rMin, rMax, phiMin, phiMax, phiLmin, phiLmax;
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
					rMax = getR(etaCkov, phiLmax, lMax);
					rMin = getR(etaCkov, phiLmin, lMax); // TODO: should it be lMin here?
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
					rMax = getR(etaCkov, phiLmax, lMin);
					rMin = getR(etaCkov, phiLmin, lMin);
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
		const TVector2& rPosLORS = tracePhot(etaTRS, phiTRS, L); // pos of estimated rPos at LORS PC
		const auto R = (rPosLORS-getPcImp()).Mod(); 
		return R;
	}
		
	void splitPhi(double& phiMax, double& phiLmax, double& phiMin, double& phiLmin, const double& etaCkov, const double& phiPhoton)
	{
		const auto& phiNewL = (phiLmax + phiLmin)/2.0;	// phi in TRS


		// Get phi in LORS 
		TVector3 dirNewTrs;
		dirNewTrs.SetMagThetaPhi(1, etaCkov, phiNewL); 
		
		
		trs2Lors(dirNewTrs, thetaP, phiP); 
		
		const auto& phiNew = dirNewTrs.Phi();

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


		// track pos in pc 2 rad
    // void pc2rad(const TVector& pc) */ 
    ClassDefNV(Populate2, 0);
};

#endif
