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


class Populate {


   using Polar3D = ROOT::Math::Polar3D<double>;
   using Rotation3D = ROOT::Math::Rotation3D;
   using RotationY = ROOT::Math::RotationY;
   using RotationZ = ROOT::Math::RotationZ;


private:

   
    TVector2 fTrkPos; // track pos in LORS at RAD // xRad, yRad
    TVector3 fTrkDir; // track position in LORS at RAD // setMagThetaPhi(1, thetaP, phiP)
    TVector2 fPc; // track pos at PC

    double nF;	     // refIdnex of freon

    const double winThick = 0.5, radThick = 1.5; const int gapThick = 8;
    const double getRefIdx = nF,  gapIdx = 1.0005, winIdx = 1.5787;

    double cosPhiRa, sinPhiRa;
    double tanThetaRa;
    double deltaX, deltaY;

public:

    Populate(TVector2 trkPos, TVector3 trkDir, double _nF) : fTrkPos(trkPos),  fTrkDir(trkDir), nF(_nF) 
    {

       
      //fPc.setX()
      Printf("init Populate class");

      Printf("Track pos at RAD : x %.3f y %.3f ", trkPos.X(), trkPos.Y());

      

    tanThetaRa = TMath::Tan(trkDir.Theta());
			cosPhiRa = TMath::Cos(trkDir.Phi());		
			sinPhiRa = TMath::Sin(trkDir.Phi());	
		

			// xRa = xPC - deltaX
      deltaX = (radThick + winThick + gapThick) * tanThetaRa * cosPhiRa;
      deltaY = (radThick + winThick + gapThick) * tanThetaRa * sinPhiRa;		
			
			setPcImp(fTrkPos.X() + deltaX, fTrkPos.Y() + deltaY);
      Printf("Track pos at PC : x %.3f y %.3f ", fPc.X(), fPc.Y());

      Printf("Track dir at RAD : theta %.3f phi %.3f ", trkDir.Theta(), trkDir.Phi());
				
    }


    // L here is "simulated" within 0..1.5 range
    ROOT::Math::Z tracePhot(const double& ckovThe, const double& ckovPhi, const double & L) const {
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
	RotationY rotY(fTrkDir.Theta());

        TRotation mphi;
        mphi.RotateZ(fTrkDir.Phi());
	RotationZ rotZ(fTrkDir.Phi());

        TRotation mrot = mphi * mtheta;
        Rotation3D mrot2 = rotZ * rotY;

        TVector3 dirCkovLORS;
        dirCkovLORS = mrot * dirCkov;
	
	Polar3D dirCkov2(dirCkov.Mag(), dirCkov.Theta(), dirCkov.Phi());	

        Polar3D dirCkovLORS2;
        dirCkovLORS2 = mrot2 * dirCkov2;

        phiCer = dirCkovLORS.Phi();     // actual value of the phi of the photon
        thetaCer = dirCkovLORS.Theta(); // actual value of thetaCerenkov of the photon


	
	Printf("trs2Lors");

	Printf("	old : rotZ %.3f, rotY %.3f", mrot.);
	Printf("	new : rotZ %.3f, rotY %.3f", mrot2);

	Printf("	old : phi %.3f, theta %.3f", phiCer, thetaCer);
	Printf("	new : phi %.3f, theta %.3f", dirCkovLORS2.Phi(), dirCkovLORS2.Theta());
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


		/*
		// track pos in rad 2 pc
    void rad2pc(const TVector& rad)
    {
			
    }

		// track pos in pc 2 rad
    void pc2rad(const TVector& pc) */ 
    
};
