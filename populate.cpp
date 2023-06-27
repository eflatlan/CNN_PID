// populate.cpp

#include <TVector2.h>
#include <TVector3.h>
#include <TRotation.h>
#include <TRandom.h>
#include <vector>

class Populate {
private:


    TVector2 fTrkPos; // track pos in LORS at RAD // xRad, yRad
    TVector3 fTrkDir; // track position in LORS at RAD // setMagThetaPhi(1, thetaP, phiP)


    double nF;	     // refIdnex of freon

public:
    TVector2 fPc; // track pos at PC
    Populate(TVector2 trkPos, TVector3 trkDir, double _nF) : fTrkPos(trkPos),  fTrkDir(trkDir), nF(_nF) 
    {

       
      //fPc.setX()
      Printf("init Populate class");

      Printf("Track pos at RAD : x %.3f y %.3f ", trkPos.X(), trkPos.Y());
      Printf("Track dir at RAD : theta %.3f phi %.3f ", trkDir.Theta(), trkDir.Phi());
      //fPc(trkPos.X())
    }

    TVector2 tracePhot(const double& ckovThe, const double& ckovPhi) const {
        double theta, phi;
        TVector3 dirTRS, dirLORS;
        dirTRS.SetMagThetaPhi(1, ckovThe, ckovPhi); // photon in TRS
        trs2Lors(dirTRS, theta, phi);
        dirLORS.SetMagThetaPhi(1, theta, phi); // photon in LORS
        return traceForward(dirLORS);          // now foward tracing
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




    TVector2 traceForward(TVector3& dirCkov) const {


	auto winThick = 0.5, radThick = 1.5; int gapThick = 8;
	auto getRefIdx = nF,  gapIdx = 1.0005, winIdx = 1.5787;

        TVector2 pos(-999, -999);
        double thetaCer = dirCkov.Theta();
        if (thetaCer > TMath::ASin(1. / getRefIdx)) {
            return pos;
        }

        double zRad = -0.5 * radThick - 0.5 * winThick;
        TVector3 posCkov(fTrkPos.X(), fTrkPos.Y(), zRad);
        propagate(dirCkov, posCkov, -0.5 * winThick);
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

        phiCer = dirCkovLORS.Phi();     // actual value of the phi of the photon
        thetaCer = dirCkovLORS.Theta(); // actual value of thetaCerenkov of the photon
    }
};
