// populate.cpp

#include <TVector2.h>
#include <TVector3.h>
#include <TRotation.h>
#include <TRandom.h>
#include <vector>

class Recon {
private:
    /*TVector3 fTrkPos; // not defined in the provided code, declaration added for compilation
    TVector3 fTrkDir; // not defined in the provided code, declaration added for compilation
    Param *fParam;    // not defined in the provided code, declaration added for compilation
    */ 
   
public:
    TVector2 tracePhot(double ckovThe, double ckovPhi) const {
        double theta, phi;
        TVector3 dirTRS, dirLORS;
        dirTRS.SetMagThetaPhi(1, ckovThe, ckovPhi); // photon in TRS
        trs2Lors(dirTRS, theta, phi);
        dirLORS.SetMagThetaPhi(1, theta, phi); // photon in LORS
        return traceForward(dirLORS);          // now foward tracing
    }

    void propagate(const TVector3 dir, TVector3& pos, double z) const {
        static TVector3 nrm(0, 0, 1);
        TVector3 pnt(0, 0, z);

        TVector3 diff = pnt - pos;
        double sint = (nrm * diff) / (nrm * dir);
        pos += sint * dir;
    }

    void refract(TVector3& dir, double n1, double n2) const {
        double sinref = (n1 / n2) * TMath::Sin(dir.Theta());
        if (TMath::Abs(sinref) > 1.) {
            dir.SetXYZ(-999, -999, -999);
        } else {
            dir.SetTheta(TMath::ASin(sinref));
        }
    }

    TVector2 traceForward(TVector3 dirCkov) const {
        TVector2 pos(-999, -999);
        double thetaCer = dirCkov.Theta();
        if (thetaCer > TMath::ASin(1. / fParam->getRefIdx())) {
            return pos;
        }
        double zRad = -0.5 * fParam->radThick() - 0.5 * fParam->winThick();
        TVector3 posCkov(fTrkPos.X(), fTrkPos.Y(), zRad);
        propagate(dirCkov, posCkov, -0.5 * fParam->winThick());
        refract(dirCkov, fParam->getRefIdx(), fParam->winIdx());
        propagate(dirCkov, posCkov, 0.5 * fParam->winThick());
        refract(dirCkov, fParam->winIdx(), fParam->gapIdx());
        propagate(dirCkov, posCkov, 0.5 * fParam->winThick() + fParam->gapThick());
        pos.Set(posCkov.X(), posCkov.Y());
        return pos;
    }

    void lors2Trs(TVector3 dirCkov, double& thetaCer, double& phiCer) const {
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

    void trs2Lors(TVector3 dirCkov, double& thetaCer, double& phiCer) const {
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
