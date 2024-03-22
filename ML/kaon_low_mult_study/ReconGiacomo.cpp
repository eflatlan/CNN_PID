#include <iostream>
#include <cmath>
#include <random>

#include <TNamed.h> //base class


// ef: new includes to replace TVector2/3


#include <vector>
// ef: using vectors instead of TClonesArray

// class AliESDtrack;  //CkovAngle() ef: commented out
// ef: what is eq in O2?

#include "HMPIDBase/Param.h"
#include "DataFormatsHMP/Cluster.h"
#include "ReconstructionDataFormats/Track.h"

#include "HMPIDBase/Param.h"
#include "HMPIDReconstruction/Recon.h" //class header
// #include "ReconstructionDataFormats/MatchInfoHMP.h"

#include <TRotation.h> //TracePhot()
#include <TH1D.h>      //HoughResponse()
#include <TRandom.h>   //HoughResponse()

#include "ReconstructionDataFormats/MatchInfoHMP.h"
#include "ReconstructionDataFormats/Track.h"
#include <TNamed.h> //base class


#include <TVector2.h>
#include <TVector3.h>

#include <vector>

#include "HMPIDBase/Param.h"
#include "DataFormatsHMP/Cluster.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/MatchInfoHMP.h"
using MatchInfo = o2::dataformats::MatchInfoHMP;

using namespace o2::hmpid;
// ClassImp(Recon);
ClassImp(o2::hmpid::Recon);

class ReconG 
{

 protected:
  int fPhotCnt; // counter of photons candidate

  //  ef : changed to smart-pointer arrays
  std::unique_ptr<int[]> fPhotFlag;      // flags of photon candidates
  std::unique_ptr<int[]> fPhotClusIndex; // cluster index of photon candidates
  std::unique_ptr<double[]> fPhotCkov;   // Ckov angles of photon candidates, [rad]
  std::unique_ptr<double[]> fPhotPhi;    // phis of photons candidates, [rad]
  std::unique_ptr<double[]> fPhotWei;    // weigths of photon candidates

  // int    *fPhotClusIndex;                     // cluster index of photon candidates

  double fCkovSigma2; // sigma2 of the reconstructed ring

  bool fIsWEIGHT;     // flag to consider weight procedure
  float fDTheta;      // Step for sliding window
  float fWindowWidth; // Hough width of sliding window

  double fRingArea; // area of a given ring
  double fRingAcc;  // fraction of the ring accepted by geometry

  TVector3 fTrkDir; // track direction in LORS at RAD

  TVector2 fTrkPos; // track positon in LORS at RAD   // XY mag
  TVector2 fMipPos; // mip positon for a given trackf // XY
  TVector2 fPc;     // track position at PC           // XY

  o2::hmpid::Param* fParam = o2::hmpid::Param::instance(); // Pointer to HMPIDParam

 private:
  ReconG(const ReconG& r);            // dummy copy constructor
  ReconG& operator=(const ReconG& r); // dummy assignment operator
  //
  ClassDef(ReconG, 3);

public:

  ReconG(double ) : fParam(o2::hmpid::Param::instance())

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  bool findPhotCkov(double cluX, double cluY, double& thetaCer, double& phiCer)
  {
    // Finds Cerenkov angle  for this photon candidate
    // Arguments: cluX,cluY - position of cadidate's cluster
    // Returns: Cerenkov angle

    TVector3 dirCkov;

    double zRad = -0.5 * fParam->radThick() - 0.5 * fParam->winThick();   // z position of middle of RAD
    TVector3 rad(fTrkPos.X(), fTrkPos.Y(), zRad);                         // impact point at middle of RAD
    TVector3 pc(cluX, cluY, 0.5 * fParam->winThick() + fParam->gapIdx()); // mip at PC
    double cluR = TMath::Sqrt((cluX - fPc.X()) * (cluX - fPc.X()) +
                                (cluY - fPc.Y()) * (cluY - fPc.Y())); // ref. distance impact RAD-CLUSTER
    double phi = (pc - rad).Phi();                                  // phi of photon

    double ckov1 = 0;
    double ckov2 = 0.75 + fTrkDir.Theta(); // start to find theta cerenkov in DRS
    const double kTol = 0.01;
    Int_t iIterCnt = 0;
    while (1) {
        if (iIterCnt >= 50) {
        return kFALSE;
        }
        double ckov = 0.5 * (ckov1 + ckov2);
        dirCkov.SetMagThetaPhi(1, ckov, phi);
        TVector2 posC = traceForward(dirCkov);   // trace photon with actual angles
        double dist = cluR - (posC - fPc).Mod(); // get distance between trial point and cluster position
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
        lors2Trs(dirCkov, thetaCer, phiCer);  // find ckov (in TRS:the effective Cherenkov angle!)
        return kTRUE;
        }
    }

    } // FindPhotTheta()
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    bool findPhotCkov2(double cluX, double cluY, double& thetaCer, double& phiCer)
    {

    // TVector3  emissionV;
    /** set emission point **/
    // emissionV.SetXYZ(xEm, yEm, zEm);

    // TVector3 directionV;
    /** set track direction vector **/
    // directionV.SetXYZ(trkPx, trkPy, trkPz);

    double zRad = -0.5 * fParam->radThick() - 0.5 * fParam->winThick(); // z position of middle of RAD

    
    TVector3 emissionV(fTrkPos.X(), fTrkPos.Y(), zRad);                 // impact point at middle of RAD

    TVector3 photonHitV, apparentV, surfaceV;
    photonHitV.SetXYZ(cluX, cluY, 0.5 * fParam->winThick() + fParam->gapIdx());
    apparentV = photonHitV - emissionV;
    surfaceV = emissionV;
    // surfaceV.SetZ(0);

    Double_t n1 = fParam->getRefIdx();
    Double_t n2 = 1.;
    Double_t apparentTheta = apparentV.Theta();
    Double_t correctedTheta = asin(n2 / n1 * sin(apparentTheta));
    Double_t deltaTheta = apparentTheta - correctedTheta;

    TVector3 perpV = apparentV.Cross(surfaceV);
    TVector3 cherenkovV = apparentV;
    // cherenkovV.Rotate(deltaTheta, perpV);

    lors2Trs(cherenkovV, thetaCer, phiCer);

    // thetaCer = cherenkovV.Angle(fTrkDir);
    // phiCer   = cherenkovV.Phi();

    return kTRUE;
    }
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    TVector2 traceForward(TVector3 dirCkov) const
    {
    // Trace forward a photon from (x,y) up to PC
    //  Arguments: dirCkov photon vector in LORS
    //    Returns: pos of traced photon at PC

    TVector2 pos(-999, -999);
    double thetaCer = dirCkov.Theta();
    if (thetaCer > TMath::ASin(1. / fParam->getRefIdx())) {
        return pos;
    }                                                                           // total refraction on WIN-GAP boundary
    double zRad = -0.5 * fParam->radThick() - 0.5 * fParam->winThick();         // z position of middle of RAD
    TVector3 posCkov(fTrkPos.X(), fTrkPos.Y(), zRad);                           // RAD: photon position is track position @ middle of RAD
    propagate(dirCkov, posCkov, -0.5 * fParam->winThick());                     // go to RAD-WIN boundary
    refract(dirCkov, fParam->getRefIdx(), fParam->winIdx());                    // RAD-WIN refraction
    propagate(dirCkov, posCkov, 0.5 * fParam->winThick());                      // go to WIN-GAP boundary
    refract(dirCkov, fParam->winIdx(), fParam->gapIdx());                       // WIN-GAP refraction
    propagate(dirCkov, posCkov, 0.5 * fParam->winThick() + fParam->gapThick()); // go to PC
    pos.Set(posCkov.X(), posCkov.Y());
    return pos;

    } // TraceForward()

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    void lors2Trs(TVector3 dirCkov, double& thetaCer, double& phiCer) const
    {
    // Theta Cerenkov reconstruction
    //  Arguments: dirCkov photon vector in LORS
    //    Returns: thetaCer of photon in TRS
    //               phiCer of photon in TRS
    //  TVector3 dirTrk;
    //  dirTrk.SetMagThetaPhi(1,fTrkDir.Theta(),fTrkDir.Phi()); -> dirTrk.SetCoordinates(1,fTrkDir.Theta(),fTrkDir.Phi())
    //  double  thetaCer = TMath::ACos(dirCkov*dirTrk);

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
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    void trs2Lors(TVector3 dirCkov, double& thetaCer, double& phiCer) const
    {
    // Theta Cerenkov reconstruction
    //  Arguments: dirCkov photon vector in TRS
    //    Returns: thetaCer of photon in LORS
    //               phiCer of photon in LORS

    // TRotation mtheta;
    // mtheta.RotateY(fTrkDir.Theta()); ef : changed to :

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

    void findRingGeom(double ckovAng, int level)
    {
    // Find area covered in the PC acceptance
    // Arguments: ckovAng - cerenkov angle
    //            level   - precision in finding area and portion of ring accepted (multiple of 50)
    //   Returns: area of the ring in cm^2 for given theta ckov

    Int_t kN = 50 * level;
    Int_t nPoints = 0;
    Double_t area = 0;

    Bool_t first = kFALSE;
    TVector2 pos1;

    for (Int_t i = 0; i < kN; i++) {
        if (!first) {
        pos1 = tracePhot(ckovAng, Double_t(TMath::TwoPi() * (i + 1) / kN)); // find a good trace for the first photon
        if (pos1.X() == -999) {
            continue;
        } // no area: open ring
        if (!fParam->isInside(pos1.X(), pos1.Y(), 0)) {
            pos1 = intWithEdge(fMipPos, pos1); // find the very first intersection...
        } else {
            if (!fParam->isInDead(pos1.X(), pos1.Y())) {
            nPoints++;
            } // photon is accepted if not in dead zone
        }
        first = kTRUE;
        continue;
        }
        TVector2 pos2 = tracePhot(ckovAng, Double_t(TMath::TwoPi() * (i + 1) / kN)); // trace the next photon
        if (pos2.X() == -999) {
        {
            continue;
        }
        } // no area: open ring
        if (!fParam->isInside(pos2.X(), pos2.Y(), 0)) {
        pos2 = intWithEdge(fMipPos, pos2);
        } else {
        if (!fParam->isInDead(pos2.X(), pos2.Y())) {
            nPoints++;
        } // photon is accepted if not in dead zone
        }
        area += TMath::Abs((pos1 - fMipPos).X() * (pos2 - fMipPos).Y() - (pos1 - fMipPos).Y() * (pos2 - fMipPos).X()); // add area of the triangle...
        pos1 = pos2;
    }
    //---  find area and length of the ring;
    fRingAcc = (Double_t)nPoints / (Double_t)kN;
    area *= 0.5;
    fRingArea = area;

    } // FindRingGeom()
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const TVector2 intWithEdge(TVector2 p1, TVector2 p2)
    {
    // It finds the intersection of the line for 2 points traced as photons
    // and the edge of a given PC
    // Arguments: 2 points obtained tracing the photons
    //   Returns: intersection point with detector (PC) edges

    double xmin = (p1.X() < p2.X()) ? p1.X() : p2.X();
    double xmax = (p1.X() < p2.X()) ? p2.X() : p1.X();
    double ymin = (p1.Y() < p2.Y()) ? p1.Y() : p2.Y();
    double ymax = (p1.Y() < p2.Y()) ? p2.Y() : p1.Y();

    double m = TMath::Tan((p2 - p1).Phi());
    TVector2 pint;
    // intersection with low  X
    pint.Set((double)(p1.X() + (0 - p1.Y()) / m), 0.);
    if (pint.X() >= 0 && pint.X() <= fParam->sizeAllX() &&
        pint.X() >= xmin && pint.X() <= xmax &&
        pint.Y() >= ymin && pint.Y() <= ymax) {
        return pint;
    }
    // intersection with high X
    pint.Set((double)(p1.X() + (fParam->sizeAllY() - p1.Y()) / m), (double)(fParam->sizeAllY()));
    if (pint.X() >= 0 && pint.X() <= fParam->sizeAllX() &&
        pint.X() >= xmin && pint.X() <= xmax &&
        pint.Y() >= ymin && pint.Y() <= ymax) {
        return pint;
    }
    // intersection with left Y
    pint.Set(0., (double)(p1.Y() + m * (0 - p1.X())));
    if (pint.Y() >= 0 && pint.Y() <= fParam->sizeAllY() &&
        pint.Y() >= ymin && pint.Y() <= ymax &&
        pint.X() >= xmin && pint.X() <= xmax) {
        return pint;
    }
    // intersection with righ Y
    pint.Set((double)(fParam->sizeAllX()), (double)(p1.Y() + m * (fParam->sizeAllX() - p1.X()))); // ef: Set->SetCoordinates
    if (pint.Y() >= 0 && pint.Y() <= fParam->sizeAllY() &&
        pint.Y() >= ymin && pint.Y() <= ymax &&
        pint.X() >= xmin && pint.X() <= xmax) {
        return pint;
    }
    return p1;
    } // IntWithEdge()


      TVector2 getMip() const
  {
    return fMipPos;
  } // mip coordinates

  double getRingArea() const
  {
    return fRingArea;
  } // area of the current ring in cm^2

  double getRingAcc() const
  {
    return fRingAcc;
  }                                                                                          // portion of the ring ([0,1]) accepted by geometry.To scale n. of photons
  double findRingExt(double ckov, int ch, double xPc, double yPc, double thRa, double phRa); // find ring acceptance by external parameters

  void setTrack(double xRad, double yRad, double theta, double phi)
  {
    fTrkDir.SetMagThetaPhi(1., theta, phi);
    fTrkPos.Set(xRad, yRad);
  } // set track parameter at RAD

  void setImpPC(double xPc, double yPc)
  {
    fPc.Set(xPc, yPc);
  } // set track impact to PC

  void setMip(double xmip, double ymip)
  {
    fMipPos.Set(xmip, ymip);
  } // set track impact to PC

};
