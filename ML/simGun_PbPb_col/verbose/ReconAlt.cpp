#ifndef Recon2_H

#define Recon2_H

// import ROOT

#include "SimpleCluster.cpp"

// #include <ROOT/XYZVector.h>

#include <Math/GenVector/Rotation3D.h>

#include "Math/GenVector/RotationX.h"

#include "Math/GenVector/RotationY.h"

#include "Math/GenVector/RotationZ.h"

#include "Math/Vector2D.h"
#include "Math/Vector3D.h"

#include "Math/Vector4D.h"
// #include "math_utils.h" // Include the necessary header file for math_utils

// Add the necessary using directive for the 'o2' namespace
using namespace o2;
using namespace ROOT::Math;

class Recon2 {

private:
  Recon2(const Recon2 &r);            // dummy copy constructor
  Recon2 &operator=(const Recon2 &r); // dummy assignment operator

public:
  // Additional methods and private members...
  Recon2(double thetaP, double phiP, double xMip, double yMip, double xRad,
         double yRad, const std::vector<SimpleCluster> &clustersIn) {

    fTrkDir.SetCoordinates(1., thetaP, phiP);

    fTrkPos.SetXY(xRad, yRad);

    fMipPos.SetXY(xMip, yMip);

    fPc.SetXY(xMip, yMip);

    clusters = clustersIn;

    mThetaP = thetaP;

    mPhiP = phiP;
  }

  /* TVector3 fTrkDir; // track direction in LORS at RAD

TVector2 fTrkPos; // track positon in LORS at RAD   // XY mag
TVector2 fMipPos; // mip positon for a given trackf // XY
TVector2 fPc;     // track position at PC           // XY*/

  //  math_utils::::Vector3D<double> posCkov(fTrkPos.X(), fTrkPos.Y(), zRad);

  Polar3DVector fTrkDir; // track direction in LORS at RAD

  XYVector fTrkPos; // track positon in LORS at RAD   // XY mag
  XYVector fMipPos; // mip positon for a given trackf // XY
  XYVector fPc;     // track position at PC           // XY

  double mThetaP, mPhiP;

  std::vector<SimpleCluster> clusters;

  double radThick() const { return 1.5; }

  double winThick() const { return 0.5; }

  double gapThick() const { return 8.0; }

  double winIdx() const { return 1.583; }

  double gapIdx() const { return 1.0005; }

  double getRefIdx() const { return 1.2904; }

  // template <typename T>

  bool findPhotCkov(double cluX, double cluY, double &thetaCer, double &phiCer)

  {

    // Finds Cerenkov angle  for this photon candidate

    // Arguments: cluX,cluY - position of cadidate's cluster

    // Returns: Cerenkov angle

    Polar3DVector dirCkov;

    double zRad =
        -0.5 * radThick() - 0.5 * winThick(); // z position of middle of RAD

    XYZVector rad(fTrkPos.X(), fTrkPos.Y(),
                  zRad); // impact point at middle of RAD

    XYZVector pc(cluX, cluY, 0.5 * winThick() + gapIdx()); // mip at PC

    double cluR = TMath::Sqrt(
        (cluX - fTrkPos.X()) * (cluX - fTrkPos.X()) +

        (cluY - fTrkPos.Y()) *
            (cluY - fTrkPos.Y())); // ref. distance impact RAD-CLUSTER

    double phi = (pc - rad).Phi(); // phi of photon

    double ckov1 = 0;

    double ckov2 =
        0.75 + fTrkDir.Theta(); // start to find theta cerenkov in DRS

    const double kTol = 0.01;

    int iIterCnt = 0;

    while (1) {

      if (iIterCnt >= 50) {

        return kFALSE;
      }

      double ckov = 0.5 * (ckov1 + ckov2);

      dirCkov.SetCoordinates(1, ckov, phi);

      XYVector posC = traceForward(dirCkov); // trace photon with actual angles

      double dist =
          cluR -
          (posC - fTrkPos)
              .Mag2(); // get distance between trial point and cluster position

      if (posC.X() == -999) {

        dist = -999; // total reflection problem
      }

      iIterCnt++; // counter step

      if (dist > kTol) {

        ckov1 = ckov; // cluster @ larger ckov

      } else if (dist < -kTol) {

        ckov2 = ckov; // cluster @ smaller ckov

      } else { // precision achived: ckov in DRS found

        dirCkov.SetCoordinates(1, ckov, phi);

        lors2Trs(dirCkov, thetaCer,
                 phiCer); // find ckov (in TRS:the effective Cherenkov angle!)

        return kTRUE;
      }
    }

  } // FindPhotTheta()

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // template <typename T> // typename

  XYVector traceForward(Polar3DVector dirCkov) const

  {

    // Trace forward a photon from (x,y) up to PC

    //  Arguments: dirCkov photon vector in LORS

    //    Returns: pos of traced photon at PC

    XYVector pos(-999, -999);

    double thetaCer = dirCkov.Theta();

    if (thetaCer > TMath::ASin(1. / getRefIdx())) {

      return pos; // total refraction on WIN-GAP boundary
    }

    double zRad =
        -0.5 * radThick() - 0.5 * winThick(); // z position of middle of RAD

    XYZVector posCkov(fTrkPos.X(), fTrkPos.Y(), zRad);

    // RAD: photon position is track position @ middle of RAD

    propagate(dirCkov, posCkov, -0.5 * winThick()); // go to RAD-WIN boundary

    refract(dirCkov, getRefIdx(), winIdx()); // RAD-WIN refraction

    propagate(dirCkov, posCkov, 0.5 * winThick()); // go to WIN-GAP boundary

    refract(dirCkov, winIdx(), gapIdx()); // WIN-GAP refraction

    propagate(dirCkov, posCkov, 0.5 * winThick() + gapThick()); // go to PC

    pos.SetCoordinates(posCkov.X(), posCkov.Y());

    return pos;

  } // TraceForward()

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  void lors2Trs(Polar3DVector dirCkov, double &thetaCer, double &phiCer) const

  {

    // Theta Cerenkov reconstruction

    //  Arguments: dirCkov photon vector in LORS

    //    Returns: thetaCer of photon in TRS

    //               phiCer of photon in TRS

    //  TVector3 dirTrk;

    //  dirTrk.SetMagThetaPhi(1,fTrkDir.Theta(),fTrkDir.Phi()); ->
    //  dirTrk.SetCoordinates(1,fTrkDir.Theta(),fTrkDir.Phi())

    //  double  thetaCer = TMath::ACos(dirCkov*dirTrk);

    ROOT::Math::Rotation3D mtheta(ROOT::Math::RotationY(
        -fTrkDir
             .Theta())); // TRotation mtheta;  mtheta.RotateY(-fTrkDir.Theta());

    ROOT::Math::Rotation3D mphi(
        ROOT::Math::RotationZ(-fTrkDir.Phi())); // mphi.RotateZ(-fTrkDir.Phi());

    ROOT::Math::Rotation3D mrot = mtheta * mphi;

    Polar3DVector dirCkovTRS;

    dirCkovTRS = mrot * dirCkov;

    phiCer = dirCkovTRS.Phi(); // actual value of the phi of the photon

    thetaCer =
        dirCkovTRS.Theta(); // actual value of thetaCerenkov of the photon
  }

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  void trs2Lors(Polar3DVector dirCkov, double &thetaCer, double &phiCer) const

  {

    // Theta Cerenkov reconstruction

    //  Arguments: dirCkov photon vector in TRS

    //    Returns: thetaCer of photon in LORS

    //               phiCer of photon in LORS

    ROOT::Math::Rotation3D mtheta(ROOT::Math::RotationY(
        fTrkDir.Theta())); // TRotation mtheta; mtheta.RotateY(fTrkDir.Theta());

    ROOT::Math::Rotation3D mphi(ROOT::Math::RotationZ(fTrkDir.Phi()));

    ROOT::Math::Rotation3D mrot = mphi * mtheta;

    Polar3DVector dirCkovLORS = mrot * dirCkov;

    phiCer = dirCkovLORS.Phi(); // actual value of the phi of the photon

    thetaCer =
        dirCkovLORS.Theta(); // actual value of thetaCerenkov of the photon
  }

  void propagate(const Polar3DVector &dir, XYZVector &pos, double z) const

  {

    // Finds an intersection point between a line and XY plane shifted along Z.

    // Arguments:  dir,pos   - vector along the line and any point of the line

    //             z         - z coordinate of plain

    //   Returns:  none

    //   On exit:  pos is the position if this intesection if any

    XYZVector nrm(0, 0, 1);

    XYZVector pnt(0, 0, z);

    XYZVector diff = pnt - pos;

    double sint = 0; //(nrm * diff) / (nrm * dir);

    pos += sint * dir;

  } // Propagate()

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // template <typename T> // typename

  void refract(Polar3DVector &dir, double n1, double n2) const

  {

    // Refract direction vector according to Snell law

    // Arguments:

    //            n1 - ref idx of first substance

    //            n2 - ref idx of second substance

    //   Returns: none

    //   On exit: dir is new direction

    double sinref = (n1 / n2) * TMath::Sin(dir.Theta());

    if (TMath::Abs(sinref) > 1.) {

      dir.SetXYZ(-999, -999, -999);

    } else {

      dir.SetTheta(TMath::ASin(sinref));
    }
  }
};

#endif
