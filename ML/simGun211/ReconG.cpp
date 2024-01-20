#include <iostream>
#include <cmath>
#include <random>

#include <TNamed.h> //base class

//##include <Vector3D.h>#include <Vector2D.h>
#include <TVector2.h>
#include <TVector3.h>
#include <vector>
// ef: using vectors instead of TClonesArray

// class AliESDtrack;  //CkovAngle() ef: commented out
// ef: what is eq in O2?

#include "HMPIDBase/Param.h"
//#include "DataFormatsHMP/Cluster.h"
//#include "ReconstructionDataFormats/Track.h"
using namespace o2;
#include "MathUtils/Cartesian.h"

//namespace ParticleUtils

class ReconG {

 protected:
  int fPhotCnt; // counter of photons candidate

  //  ef: changed to smart-pointer arrays
  // int *fPhotFlag;
  std::unique_ptr<int[]> fPhotFlag;      // flags of photon candidates
  std::unique_ptr<int[]> fPhotClusIndex; // cluster index of photon candidates
  std::unique_ptr<double[]> fPhotCkov;   // Ckov angles of photon candidates, [rad]
  std::unique_ptr<double[]> fPhotPhi;    // phis of photons candidates, [rad]
  std::unique_ptr<double[]> fPhotWei;    // weigths of photon candidates

  double fCkovSigma2; // sigma2 of the reconstructed ring

  bool fIsWEIGHT;     // flag to consider weight procedure
  float fDTheta;      // Step for sliding window
  float fWindowWidth; // Hough width of sliding window

  double fRingArea; // area of a given ring
  double fRingAcc;  // fraction of the ring accepted by geometry


  /*
  math_utils::Vector3D<double> fTrkDir; // track direction in LORS at RAD

  math_utils::Vector2D<double> fTrkPos; // track positon in LORS at RAD   // XY mag
  math_utils::Vector2D<double> fMipPos; // mip positon for a given trackf // XY
  math_utils::Vector2D<double> fPc;     // track position at PC           // XY
  */
  TVector3 fTrkDir; // track direction in LORS at RAD

  TVector2 fTrkPos; // track positon in LORS at RAD   // XY mag
  TVector2 fMipPos; // mip positon for a given trackf // XY
  TVector2 fPc;     // track position at PC           // XY

  std::unique_ptr<o2::hmpid::Param> fParam; // Pointer to HMPIDParam
  double refIdx;
  
    const double radThick = 1.5;
    const double winThick = 0.5;
    const double gapThick = 8;
    const double winIdx = 1.5787;
    const double gapIdx = 1.0005;
    double etaC;
 private:
  ReconG(const ReconG& r);            // dummy copy constructor
  ReconG& operator=(const ReconG& r); // dummy assignment operator

 public: //   ReconG reconG(thetaP, phiP, xP, yP, xPC, yPC, nF);
    ReconG(double _theta, double _phi, double _xRad, double _yRad, double _xPC, double _yPC, double n, double _etaC) : refIdx(n), etaC(_etaC)
    { 
      setTrack(_xRad, _yRad, _theta, _phi);
      setImpPC(_xPC, _yPC);



     Printf("fTrkPos: x %f y %f ", fTrkPos.X(), fTrkPos.Y());
     Printf("fPC : x %f y %f ", fPc.X(), fPc.Y());
    }

    ReconG(double _theta, double _phi, double _xRad, double _yRad, double _xPC, double _yPC, double n) : refIdx(n)
    { 
      setTrack(_xRad, _yRad, _theta, _phi);
      setImpPC(_xPC, _yPC);



     Printf("fTrkPos: x %f y %f ", fTrkPos.X(), fTrkPos.Y());
     Printf("fPC : x %f y %f ", fPc.X(), fPc.Y());
    }


/*    ReconG(double _xRad, double _yRad, double _theta, double _phi, double _xPC, double _yPC) : / *fParam(o2::hmpid::Param::instance()) * /
    { 
      setTrack(_xRad, _yRad, _theta, _phi);
      setImpPC(_xPC, _yPC);
    }*/

bool findPhotCkov(double cluX, double cluY, double& thetaCer, double& phiCer, double etaC)
  {

    bool found = false;
    // Finds Cerenkov angle  for this photon candidate
    // Arguments: cluX,cluY - position of cadidate's cluster
    // Returns: Cerenkov angle

    TVector3 dirCkov;
        

    double zRad = -0.5 * radThick;   // z position of middle of RAD
    //double zRad = -0.5 * radThick - 0.5 * winThick;   // z position of middle of RAD
    TVector3 rad(fTrkPos.X(), fTrkPos.Y(), zRad);                         // impact point at middle of RAD


    // TODO : this was changed:
    TVector3 pc(cluX, cluY, 0.5 * radThick + winThick + gapThick); // mip at PC
    // TVector3 pc(cluX, cluY, 0.5 * winThick + gapIdx); // mip at PC
    double cluR = TMath::Sqrt((cluX - fPc.X()) * (cluX - fPc.X()) +
                                (cluY - fPc.Y()) * (cluY - fPc.Y())); // ref. distance impact RAD-CLUSTER
    double phi = (pc - rad).Phi();                                  // phi of photon
    Printf("ReconG findCkov: cluR  %.3f " , cluR);

    Printf("ReconG findCkov: cluX %.3f fPcX %.3f" ,cluX, fPc.X());
    Printf("ReconG findCkov: cluY %.3f fPcY %.3f" ,cluY, fPc.Y());

    double ckov1 = 0;
    double ckov2 = 0.75 + fTrkDir.Theta(); // start to find theta cerenkov in DRS
    const double kTol = 0.01;
    Int_t iIterCnt = 0;

    TGraph* tCkovReconGGraph = new TGraph(50);	
    
    while (1) {

        if (iIterCnt >= 50) {
          found = false;
	  break;
        }
        double ckov = 0.5 * (ckov1 + ckov2);
	
        dirCkov.SetMagThetaPhi(1, ckov, phi);
        TVector2 posC = traceForward(dirCkov);   // trace photon with actual angles
        double dist = cluR - (posC - fPc).Mod(); // get distance between trial point and cluster position


	tCkovReconGGraph->SetPoint(iIterCnt, iIterCnt, ckov);
        
	/*
	auto sinThetaP = TMath::Sin(fTrkDir.Theta());

        auto cosDiffPhi = TMath::Cos(phi-fTrkDir.Phi());

        auto cosThetaProd =  TMath::Cos(theta)*TMath::Cos(fTrkDir.Theta());

	auto etaC = TMath::Acos(sinThetaP + cosDiffPhi + cosThetaProd);
	Printf("etaC = %.4f",etaC);*/

	// TODO: nb!! cant really see where cosEtaC = sinThetaP *cos(phi-phiP).. is done??
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

          found = true;
        }



//	Printf("ReconG findCkov: cnt %d| ckov %.4f : ckov1 %.4f : ckov2 %.4f | dist %.3f" , iIterCnt, ckov, ckov1, ckov2, dist);


/*
	Printf("ReconG findCkov: cluX %.3f fPcX %.3f" ,cluX, fPc.X());
	Printf("ReconG findCkov: cluY %.3f fPcY %.3f" ,cluY, fPc.Y());*/
    }

const auto infString4 = Form("ReconG meth: | Actual Ckov : %.3f | Reconstructed Ckov1 = %.3f",etaC, thetaCer);
TCanvas* tCkovGraphG = new TCanvas("T ReconG","etaC Graph ReconG", 1600, 1600);		tCkovGraphG->cd();	
tCkovReconGGraph->SetTitle(infString4);
tCkovReconGGraph->Draw("AP");
      return found;

    } // FindPhotTheta()
   
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    TVector2 traceForward(TVector3 dirCkov) const
    {
    // Trace forward a photon from (x,y) up to PC
    //  Arguments: dirCkov photon vector in LORS
    //    Returns: pos of traced photon at PC

    TVector2 pos(-999, -999);
    double thetaCer = dirCkov.Theta();
    if (thetaCer > TMath::ASin(1. / refIdx)) {
        return pos;
    }                                                       


    // radThick value is value to change for varying L
    // total refraction on WIN-GAP boundary


    // radThick' = (radThick - L)
    // double zRad = - radThick' - 0.5 * winThick;
    double zRad = -0.5 * radThick - 0.5 * winThick;         // z position of middle of RAD
    TVector3 posCkov(fTrkPos.X(), fTrkPos.Y(), zRad);                           // RAD: photon position is track position @ middle of RAD
    propagate(dirCkov, posCkov, -0.5 * winThick);                     // go to RAD-WIN boundary
    refract(dirCkov, refIdx, winIdx);                    // RAD-WIN refraction
    propagate(dirCkov, posCkov, 0.5 * winThick);                      // go to WIN-GAP boundary
    refract(dirCkov, winIdx, gapIdx);                       // WIN-GAP refraction
    propagate(dirCkov, posCkov, 0.5 * winThick + gapThick); // go to PC
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


    void propagate(const TVector3& dir, TVector3& pos, double z) const
    {
    // Finds an intersection point between a line and XY plane shifted along Z.
    // Arguments:  dir,pos   - vector along the line and any point of the line
    //             z         - z coordinate of plain
    //   Returns:  none
    //   On exit:  pos is the position if this intesection if any
    static TVector3 nrm(0, 0, 1);
    TVector3 pnt(0, 0, z);

    TVector3 diff = pnt - pos;
    double sint = 0; //(nrm * diff) / (nrm * dir);
    pos += sint * dir;
    } // Propagate()

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    void refract(TVector3& dir, const double n1, double n2) const
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
  } 
 
  void setTrack(double xRad, double yRad, double theta, double phi)
  {
    fTrkDir.SetMagThetaPhi(1., theta, phi);
    fTrkPos.Set(xRad, yRad);
    Printf("setTrack xRad %.3f yRad %.3f theta %.3f phi %.3f", xRad, yRad, theta, phi);
  } // set track parameter at RAD


  void setImpPC(double xPc, double yPc)
  {
    Printf("setImpPC xPc %.3f yPc %.3f", xPc, yPc);
    fPc.Set(xPc, yPc);
  } // set track impact to PC
  void setMip(double xmip, double ymip)
  {
    fMipPos.Set(xmip, ymip);
  } // set track impact to PC
};
