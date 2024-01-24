

!pip install cppyy
import cppyy

Alisigma2="""
#ifndef Alisigma2_H
#define Alisigma2_H
class Alisigma2_ {
public:
    Alisigma2_(double refIdx) {
        mRefIdx = refIdx;
    }

    // Member function declarations
    double sigma2(double trkTheta, double trkPhi, double ckovTh, double ckovPh);
    double sigLoc(double trkTheta, double trkPhi, double thetaC, double phiC, double betaM);
    double sigCrom(double trkTheta, double trkPhi, double thetaC, double phiC, double betaM);
    double sigGeom(double trkTheta, double trkPhi, double thetaC, double phiC, double betaM);
    double sigmaCorrFact(int iPart, double occupancy);

    double getRefIdx() const { return mRefIdx; } // running refractive index
    double mRefIdx = 1.2904; // running refractive index

    double radThick() const { return 1.5; }  //<--TEMPORAR--> to be removed in future. Radiator thickness
    double winThick() const { return 0.5; }  //<--TEMPORAR--> to be removed in future. Window thickness
    double gapThick() const { return 8.0; }  //<--TEMPORAR--> to be removed in future. Proximity gap thickness

    double winIdx() const { return 1.5833; } //<--TEMPORAR--> to be removed in future. Mean refractive index of WIN material (SiO2)

    //double winIdx() const { return 1.5787; } //<--TEMPORAR--> to be removed in future. Mean refractive index of WIN material (SiO2)
    double gapIdx() const { return 1.0005; } //<--TEMPORAR--> to be removed in future. Mean refractive index of GAP material (CH4)


};

double Alisigma2_::sigma2(double trkTheta, double trkPhi, double ckovTh, double ckovPh)
{
  // Analithical calculation of total error (as a sum of localization, geometrical and chromatic errors)
  // on Cerenkov angle for a given Cerenkov photon
  // created by a given MIP. Fromulae according to CERN-EP-2000-058
  // Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]
  //            dip and azimuthal angles for MIP taken at the entrance to radiator, [radians]
  //            MIP beta
  //   Returns: absolute error on Cerenkov angle, [radians]

  double trkBeta = 1. / (std::cos(ckovTh) * getRefIdx());

  if (trkBeta > 1) {
    trkBeta = 1; // protection against bad measured thetaCer
  }
  if (trkBeta < 0) {
    trkBeta = 0.0001; //
  }

  double sigLocVar = sigLoc(trkTheta, trkPhi, ckovTh, ckovPh, trkBeta);
  double sigGeomVar =sigGeom(trkTheta, trkPhi, ckovTh, ckovPh, trkBeta);
  double sigCromVar = sigCrom(trkTheta, trkPhi, ckovTh, ckovPh, trkBeta);
  double sigmaRes2 = (sigLocVar*sigLocVar + sigGeomVar*sigGeomVar + sigCromVar*sigCromVar);
  return sigmaRes2;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double Alisigma2_::sigLoc(double trkTheta, double trkPhi, double thetaC, double phiC, double betaM)
{
  // Analitical calculation of localization error (due to finite segmentation of PC) on Cerenkov angle for a given
  // Cerenkov photon
  // created by a given MIP. Fromulae according to CERN-EP-2000-058
  // Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]
  //            dip and azimuthal angles for MIP taken at the entrance to radiator, [radians]
  //            MIP beta
  //   Returns: absolute error on Cerenkov angle, [radians]

  double phiDelta = phiC - trkPhi;

  double sint = std::sin(trkTheta);
  double cost = std::cos(trkTheta);
  double sinf = std::sin(trkPhi);
  double cosf = std::cos(trkPhi);
  double sinfd = std::sin(phiDelta);
  double cosfd = std::cos(phiDelta);
  double tantheta = std::tan(thetaC);

  double alpha = cost - tantheta * cosfd * sint;                               // formula (11)
  double k = 1. - getRefIdx() * getRefIdx() + alpha * alpha / (betaM * betaM); // formula (after 8 in the text)
  if (k < 0) {
    return 1e10;
  }
  double mu = sint * sinf + tantheta * (cost * cosfd * sinf + sinfd * cosf); // formula (10)
  double e = sint * cosf + tantheta * (cost * cosfd * cosf - sinfd * sinf);  // formula (9)

  double kk = betaM * std::sqrt(k) / (gapThick() * alpha); // formula (6) and (7)
  // formula (6)
  double dtdxc = kk * (k * (cosfd * cosf - cost * sinfd * sinf) - (alpha * mu / (betaM * betaM)) * sint * sinfd);
  // formula (7)            pag.4
  double dtdyc = kk * (k * (cosfd * sinf + cost * sinfd * cosf) + (alpha * e / (betaM * betaM)) * sint * sinfd);

  double errX = 0.2, errY = 0.25; // end of page 7
  return std::sqrt(errX * errX * dtdxc * dtdxc + errY * errY * dtdyc * dtdyc);
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double Alisigma2_::sigCrom(double trkTheta, double trkPhi, double thetaC, double phiC, double betaM)
{
  // Analitical calculation of chromatic error (due to lack of knowledge of Cerenkov photon energy)
  // on Cerenkov angle for a given Cerenkov photon
  // created by a given MIP. Fromulae according to CERN-EP-2000-058
  // Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]
  //            dip and azimuthal angles for MIP taken at the entrance to radiator, [radians]
  //            MIP beta
  //   Returns: absolute error on Cerenkov angle, [radians]

  double phiDelta = phiC - trkPhi;

  double sint = std::sin(trkTheta);
  double cost = std::cos(trkTheta);
  double cosfd = std::cos(phiDelta);
  double tantheta = std::tan(thetaC);

  double alpha = cost - tantheta * cosfd * sint;                         // formula (11)
  double dtdn = cost * getRefIdx() * betaM * betaM / (alpha * tantheta); // formula (12)

  //  double f = 0.00928*(7.75-5.635)/std::sqrt(12.);
  double f = 0.0172 * (7.75 - 5.635) / std::sqrt(24.);

  return f * dtdn;
} // SigCrom()
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double Alisigma2_::sigGeom(double trkTheta, double trkPhi, double thetaC, double phiC, double betaM)
{
  // Analitical calculation of geometric error (due to lack of knowledge of creation point in radiator)
  // on Cerenkov angle for a given Cerenkov photon
  // created by a given MIP. Formulae according to CERN-EP-2000-058
  // Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]
  //            dip and azimuthal angles for MIP taken at the entrance to radiator, [radians]
  //            MIP beta
  //   Returns: absolute error on Cerenkov angle, [radians]

  double phiDelta = phiC - trkPhi;

  double sint = std::sin(trkTheta);
  double cost = std::cos(trkTheta);
  double sinf = std::sin(trkPhi);
  double cosfd = std::cos(phiDelta);
  double costheta = std::cos(thetaC);
  double tantheta = std::tan(thetaC);

  double alpha = cost - tantheta * cosfd * sint; // formula (11)

  double k = 1. - getRefIdx() * getRefIdx() + alpha * alpha / (betaM * betaM); // formula (after 8 in the text)
  if (k < 0) {
    return 1e10;
  }

  double eTr = 0.5 * radThick() * betaM * std::sqrt(k) / (gapThick() * alpha); // formula (14)
  double lambda = (1. - sint * sinf) * (1. + sint * sinf);                       // formula (15)

  double c1 = 1. / (1. + eTr * k / (alpha * alpha * costheta * costheta));                     // formula (13.a)
  double c2 = betaM * std::pow(k, 1.5) * tantheta * lambda / (gapThick() * alpha * alpha); // formula (13.b)



  double c3 = (1. + eTr * k * betaM * betaM) / ((1 + eTr) * alpha * alpha);                    // formula (13.c)
  double c4 = std::sqrt(k) * tantheta * (1 - lambda) / (gapThick() * betaM);                 // formula (13.d)
  double dtdT = c1 * (c2 + c3 * c4);
  double trErr = radThick() / (std::sqrt(12.) * cost);

  return trErr * dtdT;
} // SigGeom()
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double Alisigma2_::sigmaCorrFact(int iPart, double occupancy)
{
  double corr = 1.0;

  switch (iPart) {
    case 0:
      corr = 0.115 * occupancy + 1.166;
      break;
    case 1:
      corr = 0.115 * occupancy + 1.166;
      break;
    case 2:
      corr = 0.115 * occupancy + 1.166;
      break;
    case 3:
      corr = 0.065 * occupancy + 1.137;
      break;
    case 4:
      corr = 0.048 * occupancy + 1.202;
      break;
  }

  return corr;
}
#endif // Alisigma2_H
"""
cpp_recon = """
#include <Eigen/Dense>

#include <math.h>

class Recon {
    public:

        Recon () {}
        double pi = std::atan(1)*4;




        /*inline double asin(double x)
        { if (x < -1.) return -pi/2;
            if (x >  1.) return  pi/2;
            return std::asin(x);
        }*/


        double mTrackPhi = 0;
        double mTrackTheta = 0;
        double mTrackX = 0;
        double mTrackY = 0;

        Eigen::Vector3d fTrkDir; // Track direction as a 3D vector
        Eigen::Vector2d fTrkPos; // Track position as a 2D vector
        Eigen::Vector2d fMipPos; // Mip position as a 2D vector
        Eigen::Vector2d fPc;     // PC position as a 2D vector



        TVector3_ fTrkDir; // Just for test
        TVector2_ fTrkPos;  // Just for test
        TVector2_ fMipPos;
        TVector2_ fPc;
        double getRefIdx() const { return mRefIdx; } // running refractive index
        double mRefIdx = 1.2904; // running refractive index

        double radThick() const { return 1.5; }  //<--TEMPORAR--> to be removed in future. Radiator thickness
        double winThick() const { return 0.5; }  //<--TEMPORAR--> to be removed in future. Window thickness
        double gapThick() const { return 8.0; }  //<--TEMPORAR--> to be removed in future. Proximity gap thickness

        double winIdx() const { return 1.5833; } //<--TEMPORAR--> to be removed in future. Mean refractive index of WIN material (SiO2)

        //double winIdx() const { return 1.5787; } //<--TEMPORAR--> to be removed in future. Mean refractive index of WIN material (SiO2)
        double gapIdx() const { return 1.0005; } //<--TEMPORAR--> to be removed in future. Mean refractive index of GAP material (CH4)
/*
using MatchInfo = o2::dataformats::MatchInfoHMP;

using namespace o2::hmpid;
// ClassImp(Recon);


class Recon {


    Recon () {}
    double getRefIdx() const { return mRefIdx; } // running refractive index
    double mRefIdx = 1.2904; // running refractive index

    double radThick() const { return 1.5; }  //<--TEMPORAR--> to be removed in future. Radiator thickness
    double winThick() const { return 0.5; }  //<--TEMPORAR--> to be removed in future. Window thickness
    double gapThick() const { return 8.0; }  //<--TEMPORAR--> to be removed in future. Proximity gap thickness

    double winIdx() const { return 1.5833; } //<--TEMPORAR--> to be removed in future. Mean refractive index of WIN material (SiO2)

    //double winIdx() const { return 1.5787; } //<--TEMPORAR--> to be removed in future. Mean refractive index of WIN material (SiO2)
    double gapIdx() const { return 1.0005; } //<--TEMPORAR--> to be removed in future. Mean refractive index of GAP material (CH4)



    ClassImp(o2::hmpid::Recon);
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    void initVars(int n)
    {
    //..
    // Init some variables
    //..
    if (n <= 0) {
        return;
    }

    // ef : changed to smart-pointer Array
    // fPhotFlag = new int[n];
    fPhotFlag = std::unique_ptr<int[]>(new int[n]);
    fPhotClusIndex = std::unique_ptr<int[]>(new int[n]);

    fPhotCkov = std::unique_ptr<double[]>(new double[n]);
    fPhotPhi = std::unique_ptr<double[]>(new double[n]);
    fPhotWei = std::unique_ptr<double[]>(new double[n]);
    //
    }

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    void ckovAngle(o2::dataformats::MatchInfoHMP* match, const std::vector<o2::hmpid::Cluster> clusters, int index, double nmean, float xRa, float yRa)
    {
    // Pattern recognition method based on Hough transform
    // Arguments:   pTrk     - track for which Ckov angle is to be found
    //              pCluLst  - list of clusters for this chamber
    //   Returns:            - track ckov angle, [rad],

    const int nMinPhotAcc = 3; // Minimum number of photons required to perform the pattern recognition

    int nClusTot = clusters.size();

    initVars(nClusTot);

    float xPc, yPc, th, ph;

    match->getHMPIDtrk(xPc, yPc, th, ph); // initialize this track: th and ph angles at middle of RAD

    setTrack(xRa, yRa, th, ph);

    setRefIdx(nmean);

    float mipQ = -1, mipX = -1, mipY = -1;
    int chId = -1, sizeClu = -1;

    fPhotCnt = 0;

    int nPads = 0;

    for (int iClu = 0; iClu < clusters.size(); iClu++) { // clusters loop

        o2::hmpid::Cluster cluster = clusters.at(iClu);
        nPads += cluster.size();
        if (iClu == index) { // this is the MIP! not a photon candidate: just store mip info
        mipX = cluster.x();
        mipY = cluster.y();
        mipQ = cluster.q();
        sizeClu = cluster.size();
        continue;
        }

        chId = cluster.ch();
        if (cluster.q() > 2 * qCut() || cluster.size() > 4) {
        continue;
        }
        double thetaCer, phiCer;
        if (findPhotCkov(cluster.x(), cluster.y(), thetaCer, phiCer)) { // find ckov angle for this  photon candidate
        fPhotCkov[fPhotCnt] = thetaCer;                               // actual theta Cerenkov (in TRS)
        fPhotPhi[fPhotCnt] = phiCer;
        fPhotClusIndex[fPhotCnt] = iClu; // actual phi   Cerenkov (in TRS): -pi to come back to "unusual" ref system (X,Y,-Z)
        fPhotCnt++;                      // increment counter of photon candidates
        }
    } // clusters loop

    match->setHMPIDmip(mipX, mipY, mipQ, fPhotCnt);     // store mip info in any case
    match->setIdxHMPClus(chId, index + 1000 * sizeClu); // set index of cluster
    match->setMipClusSize(sizeClu);

    if (fPhotCnt < nMinPhotAcc) {         // no reconstruction with <=3 photon candidates
        match->setHMPsignal(kNoPhotAccept); // set the appropriate flag
        return;
    }

    fMipPos.Set(mipX, mipY);

    // PATTERN RECOGNITION STARTED:
    if (fPhotCnt > multCut()) {
        fIsWEIGHT = kTRUE;
    } // offset to take into account bkg in reconstruction
    else {
        fIsWEIGHT = kFALSE;
    }

    float photCharge[10] = {0x0};

    int iNrec = flagPhot(houghResponse(), clusters, photCharge); // flag photons according to individual theta ckov with respect to most probable
    // int iNrec = flagPhot(houghResponse(), clusters); // flag photons according to individual theta ckov with respect to most probable

    match->setPhotCharge(photCharge);
    match->setHMPIDmip(mipX, mipY, mipQ, iNrec); // store mip info

    if (iNrec < nMinPhotAcc) {
        match->setHMPsignal(kNoPhotAccept); // no photon candidates are accepted
        return;
    }

    int occupancy = (int)(1000 * (nPads / (6. * 80. * 48.)));

    double thetaC = findRingCkov(clusters.size()); // find the best reconstructed theta Cherenkov
    findRingGeom(thetaC, 2);

    match->setHMPsignal(thetaC + occupancy); // store theta Cherenkov and chmaber occupancy
    // match->SetHMPIDchi2(fCkovSigma2);                                                        //store experimental ring angular resolution squared

    // deleteVars(); ef : in case of smart-pointers, should not be necessary?
    } // CkovAngle() */
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    bool findPhotCkov(double cluX, double cluY, double& thetaCer, double& phiCer)
    {
        // Finds Cerenkov angle  for this photon candidate
        // Arguments: cluX,cluY - position of cadidate's cluster
        // Returns: Cerenkov angle

        TVector3_ dirCkov;

        double zRad = -0.5 * radThick() - 0.5 * winThick();     // z position of middle of RAD
        TVector3_ rad(fTrkPos.X(), fTrkPos.Y(), zRad);                           // impact point at middle of RAD
        TVector3_ pc(cluX, cluY, 0.5 * winThick() + gapThick()); // mip at PC
        double cluR = std::sqrt((cluX - fPc.X()) * (cluX - fPc.X()) +
                                    (cluY - fPc.Y()) * (cluY - fPc.Y())); // ref. distance impact RAD-CLUSTER
        double phi = (pc - rad).Phi();                                  // phi of photon

        double ckov1 = 0;
        double ckov2 = 0.75 + fTrkDir.Theta(); // start to find theta cerenkov in DRS
        const double kTol = 0.01;
        int iIterCnt = 0;
        while (1) {
            if (iIterCnt >= 50) {
                return false;
            }
            double ckov = 0.5 * (ckov1 + ckov2);
            dirCkov.SetMagThetaPhi(1, ckov, phi);
            TVector2_ posC = traceForward(dirCkov);   // trace photon with actual angles
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
                return true;
            }
        }
    } // FindPhotTheta()
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    TVector2_ traceForward(TVector3_ dirCkov) const
    {
        // Trace forward a photon from (x,y) up to PC
        //  Arguments: dirCkov photon vector in LORS
        //    Returns: pos of traced photon at PC

        TVector2_ pos(-999, -999);
        double thetaCer = dirCkov.Theta();
        if (thetaCer > std::asin(1. / getRefIdx())) {
            return pos;
        }                                                                           // total refraction on WIN-GAP boundary
        double zRad = -0.5 * radThick() - 0.5 * winThick();         // z position of middle of RAD
        TVector3_ posCkov(fTrkPos.X(), fTrkPos.Y(), zRad);                           // RAD: photon position is track position @ middle of RAD
        propagate(dirCkov, posCkov, -0.5 * winThick());                     // go to RAD-WIN boundary
        refract(dirCkov, getRefIdx(), winIdx());                    // RAD-WIN refraction
        propagate(dirCkov, posCkov, 0.5 * winThick());                      // go to WIN-GAP boundary
        refract(dirCkov, winIdx(), gapIdx());                       // WIN-GAP refraction
        propagate(dirCkov, posCkov, 0.5 * winThick() + gapThick()); // go to PC
        pos.Set(posCkov.X(), posCkov.Y());
        return pos;

    } // TraceForward()

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    void lors2Trs(TVector3_ dirCkov, double& thetaCer, double& phiCer) const
    {
        // Theta Cerenkov reconstruction
        //  Arguments: dirCkov photon vector in LORS
        //    Returns: thetaCer of photon in TRS
        //               phiCer of photon in TRS
        //  TVector3_ dirTrk;
        //  dirTrk.SetMagThetaPhi(1,fTrkDir.Theta(),fTrkDir.Phi()); -> dirTrk.SetCoordinates(1,fTrkDir.Theta(),fTrkDir.Phi())
        //  double  thetaCer = TMath::ACos(dirCkov*dirTrk);

        TRotation mtheta;
        mtheta.RotateY(-fTrkDir.Theta());

        TRotation mphi;
        mphi.RotateZ(-fTrkDir.Phi());

        TRotation mrot = mtheta * mphi;

        TVector3_ dirCkovTRS;
        dirCkovTRS = mrot * dirCkov;
        phiCer = dirCkovTRS.Phi();     // actual value of the phi of the photon
        thetaCer = dirCkovTRS.Theta(); // actual value of thetaCerenkov of the photon
    }
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    void trs2Lors(TVector3_ dirCkov, double& thetaCer, double& phiCer) const
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

        TVector3_ dirCkovLORS;
        dirCkovLORS = mrot * dirCkov;

        phiCer = dirCkovLORS.Phi();     // actual value of the phi of the photon
        thetaCer = dirCkovLORS.Theta(); // actual value of thetaCerenkov of the photon
    }


    /*
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    void findRingGeom(double ckovAng, int level)
    {
        // Find area covered in the PC acceptance
        // Arguments: ckovAng - cerenkov angle
        //            level   - precision in finding area and portion of ring accepted (multiple of 50)
        //   Returns: area of the ring in cm^2 for given theta ckov

        Int_t kN = 50 * level;
        Int_t nPoints = 0;
        doublearea = 0;

        Bool_t first = kFALSE;
        TVector2_ pos1;

        for (Int_t i = 0; i < kN; i++) {
            if (!first) {
            pos1 = tracePhot(ckovAng, double(pi * 2 * (i + 1) / kN)); // find a good trace for the first photon
            if (pos1.X() == -999) {
                continue;
            } // no area: open ring
            if (!isInside(pos1.X(), pos1.Y(), 0)) {
                pos1 = intWithEdge(fMipPos, pos1); // find the very first intersection...
            } else {
                if (!isInDead(pos1.X(), pos1.Y())) {
                nPoints++;
                } // photon is accepted if not in dead zone
            }
            first = kTRUE;
            continue;
            }
            TVector2_ pos2 = tracePhot(ckovAng, double(pi * 2 * (i + 1) / kN)); // trace the next photon
            if (pos2.X() == -999) {
            {
                continue;
            }
            } // no area: open ring
            if (!isInside(pos2.X(), pos2.Y(), 0)) {
            pos2 = intWithEdge(fMipPos, pos2);
            } else {
            if (!isInDead(pos2.X(), pos2.Y())) {
                nPoints++;
            } // photon is accepted if not in dead zone
            }
            area += std::abs((pos1 - fMipPos).X() * (pos2 - fMipPos).Y() - (pos1 - fMipPos).Y() * (pos2 - fMipPos).X()); // add area of the triangle...
            pos1 = pos2;
        }
        //---  find area and length of the ring;
        fRingAcc = (Double_t)nPoints / (Double_t)kN;
        area *= 0.5;
        fRingArea = area;

    } // FindRingGeom()
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const TVector2_ intWithEdge(TVector2_ p1, TVector2_ p2)
    {
        // It finds the intersection of the line for 2 points traced as photons
        // and the edge of a given PC
        // Arguments: 2 points obtained tracing the photons
        //   Returns: intersection point with detector (PC) edges

        double xmin = (p1.X() < p2.X()) ? p1.X() : p2.X();
        double xmax = (p1.X() < p2.X()) ? p2.X() : p1.X();
        double ymin = (p1.Y() < p2.Y()) ? p1.Y() : p2.Y();
        double ymax = (p1.Y() < p2.Y()) ? p2.Y() : p1.Y();

        double m = std::tan((p2 - p1).Phi());
        TVector2_ pint;
        // intersection with low  X
        pint.Set((double)(p1.X() + (0 - p1.Y()) / m), 0.);
        if (pint.X() >= 0 && pint.X() <= sizeAllX() &&
            pint.X() >= xmin && pint.X() <= xmax &&
            pint.Y() >= ymin && pint.Y() <= ymax) {
            return pint;
        }
        // intersection with high X
        pint.Set((double)(p1.X() + (sizeAllY() - p1.Y()) / m), (double)(sizeAllY()));
        if (pint.X() >= 0 && pint.X() <= sizeAllX() &&
            pint.X() >= xmin && pint.X() <= xmax &&
            pint.Y() >= ymin && pint.Y() <= ymax) {
            return pint;
        }
        // intersection with left Y
        pint.Set(0., (double)(p1.Y() + m * (0 - p1.X())));
        if (pint.Y() >= 0 && pint.Y() <= sizeAllY() &&
            pint.Y() >= ymin && pint.Y() <= ymax &&
            pint.X() >= xmin && pint.X() <= xmax) {
            return pint;
        }
        // intersection with righ Y
        pint.Set((double)(sizeAllX()), (double)(p1.Y() + m * (sizeAllX() - p1.X()))); // ef: Set->SetCoordinates
        if (pint.Y() >= 0 && pint.Y() <= sizeAllY() &&
            pint.Y() >= ymin && pint.Y() <= ymax &&
            pint.X() >= xmin && pint.X() <= xmax) {
            return pint;
        }
        return p1;
    } // IntWithEdge()
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    double findRingCkov(int)
    {
        // Loops on all Ckov candidates and estimates the best Theta Ckov for a ring formed by those candidates. Also estimates an error for that Theat Ckov
        // collecting errors for all single Ckov candidates thetas. (Assuming they are independent)
        // Arguments: iNclus- total number of clusters in chamber for background estimation
        //    Return: best estimation of track Theta ckov

        doublewei = 0.;
        doubleweightThetaCerenkov = 0.;

        doubleckovMin = 9999., ckovMax = 0.;
        doublesigma2 = 0; // to collect error squared for this ring

        for (Int_t i = 0; i < fPhotCnt; i++) { // candidates loop
            if (fPhotFlag[i] == 2) {
            if (fPhotCkov[i] < ckovMin) {
                ckovMin = fPhotCkov[i];
            } // find max and min Theta ckov from all candidates within probable window
            if (fPhotCkov[i] > ckovMax) {
                ckovMax = fPhotCkov[i];
            }
            weightThetaCerenkov += fPhotCkov[i] * fPhotWei[i];
            wei += fPhotWei[i]; // collect weight as sum of all candidate weghts

            sigma2 += 1. / sigma2(fTrkDir.Theta(), fTrkDir.Phi(), fPhotCkov[i], fPhotPhi[i]);
            }
        } // candidates loop

        if (sigma2 > 0) {
            fCkovSigma2 = 1. / sigma2;
        } else {
            fCkovSigma2 = 1e10;
        }

        if (wei != 0.) {
            weightThetaCerenkov /= wei;
        } else {
            weightThetaCerenkov = 0.;
        }
        return weightThetaCerenkov;

    } // FindCkovRing()
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    int flagPhot(double ckov, const std::vector<o2::hmpid::Cluster> clusters, float* photChargeVec)
    // int flagPhot(double ckov, const std::vector<o2::hmpid::Cluster> clusters)
    {
        // Flag photon candidates if their individual ckov angle is inside the window around ckov angle returned by  HoughResponse()
        // Arguments: ckov- value of most probable ckov angle for track as returned by HoughResponse()
        //   Returns: number of photon candidates happened to be inside the window

        // Photon Flag:  Flag = 0 initial set;
        //               Flag = 1 good candidate (charge compatible with photon);
        //               Flag = 2 photon used for the ring;

        Int_t steps = (Int_t)((ckov) / fDTheta); // how many times we need to have fDTheta to fill the distance between 0  and thetaCkovHough

        doubletmin = (Double_t)(steps - 1) * fDTheta;
        doubletmax = (Double_t)(steps)*fDTheta;
        doubletavg = 0.5 * (tmin + tmax);

        tmin = tavg - 0.5 * fWindowWidth;
        tmax = tavg + 0.5 * fWindowWidth;

        Int_t iInsideCnt = 0;                  // count photons which Theta ckov inside the window
        for (Int_t i = 0; i < fPhotCnt; i++) { // photon candidates loop
            fPhotFlag[i] = 0;
            if (fPhotCkov[i] >= tmin && fPhotCkov[i] <= tmax) {
            fPhotFlag[i] = 2;
            o2::hmpid::Cluster cluster = clusters.at(fPhotClusIndex[i]);
            float charge = cluster.q();
            if (iInsideCnt < 10) {
                photChargeVec[iInsideCnt] = charge;
            } // AddObjectToFriends(pCluLst,i,pTrk);
            iInsideCnt++;
            }
        }

        return iInsideCnt;

    } // FlagPhot() */
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    TVector2_ tracePhot(double ckovThe, double ckovPhi) const
    {
        // Trace a single Ckov photon from emission point somewhere in radiator up to photocathode taking into account ref indexes of materials it travereses
        // Arguments: ckovThe,ckovPhi- photon ckov angles in TRS, [rad]
        //   Returns: distance between photon point on PC and track projection

        double theta, phi;
        TVector3_ dirTRS, dirLORS;
        dirTRS.SetMagThetaPhi(1, ckovThe, ckovPhi); // photon in TRS
        trs2Lors(dirTRS, theta, phi);
        dirLORS.SetMagThetaPhi(1, theta, phi); // photon in LORS
        return traceForward(dirLORS);          // now foward tracing

    } // tracePhot()
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    void propagate(const TVector3_ dir, TVector3_& pos, double z) const
    {
        // Finds an intersection point between a line and XY plane shifted along Z.
        // Arguments:  dir,pos   - vector along the line and any point of the line
        //             z         - z coordinate of plain
        //   Returns:  none
        //   On exit:  pos is the position if this intesection if any
        static TVector3_ nrm(0, 0, 1);
        TVector3_ pnt(0, 0, z);

        TVector3_ diff = pnt - pos;
        double sint = (nrm * diff) / (nrm * dir);
        pos += sint * dir;
    } // Propagate()
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    void refract(TVector3_& dir, double n1, double n2) const
    {
        // Refract direction vector according to Snell law
        // Arguments:
        //            n1 - ref idx of first substance
        //            n2 - ref idx of second substance
        //   Returns: none
        //   On exit: dir is new direction

        double sinref = (n1 / n2) * std::sin(dir.Theta());
        if (std::abs(sinref) > 1.) {
            dir.SetXYZ(-999, -999, -999);
        } else {
            dir.SetTheta(std::asin(sinref)); //            dir.SetTheta(TMath::ASin(sinref));
        }
    }
};
"""


  # cppyy.cppdef(Alisigma2)

# # Create an instance of the Param class
# Alisigma2_ = cppyy.gbl.Alisigma2_(1.2904)
# # Particle masses in GeV/c^2
# import multiprocessing
# def process_track(params):
#     recon, cluX, cluY, thetaP, phiP = params
#     # Call findPhotCkov and receive updated thetaCer and phiCer
#     updated_thetaCer, updated_phiCer = recon.findPhotCkov(cluX, cluY)



#     Alisigma2_ = cppyy.gbl.Alisigma2_(1.2904)

#     #Double_t aliSigma2(Double_t trkTheta,Double_t trkPhi, Double_t ckovTh, Double_t ckovPh)
#     # Call the sigLoc method

#     print(f"updated_thetaCer {updated_thetaCer} updated_phiCer {updated_phiCer}")

#     sigma2 = Alisigma2_.sigma2(thetaP, phiP,  updated_thetaCer, updated_phiCer)
#     print(f"sigma2 {sigma2}")

#     return updated_thetaCer, updated_phiCer, sigma2


# #theta_cer_padded, phi_cer_padded = parallel_process(recon, x_padded, y_padded, self.ThetaP[0], self.PhiP[0])

# def parallel_process(recon, x_padded, y_padded, thetaP, phiP):
#     params_list = [(recon, cluX, cluY, thetaP, phiP) for cluX, cluY in zip(x_padded, y_padded)]

#     with multiprocessing.Pool() as pool:
#         updated_values = pool.map(process_track, params_list)

#     # Unpack and update the theta_cer_padded and phi_cer_padded arrays
#     theta_cer_padded = []
#     phi_cer_padded = []
#     for updated_thetaCer, updated_phiCer in updated_values:
#         theta_cer_padded.append(updated_thetaCer)
#         phi_cer_padded.append(updated_phiCer)

#     return theta_cer_padded, phi_cer_padded

# def get_other_tracks(event_data_dict, x_padded, y_padded):

#     print(f"called get_other_tracks(event_data_dict, x_padded, y_padded)")



#     class ClassForPhotonProbability:

#         def __init__(self):
#             dict = ["Momentum", "xRad", "yRad", "xMip", "yMip", "ThetaP", "PhiP"]

#             # Initialize all attributes to None or a default value
#             self.Momentum = None
#             self.xRad = None
#             self.yRad = None
#             self.xMip = None
#             self.yMip = None
#             self.ThetaP = None
#             self.PhiP = None

#             self.fTrkPos = [0, 0]
#             self.fMipPos = [0, 0]
#             #self.fTrkDir = [self.xMip, self.yMip]



#         def populate_attributes(self):
#             # Iterate over each key and fill the class attributes with corresponding values

#             setattr(self, "Momentum", event_data_dict["Momentum"])
#             setattr(self, "xRad", event_data_dict["xRad"])
#             setattr(self, "yRad", event_data_dict["yRad"])
#             setattr(self, "xMip", event_data_dict["xMip"])
#             setattr(self, "yMip", event_data_dict["yMip"])
#             setattr(self, "ThetaP", event_data_dict["ThetaP"])
#             setattr(self, "PhiP", event_data_dict["PhiP"])




#             self.fTrkPos = [self.xRad[0], self.yRad[1]]
#             self.fMipPos = [self.xMip[0], self.yMip[1]]

#             recon = Recon(self.fTrkPos, self.fMipPos, self.ThetaP[0])

#             phiCer, thetaCer = 0, 0



#             # Example usage
#             theta_cer_padded, phi_cer_padded = parallel_process(recon, x_padded, y_padded, self.ThetaP[0], self.PhiP[0])


#             # x_padded, y_padded

#             #Double_t aliSigma2(Double_t trkTheta,Double_t trkPhi, Double_t ckovTh, Double_t ckovPh)
#             # Call the sigLoc method


#             #result = Alisigma2_.sigLoc(1.0, 0.5, 0.3, 0.4, 0.9)
#             print(f"shape theta_cer_padded {theta_cer_padded.shape}")


#     photon_prob = ClassForPhotonProbability()
#     print("photon_prob = ClassForPhotonProbability()")

#     # Populate the attributes using event_data_dict
#     photon_prob.populate_attributes()



#     #return theta_c_hyps



# cppyy.cppdef(Alisigma2)

# # Create an instance of the Param class
# Alisigma2_ = cppyy.gbl.Alisigma2_(1.2904)
# # Particle masses in GeV/c^2
