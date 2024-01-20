#include "SimpleCluster.cpp"



#ifndef Recon_H
#define Recon_H


class Recon 
{
    
    
    public : 
    
     bool findPhotCkov(double cluX, double cluY, double& thetaCer, double& phiCer)
    {
      // Finds Cerenkov angle  for this photon candidate
      // Arguments: cluX,cluY - position of cadidate's cluster
      // Returns: Cerenkov angle

      TVector3 dirCkov;

      double zRad = -0.5 * radThick() - 0.5 * winThick();     // z position of middle of RAD
      TVector3 rad(fTrkPos.X(), fTrkPos.Y(), zRad);                           // impact point at middle of RAD
      TVector3 pc(cluX, cluY, 0.5 * winThick() + gapThick()); // mip at PC
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
    std::vector<SimpleCluster> clusters = {};
    
    Recon(double thetaP, double phiP, double xMip, double yMip, double xRad, double yRad, const std::vector<SimpleCluster>& clustersIn) {
        fTrkDir.SetMagThetaPhi(1., thetaP, phiP);
        fTrkPos.Set(xRad, yRad);
        fMipPos.Set(xMip, yMip);
        fPc.Set(xMip, yMip);
        clusters = clustersIn;
        mThetaP = thetaP;
        mPhiP = phiP;
        
        
    }
    
    double thetaP() const { return mThetaP; }  //<--TEMPORAR--> to be removed in future. Radiator thickness
    double phiP() const { return mPhiP; }  //<--TEMPORAR--> to be removed in future. Window thickness

 
    
    
    
    private :
    TVector3 fTrkDir; // track direction in LORS at RAD

    TVector2 fTrkPos; // track positon in LORS at RAD   // XY mag
    TVector2 fMipPos; // mip positon for a given trackf // XY
    TVector2 fPc;     // track position at PC           // XY   
    double mThetaP,  mPhiP;
    

    double radThick() const { return 1.5; }  //<--TEMPORAR--> to be removed in future. Radiator thickness
    double winThick() const { return 0.5; }  //<--TEMPORAR--> to be removed in future. Window thickness
    double gapThick() const { return 8.0; }  //<--TEMPORAR--> to be removed in future. Proximity gap thickness
    double winIdx() const { return 1.583; } //<--TEMPORAR--> to be removed in future. Mean refractive index of WIN material (SiO2)
    double gapIdx() const { return 1.0005; } //<--TEMPORAR--> to be removed in future. Mean refractive index of GAP material (CH4)
    double getRefIdx() const { return 1.2904; } //<--TEMPORAR--> to be removed in future. Mean refractive index of GAP material (CH4)
    

    /*
    def winIdx(self):
        return np.sqrt(1 + 46.411 / (10.666 * 10.666 - self.eV * self.eV) + 228.71 / (18.125 * 18.125 - self.eV * self.eV))

    def gapIdx(self):
        return 1 + 0.12489e-6 / (2.62e-4 - self.eV * self.eV / 1239.84 / 1239.84)

    # Assuming std::sqrt is the square root function
    def nIdxRad(self, eV, temp=25):  # Default temperature is 20 unless provided otherwise
        eV_term = 1239.84 / self.eV
        return np.sqrt(1 + 0.554 * eV_term * eV_term / (eV_term * eV_term - 5769)) - 0.0005 * (temp - 20) */
    
    
    
    void propagate(const TVector3 dir, TVector3& pos, double z) const
    {
      // Finds an intersection point between a line and XY plane shifted along Z.
      // Arguments:  dir,pos   - vector along the line and any point of the line
      //             z         - z coordinate of plain
      //   Returns:  none
      //   On exit:  pos is the position if this intesection if any
      static TVector3 nrm(0, 0, 1);
      TVector3 pnt(0, 0, z);

      TVector3 diff = pnt - pos;
      double sint = (nrm * diff) / (nrm * dir);
      pos += sint * dir;
    } // Propagate()
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    void refract(TVector3& dir, double n1, double n2) const
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

   

    TVector2 traceForward(TVector3 dirCkov) const
    {
      // Trace forward a photon from (x,y) up to PC
      //  Arguments: dirCkov photon vector in LORS
      //    Returns: pos of traced photon at PC

      TVector2 pos(-999, -999);
      double thetaCer = dirCkov.Theta();
      if (thetaCer > TMath::ASin(1. / getRefIdx())) {
        return pos;
      }                                                                           // total refraction on WIN-GAP boundary
      double zRad = -0.5 * radThick() - 0.5 * winThick();         // z position of middle of RAD
      TVector3 posCkov(fTrkPos.X(), fTrkPos.Y(), zRad);                           // RAD: photon position is track position @ middle of RAD
      propagate(dirCkov, posCkov, -0.5 * winThick());                     // go to RAD-WIN boundary
      refract(dirCkov, getRefIdx(), winIdx());                    // RAD-WIN refraction
      propagate(dirCkov, posCkov, 0.5 * winThick());                      // go to WIN-GAP boundary
      refract(dirCkov, winIdx(), gapIdx());                       // WIN-GAP refraction
      propagate(dirCkov, posCkov, 0.5 * winThick() + gapThick()); // go to PC
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


  /*
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

    fParam->setRefIdx(nmean);

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
      if (cluster.q() > 2 * fParam->qCut() || cluster.size() > 4) {
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
    if (fPhotCnt > fParam->multCut()) {
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
  } // CkovAngle() 
  */

};    
   
#endif
