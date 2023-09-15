
template <typename T>
bool Recon::findPhotCkov(double cluX, double cluY, double& thetaCer, double& phiCer)
{
  // Finds Cerenkov angle  for this photon candidate
  // Arguments: cluX,cluY - position of cadidate's cluster
  // Returns: Cerenkov angle

  Polar3DVector dirCkov;

  double zRad = -0.5 * fParam->radThick() - 0.5 * fParam->winThick(); // z position of middle of RAD

  math_utils::Vector3D<double> rad(fTrkPos.X(), fTrkPos.Y(), zRad);                         // impact point at middle of RAD
  math_utils::Vector3D<double> pc(cluX, cluY, 0.5 * fParam->winThick() + fParam->gapIdx()); // mip at PC

  double cluR = TMath::Sqrt((cluX - fTrkPos.X()) * (cluX - fTrkPos.X()) +
                            (cluY - fTrkPos.Y()) * (cluY - fTrkPos.Y())); // ref. distance impact RAD-CLUSTER

  double phi = (pc - rad).Phi(); // phi of photon

  double ckov1 = 0;
  double ckov2 = 0.75 + fTrkDir.Theta(); // start to find theta cerenkov in DRS
  const double kTol = 0.01;
  int iIterCnt = 0;

  while (1) {
    if (iIterCnt >= 50) {
      return kFALSE;
    }
    double ckov = 0.5 * (ckov1 + ckov2);

    dirCkov.SetCoordinates(1, ckov, phi);
    o2::math_utils::Vector2D<double> posC = traceForward(dirCkov); // trace photon with actual angles
    double dist = cluR - (posC - fTrkPos).Mag2();                  // get distance between trial point and cluster position

    if (posC.X() == -999) {
      dist = -999; // total reflection problem
    }
    iIterCnt++;    // counter step
    if (dist > kTol) {
      ckov1 = ckov; // cluster @ larger ckov
    } else if (dist < -kTol) {
      ckov2 = ckov; // cluster @ smaller ckov
    } else {          // precision achived: ckov in DRS found

      dirCkov.SetCoordinates(1, ckov, phi); 
      lors2Trs(dirCkov, thetaCer, phiCer);  // find ckov (in TRS:the effective Cherenkov angle!)
      return kTRUE;
    }
  }
} // FindPhotTheta()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


template <typename T> // typename
o2::math_utils::Vector2D<double> Recon::traceForward(Polar3DVector dirCkov) const
{
  // Trace forward a photon from (x,y) up to PC
  //  Arguments: dirCkov photon vector in LORS
  //    Returns: pos of traced photon at PC

  math_utils::Vector2D<double> pos(-999, -999);
  double thetaCer = dirCkov.Theta();
  if (thetaCer > TMath::ASin(1. / fParam->getRefIdx())) {
    return pos;                                                       // total refraction on WIN-GAP boundary
  }
  double zRad = -0.5 * fParam->radThick() - 0.5 * fParam->winThick(); // z position of middle of RAD

  math_utils::Vector3D<double> posCkov(fTrkPos.X(), fTrkPos.Y(), zRad);
									      // RAD: photon position is track position @ middle of RAD

  propagate(dirCkov, posCkov, -0.5 * fParam->winThick());                     // go to RAD-WIN boundary
  refract(dirCkov, fParam->getRefIdx(), fParam->winIdx());                    // RAD-WIN refraction
  propagate(dirCkov, posCkov, 0.5 * fParam->winThick());                      // go to WIN-GAP boundary
  refract(dirCkov, fParam->winIdx(), fParam->gapIdx());                       // WIN-GAP refraction
  propagate(dirCkov, posCkov, 0.5 * fParam->winThick() + fParam->gapThick()); // go to PC
  pos.SetCoordinates(posCkov.X(), posCkov.Y());
  return pos;
} // TraceForward()

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Recon::lors2Trs(Polar3DVector dirCkov, double& thetaCer, double& phiCer) const
{
  // Theta Cerenkov reconstruction
  //  Arguments: dirCkov photon vector in LORS
  //    Returns: thetaCer of photon in TRS
  //               phiCer of photon in TRS
  //  TVector3 dirTrk;
  //  dirTrk.SetMagThetaPhi(1,fTrkDir.Theta(),fTrkDir.Phi()); -> dirTrk.SetCoordinates(1,fTrkDir.Theta(),fTrkDir.Phi())
  //  double  thetaCer = TMath::ACos(dirCkov*dirTrk);
  
  ROOT::Math::Rotation3D mtheta(ROOT::Math::RotationY(-fTrkDir.Theta()));     // TRotation mtheta;  mtheta.RotateY(-fTrkDir.Theta());


  

  ROOT::Math::Rotation3D mphi(ROOT::Math::RotationZ(-fTrkDir.Phi()));         // mphi.RotateZ(-fTrkDir.Phi());

  ROOT::Math::Rotation3D mrot = mtheta * mphi;

  math_utils::Vector3D<double> dirCkovTRS;
  dirCkovTRS = mrot * dirCkov;
  phiCer = dirCkovTRS.Phi();     // actual value of the phi of the photon
  thetaCer = dirCkovTRS.Theta(); // actual value of thetaCerenkov of the photon
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Recon::trs2Lors(math_utils::Vector3D<double> dirCkov, double& thetaCer, double& phiCer) const
{
  // Theta Cerenkov reconstruction
  //  Arguments: dirCkov photon vector in TRS
  //    Returns: thetaCer of photon in LORS
  //               phiCer of photon in LORS



  ROOT::Math::Rotation3D mtheta(ROOT::Math::RotationY(fTrkDir.Theta()));      // TRotation mtheta; mtheta.RotateY(fTrkDir.Theta());

  ROOT::Math::Rotation3D mphi(ROOT::Math::RotationZ(fTrkDir.Phi()));

  ROOT::Math::Rotation3D mrot = mphi * mtheta;

  math_utils::Vector3D<double> dirCkovLORS = mrot * dirCkov;
  phiCer = dirCkovLORS.Phi();     // actual value of the phi of the photon
  thetaCer = dirCkovLORS.Theta(); // actual value of thetaCerenkov of the photon
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Recon::findRingGeom(double ckovAng, int level)
{
  // Find area covered in the PC acceptance
  // Arguments: ckovAng - cerenkov angle
  //            level   - precision in finding area and portion of ring accepted (multiple of 50)
  //   Returns: area of the ring in cm^2 for given theta ckov

  int kN = 50 * level;
  int nPoints = 0;
  double area = 0;

  bool first = kFALSE;

  math_utils::Vector2D<double> pos1;

  for (int i = 0; i < kN; i++) {
    if (!first) {
      pos1 = tracePhot(ckovAng, double(TMath::TwoPi() * (i + 1) / kN)); // find a good trace for the first photon
      if (pos1.X() == -999) {
        continue; // no area: open ring
      }
      if (!fParam->isInside(pos1.X(), pos1.Y(), 0)) {
        pos1 = intWithEdge(fMipPos, pos1); // find the very first intersection...
      } else {
        if (!Param::isInDead(1.0f, 1.0f)) {// ef: moved method from Param.cxx to h
          nPoints++;                       // photon is accepted if not in dead zone
        }
      }
      first = kTRUE;
      continue;
    }
    math_utils::Vector2D<double> pos2 = tracePhot(ckovAng, double(TMath::TwoPi() * (i + 1) / kN)); // trace the next photon
    if (pos2.X() == -999) {
      continue; // no area: open ring
    } 
    if (!fParam->isInside(pos2.X(), pos2.Y(), 0)) {
      pos2 = intWithEdge(fMipPos, pos2);
    } else {
      if (!Param::isInDead(pos2.X(), pos2.Y())) {
        nPoints++; // photon is accepted if not in dead zone
      }
    }

    area += TMath::Abs((pos1 - fMipPos).X() * (pos2 - fMipPos).Y() - (pos1 - fMipPos).Y() * (pos2 - fMipPos).X()); // add area of the triangle...
    pos1 = pos2;
  }
  //---  find area and length of the ring;
  fRingAcc = (double)nPoints / (double)kN;
  area *= 0.5;
  fRingArea = area;
} // FindRingGeom()

// template <typename T>
void Recon::propagate(const Polar3DVector& dir, math_utils::Vector3D<double>& pos, double z) const
{
  // Finds an intersection point between a line and XY plane shifted along Z.
  // Arguments:  dir,pos   - vector along the line and any point of the line
  //             z         - z coordinate of plain
  //   Returns:  none
  //   On exit:  pos is the position if this intesection if any
  static math_utils::Vector3D<double> nrm(0, 0, 1);
  math_utils::Vector3D<double> pnt(0, 0, z);

  math_utils::Vector3D<double> diff = pnt - pos;
  double sint = 0; //(nrm * diff) / (nrm * dir);
  pos += sint * dir;
} // Propagate()

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T> // typename
void Recon::refract(Polar3DVector& dir, double n1, double n2) const
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