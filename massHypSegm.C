#include <stdio.h>
#include <TROOT.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TPad.h>
#include <TH2F.h>
#include <Riostream.h>
#include <TRandom.h>
#include <TVector3.h>
#include <TRotation.h>
#include <TVector2.h>
#include <TH1D.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TF1.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <memory>

#include "SaveData.cpp"
#include "RandomValues.cpp"
#include "ParticleUtils.cpp"
// sudo yum install hdf5-devel


#include <Math/Vector3D.h>
#include <Math/Vector2D.h>
#include <Math/GenVector/Rotation3D.h>
#include <Math/GenVector/RotationX.h>
#include <Math/GenVector/RotationY.h>
#include <Math/GenVector/RotationZ.h>

#include "Math/Vector3D.h"
#include "Math/Vector2D.h"
#include "Math/GenVector/Rotation3D.h"
#include "Math/GenVector/RotationX.h"
#include "Math/GenVector/RotationY.h"
#include "Math/GenVector/RotationZ.h"



static constexpr float arrWaveLenDefault[30] = {
  162, 164, 166, 168, 170, 172, 174, 176, 178, 180,
  182, 184, 186, 188, 190, 192, 194, 196, 198, 200,
  202, 204, 206, 208, 210, 212, 214, 216, 218, 220};
static constexpr float nm2eV = 1239.842609;


class ParticleInfo : public TObject {
public:
    // Your data members here
    float momentum;
    float mass;
    float energy;
    float refractiveIndex;
    float ckov;
    float xP;
    float yP;
    float thetaP;
    float phiP;
    TH2F* map;
    std::vector<Bin> filledBins;
    ParticleInfo() {
        // Initialize your data members
    }

    virtual ~ParticleInfo() {}

    ClassDef(ParticleInfo,1) // For ROOT dictionary
};


/*
struct ParticleInfo {
    float momentum;
    float mass;
    float energy;
    float refractiveIndex;
    float ckov;
    TH2F* map;
}; */ 


float arrW[750]= {0.};


using namespace o2;
using namespace o2::hmpid;

#include <HMPIDBase/Param.h>
void setStyleInd(TH2* th1f, float ratio = 1.2);
void setStyleInd(TH1* th1f, float ratio = 1.2);



const float mass_Pion = 0.1396, mass_Kaon = 0.4937, mass_Proton = 0.938; // masses in
std::array<float, 3> masses = {mass_Pion, mass_Kaon, mass_Proton};
// get Ring-radius from Cherenkov-angle
float getRadiusFromCkov(float ckovAngle);
float calcRingGeom(float ckovAng, int level);
float GetFreonIndexOfRefraction(float x);
float GetQuartzIndexOfRefraction(float x);
float BackgroundFunc(float *x, float *par);

const float fDTheta = 0.001;  // increment
float kThetaMax = 0.75;
int nChannels = (int)(kThetaMax / fDTheta + 0.5);


void setStyle();


auto phots = new TH1D("Photon Candidates", "Photon Candidates;angle [rad]; counts/1 mrad", nChannels, 0, kThetaMax);
auto photsw = new TH1D("Photon Weights", "Photon Weights;angle [rad]; counts/1 mrad", nChannels, 0, kThetaMax);
auto resultw = new TH1D("Sum of Weights in Window at Bin", "Sum of Weights in Window at bin;angle [rad]; counts/1 mrad", nChannels, 0, kThetaMax);
TH1F *hTheta = new TH1F("Background","Background; angle [rad]; counts/1 mrad",750,0.,0.75); 

TH1F *hThetaCh = new TH1F("Cherenkov Photons","Cherenkov Photons; angle [rad]; counts/1 mrad",750,0.,0.75); 

float meanCherenkovAngle;


std::vector<float> photonCandidates;

float ckovTrackOut = 0;
float /*std::array<TH1D*, 3>*/ houghResponse(std::vector<float>& photonCandidates, float fWindowWidth);



const float defaultPhotonEnergy = 6.75; 
const float refIndexFreon = GetFreonIndexOfRefraction(defaultPhotonEnergy);
const float refIndexQuartz = GetQuartzIndexOfRefraction(defaultPhotonEnergy);
const float  refIndexCH4 = 1.00; 

double getCkovFromCoords(double xP, double yP, double x, double y, double phiP, double thetaP, float nF, float nQ, float nG);
std::vector<std::pair<double, double>>  backgroundStudy(std::vector<Bin>& mapBins, float ckovActual = 0.5, float occupancy = 0.01, float thetaTrack = 0, double photonEnergy = 6.75, float& xP, float& yP, float& thetaP, float& phiP, std::vector<double> ckovHyps);

const float CH4GapWidth = 8;
const float  RadiatorWidth = 1.;
const float  QuartzWindowWidth = 0.5;
const float  EmissionLenght = RadiatorWidth/2;


const float tGap = 8;
const float  rW = 1.;
const float  qW = 0.5;
const float  L = rW/2;


TH1* getMaxInRange(TH1* th1, float& up, float mid, float width);
float getMaxInRange(TH1* th1, int start, int width);
float getMaxInRange(TH1* th1, float mid, float width);


// mass_Pion_sq mass_Kaon_sq mass_Proton_sq GeV/c^2
const float mass_Pion_sq = mass_Pion*mass_Pion, mass_Kaon_sq = mass_Kaon*mass_Kaon,  mass_Proton_sq = mass_Proton*mass_Proton;


TRandom2* rndInt = new TRandom2(1); 

float calcCkovFromMass(float p, float n, float m);
std::array<float, 3> calcCherenkovHyp(float p, float n);


float randomMass(); float randomEnergy();

std::vector<Bin> fillMapVector(TH2F* map)
{
    std::vector<Bin> filledBins;

    int nBinsX = map->GetNbinsX();
    int nBinsY = map->GetNbinsY();

    for (int i = 1; i <= nBinsX; ++i) {
        for (int j = 1; j <= nBinsY; ++j) {
            float binContent = map->GetBinContent(i, j);
            if (binContent > 0) {
                float x = map->GetXaxis()->GetBinCenter(i);
                float y = map->GetYaxis()->GetBinCenter(j);
                filledBins.push_back(Bin{x, y});
            }
        }
    }

    return filledBins;
}


float randomMomentum()
{
  return 1+4*rndInt->Gaus(0.5, 0.25);
}


// create a number of random particles 
/*float[] randomParticles(int numParticles)
{
  DataSaver dataSaver(Form("RandomParticles%3d.root",numParticles));

  float[] randomMomentum = 

}*/

TH2F* tHistMass = new TH2F("test", "test; Momentum (GeV/c); Cherenkov Angle, #theta_{ch} (rad)", 5000, 0., 5., 800, 0., 0.8);


void testHyp()
{  
  TCanvas *tCkov = new TCanvas("ckov","ckov",800,800);  
//TH2F *hClusterMap = new TH2F("Cluster Map", "Cluster Map; x [cm]; y [cm]",1000,-10.,10.,1000,-10.,10.);
  rndInt->SetSeed(0);
  for(float p = 0.; p < 5; p+= 0.001)
  { 

    
    auto photonEnergy = randomEnergy();
    auto n = GetFreonIndexOfRefraction(photonEnergy); // refractive index
    Printf("P =  %f  || n = %f", p, n);
    auto ckovAngles = calcCherenkovHyp(p, n);
    for(auto& ckovAngle:ckovAngles){
      if(!TMath::IsNaN(ckovAngle)){tHistMass->Fill(p, ckovAngle);}
    }
  }
  tCkov->cd();
  tHistMass->Draw();
  //calcCherenkovHyp(1.5, 1.289);

}

std::array<float, 3> ckovMassHyp;
const float ckovConstraintWhidth = 0.015; // 15 mrad per side for ckovangle constraint 
					   // from mass hypotheses


float ckovActual = 0.;
float massActual = 0.;
float pActual = 0.;

void saveDataInst();
/*std::shared_ptr<TFile>*/void saveParticleInfoToROOT(const std::vector<ParticleInfo>& particleVector);




void testRandomMomentum(int numObjects = 10, float thetaTrackInclination = 0, double occupancy = 0.03)
{  


// create random mass, energy (and from this ref-index) and momentum :


  
  std::vector<RandomValues> randomObjects(numObjects);
  rndInt->SetSeed(0);

  std::vector<ParticleInfo> particleVector;

  int i = 0;
  for(const auto& randomValue : randomObjects){
     // get cherenkov angle from mass momentum and refindex
     const auto& ckov = calcCkovFromMass(randomValue.momentum, randomValue.refractiveIndex, randomValue.mass); //  calcCkovFromMass(momentum, n, mass)
     
     // theoretical ckov angles : 
     const auto& ckovHyps = calcCherenkovHyp(randomValue.momentum, randomValue.refractiveIndex); 
    
     // get the map with a given occupancy and ckov angle calculated 
     double photonEnergy = randomValue.energy;
     std::vector<Bin> mapBins;
    
     /* old : 
     const auto& map = backgroundStudy(mapBins, ckov, occupancy, thetaTrackInclination, photonEnergy); // cherenkov angle mean / occupancy / theta track inclination (perpendicular =)
     auto filledBins = fillMapVector(map); */
    
     // new : fill instead a vector with x, y pairs
     float xP, yP, thetaP, phiP; 
     const auto& filledBins = backgroundStudy(mapBins, ckov, occupancy, thetaTrackInclination, photonEnergy, xP, yP, thetaP, phiP, ckovHyps); // cherenkov angle mean / occupancy / theta track inclination (perpendicular =)
     
     

    
      // TODO : hSignalAndNoiseMap just placeholder, make instead a TH2F from the filledBins?
     const auto& map = hSignalAndNoiseMap;

     // make sure the momentum is valid for the given particle (e.g., no - in the square-root in calcCkovFromMass and acos [-1..1])
    if (ckov == 0) {
      continue;
    }
     i++;

    
     // ParticleInfo Has to be extended to also contain xP, yP, (xMIP, yMIP?), thetaP, phiP 
     ParticleInfo particle;
     particle.filledBins = mapBins;
     particle.momentum = randomValue.momentum;
     particle.mass = randomValue.mass;
     particle.energy = randomValue.energy;
     particle.refractiveIndex = randomValue.refractiveIndex;
     particle.ckov = ckov;
     particle.xP = xP;
     particle.yP = yP;
     particle.thetaP = thetaP;
     particle.phiP = phiP;

     particle.map = map;  
     particleVector.emplace_back(particle);

     Printf("CkovAngle %f Mass %f RefIndex %f Momentum %f | Num Entries in Map : %d", ckov, particle.mass, particle.refractiveIndex, particle.momentum, particle.filledBins.size()); 
     //map->SaveAs(Form("map%d.root", i));
  }

  // save object
  
  saveParticleInfoToROOT(particleVector);
}


//(mapBins, ckov, occupancy, thetaTrackInclination, photonEnergy, 

std::vector<std::pair<double, double>>  backgroundStudy(std::vector<Bin>& mapBins, float ckovActual = 0.5, float occupancy = 0.01, float thetaTrack = 0, double photonEnergy = 6.75, float& xP, float& yP, float& thetaP, float& phiP, std::vector<double> ckovHyps)   
{

  auto ckovAngle = ckovActual;

  Int_t NumberOfEvents = 1; Int_t NumberOfClusters = 13; float Hwidth = 15.;
  //testRandomMomentum();
  gStyle->SetOptStat("ei");

  TRandom2* rndP = new TRandom2(1); 
  rndP->SetSeed(0);


  // number of cherenkov photons in the cherenkov ring:
  const auto numberOfCkovPhotons = rndP->Poisson(13);

  photonCandidates.clear();
  float ThetaP = 0; // [rad]  // endre denne
  float PhiP=0,PhiF=0,DegPhiP=0;

  float DegThetaP = thetaTrack;

  //float /*RadiatorWidth,*/ QuartzWindowWidth,CH4GapWidth,EmissionLenght;
  float FreonIndexOfRefraction,QuartzIndexOfRefraction,CH4IndexOfRefraction;

  float PhotonEnergy = photonEnergy; 
  
  FreonIndexOfRefraction = GetFreonIndexOfRefraction(photonEnergy);
  QuartzIndexOfRefraction = GetQuartzIndexOfRefraction(photonEnergy);
  CH4IndexOfRefraction = 1.00;
  
  
  // TODO: this must be changed??
  float nF = FreonIndexOfRefraction;
  float nQ = QuartzIndexOfRefraction;
  float nG = CH4IndexOfRefraction;

  
  setStyle();

  //TH2F *hSignalAndNoiseMap = new TH2F("Signal and Noise ", "Signal and Noise ; x [cm]; y [cm]",1000,-25.,25.,1000,-25.,25.);

  TH2F *hSignalAndNoiseMap = new TH2F("Signal and Noise ", "Signal and Noise ; x [cm]; y [cm]",160,0.,159.,144,0,143);

  float mapArray[40][40]{};
      

  /* old code 
  for(Int_t iEvt = 0; iEvt<NumberOfEvents; iEvt++){
    
    //Printf("event number = %i",iEvt);
    
    float Xcen[100000],Ycen[100000];
     
    //DegThetaP = 4.;//0.*(1 - 2*gRandom->Rndm(iEvt));
    gRandom->SetSeed(0);


    DegPhiP   = 7.5;//360*gRandom->Rndm(iEvt); 
        
    // hvordan endre denne:?
    ThetaP = TMath::Pi()*DegThetaP/180;
    PhiP = TMath::Pi()*DegPhiP/180;  
     
    
    const int numBackGround = occupancy*6*80*48; //occupancy = 0.03
    NumberOfClusters = numBackGround;

    
    for(Int_t n1=0;n1<NumberOfClusters; n1++) {// clusters loop
      
    //  Printf("cluster = %i",n1);
      
      Xcen[n1] = 160*(gRandom->Rndm(n1));
      Ycen[n1] = 144*(gRandom->Rndm(n1));

      //noiseMap->Fill(Xcen[n1], Ycen[n1]);
      hSignalAndNoiseMap->Fill(Xcen[n1], Ycen[n1]);

      mapBins.push_back(Bin{Xcen[n1], Ycen[n1]});
      //mapArray[Xcen[n1]+20][Ycen[n1]+20] = 1;
      //hSignalAndNoiseMap2->Fill(Xcen[n1], Ycen[n1]);      
   } 
	  
  Int_t NphotTot = 0;
  float MeanTot = 0;
     
 }*/

   
 // MIP azimuthal angle:
 phiP = static_cast<float>((3.14159)*(1-2*gRandom->Rndm(1)));

 // MIP polar angle:
 // a.o.n; only between +- 5 deg
 thetaP = static_cast<float>((5/180*3.149265)*(1-2*gRandom->Rndm(1)));

 // rndm value in ranfe 0.4, 0.7?
 TRandom2* rnd = new TRandom2(1);
 rnd->SetSeed(0);

  const auto cosThetaP = TMath::Cos(thetaP);
  const auto sinThetaP = TMath::Sin(thetaP);
  const auto tanThetaP = TMath::Tan(thetaP);
  const auto cosPhiP = TMath::Cos(phiP);
  const auto sinPhiP = TMath::Sin(phiP);

 // place the impact point in x[10..150] and y[10..134]
 xP = static_cast<float>((160-10)*(1*gRandom->Rndm())-10);
 yP = static_cast<float>((144-10)*(1*gRandom->Rndm())-10);

 const auto xMipLocal = tanThetaP*cosPhiP;
 const auto yMipLocal = tanThetaP*sinPhiP;
	
 std::vector<std::pair<double, double>> cherenkovPhotons(numberOfCkovPhotons);
 for(Int_t i=0; i < numberOfCkovPhotons; i++) {
   
   // endre std-dev her til å følge prob-dist?!
   float ckovAnglePhot = rnd->Gaus(ckovAngle, 0.012);		    // random CkovAngle, with 0.012 std-dev
   float alpha = static_cast<float>((3.14159)*(1-2*gRandom->Rndm(1)));    // angle in photon-map (-pi to pi)
    
   // alpha, ckovAnglePhot will temp be placeholders..
   float phiL = alpha;
   float etaC = ckovAnglePhot;
	 
   // populate map with cherenkov photon	 
   const auto& ckovPhotonCoordinates = makeCkovPhoton(thetaP, phiP, phiL, etaC, nF, nG, nQ); // nF, nG, nQ to be added here? 

   cherenkovPhotons.push_back(ckovPhotonCoordinates);
	 
	 
   // old code, just a circle:
   // hThetaCh->Fill(ckovAngle)
   // trenger eliptisk skewing her:!
	  
	 
   /*	 
   float ringRadius = getRadiusFromCkov(ckovAnglePhot); // R in photon-map
   Printf("Cherenkov Photon : Angle = %f Radius = %f", ckovAnglePhot, ringRadius);	
   // populate the photon maps randomly radially

   // get x and y values of Photon-candidate:
   float x = xMip+ static_cast<float>(TMath::Cos(alpha)*ringRadius);
   float y = yMip+ static_cast<float>(TMath::Sin(alpha)*ringRadius);  

   // populating the pad
   hSignalAndNoiseMap->Fill(x,y);
   mapBins.push_back(Bin{x, y});
   //mapArray[Xcen[n1]+20][Ycen[n1]+20] = 1;

   // add Actual Cherenkov Photon to Candidates
   photonCandidates.emplace_back(ckovAnglePhot);*/ 

  } 
	
  // local coordinates are transformed to global here : 
  // also population of noise is done here
  

  const auto& photonCandidatesCoords = segment(xP, yP, xMipLocal, yMipLocal, cherenkovPhotons, ckovHyps, nF, nQ, nG);	
	
  int cnt = 0;
  for(auto& b : mapBins) {/*Printf("xval %f", b.x);*/ cnt++;}
  Printf(" mapBins size = %d", mapBins.size());
  Printf(" mapBins cnt size = %d",   cnt);
 //Printf("Hough Window size = %f", Hwidth);
 
 /*auto ckovAnglePredicted = houghResponse(photonCandidates,  Hwidth); */

 return photonCandidatesCoords;
 
}
//**********************************************************************************************************************************************************************************************************
float GetFreonIndexOfRefraction(float photonEnergy)
{
  auto x = photonEnergy;
  float k = 1.177 + (0.0172)*x;
  return k;
}
//**********************************************************************************************************************************************************************************************************
float GetQuartzIndexOfRefraction(float x)
  
{
  float k = TMath::Sqrt(1 + 46.411/(113.763556 - x) + 228.71/(328.51563 - x));
  return k;
}

//*********************************************************************************************************************************************************************************************************
float BackgroundFunc(float *x, float *par)
{
 float xx = x[0];
  
 float f = par[0]*TMath::Tan(TMath::ASin(1.2903*TMath::Sin(xx)))*(1+TMath::Tan(TMath::ASin(1.2903*TMath::Sin(xx)))*TMath::Tan(TMath::ASin(1.2903*TMath::Sin(xx)))*1.2903*TMath::Cos(xx)/cos(asin(1.2903*TMath::Sin(xx))));
  
 return f;
}       
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


// get Ring-radius from Cherenkov-angle
float getRadiusFromCkov(float ckovAngle)
{

  //// refIndexFreon refIndexQuartz refIndexCH4
  float sin_ckov = static_cast<float>(TMath::Sin(ckovAngle));
  float sin_qz = static_cast<float>(sin_ckov*(refIndexFreon/refIndexQuartz));
  float sin_theta0 = static_cast<float>(sin_qz*(refIndexQuartz/refIndexCH4));

  float R_ckov = sin_ckov*(RadiatorWidth - EmissionLenght);
  float R_qz = sin_qz * QuartzWindowWidth;
  float R_0 = sin_theta0*CH4GapWidth;
  //Printf("Radiuses  , R_ckov  % f + R_qz %f + R_0 %f", R_ckov, R_qz, R_0);
  float R = static_cast<float>(R_ckov + R_qz + R_0);
  return R;
} 


/*
*/

float /*std::array<TH1D*, 3>*/ houghResponse(std::vector<float>& photonCandidates, float fWindowWidth)
{

  // ckovMassHyp applies constraints!!
  
  auto l0 = ckovMassHyp[0] - ckovConstraintWhidth;
  auto u0 = ckovMassHyp[0] + ckovConstraintWhidth;
  auto l1 = ckovMassHyp[1] - ckovConstraintWhidth;
  auto u1 = ckovMassHyp[1] + ckovConstraintWhidth;
  auto l2 = ckovMassHyp[2] - ckovConstraintWhidth;
  auto u2 = ckovMassHyp[2] + ckovConstraintWhidth;
  Printf("Regions %f %f || %f %f || %f %f  ", l0, u0, l1, u1, l2, u2);

  
  
  /* ef : changed from this:
  // TH1D *resultw = new TH1D("resultw","resultw"       ,nChannels,0,kThetaMax);
  // TH1D *phots   = new TH1D("Rphot"  ,"phots"         ,nChannels,0,kThetaMax);
  // TH1D *photsw  = new TH1D("RphotWeighted" ,"photsw" ,nChannels,0,kThetaMax); */

  int nBin = (int)(kThetaMax / fDTheta);


  // ef : nCorrBand w the previous setting led to the bin1 = 1 bin2 = 750
  int nCorrBand = (int)(fWindowWidth / (2 /** fDTheta*/));

  Printf("nBin %d nCorrBand %d", nBin, nCorrBand);
  int binMax, sumMax = 0;
  std::vector<float> okAngles;
  okAngles.clear();
  for (const auto& angle : photonCandidates) { // photon cadidates loop

    if (angle < 0 || angle > kThetaMax)
      continue;

    // ef : check if angle in ckovMassHyp:
    if (angle < l2 || angle > u0)
      continue;


    phots->Fill(angle);

    int bin = (int)(0.5 + angle / (fDTheta));
    float weight = 1.;
    if (true) {
      float lowerlimit = ((float)bin) * fDTheta - 0.5 * fDTheta;
      float upperlimit = ((float)bin) * fDTheta + 0.5 * fDTheta;


      float rLow = getRadiusFromCkov(lowerlimit);
      float areaLow =  0.5*TMath::Pi()*TMath::Sq(rLow);// calcRingGeom(lowerlimit, 2);
 
      float rHigh = getRadiusFromCkov(upperlimit);
      float areaHigh =  0.5*TMath::Pi()*TMath::Sq(rHigh);// calcRingGeom(lowerlimit, 2);
      //Printf("Areas : areaLow %f, areaHigh %f ", areaLow, areaHigh);

      float diffArea = areaHigh - areaLow;
      
      if (diffArea > 0)
        weight = 1. / diffArea;
    }
    okAngles.emplace_back(angle);
    photsw->Fill(angle, weight);

    int nnn = static_cast<int>(angle*1000);
    arrW[nnn] += weight;
    //fPhotWei.emplace_back(weight); ef: do i need this?
  } // photon candidates loop

  for (int i = 1; i <= nBin; i++) {
    int bin1 = i - nCorrBand;
    int bin2 = i + nCorrBand;
    if (bin1 < 1)
      bin1 = 1;
    if (bin2 > nBin)
      bin2 = nBin;
    float sumPhots = phots->Integral(bin1, bin2);

    /*Printf("bin1 %d ; bin2 %d; sumPhots %f ", bin1, bin2, sumPhots);
    if (sumPhots < 3)
      continue; // if less then 3 photons don't trust to this ring*/
    float sumPhotsw = photsw->Integral(bin1, bin2);
    if ((float)((i /*+ 0.5*/) * fDTheta) > 0.7)
      continue;

    if (sumPhotsw > sumMax){
      binMax = i;
      sumMax = sumPhotsw;
      
    }
    resultw->Fill((float)((i /*+ 0.5*/) * fDTheta), sumPhotsw);
  }
  // evaluate the "BEST" theta ckov as the maximum value of histogramm

  float* pVec = (float*)resultw->GetArray();
  int locMax = TMath::LocMax(nBin, pVec);


  float smtest = 0; int ent=0;

  for(const auto& ok:okAngles){
    if (TMath::Abs(ok*1000-locMax) > nCorrBand)
      continue;     
    smtest+=ok; ent++;	
  }
  auto avgTest = smtest/ent;
  Printf("avgTest %f ent %d", avgTest, ent);
  



  Printf("pVec %f locMax %d", *pVec, locMax);
  Printf("sumMax %d, binMax %d", sumMax, binMax);
//photsw
  Printf("Entries : resultw %f photsw %f phots %f", resultw->GetEntries(), photsw->GetEntries(), phots->GetEntries());
  Printf("Max resultw %f photsw %f phots %f",  resultw->GetMaximum(10000.), photsw->GetMaximum(10000.), phots->GetMaximum(10000.));
  Printf("Min resultw %f photsw %f phots %f",  resultw->GetMinimum(-1.), photsw->GetMinimum(-1.), phots->GetMinimum(-1.));

  // ef: not this method, raw-pointers should not be used with new/delete-keywords
  //     smart-pointers are deleted when the fcuntion exits scope :
  // delete phots;delete photsw;delete resultw; // Reset and delete objects

  
  float ckovTrack = static_cast<float>(locMax * fDTheta + 0.5 * fDTheta); // final most probable track theta ckov
  ckovTrackOut = ckovTrack;


  float sumCkov = 0.;
  int entries = 0;
  float ckovSliding = phots->Integral(locMax-nCorrBand, locMax+nCorrBand);

  //binMax
  //locMax
  for(int i = locMax-nCorrBand; i < locMax+nCorrBand; i++)
  {

     auto val = phots->GetBinContent(i);
     auto w = photsw->GetBinContent(i);
     //Printf("Ckov Photons : Bin %d, Value %f; Weight %f", i, val,w);

     if(val > 0.){
       sumCkov += fDTheta*i;//+0.5*fDTheta;
       entries++;
     }
  }

  const float avgCkov = sumCkov/static_cast<float>(entries);
  ckovTrack = avgCkov;
  Printf("Avg Ckov = %f", avgCkov);

  Printf("sumCkov = %f || ckovSliding = %f", sumCkov,ckovSliding);

  TCanvas *houghCanvas = new TCanvas("Hough Canvas","Hough Canvas", 800, 800);
  houghCanvas->Divide(2,2);

  houghCanvas->cd(1);
  auto hThetaCl2 = static_cast<TH1*>(hTheta->Clone());
  setStyleInd(hThetaCl2);
  hThetaCl2->Draw();

  houghCanvas->cd(2);
  setStyleInd(phots);
  phots->Draw();
  TLatex lt2;
  lt2.DrawLatexNDC(.15, .85, Form("Window #eta_{c} :"));
  lt2.DrawLatexNDC(.16, .775, Form("Width = %.0f [mRad]", fWindowWidth));
  lt2.DrawLatexNDC(.16, .7, Form("Entries = %d", entries));

  auto up = getMaxInRange(phots, ckovTrack, fWindowWidth);
  TLine* l3 = new TLine(ckovTrack-(fWindowWidth/2)*0.001, 0, ckovTrack-(fWindowWidth/2)*0.001, up);
  TLine* l4 = new TLine(ckovTrack+(fWindowWidth/2)*0.001, 0, ckovTrack+(fWindowWidth/2)*0.001, up);
  TLine* l5 = new TLine(ckovTrack-(fWindowWidth/2)*0.001,up, ckovTrack+(fWindowWidth/2)*0.001, up);
  l3->SetLineColor(kGreen);  l4->SetLineColor(kGreen);l5->SetLineColor(kGreen);
  l3->SetLineWidth(2);    l4->SetLineWidth(2);l5->SetLineWidth(2);
  l3->Draw();  l4->Draw();l5->Draw();


  Printf("TBox Max %f ", photsw->GetBinContent(phots->GetMaximumBin()));
  Printf("TBox MaxBin %f ", phots->GetMaximumBin());


  int binBox = phots->GetYaxis()->GetLast();
  int yMax = phots->GetYaxis()->GetXmax();
  Printf("TBox binY %d ", binBox);
    Printf("TBox maxY %d ", yMax);

  TBox* box = new TBox(/*ckovTrack*/ckovTrack-nCorrBand*0.001,0.,/*ckovTrack*/ckovTrack+nCorrBand*0.001,phots->GetBinContent(ckovTrack*1000/*photsw->GetMaximumBin()*/));
  box->SetLineColor(kGreen);
  box->SetFillStyle(0);  
  box->SetLineWidth(2);  


  houghCanvas->cd(3);
  setStyleInd(photsw);
  float up_ = 0;
  auto th = getMaxInRange(photsw, up_, ckovTrack, fWindowWidth);
  th->Draw();




  auto l3c = static_cast<TLine*>(l3->Clone());  auto l4c = static_cast<TLine*>(l4->Clone());  auto l5c = static_cast<TLine*>(l5->Clone());
  l5c->SetY1(up_);
  l3c->SetY2(up_); l4c->SetY2(up_);l5c->SetY2(up_);
  l3c->Draw();  l4c->Draw();l5c->Draw();



  auto pad4 = static_cast<TPad*>(houghCanvas->cd(4));
  pad4->PaintText(.2, .9, Form("Track Cherenkov"));
  pad4->PaintText(.2, .8, Form("Actual Value = %.3f", meanCherenkovAngle));
  pad4->PaintText(.2, .7, Form("Predicted Value = %.3f", ckovTrack));
  pad4->PaintTextNDC(.2, .9, Form("Track Cherenkov"));
  pad4->PaintTextNDC(.2, .8, Form("Actual Value = %.3f", meanCherenkovAngle));
  pad4->PaintTextNDC(.2, .7, Form("Predicted Value = %.3f", ckovTrack));
  
  gStyle->SetPaintTextFormat("1.3f");
  setStyleInd(resultw);
  resultw->Draw();


  TLine* l = new TLine(/*ckovTrack*/ckovTrack, 0, /*ckovTrack*/ckovTrack, resultw->GetBinContent(locMax)*2);

  //TLine* l = new TLine(resultw->GetMaximumBin(), 0., resultw->GetMaximumBin(), resultw->GetBinContent(resultw->GetMaximumBin()));
  l->SetLineColor(kGreen);
  l->SetLineWidth(2);  
  l->Draw();

   
  Printf("ckovTrack = %f", ckovTrack);
  return ckovTrack;

}




void setStyleInd(TH1* th1f, float ratio = 1.2)
{
  th1f->SetTitleSize((th1f->GetTitleSize("x")*ratio), "xy");
  th1f->SetLabelSize((th1f->GetLabelSize("x")*ratio), "xy");
}


/*
void setPad(TPad*, float l, float, r, float t, float b)
{
  th1f->SetPadMaring((th1f->GetTitleSize("x")*ratio), "xy");
}*/


void setStyleInd(TH2* th1f, float ratio)
{
  th1f->SetTitleSize((th1f->GetTitleSize("x")*ratio), "xy");
  th1f->SetLabelSize((th1f->GetLabelSize("x")*ratio), "xy");
}





float getMaxInRange(TH1* th1, float mid, float width)
{


  //Printf("mid %f || width %f ", mid, width);
  float max = 1.0;

  const int startBin = static_cast<int>(mid*1000);
  const int nBin2 = static_cast<int>(width/2);
  //Printf("startBin %d || endBin %d ", startBin-nBin2, startBin+nBin2);


  int start = static_cast<int>(startBin-nBin2);
  int end = static_cast<int>(startBin+nBin2);
  for(int i = start; i < end; i++){
    auto binEnt = th1->GetBinContent(i);
    //if (binEnt > max) max = binEnt;
    //Printf("ent i %d || val %f ", i, binEnt);
  }
  return max;
}

TH1* getMaxInRange(TH1* th1, float& up, float mid, float width)
{
  TH1* thOut = static_cast<TH1*>(th1);

  //Printf("mid %f || width %f ", mid, width);
  float max = 1.0;

  const int startBin = static_cast<int>(mid*1000);
  const int nBin2 = static_cast<int>(width/2);
  //Printf("startBin %d || endBin %d ", startBin-nBin2, startBin+nBin2);


  int start = static_cast<int>(startBin-nBin2);
  int end = static_cast<int>(startBin+nBin2);

  float r = 0;
  for(const auto& i : arrW){
    if (i > r) 
      r = i;
  }
  
  for(int i = start; i < end; i++){
    auto binEnt = th1->GetBinContent(i);
    auto binent = arrW[i-1];
    thOut->SetBinContent(binent, i);
    //if (binEnt > max) max = binent;
    //Printf("ent i %d || val %f || val2 %f", i, binEnt, binent);

  }
  thOut->GetXaxis()->SetRangeUser(0., r);
  up = max;
  return thOut;
}

// mass_Pion_sq mass_Kaon_sq mass_Proton_sq
std::array<float, 3> calcCherenkovHyp(float p, float n)
{
  const float p_sq = p*p;
  const float cos_ckov_denom = p*refIndexFreon;
  const auto cos_ckov_Pion = static_cast<float>(TMath::Sqrt(p_sq + mass_Pion_sq)/(cos_ckov_denom)); // n = refIndexFreon 1.289 later make it random?

  const float cos_ckov_Kaon = static_cast<float>(TMath::Sqrt(p_sq + mass_Kaon_sq)/(cos_ckov_denom)); 
  const float cos_ckov_Proton = static_cast<float>(TMath::Sqrt(p_sq + mass_Proton_sq)/(cos_ckov_denom));
  

  const float ckovAnglePion = static_cast<float>( TMath::ACos(cos_ckov_Pion)); 
  const float ckovAngleKaon = static_cast<float>(TMath::ACos(cos_ckov_Kaon)); 
  const float ckovAngleProton = static_cast<float>(TMath::ACos(cos_ckov_Proton)); 

  Printf("Pion %.3f Kaon %.3f Proton %.3f", ckovAnglePion, ckovAngleKaon, ckovAngleProton);

  return {ckovAnglePion, ckovAngleKaon, ckovAngleProton};
}



float calcCkovFromMass(float p, float n, float m)
{
  const float p_sq = p*p;
  const float cos_ckov_denom = p*refIndexFreon;

  // sanity check ;)
  if(p_sq + m*m < 0){
    return 0;
  }

  const auto cos_ckov = static_cast<float>(TMath::Sqrt(p_sq + m*m)/(cos_ckov_denom));

  // sanity check ;)
  if(cos_ckov > 1 || cos_ckov < -1)
    return 0;

  const auto ckovAngle = static_cast<float>(TMath::ACos(cos_ckov));

  return ckovAngle;
}




float randomMass() 
{  
  auto index = static_cast<int>(rndInt->Integer(3));
  Printf("randomMass indes = %d", index);
  return masses[index];
}

float randomEnergy()
{
  // random energy from the arrWaveLenDefault
  auto index = static_cast<int>(rndInt->Integer(30));
  Printf("rEn indes = %d", index);
  float photonEnergy = static_cast<float>(nm2eV/arrWaveLenDefault[index]);
  return photonEnergy;
}

void setStyle()
{
  gStyle->SetOptStat("ei");
  phots->SetTitleSize(phots->GetTitleSize("x")*1.3, "xy");
  phots->SetLabelSize(phots->GetLabelSize("x")*1.3, "xy");

  photsw->SetTitleSize(phots->GetTitleSize("x"), "xy");
  photsw->SetLabelSize(phots->GetLabelSize("x"), "xy");

  resultw->SetTitleSize(phots->GetTitleSize("x"), "xy");
  resultw->SetLabelSize(phots->GetLabelSize("x"), "xy");

  hTheta->SetTitleSize(phots->GetTitleSize("x"), "xy");
  hTheta->SetLabelSize(phots->GetLabelSize("x"), "xy");
}

/*
void saveDataVector()
{
  std::vector<DataInfo> dataVector(100);  // Assuming you have a vector of Data

    // Fill the vector with some data
    // In a real case, you would likely fill this from your actual data source
    for (int i = 0; i < 100; i++) {
        dataVector[i].momentum = i * 0.5;  // some dummy values
        dataVector[i].typeOfParticle = i % 3;
        dataVector[i].refractiveIndex = i * 0.1;
        for(int j=0; j<10; j++){
            for(int k=0; k<10; k++){
                dataVector[i].pads[j][k] = i+j+k; // some dummy values
            }
        }
    }

    DataSaver dataSaver("data.root");
    dataSaver.fillData(dataVector);
    dataSaver.save();
}


void saveDataInst()
{

  //std::vector<int> dataVector(100);  // Assuming you have a vector of Data	
  std::vector<DataInfo> dataVector(100);  // Assuming you have a vector of Data

    // Fill the vector with some data
    // In a real case, you would likely fill this from your actual data source
    for (int i = 0; i < 100; i++) {
        dataVector[i].momentum = i * 0.5;  // some dummy values
        dataVector[i].typeOfParticle = i % 3;
        dataVector[i].refractiveIndex = i * 0.1;
        for(int j=0; j<10; j++){
            for(int k=0; k<10; k++){
                dataVector[i].pads[j][k] = i+j+k; // some dummy values
            }
        }
    }

    DataSaver dataSaver("data.root");
    dataSaver.fillData(dataVector);
    dataSaver.save();
} */ 


/*
// save TH2F* in own TTree since it caused segmentation fault when writing to same TTree as the other elements
std::shared_ptr<TFile> saveParticleInfoToROOT(const std::vector<ParticleInfo>& particleVector) {
    // Create a smart pointer for the TFile
    std::shared_ptr<TFile> outputFile(new TFile("outputFile.root", "RECREATE"));

    // Create a smart pointer for the TTree
    std::shared_ptr<TTree> tree(new TTree("tree", "ParticleInfo Tree"));

    // Create variables to hold the values of ParticleInfo properties
    float momentum;
    float mass;
    float energy;
    float refractiveIndex;
    float ckov;

    // Set the branch addresses for the TTree
    tree->Branch("momentum", &momentum);
    tree->Branch("mass", &mass);
    tree->Branch("energy", &energy);
    tree->Branch("refractiveIndex", &refractiveIndex);
    tree->Branch("ckov", &ckov);

    // Loop over the ParticleInfo objects and fill the TTree
    for (const auto& particle : particleVector) {
        momentum = particle.momentum;
        mass = particle.mass;
        energy = particle.energy;
        refractiveIndex = particle.refractiveIndex;
        ckov = particle.ckov;

        tree->Fill();

        // Ensure the map is not a nullptr before writing it to the file
        if (particle.map) {
            particle.map->Write(); // This will write the TH2F* to the file
        }
    }

    // Write the TTree to the TFile
    outputFile->cd();
    tree->Write();

    return outputFile;
}*/ 

void saveParticleInfoToROOT2(const std::vector<ParticleInfo>& particleVector) {
    // Create a new TFile
    TFile* outputFile = new TFile("outputFile2.root", "RECREATE");

    // Create a TTree
    TTree* tree = new TTree("tree", "ParticleInfo Tree");

    // Create a ParticleInfo pointer to hold the values
    ParticleInfo* particleInfo = nullptr;

    // Create a branch for ParticleInfo objects
    tree->Branch("ParticleInfo", "ParticleInfo", &particleInfo, 32000, 0);

    // Loop over the ParticleInfo objects and fill the TTree
    for (const auto& particle : particleVector) {
        // Delete the old object and create a new one with the current particle's values
        delete particleInfo;
        particleInfo = new ParticleInfo(particle);

        // Fill the tree
        tree->Fill();
    }

    // Write the TTree to the TFile
    tree->Write();

    // Close the TFile
    outputFile->Close();

    // Delete the last object created
    delete particleInfo;
}


void readParticleInfoFromROOT() {
    // Open the file
    TFile* inputFile = new TFile("outputFile.root", "READ");

    // Get the TTree
    TTree* tree = (TTree*)inputFile->Get("tree");

    // Variables to hold the values of ParticleInfo properties
    float momentum;
    float mass;
    float energy;
    float refractiveIndex;
    float ckov;

    // Set the branch addresses for the TTree
    tree->SetBranchAddress("momentum", &momentum);
    tree->SetBranchAddress("mass", &mass);
    tree->SetBranchAddress("energy", &energy);
    tree->SetBranchAddress("refractiveIndex", &refractiveIndex);
    tree->SetBranchAddress("ckov", &ckov);

    // Get number of entries in the TTree
    Long64_t nEntries = tree->GetEntries();

    // Loop over the entries in the TTree
    for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {
        // Load the entry
        tree->GetEntry(iEntry);

        // Print the ParticleInfo object's values
        std::cout << "ParticleInfo object " << iEntry << ":\n";
        std::cout << "  momentum: " << momentum << "\n";
        std::cout << "  mass: " << mass << "\n";
        std::cout << "  energy: " << energy << "\n";
        std::cout << "  refractiveIndex: " << refractiveIndex << "\n";
        std::cout << "  ckov: " << ckov << "\n";

        // Go to the "maps" directory
        TDirectory* mapsDir = (TDirectory*)inputFile->Get("maps");
        mapsDir->cd();

        // Get the associated map
        TString mapName = TString::Format("hist_%lld", iEntry);
        TH2F* map = (TH2F*)gDirectory->Get(mapName);

        if (map) {
            // Print the number of entries in the map
            std::cout << "  map entries: " << map->GetEntries() << "\n";
        } else {
            std::cout << "  map not found\n";
        }

        // Go back to the top directory
        inputFile->cd();
    }

    // Close the file
    inputFile->Close();
}

/*
void readParticleInfoFromROOT2() {
    // Open the file
    TFile* inputFile = new TFile("outputFile2.root", "READ");

    // Get the TTree
    TTree* tree = (TTree*)inputFile->Get("tree");

    // Create a ParticleInfo pointer
    ParticleInfo* particleInfo = nullptr;

    // Set the branch address
    Printf("Setting ParticleInfo adress");
    tree->SetBranchAddress("ParticleInfo", &particleInfo);

    // Get number of entries in the TTree
    Long64_t nEntries = tree->GetEntries();

    Printf("Entering Loop");

    // Loop over the entries in the TTree
    for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {
        // Load the entry


        if(particleInfo) {
            delete particleInfo;
            particleInfo = nullptr;
        }
        std::cout << "tree->GetEntry(iEntry); " << iEntry << ":\n";
        tree->GetEntry(iEntry);

        if(particleInfo) {
            std::cout << "ParticleInfo object " << iEntry << ":\n";
            std::cout << "  momentum: " << particleInfo->momentum << "\n";
            std::cout << "  mass: " << particleInfo->mass << "\n";
            std::cout << "  energy: " << particleInfo->energy << "\n";
            std::cout << "  refractiveIndex: " << particleInfo->refractiveIndex << "\n";
            std::cout << "  ckov: " << particleInfo->ckov << "\n";
        } else {
            std::cout << " ParticleInfo object not found\n";
        }

        // Get the associated histogram
        TString histName = TString::Format("hist_%lld", iEntry);
        TH2F* hist = (TH2F*)inputFile->Get(histName);
        if (hist) {
            // Print the number of points in the histogram
            std::cout << "  map points: " << hist->GetEntries() << "\n";
        } else {
            std::cout << "  map not found\n";
        }
    }

    // Close the file
    inputFile->Close();
}*/ 

void saveHD5(const std::vector<ParticleInfo>& particleVectorIn)
{
	// Convert TH2F* maps to 2D array and save ParticleInfo structs to an HDF5 file.
	std::vector<ParticleUtils::ParticleInfo> convertedVector;
	for (const auto& particle : particleVectorIn) {
	    ParticleUtils::ParticleInfo newParticle;
	    newParticle.momentum = particle.momentum;
	    newParticle.mass = particle.mass;
	    newParticle.energy = particle.energy;
	    newParticle.refractiveIndex = particle.refractiveIndex;
	    newParticle.ckov = particle.ckov;
	    newParticle.filledBins = particle.filledBins;

	    
	    // Assuming the map is a TH2F*, convert it to a 2D array
	    //newParticle.map = ParticleUtils::convertTH2FTo2DArray(particle.map);

	    // Add the new ParticleInfo with the converted map to the new vector
	    convertedVector.push_back(newParticle);

	    // Print the ParticleInfo object's values
            std::cout << "HD5 reading object "  << ":\n";
            std::cout << "  momentum: " << particle.momentum << "\n";
            std::cout << "  mass: " << particle.mass << "\n";
            std::cout << "  energy: " << particle.energy << "\n";
            std::cout << "  refractiveIndex: " << particle.refractiveIndex << "\n";
            std::cout << "  ckov: " << particle.ckov << "\n";

       
	}

	// Save the new vector of ParticleInfo with the converted maps to an HDF5 file
	ParticleUtils::saveParticleInfoToHDF5(convertedVector);



	// reading the entreis...
    	ParticleUtils::loadParticleInfoFromHDF5("ParticleInfo.h5");
	// Suppose you have a std::vector<ParticleInfo> named particleVector
	//std::vector<ParticleUtils::ParticleInfo> particleVector;

	// Populate particleVector...

	// Now save the vector of ParticleInfo to an HDF5 file
	//ParticleUtils::saveParticleInfoToHDF5(particleVector);
}

void saveParticleInfoToROOT(const std::vector<ParticleInfo>& particleVector) {
    // Create a new TFile
    TFile* outputFile = new TFile("outputFile.root", "RECREATE");

    // Create a TTree
    TTree* tree = new TTree("tree", "ParticleInfo Tree");

    // Create variables to hold the values of ParticleInfo properties
    float momentum;
    float mass;
    float energy;
    float refractiveIndex;
    float ckov;

    // Set the branch addresses for the TTree
    tree->Branch("momentum", &momentum, "momentum/F");
    tree->Branch("mass", &mass, "mass/F");
    tree->Branch("energy", &energy, "energy/F");
    tree->Branch("refractiveIndex", &refractiveIndex, "refractiveIndex/F");
    tree->Branch("ckov", &ckov, "ckov/F");

    // Create a directory to store the maps
    TDirectory* mapsDir = outputFile->mkdir("maps");
    mapsDir->cd();

    // Counter for unique histogram names
    int histCounter = 0;

    // Loop over the ParticleInfo objects and fill the TTree
    for (const auto& particle : particleVector) {
        momentum = particle.momentum;
        mass = particle.mass;
        energy = particle.energy;
        refractiveIndex = particle.refractiveIndex;
        ckov = particle.ckov;
	//TCanvas* canvas = new TCanvas("canvas", "Map Canvas", 800, 800);
        // Write each histogram to the maps directory with a unique name
        TString histName = TString::Format("hist_%d", histCounter++);
        TH2F* histCopy = new TH2F(*particle.map);
        histCopy->SetName(histName);

	// Draw the map on the canvas with a 1:1 aspect ratio
	//canvas->SetCanvasSize(800, 800);
	//canvas->SetFixedAspectRatio();
        histCopy->Write();
        tree->Fill();
    }

    // Go back to the top directory
    outputFile->cd();

    // Write the TTree to the TFile
    tree->Write();

    // Close the TFile
    outputFile->Close();


    // save by class: some problems
    //saveParticleInfoToROOT2(particleVector);

   
    saveHD5(particleVector);
    

    
    // works good
    Printf("\n\n Reading from file now...");
    readParticleInfoFromROOT();
}


// TODO : add nQ
// only consider photons in the correct range:
std::vector<std::pair<double, double>> segment(double yP, double yP, double xMipLocal, double yMipLocal, std::vector<std::pair<double, double>>& cherenkovPhotons, std::vector<double> ckovHyps, float nF, float nQ, float nG))
{ 
  // TODO: ckovHyps : get std-dev for Theta_ckov of pion kaon and proton from the values theta_i
  
  // placeholders for the above : 
  double stdDevPion = 0.08; double stdDevKaon = 0.08; double stdDevProton = 0.08;
  double ckovPionMin = ckovHyps[0] - 3 * stdDevPion;
  double ckovPionMax = ckovHyps[0] + 3 * stdDevPion;

  double ckovKaonMin = ckovHyps[1] - 3 * stdDevKaon;
  double ckovKaonMax = ckovHyps[1] + 3 * stdDevKaon;

  double ckovProtonMin = ckovHyps[2] - 3 * stdDevProton;
  double ckovProtonMax = ckovHyps[2] + 3 * stdDevProton;


  // these are in local coordinate system
  // create bbox of minimal possible ckov
  // Proton with eta_c = theta_c_proton - 3*std_dev_proton
  const auto& coordsMinPhi0 = populateMap(thetaP, phiP, 0, ckovProtonMin);
  xMinPhi0 = coordsMinPhi0.first;
  yMinPhi0 = coordsMinPhi0.first;

  // these are in local coordinate system
  const auto& coordsMinPhiPi = populateMap(thetaP, phiP, TMath::Pi, ckovProtonMin);
  xMinPhiPi = coordsMinPhiPi.first;
  yMinPhiPi = coordsMinPhiPi.first;

  const auto l1Min = TMath::Sqrt((xMinPhiPi-  xMipLocal)*(xMinPhiPi-xMipLocal) + (yMinPhiPi-yMipLocal)*(yMinPhiPi-yMipLocal));

  const auto l2Min = TMath::Sqrt((xMinPhi0-xMipLocal)*(xMinPhi0-xMipLocal) + (yMinPhi0-yMipLocal)*(yMinPhi0-yMipLocal));

 const auto rMin = getR(thetaP, phiP, etaCMin);


  // these are in local coordinate system
  // create bbox of maximal possible ckov
  // Proton with eta_c = theta_c_pion + 3*std_dev_pion

  const auto& coordsMinPhi0 = populateMap(thetaP, phiP, 0, ckovPionMax)
  xMinPhi0 = coordsMinPhi0.first;
  yMinPhi0 = coordsMinPhi0.first;

  // these are in local coordinate system
  const auto& coordsMinPhiPi = populateMap(thetaP, phiP, TMath::Pi, ckovPionMax)
  xMaxPhiPi = coordsMaxPhiPi.first;
  yMaxPhiPi = coordsMaxPhiPi.first;

  const auto l1Max = TMath::Sqrt((xMaxPhiPi-xMipLocal)*(xMaxPhiPi-xMipLocal) + (yMaxPhiPi-yMipLocal)*(yMaxPhiPi-yMipLocal));

  const auto l2Max = TMath::Sqrt((xMaxPhi0-xMipLocal)*(xMaxPhi0-xMipLocal) + (yMaxPhi0-yMipLocal)*(yMaxPhi0-yMipLocal));

  const auto rMax = getR(thetaP, phiP, etaCMax);


  // populate with background:
  const auto area = rMax*2*(l1Max+l2Max);
  const auto numBackgroundPhotons = area*occupancy; 
  std::vector<std::pair<double, double>> backGroundPhotons(numBackgroundPhotons);
  
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis1(0, l1Max + l2Max);
std::uniform_real_distribution<> dis2(-rMax, rMax);


for (auto& pair : backGroundPhotons) {
  pair.first = dis1(gen);
  pair.second = dis2(gen);
}

 std::vector<std::pair<double, double>> photonCandidates;

// here : loop through all backGroundPhotons and cherenkovPhotons

  for(const auto& photons : cherenkovPhotons)
  {  
    const auto& x = photons.first;
    const auto& y = photons.second;
    if(x > xMin && x < xMax && y > yMin && y < yMax){
      
      bool withinRange = false; 
      // check if the coordinates also corresponds to one of the possible cherenkov candidates
      const auto& ckov = getCkovFromCoords(xP, yP, x, y, phiP, thetaP, nF, nQ, nG);      
      
      // TODO: later, also add to candidates (i.e., pionCandidates, kaonCandidates...)
      if( ckovPionMin <  ckov & ckov < ckovPionMax ){
        withinRange = true;
      }
      
      if( ckovKaonMin <  ckov & ckov < ckovKaonMax ){
        withinRange = true;
      }
      
      if( ckovProtonMin <  ckov & ckov < ckovProtonMax ){
        withinRange = true;
      }
      
      if(withinRange){
        // transform to global coordinates:
        const auto& coords = local2Global(phi, xP, yP, x, y);

        filledBins.push_back(coords);
      }
        
        // Fill map
      //hSignalAndNoiseMap->Fill(Xcen[n1], Ycen[n1]);
      
    }    

  }

  for(const auto& photons : backGroundPhotons)
  {  
    const auto& x = photons.first;
    const auto& y = photons.second;
    if(x > xMin && x < xMax && y > yMin && y < yMax){

      bool withinRange = false; 
      // check if the coordinates also corresponds to one of the possible cherenkov candidates
      const auto& ckov = getCkovFromCoords(xP, yP, x, y, phiP, thetaP, nF, nQ, nG); 
           
      // TODO: later, also add to candidates (i.e., pionCandidates, kaonCandidates...)
      if( ckovPionMin <  ckov & ckov < ckovPionMax ){
        withinRange = true;
      }
      
      if( ckovKaonMin <  ckov & ckov < ckovKaonMax ){
        withinRange = true;
      }
      
      if( ckovProtonMin <  ckov & ckov < ckovProtonMax ){
        withinRange = true;
      }
      
      if(withinRange){
        // transform to global coordinates:
        const auto& coords = local2Global(phi, xP, yP, x, y);

        filledBins.push_back(coords);
      }      
    }    
  }
}

double getCkovFromCoords(double xP, double yP, double x, double y, double phiP, double thetaP, float nF, float nQ, float nG)
{    
  double phiF = 0;
  double ThetaF1,ThetaF2,ThetaF=0,ThetaLimite;
  float Xpi=0,Ypi=0,Xf=0,Yf=0,Xf1=0,Yf1=0,Xf2=0,Yf2=0; 
  float ThetaCherenkov, PhiCherenkov, DegPhiCherenkov;

  
  float deltaX = (rW+qW+tGap-L)*TMath::Tan(thetaP)*TMath::Cos(phiP);
  float deltaY = (rW+qW+tGap-L)*TMath::Tan(thetaP)*TMath::Sin(phiP);
	  
  Xpi = xP - deltaX;
  Ypi = yP - deltaY;

  TVector3 v2(Xcen[n1]-Xpi-EmissionLenght*TMath::Tan(ThetaP)*TMath::Cos(PhiP),Ycen[n1]-Ypi-EmissionLenght*TMath::Tan(ThetaP)*TMath::Sin(PhiP),RadiatorWidth+QuartzWindowWidth+CH4GapWidth-EmissionLenght); 

  PhiF = v2.Phi();      

  ThetaLimite = TMath::ASin(nG/nQ);

  double ThetaF0 = TMath::ASin(nQ/nF*TMath::Sin(ThetaLimite))-0.00001;

  //  Printf("ThetaF0 = %f",ThetaF0*TMath::RadToDeg());

  double ThetaF01 = TMath::ASin((nF/nQ)*(TMath::Sin(ThetaF0)));      

  double ThetaF02 = TMath::ASin((nQ/nG)*(TMath::Sin(ThetaF01)));

  float X01 = EmissionLenght*TMath::Tan(ThetaP)*TMath::Cos(PhiP);

  float Y01 =  EmissionLenght*TMath::Tan(ThetaP)*TMath::Sin(PhiP);

  float X02 = (RadiatorWidth - EmissionLenght)*TMath::Tan(ThetaF0)*TMath::Cos(PhiF)+QuartzWindowWidth*TMath::Tan(ThetaF01)*TMath::Cos(PhiF)+CH4GapWidth*TMath::Tan(ThetaF02)*TMath::Cos(PhiF);

  float Y02 = (RadiatorWidth - EmissionLenght)*TMath::Tan(ThetaF0)*TMath::Sin(PhiF) + QuartzWindowWidth*TMath::Tan(ThetaF01)*TMath::Sin(PhiF) + CH4GapWidth*TMath::Tan(ThetaF02)*TMath::Sin(PhiF);  

  float X0 = X01 + X02;
  float Y0 = Y01 + Y02;

  double ThetaMin = 0;
  //  Double_t ThetaMax = 0.75+ThetaP;
  double ThetaMax = ThetaF0;

  Xf = 999;
  Yf = 999;

  Int_t nWhile = 0;

  while(TMath::Sqrt((Xf-Xcen[n1]+Xpi)*(Xf-Xcen[n1]+Xpi)+(Yf-Ycen[n1]+Ypi)*(Yf-Ycen[n1]+Ypi))>0.0001)
  { 
    nWhile++;

    ThetaF = static_cast<double>( (0.5*(ThetaMax - ThetaMin) + ThetaMin));

    ThetaF1 = TMath::ASin((nF/nQ)*(TMath::Sin(ThetaF)));     
    ThetaF2 = TMath::ASin((nQ/nG)*(TMath::Sin(ThetaF1)));

    Xf1 = EmissionLenght*TMath::Tan(ThetaP)*TMath::Cos(PhiP);
    Yf1 =  EmissionLenght*TMath::Tan(ThetaP)*TMath::Sin(PhiP);

    Xf2 = (RadiatorWidth - EmissionLenght)*TMath::Tan(ThetaF)*TMath::Cos(PhiF)+QuartzWindowWidth*TMath::Tan(ThetaF1)*TMath::Cos(PhiF)+CH4GapWidth*TMath::Tan(ThetaF2)*TMath::Cos(PhiF);
    Yf2 = (RadiatorWidth - EmissionLenght)*TMath::Tan(ThetaF)*TMath::Sin(PhiF) + QuartzWindowWidth*TMath::Tan(ThetaF1)*TMath::Sin(PhiF) + CH4GapWidth*TMath::Tan(ThetaF2)*TMath::Sin(PhiF);  

    Xf = Xf1 + Xf2;
    Yf = Yf1 + Yf2;
    
    if(TMath::Sqrt((Xf-X0)*(Xf-X0)+(Yf-Y0)*(Yf-Y0))>TMath::Sqrt((Xcen[n1]-Xpi-X0)*(Xcen[n1]-Xpi-X0)+(Ycen[n1]-Ypi-Y0)*(Ycen[n1]-Ypi-Y0)))
    {
      ThetaMin = ThetaF;
    }

    else 
    {
      ThetaMax = ThetaF;   
    }  

    if(nWhile>30) break;

  } // while 
        
	      
  TVector3 vP((TMath::Sin(ThetaP))*(TMath::Cos(PhiP)),(TMath::Sin(ThetaP))*(TMath::Sin(PhiP)),(TMath::Cos(ThetaP)));
  TVector3 vz(0.,0.,1.);

  TVector3 v1 = vP.Cross(vz);

  TVector3 vF((TMath::Sin(ThetaF))*(TMath::Cos(PhiF)),(TMath::Sin(ThetaF))*(TMath::Sin(PhiF)),(TMath::Cos(ThetaF)));

  
  
  if(ThetaP==0)
  {	      
    ThetaCherenkov = ThetaF;		  
    PhiCherenkov = PhiF;	      
  }	  
  else		
  {
    vF.Rotate(ThetaP,v1);
    ThetaCherenkov = vF.Theta();
    PhiCherenkov = vF.Phi();
  }

  DegPhiCherenkov = 180*PhiCherenkov/(TMath::Pi());

  if(DegPhiCherenkov<0) DegPhiCherenkov+=360;	
}


std::pair<double, double> local2Global(double phi, double xP, double yP, double xL, double yL)
{
  const auto cosPhiP = TMath::Cos(phiP);
  const auto sinPhiP = TMath::Sin(phiP);
  const auto x = cosPhiP*xL + sinPhiP*yL - xP;
  const auto y = sinPhiP*xL - cosPhiP*yL - yP;
  return std::make_pair(x, y);
}

std::pair<double, double> global2Local(double phi, double xP, double yP, double xG, double yG)
{
  const auto cosPhiP = TMath::Cos(phiP);
  const auto sinPhiP = TMath::Sin(phiP);
  const auto x = cosPhiP*xG - sinPhiP*yG + xP;
  const auto y = sinPhiP*xG + cosPhiP*yG + yP;
  return std::make_pair(x, y);
}


// populate map in local reference system
// based on one etaC value (i.e., one photon at a time)
std::pair<double, double> populateMap(double thetaP, double phiP, double phiL, double etaC)
{
  const auto cosThetaP = TMath::Cos(thetaP);
  const auto sinThetaP = TMath::Sin(thetaP);
  const auto tanThetaP = TMath::Tan(thetaP);

  const auto cosPhiP = TMath::Cos(phiP);
  const auto sinPhiP = TMath::Sin(phiP);

  const auto cosPhiL = TMath::Cos(phiL);
  const auto sinPhiL = TMath::Sin(phiL);
  
  const auto cosEtaC = TMath::Cos(etaC);
  const auto sinEtaC = TMath::Sin(etaC);


  const auto rwlGap = (rW - L)/(TMath::Sqrt(1-sinEtaC*sinEtaC));

  const auto qwGap = (qW*nF)/(TMath::Sqrt(nQ*nQ-sinEtaC*sinEtaC*nF*nF));

  numerator = (tGap + tanThetaP*cosPhiL*sinPhiP * (rwlGap + qwGap));

  const auto denominator = 1- (tanThetaP*cosPhiL*sinPhiP*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));

  const auto tGapZ = numerator/denominator;


  const auto Lz = (rW-L) + qW + tGapZ;


  const auto tGapGap = (tGapZ*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF);

  const auto T = sinEtaC*(rwlGap+qwGap+tGapGap);

  const auto x = T*(cosPhiP*cosPhiL - sinPhiP*sinPhiL/cosThetaP) + Lz*tanThetaP*cosPhiP;

  const auto y = T*(sinPhiP*cosPhiL - cosPhiP*sinPhiL/cosThetaP) + Lz*tanThetaP*sinPhiP;

  return std::make_pair(x, y);
}




// get R at phiLocal = pi/2 V = 3pi/2
double getR(double thetaP, double phiP, double etaC)
{
  const auto cosThetaP = TMath::Cos(thetaP);
  const auto sinThetaP = TMath::Sin(thetaP);
  const auto tanThetaP = TMath::Tan(thetaP);

  const auto cosPhiP = TMath::Cos(phiP);
  const auto sinPhiP = TMath::Sin(phiP);

  const auto cosPhiL = 0; // phiLocal = pi/2
  const auto sinPhiL = 1; // --||-- 
  
  const auto cosEtaC = TMath::Cos(etaC);
  const auto sinEtaC = TMath::Sin(etaC);

  const auto rwlGap = (rW - L)/(TMath::Sqrt(1-sinEtaC*sinEtaC));

  const auto qwGap = (qW*nF)/(TMath::Sqrt(nQ*nQ-sinEtaC*sinEtaC*nF*nF));

  numerator = (tGap + tanThetaP*cosPhiL*sinPhiP * (rwlGap + qwGap));

  const auto denominator = 1- (tanThetaP*cosPhiL*sinPhiP*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));

  const auto tGapZ = numerator/denominator;


  const auto Lz = (rW-L) + qW + tGapZ;


  const auto tGapGap = (tGapZ*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF);

  const auto R = sinEtaC*(rwlGap+qwGap+tGapGap)/cos_theta_P;

  return R;
}


