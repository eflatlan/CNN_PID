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
#include "CkovTools.cpp"
#include "populate.cpp"
// sudo yum install hdf5-devel
#include <HMPIDBase/Param.h>

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



using namespace o2;
using namespace o2::hmpid;



class ParticleInfo : public TObject {
public:
    // Your data members here
    float momentum;
    float mass;
    float energy;
    float refractiveIndex;
    float ckov;
    float xRad;
    float yRad;
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





void populateRegions(Populate* populate, std::vector<std::pair<double, double>>& vecArr, TH2F* map, const double& eta, const double& l) {
	 
	 const int kN = vecArr.size();
	 for(int i = 0; i < vecArr.size(); i++){
		  const auto& value = populate->tracePhot(eta, Double_t(TMath::TwoPi()*(i+1)/kN), l);
		  if(value.X() > 0 && value.X() < 156.0 && value.Y() > 0 && value.Y() < 144) {
		    map->Fill(value.X(), value.Y());
		    vecArr[i] = std::make_pair(value.X(), value.Y());
		  }
	 }
 		
 }

void setStyleInd(TH2* th1f, float ratio = 1.2);


// TODO : add std-dev or similar here :?
double getThetaP(double momentum)
{
  double degThetaP;
  float pH, pL, degThetaH, degThetaL;
  if (momentum < 0.5){
		pH = 0.5; pL = 0.4; degThetaH = 42.5; degThetaL = 50;
  } else if (momentum >= 0.5 && momentum < 0.7){
		pH = 0.7; pL = 0.5; degThetaH = 27.5; degThetaL = 42.5;
  } else if (momentum >= 0.7 && momentum < 1){
		pH = 1; pL = 0.7; degThetaH = 22.5; degThetaL = 27.5;
  } else if (momentum >= 1 && momentum < 1.5){
		pH = 1.5; pL = 1; degThetaH = 15; degThetaL = 22.5;
  } else if (momentum >= 1.5 && momentum < 2.5){
		pH = 2.5; pL = 1.5; degThetaH = 10; degThetaL = 15;
  } else {
		pH = 5; pL = 2.5; degThetaH = 8; degThetaL = 10;
  }
  
  degThetaP = degThetaL + (degThetaH-degThetaL)/(pH-pL) * (momentum - pL);

  return degThetaP*3.1415/180;
}



void setStyleInd(TH1* th1f, float ratio = 1.2);



const float mass_Pion = 0.1396, mass_Kaon = 0.4937, mass_Proton = 0.938; // masses in
std::array<float, 3> masses = {mass_Pion, mass_Kaon, mass_Proton};
// get Ring-radius from Cherenkov-angle

float calcRingGeom(float ckovAng, int level);
float GetFreonIndexOfRefraction(float x);

//float GetQuartzIndexOfRefraction(float x);
//float BackgroundFunc(float *x, float *par);

const float fDTheta = 0.001;  // increment
float kThetaMax = 0.75;
int nChannels = (int)(kThetaMax / fDTheta + 0.5);


void setStyle();




TH1D* phots = new TH1D("Photon Candidates", "Photon Candidates;angle [rad]; counts/1 mrad", nChannels, 0, kThetaMax);
TH1D* photsw = new TH1D("Photon Weights", "Photon Weights;angle [rad]; counts/1 mrad", nChannels, 0, kThetaMax);
TH1D* resultw = new TH1D("Sum of Weights in Window at Bin", "Sum of Weights in Window at bin;angle [rad]; counts/1 mrad", nChannels, 0, kThetaMax);
TH1F *hTheta = new TH1F("Background","Background; angle [rad]; counts/1 mrad",750,0.,0.75); 

TH1F *hThetaCh = new TH1F("Cherenkov Photons","Cherenkov Photons; angle [rad]; counts/1 mrad",750,0.,0.75); 

float meanCherenkovAngle;


std::vector<float> photonCandidates;

float ckovTrackOut = 0;
float /*std::array<TH1D*, 3>*/ houghResponse(std::vector<float>& photonCandidates, float fWindowWidth);




// make event number i:
std::vector<std::pair<double, double>>  makeEvent(std::vector<Bin>& mapBins, float occupancy, RandomValues& randomValue, ParticleInfo& particle, int cntEvent);



const float defaultPhotonEnergy = 6.75; 
const float refIndexFreon = GetFreonIndexOfRefraction(defaultPhotonEnergy);
const float refIndexQuartz = 1.5875;//GetQuartzIndexOfRefraction(defaultPhotonEnergy);
const float  refIndexCH4 = 1.00; 
const float CH4GapWidth = 8;
const float  RadiatorWidth = 1.5; // ? was 1
const float  QuartzWindowWidth = 0.5;
const float  EmissionLenght = RadiatorWidth/2;


const float tGap = 8;
const float  rW = 1.5; // ? was 1
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


double scaleThetaP;

void testRandomMomentum(int numObjects = 10, float thetaTrackInclination = 0, double occupancy = 0.03, double _scaleThetaP = 1)
{  

  scaleThetaP = _scaleThetaP;


// create random mass, energy (and from this ref-index) and momentum :


  
  std::vector<RandomValues> randomObjects(numObjects);
  rndInt->SetSeed(0);

  std::vector<ParticleInfo> particleVector;

  int i = 0;

  int eventCnt = 0;
  for(auto& randomValue : randomObjects){
     // get cherenkov angle from mass momentum and refindex

    
     // get the map with a given occupancy and ckov angle calculated 
     std::vector<Bin> mapBins;


     ParticleInfo particle;

     Printf(" enter makeEvent %d", eventCnt); 
     const auto filledBins = makeEvent(mapBins, occupancy, randomValue, particle, eventCnt); 

     eventCnt++;
     
     Printf(" exit makeEvent %d", eventCnt); 
     Printf("testRandomMomentum():  mapBins.size = %zu", mapBins.size()); 
    
      // TODO : hSignalAndNoiseMap just placeholder, make instead a TH2F from the filledBins?

     /*
     const auto& map =  new TH2F("Signal and Noise2 ", "Signal and Noise2 ; x [cm]; y [cm]",160,0.,159.,144,0,143); */ 

     // make sure the momentum is valid for the given particle (e.g., no - in the square-root in calcCkovFromMass and acos [-1..1])
    if (particle.ckov == 0) {
      continue;
    }
    i++;

    
     // ParticleInfo Has to be extended to also contain xRad, yP, (xMIP, yMIP?), thetaP, phiP 

     particleVector.emplace_back(particle);

     Printf("CkovAngle %f Mass %f RefIndex %f Momentum %f | Num Entries in Map : %zu", particle.ckov, particle.mass, particle.refractiveIndex, particle.momentum, particle.filledBins.size()); 
     //map->SaveAs(Form("map%d.root", i));
  }

  // save object
  
  saveParticleInfoToROOT(particleVector);
}


//(mapBins, ckov, occupancy, thetaTrackInclination, photonEnergy, 



std::vector<std::pair<double, double>>  makeEvent(std::vector<Bin>& mapBins, float occupancy, RandomValues& randomValue, ParticleInfo& particle, int eventCnt)  
{


   // randomValue.momentum = 0.6;
   const auto& ckov = calcCkovFromMass(randomValue.momentum, randomValue.refractiveIndex, randomValue.mass); //  calcCkovFromMass(momentum, n, mass)
     
   // theoretical ckov angles : 

  // TODO: fjern denne ighen
  // auto momentum = 1.5; // 
  // const auto& ckovHyps2 = calcCherenkovHyp(momentum, randomValue.refractiveIndex); 
  

  // funker ikke segm på p = .5? / jo, kun at momentum settes til annen verdi enn rV.momentum
  const auto& ckovHyps = calcCherenkovHyp(randomValue.momentum, randomValue.refractiveIndex); 




  if(particle.energy < 5) particle.energy = 6.75; 
  auto ckovAngle = ckov;

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





  //float /*RadiatorWidth,*/ QuartzWindowWidth,CH4GapWidth,EmissionLenght;
  float FreonIndexOfRefraction,QuartzIndexOfRefraction,CH4IndexOfRefraction;


  // randomly created in RandomValues.cpp

  
  FreonIndexOfRefraction = GetFreonIndexOfRefraction(particle.energy);
  QuartzIndexOfRefraction = 1.5787;// GetQuartzIndexOfRefraction(particle.energy);
  CH4IndexOfRefraction = 1.0005;
  
  
  // TODO: this must be changed??
  float nF = FreonIndexOfRefraction;
  float nQ = QuartzIndexOfRefraction;
  float nG = CH4IndexOfRefraction;

  
  setStyle();

  //TH2F *hSignalAndNoiseMap = new TH2F("Signal and Noise ", "Signal and Noise ; x [cm]; y [cm]",1000,-25.,25.,1000,-25.,25.);



  float mapArray[40][40]{};
      

  // rndm value in ranfe 0.4, 0.7?
 TRandom2* rnd = new TRandom2(1);
 rnd->SetSeed(0);
 gRandom->SetSeed(0);
 // MIP azimuthal angle:

 // TODO: change phiP back again
 // setting phiP of track
 double phiP = TMath::Pi()/4;// static_cast<float>((3.14159)*(1-2*gRandom->Rndm(1)));
 phiP = 1*static_cast<float>((3.14159)*(1-2*gRandom->Rndm(1)));



 // setting thetaP of track
 // MIP polar angle:
 double thetaP = 1*scaleThetaP*getThetaP(randomValue.momentum);// = static_cast<float>(0.1*(1-2*gRandom->Rndm(1)));
  
         


 // place the impact point on rad : in x[10..150] and y[10..134]

  auto diff = 10;


  // TODO: change back
  // setting radiator impact point of track
  double xRad = 1*static_cast<float>((160-diff)*(1*gRandom->Rndm())+diff);
  double yRad = 1*static_cast<float>((144-diff)*(1*gRandom->Rndm())+diff);
 
  auto winThick = 0.5, radThick = 1.5; int gapThick = 8;
  auto getRefIdx = static_cast<double>(nF),  gapIdx = 1.0005, winIdx = 1.5787;


  // assumming L = 0.5*rW
  auto delta = (radThick*0.5 + winThick + gapThick)*TMath::Tan(thetaP);
  auto xPC =  xRad + delta*TMath::Cos(phiP);
  auto yPC =  yRad + delta*TMath::Sin(phiP);


  

  /*hSignalMIP->Fill(xRad, yRad);
  hSignalMIPpc->Fill(xPC, yPC);*/
 

 const auto lMax = 1.5, lMin = 0.0;



  // double L = 0;
  double L = static_cast<float>((rW)*(gRandom->Rndm(1)));
  
  // TODO; change back to:
  // double radParams[6] = {xRad,yRad,L,thetaP, phiP, randomValue.momentum};
  double radParams[7] = {xRad,yRad,L,thetaP, phiP, randomValue.momentum, randomValue.mass};
  double refIndexes[3] = {nF, nQ, nG};



 
 TVector2 trkPos(xRad, yRad);
 TVector3 trkDir; trkDir.SetMagThetaPhi(1, thetaP, phiP);
 
 // getR* getr = new getR(trkPos, trkDir, nF);
 


 Populate* populate = new Populate(trkPos, trkDir, nF);

 CkovTools ckovTools(radParams, refIndexes, ckovHyps, occupancy, ckovAngle, eventCnt);



 const TVector2& posMaxProton = populate->tracePhot(ckovTools.getMaxCkovProton(), 0, lMax);





  Printf("bgstudy segment : phiP %f thetaP %f xRad %f yP %f ", phiP, thetaP, xRad, yRad);


  Printf(" makeEvent%d : enter  numberOfCkovPhotons loop",eventCnt); 
  Printf(" makeEvent%d : numberOfCkovPhotons %d", eventCnt, numberOfCkovPhotons); 
 std::vector<std::array<double, 3>> cherenkovPhotons(numberOfCkovPhotons);

 //Populate populate()


 
 
 // maximum possible radius between points
 // TODO : change to be restricted by ckovHyps! (exceeding p-threshold)
 const auto rMax = (populate->getPcImp() - posMaxProton).Mod();
 Printf("posMaxProton --> rMax %.2f", rMax);


 //auto sLine = std::make_unique<std::unique_ptr<TLine>[]>(37); 
 std::vector<std::pair<double, double>> tmp;
 //std::vector<std::unique_ptr<TLine>> sLine;
 



 const int kN = 200;

 

 Printf(" makeEvent%d : ckovHyps = <%.3f, %.3f> | <%.3f, %.3f> | <%.3f, %.3f>", eventCnt, ckovTools.getMinCkovPion(),ckovTools.getMaxCkovPion(),ckovTools.getMinCkovKaon(),
ckovTools.getMaxCkovKaon(),ckovTools.getMinCkovProton(), ckovTools.getMaxCkovProton(), eventCnt); 

// void populateRegions(std::vector<std::array<double, 3>>& vecArr, TH2F* map, const double& eta, const double& l);
 
 /*
 std::vector<std::pair<double, double>> maxPionVec, maxKaonVec, maxProtonVec;
  
 std::vector<std::array<double, 3>> arrMaxPion, arrMaxKaon, arrMaxProton;
 
 maxPionVec.reserve(kN); maxKaonVec.reserve(kN); maxProtonVec.reserve(kN); 
 maxPionVec.resize(kN); maxKaonVec.resize(kN); maxProtonVec.resize(kN);
 Printf(" makeEvent%d : populating loop",eventCnt); 

  std::vector<std::pair<double, double>> minPionVec, minKaonVec, minProtonVec;
  minPionVec.reserve(kN); minKaonVec.reserve(kN); minProtonVec.reserve(kN); 
  minPionVec.resize(kN); minKaonVec.resize(kN); minProtonVec.resize(kN);
  Printf(" makeEvent%d : populating loop", eventCnt); 

  populateRegions(populate, maxPionVec, hMaxPion, ckovTools.getMaxCkovPion(), lMin);
  populateRegions(populate, maxKaonVec, hMaxKaon, ckovTools.getMaxCkovKaon(), lMin);
  populateRegions(populate, maxProtonVec, hMaxProton, ckovTools.getMaxCkovProton(), lMin);
  populateRegions(populate, minPionVec, hMinPion, ckovTools.getMinCkovPion(), lMax);
  populateRegions(populate, minKaonVec, hMinKaon, ckovTools.getMinCkovKaon(), lMax);
  populateRegions(populate, minProtonVec, hMinProton, ckovTools.getMinCkovProton(), lMax);
*/ 
 // track hit at PC in LORS : 
 const auto posPC = populate->getPcImp();


 Printf(" makeEvent%d : exit populating loop ", eventCnt); 


 
 // Printf("		sizes = kN %d | %zu %zu %zu", kN, maxPionVec.size(),    maxKaonVec.size(), maxProtonVec.size()); 

 auto deltaPcRad = (populate->getPcImp() - populate->getTrackPos());
 //else diff = vec.Phi() - phiL;
 Printf("deltaPhi (rad->pc) %.4f | phiP  %.4f ", deltaPcRad.Phi(), phiP);

 int photonCount = 0;
 for(photonCount=0; photonCount < numberOfCkovPhotons; photonCount++) {
   
   // TODO: endre std-dev her til å følge prob-dist?!
   // TODO: change back to 0.008
   float etaC = rnd->Gaus(ckovAngle, 0.008);		    // random CkovAngle, with 0.012 std-dev

   float phiL = static_cast<float>((TMath::Pi())*(1-2*gRandom->Rndm(1)));
   // angle around ckov Cone of photon



   TVector2 phot = populate->tracePhot(etaC, phiL, L); // tracePhot(double ckovThe, double ckovPhi)
   //Printf("made photon using populate->tracePhot | L = %.3f, etaC = %.4f, ckovAngle = %.4f | x %.3f, y %.3f ", L, etaC, ckovAngle, phot.X(), phot.Y());
	 
  //auto vec = phot - populate->getImpPc();     
   auto vec = phot - populate->getTrackPos();


   float diff;
   if(phiL <  0) diff = TMath::TwoPi() - (vec.Phi() - phiL);
   else diff = vec.Phi() - phiL;
   //Printf("phiL %.4f | phiRingAngle  %.4f | diff %.4f | phiP %.4f ", phiL, vec.Phi(), diff,  phiP);

   // populate map with cherenkov photon

   // later change below fcn to take valeus from setPhoton:
   //Printf(" makeEvent : enter  ckovTools.makeCkovPhoton"); 	 
   //const auto& ckovPhotonCoordinates = ckovTools.makeCkovPhoton(phiL, etaC);
   
   //std::array<double, 3> cand = {ckovPhotonCoordinates.first,ckovPhotonCoordinates.second , etaC};


   // photX photY are "ideal" values, should add some uncertainty to them...
   std::array<double, 3> cand = {phot.X(), phot.Y() , etaC};
  

   //cherenkovPhotons.emplace_back(cand);
   cherenkovPhotons[photonCount] = cand;

   //Printf(" makeEvent : ckovTools.makeCkovPhoton returned x %f y%f", ckovPhotonCoordinates.first + xRad, ckovPhotonCoordinates.second + yRad); 	    
   //hSignalAndNoiseMap->Fill(phot.X(), phot.Y());

  } 


   Printf(" makeEvent%d : Number of photons created :  %d ",eventCnt, photonCount); 

	
  // local coordinates are transformed to global here : 
  // also population of noise is done here

  //Printf(" makeEvent : num cherenkovPhotons %f", cherenkovPhotons.size()); 
  typedef std::vector<std::pair<double, double>> MapType;
  MapType temp;

 

  Printf(" makeEvent%d : enter  ckovTools.segment",eventCnt);
  const auto photonCandidatesCoords = ckovTools.segment(cherenkovPhotons, temp); // temp --> mapBins



  Printf(" makeEvent%d : num photonCandidatesCoords %zu",eventCnt, photonCandidatesCoords.size()); 	 
  for(const auto& photons : photonCandidatesCoords){
    //Printf("photon x %f y %f", photons.first, photons.second);
    //hSignalAndNoiseMap->Fill(photons.first, photons.second);
  }
    Printf("makeEvent%d() exit const auto& photons : photonCandidatesCoords", eventCnt);
  
  

  /*
  TCanvas *thSignalAndNoiseMap = new TCanvas("hSignalAndNoiseMap","hSignalAndNoiseMap",800,800);  
  thSignalAndNoiseMap->cd();
  hSignalAndNoiseMap->Draw();
  hSignalMIP->Draw("same");
  hSignalMIPpc->Draw("same"); */ 


  

  /*
  for(const auto& s : sLine) 
  { 
    // s->Draw();
  } */
  /*
  hMaxPion->Draw("same");     //  Printf("makeEvent()  hMaxPion->Draw");
  hMinPion->Draw("same");
  hMaxProton->Draw("same");   //Printf("makeEvent()  hMaxProton->Draw");
  hMinProton->Draw("same");
  hMinKaon->Draw("same");
  hMaxKaon->Draw("same"); */



/*
incL->SetMarkerColor(kRed);
incL2->SetMarkerColor(kRed);
incL2->Draw("same"); incL2->SetMarkerStyle(3);
incL->Draw("same");incL->SetMarkerStyle(2);
*/
  /*
	hMinPionMaxL->Draw("same");
	hMaxPionMinL->Draw("same");*/



  //l4->Draw("same");l3->Draw("same");
  //l2->Draw("same");//l1->Draw("same");
  //Printf("makeEvent()  l1->Draw()");
  //hSignalAndNoiseMap->Show();
  /*thSignalAndNoiseMap->SaveAs("thSignalAndNoiseMap.png");
  thSignalAndNoiseMap->Show();
  gPad->Update();*/ 
  int cnt = 0;
  for(auto& b : mapBins) {/*Printf("xval %f", b.x);*/ cnt++;}
  Printf(" mapBins size = %zu", mapBins.size());
  Printf(" mapBins cnt size = %d",   cnt);
 //Printf("Hough Window size = %f", Hwidth);
 
 /*auto ckovAnglePredicted = houghResponse(photonCandidates,  Hwidth); */

  particle.filledBins = mapBins;
  particle.momentum = randomValue.momentum;
  particle.mass = randomValue.mass;
  particle.energy = randomValue.energy;
  particle.refractiveIndex = randomValue.refractiveIndex;
  particle.ckov = ckov;
  particle.xRad = xRad;
  particle.yRad = yRad;
  particle.thetaP = thetaP;
  particle.phiP = phiP;
  //particle.map = hSignalAndNoiseMap;  


  Printf("end makeEvent%d", eventCnt);
  return photonCandidatesCoords;
 
} // 
//**********************************************************************************************************************************************************************************************************
float GetFreonIndexOfRefraction(float photonEnergy)
{
  auto x = photonEnergy;
  float k = 1.177 + (0.0172)*x;
  return k;
}
//**********************************************************************************************************************************************************************************************************


/*float GetQuartzIndexOfRefraction(float x)
{
  float k = TMath::Sqrt(1 + 46.411/(113.763556 - x) + 228.71/(328.51563 - x));
  return k;
}*/ 






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
  const float cos_ckov_denom = p*n;
  const auto cos_ckov_Pion = static_cast<float>(TMath::Sqrt(p_sq + mass_Pion_sq)/(cos_ckov_denom)); // n = refIndexFreon 1.289 later make it random?

  const float cos_ckov_Kaon = static_cast<float>(TMath::Sqrt(p_sq + mass_Kaon_sq)/(cos_ckov_denom)); 
  const float cos_ckov_Proton = static_cast<float>(TMath::Sqrt(p_sq + mass_Proton_sq)/(cos_ckov_denom));
  

  const float ckovAnglePion = static_cast<float>( TMath::ACos(cos_ckov_Pion)); 
  const float ckovAngleKaon = static_cast<float>(TMath::ACos(cos_ckov_Kaon)); 
  const float ckovAngleProton = static_cast<float>(TMath::ACos(cos_ckov_Proton)); 

  //Printf("Pion %.3f Kaon %.3f Proton %.3f", ckovAnglePion, ckovAngleKaon, ckovAngleProton);

  return {ckovAnglePion, ckovAngleKaon, ckovAngleProton};
}



float calcCkovFromMass(float p, float n, float m)
{
  const float p_sq = p*p;
  const float cos_ckov_denom = p*n;

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
        //TH2F* histCopy = new TH2F(*particle.map);
        //histCopy->SetName(histName);

	// Draw the map on the canvas with a 1:1 aspect ratio
	//canvas->SetCanvasSize(800, 800);
	//canvas->SetFixedAspectRatio();
        //histCopy->Write();
        //tree->Fill();
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



