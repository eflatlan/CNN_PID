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





float arrW[750]= {0.};


//using namespace o2;
//using namespace o2::hmpid;

#include <HMPIDBase/Param.h>




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

     
void backgroundStudy(std::vector<Bin>& mapBins, std::pair<float, float>& mipPos, float ckovActual = 0.5, float occupancy = 0, float thetaTrack = 0);

const float CH4GapWidth = 8;
const float  RadiatorWidth = 1.;
const float  QuartzWindowWidth = 0.5;
const float  EmissionLenght = RadiatorWidth/2;

struct ParticleInfo {
	float momentum;
	float mass;
	float energy;
	float refractiveIndex;
	float ckov;
	std::vector<Bin> filledBins;
	    std::pair<float, float> mipPos;
	//TH2F* map;
};


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




TH2F* tHistMass = new TH2F("test", "test; Momentum (GeV/c); Cherenkov Angle, #theta_{ch} (rad)", 5000, 0., 5., 800, 0., 0.8);


void testRingRadius()
{
  // Creating a canvas for the histograms
  TCanvas *tRingRadius = new TCanvas("ringRadius", "ringRadius", 1200, 800);
  tRingRadius->Divide(2, 2);

  gStyle->SetOptStat("ei");
  // Creating TH2F histograms for each particle class
  TH2F *hRingRadiusPion = new TH2F("Ring Radii - Pion", "Ring Radii - Pion; Momentum; Ring Radius (cm)",
                                   5000, 0., 5., 800, 0., 8);
  TH2F *hRingRadiusKaon = new TH2F("Ring Radii - Kaon", "Ring Radii - Kaon; Momentum; Ring Radius (cm)",
                                    5000, 0., 5.,800, 0., 8);
  TH2F *hRingRadiusProton = new TH2F("Ring Radii - Proton", "Ring Radii - Proton; Momentum; Ring Radius (cm)",
                                     5000, 0., 5.,800, 0., 8);
  TH2F *hRingRadiusCombined = new TH2F("Ring Radii - Combined", "Ring Radii - Combined; Momentum; Ring Radius (cm)",
                                     5000, 0., 5., 800, 0., 8);

  rndInt->SetSeed(0);

  for(float p = 0.; p < 5; p+= 0.01){
    // Calculate refractive index for random photon energy
    auto photonEnergy = randomEnergy();
    const auto n = 1.28;

    // Calculate Cherenkov angle for each mass
    std::array<float, 3> ckovAngles = calcCherenkovHyp(p, n);

    int cntx = 0;
    for (auto& ckovAngle : ckovAngles) {
      // Ignore NaN Cherenkov angles
      if (TMath::IsNaN(ckovAngle)) {
        Printf("fill values IsNaN ckov %f ", ckovAngle);
        continue;
      }

      // Get ring radius for current Cherenkov angle
      float ringRadius = getRadiusFromCkov(ckovAngle);

      // Add ring radius to respective histograms based on particle class
      if (cntx == 0) {
        hRingRadiusPion->Fill(p, ringRadius);
        hRingRadiusCombined->Fill(p, ringRadius);
      } else if (cntx == 1) {
        hRingRadiusKaon->Fill(p, ringRadius);
        hRingRadiusCombined->Fill(p, ringRadius);
      } else if (cntx == 2) {
        hRingRadiusProton->Fill(p, ringRadius);
        hRingRadiusCombined->Fill(p, ringRadius);
      }

      Printf("fill values : ckov %f  ringRadius %f | n = %d", ckovAngle, ringRadius, cntx);
      cntx++;
    }
  }

  // Set the statistics options for each histogram
  hRingRadiusPion->SetStats(kTRUE);
  hRingRadiusKaon->SetStats(kTRUE);
  hRingRadiusProton->SetStats(kTRUE);
  hRingRadiusCombined->SetStats(kTRUE);

  // Draw histograms on the canvas
  tRingRadius->cd(1);
  hRingRadiusPion->Draw();
  tRingRadius->cd(2);
  hRingRadiusKaon->Draw();
  tRingRadius->cd(3);
  hRingRadiusProton->Draw();
  tRingRadius->cd(4);
  hRingRadiusCombined->Draw();
  gStyle->SetOptStat("ei");
}




void testHyp()
{  
  TCanvas *tCkov = new TCanvas("ckov","ckov",800,800);  
//TH2F *hClusterMap = new TH2F("Cluster Map", "Cluster Map; x [cm]; y [cm]",1000,-10.,10.,1000,-10.,10.);
  rndInt->SetSeed(0);
  for(float p = 0.; p < 5; p+= 0.001)
  { 

    
    auto photonEnergy = randomEnergy();
    auto n = 1.28;//GetFreonIndexOfRefraction(photonEnergy); // refractive index
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


void saveHD5(const std::vector<ParticleInfo>& particleVectorIn);

void testRandomMomentum(int numObjects = 10, float thetaTrackInclination = 0)
{  


// create random mass, energy (and from this ref-index) and momentum :


  
  std::vector<RandomValues> randomObjects(numObjects);
  rndInt->SetSeed(0);

  std::vector<ParticleInfo> particleVector;

  int i = 0;
  for(const auto& randomValue : randomObjects){



     // get cherenkov angle from mass momentum and refindex
     const auto& ckov = calcCkovFromMass(randomValue.momentum, randomValue.refractiveIndex, randomValue.mass); //  calcCkovFromMass(momentum, n, mass)

     // get the map with a given occupancy and ckov angle calculated 

     std::vector<Bin> mapBins;
     std::pair<float, float> mipPos;


     if (ckov < 0.1) {
       continue;
     }

     backgroundStudy(mapBins, mipPos, ckov, 0.001, thetaTrackInclination); // cherenkov angle mean / occupancy / theta track inclination (perpendicular =)

     //std::cout << "study ret&:particleInfo.mipPos " << mipPos.first << " " << mipPos.second << std::endl; 

     //auto filledBins = fillMapVector(map);

  

     // make sure the momentum is valid for the given particle (e.g., no - in the square-root in calcCkovFromMass and acos [-1..1])

     i++;

     ParticleInfo particle;
     particle.filledBins = mapBins;
     particle.momentum = randomValue.momentum;
     particle.mass = randomValue.mass;
     particle.energy = randomValue.energy;
     particle.refractiveIndex = randomValue.refractiveIndex;
     particle.ckov = ckov;
     //particle.map = map; 
     particle.mipPos = mipPos;
     particleVector.emplace_back(particle);

     Printf("CkovAngle %f Mass %f RefIndex %f Momentum %f | Num Entries in Map : %d", ckov, particle.mass, particle.refractiveIndex, particle.momentum, particle.filledBins.size()); 
     //map->SaveAs(Form("map%d.root", i));
  }

  saveHD5(particleVector);
  //saveParticleInfoToROOT(particleVector);
}




void backgroundStudy(std::vector<Bin>& mapBins, std::pair<float, float>& mipPos, float ckovActual = 0.5, float occupancy = 0, float thetaTrack = 0)   
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
  float ThetaF1,ThetaF2,ThetaF=0,ThetaLimite;
  float Xpi=0,Ypi=0,Xf=0,Yf=0,Xf1=0,Yf1=0,Xf2=0,Yf2=0,Xp=0,Yp=0; 

  float PhotonEnergy = 6.75; 
  
  FreonIndexOfRefraction = GetFreonIndexOfRefraction(PhotonEnergy);
  QuartzIndexOfRefraction = GetQuartzIndexOfRefraction(PhotonEnergy);
  CH4IndexOfRefraction = 1.00;

  setStyle();

  //TH2F *hSignalAndNoiseMap = new TH2F("Signal and Noise ", "Signal and Noise ; x [cm]; y [cm]",1000,-25.,25.,1000,-25.,25.);

  //TH2F *hSignalAndNoiseMap = new TH2F("Signal and Noise ", "Signal and Noise ; x [cm]; y [cm]",160*10,0.,159.,144*10,0,143);

  //float mapArray[40][40]{};
  

  float Deltax = (RadiatorWidth+QuartzWindowWidth+CH4GapWidth-EmissionLenght)*TMath::Tan(ThetaP)*TMath::Cos(PhiP);
  float Deltay = (RadiatorWidth+QuartzWindowWidth+CH4GapWidth-EmissionLenght)*TMath::Tan(ThetaP)*TMath::Sin(PhiP);
	  
  Xpi = Xp - Deltax;
  Ypi = Yp - Deltay;
	  
  float ThetaCherenkov[100] = {0x0}, PhiCherenkov[100] = {0x0}, DegPhiCherenkov[100] = {0x0};
    gRandom->SetSeed(0);

  
  for(Int_t iEvt = 0; iEvt<NumberOfEvents; iEvt++){
    
    //Printf("event number = %i",iEvt);
    
    float Xcen[100],Ycen[100];
     
    //DegThetaP = 4.;//0.*(1 - 2*gRandom->Rndm(iEvt));



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
      //hSignalAndNoiseMap->Fill(Xcen[n1], Ycen[n1]);

      mapBins.push_back(Bin{Xcen[n1], Ycen[n1]});
      //mapArray[Xcen[n1]+20][Ycen[n1]+20] = 1;
      //hSignalAndNoiseMap2->Fill(Xcen[n1], Ycen[n1]);
      
   }
	  
  Int_t NphotTot = 0;
  float MeanTot = 0;
     
 }


 // rndm value in ranfe 0.4, 0.7?
 TRandom2* rnd = new TRandom2(1);
 rnd->SetSeed(0);



 // place the MIP in x[10..150] and y[10..134]
 const float xMip = static_cast<float>((160-20)*(1*gRandom->Rndm(1))+10);
 const float yMip = static_cast<float>((144-20)*(1*gRandom->Rndm(1))+10);


 mipPos = {xMip, yMip};



 Printf("Cherenkov Theta = %f numberOfCkovPhotons= %d", ckovAngle, numberOfCkovPhotons);	
 for(Int_t i=0; i < numberOfCkovPhotons; i++) {
   
   // endre std-dev her til å følge prob-dist?!
   float ckovAnglePhot = rnd->Gaus(ckovAngle, 0.02);		    // random CkovAngle, with 0.012 std-dev
   


   // hThetaCh->Fill(ckovAngle);


   // trenger eliptisk skewing her:!
   float ringRadius = getRadiusFromCkov(ckovAnglePhot); // R in photon-map


   if(ringRadius < 0 ) {
   } 
   Printf("Cherenkov Photon : Angle = %f Radius = %f", ckovAnglePhot, ringRadius);	



   // populate the photon maps randomly radially
   float alpha = static_cast<float>((3.14159)*(1-2*gRandom->Rndm(1)));    // angle in photon-map (-pi to pi)



   // get x and y values of Photon-candidate:
   float x = xMip+ static_cast<float>(TMath::Cos(alpha)*ringRadius);
   float y = yMip+ static_cast<float>(TMath::Sin(alpha)*ringRadius);  



   // populating the pad
   //hSignalAndNoiseMap->Fill(x,y);
   mapBins.push_back(Bin{x, y});
   //mapArray[Xcen[n1]+20][Ycen[n1]+20] = 1;

   // add Actual Cherenkov Photon to Candidates
   photonCandidates.emplace_back(ckovAnglePhot);

  } 
  int cnt = 0;
  for(auto& b : mapBins) {/*Printf("xval %f", b.x);*/ cnt++;}
  Printf(" mapBins size = %d", mapBins.size());
  Printf(" mapBins cnt size = %d",   cnt);
 //Printf("Hough Window size = %f", Hwidth);
 
 /*auto ckovAnglePredicted = houghResponse(photonCandidates,  Hwidth); */

 //return hSignalAndNoiseMap;
 
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
  const double p_sq = p*p;
  const double cos_ckov_denom = p*n;

  // sanity check ;)
  if(p_sq + m*m < 0){
    return 0;
  }

  const auto cos_ckov = static_cast<double>(TMath::Sqrt(p_sq + m*m)/(cos_ckov_denom));

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
	    newParticle.refractiveIndex = particle.refractiveIndex;
	    newParticle.ckov = particle.ckov;
	    newParticle.filledBins = particle.filledBins;
	    newParticle.mipPos = particle.mipPos;
	    
	    // Assuming the map is a TH2F*, convert it to a 2D array
	    //newParticle.map = ParticleUtils::convertTH2FTo2DArray(particle.map);

	    // Add the new ParticleInfo with the converted map to the new vector
	    convertedVector.push_back(newParticle);

	   
       
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


    /*
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
    //tree->Write();

    // Close the TFile
    outputFile->Close();


    // save by class: some problems
    //saveParticleInfoToROOT2(particleVector);

    */   
    saveHD5(particleVector);
    

    
    // works good
    Printf("\n\n Reading from file now...");
    //readParticleInfoFromROOT();
}




