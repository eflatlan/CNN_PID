#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <random>
#include "populate.cpp"

#include "populate2.cpp" // TODO: change name of class and file here 

#include <math.h>
//#include "ReconE.cpp"
#include "ReconG.cpp"

//namespace ParticleUtils
using namespace o2;




class ArrAndMap
{

  public :

  int scale = 10;

	TH2F* ckovCandMapRange = nullptr;

	TH2F* bgCandMapRange  = nullptr;

	TH2F* ckovCandMapOutRange = nullptr;

	TH2F* bgCandMapOutRange = nullptr;


  TH2F *hSignalAndNoiseMap =  nullptr;

		TH2F *hSignalMIP =  nullptr;
		TH2F *hSignalMIPpc =  nullptr;

		TH2F* arrayMaxMin[8];
		TH2F *hMaxProton =  nullptr;
		TH2F *hMaxPion =  nullptr;
		TH2F *hMaxPionMinL =  nullptr;
		TH2F *hMinPionMaxL =  nullptr;
		TH2F *hMaxKaon =  nullptr;


		TH2F *hMinProton =  nullptr;
		TH2F *hMinPion =  nullptr;
		TH2F *hMinKaon = nullptr;

    Populate* populatePtr = nullptr;


		int eventCnt;
  ArrAndMap(int eventCount, Populate* populatePtr) : eventCnt(eventCount),  populatePtr(populatePtr) {
		ckovCandMapRange = new TH2F("ckovCandMapRange", "ckovCandMapRange", 160*scale, 0, 159, 144*scale, 0, 143);

		bgCandMapRange  = new TH2F("bgCandMapRange", "bgCandMapRange", 160*scale, 0, 159, 144*scale, 0, 143);

		ckovCandMapOutRange = new TH2F("ckovCandMapOutRange", "ckovCandMapOutRange", 160*scale, 0, 159, 144*scale, 0, 143);

		bgCandMapOutRange = new TH2F("bgCandMapOutRange", "bgCandMapOutRange", 160*scale, 0, 159, 144*scale, 0, 143);


	

		hSignalAndNoiseMap = new TH2F("Signal and Noise ", "Signal and Noise ; x [cm]; y [cm]",160*scale,0.,159.,144*scale,0,143);

		hSignalMIP = new TH2F("hmip ", "hmip ; x [cm]; y [cm]",160*scale,0.,159.,144*scale,0,143);
		hSignalMIPpc = new TH2F("hmip pc", "hmip pc; x [cm]; y [cm]",160*scale,0.,159.,144*scale,0,143);


		hMaxProton = new TH2F("maxPoss Ckov Proton", "maxPoss Ckov Proton; x [cm]; y [cm]",160*scale,0.,159.,144*scale,0,143);
		hMaxPion = new TH2F("maxPoss Ckov Pion", "maxPoss Ckov Pion; x [cm]; y [cm]",160*scale,0.,159.,144*scale,0,143);
		hMaxPionMinL = new TH2F("maxPoss Ckov Pion min L", "maxPoss Ckov Pion min L; x [cm]; y [cm]",160*scale,0.,159.,144*scale,0,143);
		hMinPionMaxL = new TH2F("minPoss Ckov Pion max L", "minPoss Ckov Pion max L; x [cm]; y [cm]",160*scale,0.,159.,144*scale,0,143);
		hMaxKaon = new TH2F("maxPoss Ckov Kaon", "maxPoss Ckov Kaon; x [cm]; y [cm]",160*scale,0.,159.,144*scale,0,143);
		hMaxProton->SetMarkerColor(kGreen+3);
		hMaxKaon->SetMarkerColor(kRed);

		hMinProton = new TH2F("minPoss Ckov Proton", "minPoss Ckov Proton; x [cm]; y [cm]",160*scale,0.,159.,144*scale,0,143);
		hMinPion = new TH2F("minPoss Ckov Pion", "minPoss Ckov Pion; x [cm]; y [cm]",160*scale,0.,159.,144*scale,0,143);
		hMinKaon = new TH2F("minPoss Ckov Kaon", "minPoss Ckov Kaon; x [cm]; y [cm]",160*scale,0.,159.,144*scale,0,143);


		hMaxProton->SetMarkerColor(kGreen+3);
		hMaxKaon->SetMarkerColor(kRed);
		ckovCandMapRange->SetMarkerColor(kGreen);
		ckovCandMapOutRange->SetMarkerColor(kGreen + 2);
		bgCandMapRange->SetMarkerColor(kRed-1);
		bgCandMapOutRange->SetMarkerColor(kRed+2);

		ckovCandMapRange->SetMarkerStyle(3);
		ckovCandMapOutRange->SetMarkerStyle(2);
		bgCandMapRange->SetMarkerStyle(3);
		bgCandMapOutRange->SetMarkerStyle(2);




		hMinProton->SetMarkerColor(kGreen+3);
		hMinKaon->SetMarkerColor(kRed);
		hMaxPion->SetMarkerColor(kBlue + 4);    // max ckov max L
		hMinPion->SetMarkerColor(kBlue); 				// min ckov min L
		hMinPionMaxL->SetMarkerColor(kBlue); 		// min ckov max L
		hMaxPionMinL->SetMarkerColor(kBlue + 4);// max ckov min L
		/*
		hMaxPion->SetMarkerStyle(3);
		hMinPionMaxL->SetMarkerStyle(3);*/
		hSignalMIPpc->SetMarkerStyle(3);
		hSignalMIP->SetMarkerStyle(3);
		hSignalMIPpc->SetMarkerStyle(3);
		hSignalMIPpc->SetMarkerColor(kRed);
		hSignalAndNoiseMap->SetMarkerStyle(2);

  }



  
	void drawTotalMap()
	{
		TCanvas* tcnvRane = new TCanvas(Form("tcnvRane%d", eventCnt), Form("tcnvRane%d", eventCnt), 1600, 800);
		tcnvRane->cd();
		ckovCandMapRange->Draw();
		ckovCandMapOutRange->Draw("same");
		bgCandMapRange->Draw("same");
		bgCandMapOutRange->Draw("same");
		
	
		hMaxPion->Draw("same");     //  Printf("makeEvent()  hMaxPion->Draw");
		hMinPion->Draw("same");
		hMaxProton->Draw("same");   //Printf("makeEvent()  hMaxProton->Draw");
		hMinProton->Draw("same");
		hMinKaon->Draw("same");
		hMaxKaon->Draw("same");
	}

	void drawMaxRegions()
	{


		const auto trkPC = populatePtr->getPcImp();
		const auto trkRad = populatePtr->getTrackPos();
		TH2F* trkPCMap = new TH2F("trkPCMap ", "trkPCMap; x [cm]; y [cm]",160*20,0.,159.,144*20,0,143);
		TH2F* trkRadMap = new TH2F("trkRadMap ", "trkRadMap; x [cm]; y [cm]",160*20,0.,159.,144*20,0,143);
	  trkRadMap->Fill(trkRad.X(), trkRad.Y());
	  trkPCMap->Fill(trkPC.X(), trkPC.Y());


		trkRadMap->SetMarkerStyle(3);  trkRadMap->SetMarkerColor(kGreen+4);
		trkPCMap->SetMarkerStyle(3);  trkPCMap->SetMarkerColor(kGreen+2);

		TCanvas *thSignalAndNoiseMap = new TCanvas(Form("hSignalAndNoiseMap%d", eventCnt),Form("hSignalAndNoiseMap%d", eventCnt),800,800);  
		thSignalAndNoiseMap->cd();
		hSignalAndNoiseMap->Draw();
		hSignalMIPpc->Draw("same");

		hMaxPion->Draw("same");     //  Printf("makeEvent()  hMaxPion->Draw");
		hMinPion->Draw("same");
		hMaxProton->Draw("same");   //Printf("makeEvent()  hMaxProton->Draw");
		hMinProton->Draw("same");
		hMinKaon->Draw("same");
		hMaxKaon->Draw("same");
	}



	void drawTotalMapAndMaxRegions()
	{
		TCanvas* tcnvRane = new TCanvas(Form("TotalMap%d", eventCnt), Form("TotalMap%d", eventCnt), 1600, 800);

		

		tcnvRane->cd();
		ckovCandMapRange->Draw();
		ckovCandMapOutRange->Draw("same");
		bgCandMapRange->Draw("same");
		bgCandMapOutRange->Draw("same");

		hSignalAndNoiseMap->Draw();


		hMaxPion->Draw("same");     //  Printf("makeEvent()  hMaxPion->Draw");
		hMinPion->Draw("same");
		hMaxProton->Draw("same");   //Printf("makeEvent()  hMaxProton->Draw");
		hMinProton->Draw("same");
		hMinKaon->Draw("same");
		hMaxKaon->Draw("same");
		
				hSignalMIP->Draw("same");
		hSignalMIPpc->Draw("same");
	}

};

class CkovTools {


private: 

  
	// for storing candidates...
	struct Candidate {
		double x, y = 0.;
		double R = 0.;
		double phiL = 0., phi = 0.;
		bool isCandidate = false;
	};

	struct Candidate2 {
		double x, y = 0.;
		//double R = 0.;
		//double phiL = 0., phi = 0.;
		/*uint32_t*/int candStatus = 0;//(000, 001, 010, 011, 100, 110, 101, 111);
	};



  bool print = false; // ef: TODO : later pass this in ctor 

  Populate2* populate2Ptr = nullptr;// =  
	Populate* populatePtr = nullptr;// 

// using array = std::array;
using vecArray3 = std::vector<std::array<double,3>>;
using vecArray2 = std::vector<std::array<double,2>>;

// using arrArray3 = std::vector<std::array<double,3>>;


using segType = std::vector<std::pair<double, double>>;
segType segPionLocal = {{}};

bool kaonStatus = true, pionStatus = true, protonStatus = true;
//TLine* tlinePion;

static constexpr double PI = M_PI;
static constexpr double halfPI = M_PI/2;
static constexpr double twoPI = M_PI*2;


static constexpr double stdDevPion = 0.008; 
static constexpr double stdDevKaon = 0.008; 
static constexpr double stdDevProton = 0.008;
static constexpr float tGap = 8;
static constexpr float  rW = 1.5; // was 1?
static constexpr float  qW = 0.5;
static constexpr float lMax = 1;


static constexpr float CH4GapWidth = 8;
static constexpr float  RadiatorWidth = 1.5; // was 1?
static constexpr float  QuartzWindowWidth = 0.5;



// L value that is generated by simulation, will be changed
//
float L = rW/2;


  // L value for reconstruction
  static constexpr float  EmissionLenght = RadiatorWidth/2;

  float thetaP, phiP, xPC, yPC, xRad, yRad; 
 float nF, nQ, nG;  
 double occupancy;
 std::array<float, 3> ckovHyps;
 std::vector<std::pair<double, double>> photons;


 // instantiate ckovhyp boundaries to "invalid values" 
 // --> later set them if they exceed momentum threshold
 double ckovPionMin = 999, ckovPionMax = 0, ckovKaonMin = 999, ckovKaonMax = 0, ckovProtonMin = 999,ckovProtonMax = 0;

 // double mRMax, mL1Max, mL2Max, mRMin, mL1Min, mL2Min;

 double momentum, mass;

 float cosThetaP, sinThetaP, tanThetaP;
 float cosPhiP, sinPhiP, tanPhiP;
 float trackCkov;
 double xMipLocal, yMipLocal;
 float phiLCurrent, etaCCurrent;


  ArrAndMap* mArrAndMap = nullptr;// new ArrAndMap(eventCnt);

 int eventCnt;
 TVector2 trkPC;

public:
typedef std::vector<std::pair<double, double>> MapType;

void addPhoton(double phiL, double etaC)
{
  photons.push_back(std::make_pair(phiL, etaC));
}


// set new values of phiL and etaC for new photon
void setPhoton(double phiL, double etaC)
{
  phiLCurrent = phiL;
  etaCCurrent = etaC;
}


//   double pc[5] = {xRad,yRad,L,thetaP, phiP};
double refIndexes[3] = {nF, nQ, nG};




CkovTools (double radParams[7], double refIndexes[3], 
           std::array<float, 3> ckovHyps, double occupancy, float trackCkov, int eventCnt)
  : 
    ckovHyps(ckovHyps), occupancy(occupancy) , trackCkov(trackCkov), eventCnt(eventCnt) {  
    
    
  // double radParams[6] = {xRad,yRad,L,thetaP, phiP, randomValue.momentum};
    
  xRad= radParams[0];
  yRad= radParams[1];
  L = radParams[2]; 
  thetaP = radParams[3];
  phiP = radParams[4];
  momentum = radParams[5];	 
  mass = radParams[6];

  nF = refIndexes[0];

	// set tehse to be constant?
  nQ = 1.5787; // ??? 1.5787;// TODO: check this !
  nG = 1.005;
 
  Printf(" CkovTools momentum = %.2f, refFreon = %.2f; ckovHyps : %.2f %.2f %.2f", momentum, nF, ckovHyps[0], ckovHyps[1], ckovHyps[2]);

	if(TMath::IsNaN(ckovHyps[0])){
 	  Printf("Pion CkovHyps is Nan!");
	  setPionStatus(false);
	  // setIsNan()?
	} else {
		setPionStatus(true);
 	  Printf("Pion CkovHyps %.2f", ckovHyps[0]);
  }

  if(TMath::IsNaN(ckovHyps[1])){
	  Printf("Kaon CkovHyps is Nan!");
   	setKaonStatus(false);
	  // setIsNan()?
	} 
	else { 
	  	//setMaxRadius();
  }

  if(TMath::IsNaN(ckovHyps[2])){
  	  Printf("Proton CkovHyps is Nan!");
	  	setProtonStatus(false);
	} else {
		  //setMaxRadius();
	}

	if(getPionStatus()){
	  ckovPionMin = ckovHyps[0] - 4 * stdDevPion;
	  ckovPionMax = ckovHyps[0] + 4 * stdDevPion;
 	  Printf("init CkovTools constructor : getPionStatus() true ! minPion %.2f, maxPion %.2f ", ckovPionMin, ckovPionMax);
  }	else {
 	  Printf("init CkovTools constructor : getPionStatus() was false !");
  }
  
  if(getKaonStatus()){
	  ckovKaonMin = ckovHyps[1] - 4 * stdDevKaon;
  	ckovKaonMax = ckovHyps[1] + 4 * stdDevKaon;
  } if(getProtonStatus()){
		ckovProtonMin = ckovHyps[2] - 4 * stdDevProton;
		ckovProtonMax = ckovHyps[2] + 4 * stdDevProton;
	}

 	cosThetaP = TMath::Cos(thetaP);
  sinThetaP = TMath::Sin(thetaP);
  tanThetaP = TMath::Tan(thetaP);

cosPhiP = TMath::Cos(phiP);
sinPhiP = TMath::Sin(phiP);


TRotation rotZ; rotZ.RotateZ(phiP);
TRotation rotY; rotY.RotateY(thetaP);

TVector3 ip(0,0,rW-L+tGap+qW);
TVector3 op; op = rotZ*rotY*ip;

xMipLocal = tanThetaP*cosPhiP*(rW-L + tGap + qW);
yMipLocal = tanThetaP*sinPhiP*(rW-L + tGap + qW);

xPC = xMipLocal + xRad;
yPC = yMipLocal + yRad;


// fyll : 


trkPC.Set(xPC, yPC); // TODO :check it equals populate.getPcImp();

 
TVector2 trkPos(xRad, yRad);
TVector3 trkDir; 
trkDir.SetMagThetaPhi(1, thetaP, phiP);




populatePtr = new Populate(trkPos, trkDir, nF);


  mArrAndMap = new ArrAndMap(eventCnt, populatePtr);
  (mArrAndMap->hSignalMIP)->Fill(xRad, yRad);
  (mArrAndMap->hSignalMIPpc)->Fill(xPC, yPC);


populate2Ptr = new Populate2(trkPos, trkDir, nF);

Populate populate(trkPos, trkDir, nF);



Printf("init Ckovtools \n MIP Root : %f %f %f \n MIP local %f %f",op.Px(),op.Py(),op.Pz(),xMipLocal,yMipLocal);
      // constructor body goes here, if needed

Printf("Ckovtools :: trkPC %.2f %.2f",trkPC.X(), trkPC.Y());
      // constructor body goes here, if needed

/*
mRMax = getR_Lmax(ckovPionMax, halfPI);
mL2Max = getR_Lmax(ckovPionMax, 0);
mL1Max = getR_Lmax(ckovPionMax, PI);


// needs to be changed if to be used! now uses L = lMax
mRMin = getR_Lmax(ckovProtonMin, halfPI);
      mL1Min = getR_Lmax(ckovProtonMin, PI);
      mL2Min = getR_Lmax(ckovProtonMin, 0);*/ 
		
			/*for(const auto& c : ckovHyps) {
				auto R = getR_Lmax(c, halfPI);
				auto l2 = getR_Lmax(c, 0);
				auto l1 = getR_Lmax(c, PI);
				Printf(" CkovTools : CkovHyp %f, R %f, l1 %f l2 %f", c, R, l2, l1);
      } */  
}




/*void setSegLines(TLine* l1, TLine* l2, TLine* l3, TLine* l4)
{
	    
}*/

void setSegment(const segType& segPion, segType& segPionRot)
{
  Printf("Ckovtools : enter setsegment!");
  // segType segPionRot; 
  // segPionRot.resize(segPion.length());


  /*
  auto sX = segPion[0].first, sY = segPion[0].second;
  auto eX = segPion[0].first, eY = segPion[0].second;
  local2PhiRing(sX, sY, xMipLocal, yMipLocal);

  local2PhiRing(eX, eY, xMipLocal, yMipLocal);
	tlinePion->SetX1(sX); tlinePion->SetY1(sY); 
	tlinePion->SetX1(eX); tlinePion->SetY1(eY); */ 

  int i = 0;
  segPionLocal.reserve(segPion.size());
  segPionLocal.resize(segPion.size());
  Printf("sizes = %zu 1 %zu!", segPionRot.size(), segPionLocal.size());
  for(const auto& p : segPion) 
  {
    auto x = p.first, y = p.second;
    local2PhiRing(x, y, xMipLocal, yMipLocal);
    segPionRot[i] = std::make_pair(x,y);
    segPionLocal[i] = std::make_pair(x,y);
    //Printf("Ckovtools : setsegment! x %.2f y %.2f", x, y);
		i++;
  }    
  Printf("Ckovtools : exit setsegment!");
  Printf("sizes = %zu %zu!", segPionRot.size(), segPionLocal.size());
}

void setPionStatus(bool status)
{ 
  pionStatus = status;
}

bool getPionStatus() const
{
  return pionStatus;
}

void setKaonStatus(bool status)
{ 
  kaonStatus = status;
}

bool getKaonStatus() const
{
  return kaonStatus;
}

void setProtonStatus(bool status)
{ 
  protonStatus = status;
}

bool getProtonStatus() const
{
  return protonStatus;
}


double getMinCkovKaon()
{
  return ckovKaonMin;
}

double getMaxCkovKaon()
{
  return ckovKaonMax;
}

double getMinCkovPion()
{
  return ckovPionMin;
}

double getMaxCkovPion()
{
  return ckovPionMax;
}

double getMinCkovProton()
{
  return ckovProtonMin;
}

double getMaxCkovProton()
{
  return ckovProtonMax;
}


double getPhiP()
{
  return phiP;
}

double getThetaP()
{
  return thetaP;
}

// TODO: add nQ
// only consider photons in the correct range:


// TODO: add nQ
// only consider photons in the correct range:


// x, y, etaC

std::vector<std::pair<double, double>> segment(std::vector<std::array<double, 3>>& cherenkovPhotons, MapType& bins) { 
  MapType pionCandidates, kaonCandidates, protonCandidates;
  
  Printf("ckovTools enter  ckovTools.segment"); 	 
  //const auto infString = Form("localRef #Theta_{p}  = %.4f #Phi_{p} = %.4f L = %.2f \n #Theta_{C} = %.4f maxCkov = %.4f ; x [cm]; y [cm]", thetaP,phiP, L,trackCkov,ckovPionMax); 

  //TH2F *localRefMIP = new TH2F("localRefMIP ", infString,800,-40.,-40.,900,-40.,40.);

  
  /*
  TVector2 trkPos(xRad, yRad);
  TVector3 trkDir; 
  trkDir.SetMagThetaPhi(1, thetaP, phiP);
  Populate* populatePtr = new Populate(trkPos, trkDir, nF);
  Populate populate(trkPos, trkDir, nF);*/


  // void setArrayMin(Populate* populate, double etaTRS, vecArray3 inPutVector) 


  const size_t kN = 150;

	//array<array<double, 3> ,kN> arrMaxPion; 

  

	// disse kan brukes istedet for maxKaonVec...?
	vecArray3 arrMaxPion, arrMinPion, arrMinProton, arrMaxProton, arrMinKaon, arrMaxKaon;

	vecArray2 arrMaxPionPos, arrMinPionPos, arrMinProtonPos, arrMaxProtonPos, arrMinKaonPos, arrMaxKaonPos;

	// do not reserve size? NO prob just dont resize :) 
	arrMaxPion.reserve(kN);
	arrMinPion.reserve(kN);
	//arrMaxPion.resize(kN);
	//arrMinPion.resize(kN);

	arrMaxPionPos.reserve(kN);
	arrMinPionPos.reserve(kN);
	arrMaxPionPos.resize(kN);
	arrMinPionPos.resize(kN);

	// arrMaxPion.resize(kN);
		
	// check if candidate can be proton (i.e., that it exceeds momentum threshold)
	if(getProtonStatus()) {
		arrMaxProton.reserve(kN);
		arrMinProton.reserve(kN);
		//arrMaxProton.resize(kN);
		//arrMinProton.resize(kN);

		arrMinProtonPos.reserve(kN);
		arrMaxProtonPos.reserve(kN);
		arrMinProtonPos.resize(kN);
		arrMaxProtonPos.resize(kN);
				
		// also later add max, i forste runde sjekk at ckov i {minProton, maxPion}
		Printf("calling setArrayMin w getMinCkovProton() = %.2f", getMinCkovProton());
  	setArrayMin(getMinCkovProton(), arrMinProton, arrMinProtonPos, kN);

		Printf("calling setArrayMax w getMaxCkovProton() = %.2f", getMaxCkovProton());
  	setArrayMax(getMaxCkovProton(), arrMaxProton, arrMaxProtonPos, kN);
	}
	

	// check if candidate can be Kaon (i.e., that it exceeds momentum threshold)
	if(getKaonStatus()) {
		arrMaxKaon.reserve(kN);
		arrMinKaon.reserve(kN);				

		//arrMaxKaon.resize(kN);
		//arrMinKaon.resize(kN);	

		arrMaxKaonPos.reserve(kN);
		arrMinKaonPos.reserve(kN);				

		arrMaxKaonPos.resize(kN);
		arrMinKaonPos.resize(kN);

		Printf("calling setArrayMin w getMinCkovKaon() = %.2f", getMinCkovKaon());
  	setArrayMin(getMinCkovKaon(), arrMinKaon, arrMinKaonPos, kN);

		Printf("calling setArrayMax w getMaxCkovKaon() = %.2f", getMaxCkovKaon());
  	setArrayMax(getMaxCkovKaon(), arrMaxKaon, arrMaxKaonPos, kN);
	}

	


	Printf("BF : Length of elem vectors : arrMaxPion %d", arrMaxPion.size());	
  Printf("calling setArrayMax w getMaxCkovPion() = %.2f", getMaxCkovPion());
  setArrayMax(getMaxCkovPion(), arrMaxPion, arrMaxPionPos, kN);
	Printf("AFTER : Length of elem vectors : arrMaxPion %d", arrMaxPion.size());	



	Printf(" BF : Length of elem vectors : arrMinPion %d", arrMinPion.size());
  Printf("calling setArrayMin w getMinCkovPion() = %.2f", getMinCkovPion());
  setArrayMin(getMinCkovPion(), arrMinPion, arrMinPionPos, kN);
	Printf("AFTER : Length of elem vectors : arrMinPion %d", arrMinPion.size());
	// fill all the values in the maps 
	
	Printf("Length of elem vectors : arrMinPion %d", arrMinPion.size());
	Printf("Length of elem vectors : arrMaxPion %d", arrMaxPion.size());
	
	if(true) {

		Printf("		fillMapFromVec(mArrAndMap->hMaxPion, arrMaxPion);// map, array");
		fillMapFromVec(mArrAndMap->hMaxPion, arrMaxPionPos);// map, array
		fillMapFromVec(mArrAndMap->hMaxKaon, arrMaxKaonPos);// map, array
		fillMapFromVec(mArrAndMap->hMaxProton, arrMaxProtonPos);// map, array

		fillMapFromVec(mArrAndMap->hMinPion, arrMinPionPos);// map, array
		fillMapFromVec(mArrAndMap->hMinKaon, arrMinKaonPos);// map, array
		fillMapFromVec(mArrAndMap->hMinProton, arrMinProtonPos);// map, array
	}




  // get track impacrt point at PC


	const int scale = 2;


  Printf("ckovTools segment : exit const auto& p : segPionLocal"); 	 

  // TODO: ckovHyps : get std-dev for Theta_ckov of pion kaon and proton from the values theta_i

 // initialize recon with track input params 
 

 // TODO: change this, xMipLocal just placeholder
 // not sure if xPC simply is obtained like this
 double xPC = xMipLocal, yPC = yMipLocal;
 local2GlobalRef(xPC, yPC);

// ReconG(double _theta, double _phi, double _xRad, double _yRad, double _xPC, double _yPC, double n, double _etaC) : refIdx(n), etaC(_etaC)
 //ReconG reconG(thetaP, phiP, xRad , yRad, xRad + xMipLocal,  yRad + yMipLocal, nF);

 Printf("dX %f dY %f ", xMipLocal, yMipLocal);

 // TH2F *hNoiseMap = new TH2F(" Noise ", " Noise ; x [cm]; y [cm]",160,0.,159.,144,0,143);
 int numPhotons= 0;


  const auto area = 144*156; 
 
 
 
  // number of bg-photons produced
  const auto numBackgroundPhotons = static_cast<int>(area*occupancy*0); 

  // number of correctly identified ckov photons
	int numFoundActualCkov = 0;
		
		
	// in case ckov photons produced are out of range of chamber dimensions:	
	int numActualCkov = 0;
	
	
  // number of bg-photons labellled as being withing range:
	int numBackgroundLabeledCkov = 0;

 	std::vector<std::array<double, 3>> backGroundPhotons(numBackgroundPhotons);
  // std::vector<std::pair<double, double>> backGroundPhotons(numBackgroundPhotons);
   //   Printf("CkovTools segment backGroundPhotons vector created");
  std::random_device rd;
  std::mt19937 gen(rd());

  // TODO : based on ckovMaxProton, add bg here : 
  // in  phiRing ref sys

  std::uniform_real_distribution<> dis1(0, 156);
  std::uniform_real_distribution<> dis2(0, 144);

  for (auto& pair : backGroundPhotons) {
    pair[0] = dis1(gen);
    pair[1] = dis2(gen);
    pair[2] = 0.;// ph for etaC here 
     // Printf("CkovTools segment created bg x%f y%f", pair.first, pair.second);
  }

  /*
  
  
  */ 



  /*
	// add background photons here :
  for(const auto& p :cherenkovPhotons){

    const auto& xDif = p[0] - (xMipLocal + xRad); 
    const auto& yDif = p[1] - (yMipLocal + yRad) ; 
    auto R = TMath::Sqrt(xDif*xDif+yDif*yDif);
    Printf(" Ckovtools segments cherenkovPhotons : x %f y %f R = %f", p[0], p[1], R);
    numPhotons++;
  } */ 

  //TH2F *hSignalMap = new TH2F("Signal", Form("Cherenkov Photons = %d  #Theta_{p}  = %f #Phi_{p} = %f ; x [cm]; y [cm]", numPhotons, thetaP,phiP),160,0.,159.,144,0,143);   

  Printf("ckovtools segmetns cherenkovPhotons size  = %zu", cherenkovPhotons.size()); 
 
  MapType filledBins;
  // placeholders for the above : 

  // these are in local coordinate system
  // create bbox of minimal possible ckov
  // Proton with eta_c = theta_c_proton - 3*std_dev_proton
	/*
  double xLU = -mL1Max;
  double yLU = mRMax;

  double xRU = mL2Max;
  double yRU = mRMax;

  double xLD = -mL1Max;
  double yLD = -mL1Max;

  double xRD = mL2Max;
  double yRD = -mRMax;


	phiRing2Local(xLU, yLU, xMipLocal,yMipLocal);
	phiRing2Local(xRU, yRU, xMipLocal,yMipLocal);
	phiRing2Local(xRD, yRD, xMipLocal,yMipLocal);
	phiRing2Local(xLD, yLD, xMipLocal,yMipLocal);

  
  const auto l1Max = mL1Max;//TMath::Sqrt((xMaxPhiPi-xMipRing)*(xMaxPhiPi-xMipRing) + (yMaxPhiPi-yMipRing)*(yMaxPhiPi-yMipRing));
  const auto l2Max = mL2Max;//TMath::Sqrt((xMaxPhi0-xMipRing)*(xMaxPhi0-xMipRing) + (yMaxPhi0-yMipRing)*(yMaxPhi0-yMipRing)); 
  */ 
    /*
    TBox* localBox = new TBox(xMaxPhi0, yMaxPhiPi - mRMax, xMaxPhiPi, yMaxPhiPi+mRMax);
    localBox->SetLineColor(kGreen);
    
    const auto& coords1 = local2Global(xMaxPhi0, yMaxPhi0 - mRMax);
    const auto& coords2 = local2Global(xMaxPhiPi, yMaxPhiPi - mRMax);
		const auto& coords3 = local2Global(xMaxPhi0, yMaxPhi0+mRMax);
		const auto& coords4 = local2Global(xMaxPhiPi, yMaxPhiPi+mRMax);
    
    TBox* globalBoxSignal = new TBox(coords1.first, coords1.second, coords2.first, coords2.second); */
    
		/*
		local2GlobalRef(xLU,yLU);
		local2GlobalRef(xRU,yRU);
		local2GlobalRef(xLD,yLD);
		local2GlobalRef(xRD,yRD);


    TLine *tlineLeftGlobal = new TLine(xLU, yLU, xLD, yLD);

    TLine *tlineRightGlobal = new TLine(xRU, yRU, xRD, yRD);

    TLine *tlineUpGlobal = new TLine(xLU, yLU, xRU, yRU);

    
    TLine *tlineDownGlobal = new TLine(xLD, yLD, xRD, yRD);
    
	  tlineUpGlobal->SetLineColor(kGreen);
	  tlineDownGlobal->SetLineColor(kBlue); */ 
       
    //globalBoxSignal->SetLineColor(kGreen);
    /*
    TLine *tlineUpLocal = new TLine(xMaxPhi0,yMaxPhi0 - mRMax, xMaxPhiPi, yMaxPhiPi - mRMax);
    TLine *tlineDownLocal = new TLine(xMaxPhi0,yMaxPhi0 + mRMax, xMaxPhiPi, yMaxPhiPi + mRMax);
		*/


  /*
    TLine *tlineUpLocal = new TLine(-mL1Max,mRMax, mL2Max, mRMax);

    TLine *tlineDownLocal = new TLine(-mL1Max,-mRMax, mL2Max, -mRMax);



	
	auto l1 = mL1Max; auto r1 = mRMax; auto r2 = mRMax; auto l2 = mL2Max;
	phiRing2Local(l1,r1,xMipLocal,yMipLocal);
	phiRing2Local(l2,r2,xMipLocal,yMipLocal);
  TLine *tlineUpLocalR = new TLine(-l1, r1,l2,r1);
  TLine *tlineDownLocalR = new TLine(-l1,-r2,l2,-r2);
  tlineUpLocal->SetLineColor(kGreen);
  tlineDownLocal->SetLineColor(kBlue);

  tlineUpLocalR->SetLineColor(kGreen);
  tlineDownLocalR->SetLineColor(kBlue);
    
  const auto rMax2 = r2;
 */ 
  //Printf("rMax2 = %f w ckov = %f", rMax2, ckovPionMax);
  // populate with background:

 


  //Printf("CkovTools segment : rMax %f l1Max %f l2Max %f Area %f ", rMax, l1Max, l2Max, area);


 // Printf("CkovTools segment backGroundPhotons vector numBackgroundPhotons = %d", numBackgroundPhotons);
  // std::vector<std::pair<double, double>> backGroundPhotons(numBackgroundPhotons);
  
  
   //   Printf("CkovTools segment backGroundPhotons vector created");
  /*std::random_device rd;
  std::mt19937 gen(rd());

  // TODO : based on ckovMaxProton, add bg here : 
  // in  phiRing ref sys

  std::uniform_real_distribution<> dis1(0, 156);
  std::uniform_real_distribution<> dis2(0, 144);

  for (auto& pair : backGroundPhotons) {
    pair.first = dis1(gen);
    pair.second = dis2(gen);
     // Printf("CkovTools segment created bg x%f y%f", pair.first, pair.second);
  } */ 

  std::vector<std::pair<double, double>> photonCandidates;


  // const auto infString2 = Form("globalREf #Theta_{p}  = %f #Phi_{p} = %f ; x [cm]; y [cm]", thetaP,phiP); 
  	// TH2F *globalREfMIP = new TH2F("globalREfMIP ", infString2,160,0.,159.,144,0,143);

	 // const auto& coordsMIP = local2Global(xMipLocal, yMipLocal);
	 // globalREfMIP->Fill(coordsMIP.first, coordsMIP.second);
 
    double xML = xMipLocal, yML = yMipLocal;
    //localRefMIPUnrot->Fill(xMipLocal, yMipLocal);
    //local2PhiRing(xML,yML,xML,yML);  
    //localRefMIP->Fill(xML,yML);
    //localRefMIP2->Fill(xML,yML);


  // change datatype, should hold at least x, y..+? {later + cluster-size...}

  // Make posPhoton to send to checkRange




  std::random_device rd2;
  std::mt19937 genUnc(rd2());
  const float sx = 0.8, sy = 0.84;
  const float uncX = sx/TMath::Sqrt(12), uncY = sy/TMath::Sqrt(12);
  // TODO: should this value be divided by two (four?)?
  //

  // TODO: add sqrt(uncX^2 + uncY^2) to radius of ckov hyps?



	// TODO: here the unc is - unc to + unc, in reality we then add 2* as much as s
  // it should be
  // within the ~1mm unc range, generate a value +- unc to add to 
  // ideal cluster-position
  std::uniform_real_distribution<> disX(-uncX/2, uncX/2);
  std::uniform_real_distribution<> disY(-uncY/2, uncY/2);


  // here : loop through all backGroundPhotons and cherenkovPhotons
  
  int iPhotCount = 0;



	// obv change this later to also add bg photons
	const size_t numPhotonsTemp = cherenkovPhotons.size();
	std::vector<Candidate> pionCands, kaonCands, protonCands;
	pionCands.reserve(numPhotonsTemp);
	kaonCands.reserve(numPhotonsTemp);
	protonCands.reserve(numPhotonsTemp);

	std::vector<Candidate2> candCombined;
	candCombined.reserve(numPhotonsTemp);
	
	//inPutVector.emplace_back(std::array<double, 3>{phiL, phiR, r});
	
	
	std::vector<std::array<double, 3>> ckovAndBackground;
	

	// add unc to cherenkov photons
  for(const auto& photons : cherenkovPhotons) {  
	  auto dX = disX(genUnc);
    auto dY = disY(genUnc);
    auto x = photons[0] + dX;
    auto y = photons[1] + dY;
    auto etaC = photons[2] ; 
  	ckovAndBackground.emplace_back(std::array<double, 3>{x, y, etaC});
	} 
	
  for(const auto& p : backGroundPhotons) {
  	ckovAndBackground.emplace_back(std::array<double, 3>{p});
	} 
	
	Printf("length ckovPhotons %zu length background %zu length total %zu",cherenkovPhotons.size(), backGroundPhotons.size(), ckovAndBackground.size());
	
	// gjør dette mer elegant snre
	
  // ckovAndBackground


  Printf(" enter const auto& photons : ckovAndBackground");
  for(const auto& photons : ckovAndBackground) {
		iPhotCount++;


  	Printf("%d", iPhotCount);
    // photX photY are "ideal" values, should add some uncertainty to them... 
    //auto dX = disX(genUnc);
    //auto dY = disY(genUnc);
    const auto& x = photons[0], y = photons[1];
    //auto x = photons[0] + dX;
    //auto y = photons[1] + dY;
    //Printf("ckovloop Unc | xIdeal %.2f, uncX %.2f , genX = %.2f", photons[0], dX, x);
    //Printf("ckovloop Unc | yIdeal %.2f, uncY %.2f , genY = %.2f", photons[1], dY, y);
    
    if(!(x > 0 && x < 156.0 && y > 0 && y < 144)) {
      Printf("Photon outside of chamber dimensions!");
      continue;
    }

    const auto& etaC = photons[2];


    double xL = x, yL = y;
    double xG = x, yG = y;
    /*localRefUnrot->Fill(x,y);
    // transform to phiRing ref-system
    local2PhiRing(xL, yL, xMipLocal, yMipLocal);
    localRef->Fill(xL, yL);*/ 


		/*
    Printf("ckovtools cherenkov photons x > xMaxPhiPi && x < xMaxPhi0 && y > yMaxPhiPi && y < yMaxPhi0");
    Printf("ckovtools cherenkov photons x  %f > xMaxPhiPi %f && x %f < xMaxPhi0 %f && y %f > yMaxPhiPi  %f && y %f < yMaxPhi0 %f",  x, xMaxPhiPi, x , xMaxPhi0 , y , yMaxPhiPi , y , yMaxPhi0);*/
    
	
    //Printf("\nckovtools cherenkov photons x  %f > -mL1Max %f && x %f < mL2Max %f && y %f > -mRMax  %f && y %f < mRMax %f\n",  x, -mL1Max, x , mL2Max , y , -mRMax , y , mRMax);

    double thetaCer, phiCer;
    //local2GlobalRef(xG, yG);
    // double cluX, double cluY, double& thetaCer, double& phiCer
    

    //Printf("CkovTools segment thetaCer %f phiCer %f", thetaCer, phiCer);


    auto xAbs = TMath::Abs(x);
    auto yAbs = TMath::Abs(y);

    //Printf("\nckovtools cherenkov photons xAbs  %f > mL1Max %f && x %f < mL2Max %f && yAbs %f > mRMax  %f && y %f < mRMax %f\n",  x, mL1Max, xAbs , mL2Max , yAbs , mRMax , y , mRMax);




    // not really helpful? 
    if(true){
    //if(x > -mL1Max && x < mL2Max && y > -mRMax  && y < mRMax){

    //double thetaCer, phiCer;
    //reconG.findPhotCkov(xG, yG, thetaCer, phiCer);	
    //auto ckov = thetaCer;


      bool withinRange = true; 
      // check if the coordinates also corresponds to one of the possible cherenkov candidates


      // TODO : check if this method is wrong??
      // use here instead method from Recon.cxx

      //reconG.findPhotCkov(xG, yG, thetaCer, phiCer, etaC);	
      //auto ckov = thetaCer;
			//auto ckov2 = ckov;
      // const auto& ckov2 = getCkovFromCoords(xG, yG, phiP, thetaP, nF, nQ, nG, etaC);      


      // Printf("CkovTools segment | actual ckov = %.3f | bgstdy method:  %.3f | ReconMEthods:  thetaCer %f phiCer %f ", etaC, ckov2,  thetaCer, phiCer);



const TVector2 posPhoton(x, y);

// find reference  radius to be compared to segmetation radiuses 
const auto rPhoton = (posPhoton - trkPC).Mod();

const auto phiPhoton = (posPhoton - trkPC).Phi();




const auto pc = populatePtr->getPcImp();

      
	
	/*
	***************************************************************************
	***************************************************************************
	***************************************************************************
	***************************************************************************
	***************************************************************************
	***************************************************************************
	***************************************************************************
	***************************************************************************
	*/ 
 
  // iteration phiL approach init
	// checking if inside outer pionRadius (Lmin, etaMax)

	
  bool isPhotonProtonCand = false, isPhotonKaonCand = false, isPhotonPionCand = false;

	bool isMaxProtonOk = false, isMinProtonOk = false, isMaxKaonOk = false, isMinKaonOk = false, isMaxPionOk = false, isMinPionOk = false;


	Printf("\n\n ============================================================ \n");		
	if(getProtonStatus() and getKaonStatus() and getPionStatus()) {
		
		Printf("Pion%d can be Pion Kaon and Proton from p-hyp \n", iPhotCount);	
		// verify thatrPhoton > rMinProton(@ phiEstimated = phiPhoton)

		Printf("\n\npopulate2Ptr->checkUnder(posPhoton, rPhoton, phiPhoton, above, arrMinProton, getMinCkovProton());");
  	isMinProtonOk = populate2Ptr->checkUnder(posPhoton, rPhoton, phiPhoton, arrMinProton, getMinCkovProton(), "Proton");

		// check if rPhoton > rMax(@ phiEstimated = phiPhoton)
		if(isMinProtonOk) {
			Printf("isMinProtonOk = true");
			Printf("\n\npopulate2Ptr->checkOver(posPhoton, rPhoton, phiPhoton, arrMaxPion, getMaxCkovPion());");

			isMaxPionOk = populate2Ptr->checkOver(posPhoton, rPhoton, phiPhoton, arrMaxPion, getMaxCkovPion(), "Pion");
			
			// this means rProtonMin < rPhoton < rMaxPion
			if(isMaxPionOk) {
				
				
				// this means rProtonMax > rPhoton; then also the rPh < rPionMax and rPh < rKaonMax

				isMaxProtonOk = populate2Ptr->checkOver(posPhoton, rPhoton, phiPhoton, arrMaxProton, getMaxCkovProton(), "Proton");
				Printf("===================================="); 
				if(isMaxProtonOk) {
					// we have proton-candiate
					isPhotonProtonCand = true;
					Printf("Photon%d is a Proton Candidate", iPhotCount); 
					// isMaxPionOk = true;
					isMaxKaonOk = true;
				}	else {
					Printf("Photon%d not a Proton Candidate", iPhotCount); 
				}
				Printf("====================================\n"); 
				Printf("===================================="); 
				// this means rPionMin < rPhoton --> and also rKaonMin < rPhoton



				isMinPionOk = populate2Ptr->checkUnder(posPhoton, rPhoton, phiPhoton, arrMinPion, getMinCkovPion(), "Pion");
				if(isMinPionOk) {
					// we have pion-candiate
					isPhotonPionCand = true;
					Printf("Photon%d is a Pion Candidate", iPhotCount); 
					isMinKaonOk = true;
				}	else {
					Printf("Photon%d not a Pion Candidate", iPhotCount); 
				}			
					Printf("===================================="); 				
				
				// if rPhoton < rPhotonMin, check if rPhoton > rKaonMin
				if(!isMinKaonOk) {
					// isMinKaonOk = populate2Ptr->checkCond()
				isMinKaonOk = populate2Ptr->checkUnder(posPhoton, rPhoton, phiPhoton, arrMinKaon, getMinCkovKaon(), "Kaon");
				}
				
				Printf("\n===================================="); 
				// only check isMaxKaon if isMinKaonOk
				if(isMinKaonOk) {
					if(!isMaxKaonOk) {
						isMaxKaonOk = populate2Ptr->checkOver(posPhoton, rPhoton, phiPhoton, arrMaxKaon, getMaxCkovKaon(), "Kaon");
						// isMaxKaonOk = populate2Ptr->checkCond()
					} 

					if(isMaxKaonOk and isMinKaonOk) {
						Printf("Photon%d is a Kaon Candidate", iPhotCount); 
						// we have kaon cand
						isPhotonKaonCand = true;
					} else {
						Printf("Photon%d not a Kaon Candidate", iPhotCount); 
					}
				} else {
					Printf("Photon%d not a Kaon Candidate", iPhotCount); 
				}
				Printf("====================================\n"); 
							
			}  // end if isMaxPionOk == true
						
			// else isMaxPionOk == false
			else {
				Printf("===============================================================");	
				Printf("\n Photon%d not a Hadron Candidate :\n rPhoton %.2f > rPionMax",iPhotCount,  rPhoton);
				Printf("===============================================================");	
			}
			
		}	// end if isMinProtonOk == true
		
		// else isMinProtonOk == false
		// this means  rPhoton < rProtonMin : 
		else {
			Printf("===============================================================");	
			Printf("\n Photon%d not a Hadron Candidate :\n rPhoton %.2f < rProtonMin",iPhotCount, rPhoton);
			Printf("===============================================================");	
		}

	
	} // end if getProtonStatus() and getKaonStatus and getPionStatus

	// if Proton not is possible because momentum-threshold is not exceeded:
	else if(!getProtonStatus() and getKaonStatus() and getPionStatus()) {
		Printf("\n====================================================="); 
		Printf("Photon%d can be Pion and Kaon from p-hyp", iPhotCount);	
		// verify that rPhoton > rMinKaon(@ phiEstimated = phiPhoton)
		Printf("\npopulate2Ptr->checkUnder(posPhoton, rPhoton, phiPhoton, arrMinProton, getMinCkovProton());");
  	isMinKaonOk = populate2Ptr->checkUnder(posPhoton, rPhoton, phiPhoton, arrMinKaon, getMinCkovKaon(), "Kaon");


		// this means rProtonMin < rPhoton < rMaxPion
		if(isMinKaonOk) {
			isMaxPionOk = populate2Ptr->checkOver(posPhoton, rPhoton, phiPhoton, arrMinPion, getMinCkovPion(), "Pion");

			if(isMaxPionOk) {

				isMinPionOk = populate2Ptr->checkUnder(posPhoton, rPhoton, phiPhoton, arrMinPion, getMinCkovPion(), "Pion");
				
			  //isMinPionOk = ...
				// this means rPionMin < rPhoton
				if(isMinPionOk) {
					isPhotonPionCand = true;	
					Printf("Photon%d is Pion Candiate", iPhotCount); 
					// we have pion-candiate
				}	

				isMaxKaonOk = populate2Ptr->checkOver(posPhoton, rPhoton, phiPhoton, arrMaxKaon, getMaxCkovKaon(), "Kaon");
				// isMaxKaonOk = populate2Ptr->
				if(isMaxKaonOk) {
					isPhotonKaonCand = true;
					Printf("Photon%d is Kaon Candiate", iPhotCount); 
					// we have kaon-candiate
				}							
			}	
			else {
				Printf("Photon%d could be Pion/Kaon (from p-hyp) but was out of radius-range",iPhotCount);
			}	
		} // end if isMinKaonOk
		else {
			Printf("Photon%d could be Pion/Kaon (from p-hyp) but was out of radius-range",iPhotCount);
		}
	} // end else if (!getProtonStatus() and getKaonStatus() and getPionStatus())



	// denne ser litt rar ut?

	// if neither Proton or Kaon is possible because momentum-threshold is not exceeded:
	else if(!getProtonStatus() and !getKaonStatus() and getPionStatus()) {
		Printf("Photon%d can be Pion from p-hyp", iPhotCount);	

		Printf("populate2Ptr->checkUnder(posPhoton, rPhoton, phiPhotob, arrMaxPion, getMaxCkovPion());",iPhotCount);

		isMaxPionOk = populate2Ptr->checkOver(posPhoton, rPhoton, phiPhoton, arrMaxPion, getMaxCkovPion(), "Pion");
		
		
		if(isMaxPionOk) {

			Printf("populate2Ptr->checkOver(posPhoton, rPhoton, phiPhoton, arrMinPion, getMaxCkovPion());",iPhotCount);
			
			bool isMinPionOk = populate2Ptr->checkUnder(posPhoton, rPhoton, phiPhoton, arrMinPion, getMinCkovPion(), "Pion");
			
			if(isMinPionOk) {
				// we have pion candidate
				Printf("Photon%d is Pion Candidate", iPhotCount);	
				isPhotonPionCand = true; 
			} else {
				Printf("Photon%d could only have been Pion, but didnt fall within radius-range",iPhotCount);
			}
		} 
		else {
			Printf("Photon%d could only have been Pion, but didnt fall within radius-range",iPhotCount);
		}			
		
	} // end else if(!getProtonStatus() and !getKaonStatus() and getPionStatus())

	// this should not really be "possible", but can be if we have another candidate that is not pion/kaon/proton
	else {
		Printf("Photon%d was not Pion Kaon or Proton from p-hyp",iPhotCount);
		// we got particle that should not be able to be pion, kaon or proton
	}


	
	Printf("\n===============================================================");	
	Printf("	Photon%d  Candidate: Pion = %d Kaon = %d Proton %d", iPhotCount, isPhotonPionCand, isPhotonKaonCand, isPhotonProtonCand);
	Printf("===============================================================\n");	


	// correctly identified cherenkov photons
	if((isPhotonPionCand or isPhotonKaonCand or isPhotonProtonCand) and etaC != 0) {
	  numFoundActualCkov ++;
	}

	// correctly identified cherenkov photons
	if(!(isPhotonPionCand or isPhotonKaonCand or isPhotonProtonCand) and etaC != 0) {
	  Printf("Ckov %.2f", etaC);
	  
	  print = true;
		//throw std::invalid_argument("Photon Ckov not found!??");
	}
	
	if(etaC != 0) {
		numActualCkov ++;
	}

	
	

	// number of found bg photons that are labelled as within proper range:
	if((isPhotonPionCand or isPhotonKaonCand or isPhotonProtonCand) and etaC == 0) {
	  numBackgroundLabeledCkov ++;	
	}

	// struct with x, y, phiLocal, phi, R + bool isCandidate? 
	// or just x, y, /* phiLocal,*/ phi, R; where x, y... = 0 if !isCandidate

	// better to not store rPhot and phiPhot? this can easily be obtained in ML photon
	/*pionCands.emplace_back(Candidate{x, y, rPhoton, phiPhoton, isPhotonPionCand});	
	kaonCands.emplace_back(Candidate{x, y, rPhoton, phiPhoton, isPhotonKaonCand});	
	protonCands.emplace_back(Candidate{x, y, rPhoton, phiPhoton, isPhotonProtonCand});
	*/ 


	int cStatus = 4*static_cast<int>(isPhotonPionCand) + 2*static_cast<int>(isPhotonKaonCand) + 1*static_cast<int>(isPhotonProtonCand);
	candCombined.emplace_back(Candidate2{x, y, cStatus});



	if(true) {
		// TODO: plot as same in maxPion..Range
		if( x > 0 && x < 156 && y < 144 && y > 0) {
			// this means it falls within range
			if(cStatus != 0) {

				// this means  candidate is a ckov-photon
				if(etaC != 0) {
					(mArrAndMap->ckovCandMapRange)->Fill(x,y);
				}	else { // falls within range, but is bg
					(mArrAndMap->bgCandMapRange)->Fill(x,y);
				}
			}
			// falls out of range
			else {
				// this means  candidate is a ckov-photon, but out of range
				if(etaC != 0) {
					(mArrAndMap->ckovCandMapOutRange)->Fill(x,y);
				}	else { // falls out of range, and is bg
					(mArrAndMap->bgCandMapOutRange)->Fill(x,y);			
				}
			}
		}
  }

	// or store as 1 vector, where candidate status 8 => 2^3 (000, 001, 010, 011, 100, 110, 101, 111)
	
	// phiL, phi, R of maxPionVec vectorh



	/*
	***************************************************************************
	***************************************************************************
	***************************************************************************
	***************************************************************************
	***************************************************************************
	***************************************************************************
	***************************************************************************
	***************************************************************************
	*/ 


	// *<init ClosedForm segm>//		
  // "Closed form approch" init    

  // check if inside pionMax and outside protonMin



	/*
	Printf("\n\n");
  TVector2 rPosPion;
  bool pionBelow = populate.checkRangeBelow(posPhoton, getMaxCkovPion(), rPosPion);

	TVector2 rPosPionMin;
  populate.checkRangeBelow(posPhoton, getMinCkovPion(), rPosPionMin);

  mapPhotons->Fill(posPhoton.X(), posPhoton.Y());
  mapPionMax->Fill(rPosPion.X(), rPosPion.Y());
  mapPionMin->Fill(rPosPionMin.X(), rPosPionMin.Y());

	/	*
  mapPionMaxRev->Fill(rPosPion.X(), rPosPion.Y());
  mapPionMinRev->Fill(rPosPionMin.X(), rPosPionMin.Y());
  * / 


  TVector2 rPosKaonMin, rPosKaonMax;
	populate.checkRangeBelow(posPhoton, getMinCkovKaon(), rPosKaonMin);
	populate.checkRangeBelow(posPhoton, getMaxCkovKaon(), rPosKaonMax);
  mapKaonMax->Fill(rPosKaonMax.X(), rPosKaonMax.Y());
  mapKaonMin->Fill(rPosKaonMin.X(), rPosKaonMin.Y());


  TVector2 rPosProtonMin, rPosProtonMax;
	populate.checkRangeBelow(posPhoton, getMinCkovProton(), rPosProtonMin);
	populate.checkRangeBelow(posPhoton, getMaxCkovProton(), rPosProtonMax);
  mapProtonMax->Fill(rPosProtonMax.X(), rPosProtonMax.Y());
  mapProtonMin->Fill(rPosProtonMin.X(), rPosProtonMin.Y());


  TVector2 rPosProton, temp;
  bool protonBelow = populate.checkRangeBelow(posPhoton, getMinCkovProton(), rPosProton);




      // if not inside ckovMaxPion, continue loop
      if(!pionBelow){
        Printf("!pionBelow");
        continue;
      } if(protonBelow == false && getProtonStatus() == true){
        Printf("protonBelow = 0 && getProtonStatus() = 1");
        continue;
      }
      // global bounds ok, can check candidates 
      else {
        // check Pion
        if(getPionStatus()){ // shouldt really be possible in this case but.. 

            Printf("\n\n ------------ \n\n Region: minProton < ckov < maxProton ok");
            TVector2 rPosPionB;
            bool pionBelow = populate.checkRangeBelow(posPhoton, getMaxCkovPion(), rPosPionB);

            Printf("getPionStatus Ok --> checkRangeAbove getMinCkovPion %.2f ", getMinCkovPion());
            if(populate.checkRangeAbove(posPhoton, getMinCkovPion(), rPosPionB)) {
              mapPion->Fill(xG, yG);
              // add candidate to pions 
              pionCandidates.push_back(std::make_pair(xG, yG));
              Printf("pionCand found ");
              filledBins.push_back(std::make_pair(xG, yG));
              hSignalMap->Fill(xG, yG);
              // range is ok, can check other candidates
              // check if ckov > ckovMaxPion
            }
        }
        if(getKaonStatus()){ // shouldt really be possible in this case but..   		    


          Printf("\n getKaonStatus Ok --> checkRange2 getMinCkovKaon %.2f, getMaxCkovKaon %.2f ", getMaxCkovKaon());

          if (populate.checkRange2(posPhoton, getMaxCkovKaon(), getMinCkovKaon(), temp))
          {	
            Printf("kaonCand found ");			    
            // kaon range ok:
            kaonCandidates.push_back(std::make_pair(xG, yG));
            mapKaon->Fill(xG, yG);
          }
        } 
        if (getProtonStatus()) {
    Printf("\n getProtonStatus Ok --> checkRangeAbove getMaxCkovProton %.2f ", getMinCkovPion());
          if(populate.checkRangeBelow(posPhoton, getMaxCkovProton(), temp)) {
            Printf("protonCand found ");			    
            // proton range ok:
            mapProton->Fill(xG, yG);
            protonCandidates.push_back(std::make_pair(xG, yG));
          }
        } 
      }  // end else / if pionBelow 
      Printf("\n\n");   
			//<end ClosedForm segm>//		
			*/
    }// end else ifTrue
  }  // end for ckovPhotons


	// iterate through photons to check : 



  Printf("	Identified the following  : \n");
  

  Printf("	numBackgroundPhotons %d", numBackgroundPhotons);
  Printf("	numFoundActualCkov %d", numFoundActualCkov);
  Printf("	numActualCkov %d", numActualCkov);
  Printf("	numBackgroundLabeledCkov %d", numBackgroundLabeledCkov);        
	
	
	/*
  if(numFoundActualCkov != numActualCkov) {
			throw std::invalid_argument("wtf value does numActualCkov have?");
  }*/ 

	Printf("=========================================================================");
	Printf("CkovHyps, possible candidates from Momentum | Pion = %d, Kaon = %d, Proton = %d", getPionStatus(), getKaonStatus(), getProtonStatus());


	int cntTemp = 0;
  for(const auto& c : candCombined){
		//Printf("Phot%d : x = %.2f, y = %.2f || statusCand = %d", cntTemp++, c.x, c.y, c.candStatus); 
	}
	Printf("=========================================================================");


  Printf("number of candidates : proton %d, kaon %d, pion %d", protonCandidates.size(), kaonCandidates.size(), pionCandidates.size());
 
  
//for(const auto& pair: filledBins)
//	Printf("CkovTools segment candidates: x%f y%f", pair.first, pair.second);    


// hNoiseMap->SetMarkerColor(kRed);
Printf("CkovTools segment filledBins Size %zu", filledBins.size());



// get impact points of track @ RAD and PC
//const auto trkPC = populate.getPcImp();
/*const auto trkRad = populate.getTrackPos();
TH2F* trkPCMap = new TH2F("trkPCMap ", "trkPCMap; x [cm]; y [cm]",160*20,0.,159.,144*20,0,143);
TH2F* trkRadMap = new TH2F("trkRadMap ", "trkRadMap; x [cm]; y [cm]",160*20,0.,159.,144*20,0,143);
trkRadMap->Fill(trkRad.X(), trkRad.Y());
trkPCMap->Fill(trkPC.X(), trkPC.Y());*/





/*
mapProton->SetMarkerColor(kGreen + 4);
localRefMIP2->Draw("same");  

TCanvas *segm = new TCanvas("semg","semg",800,800);  
segm->Divide(2,2);

segm->cd(1);
mapPhotons->SetMarkerStyle(2); 

mapPhotons->Draw();
localRefMIP2->Draw("same");  
trkPCMap->Draw("same"); trkRadMap->Draw("same"); 


segm->cd(2);
mapPion->SetTitle(Form("mapPion : min%d max%d photons%d", mapPionMin->GetEntries(), mapPionMax->GetEntries(), mapPhotons->GetEntries()));
mapPion->Draw();
mapPhotons->Draw("same");
localRefMIP2->Draw("same");
mapPionMin->Draw("same");
mapPionMax->Draw("same");
trkPCMap->Draw("same"); trkRadMap->Draw("same"); 

segm->cd(3);

mapKaon->SetTitle(Form("mapKaon: min%d max%d photons%d", mapKaonMin->GetEntries(), mapKaonMax->GetEntries(), mapPhotons->GetEntries()));
mapKaon->Draw();
mapPhotons->Draw("same");
mapKaonMin->Draw("same");
mapKaonMax->Draw("same");
trkPCMap->Draw("same"); trkRadMap->Draw("same"); 

segm->cd(4);
mapProton->SetTitle(Form("mapProton : min%d max%d photons%d", mapProtonMin->GetEntries(), mapProtonMax->GetEntries(), mapPhotons->GetEntries()));
mapProton->Draw();

mapPhotons->Draw("same");
mapProtonMin->Draw("same");
mapProtonMax->Draw("same");
trkPCMap->Draw("same"); trkRadMap->Draw("same"); 

gPad->Update();
segm->Show();



TCanvas *thLocal = new TCanvas("thLocal","thLocal",800,800);  
thLocal->cd();
//localBox->Draw();
//localRef->SetMarkerColor(kBlue);


localRefMIP->SetMarkerStyle(3);
localRef->SetMarkerStyle(2);
localRefUnrot->SetMarkerStyle(2);
localRefUnrot->SetMarkerColor(kBlue);
localRefBG->SetMarkerColor(kRed);
localRefBG->SetMarkerStyle(2);
localRefMIPUnrot->SetMarkerColor(kBlue);
localRefMIPUnrot->SetMarkerStyle(3);
localRefMIPUnrot->Draw("same");
 
localRef->Draw();
localRefBG->Draw("same");
localRefUnrot->Draw("same");

Printf("ckovtools segment : localPion->Draw()");
//localPion->Draw("same");    
localRefMIP->Draw("same");// *
tlineUpLocal->Draw();
tlineDownLocal->Draw();
tlineUpLocalR->Draw();
tlineDownLocalR->Draw();//  /
Printf("ckovtools segment : localPion->Draw()");
gPad->Update();
// * / 
*/


	if(print) {
			
		mArrAndMap->drawTotalMapAndMaxRegions();
		mArrAndMap->drawTotalMap();
		mArrAndMap->drawMaxRegions();		  
  
		Printf("momentum %.2f | mass %.2f | thetaP %.2f | phiP %.2f | L %.2f", momentum, mass, thetaP, phiP, L);
		
		
		throw std::invalid_argument("print invoked");
  // drawTotalMap / drawMaxRegions
  }

	return filledBins;
} // end segment


//r0G = rGL + rLP' ; rLP' = Rz(phi)[rLP]
void phiRing2Local(double &xL, double &yL, const double& xMipL, const double& yMipL)
{	  
  TRotation mThetaRot;
  mThetaRot.RotateZ(phiP);
  
  TVector3 mip(xMipL, yMipL, 0);
  TVector3 phiRingPos(xL, yL, 0);
  TVector3 op = mip + mThetaRot*phiRingPos;
  xL = op.Px();
  yL = op.Py();
}
  
//r0G = rGL + rLP' ; rLP' = Rz(phi)[rLP]
// rLP = Rz(-phi) * [r0G-rGL]
void local2PhiRing(double &xL, double &yL, const double& xMipL, const double& yMipL)
{	  
  TRotation mThetaRot;
  mThetaRot.RotateZ(-phiP);
  
  TVector3 pos(xL-xMipL, yL-yMipL, 0);
  TVector3 op = mThetaRot*pos;
  xL = op.Px();
  yL = op.Py();
}

void local2GlobalRef(double& xL, double& yL)
{	  
  xL = xL  + xRad;
  yL = yL  + yRad;	  
}


std::pair<double, double> local2Global(double xL, double yL)
{
  
  const auto x = xL  + xRad;
  const auto y = yL  + yRad;	  
  return {x, y};
}

std::pair<double, double> global2Local(double xG, double yG)
{

  //mTheta.RotateY
  
  const auto x =  xG - xRad;
  const auto y =  yG - yRad;
  return std::make_pair(x, y);
}




// NB! not uses L by simulation here!
// get radius from MIP to a specific phiL, etaC pair
double getR(double etaC, double phiL)
{

	const auto cosPhiL = TMath::Cos(phiL); // 
	const auto sinPhiL = TMath::Cos(phiL); // --||-- 
	
	const auto cosEtaC = TMath::Cos(etaC);
	const auto sinEtaC = TMath::Sin(etaC);

	const auto rwlDeltaR = (rW - L)/(TMath::Sqrt(1-sinEtaC*sinEtaC));
	const auto qwDeltaR = (qW*nF)/(TMath::Sqrt(nQ*nQ-sinEtaC*sinEtaC*nF*nF));
	const auto num = (tGap + tanThetaP*cosPhiL*sinEtaC * (rwlDeltaR + qwDeltaR));


 // ef :error was on this line :
 // 		const auto denum = 1- (tanThetaP*cosPhiL*sinPhiP*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));

	const auto denum = 1 - (tanThetaP*cosPhiL*sinEtaC*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));

	const auto tZ = num/denum;

	const auto Lz = (rW-L) + qW + tZ;

	const auto tGapDeltaR = (tZ*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));

	const auto R = sinEtaC*(rwlDeltaR+qwDeltaR+tGapDeltaR)/cosThetaP;
  //Printf("getR_Lmax : R %f |  wlGap %f qwDeltaR %f tGapDeltaR %f", R, rwlDeltaR,qwDeltaR,tGapDeltaR);
	return R;
}

	// get radius from MIP to a specific phiL, etaC pair
double getR_Lmax(double etaC, double phiL)
{ 


	const auto cosPhiL = TMath::Cos(phiL); // 
	const auto sinPhiL = TMath::Cos(phiL); // --||-- 
	
	const auto cosEtaC = TMath::Cos(etaC);
	const auto sinEtaC = TMath::Sin(etaC);

	const auto rwlDeltaR = (rW - lMax)/(TMath::Sqrt(1-sinEtaC*sinEtaC));
	const auto qwDeltaR = (qW*nF)/(TMath::Sqrt(nQ*nQ-sinEtaC*sinEtaC*nF*nF));
	const auto num = (tGap + tanThetaP*cosPhiL*sinEtaC * (rwlDeltaR + qwDeltaR));

	const auto denum = 1 - (tanThetaP*cosPhiL*sinEtaC*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));
	const auto tZ = num/denum;
	const auto Lz = (rW-lMax) + qW + tZ;

	const auto tGapDeltaR = (tZ*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));

	const auto R = sinEtaC*(rwlDeltaR+qwDeltaR+tGapDeltaR)/cosThetaP;
  //Printf("getR_Lmax : R %f |  wlGap %f qwDeltaR %f tGapDeltaR %f", R, rwlDeltaR,qwDeltaR,tGapDeltaR);
	return R;
}



// populate map in local reference system
// based on one etaC value (i.e., one photon at a time)
const std::pair<double, double> makeCkovPhoton(double phiL, double etaC)
{

  Printf("\nmakeCkovPhoton : phiL %f etaC %f", phiL, etaC);
	const auto cosPhiL = TMath::Cos(phiL);
	const auto sinPhiL = TMath::Sin(phiL);
	
  // if isNan; find etaC that corresponds to rwlDeltaR ^ qwDeltaR NOT nan
  //if(isNan){}
  
	const auto cosEtaC = TMath::Cos(etaC);
	const auto sinEtaC = TMath::Sin(etaC);


	const auto rwlDeltaR = (rW - L)/(TMath::Sqrt(1-sinEtaC*sinEtaC));

	const auto qwDeltaR = (qW*nF)/(TMath::Sqrt(nQ*nQ-sinEtaC*sinEtaC*nF*nF));


	const auto num = (tGap + tanThetaP*cosPhiL*sinEtaC*(rwlDeltaR + qwDeltaR));

	const auto denum = 1 - (tanThetaP*cosPhiL*sinEtaC*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));


  Printf("makeCkovPhoton : rwlDeltaR %f qwDeltaR %f num %f", rwlDeltaR, qwDeltaR, num);

	const auto tZ = num/denum;


	const auto Lz = (rW-L) + qW + tZ;


	const auto tGapDeltaR = (tZ*nF)/(TMath::Sqrt(nG*nG-sinEtaC*sinEtaC*nF*nF));

	const auto T = sinEtaC*(rwlDeltaR+qwDeltaR+tGapDeltaR);


	auto z = - T*tanThetaP * cosPhiL + Lz;
  TRotation mZ1;
  mZ1.RotateZ(phiP);

  TRotation my;
  my.RotateY(thetaP);

  TRotation mZ2;
  mZ2.RotateZ(phiL);

  TVector3 vz(0,0,Lz); TVector3 vx(T,0,0);
  
  TVector3 pos; pos =  vz + mZ2 * vx;
  TVector3 posOut; posOut = mZ1*my*pos*(1/cosThetaP);
  
	const auto x = T*(cosPhiP*cosPhiL - sinPhiP*sinPhiL/cosThetaP) + Lz*tanThetaP*cosPhiP;

	const auto y = T*(sinPhiP*cosPhiL + cosPhiP*sinPhiL/cosThetaP) + Lz*tanThetaP*sinPhiP;


  Printf("makeCkovPhoton : Lz %f tGapDeltaR %f T %f", Lz, tGapDeltaR, T);

	Printf("makeCkovPhoton Rot x %f, y %f z %f || z != %f", posOut.Px(),posOut.Py(),posOut.Pz(), rW-L+tGap+qW);
	Printf("makeCkovPhoton Harcoded x %f, y %f z %f || z != %f \n", x,y,z, rW-L+tGap+qW);
	/*
	const auto xR = T*(cosPhiP*cosPhiL - sinPhiP*sinPhiL/cosThetaP) + Lz*tanThetaP*cosPhiP;

	const auto y = T*(sinPhiP*cosPhiL - cosPhiP*sinPhiL/cosThetaP) + Lz*tanThetaP*sinPhiP;*/
			                            

  //Printf("makeCkovPhoton : x %f, y %f || rot : x %f y %f \n", x, y, posOut.Px(), posOut.Py());
      
  //return std::make_pair(x, y);
	return std::make_pair(posOut.Px(), posOut.Py());
}


double getCkovFromCoords(double xL, double yL, double phiP, double thetaP, float nF, float nQ, float nG, double etaC)
{       


  Printf("ReconG findCkov: cluX %.3f > fPCX %.3f >  cluY %.3f > fPCY %.3f  ", xL, xRad + xMipLocal,yL, yRad + yMipLocal);
  

Printf("ReconG findCkov: cluX %.3f > fPCX %.3f >  cluY %.3f > fPCY %.3f  ", xL, xRad + xMipLocal,yL, yRad + yMipLocal);

  const auto infString = Form("localRefCk xRad%.2f yRad%.2f",xRad,yRad);
  TH2F *localRefCk = new TH2F("localRefCk ", infString,40,-20.,20.,40,-20,20);


  const auto infString2 = Form("localRefCk2 xRad%.2f yRad%.2f",xRad,yRad);
  TH2F *localRefCk2 = new TH2F("localRefCk2 ", infString2,40,-20.,20.,40,-20,20);

  const auto infString3 = Form("localRefCk3 xRad%.2f yRad%.2f",xRad,yRad);
  TH2F *localRefCk3 = new TH2F("localRefCk3 ", infString3,40,-20.,20.,40,-20,20);

  localRefCk->Fill(xL,yL);
  //const auto& coords = local2Global(xL, yL);


  

	double phiF = 0;
	double thetaF1,thetaF2,thetaF=0,thetaLimite;
	float xPi=0,yPi=0,xF=0,yF=0,xF1=0,yF1=0,xF2=0,yF2=0; 
	float ThetaCherenkov, PhiCherenkov, DegPhiCherenkov;

	
	float deltaX = (rW-L+qW+tGap)*tanThetaP*cosPhiP;
	float deltaY = (rW-L+qW+tGap)*tanThetaP*sinPhiP;
		
	xPi = xRad;// - deltaX;
	yPi = yRad;// - deltaY;
  const auto x = xL - deltaX;//coords.first;
  const auto y = yL - deltaY;//coords.second;
	TVector3 v2(x-xPi-L*tanThetaP*cosPhiP, y-yPi-L*tanThetaP*sinPhiP,rW+qW+tGap-L); 

	phiF = v2.Phi();      

	thetaLimite = TMath::ASin(nG/nQ);

	double thetaF0 = TMath::ASin(nQ/nF*TMath::Sin(thetaLimite))-0.00001;

	//  Printf("thetaF0 = %f",thetaF0*TMath::RadToDeg());

	double thetaF01 = TMath::ASin((nF/nQ)*(TMath::Sin(thetaF0)));      

	double thetaF02 = TMath::ASin((nQ/nG)*(TMath::Sin(thetaF01)));

	float x01 = L*tanThetaP*cosPhiP;

	float y01 =  L*tanThetaP*sinPhiP;

	float x02 = (rW - L)*TMath::Tan(thetaF0)*TMath::Cos(phiF)+qW*TMath::Tan(thetaF01)*TMath::Cos(phiF)+tGap*TMath::Tan(thetaF02)*TMath::Cos(phiF);

	float y02 = (rW - L)*TMath::Tan(thetaF0)*TMath::Sin(phiF) + qW*TMath::Tan(thetaF01)*TMath::Sin(phiF) + tGap*TMath::Tan(thetaF02)*TMath::Sin(phiF);  

	float x0 = x01 + x02;
	float y0 = y01 + y02;

	double ThetaMin = 0;
	//  Double_t ThetaMax = 0.75+thetaP;
	double ThetaMax = thetaF0;

	xF = 999;
	yF = 999;

	Int_t nWhile = 0;
	
	TGraph* tCkovReconGraph = new TGraph(50);		
		
	while(TMath::Sqrt((xF-x+xPi)*(xF-x+xPi)+(yF-y+yPi)*(yF-y+yPi))>0.0001)
	{ 
	  nWhile++;

	  thetaF = static_cast<double>( (0.5*(ThetaMax - ThetaMin) + ThetaMin));

	  thetaF1 = TMath::ASin((nF/nQ)*(TMath::Sin(thetaF)));     
	  thetaF2 = TMath::ASin((nQ/nG)*(TMath::Sin(thetaF1)));

	  xF1 = EmissionLenght*tanThetaP*cosPhiP;
	  yF1 =  EmissionLenght*tanThetaP*sinPhiP;

	  xF2 = (rW-L)*TMath::Tan(thetaF)*TMath::Cos(phiF)+qW*TMath::Tan(thetaF1)*TMath::Cos(phiF)+tGap*TMath::Tan(thetaF2)*TMath::Cos(phiF);
	  yF2 = (rW-L)*TMath::Tan(thetaF)*TMath::Sin(phiF) + qW*TMath::Tan(thetaF1)*TMath::Sin(phiF) + tGap*TMath::Tan(thetaF2)*TMath::Sin(phiF);  

	  xF = xF1 + xF2;
	  yF = yF1 + yF2;


	  double X = xF, Y = yF, x1 = xF1, y1 = yF1, x2 = xF2, y2 = yF2;

	  local2PhiRing(X, Y, xMipLocal, yMipLocal);
	  local2PhiRing(x1, y1, xMipLocal, yMipLocal);
	  local2PhiRing(x2, y2, xMipLocal, yMipLocal);

		  localRefCk->Fill(X,Y);
		  localRefCk2->Fill(xF1,yF1);
		  localRefCk3->Fill(xF2,yF2);

		  tCkovReconGraph->SetPoint(nWhile, nWhile, thetaF);

	  //Printf("localRefCk xF %.2f, yF %.2f thetaF %.4f",xF,yF,thetaF );
	  if(TMath::Sqrt((xF-x0)*(xF-x0)+(yF-y0)*(yF-y0))>TMath::Sqrt((x-xPi-x0)*(x-xPi-x0)+(y-yPi-y0)*(y-yPi-y0)))
	  {
	    ThetaMin = thetaF;
	  }
	  else 
	  {
	    ThetaMax = thetaF;   
	  }  

	  if(nWhile>50) break;

	} // while 
	      


	TVector3 vP((TMath::Sin(thetaP))*(TMath::Cos(phiP)),(TMath::Sin(thetaP))*(TMath::Sin(phiP)),(TMath::Cos(thetaP)));
	TVector3 vz(0.,0.,1.);

	TVector3 v1 = vP.Cross(vz);

	TVector3 vF((TMath::Sin(thetaF))*(TMath::Cos(phiF)),(TMath::Sin(thetaF))*(TMath::Sin(phiF)),(TMath::Cos(thetaF)));

	
	
	if(thetaP==0)
	{	      
	  ThetaCherenkov = thetaF;		  
	  PhiCherenkov = phiF;	      
	}	  
	else		
	{
	  vF.Rotate(thetaP,v1);
	  ThetaCherenkov = vF.Theta();
	  PhiCherenkov = vF.Phi();
	}




	DegPhiCherenkov = 180*PhiCherenkov/(TMath::Pi());

              const auto infString4 = Form("bg meth: xRad%.2f yRad%.2f | Actual Ckov : %.3f | Reconstructed Ckov1 = %.3f, Ckov2 = %.3f, avgRecCkov = %.3f",xRad,yRad, etaC, thetaF, vF.Theta(), thetaF/2 + vF.Theta()/2);



	TCanvas* tCkovGraph = new TCanvas("t2","etaC Graph", 1600, 1600);
	tCkovGraph->cd();
	tCkovReconGraph->SetTitle(infString4);
	tCkovReconGraph->Draw("AP");
	

	TCanvas* tcanv = new TCanvas("test", "test", 1600, 1600);
	tcanv->cd();


  localRefCk->SetTitle(infString4);
  localRefCk->Draw();

  localRefCk2->SetMarkerColor(kRed);
  localRefCk2->Draw("same");
  localRefCk3->SetMarkerColor(kBlue);
  localRefCk3->Draw("same");


	if(DegPhiCherenkov<0) DegPhiCherenkov+=360;
  return static_cast<double>(ThetaCherenkov);	
}


// placeholder...
void setArrayMax(double etaTRS, vecArray3& inPutVectorAngle, vecArray2& inPutVectorPos, const size_t kN)
{
  // const size_t kN = inPutVector.size();
  const float lMin = 0.;
  const auto trkPC2 = populatePtr->getPcImp();

  Printf("\n\n setArrayMax() enter --> kN %zu, etaTRS %.2f", kN, etaTRS);
	for(int i = 0; i < kN; i++){

		const auto phiL = Double_t(TMath::TwoPi()*(i+1)/kN);
		TVector3 dirTrs, dirLORS;

		// set TRS values :
		dirTrs.SetMagThetaPhi(1, etaTRS, phiL);

		Printf("\n setArrayMax() dirTrs.SetMagThetaPhi(1, etaTRS = x %.2f, phiL = x %.2f); ", etaTRS, phiL);
		

		double thetaR, phiR; // phiR is value of phi @ estimated R in LORS

		// make this fcn in populate instead?		
		populatePtr->trs2Lors(dirTrs, thetaR, phiR);
		
		Printf("setArrayMax() called trs2Lors on dirTRS : return { thetaR = %.2f, phiR = %.2f}", thetaR, phiR);
		

		dirLORS.SetMagThetaPhi(1, thetaR, phiR);

		Printf("setArrayMax() dirLORS {x %.2f y %.2f z %.2f}", dirLORS.X(), dirLORS.Y(), dirLORS.Z());


		// temp
		// ckovThe, const double& ckovPhi, const double & L
		// this should return the same as max
		const auto t = populatePtr->tracePhot(etaTRS, phiL, lMin);
		Printf("setArrayMax() tracePhot returned TVector2 {x %.2f y %.2f}", t.X(), t.Y());
		// temp
		
		Printf("setArrayMax() called  populatePtr->traceForward(dirLORS (thetaR %.2f, phiR %.2f) lMin =  %.2f", thetaR, phiR, lMin);
		const auto& max = populatePtr->traceForward(dirLORS, lMin);
		
		Printf("setArrayMax() traceForward returned TVector2 {x %.2f y %.2f}", max.X(), max.Y());
		
		// add protection if traceForward returns 0 or -999?
		const auto r = (max - trkPC2).Mod();

 		Printf("setArrayMax() --> i %d, {x %.2f y %.2f} - MIP {x %.2f y %.2f}", i, max.X(), max.Y(), trkPC.X(), trkPC.Y());

		//Printf("setArrayMax2() --> i %d, {x %.2f y %.2f} - MIP {x %.2f y %.2f}", i, max.X(), max.Y(), trkPC2.X(), trkPC2.Y());


		// phiR in [-pi, pi]? set to 0..2pi?
		// inPutVector.emplace_back(std::array<double, 3>{phiL, phiR, r});
		if(phiR < 0) {
			phiR = TMath::TwoPi() + phiR;
    }

		// protections if r > value?
		//if(r > )


		// TODO: if it goes out of map, find intersection with chamber-edges??
		// really jsut check wether TraceForward returned x, y = -999.
 		// to set points out of the map here is fine	


		if((max.Y() == -999) or (max.X() == -999)) {
    	//inPutVector[i] = {0,0,0, 0, 0};
    	//inPutVector.emplace_back(std::array<double, 3>{0, 0, 0});
			// placeholder, find better solution? 
			if(max.Y() == -999) {
				Printf("setArrayMax() max.Y() %.1f == -999", max.Y());
			}
			if(max.X() == -999) {
				Printf("setArrayMax() max.X() %.1f == -999", max.X());
			}
    } else {
    	//inPutVector.emplace_back(std::array<double, 3>{phiL, phiR, r});
    	inPutVectorPos[i] = {max.X(), max.Y()};
    	//inPutVectorAngle[i] = {phiL, phiR, r};
    	inPutVectorAngle.emplace_back(std::array<double, 3>{phiL, phiR, r});
    	Printf("setArrayMax() emplacing element %d : phiL %.2f, phiR %.2f, r %.2f", i, phiL, phiR, r);
    }
    	//inPutVector[i] = {phiL, phiR, r};
    // Printf("setArrayMax() emplacing element %d : phiL %.2f, phiR %.2f, r %.2f", i, phiL, phiR, r);

		/*if(maxPion.X() > 0 && maxPion.X() < 156.0 && maxPion.Y() > 0 && maxPion.Y() < 144) {
			// hMaxPion->Fill(maxPion.X(), maxPion.Y());
			// maxPionVec[i] = std::make_pair(maxPion.X(), maxPion.Y());
			Printf("maxPion loop i = %d, maxSize = %zu", i, maxPionVec.size()); 
		}*/
	}
  	Printf("\n");
	for(const auto& ip : inPutVectorAngle) {
		const auto& phiL_ = ip[0];
		const auto& phiR_ = ip[1]; 
		const auto& r_ = ip[2];  
		Printf("setArrayMax() --> checking inputVector | : phiL %.2f, phiR %.2f, r %.2f", phiL_, phiR_, r_);
  } 
}




// placeholder...
void setArrayMin(double etaTRS, vecArray3& inPutVectorAngle, vecArray2& inPutVectorPos, const size_t kN)
{
  // const size_t kN = inPutVector.size();

  // higher L == > lower radius
  const auto lMax = 1.5;
  const auto trkPC2 = populatePtr->getPcImp();

  Printf("\n\n setArrayMin() enter --> kN %zu, etaTRS %.2f", kN, etaTRS);
	for(int i = 0; i < kN; i++){

		const auto phiL = Double_t(TMath::TwoPi()*(i+1)/kN);
		TVector3 dirTrs, dirLORS;

		// set TRS values :
		dirTrs.SetMagThetaPhi(1, etaTRS, phiL);

		Printf("\n setArrayMin() dirTrs.SetMagThetaPhi(1, etaTRS = x %.2f, phiL = x %.2f); ", etaTRS, phiL);
		

		double thetaR, phiR; // phiR is value of phi @ estimated R in LORS

		// make this fcn in populate instead?		
		populatePtr->trs2Lors(dirTrs, thetaR, phiR);
		
		Printf("setArrayMin() called trs2Lors on dirTRS : return { thetaR = %.2f, phiR = %.2f}", thetaR, phiR);
		

		dirLORS.SetMagThetaPhi(1, thetaR, phiR);

		Printf("setArrayMin() dirLORS {x %.2f y %.2f z %.2f}", dirLORS.X(), dirLORS.Y(), dirLORS.Z());


		// temp
		// ckovThe, const double& ckovPhi, const double & L
		// this should return the same as max
		const auto t = populatePtr->tracePhot(etaTRS, phiL, lMax);
		Printf("setArrayMin() tracePhot returned TVector2 {x %.2f y %.2f}", t.X(), t.Y());
		// temp
		
		Printf("setArrayMin() called  populatePtr->traceForward(dirLORS (thetaR %.2f, phiR %.2f) lMax =  %.2f", thetaR, phiR, lMax);
		const auto& max = populatePtr->traceForward(dirLORS, lMax);
		
		Printf("setArrayMin() traceForward returned TVector2 {x %.2f y %.2f}", max.X(), max.Y());
		
		// add protection if traceForward returns 0 or -999?
		const auto r = (max - trkPC2).Mod();

 		Printf("setArrayMin() --> i %d, {x %.2f y %.2f} - MIP {x %.2f y %.2f}", i, max.X(), max.Y(), trkPC.X(), trkPC.Y());

		Printf("setArrayMin2() --> i %d, {x %.2f y %.2f} - MIP {x %.2f y %.2f}", i, max.X(), max.Y(), trkPC2.X(), trkPC2.Y());


		// phiR in [-pi, pi]? set to 0..2pi?
		// inPutVector.emplace_back(std::array<double, 3>{phiL, phiR, r});
		if(phiR < 0) {
			phiR = TMath::TwoPi() + phiR;
    }

		// protections if r > value?
		//if(r > )


		// TODO: if it goes out of map, find intersection with chamber-edges??
		// really jsut check wether TraceForward returned x, y = -999.
 		// to set points out of the map here is fine	


		if((max.Y() == -999) or (max.X() == -999)) {

			// set element to false, make better solution later..
    	//inPutVector[i] = {0, 0, 0};
			// placeholder, find better solution? s
			if(max.Y() == -999) {
				Printf("setArrayMin() max.Y() %.1f == -999", max.Y());
			}
			if(max.X() == -999) {
				Printf("setArrayMin() max.X() %.1f == -999", max.X());
			}
    } else {
    	inPutVectorPos[i] = {max.X(), max.Y()};
    	//inPutVectorAngle[i] = {phiL, phiR, r};
    	inPutVectorAngle.emplace_back(std::array<double, 3>{phiL, phiR, r});
    	
			//boolArray[i] = true;
    	Printf("setArrayMin() emplacing element %d : phiL %.2f, phiR %.2f, r %.2f", i, phiL, phiR, r);
    }
    	//inPutVector[i] = {phiL, phiR, r};
    // Printf("setArrayMin() emplacing element %d : phiL %.2f, phiR %.2f, r %.2f", i, phiL, phiR, r);

		/*if(maxPion.X() > 0 && maxPion.X() < 156.0 && maxPion.Y() > 0 && maxPion.Y() < 144) {
			// hMaxPion->Fill(maxPion.X(), maxPion.Y());
			// maxPionVec[i] = std::make_pair(maxPion.X(), maxPion.Y());
			Printf("maxPion loop i = %d, maxSize = %zu", i, maxPionVec.size()); 
		}*/
	}
  	Printf("\n");
	for(const auto& ip : inPutVectorAngle) {
		const auto& phiL_ = ip[0];
		const auto& phiR_ = ip[1]; 
		const auto& r_ = ip[2];  
		Printf("setArrayMin() --> checking inputVector | : phiL %.2f, phiR %.2f, r %.2f", phiL_, phiR_, r_);
  } 
}


// map i.e., hMaxPion | vecArray i.e. arrMaxPion
 
void fillMapFromVec(TH2F* map, const vecArray2& vecArr)
{
	for(const auto vE : vecArr) {

		//Printf("fillMapFromVec El %.2f %.2f", vE[0], vE[1]);
		map->Fill(vE[0], vE[1]);
	}
}


// placeholder...
void setArrayMin2(Populate* populate, double etaTRS, vecArray3& inPutVector)
{
  const size_t kN = inPutVector.size();
  const auto lMax = 1.5;
		
	for(int i = 0; i < kN; i++){

		const auto phiL = Double_t(TMath::TwoPi()*(i+1)/kN);
		TVector3 dirTrs, dirLORS;
		dirTrs.SetMagThetaPhi(1, etaTRS, phiL);
		double thetaR, phiR; // phiR is value of phi @ estimated R in LORS
		populatePtr->trs2Lors(dirTrs, thetaR, phiR);
		dirLORS.SetMagThetaPhi(1, thetaR, phiR);

		const auto& max = populatePtr->traceForward(dirLORS, lMax);

		const auto r = (max - trkPC).Mod();
		inPutVector[i] = {phiL, phiR, r};


	}
}


void populateRegions(std::vector<std::pair<double, double>>& vecArr, TH2F* map, const double& eta, const double& l) {	 
	 const int kN = vecArr.size();
	 for(int i = 0; i < vecArr.size(); i++){
		  const auto& value = populatePtr->tracePhot(eta, Double_t(TMath::TwoPi()*(i+1)/kN), l);
		  if(value.X() > 0 && value.X() < 156.0 && value.Y() > 0 && value.Y() < 144) {
		    map->Fill(value.X(), value.Y());
		    vecArr[i] = std::make_pair(value.X(), value.Y());
		  }
	 }
 		
 }


}; // end class CkovTools







