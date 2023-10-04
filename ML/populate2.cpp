
#ifndef TEST_POPULATE2
#define TEST_POPULATE2
// populate.cpp
#include <TVector2.h>
#include <TVector3.h>
#include <TRotation.h>
#include <TRandom.h>
#include <vector>


//#include "Math/GenVector/Rotation3D.h"
#include "Math/GenVector/Rotation3D.h"
#include "Math/GenVector/Cartesian3D.h"
#include "Math/GenVector/DisplacementVector3D.h"

#include "Math/GenVector/RotationYfwd.h"
#include "Math/GenVector/RotationZfwd.h"


class Populate2 {

   // using array = std::array;
	using vecArray4 = std::vector<std::array<double,4>>;

   using vecArray3 = std::vector<std::array<double,3>>;
   using Polar3D = ROOT::Math::Polar3D<double>;
   using Rotation3D = ROOT::Math::Rotation3D;
   using RotationY = ROOT::Math::RotationY;
   using RotationZ = ROOT::Math::RotationZ;


private:



    TVector2 fTrkPos; // track pos in LORS at RAD // xRad, yRad
    TVector3 fTrkDir; // track position in LORS at RAD // setMagThetaPhi(1, thetaP, phiP)
    TVector2 fPc; // track pos at PC

    TVector2 fMipPos; // MIP-pos at PC

    double nF, getRefIdx;	     // refIdnex of freon

    const double winThick = 0.5, radThick = 1.5; const int gapThick = 8;
    const double gapIdx = 1.0005, winIdx = 1.5787;

    const double nQ = winIdx, nG = gapIdx;

    double thetaP, cosPhiP, sinPhiP, tanThetaP;
    double deltaX, deltaY;

    double phiP, sinThetaP, cosThetaP;


    TVector2 fTrkPos2D;


		// unc in HMPI in mm (should prob be 2*sqrt(2) or 1*sqrt(2))
		static constexpr float uncSquared = .2;

public:


    TH2F* getLimMin() {
        return limMin;
    }

    TH2F* getLimMax() {
        return limMax;
    }



    Populate2(TVector2 trkPos, TVector3 trkDir, double _nF) : fTrkPos(trkPos),  fTrkDir(trkDir), nF(_nF) 
    {

 			/*
	    limMin = new TH2F("minPos","minPos", 1600, 0, 159, 1440, 0, 143);

			limMax = new TH2F("maxPos", "maxPos", 1600, 0, 159, 1440, 0, 143);

			

				TVector3 fTrkDir; // track direction in LORS at RAD

				TVector2 fTrkPos; // track positon in LORS at RAD   // XY mag
				TVector2 fMipPos; // mip positon for a given trackf // XY
				TVector2 fPc;     // track position at PC           // XY

			*/			
			
	    fTrkPos2D.Set(trkPos.X(), trkPos.Y()); 


      //fPc.setX()
      //Printf("init Populate class");

      
		  thetaP = trkDir.Theta();
		  tanThetaP = TMath::Tan(thetaP);
		  cosThetaP = TMath::Cos(thetaP);		
		  sinThetaP = TMath::Sin(thetaP);

		  getRefIdx = nF;  

		  phiP = trkDir.Phi();
			cosPhiP = TMath::Cos(phiP);		
			sinPhiP = TMath::Sin(phiP);	
		


			// xRa = xPC - deltaX


      // should be radThick/2 here if assuming half em-length */
      deltaX = (radThick/2 + winThick + gapThick) * tanThetaP * cosPhiP;
      deltaY = (radThick/2 + winThick + gapThick) * tanThetaP * sinPhiP;		
			

      // NB! TODO: here PC impact point based on L = rW/2!!
      //setPcImp(fTrkPos.X() + deltaX, fTrkPos.Y() + deltaY);
      
      setPcImp(fTrkPos.X() + deltaX, fTrkPos.Y() + deltaY);
      //Printf("Track pos at RAD : x %.3f y %.3f ", trkPos.X(), trkPos.Y());
      //Printf("nF %.3f : x %.3f y %.3f ", nF, trkPos.X(), trkPos.Y());
      
      //Printf("Track pos at PC : x %.3f y %.3f ", fPc.X(), fPc.Y());

      //Printf("Track dir at RAD : theta %.3f phi %.3f ", trkDir.Theta(), trkDir.Phi());      
				
    }

		TH2F* limMin = nullptr; 
		TH2F* limMax = nullptr; 



  // Trace a single Ckov photon from emission point somewhere in radiator up to photocathode taking into account ref indexes of materials it travereses
  // Arguments: ckovThe,ckovPhi- photon ckov angles in TRS, [rad]
  //   Returns: distance between photon point on PC and track projection
    // L here is "simulated" within 0..1.5 range
    TVector2 tracePhot(const double& ckovThe, const double& ckovPhi, const double & L) const {
    		Printf("populate2 : tracePhot()");
        double theta, phi;
        TVector3 dirTRS, dirLORS;
        dirTRS.SetMagThetaPhi(1, ckovThe, ckovPhi); // photon in TRS
        trs2Lors(dirTRS, theta, phi);
        dirLORS.SetMagThetaPhi(1, theta, phi); // photon in LORS
        return traceForward(dirLORS, L);          // now foward tracing
    }

    void propagate(const TVector3& dir, TVector3& pos, double z) const {
        static TVector3 nrm(0, 0, 1);
        TVector3 pnt(0, 0, z);

        TVector3 diff = pnt - pos;
        double sint = (nrm * diff) / (nrm * dir);
        pos += sint * dir;
    }

    void refract(TVector3& dir, const double& n1, const double& n2) const {
        double sinref = (n1 / n2) * TMath::Sin(dir.Theta());
        if (TMath::Abs(sinref) > 1.) {
            dir.SetXYZ(-999, -999, -999);
        } else {
            dir.SetTheta(TMath::ASin(sinref));
        }
    }



		// TODO :this value is not set? (getRefIdx)
    TVector2 traceForward(TVector3& dirCkov, const double& L) const {
			
				// getrefIdx set in ctor?
	      //  auto getRefIdx = nF;
		
        TVector2 pos(-999, -999);
        double thetaCer = dirCkov.Theta();
        
        // if called to make the boundaries from hyp: since were 
        // using n x std_dev to make the boundaries, this value can be exceeded, 
        // therefore instead "Threshold" the value
        
        
        // ef :TODO: ask Giacomo, is this valid?
        // the ckov photon can be above this value, just not the track cherenkov value?
        
        // its the dispersion of chromacity that makes the normal distribution around the track Ckov value?
        if (thetaCer > TMath::ASin(1. / getRefIdx)) {
        		LOGP(debug, "populate2: traceForward() INVOKED thetaCer {} > TMath::ASin(1. / getRefIdx{}))", getRefIdx, thetaCer);
        		
        		
        		// instead, set the threshold ? : 
        		//thetaCer = TMath::ASin(1. / getRefIdx);
            return pos; // ef: TODO this was changed to the above
        }
        


				// change radThick to other value to change L
 				// auto radThick' = (radThick - L);
        // double zRad = - radThick' - 0.5 * winThick;
        double zRad = - (radThick - L) - 0.5 * winThick; 
  
				// TODO: which value should be changed??

        TVector3 posCkov(fTrkPos.X(), fTrkPos.Y(), zRad);
        propagate(dirCkov, posCkov, -0.5 * winThick); // TODO giacomo spm :er ikke dette også feil!
        refract(dirCkov, getRefIdx, winIdx);
        propagate(dirCkov, posCkov, 0.5 * winThick);
        refract(dirCkov, winIdx, gapIdx);
        propagate(dirCkov, posCkov, 0.5 * winThick + gapThick);
        pos.Set(posCkov.X(), posCkov.Y());


				if(pos.X() == -999) {
					Printf("	traceForward() : pos.X() = -999!");
				}	if(pos.Y() == -999) {
					Printf("	traceForward() : pos.Y() = -999!");
				}

        return pos;
    }

    void lors2Trs(const TVector3& dirCkov, double& thetaCer, double& phiCer) const {
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

    void trs2Lors(const TVector3& dirCkov, double& thetaCer, double& phiCer) const {

        TRotation mtheta;
        mtheta.RotateY(fTrkDir.Theta());

        TRotation mphi;
        mphi.RotateZ(fTrkDir.Phi());

        TRotation mrot = mphi * mtheta;


        TVector3 dirCkovLORS;
        dirCkovLORS = mrot * dirCkov;
	

        //Polar3D dirCkovLORS2;
        //dirCkovLORS2 = mrot * dirCkov;

        phiCer = dirCkovLORS.Phi();     // actual value of the phi of the photon
        thetaCer = dirCkovLORS.Theta(); // actual value of thetaCerenkov of the photon

	/*Printf("trs2Lors");
	Printf("	old : phi %.3f, theta %.3f", phiCer, thetaCer);
	Printf("	new : phi %.3f, theta %.3f", dirCkovLORS2.Phi(), dirCkovLORS2.Theta());*/
    }


    // getter and setter functions

    void setTrackPos(double x, double y)
    {
			fTrkPos.SetX(x); 
 			fTrkPos.SetY(y); 
    }

    void setTrackPos(const TVector2& fTrkIn)
    {
			fTrkPos = fTrkIn;
    }
    
    TVector2 getTrackPos() const
    {
			return fTrkPos;
    }


    void setPcImp(double x, double y)
    {
			fPc.SetX(x); 
 			fPc.SetY(y); 
    }

    void setPcImp(const TVector2& fPcImpIn)
    {
			fPc = fPcImpIn;
    }
    
    TVector2 getPcImp() const
    {
			return fPc;
    }


    TVector3 getTrkDir() const
    {
			return fTrkDir;
    }


	bool checkOver(const TVector2& posPhoton, const double& rPhoton, const double& phiPhoton, vecArray4& vec, const double& etaCkov, const char* hadronType) {
		// vec : contains phiL, phi, R w etaMin/etaMax for hadron species
		
		const auto sizeVec = vec.size();	


		const double lMin = 0.;
    bool decisionTaken = false, condition = false;

    //const std::array<double, 3>	fst = vec[index];	

    const size_t kN = vec.size();
    int initValue;
		
		auto phiC = phiPhoton - phiP;
		
		if(phiC > TMath::TwoPi()){
			phiC = phiC - TMath::TwoPi();
    }
		

		Printf("\n\n enter checkOver(%s) phiPhoton %.4f, phiP %.4f", hadronType, phiPhoton, phiP);


		//LOGP(info, "Some radius checks checkOver rPhoton {} : [sizeVec - 2]{} [sizeVec - 1]{} [0]{} [1]{}", rPhoton, vec[sizeVec - 2][3],vec[sizeVec - 1][3],vec[0][3],vec[1][3]);
		
		if (phiC < TMath::Pi()/2)
			initValue = 0;
		else if (phiC < TMath::Pi())
			initValue = static_cast<int>(kN/4);
		else if (phiC < (3/2)*TMath::Pi())
			initValue = static_cast<int>(kN/2);
		else if (phiC < TMath::TwoPi())
			initValue = static_cast<int>((3/2)*kN);
		else {
			initValue = 0;
			Printf("phiC %.2f", phiC);
			throw std::invalid_argument("wtf value does phiC have?");
		} 
		// TODO: set iv back to reasonable"" value:)
 		initValue = 0;

		int iCnt = 0;
		auto phi2 = vec[initValue][2]; // 


	  //Printf(" Enter checkOver(%s) \n initValue = %d \n  vec[initValue-1][1] = %.3f || phi2 = vec[initValue][1] =  %.3f |  vec[initValue+1][1] = %.3f;", hadronType,initValue, vec[sizeVec-1][1], vec[initValue][1], vec[initValue+1][1]);


	//Printf(" phi::: %.3f | %.3f | %.3f | %.3f | %.3f | ;", vec[sizeVec-1][1],vec[0][1], vec[1][1], vec[2][1], vec[3][1]);

	//Printf(" r::: %.3f | %.3f | %.3f | %.3f | %.3f | ;", vec[sizeVec-1][2],vec[0][2], vec[1][2], vec[2][2], vec[3][2]);

		// set increment opposite way if phi1 > phiPhoton
		int inc = 1;



    // vec[initValue-1][1] = 3.873 || phi2 = vec[initValue][1] =  3.908 

	  // 1st check if it is between 1st and last: ?
		if(vec[sizeVec - 1][2] < phiPhoton && phiPhoton < vec[0][2]) {
		}


		/* iterate backwards?
		if (phi2 > phiPhoton)
			inc = -1;		*/ 

		// accesing the correct phi


		// TODO: change back to being 0?
		// double prev = 0.;// vec[initValue-1][1];

		double prev = vec[initValue][2];
		while(/*phi2 < phiPhoton*/ true) {



      // what if  initValue + iCnt + 1  > vec size?
			if(initValue + iCnt + 1 >= kN) { // do vec size instead? 
				Printf("initValue %d + iCnt %d + 1 >= kN%d", initValue, iCnt, kN);
				Printf("\n prev %.6f  phi2 %.6f < phiPhoton  %.6f \n", prev, phi2,  phiPhoton);

				//throw std::invalid_argument("Problem in while-loop???");
				return true; // false
				// We have open contour, or contour goes out of map --> The particle should be l
				// labelled as being a candidate			

				// return false;
				// ef : returner, og tegn kontur?
				// ef : dette betyr at det er en open contour?
				// throw std::invalid_argument("Problem in while-loop???");
			}

			phi2 = vec[initValue + iCnt + 1][2]; // waas not +1_
			//Printf("	while(phi12< phiPhoton) {  || phi2 %.4f, prev %.4f phiPhoton %.4f, iCnt %d, initValue + iCnt = %d", phi2, prev, phiPhoton, iCnt, initValue + iCnt);



			// these 2 ifs deals with "going around"
			// this means phi2 went "around" and went from 6.xx to 0.xx
			// then phi2  = 0.xx, prev = 6.xx and phiPhoton is in this range
			if(prev > phi2 && prev < phiPhoton && phi2 < phiPhoton && phi2 != 0) {
				Printf("\n prev %.6f > phi2 %.6f && phi2 %.6f < phiPhoton  %.3fexit! \n", prev, phi2, phi2, phiPhoton);
				break;
			} 
			if(prev > phi2 && phi2 > phiPhoton && phi2 != 0) {
				Printf("\n prev %.6f > phi2 %.6f && phi2 %.6f > phiPhoton %.6f, exit! \n", prev, phi2,phi2, phiPhoton);
				break;
			}


			/*denne stmmer vel ikke?:
			if(prev > phiPhoton && phiPhoton > phi2) {
				Printf("\n prev %.6f > phi2 %.6f, exit! \n", prev, phi2);
				break;
			} */  

			// not break if prev = 0, phi2 > phiPhoton at 1st element
			if(phi2 > phiPhoton && prev < phiPhoton && prev != 0) {
				Printf("\n phi2 %.6f > phiPhoton %.6f && prev %.6f < phiPhoton %.6f: exit!\n",  phi2, phiPhoton, prev, phiPhoton);
			 	break; 	
			}
			iCnt += inc;
			prev = phi2;
		}



		// hvis det går rundt her? 
    // vil det bli break Segmenetation fault siden vi prøver å finne element som er out of range

		// TODO: this has to be changed if inc = -1?

		double phiL1, phiL2, phi1, r1, r2;

		int index = iCnt - 1;


		// hvis ikke har gått rundt:
		if(initValue + iCnt+1 < sizeVec) {
			//double phiL1 = vec.at(index).at(0);
			phiL1 = vec[initValue + iCnt][0];
			phiL2 = vec[initValue + iCnt+1][0];
			

			// was index
			phi1 = vec[initValue + iCnt][2];
			phi2 = vec[initValue + iCnt + 1][2]; // iCnt

			// get radiuses
			r1 = vec[initValue + iCnt][3];
			r2 = vec[initValue + iCnt+1][3];

		} else {
			//double phiL1 = vec.at(index).at(0);
			phiL1 = vec[sizeVec - 1][0];
			phiL2 = vec[0][0];
			

			// was index
			phi1 = vec[sizeVec - 1][2];
			phi2 = vec[0][2]; // iCnt

			// get radiuses
			r1 = vec[sizeVec - 1][3];
			r2 = vec[0][3];
		} 

		if(phi2 == 0) {
			
			Printf("	checkOver() -->  phi2 == 0, setting back to 1st elem in arr");
			Printf("	Old vals : phi2 %.2f, phiL2 %.2f, r2 %.2f", phi2, phiL2, r2);			
			phiL2 = vec[0][0];
			phi2 = vec[0][2];
			r2 = vec[0][3];
			Printf("	New vals : phi2 %.2f, phiL2 %.2f, r2 %.2f", phi2, phiL2, r2);			
		}

		// correct indexes are found 

		Printf("	checkOver(%s) -->  kN %d | initValue %d | phiPhoton %.6f, phiP %.3f, phiC %.3f, phi2 %.6f, phi1 %.3f,", hadronType, kN, initValue, phiPhoton, phiP, phiC, phi2, phi1);

		//Printf("phiPhoton %.2f| vec[initValue + iCnt - 1] %.2f,  vec[initValue + iCnt] %.2f, vec[initValue + iCnt + 1] %.2f", phiPhoton, vec[initValue + iCnt-1][1],  vec[initValue + iCnt][1], vec[initValue + iCnt+ 1][1]);
		
				

		// eta is min/max ckovHyps +- 3 std-dev
		


		// TODO: change also logic behind?

	  double rMin, rMax, phiMin, phiMax, phiLmin, phiLmax;
		if(r2 > r1) {
			rMin = r1;
			rMax = r2;

			phiMin = phi1;
			phiMax = phi2;

			phiLmin = phiL1;
			phiLmax = phiL2;
		}	

		else {

			rMax = r1;
			rMin = r2;

			phiMax = phi1;
			phiMin = phi2;

			phiLmax = phiL1;
			phiLmin = phiL2;
		}


		{
				TVector2 newPos1, newPos2; 

        float phiLmaxCp = phiLmax;
        float phiLminCp = phiLmin;

				float rMaxCp = getR(newPos1, etaCkov, phiLmaxCp, lMin);
				float rMinCp = getR(newPos2, etaCkov, phiLminCp, lMin);
				Printf("	newPos1 = %.2f,  %.2f ||newPos2 = %.2f, %.2f", newPos1.X(), newPos1.Y(), newPos2.X(), newPos2.Y());
				limMax->Fill(newPos1.X(), newPos1.Y());
				limMin->Fill(newPos2.X(), newPos2.Y());
				Printf("	rPhoton = %.2f, rMax = %.2f, rMin = %.2f || from getR :  rMax = %.2f, rMin = %.2f", rPhoton, rMax, rMin, rMaxCp, rMinCp);

		}

		Printf("	rPhoton = %.2f, rMax = %.2f, rMin = %.2f", rPhoton, rMax, rMin);

		// check if under a radius 
		// 
		while(!decisionTaken) {

			Printf("	| rPhoton = %.2f, rMax = %.2f, rMin = %.2f",  rPhoton, rMax, rMin);

			Printf("	phiMax = %.4f, phiLmax = %.4f, phiMin = %.4f, phiLmin = %.4f, etaCkov %.4f, phiPhoton %.4f",phiMax, phiLmax, phiMin, phiLmin, etaCkov, phiPhoton);



			/*
			if((rPhoton < rMax + uncSquared)) {
				// stop iterating, condition is ok!
				decisionTaken = true;
				condition = true;
				break; // just tbs
			}
			else if((rPhoton >= rMax + uncSquared)) {
				// stop iterating, condition is false
				decisionTaken = true;
				condition = false;
				break; // just tbs
			} */				
			
			if((rPhoton < rMin + uncSquared)) {

			   Printf("rPhoton %.2f < rMin + uncSquared %.2f",  rPhoton, rMin + uncSquared);
				// stop iterating, condition is ok!
				decisionTaken = true;
				condition = true;
				break; // just tbs
			}
			else if((rPhoton > rMax + uncSquared)) {
			   Printf("rPhoton %.2f > rMax + uncSquared %.2f",  rPhoton, rMax + uncSquared);
				// stop iterating, condition is false
				decisionTaken = true;
				condition = false;
				break; // just tbs
			} else if((rPhoton > rMin + (uncSquared) && rPhoton < rMax + (uncSquared))) {

				// iterate by splitting r1, r2 -> phi1 phi2
			
				// find new value for eta1, eta2; (passed by ref)
				// phiL1, phiL2 also passed by ref

				splitPhi(phiMax, phiLmax, phiMin, phiLmin, etaCkov, phiPhoton); // 

				TVector2 newPos1, newPos2;
				rMax = getR(newPos1, etaCkov, phiLmax, lMin);
				rMin = getR(newPos2, etaCkov, phiLmin, lMin); // TODO: should it be lMin here?
				//limMax->Fill(newPos1.X(), newPos1.Y());
				//limMin->Fill(newPos2.X(), newPos2.Y());

					Printf("limMax num BinEntries = %f ", limMax->GetEntries());
					Printf("limMin num BinEntries = %f ", limMin->GetEntries());
				Printf("	checkOver(%s) | rMin = %.2f , rMax = %.2f,  \n",hadronType,  rMin, rMax);
				Printf("	limMin %.2f %.2f,  \n", newPos1.X(), newPos1.Y());
			} else {
				

				if(iCnt + 1 < sizeVec) {
					Printf("\n iCnt + 1 < sizeVec");
					Printf(" ??? r ==> -1 : %.3f 0 : %.3f 1: %.3f  " , vec[iCnt -1][2], vec[iCnt][2], vec[iCnt+1][2]);
					Printf("??? phi ==> -1 : %.3f 0 : %.3f 1: %.3f  " , vec[iCnt -1][1], vec[iCnt][1], vec[iCnt+1][1]);
					Printf("??? phiL ==> -1 : %.3f 0 : %.3f 1: %.3f  " , vec[iCnt -1][0], vec[iCnt][0], vec[iCnt+1][0]);
				} else {
					Printf("\n iCnt + 1 == sizeVec");
					Printf(" ??? r ==> -1 : %.3f 0 : %.3f 1: %.3f  " , vec[iCnt -1][2], vec[iCnt][2], vec[0][2]);
					Printf("??? phi ==> -1 : %.3f 0 : %.3f 1: %.3f  " , vec[iCnt -1][1], vec[iCnt][1], vec[0][1]);
					Printf("??? phiL ==> -1 : %.3f 0 : %.3f 1: %.3f  " , vec[iCnt -1][0], vec[iCnt][0], vec[0][0]);
				}
				Printf("r2 %.3f r1 %.3f ", r2, r1);
				Printf("phi2 %.6f phi1 %.3f ", phi2, phi1);
				Printf("phiL2 %.3f phiL2 %.3f ", phiL2, phiL1);

				throw std::invalid_argument("decisionTaken???");
			}
		}
	

		Printf("exit checkOver(%s) with condition = %d|  rMin = %.2f , rMax = %.2f, rPhoton = %.2f | L used for lim : %.2f, etaC = %.2f", hadronType, condition, rMin, rMax, rPhoton, lMin, etaCkov);

Printf("	phiMin = %.4f <  phiPhoton %.4f <  phiMax = %.4f, ",phiMin, phiPhoton, phiMax);

		if(rMin > rMax ) {Printf("NB!!!!!!!! exit checkOver(%s) :: rMin %.2f > rMax %.2f ", hadronType, rMin, rMax);}


		// the condition to be evaluated {false/true}
		return condition;
	}




  
  /*  vecArray4& vec {phiL, phiR, phiPC, r}
			// phiL : photon_phi i TRS system 
    	// phiR : photon_phi i LORS system
	  	// phiPC : (photon - MIP).Phi(); --> vector phiL "on paper"
    	 
			const auto r = (max - trkPC2).Mod(); // trkPC2 : track impact @ PC
			// from alinot_paattrec : this is MIP2photon distance?
			// const auto r = (max - trkPC2).Mod(); // trkPC2 : track impact @ PC
	*/

	bool checkUnder(const TVector2& posPhoton, const double& rPhoton, const double& phiPhoton, vecArray4& vec, const double& etaCkov, const char* hadronType, bool& shutDownOnOpen) {
		// vec : contains phiL, phi, R w etaMin/etaMax for hadron species
		
		
		const auto sizeVec = vec.size();
				
				

		
		//const double lMin = 0., 
		const double lMax = 1.5;
    bool decisionTaken = false, condition = false;

    //const std::array<double, 3>	fst = vec[index];	

    const size_t kN = vec.size();
    int initValue;
		
		auto phiC = phiPhoton - phiP; // phiC : value to be "approximated"		
																	// phiPhoton : phiValue of photon in LORS
																	// phiP track phi value ()
		
		if(phiC > TMath::TwoPi()){
			phiC = phiC - TMath::TwoPi();
    }
		

		Printf("\n\n enter checkUnder(%s) phiPhoton %.4f, phiP %.4f",  hadronType, phiPhoton, phiP);

		//LOGP(info, "Some radius checks checkUnder rPhoton {} : [sizeVec - 2]{} [sizeVec - 1]{} [0]{} [1]{}", rPhoton,vec[sizeVec - 2][3],vec[sizeVec - 1][3],vec[0][3],vec[1][3]);
		


		if (phiC < TMath::Pi()/2)
			initValue = 0;
		else if (phiC < TMath::Pi())
			initValue = static_cast<int>(kN/4);
		else if (phiC < (3/2)*TMath::Pi())
			initValue = static_cast<int>(kN/2);
		else if (phiC < TMath::TwoPi())
			initValue = static_cast<int>((3/2)*kN);
		else {
			initValue = 0;
			Printf("phiC %.2f", phiC);
			throw std::invalid_argument("wtf value does phiC have?");
		} 
		// TODO: set iv back to reasonable"" value:)
 		initValue = 0;

		int iCnt = 0;
		auto phi2 = vec[initValue][2]; // phiPC : (photon - MIP).Phi(); --> vector



	  //Printf(" Enter checkUnder(%s) \n initValue = %d \n  vec[initValue-1][1] = %.3f || phi2 = vec[initValue][1] =  %.3f |  vec[initValue+1][1] = %.3f;", hadronType,initValue, vec[sizeVec-1][1], vec[initValue][1], vec[initValue+1][1]);

	//Printf(" phi::: %.3f | %.3f | %.3f | %.3f | %.3f | ;", vec[sizeVec-1][1],vec[0][1], vec[1][1], vec[2][1], vec[3][1]);

	//Printf(" r::: %.3f | %.3f | %.3f | %.3f | %.3f | ;", vec[sizeVec-1][2],vec[0][2], vec[1][2], vec[2][2], vec[3][2]);

		// set increment opposite way if phi1 > phiPhoton
		int inc = 1;



		/* iterate backwards?
		if (phi2 > phiPhoton)
			inc = -1;		*/ 

		// accesing the correct phi


		// TODO: change back to being 0?
		// double prev = 0.;// vec[initValue-1][1];

		double prev = vec[initValue][2];
		Printf("\n Enter while loop with prev %.6f  \n", prev);
		while(/*phi2 < phiPhoton*/ true) {


      // what if  initValue + iCnt + 1  > vec size?
			if(initValue + iCnt + 1 >= kN) { // do vec size instead? 
				shutDownOnOpen = true;
				Printf("initValue %d + iCnt %d + 1 >= kN%d", initValue, iCnt, kN);
				return false; // true
				// We have open contour, or contour goes out of map --> The particle should be l
				// labelled as being a candidate			

				// return false;
				// ef : returner, og tegn kontur?
				// ef : dette betyr at det er en open contour?
				// throw std::invalid_argument("Problem in while-loop???");
			}
			phi2 = vec[initValue + iCnt + 1][2]; // waas not +1_
			//Printf("	while(phi12< phiPhoton) {  || phi2 %.4f, prev %.4f phiPhoton %.4f, iCnt %d, initValue + iCnt = %d", phi2, prev, phiPhoton, iCnt, initValue + iCnt);

			//Printf("\n prev %.6f  phi2 %.6f phiPhoton %.6f \n", prev, phi2, phiPhoton);

			// these 2 ifs deals with "going around"
			// this means phi2 went "around" and went from 6.xx to 0.xx
			// then phi2  = 0.xx, prev = 6.xx and phiPhoton is in this range
			if(prev > phi2 && prev < phiPhoton && phi2 < phiPhoton && prev != 0 ) {
				Printf("\n prev %.6f > phi2 %.6f && phi2 %.6f < phiPhoton  %.3fexit! \n", prev, phi2,phi2, phiPhoton);
				break;
			} 
			if(prev > phi2 && phi2 > phiPhoton && prev != 0) {
				Printf("\n prev %.6f > phi2 %.6f && phi2 %.6f > phiPhoton %.6f, exit! \n", prev, phi2,phi2, phiPhoton);
				break;
			}


			/*denne stmmer vel ikke?:
			if(prev > phiPhoton && phiPhoton > phi2) {
				Printf("\n prev %.6f > phi2 %.6f, exit! \n", prev, phi2);
				break;
			} */  

			// not break if prev = 0, phi2 > phiPhoton at 1st element
			if(phi2 > phiPhoton && prev < phiPhoton && prev != 0) {
				Printf("\n phi2 %.6f > phiPhoton %.6f && prev %.6f < phiPhoton %.6f: exit!\n",  phi2, phiPhoton, prev, phiPhoton);
			 	break; 	
			}
			iCnt += inc;
			prev = phi2;
		}

		// TODO: this has to be changed if inc = -1?
		int index = iCnt - 1;

		double phiL1, phiL2, phi1, r1, r2;


		// hvis ikke har gått rundt:
		if(initValue + iCnt+ 1 < sizeVec) {

		Printf("initValue %d + iCnt %d + 1 < sizeVec%d", initValue, iCnt, sizeVec);
			//double phiL1 = vec.at(index).at(0);
			phiL1 = vec[initValue + iCnt][0];
			phiL2 = vec[initValue + iCnt+1][0];
			

			// was index
			phi1 = vec[initValue + iCnt][2];
			phi2 = vec[initValue + iCnt + 1][2]; // iCnt

 
			// get radiuses
			r1 = vec[initValue + iCnt][3];
			r2 = vec[initValue + iCnt+1][3];



		} else {

			//double phiL1 = vec.at(index).at(0);
			phiL1 = vec[initValue + iCnt][0];
			phiL2 = vec[0][0];
			

			// was index
			phi1 = vec[initValue + iCnt][2];
			phi2 = vec[0][2]; // iCnt

			// get radiuses
			r1 = vec[initValue + iCnt][3];
			r2 = vec[0][3];
		} 


		if(phi2 == 0) {
			
			Printf("	checkUnder(%s) -->  phi2 == 0, setting back to 1st elem in arr");
			Printf("	Old vals : phi2 %.2f, phiL2 %.2f, r2 %.2f", phi2, phiL2, r2);			
			phiL2 = vec[0][0];
			phi2 = vec[0][2];
			r2 = vec[0][3];
			Printf("	New vals : phi2 %.2f, phiL2 %.2f, r2 %.2f", phi2, phiL2, r2);			
		}


		Printf("	checkUnder(%s) -->  kN %d | initValue %d | phiPhoton %.6f, phiP %.3f, phiC %.3f, phi2 %.6f, phi1 %.3f,", hadronType, kN, initValue, phiPhoton, phiP, phiC, phi2, phi1);

		//Printf("phiPhoton %.2f| vec[initValue + iCnt - 1] %.2f,  vec[initValue + iCnt] %.2f, vec[initValue + iCnt + 1] %.2f", phiPhoton, vec[initValue + iCnt-1][1],  vec[initValue + iCnt][1], vec[initValue + iCnt+ 1][1]);
		
				

		// eta is min/max ckovHyps +- 3 std-dev
		


		// TODO: change also logic behind?

	  double rMin, rMax, phiMin, phiMax, phiLmin, phiLmax;
		if(r2 > r1) {
			rMin = r1;
			rMax = r2;
			phiMax = phi2;
			phiMin = phi1;
			phiLmax = phiL2;
			phiLmin = phiL1;
		}	

		else {
			rMin = r2;
			rMax = r1;
			phiMax = phi1;
			phiMin = phi2;
			phiLmax = phiL1;
			phiLmin = phiL2;
		}


		{
				TVector2 newPos1, newPos2; 

        float phiLmaxCp = phiLmax;
        float phiLminCp = phiLmin;

				float rMaxCp = getR(newPos1, etaCkov, phiLmaxCp, lMax);
				float rMinCp = getR(newPos2, etaCkov, phiLminCp, lMax);
				limMax->Fill(newPos1.X(), newPos1.Y());
				limMin->Fill(newPos2.X(), newPos2.Y());

				Printf("	newPos1 = %.2f,  %.2f ||newPos2 = %.2f, %.2f", newPos1.X(), newPos1.Y(), newPos2.X(), newPos2.Y());

				Printf("	rPhoton = %.2f, rMax = %.2f, rMin = %.2f || from getR :  rMax = %.2f, rMin = %.2f", rPhoton, rMax, rMin, rMaxCp, rMinCp);

		}

		Printf("	rPhoton = %.2f, rMax = %.2f, rMin = %.2f", rPhoton, rMax, rMin);

		// check if outside of radius


		// rA = rPhoton > rMax - uncSq
		// rB = rPhoton < rMin - uncSq
		// rC rPhoton > rMin && rPohoton < rMax
		
		// rPhoton > r 


		while(!decisionTaken) {
			Printf("	rPhoton = %.2f, rMax = %.2f, rMin = %.2f", rPhoton, rMax, rMin);
			/*
			if((rPhoton > rMin - uncSquared)) {
				// stop iterating, condition is true
				decisionTaken = true;
				condition = true;
			} else if((rPhoton  <= rMin - uncSquared)) {
				// stop iterating, condition is false!
				decisionTaken = true;
				condition = false;
			} */
			 
			if((rPhoton  >  rMax - uncSquared)) { // was rMin!!
			   Printf("rPhoton %.2f >  rMax - uncSquared %.2f",  rPhoton, rMax - uncSquared);
				// stop iterating, condition is true
				decisionTaken = true;
				condition = true;
			} else if((rPhoton < rMin - uncSquared)) {
			   Printf("rPhoton %.2f  < rMin - uncSquared %.2f",  rPhoton, rMin - uncSquared);
				// stop iterating, condition is false!
				decisionTaken = true;
				condition = false;
			} else if((rPhoton + (uncSquared)> rMin && rPhoton + (uncSquared)< rMax)) {
				// iterate by splitting r1, r2 -> phi1 phi2



				splitPhi(phiMax, phiLmax, phiMin, phiLmin, etaCkov, phiPhoton); // 
				TVector2 newPos1, newPos2;
				rMax = getR(newPos1, etaCkov, phiLmax, lMax);
				rMin = getR(newPos2, etaCkov, phiLmin, lMax); 
				limMax->Fill(newPos1.X(), newPos1.Y());
				limMin->Fill(newPos2.X(), newPos2.Y());

				Printf("	checkUnder(%s) | rMin = %.2f , rMax = %.2f,  ",hadronType,  rMin, rMax);
				Printf("	limMin %.2f %.2f,  \n", newPos1.X(), newPos1.Y());
			} else {

				if(iCnt + 1 < sizeVec) {
					Printf("\n iCnt + 1 < sizeVec");
					Printf(" ??? r ==> -1 : %.3f 0 : %.3f 1: %.3f  " , vec[iCnt -1][2], vec[iCnt][2], vec[iCnt+1][2]);
					Printf("??? phi ==> -1 : %.3f 0 : %.3f 1: %.3f  " , vec[iCnt -1][1], vec[iCnt][1], vec[iCnt+1][1]);
					Printf("??? phiL ==> -1 : %.3f 0 : %.3f 1: %.3f  " , vec[iCnt -1][0], vec[iCnt][0], vec[iCnt+1][0]);
				} else {
					Printf("\n iCnt + 1 == sizeVec");
					Printf(" ??? r ==> -1 : %.3f 0 : %.3f 1: %.3f  " , vec[iCnt -1][2], vec[iCnt][2], vec[0][2]);
					Printf("??? phi ==> -1 : %.3f 0 : %.3f 1: %.3f  " , vec[iCnt -1][1], vec[iCnt][1], vec[0][1]);
					Printf("??? phiL ==> -1 : %.3f 0 : %.3f 1: %.3f  " , vec[iCnt -1][0], vec[iCnt][0], vec[0][0]);
				}
				Printf("r2 %.3f r1 %.3f ", r2, r1);
				Printf("phi2 %.6f phi1 %.3f ", phi2, phi1);
				Printf("phiL2 %.3f phiL2 %.3f ", phiL2, phiL1);

				//Printf("	checkUnder(%s) | rMin = %.2f , rMax = %.2f,  ",hadronType,  rMin,
				throw std::invalid_argument("decisionTaken???");
			}
		}

		Printf("exit checkUnder(%s) with condition = %d|  rMin = %.2f , rMax = %.2f, rPhoton = %.2f | L used for lim : %.2f, etaC = %.2f", hadronType, condition, rMin, rMax, rPhoton, lMax, etaCkov);


Printf("	phiMin = %.4f <  phiPhoton %.4f <  phiMax = %.4f \n, ",phiMin, phiPhoton, phiMax);

		if(rMin > rMax ) {Printf("NB!!!!!!!! exit checkUnder(%s) :: rMin %.2f > rMax %.2f ", hadronType, rMin, rMax);}


		// the condition to be evaluated {false/true}
		return condition;
	}

		





	double getR(TVector2& rPosLORS, const double& etaTRS, const double& phiTRS, const double& L)
	{

	
		//TVector3 dirCkov;
		//dirCkov.SetMagThetaPhi(1, etaTRS, phiTRS);
		rPosLORS = tracePhot(etaTRS, phiTRS, L); // pos of estimated rPos at LORS PC
Printf("\n enter  getR()  rPosLORS {x %.2f y %.2f} - MIP {x %.2f y %.2f}",  rPosLORS.X(), rPosLORS.Y(), (getPcImp()).X(), (getPcImp()).Y());

		const auto R = (rPosLORS-getPcImp()).Mod(); 
   Printf("\n exit : getR() : etaTRS %.2f | phiTRS %.2f | L %.2f | R %.2f ", etaTRS, phiTRS, L, R);
		return R;
	}
		
	void splitPhi(double& phiMax, double& phiLmax, double& phiMin, double& phiLmin, const double& etaCkov, const double& phiPhoton)
	{


		// phiPhoton == phi in LORS

		// split the TRS value 
		const auto& phiNewL = (phiLmax + phiLmin)/2.0;	// phi in TRS


		// Get phi in LORS 
		TVector3 dirNewTrs;
		dirNewTrs.SetMagThetaPhi(1, etaCkov, phiNewL); // phi in TRS



    double theta, phi;
		trs2Lors(dirNewTrs, theta, phi); 
		// double theta, phi;


		// nei!!
		// auto phiNew = dirNewTrs.Phi();

		auto phiNew = phi;
		if (phiNew < 0) {
			phiNew = TMath::TwoPi() + phiNew;
		}

    Printf("\n enter : splitPhi() \n phiMax %.2f , phiMin  %.2f ", phiMax, phiMin); 

   Printf("	splitPhi() etaCkov %.2f | phiNew %.2f | phiP %.2f | thetaP %.2f | phiPhoton %.2f ", etaCkov, phiNew, phiP, thetaP, phiPhoton);

    Printf("	splitPhi : phiLmax %.2f | phiLmin %.2f ", phiLmax, phiLmin);
		if((phiPhoton) < TMath::Pi() + phiP){
			
			if((phiNew-phiP) > (phiPhoton-phiP)) {
				phiMin = phiNew; 
				phiLmin = phiNewL;		 
				Printf("	splitPhi : if((phiPhoton) < TMath::Pi() + phiP --> (phiNew-phiP) > (phiPhoton-phiP)");
			} else {
				phiMax = phiNew;
				phiLmax = phiNewL;
Printf("	splitPhi : if((phiPhoton) < TMath::Pi() + phiP --> else()");
			}
		} 
		else {
			if((phiNew-phiP) < (phiPhoton-phiP)) {
        Printf("	splitPhi : if((phiPhoton) > TMath::Pi() + phiP --> (phiNew-phiP) < (phiPhoton-phiP)");
				phiMin = phiNew; 
				phiLmin = phiNewL; 
			} else {
				Printf("	splitPhi : if((phiPhoton) > TMath::Pi() + phiP --> else");
				phiMax = phiNew;
				phiLmax = phiNewL; 
			}
		}		

		Printf("	splitPhi : ");
    Printf("\n splitPhi Exit : phiLmax %.2f | phiLmin %.2f ", phiLmax, phiLmin);

	} // end splitPhi()


		// track pos in pc 2 rad
    // void pc2rad(const TVector& pc) */ 
    ClassDefNV(Populate2, 1);
};

#endif
