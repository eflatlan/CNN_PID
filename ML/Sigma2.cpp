
#include <math.h>
#include <cmath>
#include <cassert>

class Sigma2 
{
private :

	static constexpr double sigmaThetaP = 0.002; // given in Rad
	double trkBeta; 
	// static constexpr double f = 0.0172*(7.75-5.635)/TMath::Sqrt(24.);
	static constexpr double sq6 = 2.44948974278;
	const  double f =  0.0172*(7.75-5.635)/TMath::Sqrt(24.);
	double getTrkBeta() const { return trkBeta;}
	const double winThick = 0.5, radThick = 1.5; const int gapThick = 8;
	const double gapIdx = 1.0005, winIdx = 1.583; // inIdx = 1.5787;
	double nF;
	double cosThetaP, tanThetaP, sinThetaP;
	double cosPhiP, sinPhiP, tanPhiP;
	double mPhiP, mThetaP,mThetaC; 
	double cosThetaC, tanThetaC, sinThetaC;

public : 


	Sigma2(double thetaC, double phiP, double thetaP, double _nF) : mThetaC(thetaC), nF(_nF) 
	{ 

		mThetaC = thetaC;
		mPhiP = phiP;
		mThetaP = thetaP;
		sinThetaP     = TMath::Sin(mThetaP);
		cosThetaP     = TMath::Cos(mThetaP);


		cosThetaC = TMath::Cos(mThetaC);
		tanThetaC = TMath::Tan(mThetaC);
		trkBeta= 1./(TMath::Cos(mThetaC)*nF);
		Printf("mThetaC %.5f | tanThetaC %.5f " , mThetaC,tanThetaC);
		if(trkBeta > 1) trkBeta = 1;       //protection against bad measured mThetaCer  
		if(trkBeta < 0) trkBeta = 0.0001;  //

		sinThetaP     = TMath::Sin(mThetaP);
		cosThetaP     = TMath::Cos(mThetaP);
		sinPhiP     = TMath::Sin(mPhiP);
		cosPhiP     = TMath::Cos(mPhiP);

	}


	double sigma2(double phiC) const
	{
	// Analithical calculation of total error (as a sum of localization, geometrical and chromatic errors) on Cerenkov angle for a given Cerenkov photon 
	// created by a given MIP. Fromulae according to CERN-EP-2000-058 
	// Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]
	//            dip and azimuthal angles for MIP taken at the entrance to radiator, [radians]        
	//            MIP beta
	//   Returns: absolute error on Cerenkov angle, [radians]    
		


		


		double sigmaLoc = SigLoc (phiC);
		double sigGeom = SigGeom (phiC);
		double sigCrom = SigCrom (phiC);
		


		Printf("sigmaLoc %.5f | sigGeom %.5f | sigCrom %.5f" ,sigmaLoc,sigGeom, sigCrom);
		return  TMath::Sqrt(sigmaThetaP*sigmaThetaP + sigmaLoc*sigmaLoc+ sigGeom*sigGeom+ sigCrom*sigCrom);
	}
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	double SigLoc(double phiC) const
	{
	// Analitical calculation of localization error (due to finite segmentation of PC) on Cerenkov angle for a given Cerenkov photon 
	// created by a given MIP. Fromulae according to CERN-EP-2000-058 
	// Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]
	//            dip and azimuthal angles for MIP taken at the entrance to radiator, [radians]        
	//            MIP beta
	//   Returns: absolute error on Cerenkov angle, [radians]    
		
		double phiDelta = phiC - mPhiP;

		double sinPhiPd    = TMath::Sin(phiDelta);
		double cosPhiPd    = TMath::Cos(phiDelta);
		
		double alpha =cosThetaP-tanThetaC*cosPhiPd*sinThetaP; 

		                                              // formula (11)
		double k = 1.-nF*nF+alpha*alpha/(trkBeta*trkBeta);        // formula (after 8 in the text)
		if (k<0) return 1e10;
		double mu =sinThetaP*sinPhiP+tanThetaC*(cosThetaP*cosPhiPd*sinPhiP+sinPhiPd*cosPhiP);                             // formula (10)
		double e  =sinThetaP*cosPhiP+tanThetaC*(cosThetaP*cosPhiPd*cosPhiP-sinPhiPd*sinPhiP);                             // formula (9)

		double kk = trkBeta*TMath::Sqrt(k)/(gapThick*alpha);                            // formula (6) and (7)
		double dtdxc = kk*(k*(cosPhiPd*cosPhiP-cosThetaP*sinPhiPd*sinPhiP)-(alpha*mu/(trkBeta*trkBeta))*sinThetaP*sinPhiPd); // formula (6)           
		double dtdyc = kk*(k*(cosPhiPd*sinPhiP+cosThetaP*sinPhiPd*cosPhiP)+(alpha* e/(trkBeta*trkBeta))*sinThetaP*sinPhiPd); // formula (7)            pag.4

		double errX = 0.2,errY=0.25;                                                            //end of page 7
		return  TMath::Sqrt(errX*errX*dtdxc*dtdxc + errY*errY*dtdyc*dtdyc);
	}
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	double SigCrom(double phiC) const
	{
	// Analitical calculation of chromatic error (due to lack of knowledge of Cerenkov photon energy) on Cerenkov angle for a given Cerenkov photon 
	// created by a given MIP. Fromulae according to CERN-EP-2000-058 
	// Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]
	//            dip and azimuthal angles for MIP taken at the entrance to radiator, [radians]        
	//            MIP beta
	//   Returns: absolute error on Cerenkov angle, [radians]    
		
		double phiDelta = phiC - mPhiP;

		double sinPhiPd    = TMath::Sin(phiDelta);
		double cosPhiPd    = TMath::Cos(phiDelta);
		
		double alpha =cosThetaP-tanThetaC*cosPhiPd*sinThetaP;                                                 // formula (11)
		double dtdn = cosThetaP*nF*trkBeta*trkBeta/(alpha*tanThetaC);                    // formula (12)
		          
	//  double f = 0.00928*(7.75-5.635)/TMath::Sqrt(12.);
		double f = 0.0172*(7.75-5.635)/TMath::Sqrt(24.);

		return f*dtdn;
	}//SigCrom()
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	double SigGeom(double phiC) const
	{
	// Analitical calculation of geometric error (due to lack of knowledge of creation point in radiator) on Cerenkov angle for a given Cerenkov photon 
	// created by a given MIP. Formulae according to CERN-EP-2000-058 
	// Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]
	//            dip and azimuthal angles for MIP taken at the entrance to radiator, [radians]        
	//            MIP beta
	//   Returns: absolute error on Cerenkov angle, [radians]    

		double phiDelta = phiC - mPhiP;

		double sinPhiPd    = TMath::Sin(phiDelta);
		double cosPhiPd    = TMath::Cos(phiDelta);
		
		double alpha =cosThetaP-tanThetaC*cosPhiPd*sinThetaP;                                                  // formula (11)
		
		double k = 1.-nF*nF+alpha*alpha/(trkBeta*trkBeta);         // formula (after 8 in the text)
		if (k<0) return 1e10;

		double eTr = 0.5*radThick*trkBeta*TMath::Sqrt(k)/(gapThick*alpha);     // formula (14)
		double lambda = (1.-sinThetaP*sinPhiP)*(1.+sinThetaP*sinPhiP);                                                  // formula (15)

		double c1 = 1./(1.+ eTr*k/(alpha*alpha*cosThetaC*cosThetaC));                              // formula (13.a)
		double c2 = trkBeta*TMath::Power(k,1.5)*tanThetaC*lambda/(gapThick*alpha*alpha);  // formula (13.b)
		double c3 = (1.+eTr*k*trkBeta*trkBeta)/((1+eTr)*alpha*alpha);                                // formula (13.c)
		double c4 = TMath::Sqrt(k)*tanThetaC*(1-lambda)/(gapThick*trkBeta);               // formula (13.d)
		double dtdT = c1 * (c2+c3*c4);
		double trErr = radThick/(TMath::Sqrt(12.)*cosThetaP);

		return trErr*dtdT;
	}//SigGeom()

};
