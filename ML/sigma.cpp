
#include <math.h>
#include <cmath>
#include <cassert>



class Sigma {



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
 		double mPhiP, mThetaP;



 		double cosThetaC, tanThetaC, sinThetaC;
 		double mThetaC;

		double SigLoc(double phiC) const
		{
		// Analitical calculation of localization error (due to finite segmentation of PC) on Cerenkov angle for a given Cerenkov photon 
		// created by a given MIP. Fromulae according to CERN-EP-2000-058 
		// Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]
		//            dip and azimuthal angles for MIP taken at the entrance to radiator, [radians]        
		//            MIP beta
		//   Returns: absolute error on Cerenkov angle, [radians]    
			
			double phiDelta = phiC - mPhiP;

			double sint     = sinThetaP;
			double cost     = cosThetaP;
			double sinf     = sinPhiP;
			double cosf     = cosPhiP;
			double sinfd    = TMath::Sin(phiDelta);
			double cosfd    = TMath::Cos(phiDelta);

			
			double alpha =cost-tanThetaC*cosfd*sint;                                                 // formula (11)
			double k = 1.-nF*nF+alpha*alpha/(getTrkBeta()*getTrkBeta());        // formula (after 8 in the text)
			if (k<0) return 1e10;
			double mu =sint*sinf+tanThetaC*(cost*cosfd*sinf+sinfd*cosf);                             // formula (10)
			double e  =sint*cosf+tanThetaC*(cost*cosfd*cosf-sinfd*sinf);                             // formula (9)

			double kk = getTrkBeta()*TMath::Sqrt(k)/(gapThick*alpha);                            // formula (6) and (7)
			double dtdxc = kk*(k*(cosfd*cosf-cost*sinfd*sinf)-(alpha*mu/(getTrkBeta()*getTrkBeta()))*sint*sinfd); // formula (6)           
			double dtdyc = kk*(k*(cosfd*sinf+cost*sinfd*cosf)+(alpha* e/(getTrkBeta()*getTrkBeta()))*sint*sinfd); // formula (7)            pag.4

			double errX = 0.2,errY=0.25;                                                            //end of page 7

			return  errX*errX*dtdxc*dtdxc + errY*errY*dtdyc*dtdyc;
			// return  TMath::Sqrt(errX*errX*dtdxc*dtdxc + errY*errY*dtdyc*dtdyc);
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

	   //Printf("phiC %.2f, mPhiP %.2f", phiC, mPhiP);
			double phiDelta = phiC - mPhiP;

			double sint     = sinThetaP;
			double cost     = cosThetaP;
			double cosfd    = TMath::Cos(phiDelta);

			
			double alpha = cost - tanThetaC*cosfd*sint; 

	    // Printf("tanThetaC %.2f, alpha %.2f, getTrkBeta %.2f", tanThetaC, alpha, getTrkBeta());
			double dtdn = cost*nF*getTrkBeta()*getTrkBeta()/(alpha*tanThetaC);                    // formula (12)
				  

			//return f*dtdn;
			return f*f*dtdn*dtdn;
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

			double sint     = sinThetaP;
			double cost     = cosThetaP;
			double sinf     = sinPhiP;
			double cosfd    = TMath::Cos(phiDelta);


			
			double alpha =cost - tanThetaC*cosfd*sint;                                                  // formula (11)
			
			double k = 1.-nF*nF+alpha*alpha/(getTrkBeta()*getTrkBeta());         // formula (after 8 in the text)
			if (k<0) return 1e10;

			double eTr = 0.5*radThick*getTrkBeta()*TMath::Sqrt(k)/(gapThick*alpha);     // formula (14)
			double lambda = (1.-sint*sinf)*(1.+sint*sinf);                                                  // formula (15)

			double c1 = 1./(1.+ eTr*k/(alpha*alpha*cosThetaC*cosThetaC));                              // formula (13.a)
			double c2 = getTrkBeta()*TMath::Power(k,1.5)*tanThetaC*lambda/(gapThick*alpha*alpha);  // formula (13.b)
			double c3 = (1.+eTr*k*getTrkBeta()*getTrkBeta())/((1+eTr)*alpha*alpha);                                // formula (13.c)
			double c4 = TMath::Sqrt(k)*tanThetaC*(1-lambda)/(gapThick*getTrkBeta());               // formula (13.d)
			double dtdT = c1 * (c2+c3*c4);
			double trErr = radThick/(TMath::Sqrt(12.)*cost);

			// return trErr*dtdT;
			return trErr*trErr*dtdT*dtdT;
		}//SigGeom()

	public : 


	Sigma(double thetaC, double phiP, double thetaP, double _nF) : mThetaC(thetaC), nF(_nF) {

		mThetaC = thetaC;
		mPhiP = phiP;
		mThetaP = thetaP;
		cosThetaP = TMath::Cos(thetaP);
		sinThetaP = TMath::Sin(thetaP);
		tanThetaP = TMath::Tan(thetaP);

		cosPhiP = TMath::Cos(phiP);
		sinPhiP = TMath::Sin(phiP);
		cosThetaC = TMath::Cos(mThetaC);
		tanThetaC = TMath::Tan(mThetaC);
		sinThetaC = TMath::Sin(mThetaC);


	  trkBeta= 1./(TMath::Cos(mThetaC)*nF);
			Printf("mThetaC %.5f | tanThetaC %.5f " , mThetaC,tanThetaC);
		if(trkBeta > 1) trkBeta = 1;       //protection against bad measured mThetaCer  
		if(trkBeta < 0) trkBeta = 0.0001;  //

  }

		double sigma2(double ckovPh) const 
		{
		// Analithical calculation of total error (as a sum of localization, geometrical and chromatic errors) on Cerenkov angle for a given Cerenkov photon 
		// created by a given MIP. Fromulae according to CERN-EP-2000-058 
		// Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]
		//            dip and azimuthal angles for MIP taken at the entrance to radiator, [radians]        
		//            MIP beta
		//   Returns: absolute error on Cerenkov angle, [radians]    
		

			double sigmaLoc2 = SigLoc (ckovPh);
			double sigGeom2 = SigGeom (ckovPh);
			double sigCrom2 = SigCrom (ckovPh);

			Printf("sigmaLoc2 %.5f | sigGeom2 %.5f | sigCrom2 %.5f" , TMath::Sqrt(sigmaLoc2),TMath::Sqrt(sigGeom2), TMath::Sqrt(sigCrom2));
			return  TMath::Sqrt(sigmaThetaP*sigmaThetaP + sigmaLoc2 + sigGeom2 + sigCrom2);

		}


		double sigma2(double ckovPh)
		{
		// Analithical calculation of total error (as a sum of localization, geometrical and chromatic errors) on Cerenkov angle for a given Cerenkov photon 
		// created by a given MIP. Fromulae according to CERN-EP-2000-058 
		// Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]
		//            dip and azimuthal angles for MIP taken at the entrance to radiator, [radians]        
		//            MIP beta
		//   Returns: absolute error on Cerenkov angle, [radians]    
		

			double sigmaLoc2 = SigLoc (ckovPh);
			double sigGeom2 = SigGeom (ckovPh);
			double sigCrom2 = SigCrom (ckovPh);
			


			Printf("sigmaLoc2 %.5f | sigGeom2 %.5f | sigCrom2 %.5f" , TMath::Sqrt(sigmaLoc2),TMath::Sqrt(sigGeom2), TMath::Sqrt(sigCrom2));
			return  TMath::Sqrt(sigmaThetaP*sigmaThetaP + sigmaLoc2+ sigGeom2+ sigCrom2);

		}


}; // end class sigma
