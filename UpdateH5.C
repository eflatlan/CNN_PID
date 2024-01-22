#include "H5Cpp.h"

#include <iostream>

#include <vector>

#include <array>




#include "Alisigma2_.cpp"
#include "Recon.cpp"

/*
n = 100;
auto fPhotCkov = std::unique_ptr<double[]>(new double[n]);
auto fPhotPhi = std::unique_ptr<double[]>(new double[n]);
auto fPhotWei = std::unique_ptr<double[]>(new double[n]);

const float mass_Pion = 0.1396, mass_Kaon = 0.4937, mass_Proton = 0.938; // masses in
std::array<float, 3> masses = {mass_Pion, mass_Kaon, mass_Proton};
// mass_Pion_sq mass_Kaon_sq mass_Proton_sq GeV/c^2
const float mass_Pion_sq = mass_Pion*mass_Pion, mass_Kaon_sq = mass_Kaon*mass_Kaon,  mass_Proton_sq = mass_Proton*mass_Proton;

std::array<float, 3> massesSQ = {mass_Pion_sq, mass_Kaon_sq, mass_Proton_sq};


double calcCkovFromMass(float p, float n, float m)
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


}




*/ 
/*

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double Recon::findRingCkov(int)
{
  // Loops on all Ckov candidates and estimates the best Theta Ckov for a ring formed by those candidates. Also estimates an error for that Theat Ckov
  // collecting errors for all single Ckov candidates thetas. (Assuming they are independent)
  // Arguments: iNclus- total number of clusters in chamber for background estimation
  //    Return: best estimation of track Theta ckov

  Double_t wei = 0.;
  Double_t weightThetaCerenkov = 0.;

  Double_t ckovMin = 9999., ckovMax = 0.;
  Double_t sigma2 = 0; // to collect error squared for this ring

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

      sigma2 += 1. / fParam->sigma2(fTrkDir.Theta(), fTrkDir.Phi(), fPhotCkov[i], fPhotPhi[i]);
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


double findRingCkov()
{
  // Loops on all Ckov candidates and estimates the best Theta Ckov for a ring formed by those candidates. Also estimates an error for that Theat Ckov
  // collecting errors for all single Ckov candidates thetas. (Assuming they are independent)
  // Arguments: iNclus- total number of clusters in chamber for background estimation
  //    Return: best estimation of track Theta ckov

  Double_t wei = 0.;
  Double_t weightThetaCerenkov = 0.;

  Double_t ckovMin = 9999., ckovMax = 0.;
  Double_t sigma2 = 0; // to collect error squared for this ring

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

      // sigma2 += 1. / fParam->sigma2(fTrkDir.Theta(), fTrkDir.Phi(), fPhotCkov[i], fPhotPhi[i]);
    }
  } // candidates loop


	/ *
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
} 
/ *
double houghResponse()
{
  //    fIdxMip = mipId;

  Double_t kThetaMax = 0.75;
  Int_t nChannels = (Int_t)(kThetaMax / fDTheta + 0.5);
  TH1D* phots = new TH1D("Rphot", "phots", nChannels, 0, kThetaMax);
  TH1D* photsw = new TH1D("RphotWeighted", "photsw", nChannels, 0, kThetaMax);
  TH1D* resultw = new TH1D("resultw", "resultw", nChannels, 0, kThetaMax);
  Int_t nBin = (Int_t)(kThetaMax / fDTheta);
  Int_t nCorrBand = (Int_t)(fWindowWidth / (2 * fDTheta));

  for (Int_t i = 0; i < fPhotCnt; i++) { // photon cadidates loop
    Double_t angle = fPhotCkov[i];
    if (angle < 0 || angle > kThetaMax) {
      continue;
    }
    phots->Fill(angle);
    Int_t bin = (Int_t)(0.5 + angle / (fDTheta));
    Double_t weight = 1.;
    if (fIsWEIGHT) {
      Double_t lowerlimit = ((Double_t)bin) * fDTheta - 0.5 * fDTheta;
      Double_t upperlimit = ((Double_t)bin) * fDTheta + 0.5 * fDTheta;
      findRingGeom(lowerlimit);
      Double_t areaLow = getRingArea();
      findRingGeom(upperlimit);
      Double_t areaHigh = getRingArea();
      Double_t diffArea = areaHigh - areaLow;
      if (diffArea > 0) {
        weight = 1. / diffArea;
      }
    }
    photsw->Fill(angle, weight);
    fPhotWei[i] = weight;
  } // photon candidates loop

  for (Int_t i = 1; i <= nBin; i++) {
    Int_t bin1 = i - nCorrBand;
    Int_t bin2 = i + nCorrBand;
    if (bin1 < 1) {
      bin1 = 1;
    }
    if (bin2 > nBin) {
      bin2 = nBin;
    }
    Double_t sumPhots = phots->Integral(bin1, bin2);
    if (sumPhots < 3) {
      continue;
    } // if less then 3 photons don't trust to this ring
    Double_t sumPhotsw = photsw->Integral(bin1, bin2);
    if ((Double_t)((i + 0.5) * fDTheta) > 0.7) {
      continue;
    }
    resultw->Fill((Double_t)((i + 0.5) * fDTheta), sumPhotsw);
  }
  // evaluate the "BEST" theta ckov as the maximum value of histogramm
  Double_t* pVec = resultw->GetArray();
  Int_t locMax = TMath::LocMax(nBin, pVec);
  delete phots;
  delete photsw;
  delete resultw; // Reset and delete objects

  return (Double_t)(locMax * fDTheta + 0.5 * fDTheta); // final most probable track theta ckov

} // HoughResponse()
*/

/*
std::array<float, 3> ckovHypCalc(float p, float n)
{
	std::array<float, 3> ckovAngles = {0, 0, 0};
	
	int i = 0;
	for(const auto& m : massesSQ) {
		ckovAngles[i] = calcCkovFromMass(p, n, m);
		i++;
	}
	return ckovAngles;
}

bool isPhotonHadronCand(double thetaCer, double sigma, double p, double n)
{
	std::array<float, 3> ckovHyps = ckovHypCalc(p, n);
	
	if (std::abs(thetaCer - ckovHyps[0]) < 2*sigma) {
		return true;
	} else if (std::abs(thetaCer - ckovHyps[1]) < 2*sigma) {
		return true;
	} else if (std::abs(thetaCer - ckovHyps[2]) < 2*sigma) {
		return true;
	}
	
	return false;
}*/




using namespace H5;



class Candidate2 {

public:

    double mX, mY, mChi2, mXe, mYe, mPhiCer, mThetaCer, mSigmaRing;

    int mQ, mSize, mCandidateStatus, mCandidateStatusCkov;



    Candidate2(double x = 0, double y = 0, int q = 0, double chi2 = 0, double xe = 0, double ye = 0, int size = 0, int candidateStatus = 0, int candidateStatusCkov = 0, double phiCer = -13., double thetaCer= -13., double sigmaRing= -13.)

        : mX(x), mY(y), mQ(q), mChi2(chi2), mXe(xe), mYe(ye), mSize(size), mCandidateStatus(candidateStatus), mCandidateStatusCkov(candidateStatusCkov), mPhiCer(phiCer), mThetaCer(thetaCer), mSigmaRing(sigmaRing) {}

};



class ParticleInfo2 {

public:

    std::vector<Candidate2> candsCombined;

    double mxRad, myRad, mxMip, myMip, mThetaP, mPhiP, mRefIndex, mMomentum, ckovReconstructed, ckovReconstructedMassHyp;

    int mCluCharge, mCluSize, mTrackPdg;



    ParticleInfo2(float xRad, float yRad, float xMip, float yMip, float th, float ph, float refIndex, int cluCharge, int cluSize, float p, int mcTrackPdg, float ckovRecon)

        : mxRad(xRad), myRad(yRad), mxMip(xMip), myMip(yMip), mThetaP(th), mPhiP(ph), mRefIndex(refIndex), mCluCharge(cluCharge), mCluSize(cluSize), mMomentum(p), mTrackPdg(mcTrackPdg), ckovReconstructed(ckovRecon), ckovReconstructedMassHyp(0) {}



    void updateValues(double sigmaRing, double thetaCer, double phiCer, double ckovReconMass) {

        for (auto& cand : candsCombined) {

            /*cand.mSigmaRing = sigmaRing;

            cand.mThetaCer = thetaCer;

            cand.mPhiCer = phiCer;*/ auto a = 1;

            std::cout<<" cand.mSigmaRing" << cand.mSigmaRing << std::endl; 

        }

        ckovReconstructedMassHyp = ckovReconMass;

    }

};


void createOrUpdateDataset(H5::Group& group, const std::string& name, const std::vector<double>& values) {
    hsize_t dims[1] = {values.size()}; // Define the dimensions of the dataset
    H5::DataSpace dataspace(1, dims); // Create a new dataspace

    // Check if the dataset exists
    if (H5Lexists(group.getId(), name.c_str(), H5P_DEFAULT)) {
        // If the dataset exists, open and update it
        H5::DataSet dataset = group.openDataSet(name);
        dataset.write(values.data(), H5::PredType::NATIVE_DOUBLE);
    } else {
    		std::cout << "Field did not exist, creating " << std::endl;
        // If the dataset does not exist, create it
        H5::DataSet dataset = group.createDataSet(name, H5::PredType::NATIVE_DOUBLE, dataspace);
        dataset.write(values.data(), H5::PredType::NATIVE_DOUBLE);
    }
}


void updateDataset(Group& group, const std::string& name, const std::vector<double>& newValues) {
    DataSet dataset = group.openDataSet(name);
    dataset.write(newValues.data(), PredType::NATIVE_DOUBLE);
}

void updateParams(H5::Group& particleGroup, double momentum, double thetaP, double phiP, double xMip, double yMip, double xRad, double yRad) {


		/*
		auto fPhotFlag = std::unique_ptr<int[]>(new int[n]);
		auto fPhotClusIndex = std::unique_ptr<int[]>(new int[n]);
		*/

		int fPhotCnt = 0;
		

    const DataSet datasetX = particleGroup.openDataSet("x_values");
    const DataSet datasetY = particleGroup.openDataSet("y_values");
    const DataSet datasetQ = particleGroup.openDataSet("q_values");
    const DataSet datasetSize = particleGroup.openDataSet("mSize_values");


    /*DataSet datasetSigmaRing = particleGroup.openDataSet("sigmaRingValues");
    DataSet datasetThetaCer = particleGroup.openDataSet("thetaCerValues");
    DataSet datasetPhiCer = particleGroup.openDataSet("phiCerValues");*/

    
    


    DataSpace spaceX = datasetX.getSpace();
    hsize_t num_elements;
    spaceX.getSimpleExtentDims(&num_elements, NULL);

    std::vector<double> x_values(num_elements), y_values(num_elements), q_values(num_elements), size_values(num_elements), sigma_ring_values(num_elements), phi_cer_values(num_elements),theta_cer_values(num_elements);
    
    
    datasetX.read(x_values.data(), PredType::NATIVE_DOUBLE);
    datasetY.read(y_values.data(), PredType::NATIVE_DOUBLE);
    datasetQ.read(q_values.data(), PredType::NATIVE_DOUBLE);
    datasetSize.read(size_values.data(), PredType::NATIVE_DOUBLE);

    //datasetSigmaRing datasetThetaCer datasetPhiCer



		/*
    datasetSigmaRing.read(sigma_ring_values.data(), PredType::NATIVE_DOUBLE);
    
    datasetThetaCer.read(theta_cer_values.data(), PredType::NATIVE_DOUBLE);
    
    datasetPhiCer.read(phi_cer_values.data(), PredType::NATIVE_DOUBLE);*/


	  Alisigma2_ alisSigma2(1.2904);



		//std::cout<< "ze_t i = 0; i < num_elements; " << std::endl;
    std::vector<SimpleCluster> clusters; 
    for (size_t i = 0; i < num_elements; ++i) {
        clusters.emplace_back(x_values[i], y_values[i], q_values[i], size_values[i]);
        
    }

		Recon reconObj(thetaP,  phiP,  xMip,  yMip,  xRad,  yRad,  clusters);




		//std::cout<< "Correctig phi Theta Sigma" << std::endl;

    for (size_t i = 0; i < num_elements; ++i) {

        
  		   double thetaCer, phiCer, sigma2;

		     //std::cout << "thetaCer"  << thetaCer  << std::endl;


 	 			 //std::cout<< " "<< i <<"/" << num_elements ;

			    phi_cer_values[i] = -13.;
			    theta_cer_values[i] = -13.;        
			    sigma_ring_values[i] = -13.;
		     
		     if (reconObj.findPhotCkov(x_values[i], y_values[i], thetaCer, phiCer)) {
	          double sigma2 = alisSigma2.sigma2(thetaP, phiP, thetaCer, phiCer); // double trkTheta, double trkPhi, double ckovTh, double ckovPh
	          //std::cout << "sigma " << sigma2 << std::endl;
		          
       		          
	          double pionVal = 0.75;
	          	
	          const double refIndexTmp = 1.2904;

					  phi_cer_values[i] = phiCer;
					  theta_cer_values[i] = thetaCer;        
					  sigma_ring_values[i] = sigma2;
						
						
						// bool isCand = isPhotonHadronCand(thetaCer, sigma2, momentum, refIndexTmp);	          
						/*if(isCand) {
						
				      fPhotCkov[fPhotCnt] = thetaCer;                               // actual theta Cerenkov (in TRS)
      				fPhotPhi[fPhotCnt] = phiCer;
      				fPhotClusIndex[fPhotCnt] = iClu; // actual phi   Cerenkov (in TRS): -pi to come back to "unusual" ref system (X,Y,-Z)
      				fPhotCnt++;                      // increment counter of photon candidates      				
      				
						}*/
	          //isPhotonHadronCand(double thetaCer, double sigma, double p, double n)
	          /*if (std::abs(thetaCer-pionVal) < 2*sigma2) {
	          	std::cout << "std::abs(thetaCer-pionVal) < 2*sigma_ring_values[i]" << std::endl;
	          }*/	
	          //std::cout<< "found With sigma" << sigma2 << std::endl;	          
		      } // end if found Recon
		      

			    			     		      

        
        
    }
    
    // datasetSigmaRing datasetThetaCer datasetPhiCer

		createOrUpdateDataset(particleGroup, "sigmaRingValues", sigma_ring_values);
		createOrUpdateDataset(particleGroup, "thetaCerValues", theta_cer_values);
		createOrUpdateDataset(particleGroup, "phiCerValues", phi_cer_values);
								
		//updateDataset(particleGroup, "phiCerValues", phi_cer_values);		
		//std::cout << "exit updateParams" << std::endl;
    // updateParams();
}




// Function to read data from H5 file




// Function to write data to H5 file

void writeH5File(const std::string& filename, const std::vector<ParticleInfo2>& particleInfos) {


		H5File file(filename, H5F_ACC_RDWR); // Open file in read-write mode

    //H5File file(filename, H5F_ACC_TRUNC);

    // Write data from ParticleInfo2 objects to H5 file

    // ...

}

std::vector<ParticleInfo2> readH5File(const std::string& filename) {
    //H5File file(filename, H5F_ACC_RDONLY);



    try {
      H5::H5File file(filename, H5F_ACC_RDWR); // Open file in read-write mode    
      const int num_particles = file.getNumObjs(); // Assuming this gets the number of groups		    
		    
		 

		  for (int i = 0; i < num_particles; ++i) {
		    std::string groupName = "Particle" + std::to_string(i);
		    Group particleGroup = file.openGroup(groupName);


				//std::cout << "Reading group " << groupName << std::endl;
		    //ParticleInfo2 pInfo;
				//std::cout << "Reading "<< i << " of  " << num_particles << std::endl;

		   //ParticleInfo2(float xRad, float yRad, float xMip, float yMip, float th, float ph, float refIndex, int cluCharge, int cluSize, float p, int mcTrackPdg, float ckovRecon)

				double momentum, thetaP, phiP, xMip, yMip, xRad, yRad = 0;

				particleGroup.openAttribute("Momentum").read(PredType::NATIVE_DOUBLE, &momentum);

				particleGroup.openAttribute("ThetaP").read(PredType::NATIVE_DOUBLE, &thetaP);
				particleGroup.openAttribute("PhiP").read(PredType::NATIVE_DOUBLE, &phiP);
				particleGroup.openAttribute("xMip").read(PredType::NATIVE_DOUBLE, &xMip);
				particleGroup.openAttribute("yMip").read(PredType::NATIVE_DOUBLE, &yMip);
				particleGroup.openAttribute("xRad").read(PredType::NATIVE_DOUBLE, &xRad);
				particleGroup.openAttribute("yRad").read(PredType::NATIVE_DOUBLE, &yRad);

				//std::cout<< "xRad " <<xRad << std::endl;
				updateParams(particleGroup, momentum, thetaP,  phiP,  xMip,  yMip,  xRad,  yRad);
				// class SimpleCluster
				// SimpleCluster(double x, double y, double q, double size)

				//std::cout<< "updateParams finished "  << std::endl;

				//Recon reconObj(thetaP,  phiP,  xMip,  yMip,  xRad,  yRad,  clusters);


				//Recon reconObj(double thetaP, double phiP, double xMip, double yMip, double xRad, double yRad, std::vector<SimpleCluster> clustersIn) {


				//Alisigma2_ alisSigma2(1.2904);



				// Process each cluster in parallel

				    

			} // end for num_particles
			
		} catch (const H5::FileIException& e) {
        // Handle file access errors
        std::cerr << "File access error: " << e.getCDetailMsg() << std::endl;
    } catch (const H5::GroupIException& e) {
        // Handle group access errors
        std::cerr << "Group access error: " << e.getCDetailMsg() << std::endl;
    } catch (const H5::Exception& e) {
        // Handle other HDF5 related errors
        std::cerr << "HDF5 error: " << e.getCDetailMsg() << std::endl;
    } catch (const std::exception& e) {
        // Handle other standard exceptions
        std::cerr << "Standard exception: " << e.what() << std::endl;
    } catch (...) {
        // Handle all other types of exceptions
        std::cerr << "Unknown exception caught" << std::endl;
    }
			
			

    //void updateDataset(Group& group, const std::string& name, const     std::vector<double>& newValues)
		
		return std::vector<ParticleInfo2>(); // Placeholder return

}


int UpdateH5( std::string filename = "Particle321_Ckov998_.h5") {

    auto particleInfos = readH5File(filename);




    // writeH5File(filename, particleInfos);



    return 0;

}


