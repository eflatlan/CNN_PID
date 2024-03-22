#include "H5Cpp.h"
#include <iostream>
#include <vector>
#include "TH2F.h"
using vecArray2 = std::vector<std::array<double,2>>;
using namespace H5;

class ParticleUtils2 {
public:




		struct Candidate2 {
		  double mX, mY;
		  int mQ;
		  double mChi2;
		  double mXe, mYe;


    double mPhiCer = -13., mThetaCer = -13.; // initialzie to unvalid values 

    double mSigmaRing = -13.; // resolution of single photon ckov angle


		  int mSize;
		  int mCandidateStatus;
		  int mCandidateStatusCkov;

			Candidate2() : mX(0), mY(0), mQ(0), mChi2(0), mXe(0), mYe(0), mSize(0), mCandidateStatus(0), mCandidateStatusCkov(0) {}
 
 
		  Candidate2(double x = 0, double y = 0, int q = 0, double chi2 = 0, double xe = 0, double ye = 0, int size = 0, int candidateStatus = 0, int candidateStatusCkov = 0, double phiCer = -13., double thetaCer= -13., double sigmaRing= -13.)
		      : mX(x), mY(y), mQ(q), mChi2(chi2), mXe(xe), mYe(ye), mSize(size), mCandidateStatus(candidateStatus) , mCandidateStatusCkov(candidateStatusCkov), mPhiCer(phiCer), mThetaCer(thetaCer), mSigmaRing(sigmaRing)
		  {
		  }
		};

			
		 
		 
    
    struct ParticleInfo2 {

				 
				std::vector<Candidate2> candsCombined;
				double mxRad, myRad, mxMip, myMip, mThetaP, mPhiP, mRefIndex;
				int mCluCharge, mCluSize, mTrackPdg;

				double mMomentum;

				double ckovReconstructed;

				ParticleInfo2(float xRad, float yRad, float xMip, float yMip, float th, float ph, float refIndex, int cluCharge, int cluSize, float p, int mcTrackPdg, float ckovRecon)
				    : mxRad(xRad), myRad(yRad), mxMip(xMip), myMip(yMip), mThetaP(th), mPhiP(ph), mRefIndex(refIndex), mCluCharge(cluCharge), mCluSize(cluSize), mMomentum(p), mTrackPdg(mcTrackPdg), ckovReconstructed(ckovRecon)

				{
				}
        
        void setVector(const std::vector<o2::hmpid::ClusterCandidate>& clusterPerChamber) {


          for(const auto& clu : clusterPerChamber) {

		
	 //mX(x), mY(y), mQ(q), mChi2(chi2), mXe(xe), mYe(ye), mSize(size), mCandidateStatus(candidateStatus) , mCandidateStatusCkov(candidateStatusCkov), mPhiCer(phiCer), mThetaCer(thetaCer), mSigmaRig(sigmaRig)


	candsCombined.emplace_back(clu.mX,clu.mY,clu.mQ,clu.mChi2,clu.mXe,clu.mYe,clu.mSize, clu.mCandidateStatus, clu.mCandidateStatusCkov, clu.mPhiCer, clu.mThetaCer, clu.mSigmaRing);
          }
        }
    };
        

    struct ClusterCandidate {
        double x, y, chi2, q, xe, ye, thetaCer, phiCer;
        int candStatus, candStatusCkov, ch, mSize;
    };

    std::vector<ParticleInfo2> mParticleInfoVector;
        
        
    void fillCandidate(const std::vector<o2::hmpid::ClusterCandidate>& clusterPerChamber, const o2::dataformats::MatchInfoHMP& track, int mcTrackPdg) {
    
        
        float xRad,  yRad,  xPc,  yPc,  th,  ph; // make for these 
        float xMip = track.getMipX(), yMip = track.getMipY(); // and thse 
        track.getHMPIDtrk(xRad,  yRad,  xPc,  yPc,  th,  ph);

        float ckovReconstructed = track.getHMPsignal(); 


        float refIndex = track.getRefIndex(); 
        float cluCharge = track.getMipClusCharge();
        float cluSize = track.getMipClusSize();
		              
				float p = track.getHmpMom(); 
				
				 									
	// ef : this seems to get the PDG correctly?
	int pdg = track.getMipClusEventPDG();
		                      
        ParticleInfo2 particleInfo(xRad, yRad, xMip, yMip, th,  ph, refIndex, cluCharge, cluSize, p, pdg, ckovReconstructed);
        
        particleInfo.setVector(clusterPerChamber);
        
        mParticleInfoVector.push_back(particleInfo);
            	 
    }	  
    	    
    	    
    	       	  
    void writeH5() {	
    
        H5File file("ParticleInfo2.h5", H5F_ACC_TRUNC);



        H5::CompType mtype(sizeof(ClusterCandidate));
        mtype.insertMember("x", HOFFSET(ClusterCandidate, x), H5::PredType::NATIVE_DOUBLE);
        mtype.insertMember("y", HOFFSET(ClusterCandidate, y), H5::PredType::NATIVE_DOUBLE);
        
        mtype.insertMember("chi2", HOFFSET(ClusterCandidate, chi2), H5::PredType::NATIVE_DOUBLE);
        
        mtype.insertMember("q", HOFFSET(ClusterCandidate, q), H5::PredType::NATIVE_INT);
        mtype.insertMember("xe", HOFFSET(ClusterCandidate, xe), H5::PredType::NATIVE_DOUBLE);
        mtype.insertMember("ye", HOFFSET(ClusterCandidate, ye), H5::PredType::NATIVE_DOUBLE);                                
        
	

				mtype.insertMember("phiCer", HOFFSET(ClusterCandidate, phiCer), H5::PredType::NATIVE_DOUBLE);
				mtype.insertMember("thetaCer", HOFFSET(ClusterCandidate, thetaCer), H5::PredType::NATIVE_DOUBLE);     



  	mtype.insertMember("candStatusCkov", HOFFSET(ClusterCandidate, candStatusCkov),  H5::PredType::NATIVE_INT); 
			
	mtype.insertMember("candStatus", HOFFSET(ClusterCandidate, candStatus),  H5::PredType::NATIVE_INT); 

	mtype.insertMember("ch", HOFFSET(ClusterCandidate, ch),  H5::PredType::NATIVE_INT); 
				
				
    for (size_t i = 0; i < mParticleInfoVector.size(); ++i) {
    
    		        
        std::string groupName = "Particle" + std::to_string(i);

        Group particleGroup = file.createGroup(groupName);
        
        
    
        // do padding w zeroes here ? 
        auto& particle = mParticleInfoVector[i];
        std::vector<double> x_values, y_values, chi2_values, xe_values, ye_values, phiCerValues, thetaCerValues, sigmaRingValues;


        std::vector<int> candStatus_values, candStatusCkov_values, ch_values, mSize_values, q_values;  // Add a vector for mSize
        

				//ParticleInfo(floxRad, yRad, xMip, yMip, th,  ph, refIndex, cluCharge, cluSize, p, mcTrackPdg) : mxRad(xRad), myRad(yRad), mxMip(xMip), myMip(yMip), mThetaP(th),  mPhiP(ph), mRefIndex(refIndex), mCluCharge(cluCharge), mCluSize(cluSize), mMomentum(p), mTrackPdg(mcTrackPdg)
	Printf("Writing to H5 : ");
        for (const auto& cand : particle.candsCombined) {

    	    //Printf("cand.mQ ffloat %.2f", cand.mQ );
    	    // Printf("cand.mQ dig %d", cand.mQ );
            x_values.push_back(cand.mX);
            y_values.push_back(cand.mY);
            chi2_values.push_back(cand.mChi2);
            q_values.push_back(cand.mQ);
            xe_values.push_back(cand.mXe);
            ye_values.push_back(cand.mYe);
            candStatus_values.push_back(cand.mCandidateStatus);
  				  candStatusCkov_values.push_back(cand.mCandidateStatusCkov);
  				  sigmaRingValues.push_back(cand.mSigmaRing);


  				  phiCerValues.push_back(cand.mPhiCer);
  				  thetaCerValues.push_back(cand.mThetaCer);
            mSize_values.push_back(cand.mSize);  
        }


				// ef: fix this to hte dimensio that is highest?
        hsize_t dims[1] = { x_values.size() };
        DataSpace dataspace(1, dims);
        
        

        // Here I'm assuming particleGroup is defined and is of type H5::Group
        DataSet dataset = particleGroup.createDataSet("x_values", PredType::NATIVE_DOUBLE, dataspace);
        dataset.write(x_values.data(), PredType::NATIVE_DOUBLE);

        dataset = particleGroup.createDataSet("y_values", PredType::NATIVE_DOUBLE, dataspace);
        dataset.write(y_values.data(), PredType::NATIVE_DOUBLE);

        dataset = particleGroup.createDataSet("q_values", PredType::NATIVE_INT, dataspace);
        dataset.write(q_values.data(), PredType::NATIVE_INT);

        dataset = particleGroup.createDataSet("mSize_values", PredType::NATIVE_INT, dataspace);
        dataset.write(mSize_values.data(), PredType::NATIVE_INT);

        dataset = particleGroup.createDataSet("thetaCerValues", PredType::NATIVE_DOUBLE, dataspace);
        dataset.write(thetaCerValues.data(), PredType::NATIVE_DOUBLE);

        dataset = particleGroup.createDataSet("sigmaRingValues", PredType::NATIVE_DOUBLE, dataspace);
        dataset.write(sigmaRingValues.data(), PredType::NATIVE_DOUBLE);

        dataset = particleGroup.createDataSet("phiCerValues", PredType::NATIVE_DOUBLE, dataspace);
        dataset.write(phiCerValues.data(), PredType::NATIVE_DOUBLE);

        dataset = particleGroup.createDataSet("chi2_values", PredType::NATIVE_DOUBLE, dataspace);
        dataset.write(chi2_values.data(), PredType::NATIVE_DOUBLE);

        dataset = particleGroup.createDataSet("xe_values", PredType::NATIVE_DOUBLE, dataspace);
        dataset.write(xe_values.data(), PredType::NATIVE_DOUBLE);

        dataset = particleGroup.createDataSet("ye_values", PredType::NATIVE_DOUBLE, dataspace);
        dataset.write(ye_values.data(), PredType::NATIVE_DOUBLE);

        dataset = particleGroup.createDataSet("candStatus_values", PredType::NATIVE_INT, dataspace);
        dataset.write(candStatus_values.data(), PredType::NATIVE_INT);






        dataset = particleGroup.createDataSet("candStatusCkov_values", PredType::NATIVE_INT, dataspace);
        dataset.write(candStatusCkov_values.data(), PredType::NATIVE_INT);


        //dataset = particleGroup.createDataSet("ch_values", PredType::NATIVE_INT, dataspace);
        //dataset.write(ch_values.data(), PredType::NATIVE_INT);

        // Add a dataset for mSize


						  


				std::vector<Candidate2> candsCombined;


				DataSpace attr_dataspace = DataSpace(H5S_SCALAR);

				// Implementing for all the required attributes
				Attribute attribute = particleGroup.createAttribute("xRad", PredType::NATIVE_DOUBLE, attr_dataspace);
				attribute.write(PredType::NATIVE_DOUBLE, &particle.mxRad);



				attribute = particleGroup.createAttribute("yRad", PredType::NATIVE_DOUBLE, attr_dataspace);
				attribute.write(PredType::NATIVE_DOUBLE, &particle.myRad);

				attribute = particleGroup.createAttribute("xMip", PredType::NATIVE_DOUBLE, attr_dataspace);
				attribute.write(PredType::NATIVE_DOUBLE, &particle.mxMip);

				attribute = particleGroup.createAttribute("yMip", PredType::NATIVE_DOUBLE, attr_dataspace);
				attribute.write(PredType::NATIVE_DOUBLE, &particle.myMip);

				attribute = particleGroup.createAttribute("ThetaP", PredType::NATIVE_DOUBLE, attr_dataspace);
				attribute.write(PredType::NATIVE_DOUBLE, &particle.mThetaP);

				attribute = particleGroup.createAttribute("PhiP", PredType::NATIVE_DOUBLE, attr_dataspace);
				attribute.write(PredType::NATIVE_DOUBLE, &particle.mPhiP);
				
				
				attribute = particleGroup.createAttribute("RefractiveIndex", PredType::NATIVE_DOUBLE, attr_dataspace);
				attribute.write(PredType::NATIVE_DOUBLE, &particle.mRefIndex);

				attribute = particleGroup.createAttribute("CluCharge", PredType::NATIVE_INT, attr_dataspace);
				attribute.write(PredType::NATIVE_INT, &particle.mCluCharge);

				attribute = particleGroup.createAttribute("CluSize", PredType::NATIVE_INT, attr_dataspace);
				attribute.write(PredType::NATIVE_INT, &particle.mCluSize);

				attribute = particleGroup.createAttribute("Momentum", PredType::NATIVE_DOUBLE, attr_dataspace);
				attribute.write(PredType::NATIVE_DOUBLE, &particle.mMomentum);

				attribute = particleGroup.createAttribute("TrackPdg", PredType::NATIVE_INT, attr_dataspace);
				attribute.write(PredType::NATIVE_INT, &particle.mTrackPdg);

				attribute = particleGroup.createAttribute("ckovReconstructed", PredType::NATIVE_DOUBLE, attr_dataspace);
				attribute.write(PredType::NATIVE_DOUBLE, &particle.ckovReconstructed);

				
// After writing all attributes
Printf("Wrote particle with PDG %d : Ckov Reconstructed : %.3f xRad %.3f, yRad %.3f, xMip %.3f, yMip %.3f, ThetaP %.3f, PhiP %.3f, RefractiveIndex %.3f, CluCharge %d, CluSize %d, Momentum %.3f ",
       particle.mTrackPdg, particle.ckovReconstructed, particle.mxRad, particle.myRad, particle.mxMip, particle.myMip, particle.mThetaP, particle.mPhiP, particle.mRefIndex, particle.mCluCharge, particle.mCluSize, particle.mMomentum);



    	}
				
        
    }
        
};
