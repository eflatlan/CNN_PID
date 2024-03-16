#include "H5Cpp.h"

#include <iostream>

#include <vector>

#include <array>

#include "Alisigma2_.cpp"

#include "Recon.cpp"

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

  Candidate2(double x = 0, double y = 0, int q = 0, double chi2 = 0,
             double xe = 0, double ye = 0, int size = 0,
             int candidateStatus = 0, int candidateStatusCkov = 0,
             double phiCer = -13., double thetaCer = -13.,
             double sigmaRing = -13.)

      : mX(x), mY(y), mQ(q), mChi2(chi2), mXe(xe), mYe(ye), mSize(size),
        mCandidateStatus(candidateStatus),
        mCandidateStatusCkov(candidateStatusCkov), mPhiCer(phiCer),
        mThetaCer(thetaCer), mSigmaRing(sigmaRing) {}
};

class ParticleInfo2 {

public:
  std::vector<Candidate2> candsCombined;

  double mxRad, myRad, mxMip, myMip, mThetaP, mPhiP, mRefIndex, mMomentum,
      ckovReconstructed, ckovReconstructedMassHyp;

  int mCluCharge, mCluSize, mTrackPdg;

  ParticleInfo2(float xRad, float yRad, float xMip, float yMip, float th,
                float ph, float refIndex, int cluCharge, int cluSize, float p,
                int mcTrackPdg, float ckovRecon)

      : mxRad(xRad), myRad(yRad), mxMip(xMip), myMip(yMip), mThetaP(th),
        mPhiP(ph), mRefIndex(refIndex), mCluCharge(cluCharge),
        mCluSize(cluSize), mMomentum(p), mTrackPdg(mcTrackPdg),
        ckovReconstructed(ckovRecon), ckovReconstructedMassHyp(0) {}

  void updateValues(double sigmaRing, double thetaCer, double phiCer,
                    double ckovReconMass) {

    for (auto &cand : candsCombined) {

      /*cand.mSigmaRing = sigmaRing;



      cand.mThetaCer = thetaCer;



      cand.mPhiCer = phiCer;*/
      auto a = 1;

      // std::cout<<" cand.mSigmaRing" << cand.mSigmaRing << std::endl;
    }

    ckovReconstructedMassHyp = ckovReconMass;
  }
};

void readAndPrintDataset(const std::string &filename,
                         const std::string &datasetName) {
  try {
    H5::H5File file(filename, H5F_ACC_RDONLY); // Open file in read-only mode
    H5::DataSet dataset = file.openDataSet(datasetName);
    H5::DataSpace dataspace = dataset.getSpace();

    hsize_t dims_out[1];
    dataspace.getSimpleExtentDims(dims_out, NULL);
    std::vector<double> data(dims_out[0]);

    dataset.read(data.data(), H5::PredType::NATIVE_DOUBLE);

    // Print the data
    /*
    for (const auto& value : data) {
        //std::cout << value << " ";
    }
    //std::cout << std::endl; */
  } catch (const H5::Exception &e) {
    std::cerr << "HDF5 error: " << e.getCDetailMsg() << std::endl;
  }
}

void createOrUpdateDataset(H5::Group &group, const std::string &name,
                           const std::vector<double> &values) {

  hsize_t dims[1] = {values.size()}; // Define the dimensions of the dataset

  H5::DataSpace dataspace(1, dims); // Create a new dataspace

  /*//std::cout << std::fixed << std::setprecision(4);

  //std::cout << "{";
  int maxValues = std::min(static_cast<int>(values.size()), 50);
  for (int i = 0; i < maxValues; i++) {

      if (values[i] > 1000) continue;

    //std::cout << values[i];
    if (i < maxValues - 1) {
      //std::cout << ", ";
    }
  }
  //std::cout << "}" << std::endl; */

  // Check if the dataset exists

  if (H5Lexists(group.getId(), name.c_str(), H5P_DEFAULT)) {

    ////std::cout << "Field did exist | OLD : " << std::endl;

    H5::DataSet dataset = group.openDataSet(name);

    DataSpace spaceX = dataset.getSpace();

    hsize_t num_elements;

    spaceX.getSimpleExtentDims(&num_elements, NULL);

    std::vector<double> x_values(num_elements);

    dataset.read(x_values.data(), PredType::NATIVE_DOUBLE);

    /*
    //std::cout << std::fixed << std::setprecision(4);

    //std::cout << "{";
    int maxValues2 = std::min(static_cast<int>(x_values.size()), 50);
    for (int i = 0; i < maxValues2; i++) {

            if (x_values[i] > 1000) continue;

      //std::cout << x_values[i];
      if (i < maxValues2 - 1) {
        //std::cout << ", ";
      }
    }
    //std::cout << "}" << std::endl; */

    // If the dataset exists, open and update it

    dataset.write(values.data(), H5::PredType::NATIVE_DOUBLE);

  } else {

    ////std::cout << "Field did not exist, creating " << std::endl;

    // If the dataset does not exist, create it

    H5::DataSet dataset =
        group.createDataSet(name, H5::PredType::NATIVE_DOUBLE, dataspace);

    dataset.write(values.data(), H5::PredType::NATIVE_DOUBLE);
  }
}

void updateDataset(Group &group, const std::string &name,
                   const std::vector<double> &newValues) {

  DataSet dataset = group.openDataSet(name);

  dataset.write(newValues.data(), PredType::NATIVE_DOUBLE);
}

void updateParams(H5::Group &particleGroup, double momentum, double thetaP,
                  double phiP, double xMip, double yMip, double xRad,
                  double yRad, const std::string fileName) {

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

  std::vector<double> x_values(num_elements), y_values(num_elements),
      q_values(num_elements), size_values(num_elements),
      sigma_ring_values(num_elements), phi_cer_values(num_elements),
      theta_cer_values(num_elements);

  datasetX.read(x_values.data(), PredType::NATIVE_DOUBLE);

  datasetY.read(y_values.data(), PredType::NATIVE_DOUBLE);

  datasetQ.read(q_values.data(), PredType::NATIVE_DOUBLE);

  datasetSize.read(size_values.data(), PredType::NATIVE_DOUBLE);

  // datasetSigmaRing datasetThetaCer datasetPhiCer
  /*
  datasetSigmaRing.read(sigma_ring_values.data(), PredType::NATIVE_DOUBLE);
  datasetThetaCer.read(theta_cer_values.data(), PredType::NATIVE_DOUBLE);
  datasetPhiCer.read(phi_cer_values.data(), PredType::NATIVE_DOUBLE);*/

  const double refIndexTmp = 1.2904;

  Alisigma2_ alisSigma2(refIndexTmp);

  ////std::cout<< "ze_t i = 0; i < num_elements; " << std::endl;

  std::vector<SimpleCluster> clusters;

  for (size_t i = 0; i < num_elements; ++i) {
    clusters.emplace_back(x_values[i], y_values[i], q_values[i],
                          size_values[i]);
  }

  // Recon(double thetaP, double phiP, double xMip, double yMip, double xRad,
  // double yRad, const std::vector<SimpleCluster>& clustersIn) {

  Recon reconObj(thetaP, phiP, xMip, yMip, xRad, yRad, clusters);

  for (size_t i = 0; i < num_elements; ++i) {

    double thetaCer, phiCer, sigma2;

    phi_cer_values[i] = -13.;
    theta_cer_values[i] = -13.;
    sigma_ring_values[i] = -13.;

    // bool findPhotCkov(double cluX, double cluY, double& thetaCer, double&
    // phiCer)

    if (reconObj.findPhotCkov(x_values[i], y_values[i], thetaCer, phiCer)) {
      double sigma2 = alisSigma2.sigma2(thetaP, phiP, thetaCer, phiCer);
      phi_cer_values[i] = phiCer;
      theta_cer_values[i] = thetaCer;
      sigma_ring_values[i] = sigma2;
      ////std::cout << "Values: phiCer " << phiCer << " thetaCer " << thetaCer
      ///<< " sigma2 " << sigma2 << std::endl;
    } /*else {
      //std::cout << "Values not found." << std::endl;
    }*/
  }

  // datasetSigmaRing datasetThetaCer datasetPhiCer

  ////std::cout << "		createOrUpdateDataset(particleGroup, sigmaRingValues,
  ///sigma_ring_values);" << std::endl;
  createOrUpdateDataset(particleGroup, "sigmaRingValues", sigma_ring_values);

  // std::cout << "		createOrUpdateDataset(particleGroup, thetaCerValues,
  // thetaCerValues);" << std::endl;
  createOrUpdateDataset(particleGroup, "thetaCerValues", theta_cer_values);

  // std::cout << "		createOrUpdateDataset(particleGroup, phiCerValues,
  // phi_cer_values);" << std::endl;
  createOrUpdateDataset(particleGroup, "phiCerValues", phi_cer_values);

  // updateDataset(particleGroup, "phiCerValues", phi_cer_values);

  ////std::cout << "exit updateParams" << std::endl;

  // updateParams();
}

// Function to read data from H5 file

// Function to write data to H5 file

void writeH5File(const std::string &filename,
                 const std::vector<ParticleInfo2> &particleInfos) {

  H5File file(filename, H5F_ACC_RDWR); // Open file in read-write mode

  // H5File file(filename, H5F_ACC_TRUNC);

  // Write data from ParticleInfo2 objects to H5 file

  // ...
}

std::vector<ParticleInfo2> readH5File(const std::string &fileName) {

  // H5File file(filename, H5F_ACC_RDONLY);

  try {

    H5::H5File file(fileName, H5F_ACC_RDWR); // Open file in read-write mode

    const int num_particles =
        file.getNumObjs(); // Assuming this gets the number of groups

    for (int i = 0; i < num_particles; ++i) {

      std::string groupName = "Particle" + std::to_string(i);

      Group particleGroup = file.openGroup(groupName);

      ////std::cout << "Reading group " << groupName << std::endl;

      // ParticleInfo2 pInfo;

      ////std::cout << "Reading "<< i << " of  " << num_particles << std::endl;

      // ParticleInfo2(float xRad, float yRad, float xMip, float yMip, float th,
      // float ph, float refIndex, int cluCharge, int cluSize, float p, int
      // mcTrackPdg, float ckovRecon)

      double momentum, thetaP, phiP, xMip, yMip, xRad, yRad = 0;

      particleGroup.openAttribute("Momentum")
          .read(PredType::NATIVE_DOUBLE, &momentum);

      particleGroup.openAttribute("ThetaP").read(PredType::NATIVE_DOUBLE,
                                                 &thetaP);

      particleGroup.openAttribute("PhiP").read(PredType::NATIVE_DOUBLE, &phiP);

      particleGroup.openAttribute("xMip").read(PredType::NATIVE_DOUBLE, &xMip);

      particleGroup.openAttribute("yMip").read(PredType::NATIVE_DOUBLE, &yMip);

      particleGroup.openAttribute("xRad").read(PredType::NATIVE_DOUBLE, &xRad);

      particleGroup.openAttribute("yRad").read(PredType::NATIVE_DOUBLE, &yRad);

      ////std::cout<< "xRad " <<xRad << std::endl;

      updateParams(particleGroup, momentum, thetaP, phiP, xMip, yMip, xRad,
                   yRad, fileName);

      // class SimpleCluster

      // SimpleCluster(double x, double y, double q, double size)

      ////std::cout<< "updateParams finished "  << std::endl;

      // Recon reconObj(thetaP,  phiP,  xMip,  yMip,  xRad,  yRad,  clusters);

      // Recon reconObj(double thetaP, double phiP, double xMip, double yMip,
      // double xRad, double yRad, std::vector<SimpleCluster> clustersIn) {

      // Alisigma2_ alisSigma2(1.2904);

      // Process each cluster in parallel

    } // end for num_particles

  } catch (const H5::FileIException &e) {

    // Handle file access errors

    std::cerr << "File access error: " << e.getCDetailMsg() << std::endl;

  } catch (const H5::GroupIException &e) {

    // Handle group access errors

    std::cerr << "Group access error: " << e.getCDetailMsg() << std::endl;

  } catch (const H5::Exception &e) {

    // Handle other HDF5 related errors

    std::cerr << "HDF5 error: " << e.getCDetailMsg() << std::endl;

  } catch (const std::exception &e) {

    // Handle other standard exceptions

    std::cerr << "Standard exception: " << e.what() << std::endl;

  } catch (...) {

    // Handle all other types of exceptions

    std::cerr << "Unknown exception caught" << std::endl;
  }

  // void updateDataset(Group& group, const std::string& name, const
  // std::vector<double>& newValues)

  return std::vector<ParticleInfo2>(); // Placeholder return
}

int UpdateH5(std::string fileName = "Particle2212_Ckov10_.h5") {

  auto particleInfos = readH5File(fileName);

  // readAndPrintDataset(fileName, "thetaCerValues"/*, theta_cer_values*/);
  // readAndPrintDataset(fileName, "sigmaRingValues"/*, sigma_ring_values*/);
  // readAndPrintDataset(fileName, "phiCerValues"/*, phiCerValues*/);

  // writeH5File(filename, particleInfos);

  return 0;
}
