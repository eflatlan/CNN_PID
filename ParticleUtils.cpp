#include "H5Cpp.h"
#include <iostream>
#include <vector>
#include "TH2F.h"
using vecArray2 = std::vector<std::array<double,2>>;
using namespace H5;

class ParticleUtils {
public:
    struct Candidate2 {
        double x, y = 0.;
        int candStatus = 0;
    };

    struct ParticleInfo {

  			vecArray2 pionCandidates, kaonCandidates, protonCandidates;
        std::vector<Candidate2> candsCombined;
        std::array<int, 4> arrayInfo;
        float momentum;
        float mass;
        float energy;
        float refractiveIndex;
        float ckov;
        float xRad;
        float yRad;
        float xPC;
        float yPC;
        float thetaP;
        float phiP;
    };

    static void saveParticleInfoToHDF5(std::vector<ParticleInfo>& particleVector) {
        H5File file("ParticleInfo.h5", H5F_ACC_TRUNC);

        for (size_t i = 0; i < particleVector.size(); ++i) {
            auto& particle = particleVector[i];
            std::string groupName = "Particle" + std::to_string(i);

	    // std::cout << groupName << std::endl;

            Group particleGroup = file.createGroup(groupName);

            DataSpace attr_dataspace = DataSpace(H5S_SCALAR);

            // Store scalar values
            Attribute attribute = particleGroup.createAttribute("Momentum", PredType::NATIVE_FLOAT, attr_dataspace);
            attribute.write(PredType::NATIVE_FLOAT, &particle.momentum);

            attribute = particleGroup.createAttribute("Mass", PredType::NATIVE_FLOAT, attr_dataspace);
            attribute.write(PredType::NATIVE_FLOAT, &particle.mass);

            attribute = particleGroup.createAttribute("Energy", PredType::NATIVE_FLOAT, attr_dataspace);
            attribute.write(PredType::NATIVE_FLOAT, &particle.energy);

            attribute = particleGroup.createAttribute("RefractiveIndex", PredType::NATIVE_FLOAT, attr_dataspace);
            attribute.write(PredType::NATIVE_FLOAT, &particle.refractiveIndex);

            attribute = particleGroup.createAttribute("Ckov", PredType::NATIVE_FLOAT, attr_dataspace);
            attribute.write(PredType::NATIVE_FLOAT, &particle.ckov);

            attribute = particleGroup.createAttribute("xRad", PredType::NATIVE_FLOAT, attr_dataspace);
            attribute.write(PredType::NATIVE_FLOAT, &particle.xRad);

            attribute = particleGroup.createAttribute("yRad", PredType::NATIVE_FLOAT, attr_dataspace);
            attribute.write(PredType::NATIVE_FLOAT, &particle.yRad);

            attribute = particleGroup.createAttribute("xPC", PredType::NATIVE_FLOAT, attr_dataspace);
            attribute.write(PredType::NATIVE_FLOAT, &particle.xPC);

            attribute = particleGroup.createAttribute("yPC", PredType::NATIVE_FLOAT, attr_dataspace);
            attribute.write(PredType::NATIVE_FLOAT, &particle.yPC);

            attribute = particleGroup.createAttribute("ThetaP", PredType::NATIVE_FLOAT, attr_dataspace);
            attribute.write(PredType::NATIVE_FLOAT, &particle.thetaP);

            attribute = particleGroup.createAttribute("PhiP", PredType::NATIVE_FLOAT, attr_dataspace);
            attribute.write(PredType::NATIVE_FLOAT, &particle.phiP);

            // Write the arrayInfo
            hsize_t array_dims[1] = {4};
            DataSpace array_space(1, array_dims);
            DataSet array_dataset = particleGroup.createDataSet("ArrayInfo", PredType::NATIVE_INT, array_space);
            array_dataset.write(particle.arrayInfo.data(), PredType::NATIVE_INT);

            // Write the candsCombined
            H5::CompType mtype(sizeof(Candidate2));
            mtype.insertMember("x", HOFFSET(Candidate2, x), H5::PredType::NATIVE_DOUBLE);
            mtype.insertMember("y", HOFFSET(Candidate2, y), H5::PredType::NATIVE_DOUBLE);
            mtype.insertMember("candStatus", HOFFSET(Candidate2, candStatus),  H5::PredType::NATIVE_INT);

            hsize_t dims[1] = { particle.candsCombined.size() };
            H5::DataSpace dataspace(1, dims);

            H5::DataSet dataset = particleGroup.createDataSet("candsCombined", mtype, dataspace);
            dataset.write(particle.candsCombined.data(), mtype);
        }
    }

    static std::vector<ParticleInfo> loadParticleInfoFromHDF5(const std::string& filename) {
        std::vector<ParticleInfo> particleVector;

        H5File file(filename, H5F_ACC_RDONLY);

        for (hsize_t i = 0; i < file.getNumObjs(); ++i) {
            std::string groupName = file.getObjnameByIdx(i);
            Group particleGroup = file.openGroup(groupName);

            // Read simple attributes
            // (The rest of the code is the same as before without the removed parts)

            // Read candsCombined
            H5::CompType mtype(sizeof(Candidate2));
            mtype.insertMember("x", HOFFSET(Candidate2, x), H5::PredType::NATIVE_DOUBLE);
            mtype.insertMember("y", HOFFSET(Candidate2, y), H5::PredType::NATIVE_DOUBLE);
            mtype.insertMember("candStatus", HOFFSET(Candidate2, candStatus), H5::PredType::NATIVE_INT);

            ParticleInfo particleInfo;

            hsize_t dims[1] = { particleInfo.candsCombined.size() };
            H5::DataSpace dataspace(1, dims);
            DataSet dataset = particleGroup.openDataSet("candsCombined");
            dataspace = dataset.getSpace();

            hsize_t dims_out[1];
            dataspace.getSimpleExtentDims(dims_out, NULL);

            std::vector<Candidate2> candsCombined(dims_out[0]);
            dataset.read(candsCombined.data(), mtype);

            // (The rest of the code is the same as before without the removed parts)
            
            particleVector.push_back(particleInfo);
        }

        return particleVector;
    }
};

