#include "H5Cpp.h"
#include <iostream>
#include <vector>
#include "TH2F.h"

using namespace H5;

using vecArray2 = std::vector<std::array<double,2>>;

struct Bin {
    float x;
    float y;
};

class ParticleUtils {
public:
    struct ParticleInfo {
	
      vecArray2 pionCandidates, kaonCandidates, protonCandidates;
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
	    // TH2F* map;
      // std::vector<Bin> filledBins;
    };

    static void saveParticleInfoToHDF5(std::vector<ParticleInfo>& particleVector) {
        H5File file("ParticleInfo.h5", H5F_ACC_TRUNC);

        // Create a compound datatype for Bin
        CompType binType(sizeof(Bin));
        binType.insertMember("X", HOFFSET(Bin, x), PredType::NATIVE_FLOAT);
        binType.insertMember("Y", HOFFSET(Bin, y), PredType::NATIVE_FLOAT);

        for (size_t i = 0; i < particleVector.size(); ++i) {
            const auto& particle = particleVector[i];

            // Now let's create a group for each particle
            std::string groupName = "Particle" + std::to_string(i);
            Group particleGroup = file.createGroup(groupName);

            // Store scalar values
            DataSpace attr_dataspace = DataSpace(H5S_SCALAR);
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

            // Write filledBins to HDF5 file
            /*hsize_t binDims[1] = {particle.filledBins.size()};

            DataSpace binspace(1, binDims);
            DataSet binDataset = particleGroup.createDataSet("FilledBins", binType, binspace);
            binDataset.write(&particle.filledBins[0], binType);*/ 
        }
    }

    static std::vector<ParticleInfo> loadParticleInfoFromHDF5(const std::string& filename) {
        std::vector<ParticleInfo> particleVector;

        // Create a compound datatype for Bin
        CompType binType(sizeof(Bin));
        binType.insertMember("X", HOFFSET(Bin, x), PredType::NATIVE_FLOAT);
        binType.insertMember("Y", HOFFSET(Bin, y), PredType::NATIVE_FLOAT);

        // Open the file
        H5File file(filename, H5F_ACC_RDONLY);

        // Iterate over all groups in the file
        for (hsize_t i = 0; i < file.getNumObjs(); ++i) {
            std::string groupName = file.getObjnameByIdx(i);

            // Open the group
            Group particleGroup = file.openGroup(groupName);

            // Read momentum
            Attribute attribute = particleGroup.openAttribute("Momentum");
            float momentum;
            attribute.read(PredType::NATIVE_FLOAT, &momentum);

            std::cout << "Particle " << i << " Momentum: " << momentum << std::endl;

            attribute = particleGroup.openAttribute("Mass");
            float mass;
            attribute.read(PredType::NATIVE_FLOAT, &mass);

            attribute = particleGroup.openAttribute("Energy");
            float energy;
            attribute.read(PredType::NATIVE_FLOAT, &energy);

            attribute = particleGroup.openAttribute("RefractiveIndex");
            float refractiveIndex;
            attribute.read(PredType::NATIVE_FLOAT, &refractiveIndex);

            attribute = particleGroup.openAttribute("Ckov");
            float ckov;
            attribute.read(PredType::NATIVE_FLOAT, &ckov);

            /*/ Read filledBins
            DataSet binDataset = particleGroup.openDataSet("FilledBins");
            DataSpace binDataSpace = binDataset.getSpace();

            hsize_t binDims[1];
            binDataSpace.getSimpleExtentDims(binDims, NULL);
            std::vector<Bin> filledBins(binDims[0]);
            binDataset.read(&filledBins[0], binType);

            std::cout << "Particle " << i << " FilledBins: \n";
            for (auto& bin : filledBins) {
                std::cout << "X: " << bin.x << ", Y: " << bin.y << '\n';
            } */

            // Construct a ParticleInfo and add it to the vector
            ParticleInfo particleInfo;
            particleInfo.momentum = momentum;
            particleInfo.mass = mass;
            particleInfo.energy = energy;
            particleInfo.refractiveIndex = refractiveIndex;
            particleInfo.ckov = ckov;
            // particleInfo.filledBins = filledBins;
            particleVector.push_back(particleInfo);
        }

        return particleVector;
    }
};
