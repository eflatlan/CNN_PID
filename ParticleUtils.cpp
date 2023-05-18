#include "H5Cpp.h"
// sudo yum install hdf5-devel
#include "H5Cpp.h"
#include <iostream>
#include <vector>
#include "TH2F.h"

using namespace H5;

struct Bin {
    double x;
    double y;
};

class ParticleUtils {
    public:
        struct ParticleInfo {
            double momentum;
            double mass;
            double energy;
            double refractiveIndex;
            double ckov;
            // std::vector<std::vector<double>> map; // Assuming the map is a 2D array of doubles
            std::vector<Bin> filledBins;
        };

        static void saveParticleInfoToHDF5(std::vector<ParticleInfo>& particleVector) {
            H5File file("ParticleInfo.h5", H5F_ACC_TRUNC);

            // Create a compound datatype for Bin
            CompType binType(sizeof(Bin));
            binType.insertMember("X", HOFFSET(Bin, x), PredType::NATIVE_DOUBLE);
            binType.insertMember("Y", HOFFSET(Bin, y), PredType::NATIVE_DOUBLE);

            for (size_t i = 0; i < particleVector.size(); ++i) {
                const auto& particle = particleVector[i];

                // Now let's create a group for each particle
                std::string groupName = "Particle" + std::to_string(i);
                Group particleGroup = file.createGroup(groupName);

                // Store scalar values
                DataSpace attr_dataspace = DataSpace(H5S_SCALAR);
                Attribute attribute = particleGroup.createAttribute("Momentum", PredType::NATIVE_DOUBLE, attr_dataspace);
                attribute.write(PredType::NATIVE_DOUBLE, &particle.momentum);

                attribute = particleGroup.createAttribute("Mass", PredType::NATIVE_DOUBLE, attr_dataspace);
                attribute.write(PredType::NATIVE_DOUBLE, &particle.mass);

                attribute = particleGroup.createAttribute("Energy", PredType::NATIVE_DOUBLE, attr_dataspace);
                attribute.write(PredType::NATIVE_DOUBLE, &particle.energy);

                attribute = particleGroup.createAttribute("RefractiveIndex", PredType::NATIVE_DOUBLE, attr_dataspace);
                attribute.write(PredType::NATIVE_DOUBLE, &particle.refractiveIndex);

                attribute = particleGroup.createAttribute("Ckov", PredType::NATIVE_DOUBLE, attr_dataspace);
                attribute.write(PredType::NATIVE_DOUBLE, &particle.ckov);

                // Write filledBins to HDF5 file
                hsize_t binDims[1] = {particle.filledBins.size()};

                DataSpace binspace(1, binDims);
                DataSet binDataset = particleGroup.createDataSet("FilledBins", binType, binspace);
                binDataset.write(&particle.filledBins[0], binType);
            }
        }

        static std::vector<std::vector<double>> convertTH2FTo2DArray(TH2F* histogram) {
            // Create a 2D vector with the same size as the histogram
            int nBinsX = histogram->GetNbinsX();
            int nBinsY = histogram->GetNbinsY();
            std::vector<std::vector<double>> array(nBinsX, std::vector<double>(nBinsY));

            // Loop over the bins in the x and y direction
            for (int i = 1; i <= nBinsX; ++i) {
                for (int j = 1; j <= nBinsY; ++j) {
                    // Store the bin content in the array
                    array[i-1][j-1] = histogram->GetBinContent(i, j);
                }
            }

            return array;
        }

        static std::vector<ParticleInfo> loadParticleInfoFromHDF5(const std::string& filename) {
            std::vector<ParticleInfo> particleVector;

            // Create a compound datatype for Bin
            CompType binType(sizeof(Bin));
            binType.insertMember("X", HOFFSET(Bin, x), PredType::NATIVE_DOUBLE);
            binType.insertMember("Y", HOFFSET(Bin, y), PredType::NATIVE_DOUBLE);

            // Open the file
            H5File file(filename, H5F_ACC_RDONLY);

            // Iterate over all groups in the file
            for (hsize_t i = 0; i < file.getNumObjs(); ++i) {
                std::string groupName = file.getObjnameByIdx(i);

                // Open the group
                Group particleGroup = file.openGroup(groupName);

                // Read momentum
                Attribute attribute = particleGroup.openAttribute("Momentum");
                double momentum;
                attribute.read(PredType::NATIVE_DOUBLE, &momentum);

                std::cout << "Particle " << i << " Momentum: " << momentum << std::endl;

                attribute = particleGroup.openAttribute("Mass");
                double mass;
                attribute.read(PredType::NATIVE_DOUBLE, &mass);

                attribute = particleGroup.openAttribute("Energy");
                double energy;
                attribute.read(PredType::NATIVE_DOUBLE, &energy);

                attribute = particleGroup.openAttribute("RefractiveIndex");
                double refractiveIndex;
                attribute.read(PredType::NATIVE_DOUBLE, &refractiveIndex);

                attribute = particleGroup.openAttribute("Ckov");
                double ckov;
                attribute.read(PredType::NATIVE_DOUBLE, &ckov);

                // Read filledBins
                DataSet binDataset = particleGroup.openDataSet("FilledBins");
                DataSpace binDataSpace = binDataset.getSpace();

                hsize_t binDims[1];
                binDataSpace.getSimpleExtentDims(binDims, NULL);
                std::vector<Bin> filledBins(binDims[0]);
                binDataset.read(&filledBins[0], binType);

                std::cout << "Particle " << i << " FilledBins: \n";
                for (auto &bin : filledBins) {
                    std::cout << "X: " << bin.x << ", Y: " << bin.y << '\n';
                }

                // Construct a ParticleInfo and add it to the vector
                ParticleInfo particleInfo;
                particleInfo.momentum = momentum;
                particleInfo.mass = mass;
                particleInfo.energy = energy;
                particleInfo.refractiveIndex = refractiveIndex;
                particleInfo.ckov = ckov;
                particleInfo.filledBins = filledBins;
                particleVector.push_back(particleInfo);
            }

            return particleVector;
        }
    };

