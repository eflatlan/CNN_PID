#include "H5Cpp.h"
// sudo yum install hdf5-devel

#include <iostream>
#include <vector>
#include "TH2F.h"

using namespace H5;

namespace ParticleUtils {

    class ParticleDataUtils {
    public:
        struct ParticleInfo {
            double momentum;
            double mass;
            double energy;
            double refractiveIndex;
            double ckov;
            std::vector<std::vector<double>> map; // Assuming the map is a 2D array of doubles
        };

        static void saveParticleInfoToHDF5(std::vector<ParticleInfo>& particleVector) {
            H5File file("ParticleInfo.h5", H5F_ACC_TRUNC);

            for (size_t i = 0; i < particleVector.size(); ++i) {
                const auto& particle = particleVector[i];

                // Now let's create a group for each particle
                std::string groupName = "Particle" + std::to_string(i);
                Group particleGroup = file.createGroup(groupName);

                // Store scalar values
                // For the sake of brevity, here I am showing it for momentum only. You can do it similarly for other scalar attributes.
                DataSpace attr_dataspace = DataSpace(H5S_SCALAR);
                Attribute attribute = particleGroup.createAttribute("Momentum", PredType::NATIVE_DOUBLE, attr_dataspace);
                attribute.write(PredType::NATIVE_DOUBLE, &particle.momentum);

                // Store 2D array "map"
                hsize_t dims[2] = {particle.map.size(), particle.map[0].size()};
                DataSpace dataspace(2, dims);

                DataSet dataset = particleGroup.createDataSet("Map", PredType::NATIVE_DOUBLE, dataspace);
                dataset.write(&particle.map[0][0], PredType::NATIVE_DOUBLE);
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
    };
}

