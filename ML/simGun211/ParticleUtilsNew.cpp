class ParticleUtils {
public:
    struct ParticleInfo {
      std::vector<double> pionCandidatesX, pionCandidatesY;
      std::vector<double> kaonCandidatesX, kaonCandidatesY;
      std::vector<double> protonCandidatesX, protonCandidatesY;
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
            const auto& particle = particleVector[i];

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
            
            Attribute attribute = particleGroup.createAttribute("ThetaP", PredType::NATIVE_FLOAT, attr_dataspace);
            attribute.write(PredType::NATIVE_FLOAT, &particle.thetaP);

            attribute = particleGroup.createAttribute("PhiP", PredType::NATIVE_FLOAT, attr_dataspace);
            attribute.write(PredType::NATIVE_FLOAT, &particle.phiP);

            // Write the arrayInfo
            hsize_t array_dims[1] = {4};
            DataSpace array_space(1, array_dims);
            DataSet array_dataset = particleGroup.createDataSet("ArrayInfo", PredType::NATIVE_INT, array_space);
            array_dataset.write(particle.arrayInfo.data(), PredType::NATIVE_INT);

            // Write the vectors data
            hsize_t vector_dims[1] = {particle.pionCandidatesX.size()};
            DataSpace vector_space(1, vector_dims);

            // PionCandidates
            DataSet vector_dataset = particleGroup.createDataSet("PionCandidatesX", PredType::NATIVE_DOUBLE, vector_space);
            vector_dataset.write(particle.pionCandidatesX.data(), PredType::NATIVE_DOUBLE);

            vector_dataset = particleGroup.createDataSet("PionCandidatesY", PredType::NATIVE_DOUBLE, vector_space);
            vector_dataset.write(particle.pionCandidatesY.data(), PredType::NATIVE_DOUBLE);

            // KaonCandidates
            vector_dims[0] = particle.kaonCandidatesX.size();
            vector_space = DataSpace(1, vector_dims);

            vector_dataset = particleGroup.createDataSet("KaonCandidatesX", PredType::NATIVE_DOUBLE, vector_space);
            vector_dataset.write(particle.kaonCandidatesX.data(), PredType::NATIVE_DOUBLE);

            vector_dataset = particleGroup.createDataSet("KaonCandidatesY", PredType::NATIVE_DOUBLE, vector_space);
            vector_dataset.write(particle.kaonCandidatesY.data(), PredType::NATIVE_DOUBLE);

            // ProtonCandidates
            vector_dims[0] = particle.protonCandidatesX.size();
            vector_space = DataSpace(1, vector_dims);

            vector_dataset = particleGroup.createDataSet("ProtonCandidatesX", PredType::NATIVE_DOUBLE, vector_space);
            vector_dataset.write(particle.protonCandidatesX.data(), PredType::NATIVE_DOUBLE);

            vector_dataset = particleGroup.createDataSet("ProtonCandidatesY", PredType::NATIVE_DOUBLE, vector_space);
            vector_dataset.write(particle.protonCandidatesY.data(), PredType::NATIVE_DOUBLE);
        }
    }

    static std::vector<ParticleInfo> loadParticleInfoFromHDF5(const std::string& filename) {
        std::vector<ParticleInfo> particleVector;

        H5File file(filename, H5F_ACC_RDONLY);

        for (hsize_t i = 0; i < file.getNumObjs(); ++i) {
            std::string groupName = file.getObjnameByIdx(i);
            Group particleGroup = file.openGroup(groupName);

            // Read simple attributes
            // Continue as before until "Ckov"

            Attribute attribute = particleGroup.openAttribute("PhiP");
            float phiP;
            attribute.read(PredType::NATIVE_FLOAT, &phiP);

            // Read arrayInfo
            DataSet array_dataset = particleGroup.openDataSet("ArrayInfo");
            std::array<int, 4> arrayInfo;
            array_dataset.read(arrayInfo.data(), PredType::NATIVE_INT);

            // Read vector data
            hsize_t dims_out[1];
            DataSet vector_dataset = particleGroup.openDataSet("PionCandidatesX");
            DataSpace vector_space = vector_dataset.getSpace();

            vector_space.getSimpleExtentDims(dims_out, NULL);
            std::vector<double> pionCandidatesX(dims_out[0]);
            vector_dataset.read(pionCandidatesX.data(), PredType::NATIVE_DOUBLE);

            vector_dataset = particleGroup.openDataSet("PionCandidatesY");
            vector_space = vector_dataset.getSpace();

            vector_space.getSimpleExtentDims(dims_out, NULL);
            std::vector<double> pionCandidatesY(dims_out[0]);
            vector_dataset.read(pionCandidatesY.data(), PredType::NATIVE_DOUBLE);

            // Repeat the same for KaonCandidatesX, KaonCandidatesY, ProtonCandidatesX and ProtonCandidatesY...

            // Construct a ParticleInfo and add it to the vector
            ParticleInfo particleInfo;
            particleInfo.momentum = momentum;
            particleInfo.mass = mass;
            particleInfo.energy = energy;
            particleInfo.refractiveIndex = refractiveIndex;
            particleInfo.ckov = ckov;
            particleInfo.xRad = xRad;
            particleInfo.yRad = yRad;
            particleInfo.xPC = xPC;
            particleInfo.yPC = yPC;
            particleInfo.thetaP = thetaP;
            particleInfo.phiP = phiP;
            particleInfo.arrayInfo = arrayInfo;
            particleInfo.pionCandidatesX = pionCandidatesX;
            particleInfo.pionCandidatesY = pionCandidatesY;
            // Same for kaonCandidates and protonCandidates...

            particleVector.push_back(particleInfo);
        }

        return particleVector;
    }
};

