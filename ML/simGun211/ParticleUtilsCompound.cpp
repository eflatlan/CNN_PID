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
	//double R = 0.;
	//double phiL = 0., phi = 0.;
	/*uint32_t */int candStatus = 0;//(000, 001, 010, 011, 100, 110, 101, 111);
};



    struct ParticleInfo {

        std::vector<Candidate2> candsCombined;


        std::vector<double> pionCandidatesX, pionCandidatesY;
        std::vector<double> kaonCandidatesX, kaonCandidatesY;
        std::vector<double> protonCandidatesX, protonCandidatesY;

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
    };

    static void saveParticleInfoToHDF5(std::vector<ParticleInfo>& particleVector) {
        H5File file("ParticleInfo.h5", H5F_ACC_TRUNC);

        for (size_t i = 0; i < particleVector.size(); ++i) {
            auto& particle = particleVector[i];

						for (const auto& candidate : particle.pionCandidates) {
    					particle.pionCandidatesX.push_back(candidate[0]);
    					particle.pionCandidatesY.push_back(candidate[1]);
						} 				


						for (const auto& candidate : particle.kaonCandidates) {
    					particle.kaonCandidatesX.push_back(candidate[0]);
    					particle.kaonCandidatesY.push_back(candidate[1]);
						} 		


						for (const auto& candidate : particle.protonCandidates) {
    					particle.protonCandidatesX.push_back(candidate[0]);
    					particle.protonCandidatesY.push_back(candidate[1]);
						} 		
																	
						// Print all scalar values

						/*
						Printf("\n\n=======================================");
						Printf("=======================================");
						Printf("	Particle %d", i);
						Printf("=======================================");
						printf("Momentum: %.2f\n", particle.momentum);
						printf("Mass: %.2f\n", particle.mass);
						printf("Energy: %.2f\n", particle.energy);
						printf("Refractive Index: %.2f\n", particle.refractiveIndex);
						printf("Ckov: %.2f\n", particle.ckov);
						printf("xRad: %.2f\n", particle.xRad);
						printf("yRad: %.2f\n", particle.yRad);
						printf("xPC: %.2f\n", particle.xPC);
						printf("yPC: %.2f\n", particle.yPC);
						printf("ThetaP: %.2f\n", particle.thetaP);
						printf("PhiP: %.2f\n", particle.phiP);
			
						// Print vector sizes
						printf("Size of pionCandidatesX: %lu\n", particle.pionCandidatesX.size());
						printf("Size of pionCandidatesY: %lu\n", particle.pionCandidatesY.size());
						printf("Size of kaonCandidatesX: %lu\n", particle.kaonCandidatesX.size());
						printf("Size of kaonCandidatesY: %lu\n", particle.kaonCandidatesY.size());
						printf("Size of protonCandidatesX: %lu\n", particle.protonCandidatesX.size());
						printf("Size of protonCandidatesY: %lu\n", particle.protonCandidatesY.size());



						printf("Size of pionCandidates: %lu\n", particle.pionCandidates.size());
						printf("Size of kaonCandidates: %lu\n", particle.kaonCandidates.size()); */ 





            std::string groupName = "Particle" + std::to_string(i);
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



						/* The cands-combined */


						H5::CompType mtype(sizeof(Candidate2));
 				   	mtype.insertMember("x", HOFFSET(Candidate2, x), H5::PredType::NATIVE_DOUBLE);
    				mtype.insertMember("y", HOFFSET(Candidate2, y), H5::PredType::NATIVE_DOUBLE);
    				mtype.insertMember("candStatus", HOFFSET(Candidate2, candStatus), 	H5::PredType::NATIVE_INT);

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
            Attribute attribute = particleGroup.openAttribute("Momentum");
            float momentum;
            attribute.read(PredType::NATIVE_FLOAT, &momentum);

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

            attribute = particleGroup.openAttribute("xRad");
            float xRad;
            attribute.read(PredType::NATIVE_FLOAT, &xRad);

            attribute = particleGroup.openAttribute("yRad");
            float yRad;
            attribute.read(PredType::NATIVE_FLOAT, &yRad);

            attribute = particleGroup.openAttribute("xPC");
            float xPC;
            attribute.read(PredType::NATIVE_FLOAT, &xPC);

            attribute = particleGroup.openAttribute("yPC");
            float yPC;
            attribute.read(PredType::NATIVE_FLOAT, &yPC);

            attribute = particleGroup.openAttribute("ThetaP");
            float thetaP;
            attribute.read(PredType::NATIVE_FLOAT, &thetaP);

            attribute = particleGroup.openAttribute("PhiP");
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

            // Read KaonCandidates
            vector_dataset = particleGroup.openDataSet("KaonCandidatesX");
            vector_space = vector_dataset.getSpace();

            vector_space.getSimpleExtentDims(dims_out, NULL);
            std::vector<double> kaonCandidatesX(dims_out[0]);
            vector_dataset.read(kaonCandidatesX.data(), PredType::NATIVE_DOUBLE);

            vector_dataset = particleGroup.openDataSet("KaonCandidatesY");
            vector_space = vector_dataset.getSpace();

            vector_space.getSimpleExtentDims(dims_out, NULL);
            std::vector<double> kaonCandidatesY(dims_out[0]);
            vector_dataset.read(kaonCandidatesY.data(), PredType::NATIVE_DOUBLE);

            // Read ProtonCandidates
            vector_dataset = particleGroup.openDataSet("ProtonCandidatesX");
            vector_space = vector_dataset.getSpace();

            vector_space.getSimpleExtentDims(dims_out, NULL);
            std::vector<double> protonCandidatesX(dims_out[0]);
            vector_dataset.read(protonCandidatesX.data(), PredType::NATIVE_DOUBLE);

            vector_dataset = particleGroup.openDataSet("ProtonCandidatesY");
            vector_space = vector_dataset.getSpace();

            vector_space.getSimpleExtentDims(dims_out, NULL);
            std::vector<double> protonCandidatesY(dims_out[0]);
            vector_dataset.read(protonCandidatesY.data(), PredType::NATIVE_DOUBLE);



						/* 
						Cands combined : 
						*/ 

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
				    

				    dataspace.getSimpleExtentDims(dims_out, NULL);
				    
				    std::vector<Candidate2> candsCombined(dims_out[0]);
				    dataset.read(candsCombined.data(), mtype);

				    for (const auto& candidate : candsCombined) {
				        std::cout << "x: " << candidate.x << ", y: " << candidate.y << ", candStatus: " << candidate.candStatus << std::endl;
				    }


            // Construct a ParticleInfo and add it to the vector

            particleInfo.momentum = momentum;
            particleInfo.candsCombined = candsCombined;
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
            particleInfo.kaonCandidatesX = kaonCandidatesX;
            particleInfo.kaonCandidatesY = kaonCandidatesY;
            particleInfo.protonCandidatesX = protonCandidatesX;
            particleInfo.protonCandidatesY = protonCandidatesY;
											
						// Print all scalar values
						Printf("\n\n=======================================");
						Printf("=======================================");
						Printf("	Particle %d", i);
						Printf("=======================================");
						printf("Momentum: %.2f\n", particleInfo.momentum);
						printf("Mass: %.2f\n", particleInfo.mass);
						printf("Energy: %.2f\n", particleInfo.energy);
						printf("Refractive Index: %.2f\n", particleInfo.refractiveIndex);
						printf("Ckov: %.2f\n", particleInfo.ckov);
						printf("xRad: %.2f\n", particleInfo.xRad);
						printf("yRad: %.2f\n", particleInfo.yRad);
						printf("xPC: %.2f\n", particleInfo.xPC);
						printf("yPC: %.2f\n", particleInfo.yPC);
						printf("ThetaP: %.2f\n", particleInfo.thetaP);
						printf("PhiP: %.2f\n", particleInfo.phiP);
			
						// Print vector sizes
						printf("Size of pionCandidatesX: %lu\n", particleInfo.pionCandidatesX.size());
						printf("Size of pionCandidatesY: %lu\n", particleInfo.pionCandidatesY.size());
						printf("Size of kaonCandidatesX: %lu\n", particleInfo.kaonCandidatesX.size());
						printf("Size of kaonCandidatesY: %lu\n", particleInfo.kaonCandidatesY.size());
						printf("Size of protonCandidatesX: %lu\n", particleInfo.protonCandidatesX.size());
						printf("Size of protonCandidatesY: %lu\n", particleInfo.protonCandidatesY.size());

						printf("particleInfo.arrayInfo %lu %lu %lu %lu \n", particleInfo.arrayInfo[0], particleInfo.arrayInfo[1], particleInfo.arrayInfo[2], particleInfo.arrayInfo[3]);

            particleVector.push_back(particleInfo);
        }

        return particleVector;
    }
};

