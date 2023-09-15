#include <iostream>
#include <random>
#include <vector>
#include <algorithm>
#include <chrono>


class RandomValues {
private:

	static constexpr float nm2eV = 1239.842609;
	const float mass_Pion = 0.1396, mass_Kaon = 0.4937, mass_Proton = 0.938; // masses in
	std::array<float, 3> masses = {mass_Pion, mass_Kaon, mass_Proton};

	static constexpr float arrWaveLenDefault[30] = {
	  162, 164, 166, 168, 170, 172, 174, 176, 178, 180,
	  182, 184, 186, 188, 190, 192, 194, 196, 198, 200,
	  202, 204, 206, 208, 210, 212, 214, 216, 218, 220};

    std::random_device rd;
    std::mt19937 gen;
    std::uniform_int_distribution<int> massDistribution;
    std::uniform_int_distribution<int> energyDistribution;
    std::normal_distribution<float> momentumDistribution;
    
public:
    float momentum;
    float mass;
    float energy;
    float refractiveIndex;


    RandomValues() : gen(rd()) {
        massDistribution = std::uniform_int_distribution<int>(0, 2);
        energyDistribution = std::uniform_int_distribution<int>(0, 29);
        momentumDistribution = std::normal_distribution<float>(0.5, 0.1);

	

        mass = getRandomMass();
        energy = getRandomEnergy();
        refractiveIndex = calculateRefractiveIndex();
        momentum = getRandomMomentum(mass, refractiveIndex);
    }
    /*
    float checkMomentumThreshold(float m, float p, float n) {
	  const float p_sq = p*p

	;
	  const float cos_ckov_denom = p*n;

	  // sanity check ;)
	  if(p_sq + m*m < 0){
	    return 0;
	  }
          const auto cos_ckov = static_cast<float>(TMath::Sqrt(p_sq + m*m)/(cos_ckov_denom));

	  // sanity check ;)
	  if(cos_ckov > 1 || cos_ckov < -1)
	    return 0;
	  }


    }*/ 



    /*

cosCkov in 0...1
(pn)^2 = p^2 + m^2
p*p*(n^2-1) = m^2

pThre = m/(sqrt(n^2-1)) 
cosCkov = sqrt(p^2+m^2)/(p*n)
0.1396, mass_Kaon = 0.4937, mass_Proton = 0.938
    */

    float getRandomMomentum(float m, float n) {
				const auto pThre = m/(TMath::Sqrt(n*n-1));

				// should imply range pThre..5 
      	auto p = pThre + (5-pThre) * momentumDistribution(gen);

				if(p > 5 ) {
					Printf("p %.4f | m %.4f | n %.4f", p, m, n); throw std::invalid_argument("vaffanculo");
				}
				else {
					return p;
				}
/*			 

								if()

				else if(m == mass_Kaon) 
        	 return pThre + 4 * momentumDistributionKaon(gen);

				else if(m == mass_Pion) 
        	 return 1 + 4 * momentumDistributionKaon(gen);*/
    }

    float getRandomMass() {
        int index = massDistribution(gen);
        return masses[index];
    }

    float getRandomEnergy() {
        int index = energyDistribution(gen);
        float photonEnergy = nm2eV / arrWaveLenDefault[index];
        return photonEnergy;
    }

   float calculateRefractiveIndex() {
        float photonEnergy = energy;
        float k = 1.177 + (0.0172) * photonEnergy;
        return k;
    }
};

/*
example usage :
    int numObjects = 5;  // Number of objects to create
    std::vector<RandomValues> randomObjects(numObjects);

*/
