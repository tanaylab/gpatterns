using namespace std;
/*
 * MethylPatternMixtureModelsCompare.cpp
 *
 *  Created on: Aug 6, 2015
 *      Author: aviezerl
 */

#include "MethylPatternMixtureModelsCompare.h"


MethylPatternMixtureModelsCompare::MethylPatternMixtureModelsCompare() {
    // TODO Auto-generated constructor stub

}

MethylPatternMixtureModelsCompare::~MethylPatternMixtureModelsCompare() {
    // TODO Auto-generated destructor stub
}

void MethylPatternMixtureModelsCompare::compare_models(const MethylPatternData& data, MethylPatternMixtureModelScanner* model1, MethylPatternMixtureModelScanner* model2, const int& K, const int& min_k, const int& num_of_sims, ostream& model1_out, ostream& model2_out){

	model1->scan(data, K, min_k, model1_out);

	MethylPatternMixtureModel* mix = model1->get_mix();
	MethylPatternData sim_data;
	for (int i=0; i < num_of_sims; ++i){
		mix->simulate(data, &sim_data);
		model2->scan(sim_data, K, min_k, model2_out);
		sim_data.clear();
	}

}


