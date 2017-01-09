using namespace std;
/*
 * UniformModel.cpp
 *
 *  Created on: Nov 17, 2014
 *      Author: aviezerl
 */

//BASE_CC_FILE
#include "UniformModel.h"
#include "MethylPatternData.h"

UniformModel::UniformModel(): m_size(0) {
	// TODO Auto-generated constructor stub

}

UniformModel::UniformModel(const UniformModel& m):m_size(m.m_size){
}

UniformModel::~UniformModel() {
	// TODO Auto-generated destructor stub
}

void UniformModel::simulate(const ModelData& original, ModelData* data) {
	for (int i=0; i<original.get_number_of_samples(); i++) {
		const Sample* os = original.get_sample(i);
		simulate_sample(os,data);
	}
}

void UniformModel::simulate_sample(const Sample* original, ModelData* data) {
	float r;
	const MethylPatternSample* os = (MethylPatternSample*) original;
	string s("");
	for (unsigned int i=0; i<os->pattern.size(); i++) {
		if (os->pattern[i] < 0)
			s += "*";
		else {
			r = ((float)rand())/RAND_MAX;
			if (r > 0.5) {
				s += "1";
			} else {
				s += "0";
			}
		}
	}
	((MethylPatternData*)data)->add_sample(s);
}

float UniformModel::get_log_prob(const Sample*) const {
	return(-m_size * log(2));
}
