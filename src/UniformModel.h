/*
 * UniformModel.h
 *
 *  Created on: Nov 17, 2014
 *      Author: aviezerl
 */

#ifndef UNIFORMMODEL_H_
#define UNIFORMMODEL_H_

#include "UniModel.h"

class UniformModel: public UniModel {
public:
	UniformModel();
	UniformModel(const UniformModel&);
	virtual ~UniformModel();
    int learn(const ModelData& data) { return(1); }
    int learn(const ModelData& data, const vector<float>& factors) { return 1; }
    void simulate(const ModelData& original, ModelData* data);
    void simulate_sample(const Sample* original, ModelData* data);
    float get_log_prob(const Sample*) const;
    void print(ostream& out) const { 	out << "Uniform-Model" << endl; return; }
    void print_summary(ostream& out) const { return; }
    void set_size(int size) { m_size = size; }

protected:
	int m_size;
};

#endif /* UNIFORMMODEL_H_ */
