/*
 * UniModel.h
 *
 *  Created on: Sep 16, 2014
 *      Author: aviezerl
 */

#ifndef UNIMODEL_H_
#define UNIMODEL_H_

#include <vector>
#include <iostream>
using namespace std;

#include "ModelData.h"


class Sample;

class UniModel {
public:
    UniModel();
    virtual ~UniModel();
    virtual int learn(const ModelData&) = 0;
    virtual void simulate(const ModelData& original, ModelData* data) = 0;
    virtual void simulate_sample(const Sample* original, ModelData* data) = 0;
    virtual int learn(const ModelData&, const vector<float>& factors) = 0;
    virtual float get_log_prob(const Sample*) const = 0;
    virtual void print(ostream&) const = 0;
    virtual void print_summary(ostream& out) const = 0;
    float get_loglikelihood() { return m_loglikelihood; }
    float loglikelihood(const ModelData& data);
protected:
    float 		m_loglikelihood;

};

#endif /* UNIMODEL_H_ */
