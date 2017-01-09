/*
 * MixtureModel.h
 *
 *  Created on: Sep 16, 2014
 *      Author: aviezerl
 */

#ifndef MIXTUREMODEL_H_
#define MIXTUREMODEL_H_

#include <vector>
#include "UniModel.h"
using namespace std;

class MixtureModel {
public:
    MixtureModel();
    MixtureModel(const MixtureModel&);
    virtual ~MixtureModel();
    void add_model(UniModel* m) {
        m_models.push_back(m);
        m_mixing_probs.push_back(0);
        m_number_of_models++;
    }
    int learn(const ModelData&);

    UniModel* get_model(int index) const { return(m_models[index]); }
    float get_loglikelihood() const{ return m_loglikelihood; }
    int get_number_of_models() const{ return m_number_of_models; }

    float get_log_prob(const Sample* s);
    float loglikelihood(const ModelData& data);
    void print(ostream& out);
    void print_summary(ostream& out);

    void simulate(const ModelData& original, ModelData* data);

public:
    virtual void init_models(const ModelData& data) = 0;
    void init_posteriors(const ModelData& data){
    	m_mixture_posteriors.resize(m_models.size(), vector<float>(data.get_number_of_samples()));
    }
    virtual void e_step(const ModelData& data);
    virtual float m_step(const ModelData& data);

    void set_mix_prob(const int& c, const float& prob) { m_mixing_probs[c] = prob; }
    float get_mix_prob(const int& c) const{ return(m_mixing_probs[c]); }
    vector<float> get_mixture_posteriors(const int& c) const{ return(m_mixture_posteriors[c]); }



protected:
    float					m_loglikelihood;
    float 					m_delta_likelihood_converged;
    unsigned int						m_number_of_models;
    vector<float>			m_mixing_probs;			//pai
    vector<UniModel*> 		m_models;
    vector< vector<float> >	m_mixture_posteriors; //gama



};

#endif /* MIXTUREMODEL_H_ */
