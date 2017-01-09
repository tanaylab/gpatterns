using namespace std;
/*
 * MixtureModel.cpp
 *
 *  Created on: Sep 16, 2014
 *      Author: aviezerl
 */

//BASE_CC_FILE

#include "MixtureModel.h"
#include "ModelData.h"
#include "util.h"
#include <chrono>
#include <random>
#include <functional>

MixtureModel::MixtureModel():
m_loglikelihood(-std::numeric_limits<double>::max()), m_delta_likelihood_converged(0.001), m_number_of_models(0) {}

MixtureModel::MixtureModel(const MixtureModel& m): m_loglikelihood(m.m_loglikelihood), m_delta_likelihood_converged(m.m_delta_likelihood_converged), m_number_of_models(0) { }

MixtureModel::~MixtureModel() {
    m_mixing_probs.clear();
    m_models.clear();
    m_mixture_posteriors.clear();
}

float MixtureModel::get_log_prob(const Sample* s) {
    float l_prob = -std::numeric_limits<double>::max();
    for (unsigned int c=0; c<m_models.size(); c++) {
        log_sum_log(l_prob, log(m_mixing_probs[c]) + m_models[c]->get_log_prob(s));
        //prob += m_mixing_probs[c] * exp(m_models[c]->get_log_prob(s));
    }
    return(l_prob);
}

float MixtureModel::loglikelihood(const ModelData& data) {
    m_loglikelihood=0;
    for (int i=0; i<data.get_number_of_samples(); i++) {
        m_loglikelihood += get_log_prob(data.get_sample(i));
    }
    return(m_loglikelihood);
}

int MixtureModel::learn(const ModelData& data) {
    //initialize models
    init_models(data);
    float old_likelihood = loglikelihood(data);
    //	print(cerr);

    init_posteriors(data);
    //    m_mixture_posteriors.resize(m_models.size(), vector<float>(data.get_number_of_samples()));
    //run EM
    int iter=0;
    float delta = 1;
    while(delta > m_delta_likelihood_converged) {
//        cerr << endl<<  "ITER :: " << iter << endl;
        e_step(data);
        m_step(data);
        loglikelihood(data);
        delta = m_loglikelihood - old_likelihood;
        old_likelihood = m_loglikelihood;
//        print(cerr);
        iter++;
    }
    return(0);
}

void MixtureModel::e_step(const ModelData& data) {
    //	cerr <<"E-STEP" << endl;
    for (int i=0; i<data.get_number_of_samples(); i++) {
        float total_sample_posteriors=0;
        for (unsigned int c=0; c<m_models.size(); c++) {
            m_mixture_posteriors[c][i] = m_mixing_probs[c] * exp(m_models[c]->get_log_prob(data.get_sample(i)));
            total_sample_posteriors += m_mixture_posteriors[c][i];
        }
        //		cerr << "\tmix_posteriors [" << i << "] :: ";
        for (unsigned int c=0; c<m_models.size(); c++) {
            m_mixture_posteriors[c][i] /= total_sample_posteriors;
            //			cerr << "\t" << m_mixture_posteriors[c][i];
        }
        //		cerr << endl;
    }
}

float MixtureModel::m_step(const ModelData& data) {
    //	cerr <<"M-STEP" << endl;
    for (unsigned int c=0; c<m_models.size(); c++) {
        m_models[c]->learn(data, m_mixture_posteriors[c]);
        //		m_models[c]->print(cerr);
        m_mixing_probs[c] = 0;
        for (int i=0; i<data.get_number_of_samples(); i++) {
            m_mixing_probs[c] += m_mixture_posteriors[c][i];
        }
        m_mixing_probs[c] /= data.get_number_of_samples();
    }
    return(1);
}

void MixtureModel::print(ostream& out) {
    out << endl;
    for (unsigned int i=0; i<m_models.size(); i++) {
        out << i << " :: mixing prob:" << m_mixing_probs[i] << "\t";
        m_models[i]->print(out);
    }
    out << "LL = " << m_loglikelihood << endl;
}

void MixtureModel::print_summary(ostream& out) {
    for (unsigned int i=0; i<m_models.size(); i++) {
        out << m_mixing_probs[i] << "\t";
        m_models[i]->print_summary(out);
        out << "\t";
    }
}

void MixtureModel::simulate(const ModelData& original, ModelData* data) {
	random_device rd;
    default_random_engine generator(rd());
    discrete_distribution<int> d(m_mixing_probs.begin(), m_mixing_probs.end());

    for (int i=0; i < (original.get_number_of_samples()); i++) {
    	m_models[d(generator)]->simulate_sample(original.get_sample(i), data);
    }

}

