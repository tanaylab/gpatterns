using namespace std;
/*
 * MethylPatternMixtureModelScanner.cpp
 *
 *  Created on: Aug 2, 2015
 *      Author: aviezerl
 */

#include "MethylPatternMixtureModelScanner.h"

#define GOOD_PVAL 0.0

MethylPatternMixtureModelScanner::MethylPatternMixtureModelScanner(const float& unimodal_uniform_mix, const int& num_of_inits): m_unimodal_uniform_mix(unimodal_uniform_mix), m_num_of_inits(num_of_inits) {
    // TODO Auto-generated constructor stub

}

MethylPatternMixtureModelScanner::~MethylPatternMixtureModelScanner() {
    // TODO Auto-generated destructor stub
}

float MethylPatternMixtureModelScanner::scan(const MethylPatternData& data, const int& K, const int& min_k, ostream& out) {
	MethylPatternData sim_data;

	create_models();
    m_uni->learn(data);
    m_mix->learn(data);

    float log_likelihood_ratio = -2*m_uni->get_loglikelihood() + 2 * m_mix->get_loglikelihood();
    if (log_likelihood_ratio == 0) { //this means that the mixture model did not add anything !
        //there is no need to run any simulations as the best we can do is with the uni model
        out << data.get_data_id() << "\t1\tNA\t" << log_likelihood_ratio << "\t"
        << "NA" << "\t" << "NA" << "\t";
        ;
        m_mix->print_summary(out);
        m_uni->print_summary(out);
        out << endl;
        return(1);
    } else {
        int outliers=0;
        //sample likelihood ratio from unimodel sampled data and compute likelihood ratio p-value
        int samples=0;
        bool continue_sampling=true;
        float sum_ratio=0;
        vector<float> ratios;

        while (continue_sampling && samples < K) {
            samples++;
            sim_data.clear();
            m_uni->simulate(data, &sim_data);

            create_sim_models();

            m_sim_uni->learn(sim_data);
            m_sim_mix->learn(sim_data);

            float sim_likelihood_ratio = -2*m_sim_uni->loglikelihood(sim_data) + 2*m_sim_mix->get_loglikelihood();
            sum_ratio += sim_likelihood_ratio;
            ratios.push_back(sim_likelihood_ratio);
            if (sim_likelihood_ratio > log_likelihood_ratio)
                outliers++;
            if (samples == min_k) {
                //this is a breaking point
                if ((float)outliers / min_k > GOOD_PVAL/min_k) {
                    continue_sampling = false; //this will never be an interesting
                }
            }
            clear_sim_models();
        }
        float mean_ratio = sum_ratio / samples;
        float sd=0;
        for (float r : ratios) {
            sd += pow(r - mean_ratio,2);
            //		if (samples == K) cerr << r << " ";
        }
        sd = sqrt(sd / (samples-1));
        float z = (log_likelihood_ratio - mean_ratio)/sd;
        //		if (samples == K) cerr << endl;
        //compute standard deviation
        //compute p-value
        float p_val = (float)outliers / samples;
        out << data.get_data_id() << "\t" << p_val  << "\t" << z << "\t" << log_likelihood_ratio << "\t"
        << mean_ratio << "\t" << sd << "\t";
        m_mix->print_summary(out);
        m_uni->print_summary(out);
        out << endl;
        return(p_val);
    }
}

void MethylPatternMixtureModelScanner::clear_sim_models(){
	if (nullptr != m_sim_uni ){
		delete m_sim_uni;
		m_sim_uni = nullptr;
	}
	if (nullptr != m_sim_mix ){
		delete m_sim_mix;
		m_sim_mix = nullptr;
	}
}

void MethylPatternMixtureModelScanner::clear_models(){
	if (nullptr != m_uni ){
		delete m_sim_uni;
		m_uni = nullptr;
	}
	if (nullptr != m_mix ){
		delete m_mix;
		m_mix = nullptr;
	}
}
