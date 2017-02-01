using namespace std;
/*
 * MethylPatternMixtureTester.cpp
 *
 *  Created on: Jan 7, 2015
 *      Author: aviezerl
 */

#include "MethylPatternMixtureTester.h"
#include "MethylPatternUniModel.h"
#include "MethylPatternMixtureModelOneNoise.h"
#include "MethylPatternMixtureModelOneNoiseFixedMix.h"
#include "SamplingUtil.h"
#define GOOD_PVAL 0.0

MethylPatternMixtureTester::MethylPatternMixtureTester(const int& K, const float& unimodal_uniform_mix, const int& num_of_inits, const int& min_k): m_K(K), m_unimodal_uniform_mix(unimodal_uniform_mix), m_num_of_inits(num_of_inits), m_min_k(min_k)  {
    // TODO Auto-generated constructor stub

}

MethylPatternMixtureTester::~MethylPatternMixtureTester() {
    // TODO Auto-generated destructor stub
}

void MethylPatternMixtureTester::remove_unimodal_groups(ModelDataMultiGroup& data, const float& qval_thresh, vector<int>& remove_inds, const bool& hard_remove) {
    vector<float> pvals;
    for (size_t i=0; i < data.get_number_of_groups(); i++) {
        pvals.push_back( get_mixture_pval(*((MethylPatternData*)data.get_group(i))) );
    }
    adjust_pvals_and_remove(data, qval_thresh, remove_inds, pvals, hard_remove);
}

void MethylPatternMixtureTester::remove_unimodal_groups(ModelDataMultiGroup& data, const float& qval_thresh, vector<int>& remove_inds, const vector<float>& mixing_probs, const bool& hard_remove) {
    vector<float> pvals;
    for (size_t i=0; i < data.get_number_of_groups(); i++) {
        pvals.push_back( get_mixture_pval(*((MethylPatternData*)data.get_group(i)), mixing_probs) );
    }
    adjust_pvals_and_remove(data, qval_thresh, remove_inds, pvals, hard_remove);
}

void MethylPatternMixtureTester::adjust_pvals_and_remove(ModelDataMultiGroup& data, const float& qval_thresh, vector<int>& remove_inds, const vector<float>& pvals, const bool& hard_remove) {
    vector<float> qvals = p_adjust(pvals);
    for (size_t i=0; i < qvals.size(); i++) {
        if (qvals[i] > qval_thresh) {
            remove_inds.push_back(i);
        }
    }
    data.remove_groups(remove_inds, hard_remove);

}


void MethylPatternMixtureTester::remove_unimodal_groups(ModelDataMultiGroup& data, const float& qval_thresh, const vector<float>& mixing_probs, const bool& hard_remove) {
    vector<int> remove_inds;
    remove_unimodal_groups(data, qval_thresh, remove_inds, mixing_probs, hard_remove);
}

void MethylPatternMixtureTester::remove_unimodal_groups(ModelDataMultiGroup& data, const float& qval_thresh, const bool& hard_remove) {
    vector<int> remove_inds;
    remove_unimodal_groups(data, qval_thresh, remove_inds, hard_remove);
}

float MethylPatternMixtureTester::get_mixture_pval(const MethylPatternData& data) {
    MethylPatternData sim_data;
    //learn unimodel
    MethylPatternMixtureModelOneNoise uni(1, m_unimodal_uniform_mix, m_num_of_inits);
    uni.learn(data);

    //learn mixture model
    MethylPatternMixtureModelOneNoise mix(2, m_unimodal_uniform_mix, m_num_of_inits);
    mix.learn(data);

    float log_likelihood_ratio = -2*uni.get_loglikelihood() + 2 * mix.get_loglikelihood();
    if (log_likelihood_ratio == 0) { //this means that the mixture model did not add anything !
        //there is no need to run any simulations as the best we can do is with the uni model
        return(1);
    } else {
        int outliers=0;
        //sample likelihood ratio from unimodel sampled data and compute likelihood ratio p-value
        int samples=0;
        bool continue_sampling=true;
        float sum_ratio=0;
        vector<float> ratios;

        while (continue_sampling && samples < m_K) {
            samples++;
            sim_data.clear();
            uni.simulate(data, &sim_data);
            MethylPatternMixtureModelOneNoise sim_uni(1, m_unimodal_uniform_mix, m_num_of_inits);
            MethylPatternMixtureModelOneNoise sim_mix(2, m_unimodal_uniform_mix, m_num_of_inits);

            sim_uni.learn(sim_data);
            sim_mix.learn(sim_data);
            float sim_likelihood_ratio = -2*sim_uni.loglikelihood(sim_data) + 2*sim_mix.get_loglikelihood();
            sum_ratio += sim_likelihood_ratio;
            ratios.push_back(sim_likelihood_ratio);
            if (sim_likelihood_ratio > log_likelihood_ratio)
                outliers++;
            if (samples == m_min_k) {
                //this is a breaking point
                if ((float)outliers / m_min_k > GOOD_PVAL/m_min_k) {
                    continue_sampling = false; //this will never be an interesting
                }
            }
        }
        float mean_ratio = sum_ratio / samples;

        //compute standard deviation
        float sd=0;
        for (float r : ratios) {
            sd += pow(r - mean_ratio,2);
            //		if (samples == K) cerr << r << " ";
        }
        sd = sqrt(sd / (samples-1));
        //        float z = (log_likelihood_ratio - mean_ratio)/sd;


        //compute p-value
        float p_val = (float)outliers / samples;
        return(p_val);
    }
}

float MethylPatternMixtureTester::get_mixture_pval(const MethylPatternData& data, const vector<float>& mixing_probs) {
    MethylPatternData sim_data;
    //learn unimodel
    MethylPatternMixtureModelOneNoiseFixedMix uni(1, m_unimodal_uniform_mix, m_num_of_inits, mixing_probs);
    uni.learn(data);

    //learn mixture model
    MethylPatternMixtureModelOneNoiseFixedMix mix(2, m_unimodal_uniform_mix, m_num_of_inits, mixing_probs);
    mix.learn(data);

    float log_likelihood_ratio = -2*uni.get_loglikelihood() + 2 * mix.get_loglikelihood();
    if (log_likelihood_ratio == 0) { //this means that the mixture model did not add anything !
        //there is no need to run any simulations as the best we can do is with the uni model
        return(1);
    } else {
        int outliers=0;
        //sample likelihood ratio from unimodel sampled data and compute likelihood ratio p-value
        int samples=0;
        bool continue_sampling=true;
        float sum_ratio=0;
        vector<float> ratios;

        while (continue_sampling && samples < m_K) {
            samples++;
            sim_data.clear();
            uni.simulate(data, &sim_data);
            MethylPatternMixtureModelOneNoiseFixedMix sim_uni(1, m_unimodal_uniform_mix, m_num_of_inits, mixing_probs);
            MethylPatternMixtureModelOneNoiseFixedMix sim_mix(2, m_unimodal_uniform_mix, m_num_of_inits, mixing_probs);

            sim_uni.learn(sim_data);
            sim_mix.learn(sim_data);
            float sim_likelihood_ratio = -2*sim_uni.loglikelihood(sim_data) + 2*sim_mix.get_loglikelihood();
            sum_ratio += sim_likelihood_ratio;
            ratios.push_back(sim_likelihood_ratio);
            if (sim_likelihood_ratio > log_likelihood_ratio)
                outliers++;
            if (samples == m_min_k) {
                //this is a breaking point
                if ((float)outliers / m_min_k > GOOD_PVAL/m_min_k) {
                    continue_sampling = false; //this will never be an interesting
                }
            }
        }
        float mean_ratio = sum_ratio / samples;

        //compute standard deviation
        float sd=0;
        for (float r : ratios) {
            sd += pow(r - mean_ratio,2);
            //		if (samples == K) cerr << r << " ";
        }
        sd = sqrt(sd / (samples-1));
        //        float z = (log_likelihood_ratio - mean_ratio)/sd;


        //compute p-value
        float p_val = (float)outliers / samples;
        return(p_val);
    }
}

//void MethylPatternMixtureTester::filter_centers(vector<vector<MethylPatternSample> >& centers, const vector<int>& remove_inds) {
//    vector<vector<MethylPatternSample> > temp_centers;
//    if (0 != remove_inds.size()) {
//        int r1=0;
//        for (int r=0; r < centers.size(); r++) {
//            if (r != remove_inds[r1]) {
//                temp_centers.push_back(centers[r]);
//            } else {
//                r1++;
//            }
//        }
//        centers.swap(temp_centers);
//    }
//}


