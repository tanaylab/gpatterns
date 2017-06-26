using namespace std;
/*
 * MethylPatternMixtureModel.cpp
 *
 *  Created on: Sep 21, 2014
 *      Author: aviezerl
 */

//BASE_CC_FILE

#include "MethylPatternMixtureModel.h"
#include "MethylPatternUniModel.h"
#include "MethylPatternData.h"
#include "UniformModel.h"
#include <algorithm>
#include <chrono>
#include <random>
#include <string>
#include <numeric>
#include <unordered_set>
#include "SamplingUtil.h"
#include "util.h"


MethylPatternMixtureModel::MethylPatternMixtureModel(): m_uniform_prior(0), m_K(1) {}

MethylPatternMixtureModel::MethylPatternMixtureModel(const MethylPatternMixtureModel& m): MixtureModel(m),
m_uniform_prior(m.m_uniform_prior), m_K(m.m_K)  {
    m_mixing_probs = m.m_mixing_probs;
    m_mixture_posteriors = m.m_mixture_posteriors;
}

void MethylPatternMixtureModel::swap(MethylPatternMixtureModel& m) {
    std::swap(m_uniform_prior, m.m_uniform_prior);
    std::swap(m_delta_likelihood_converged, m.m_delta_likelihood_converged);
    std::swap(m_number_of_models, m.m_number_of_models);
    std::swap(m_K, m.m_K);
    std::swap(m_models, m.m_models);
    std::swap(m_loglikelihood, m.m_loglikelihood);
    std::swap(m_mixing_probs, m.m_mixing_probs);
    std::swap(m_mixture_posteriors, m.m_mixture_posteriors);
    //	std::swap(m_uniform, m.m_uniform);

}

MethylPatternMixtureModel::MethylPatternMixtureModel(int number_of_mixtures, float uniform_prior): MethylPatternMixtureModel(number_of_mixtures, uniform_prior, 1) {}

MethylPatternMixtureModel::MethylPatternMixtureModel(int number_of_mixtures, float uniform_prior, const int& K):
m_uniform_prior(uniform_prior), m_K(K) {
}

MethylPatternMixtureModel::~MethylPatternMixtureModel() {
    vector<UniModel*>::iterator cur_m = m_models.begin();
    vector<UniModel*>::iterator end_m = m_models.end();
    while (cur_m != end_m) {
        delete(*cur_m);
        cur_m++;
    }
}


//initialize the models by sampling two patterns that where not already chosen (not in chosen_inds) from uniq_samples.
void MethylPatternMixtureModel::init_models_comb(const vector<Sample*>& uniq_samples, const vector<int>& uniq_sample_freqs, const ModelData& data, unordered_set<vector<int>, vector_hasher >& chosen_inds) {
    random_device rd;
    int k = (int)chosen_inds.size();

    vector<int> chosen_inds_i;

    int init_tries = 0;
    int max_tries = k * 10;
    while((int)chosen_inds.size() <= k) { //until found a sample for each model
        //TODO
        chosen_inds_i.clear();

        //create a vector with [1..number of samples]
        vector<int> freqs = uniq_sample_freqs;
        vector<int> sample_inds(freqs.size());
        iota(sample_inds.begin(), sample_inds.end(), 0);


        //select a random sample for each model
        default_random_engine generator( rd() );
        for (size_t i=0; (i <  m_number_of_models - 1) && (freqs.size() > 0); i++) {
            discrete_distribution<int> d(freqs.begin(), freqs.end());
            int c = d(generator);

            std::swap(freqs[c], freqs.back());
            freqs.pop_back();
            std::swap(sample_inds[c], sample_inds.back());

            chosen_inds_i.push_back(sample_inds.back());
            sample_inds.pop_back();
        }

        //try to put the chosen samples in the hash (if already chosen, would not change its size)
        chosen_inds.insert(chosen_inds_i);
        init_tries++;
        if (init_tries == max_tries) {
            break;
        }
    }


    //initialize the models with the chosen samples
    for (size_t c=0; c< chosen_inds_i.size(); c++) {
        ((MethylPatternUniModel*)(m_models[c]))->init_params(((MethylPatternSample*)uniq_samples[ chosen_inds_i[c] ])->pattern, data);
        adjust_params(c);
    }

    //if there are more models than samples assign the first chosen sample to the rest of the samples
    for (size_t c=chosen_inds_i.size(); c< m_number_of_models - 1; c++) {
        ((MethylPatternUniModel*)(m_models[c]))->init_params(((MethylPatternSample*)uniq_samples[ chosen_inds_i[0] ])->pattern, data);
        adjust_params(c);
    }
}

int MethylPatternMixtureModel::learn(const ModelData& data) {
    // create a vector of all the unique samples
    vector<int> uniq_sample_freqs;
    vector<Sample*> uniq_samples = data.unique_samples(uniq_sample_freqs);

    //an unordered set (hash) to store the indices of samples that where already used for initialization
    unordered_set<vector<int>, vector_hasher > chosen_inds;

    //other variable declarations
    vector<MethylPatternMixtureModel*> res_models;
    res_models.reserve(0.5 * uniq_samples.size() * (uniq_samples.size() - 1));
    int init_ind = 0, max_lhood_ind = 0;
    double max_loglikelihood = -std::numeric_limits<double>::max();

    // first init - hamming distance heuristic
    init_models_hamming(data);
    MixtureModel::learn(data);
    res_models.push_back(clone());

    max_loglikelihood = m_loglikelihood;
    max_lhood_ind = init_ind;
    init_ind++;

    //number of initializations is m_K unless there aren't enough unique samples, which can happen when:
    //a. number of unique samples < m_K
    //b. number of combinations of samples < m_K
    size_t num_of_inits = (m_number_of_models - 1) > uniq_samples.size() ? min(m_K, (int)uniq_samples.size()) : min(m_K, nchoosek(uniq_samples.size(), m_number_of_models - 1));

    // on each initialization we choose two patterns, learn the model and save it in res_models
    for (size_t k=0; k<num_of_inits; k++) {
        init_models_comb(uniq_samples, uniq_sample_freqs, data, chosen_inds);
        MixtureModel::learn(data);
        res_models.push_back(clone());
        if (m_loglikelihood > max_loglikelihood) {
            max_loglikelihood = m_loglikelihood;
            max_lhood_ind = init_ind;
        }
        init_ind++;
    }

    //copy the model with the highest likelihood to the current MethylPatternMixtureModel object
    copy(res_models[max_lhood_ind]);

    //delete other models
    vector<MethylPatternMixtureModel*>::iterator cur_m = res_models.begin();
    vector<MethylPatternMixtureModel*>::iterator end_m = res_models.end();
      while (cur_m != end_m) {
          delete(*cur_m);
          cur_m++;
      }

    return(0);
}

void MethylPatternMixtureModel::init_models(const ModelData& data) {
    //	associate each sample with one of the centers to initialize the mixing probs (pai)
    vector<int> associated_samples(m_models.size(),0);
    for (int i=0; i<data.get_number_of_samples(); i++) {
        //	        int minimal_distance = ((MethylPatternSample*)data.get_sample(0))->pattern.size();
        int preferred_model=0;
        MethylPatternSample* s = (MethylPatternSample*) data.get_sample(i);
        int minimal_distance = s->pattern.size();
        for (size_t c=0; c<m_models.size()-1; c++) {
            int dist = s->hamming_distance(((MethylPatternUniModel*)m_models[c])->get_consensus_pattern());
            if ( dist < minimal_distance) {
                minimal_distance = dist;
                preferred_model = c;
            }
        }
        associated_samples[preferred_model]++;
    }
    for (size_t c=0; c<m_models.size()-1; c++) {
        m_mixing_probs[c] = ((float) associated_samples[c] / data.get_number_of_samples()) * (1-m_uniform_prior);
    }
    m_mixing_probs[m_models.size()-1] = m_uniform_prior;
}








float MethylPatternMixtureModel::m_step(const ModelData& data) {
    MixtureModel::m_step(data);
    float uniform_mix = m_mixing_probs[m_mixing_probs.size()-1];
    for (size_t c=0; c<m_models.size()-1; c++) {
        m_mixing_probs[c] = m_mixing_probs[c] / (1-uniform_mix) * (1-m_uniform_prior);
    }
    m_mixing_probs[m_models.size()-1] = m_uniform_prior;
    return(1);
}

float MethylPatternMixtureModel::center_error(const MethylPatternSample& center, const int& model) {
    return(center.hamming_distance( ((MethylPatternUniModel*)(m_models[model]))->get_consensus_pattern()) / center.pattern.size());
}

vector<float> MethylPatternMixtureModel::center_error(const vector<MethylPatternSample>& centers) {
    vector<float> errors;
    for (size_t i=0; i < m_models.size() - 1; i++) {
        errors.push_back(center_error(centers[i], i));
    }
    return(errors);
}
