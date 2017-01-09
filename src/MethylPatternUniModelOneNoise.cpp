using namespace std;
/*
 * MethylPatternUniModelOneNoise.cpp
 *
 *  Created on: Jul 22, 2015
 *      Author: aviezerl
 */

//BASE_CC_FILE

#include "MethylPatternUniModelOneNoise.h"
#include "MethylPatternData.h"
#include "UniformModel.h"
#include <algorithm>
#include <chrono>
#include <random>
#include <string>
#include <numeric>
#include <unordered_set>
#include "SamplingUtil.h"

MethylPatternUniModelOneNoise::MethylPatternUniModelOneNoise(): MethylPatternUniModel(),
m_noise(0) {
    // TODO Auto-generated constructor stub
}

MethylPatternUniModelOneNoise::MethylPatternUniModelOneNoise(const MethylPatternUniModelOneNoise& uni) : MethylPatternUniModel(uni), m_noise(uni.m_noise) {
}

MethylPatternUniModelOneNoise::~MethylPatternUniModelOneNoise() {
    // TODO Auto-generated destructor stub
}

MethylPatternUniModelOneNoise& MethylPatternUniModelOneNoise::operator=(const MethylPatternUniModelOneNoise& src) {
    if (this == &src) {
        return *this;
    } else {
    	MethylPatternUniModel::operator=(src);
        m_noise= src.m_noise;
    }
    return(*this);
}

int MethylPatternUniModelOneNoise::learn(const ModelData& data, const vector<float>& factors) {
    m_consensus_pattern.resize(0);
    //go over samples from data to determine consensus pattern
    int numberOfSamples = data.get_number_of_samples();
    vector<float> n_zero;
    vector<float> n_one;
    for (int i=0; i<numberOfSamples; i++) {
        MethylPatternSample* s = (MethylPatternSample*)data.get_sample(i);
        if (!s) {
            cerr << "Received Sample NULL" << endl;
            return (-1);
        }
        if (n_zero.size() == 0) {
            n_zero.resize(s->pattern.size(),0);
            n_one.resize(s->pattern.size(),0);
            m_consensus_pattern.resize(s->pattern.size(), -1);
        }

        for (unsigned int a=0; a<s->pattern.size(); a++) {
            if (s->pattern[a] == 0)
                n_zero[a] += factors[i];
            if (s->pattern[a] == 1)
                n_one[a] += factors[i];
        }
    }
    for (unsigned int i=0; i<n_zero.size(); i++) {
        if (n_zero[i] + n_one[i] > 0) {
            if (n_zero[i] >= n_one[i]) {
                m_consensus_pattern[i] = 0;
            } else {
                m_consensus_pattern[i] = 1;
            }
        }
    }
    //compute the noise
    float error_bits=0;
    float total_bits=0;
    for (unsigned int index=0;index < m_consensus_pattern.size(); index++) {
        if (m_consensus_pattern[index] == 0) {
            error_bits += n_one[index];
        }
        if (m_consensus_pattern[index] == 1) {
            error_bits += n_zero[index];
        }
        total_bits+= n_zero[index] + n_one[index];
    }
    if (total_bits == 0) { //we have no samples for this model
        m_noise = 0;
    } else {
        m_noise = error_bits/total_bits;
    }

    // calc_likelihood
    if (m_noise == 0) {
        m_loglikelihood = 0;
    } else {
        m_loglikelihood = (total_bits - error_bits)*log(1-m_noise) +
                          error_bits * log(m_noise);
    }

    return(0);


}

int MethylPatternUniModelOneNoise::learn_noise(const ModelData& data) {

//count number of ones and zeros in each position
	int numberOfSamples = data.get_number_of_samples();
	    vector<float> n_zero;
	    vector<float> n_one;
	    for (int i=0; i<numberOfSamples; i++) {
	        MethylPatternSample* s = (MethylPatternSample*)data.get_sample(i);
	        if (!s) {
	            cerr << "Received Sample NULL" << endl;
	            return (-1);
	        }
	        if (n_zero.size() == 0) {
	            n_zero.resize(s->pattern.size(),0);
	            n_one.resize(s->pattern.size(),0);
	        }

	        for (unsigned int a=0; a<s->pattern.size(); a++) {
	            if (s->pattern[a] == 0)
	                n_zero[a]++;
	            if (s->pattern[a] == 1)
	                n_one[a]++;
	        }
	    }



	//calc noise
    float error_bits=0;
    float total_bits=0;
    for (unsigned int index=0;index < m_consensus_pattern.size(); index++) {
        if (m_consensus_pattern[index] == 0) {
            error_bits += n_one[index];
        }
        if (m_consensus_pattern[index] == 1) {
            error_bits += n_zero[index];
        }
        total_bits+= n_zero[index] + n_one[index];
    }
    if (total_bits == 0) { //we have no samples for this model
        m_noise = 0;
    } else {
        m_noise = error_bits/total_bits;
    }

    // calc_likelihood
    if (m_noise == 0) {
        m_loglikelihood = 0;
    } else {
        m_loglikelihood = (total_bits - error_bits)*log(1-m_noise) +
                          error_bits * log(m_noise);
    }
    return(0);
}

void MethylPatternUniModelOneNoise::simulate_sample(const Sample* original, ModelData* data) {
    float r;
    const MethylPatternSample* os = (MethylPatternSample*)original;
    string s("");
    for (unsigned int i=0; i<os->pattern.size(); i++) {
        if (os->pattern[i] < 0)
            s += "*";
        else {
            if (m_consensus_pattern[i] < 0) {
                s += "*";
            } else {
                r = ((float)rand())/RAND_MAX;
                if (m_consensus_pattern[i] == 0) {
                    if (r < m_noise) {
                        s += "1";
                    } else {
                        s += "0";
                    }
                } else {
                    if (r < m_noise) {
                        s += "0";
                    } else {
                        s += "1";
                    }
                }
            }
        }
    }
    ((MethylPatternData*)data)->add_sample(s);
}

MethylPatternSample MethylPatternUniModelOneNoise::simulate(const MethylPatternSample& original,
        const int& hamming_dist, const int& sample_num, ModelData* data) {
    float r;
    MethylPatternSample center(original, hamming_dist);
    for (int i = 0; i < sample_num; i++) {
        string s("");
        for (unsigned int i = 0; i < center.pattern.size(); i++) {
            if (center.pattern[i] < 0)
                s += "*";
            else {
                r = ((float) rand()) / RAND_MAX;
                if (center.pattern[i] == 0) {
                    if (r < m_noise) {
                        s += "1";
                    } else {
                        s += "0";
                    }
                } else if (center.pattern[i] == 1) {
                    if (r < m_noise) {
                        s += "0";
                    } else {
                        s += "1";
                    }
                }
            }
        }
        ((MethylPatternData*) data)->add_sample(s);
    }
    return(center);
}

float MethylPatternUniModelOneNoise::get_log_prob(const Sample* s) const {
    int correct=0;
    int errors=0;
    MethylPatternSample* ms = (MethylPatternSample*) s;
    for (unsigned int i=0; i<ms->pattern.size(); i++) {
        if (ms->pattern[i] == m_consensus_pattern[i]) {
            if (m_consensus_pattern[i] >= 0)
                correct++;
        } else {
            if (ms->pattern[i] >= 0)
                errors++;
        }
    }
    float l_prob;
    if (m_noise == 0) {
        if (errors > 0) {
            l_prob = -std::numeric_limits<double>::max();
        } else {
            l_prob = 0;
        }
    } else {
        l_prob = correct * log(1-m_noise)
                 + errors * log(m_noise);
    }
    return(l_prob);
}

void 		MethylPatternUniModelOneNoise::print_parameters(ostream& out) 	const {
	out << "\tnoise: " << m_noise;
}
void		MethylPatternUniModelOneNoise::print_parameters_summary(ostream& out) const{
    out << "\t" << m_noise;
}

void MethylPatternUniModelOneNoise::init_params(const vector<int>& pattern, const ModelData& data) {
	MethylPatternUniModel::init_params(pattern, data);
    learn_noise(data);
}

void MethylPatternUniModelOneNoise::init_params(const vector<int>& pattern, float noise) {
    m_consensus_pattern = pattern;
    m_noise = noise;
}


