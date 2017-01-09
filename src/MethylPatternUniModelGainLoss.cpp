using namespace std;
/*
 * MethylPatternUniModelGainLoss.cpp
 *
 *  Created on: Jul 23, 2015
 *      Author: aviezerl
 */

//BASE_CC_FILE

#include "MethylPatternUniModelGainLoss.h"
#include "MethylPatternData.h"
#include "UniformModel.h"
#include <algorithm>
#include <chrono>
#include <random>
#include <string>
#include <numeric>
#include <unordered_set>
#include "SamplingUtil.h"

MethylPatternUniModelGainLoss::MethylPatternUniModelGainLoss():  MethylPatternUniModel(), m_gain(0), m_loss(0) {
    // TODO Auto-generated constructor stub
}

MethylPatternUniModelGainLoss::MethylPatternUniModelGainLoss(const MethylPatternUniModelGainLoss& uni) : MethylPatternUniModel(uni), m_gain(uni.m_gain), m_loss(uni.m_loss) {}

MethylPatternUniModelGainLoss::~MethylPatternUniModelGainLoss() {
    // TODO Auto-generated destructor stub
}

MethylPatternUniModelGainLoss& MethylPatternUniModelGainLoss::operator=(const MethylPatternUniModelGainLoss& src) {
    if (this == &src) {
        return *this;
    } else {
        MethylPatternUniModel::operator=(src);
        m_gain = src.m_gain;
        m_loss = src.m_loss;
    }
    return(*this);
}

int MethylPatternUniModelGainLoss::learn(const ModelData& data, const vector<float>& factors) {
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

    //compute gain and loss
    float one_error_bits=0;
    float one_total_bits=0;
    float zero_error_bits=0;
    float zero_total_bits=0;
    for (unsigned int index=0;index < m_consensus_pattern.size(); index++) {
        if (m_consensus_pattern[index] == 0) {
            zero_error_bits += n_one[index];
            zero_total_bits += n_zero[index] + n_one[index];
        }
        if (m_consensus_pattern[index] == 1) {
            one_error_bits += n_zero[index];
            one_total_bits += n_zero[index] + n_one[index];
        }
    }
    if (one_total_bits == 0) { //we have no samples for this model
        m_loss = 0;
    } else {
        m_loss = one_error_bits/one_total_bits;
    }
    if (zero_total_bits == 0) { //we have no samples for this model
        m_gain = 0;
    } else {
        m_gain = zero_error_bits/zero_total_bits;
    }
    float one_correct_bits = (one_total_bits - one_error_bits);
    float zero_correct_bits = (zero_total_bits - zero_error_bits);

    // calc_likelihood
    if ((m_gain == 0 && zero_error_bits > 0) || (m_loss == 0 && one_error_bits > 0) || (m_gain == 1 && zero_correct_bits > 0) || (m_loss == 1 && one_correct_bits > 0)) {
        m_loglikelihood = -std::numeric_limits<double>::max();
        return(0);
    }

    float l_prob_gain, l_prob_loss;
    if ((m_gain == 0 && zero_error_bits == 0) || (m_gain == 1 && zero_correct_bits == 0)) {
        l_prob_gain = 0;
    } else {
        l_prob_gain = zero_correct_bits * log(1-m_gain) + zero_error_bits * log(m_gain);
    }

    if ((m_loss == 0 && one_error_bits == 0) || (m_loss == 1 && one_correct_bits == 0)) {
        l_prob_loss = 0;
    } else {
        l_prob_loss = one_correct_bits * log(1-m_loss) + one_error_bits * log(m_loss);
    }

    m_loglikelihood = l_prob_gain + l_prob_loss;


    return(0);
}

int MethylPatternUniModelGainLoss::learn_gain_loss(const ModelData& data) {
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

    //compute gain and loss
    float one_error_bits=0;
    float one_total_bits=0;
    float zero_error_bits=0;
    float zero_total_bits=0;
    for (unsigned int index=0;index < m_consensus_pattern.size(); index++) {
        if (m_consensus_pattern[index] == 0) {
            zero_error_bits += n_one[index];
            zero_total_bits+= n_zero[index] + n_one[index];
        }
        if (m_consensus_pattern[index] == 1) {
            one_error_bits += n_zero[index];
            one_total_bits+= n_zero[index] + n_one[index];
        }
    }
    if (one_total_bits == 0) { //we have no samples for this model
        m_loss = 0;
    } else {
        m_loss = one_error_bits/one_total_bits;
    }
    if (zero_total_bits == 0) { //we have no samples for this model
        m_gain = 0;
    } else {
        m_gain = zero_error_bits/zero_total_bits;
    }

    float one_correct_bits = (one_total_bits - one_error_bits);
    float zero_correct_bits = (zero_total_bits - zero_error_bits);

    // calc_likelihood
    if ((m_gain == 0 && zero_error_bits > 0) || (m_loss == 0 && one_error_bits > 0) || (m_gain == 1 && zero_correct_bits > 0) || (m_loss == 1 && one_correct_bits > 0)) {
        m_loglikelihood = -std::numeric_limits<double>::max();
        return(0);
    }

    float l_prob_gain, l_prob_loss;
    if ((m_gain == 0 && zero_error_bits == 0) || (m_gain == 1 && zero_correct_bits == 0)) {
        l_prob_gain = 0;
    } else {
        l_prob_gain = zero_correct_bits * log(1-m_gain) + zero_error_bits * log(m_gain);
    }

    if ((m_loss == 0 && one_error_bits == 0) || (m_loss == 1 && one_correct_bits == 0)) {
        l_prob_loss = 0;
    } else {
        l_prob_loss = one_correct_bits * log(1-m_loss) + one_error_bits * log(m_loss);
    }

    m_loglikelihood = l_prob_gain + l_prob_loss;



    return(0);
}

void MethylPatternUniModelGainLoss::simulate_sample(const Sample* original, ModelData* data) {
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
                    if (r < m_gain) {
                        s += "1";
                    } else {
                        s += "0";
                    }
                } else {
                    if (r < m_loss) {
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

MethylPatternSample MethylPatternUniModelGainLoss::simulate(const MethylPatternSample& original,
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
                    if (r < m_gain) {
                        s += "1";
                    } else {
                        s += "0";
                    }
                } else if (center.pattern[i] == 1) {
                    if (r < m_loss) {
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

float MethylPatternUniModelGainLoss::get_log_prob(const Sample* s) const {
    int correct_zero=0;
    int errors_zero=0;
    int correct_one=0;
    int errors_one=0;
    MethylPatternSample* ms = (MethylPatternSample*) s;
    for (unsigned int i=0; i<ms->pattern.size(); i++) {
        if (0 == m_consensus_pattern[i]) {
            if (ms->pattern[i] == m_consensus_pattern[i]) {
                correct_zero++;
            } else if (ms->pattern[i] >= 0) {
                errors_zero++;
            }
        }
        if (1 == m_consensus_pattern[i]) {
            if (ms->pattern[i] == m_consensus_pattern[i]) {
                correct_one++;
            } else if (ms->pattern[i] >= 0) {
                errors_one++;
            }
        }
    }

    float l_prob_gain, l_prob_loss, l_prob;
    if ((m_gain == 0 && errors_zero > 0) || (m_loss == 0 && errors_one > 0) || (m_gain == 1 && correct_zero > 0) || (m_loss == 1 && correct_one > 0)) {
        l_prob = -std::numeric_limits<double>::max();
        return(l_prob);
    }

    if ((m_gain == 0 && errors_zero == 0) || (m_gain == 1 && correct_zero == 0)) {
        l_prob_gain = 0;
    } else {
        l_prob_gain = correct_zero * log(1-m_gain) + errors_zero * log(m_gain);
    }

    if ((m_loss == 0 && errors_one == 0) || (m_loss == 1 && correct_one == 0) ) {
        l_prob_loss = 0;
    } else {
        l_prob_loss = correct_one * log(1-m_loss) + errors_one * log(m_loss);
    }

    l_prob = l_prob_gain + l_prob_loss;
    return(l_prob);
}

void 		MethylPatternUniModelGainLoss::print_parameters(ostream& out) 	const {
    out << "\tgain: " << m_gain << "\tloss: " << m_loss;
}
void		MethylPatternUniModelGainLoss::print_parameters_summary(ostream& out) const {
    out << "\t" << m_gain << "\t" << m_loss;
}

void MethylPatternUniModelGainLoss::init_params(const vector<int>& pattern, const ModelData& data) {
    MethylPatternUniModel::init_params(pattern, data);
    learn_gain_loss(data);
}

void MethylPatternUniModelGainLoss::init_params(const vector<int>& pattern, float gain, float loss) {
    m_consensus_pattern = pattern;
    m_gain = gain;
    m_loss = loss;
}
