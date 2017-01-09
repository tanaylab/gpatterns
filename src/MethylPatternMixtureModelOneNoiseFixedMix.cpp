using namespace std;
/*
 * MethylPatternMixtureModelOneNoiseFixedMix.cpp
 *
 *  Created on: Jul 22, 2015
 *      Author: aviezerl
 */

#include "MethylPatternMixtureModelOneNoiseFixedMix.h"

MethylPatternMixtureModelOneNoiseFixedMix::MethylPatternMixtureModelOneNoiseFixedMix() {
    // TODO Auto-generated constructor stub

}

MethylPatternMixtureModelOneNoiseFixedMix::MethylPatternMixtureModelOneNoiseFixedMix(const MethylPatternMixtureModelOneNoiseFixedMix& m): MethylPatternMixtureModelOneNoise(m) {}

MethylPatternMixtureModelOneNoiseFixedMix::MethylPatternMixtureModelOneNoiseFixedMix(int number_of_mixtures, float uniform_prior, const vector<float> mixing_probs): MethylPatternMixtureModelOneNoise(number_of_mixtures, uniform_prior) {
    for (unsigned int c=0; c<mixing_probs.size(); c++) {
        m_mixing_probs[c]=mixing_probs[c];
    }
    m_mixing_probs[m_mixing_probs.size()-1] = uniform_prior;

}

MethylPatternMixtureModelOneNoiseFixedMix::MethylPatternMixtureModelOneNoiseFixedMix(int number_of_mixtures, float uniform_prior, const int& K, const vector<float> mixing_probs): MethylPatternMixtureModelOneNoise(number_of_mixtures, uniform_prior, K) {
    for (unsigned int c=0; c<mixing_probs.size(); c++) {
        m_mixing_probs[c]=mixing_probs[c];
    }
    m_mixing_probs[m_mixing_probs.size()-1] = uniform_prior;
}



MethylPatternMixtureModelOneNoiseFixedMix::~MethylPatternMixtureModelOneNoiseFixedMix() {
    // TODO Auto-generated destructor stub
}

float MethylPatternMixtureModelOneNoiseFixedMix::m_step(const ModelData& data) {
    //	cerr <<"M-STEP" << endl;
    for (unsigned int c=0; c<m_models.size(); c++) {
        m_models[c]->learn(data, m_mixture_posteriors[c]);
        //		m_models[c]->print(cerr);
    }
    return(1);
}



