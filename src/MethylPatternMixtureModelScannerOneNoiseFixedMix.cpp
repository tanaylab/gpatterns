using namespace std;
/*
 * MethylPatternMixtureModelScannerOneNoiseFixedMix.cpp
 *
 *  Created on: Aug 2, 2015
 *      Author: aviezerl
 */

#include "MethylPatternMixtureModelScannerOneNoiseFixedMix.h"
#include "MethylPatternMixtureModelOneNoiseFixedMix.h"


MethylPatternMixtureModelScannerOneNoiseFixedMix::MethylPatternMixtureModelScannerOneNoiseFixedMix(const float& unimodal_uniform_mix, const int& num_of_inits, const vector<float>& mixing_probs): MethylPatternMixtureModelScanner(unimodal_uniform_mix, num_of_inits), m_mixing_probs(mixing_probs)  {
//	m_uni = new MethylPatternMixtureModelOneNoiseFixedMix(1, unimodal_uniform_mix, num_of_inits, mixing_probs);
//	m_mix = new MethylPatternMixtureModelOneNoiseFixedMix(2, unimodal_uniform_mix, num_of_inits, mixing_probs);
}

MethylPatternMixtureModelScannerOneNoiseFixedMix::~MethylPatternMixtureModelScannerOneNoiseFixedMix() {
    delete m_uni;
	delete m_mix;
}

void MethylPatternMixtureModelScannerOneNoiseFixedMix::create_sim_models(){
	clear_sim_models();
    m_sim_uni = new MethylPatternMixtureModelOneNoiseFixedMix(1, m_unimodal_uniform_mix, m_num_of_inits, m_mixing_probs);
    m_sim_mix = new MethylPatternMixtureModelOneNoiseFixedMix(2, m_unimodal_uniform_mix, m_num_of_inits, m_mixing_probs);
}

void MethylPatternMixtureModelScannerOneNoiseFixedMix::create_models(){
	clear_models();
    m_uni = new MethylPatternMixtureModelOneNoiseFixedMix(1, m_unimodal_uniform_mix, m_num_of_inits, m_mixing_probs);
    m_mix = new MethylPatternMixtureModelOneNoiseFixedMix(2, m_unimodal_uniform_mix, m_num_of_inits, m_mixing_probs);
}
