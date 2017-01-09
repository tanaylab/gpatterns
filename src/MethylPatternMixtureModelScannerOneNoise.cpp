using namespace std;
/*
 * MethylPatternMixtureModelScannerOneNoise.cpp
 *
 *  Created on: Aug 2, 2015
 *      Author: aviezerl
 */

#include "MethylPatternMixtureModelScannerOneNoise.h"
#include "MethylPatternMixtureModelOneNoise.h"


MethylPatternMixtureModelScannerOneNoise::MethylPatternMixtureModelScannerOneNoise(const float& unimodal_uniform_mix, const int& num_of_inits): MethylPatternMixtureModelScanner(unimodal_uniform_mix, num_of_inits)  {
//	m_uni = new MethylPatternMixtureModelOneNoise(1, unimodal_uniform_mix, num_of_inits);
//	m_mix = new MethylPatternMixtureModelOneNoise(2, unimodal_uniform_mix, num_of_inits);
}

MethylPatternMixtureModelScannerOneNoise::~MethylPatternMixtureModelScannerOneNoise() {
    delete m_uni;
	delete m_mix;
}

void MethylPatternMixtureModelScannerOneNoise::create_sim_models(){
	clear_sim_models();
    m_sim_uni = new MethylPatternMixtureModelOneNoise(1, m_unimodal_uniform_mix, m_num_of_inits);
    m_sim_mix = new MethylPatternMixtureModelOneNoise(2, m_unimodal_uniform_mix, m_num_of_inits);
}

void MethylPatternMixtureModelScannerOneNoise::create_models(){
	clear_models();
    m_uni = new MethylPatternMixtureModelOneNoise(1, m_unimodal_uniform_mix, m_num_of_inits);
    m_mix = new MethylPatternMixtureModelOneNoise(2, m_unimodal_uniform_mix, m_num_of_inits);
}


