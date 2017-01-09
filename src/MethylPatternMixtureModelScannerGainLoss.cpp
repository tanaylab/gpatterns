using namespace std;
/*
 * MethylPatternMixtureModelScannerGainLoss.cpp
 *
 *  Created on: Aug 2, 2015
 *      Author: aviezerl
 */

#include "MethylPatternMixtureModelScannerGainLoss.h"
#include "MethylPatternMixtureModelGainLoss.h"


MethylPatternMixtureModelScannerGainLoss::MethylPatternMixtureModelScannerGainLoss(const float& unimodal_uniform_mix, const int& num_of_inits): MethylPatternMixtureModelScanner(unimodal_uniform_mix, num_of_inits) {
//	m_uni = new MethylPatternMixtureModelGainLoss(1, unimodal_uniform_mix, num_of_inits);
//    m_mix = new MethylPatternMixtureModelGainLoss(2, unimodal_uniform_mix, num_of_inits);

}

MethylPatternMixtureModelScannerGainLoss::~MethylPatternMixtureModelScannerGainLoss() {
    delete m_uni;
	delete m_mix;
}

void MethylPatternMixtureModelScannerGainLoss::create_sim_models(){
	clear_sim_models();
    m_sim_uni = new MethylPatternMixtureModelGainLoss(1, m_unimodal_uniform_mix, m_num_of_inits);
    m_sim_mix = new MethylPatternMixtureModelGainLoss(2, m_unimodal_uniform_mix, m_num_of_inits);
}

void MethylPatternMixtureModelScannerGainLoss::create_models(){
	clear_models();
    m_uni = new MethylPatternMixtureModelGainLoss(1, m_unimodal_uniform_mix, m_num_of_inits);
    m_mix = new MethylPatternMixtureModelGainLoss(2, m_unimodal_uniform_mix, m_num_of_inits);
}


