/*
 * MethylPatternMixtureModelScannerOneNoiseFixedMix.h
 *
 *  Created on: Aug 2, 2015
 *      Author: aviezerl
 */

#ifndef METHYLPATTERNMIXTUREMODELSCANNERONENOISEFIXEDMIX_H_
#define METHYLPATTERNMIXTUREMODELSCANNERONENOISEFIXEDMIX_H_

#include "MethylPatternMixtureModelScanner.h"



class MethylPatternMixtureModelScannerOneNoiseFixedMix : public MethylPatternMixtureModelScanner {
public:
    MethylPatternMixtureModelScannerOneNoiseFixedMix(const float& unimodal_uniform_mix, const int& num_of_inits, const vector<float>& mixing_probs);
    virtual ~MethylPatternMixtureModelScannerOneNoiseFixedMix();

protected:
    void create_sim_models() override;
    void create_models() override;

protected:
    vector <float> m_mixing_probs;
};


#endif /* METHYLPATTERNMIXTUREMODELSCANNERONENOISEFIXEDMIX_H_ */
