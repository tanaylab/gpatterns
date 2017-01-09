/*
 * MethylPatternMixtureModelScannerOneNoise.h
 *
 *  Created on: Aug 2, 2015
 *      Author: aviezerl
 */

#ifndef METHYLPATTERNMIXTUREMODELSCANNERONENOISE_H_
#define METHYLPATTERNMIXTUREMODELSCANNERONENOISE_H_

#include "MethylPatternMixtureModelScanner.h"



class MethylPatternMixtureModelScannerOneNoise : public MethylPatternMixtureModelScanner {
public:
    MethylPatternMixtureModelScannerOneNoise(const float& unimodal_uniform_mix, const int& num_of_inits);
    virtual ~MethylPatternMixtureModelScannerOneNoise();

protected:
    void create_sim_models() override;
    void create_models() override;
};


#endif /* METHYLPATTERNMIXTUREMODELSCANNERONENOISE_H_ */
