/*
 * MethylPatternMixtureModelScannerGainLoss.h
 *
 *  Created on: Aug 2, 2015
 *      Author: aviezerl
 */

#ifndef METHYLPATTERNMIXTUREMODELSCANNERGAINLOSS_H_
#define METHYLPATTERNMIXTUREMODELSCANNERGAINLOSS_H_

#include "MethylPatternMixtureModelScanner.h"



class MethylPatternMixtureModelScannerGainLoss : public MethylPatternMixtureModelScanner {
public:
    MethylPatternMixtureModelScannerGainLoss(const float& unimodal_uniform_mix, const int& num_of_inits);
    virtual ~MethylPatternMixtureModelScannerGainLoss();

protected:
    void create_sim_models() override;
    void create_models() override;
};


#endif /* METHYLPATTERNMIXTUREMODELSCANNERGAINLOSS_H_ */
