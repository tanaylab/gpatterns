/*
 * MethylPatternMixtureModelGainLoss.h
 *
 *  Created on: Jul 23, 2015
 *      Author: aviezerl
 */

#ifndef METHYLPATTERNMIXTUREMODELGAINLOSS_H_
#define METHYLPATTERNMIXTUREMODELGAINLOSS_H_

#include "MethylPatternMixtureModel.h"



class MethylPatternMixtureModelGainLoss : public MethylPatternMixtureModel {
public:
    MethylPatternMixtureModelGainLoss();
    MethylPatternMixtureModelGainLoss(const MethylPatternMixtureModelGainLoss&);
    MethylPatternMixtureModelGainLoss(int number_of_mixtures, float uniform_prior);
    MethylPatternMixtureModelGainLoss(int number_of_mixtures, float uniform_prior, const int& K);
    virtual ~MethylPatternMixtureModelGainLoss();
    MethylPatternMixtureModelGainLoss& operator=(const MethylPatternMixtureModelGainLoss& m);
    MethylPatternMixtureModel* clone() override;
    void copy(MethylPatternMixtureModel* m) override;

    void init_models_hamming(const ModelData& data) override;

protected:
    void adjust_params(const int& c) override;
};


#endif /* METHYLPATTERNMIXTUREMODELGAINLOSS_H_ */
