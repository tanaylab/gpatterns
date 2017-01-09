/*
 * MethylPatternMixtureModelOneNoiseFixedMix.h
 *
 *  Created on: Jul 22, 2015
 *      Author: aviezerl
 */

#ifndef METHYLPATTERNMIXTUREMODELONENOISEFIXEDMIX_H_
#define METHYLPATTERNMIXTUREMODELONENOISEFIXEDMIX_H_

#include "MethylPatternMixtureModelOneNoise.h"



class MethylPatternMixtureModelOneNoiseFixedMix : public MethylPatternMixtureModelOneNoise {
public:
    MethylPatternMixtureModelOneNoiseFixedMix();
    MethylPatternMixtureModelOneNoiseFixedMix(const MethylPatternMixtureModelOneNoiseFixedMix&);
        MethylPatternMixtureModelOneNoiseFixedMix(int number_of_mixtures, float uniform_prior, const vector<float> mixing_probs);
        MethylPatternMixtureModelOneNoiseFixedMix(int number_of_mixtures, float uniform_prior, const int& K, const vector<float> mixing_probs);
        virtual ~MethylPatternMixtureModelOneNoiseFixedMix();


        void init_models(const ModelData& data) override {}
        ;
        float m_step(const ModelData& data) override;

};


#endif /* METHYLPATTERNMIXTUREMODELONENOISEFIXEDMIX_H_ */
