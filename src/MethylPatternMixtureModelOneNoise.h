/*
 * MethylPatternMixtureModelOneNoise.h
 *
 *  Created on: Jul 22, 2015
 *      Author: aviezerl
 */

#ifndef METHYLPATTERNMIXTUREMODELONENOISE_H_
#define METHYLPATTERNMIXTUREMODELONENOISE_H_

#include "MethylPatternMixtureModel.h"



class MethylPatternMixtureModelOneNoise : public MethylPatternMixtureModel {
public:
	MethylPatternMixtureModelOneNoise();
	MethylPatternMixtureModelOneNoise(const MethylPatternMixtureModelOneNoise&);
	MethylPatternMixtureModelOneNoise(int number_of_mixtures, float uniform_prior);
	MethylPatternMixtureModelOneNoise(int number_of_mixtures, float uniform_prior, const int& K);
	virtual ~MethylPatternMixtureModelOneNoise();
	MethylPatternMixtureModelOneNoise& operator=(const MethylPatternMixtureModelOneNoise& m);
	MethylPatternMixtureModel* clone() override;
	void copy(MethylPatternMixtureModel* m) override;

	void init_models_hamming(const ModelData& data) override;

protected:
	void adjust_params(const int& c) override;



};


#endif /* METHYLPATTERNMIXTUREMODELONENOISE_H_ */
