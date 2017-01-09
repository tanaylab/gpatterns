/*
 * MethylPatternUniModelOneNoise.h
 *
 *  Created on: Jul 22, 2015
 *      Author: aviezerl
 */

#ifndef METHYLPATTERNUNIMODELONENOISE_H_
#define METHYLPATTERNUNIMODELONENOISE_H_

#include "MethylPatternUniModel.h"



class MethylPatternUniModelOneNoise : public MethylPatternUniModel {
public:
    MethylPatternUniModelOneNoise();
    MethylPatternUniModelOneNoise(const MethylPatternUniModelOneNoise& uni);
    virtual ~MethylPatternUniModelOneNoise();
    MethylPatternUniModelOneNoise& operator=(const MethylPatternUniModelOneNoise&);

    int 	learn(const ModelData&, const vector<float>& factors) override;
    void 	simulate_sample(const Sample* original, ModelData* data) override;
    MethylPatternSample simulate(const MethylPatternSample& original, const int& hamming_dist, const int& sample_num, ModelData* data) override;
    float 		get_log_prob(const Sample*) const override;

    float		get_noise() const {   return m_noise;  }
    void 		set_noise(const float& noise) {    m_noise = noise;   }

    void        init_params(const vector<int>& pattern, const ModelData& data);
    void 		init_params(const vector<int>& pattern, float noise);

protected:
    int 	learn_noise(const ModelData&);
    void 		print_parameters(ostream& out) 	const override;
    void		print_parameters_summary(ostream& out) const override;

    float		m_noise;			 //the probability for noise
};


#endif /* METHYLPATTERNUNIMODELONENOISE_H_ */
