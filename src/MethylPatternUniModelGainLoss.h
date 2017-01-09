/*
 * MethylPatternUniModelGainLoss.h
 *
 *  Created on: Jul 23, 2015
 *      Author: aviezerl
 */

#ifndef METHYLPATTERNUNIMODELGAINLOSS_H_
#define METHYLPATTERNUNIMODELGAINLOSS_H_

#include "MethylPatternUniModel.h"



class MethylPatternUniModelGainLoss : public MethylPatternUniModel {
public:
    MethylPatternUniModelGainLoss();
    MethylPatternUniModelGainLoss(const MethylPatternUniModelGainLoss& uni);
    virtual ~MethylPatternUniModelGainLoss();
    MethylPatternUniModelGainLoss& operator=(const MethylPatternUniModelGainLoss&);

    int 	learn(const ModelData&, const vector<float>& factors) override;
    void 	simulate_sample(const Sample* original, ModelData* data) override;
    MethylPatternSample simulate(const MethylPatternSample& original, const int& hamming_dist, const int& sample_num, ModelData* data) override;
    float 		get_log_prob(const Sample*) const override;

    float		get_gain() const {     return m_gain;    }
    void 		set_gain(const float& gain) {        m_gain = gain;    }
    float		get_loss() const {     return m_loss;    }
    void 		set_loss(const float& loss) {        m_loss = loss;    }

    void        init_params(const vector<int>& pattern, const ModelData& data);
    void 		init_params(const vector<int>& pattern, float gain, float loss);

protected:
    int 	learn_gain_loss(const ModelData&);
    void 		print_parameters(ostream& out) 	const override;
    void		print_parameters_summary(ostream& out) const override;

    float		m_gain;			 //the probability for gaining methylation
    float		m_loss;			 //the probability for losing methylation

};


#endif /* METHYLPATTERNUNIMODELGAINLOSS_H_ */
