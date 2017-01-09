/*
 * MethylPatternUniModel.h
 *
 *  Created on: Sep 16, 2014
 *      Author: aviezerl
 */

#ifndef METHYLPATTERNUNIMODEL_H_
#define METHYLPATTERNUNIMODEL_H_

#include "UniModel.h"
#include "MethylPatternData.h"
#include <vector>
using namespace std;

class MethylPatternUniModel: public UniModel {
public:
	MethylPatternUniModel();
	MethylPatternUniModel(const MethylPatternUniModel& uni);
	virtual ~MethylPatternUniModel();
    MethylPatternUniModel& operator=(const MethylPatternUniModel&);
	virtual int 	learn(const ModelData&) override;
	virtual int 	learn(const ModelData&, const vector<float>& factors) override = 0;
	float 	loglikelihood(const ModelData& data);
	virtual void 	simulate(const ModelData& original, ModelData* data) override;
	virtual void 	simulate_sample(const Sample* original, ModelData* data) override = 0;
	virtual MethylPatternSample simulate(const MethylPatternSample& original, const int& hamming_dist, const int& sample_num, ModelData* data) = 0;
	virtual void 		print(ostream& out) 	const override;
	virtual void		print_summary(ostream& out) const override;
	virtual float 		get_log_prob(const Sample*) const override = 0;
	vector<int> get_consensus_pattern() const { return m_consensus_pattern; }
	void set_consensus_pattern(const vector<int> pattern) {m_consensus_pattern = pattern;}



	void        init_params(const vector<int>& pattern, const ModelData& data);
protected:
	 virtual void 		print_parameters(ostream& out) 	const = 0;
     virtual void		print_parameters_summary(ostream& out) const = 0;

	vector<int> m_consensus_pattern; //holds the consensus pattern for the model: 0, 1 or -1 (for na)




};

#endif /* METHYLPATTERNUNIMODEL_H_ */
