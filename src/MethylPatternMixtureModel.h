/*
 * MethylPatternMixtureModel.h
 *
 *  Created on: Sep 21, 2014
 *      Author: aviezerl
 */

#ifndef METHYLPATTERNMIXTUREMODEL_H_
#define METHYLPATTERNMIXTUREMODEL_H_

#include "MixtureModel.h"
#include "MethylPatternData.h"

struct vector_hasher {
    template <typename T>
    size_t operator()(const T p) const {
        return boost::hash_range(p.begin(), p.end());
    }

};

class MethylPatternMixtureModel: public MixtureModel {
public:
    MethylPatternMixtureModel();
    MethylPatternMixtureModel(const MethylPatternMixtureModel&);
    MethylPatternMixtureModel(int number_of_mixtures, float uniform_prior);
    MethylPatternMixtureModel(int number_of_mixtures, float uniform_prior, const int& K);
    virtual ~MethylPatternMixtureModel();
    void swap(MethylPatternMixtureModel& m);

    virtual MethylPatternMixtureModel* clone() = 0;
    virtual void copy(MethylPatternMixtureModel* m) = 0;

    int learn(const ModelData& data);

    virtual void init_models_hamming(const ModelData& data) = 0;
    void init_models_comb(const vector<Sample*>& uniq_samples, const vector<int>& uniq_sample_freqs, const ModelData& data, unordered_set<vector<int>, vector_hasher >& chosen_inds);

    vector<float> center_error(const vector<MethylPatternSample>& centers);

protected:
    void init_models(const ModelData& data) override;

    float m_step(const ModelData& data) override;


    virtual void adjust_params(const int& c) = 0;

    float center_error(const MethylPatternSample& center, const int& model);

protected:
    float m_uniform_prior;
    int m_K; // number of initializations
};

#endif /* METHYLPATTERNMIXTUREMODEL_H_ */
