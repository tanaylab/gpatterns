/*
 * MethylPatternMixtureModelScanner.h
 *
 *  Created on: Aug 2, 2015
 *      Author: aviezerl
 */

#ifndef METHYLPATTERNMIXTUREMODELSCANNER_H_
#define METHYLPATTERNMIXTUREMODELSCANNER_H_

#include "MethylPatternMixtureModel.h"


class MethylPatternMixtureModelScanner {
public:
    MethylPatternMixtureModelScanner(const float& unimodal_uniform_mix, const int& num_of_inits);
    virtual ~MethylPatternMixtureModelScanner();

    float scan(const MethylPatternData& data, const int& K, const int& min_k, ostream& out);

    MethylPatternMixtureModel* get_uni(){return(m_uni);}
    MethylPatternMixtureModel* get_mix(){return(m_mix);}
protected:
    virtual void create_sim_models() = 0;
    virtual void clear_sim_models();
    virtual void create_models() = 0;
    virtual void clear_models();


protected:
    float m_unimodal_uniform_mix;
    int m_num_of_inits;
    MethylPatternMixtureModel* m_uni = nullptr;
    MethylPatternMixtureModel* m_mix = nullptr;
    MethylPatternMixtureModel* m_sim_uni = nullptr;
    MethylPatternMixtureModel* m_sim_mix = nullptr;
};


#endif /* METHYLPATTERNMIXTUREMODELSCANNER_H_ */
