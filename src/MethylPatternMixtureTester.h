/*
 * MethylPatternMixtureTester.h
 *
 *  Created on: Jan 7, 2015
 *      Author: aviezerl
 */

#ifndef METHYLPATTERNMIXTURETESTER_H_
#define METHYLPATTERNMIXTURETESTER_H_

#include "MethylPatternData.h"
#include "ModelDataMultiGroup.h"

class MethylPatternMixtureTester {
public:
    MethylPatternMixtureTester(const int& K = 1000, const float& unimodal_uniform_mix = 0.05, const int& num_of_inits = 3, const int& min_k=100);
    virtual ~MethylPatternMixtureTester();

    float get_mixture_pval(const MethylPatternData& data);
    float get_mixture_pval(const MethylPatternData& data, const vector<float>& mixing_probs);
    void remove_unimodal_groups(ModelDataMultiGroup& data, const float& qval_thresh, vector<int>& remove_inds, const bool& hard_remove);
    void remove_unimodal_groups(ModelDataMultiGroup& data, const float& qval_thresh, vector<int>& remove_inds, const vector<float>& mixing_probs, const bool& hard_remove);
    void remove_unimodal_groups(ModelDataMultiGroup& data, const float& qval_thresh, const bool& hard_remove);
    void remove_unimodal_groups(ModelDataMultiGroup& data, const float& qval_thresh, const vector<float>& mixing_probs, const bool& hard_remove);

//    void filter_centers(vector<vector<MethylPatternSample> >& centers, const vector<int>& remove_inds);
    template<typename T>
       void filter_by_inds(vector<T>& vec,  const vector<int>& remove_inds) {
           vector<T> temp_vec;
           if (0 != remove_inds.size()) {
               int r1=0;
               for (int r=0; r < vec.size(); r++) {
                   if (r != remove_inds[r1]) {
                       temp_vec.push_back(vec[r]);
                   } else {
                       r1++;
                   }
               }
               vec.swap(temp_vec);
           }
       }

protected:
    void adjust_pvals_and_remove(ModelDataMultiGroup& data, const float& qval_thresh, vector<int>& remove_inds, const vector<float>& pvals, const bool& hard_remove);

private:
    int m_K;
    float m_unimodal_uniform_mix;
    int m_num_of_inits;
    int m_min_k;
};


#endif /* METHYLPATTERNMIXTURETESTER_H_ */
