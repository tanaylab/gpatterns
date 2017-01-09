using namespace std;
/*
 * UniModel.cpp
 *
 *  Created on: Sep 16, 2014
 *      Author: aviezerl
 */


//BASE_CC_FILE


#include "UniModel.h"

UniModel::UniModel():m_loglikelihood(0) {}

UniModel::~UniModel() {}

float UniModel::loglikelihood(const ModelData& data) {
    m_loglikelihood=0;
    for (int i=0; i<data.get_number_of_samples(); i++) {
        m_loglikelihood += get_log_prob(data.get_sample(i));
    }
    return(m_loglikelihood);
}
