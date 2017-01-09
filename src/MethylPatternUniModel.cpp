using namespace std;
/*
 * MethylPatternUniModel.cpp
 *
 *  Created on: Sep 16, 2014
 *      Author: aviezerl
 */

//BASE_CC_FILE

#include "MethylPatternUniModel.h"
#include "ModelData.h"
#include "MethylPatternData.h"
//#include <stdlib.h>
//#include <cstdlib>

MethylPatternUniModel::MethylPatternUniModel()  {
    // TODO Auto-generated constructor stub

}

MethylPatternUniModel::MethylPatternUniModel(const MethylPatternUniModel& uni) :
        m_consensus_pattern(uni.m_consensus_pattern) {
}

MethylPatternUniModel::~MethylPatternUniModel() {
    // TODO Auto-generated destructor stub
    m_consensus_pattern.clear();
}

float MethylPatternUniModel::loglikelihood(const ModelData& data) {
    m_loglikelihood=0;
    for (int i=0; i<data.get_number_of_samples(); i++) {
        m_loglikelihood += get_log_prob(data.get_sample(i));
    }
    return(m_loglikelihood);
}

MethylPatternUniModel& MethylPatternUniModel::operator=(const MethylPatternUniModel& src) {
    if (this == &src) {
        return *this;
    } else {
        m_consensus_pattern = src.m_consensus_pattern;
    }
    return(*this);
}

int MethylPatternUniModel::learn(const ModelData& data) {
    vector<float> factors(data.get_number_of_samples(), 1);
    return(learn(data, factors));
}



void MethylPatternUniModel::simulate(const ModelData& original, ModelData* data) {
    for (int i=0; i<original.get_number_of_samples(); i++) {
        const Sample* os = original.get_sample(i);
        simulate_sample(os, data);
    }
}





void MethylPatternUniModel::print(ostream& out) const {
    vector<int>::const_iterator cur_a = m_consensus_pattern.begin();
    vector<int>::const_iterator end_a = m_consensus_pattern.end();
    out << "Uni-Model :: ";
    out << "center pattern: ";
    while(cur_a != end_a) {
        if ((*cur_a) >= 0) {
            out << (*cur_a);
        } else {
            out << "*";
        }
        cur_a++;
    }
    print_parameters(out);
    out << "\tllikelihood: " << m_loglikelihood << endl;
}

void MethylPatternUniModel::print_summary(ostream& out) const {
    vector<int>::const_iterator cur_a = m_consensus_pattern.begin();
    vector<int>::const_iterator end_a = m_consensus_pattern.end();
    while(cur_a != end_a) {
        if ((*cur_a) >= 0) {
            out << (*cur_a);
        } else {
            out << "*";
        }
        cur_a++;
    }

    print_parameters_summary(out);
}

void MethylPatternUniModel::init_params(const vector<int>& pattern, const ModelData& data) {
    m_consensus_pattern = pattern;
}



