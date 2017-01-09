using namespace std;
/*
 * ModelData.cpp
 *
 *  Created on: Sep 16, 2014
 *      Author: aviezerl
 */

//BASE_CC_FILE

#include "ModelData.h"

ModelData::ModelData():m_N(0) {}

ModelData::ModelData(const vector<ModelData>& data):m_N(0) {
    for (vector<ModelData>::const_iterator it = data.begin(); it != data.end(); ++it) {
        for (int i = 0; i < (*it).get_number_of_samples(); i++) {
            add_sample((*it).m_samples[i]);
        }
    }
}

ModelData::ModelData(const vector<ModelData*>& data):m_N(0) {
    for (vector<ModelData*>::const_iterator it = data.begin(); it != data.end(); ++it) {
        for (int i = 0; i < (*it)->get_number_of_samples(); i++) {
            add_sample((*it)->m_samples[i]);
        }
    }
}

ModelData::~ModelData() {
	clear();
    // TODO Auto-generated destructor stub
}

const Sample* ModelData::get_sample(int sample_id) const {
    if (sample_id < m_N) {
        return m_samples[sample_id];
    }
    return NULL;
}

void ModelData::add_sample(Sample* s) {
    m_samples.push_back(s);
    m_N++;
}

void ModelData::print_samples(ostream& out) const {
    out << m_N << endl;

    for (int i=0; i<m_N; i++) {
        m_samples[i]->print(out);
        out << endl;
    }
}

void ModelData::clear() {
    //cerr << "Model Data :: clear" << endl;
    for (int i=0; i<m_N; i++) {
        delete(m_samples[i]);
    }
    m_samples.clear();
    m_N=0;
}
