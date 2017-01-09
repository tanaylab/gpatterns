using namespace std;
/*
 * ModelDataMultiGroup.cpp
 *
 *  Created on: Oct 7, 2014
 *      Author: aviezerl
 */

//BASE_CC_FILE

#include "ModelDataMultiGroup.h"

ModelDataMultiGroup::ModelDataMultiGroup():m_N(0) {

}

ModelDataMultiGroup::~ModelDataMultiGroup() {
//	clear();
}

const ModelData* ModelDataMultiGroup::get_group(int group_id) const {
    if (group_id < (int)m_N) {
        return m_groups[group_id];
    }
    return NULL;
}

void ModelDataMultiGroup::add_group(ModelData* d) {
    m_groups.push_back(d);
    m_N++;
}

void ModelDataMultiGroup::remove_groups(vector<int>& inds, bool hard_remove){
	if (0 == inds.size()){
		return;
	}
	vector<ModelData*> temp_groups;
	int j=0;
	for (unsigned int i=0; i < m_groups.size(); i++){
		if (i == (unsigned int)inds[j]){
			if (hard_remove){
				delete(m_groups[i]);
			}
			j++;
		} else {
			temp_groups.push_back(m_groups[i]);
		}
	}
	m_groups.swap(temp_groups);
	m_N = m_groups.size();
}

void ModelDataMultiGroup::print_groups(ostream& out) {
    out << m_N << endl;

    for (unsigned int i=0; i<m_N; i++) {
        m_groups[i]->print_samples(out);
        out << endl;
    }
}

void ModelDataMultiGroup::clear() {
    //cerr << "Model Data :: clear" << endl;
    for (unsigned int i=0; i<m_N; i++) {
        delete(m_groups[i]);
    }
    m_groups.clear();
    m_N=0;
}
