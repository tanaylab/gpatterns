/*
 * ModelDataMultiGroup.h
 *
 *  Created on: Oct 7, 2014
 *      Author: aviezerl
 */

#ifndef MODELDATAMULTIGROUP_H_
#define MODELDATAMULTIGROUP_H_

#include "ModelData.h"

class ModelDataMultiGroup {
public:
    ModelDataMultiGroup();
    virtual ~ModelDataMultiGroup();
    unsigned int get_number_of_groups() const { return(m_N); }
    virtual const ModelData* get_group(int group_id) const;
    virtual const vector<ModelData*>& get_groups() const { return(m_groups); }
    virtual void add_group(ModelData* d);
    virtual void print_groups(ostream&);
    virtual void remove_groups(vector<int>& inds, bool hard_remove = true);

    virtual void clear();

protected:
    unsigned int m_N;		//number of groups
    vector<ModelData*>	m_groups;
};

#endif /* MODELDATAMULTIGROUP_H_ */
