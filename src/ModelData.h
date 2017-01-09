/*
 * ModelData.h
 *
 *  Created on: Sep 16, 2014
 *      Author: aviezerl
 */

#ifndef MODELDATA_H_
#define MODELDATA_H_

#include <vector>
#include <iostream>
using namespace std;

struct Sample {
    Sample() {}
    virtual ~Sample() {}
    virtual void print(ostream& out) const {}
}
;

class ModelData {
public:
    ModelData();
    ModelData(const vector<ModelData>&);
    ModelData(const vector<ModelData*>&);
    virtual ~ModelData();
    int get_number_of_samples() const { return(m_N); }
    virtual const Sample* get_sample(int sample_id) const;
    virtual void add_sample(Sample* s);
    virtual void print_samples(ostream&) const;
    virtual void clear();

    virtual vector<Sample*> unique_samples() const { return(m_samples);}
    virtual vector<Sample*> unique_samples(vector<int>&) const { return(m_samples);}

protected:
    int m_N;		//number of samples
    vector<Sample*>	m_samples;
};

#endif /* MODELDATA_H_ */
