/*
 * MethylPatternData.h
 *
 *  Created on: Sep 16, 2014
 *      Author: aviezerl
 */

#ifndef METHYLPATTERNDATA_H_
#define METHYLPATTERNDATA_H_

#include "ModelData.h"
#include <boost/functional/hash.hpp>
#include <unordered_set>
#include <unordered_map>
#include <Rcpp.h>
using namespace Rcpp;


struct MethylPatternSample : public Sample {
    MethylPatternSample(int length, int val) : pattern(length,val) {}
    MethylPatternSample(const vector<int> vec) : pattern(vec) {}
    MethylPatternSample(const MethylPatternSample& samp, const int& dist);
    virtual ~MethylPatternSample() {
        pattern.clear();
    }

    virtual bool operator==(const MethylPatternSample& rhs) const{
      	return(this->pattern == rhs.pattern);
     }

    void print(ostream& out) const override{
        for (size_t j=0; j<pattern.size(); j++) {
            if (pattern[j] < 0 )
                out << "*";
            else
                out << pattern[j];
        }

    }
    int hamming_distance(const MethylPatternSample& p) const {
        return(hamming_distance(p.pattern));
    }

    int hamming_distance(const vector<int>& p) const {
        int dist=0;
        for (size_t j=0; j<pattern.size(); j++) {
            if (pattern[j] >= 0 && p[j] >= 0) {
                if (pattern[j] != p[j])
                    dist++;
            }
        }
        return dist;
    }

    MethylPatternSample* make_sample(MethylPatternSample* samp, const int& dist) const{
//    	assert((size_t)dist < pattern.size());
    	return(new MethylPatternSample(*this, dist));
    }

    vector<int> pattern;
};

struct PtrCmp{
	bool operator()(const MethylPatternSample* lhs, const MethylPatternSample* rhs) const  {
		return ((*lhs) == (*rhs));
	}
};

struct MethylPatternSamplePtrHash
{
    template <typename T>
    size_t operator()(const T* p) const
    {
        return boost::hash_range(p->pattern.begin(), p->pattern.end());
    }
};

typedef std::unordered_set<MethylPatternSample*, MethylPatternSamplePtrHash, PtrCmp> MethylPatternSamplePtrSet;
typedef std::unordered_map<MethylPatternSample*, int, MethylPatternSamplePtrHash, PtrCmp> MethylPatternSamplePtrMap;

class MethylPatternData: public ModelData {
public:
    MethylPatternData();
    MethylPatternData(const DataFrame& df);
    virtual ~MethylPatternData();

    void set_data_id(const int id) { m_id=id; }
    int get_data_id() const { return m_id; }
    int get_na() const { return m_na;  }
    int get_pattern_length() const { return m_pattern_length; }
    void add_sample(const string&);
    void print_samples(ostream&);
    void clear();
    float avg_m();
    vector<Sample*> unique_samples() const;
    vector<Sample*> unique_samples(vector<int>& sample_p) const;


protected:
    int		m_id;
    int		m_na;
    int		m_pattern_length;
};

#endif /* METHYLPATTERNDATA_H_ */
