using namespace std;
/*
 * MethylPatternData.cpp
 *
 *  Created on: Sep 16, 2014
 *      Author: aviezerl
 */

//BASE_CC_FILE

#include "MethylPatternData.h"
#include "SamplingUtil.h"
using namespace Rcpp;

MethylPatternSample::MethylPatternSample(const MethylPatternSample& samp, const int& dist) {
    pattern  = samp.pattern;
    vector<int> inds(pattern.size(), -1);
    for (size_t i=0; i < inds.size(); i++) {
        inds[i] = i;
    }

    vector<int> flip_inds = random_sampling_without_replacement(inds, dist);
    for (size_t i = 0; i < flip_inds.size(); i++) {
        if (samp.pattern[i] == 0) {
            pattern[i] = 1;
        } else if (samp.pattern[i] == 1) {
            pattern[i] = 0;
        }
    }
}

MethylPatternData::MethylPatternData():m_id(0),m_na(0),m_pattern_length(0) {}

MethylPatternData::MethylPatternData(const DataFrame& df):m_na(0),m_pattern_length(0){
    StringVector pats = df["pat"];
    IntegerVector ids = df["fid"];
    m_id = ids(0);
    for(int i=0; i < pats.size(); i++){
        add_sample(Rcpp::as< std::string >(pats(i)));
    }
}

MethylPatternData::~MethylPatternData() {}

void MethylPatternData::add_sample(const string& pattern) {
    int non_na=0;
    if (m_pattern_length == 0) {
        m_pattern_length = pattern.length();
    }
    MethylPatternSample* s = new MethylPatternSample(m_pattern_length,-1);
    for (int i=0; i<m_pattern_length; i++) {
        if (pattern[i]=='0') {
            s->pattern[i] = 0;
            non_na++;
        } else {
            if (pattern[i] == '1') {
                s->pattern[i] = 1;
                non_na++;
            }
        }
    }

    m_na += m_pattern_length - non_na;
    ModelData::add_sample(s);
}


void MethylPatternData::print_samples(ostream& out) {
    out << m_id << "\t";
    ModelData::print_samples(out);
}

void MethylPatternData::clear() {
    //	cerr << "Methyl Pattern Data :: clear" << endl;
    ModelData::clear();
    m_na = 0;
    m_pattern_length=0;
}

vector<Sample*> MethylPatternData::unique_samples() const {
    MethylPatternSamplePtrSet uniq_samples_set;
    vector<Sample*> uniq_samples;
    for (int i=0; i<m_N; i++) {
        if ( uniq_samples_set.insert( ((MethylPatternSample*)(m_samples[i])) ).second  ) {
            uniq_samples.push_back(m_samples[i]);
        }
    }
    return(uniq_samples);
}

vector<Sample*> MethylPatternData::unique_samples(vector<int>& freqs) const {
    MethylPatternSamplePtrMap uniq_samples_map;
    vector<Sample*> uniq_samples;
    for (int i=0; i<m_N; i++) {
        uniq_samples_map[ (MethylPatternSample*)(m_samples[i]) ] ++;
    }
    for (const pair<MethylPatternSample*, int>& elem : uniq_samples_map) {
        uniq_samples.push_back(elem.first);
        freqs.push_back(elem.second);
    }
    return(uniq_samples);
}

float MethylPatternData::avg_m() {
    float meth = 0;
    float unmeth = 0;
    for (Sample* sm : m_samples) {
    	MethylPatternSample* s = (MethylPatternSample*)sm;
        for (int i=0; i<m_pattern_length; i++) {
            if (s->pattern[i] == 0) {
                unmeth++;
            } else if (s->pattern[i] == 1) {
            	meth++;
            }
        }
    }
    if ((meth + unmeth) == 0){
    	return(-1);
    }
    return(meth / (meth + unmeth));
}



