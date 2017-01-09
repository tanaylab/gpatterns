using namespace std;
/*
 * MethylPatternMixtureModelOneNoise.cpp
 *
 *  Created on: Jul 22, 2015
 *      Author: aviezerl
 */

#include "MethylPatternMixtureModelOneNoise.h"
#include "MethylPatternUniModelOneNoise.h"
#include "UniformModel.h"

MethylPatternMixtureModelOneNoise::MethylPatternMixtureModelOneNoise() {
    // TODO Auto-generated constructor stub

}

MethylPatternMixtureModelOneNoise::~MethylPatternMixtureModelOneNoise() {
    // TODO Auto-generated destructor stub
}

MethylPatternMixtureModelOneNoise::MethylPatternMixtureModelOneNoise(const MethylPatternMixtureModelOneNoise& m): MethylPatternMixtureModel(m) {
    for (unsigned int i=0; i<m.m_models.size()-1; i++) {
        MethylPatternUniModelOneNoise* m_uni = (MethylPatternUniModelOneNoise*)m.m_models[i];
        MethylPatternUniModel* uni = new MethylPatternUniModelOneNoise(*m_uni);
        add_model(uni);
    }
    UniformModel* m_uniform = (UniformModel*)m.m_models[m.m_models.size()-1];
    UniformModel* uniform = new UniformModel(*m_uniform);
    add_model(uniform);
}

MethylPatternMixtureModelOneNoise::MethylPatternMixtureModelOneNoise(int number_of_mixtures, float uniform_prior): MethylPatternMixtureModel(number_of_mixtures, uniform_prior) {}

MethylPatternMixtureModelOneNoise::MethylPatternMixtureModelOneNoise(int number_of_mixtures, float uniform_prior, const int& K): MethylPatternMixtureModel(number_of_mixtures, uniform_prior, K) {
    for (int i=0; i<number_of_mixtures; i++) {
        MethylPatternUniModel* uni = new MethylPatternUniModelOneNoise();
        add_model(uni);
    }
    UniformModel* uniform = new UniformModel();
    add_model(uniform);
}

MethylPatternMixtureModelOneNoise& MethylPatternMixtureModelOneNoise::operator=(const MethylPatternMixtureModelOneNoise& m) {
    if(&m == this) {
        return(*this);
    }
    MethylPatternMixtureModelOneNoise temp(m);
    swap(temp);
    return(*this);
}

MethylPatternMixtureModel* MethylPatternMixtureModelOneNoise::clone() {
    return(new MethylPatternMixtureModelOneNoise(*this));
}

void MethylPatternMixtureModelOneNoise::copy(MethylPatternMixtureModel* m) {
    (*this) = *((MethylPatternMixtureModelOneNoise*)m);
}

void MethylPatternMixtureModelOneNoise::init_models_hamming(const ModelData& data) {
    //how to initialize the data....start by initializing one model with consensus, and
    //then selecting the pattern from the data that is farthest away (hamming distance)
    //from the existing and so on...
    //mixture coefficient is then computed by determining for each sample if it is closer to one or the other
    //noise is set for all models according to the first model.
    m_models[0]->learn(data);
    vector<int> consensus = ((MethylPatternUniModel*)(m_models[0]))->get_consensus_pattern();

    //initializing the uniform model (always the last one in the vector)
    UniformModel* uniform = (UniformModel*)m_models[m_models.size()-1];
    uniform->set_size(consensus.size());

    //initializing all other methylpattern unimodels
    float noise = ((MethylPatternUniModelOneNoise*)(m_models[0]))->get_noise();
    for (unsigned int c=1; c<m_models.size()-1; c++) {
        int distance = 0;
        //locating data sample that is the farthest away
        int selected_data = 0;
        for (int i=0; i<data.get_number_of_samples(); i++) {
            MethylPatternSample* s = (MethylPatternSample*) data.get_sample(i);
            int sample_distance=0;
            for (unsigned int prev_c=0; prev_c<c; prev_c++) {
                sample_distance +=
                    s->hamming_distance(((MethylPatternUniModel*)m_models[prev_c])->get_consensus_pattern());
            }
            if (sample_distance > distance) {
                selected_data = i;
                distance = sample_distance;
            }
        }
        //fix missing data
        //at this point we have a sample from the data for which will initialize the next model
        vector<int> next_consensus = ((MethylPatternSample*)(data.get_sample(selected_data)))->pattern;
        for (unsigned int a=0; a<consensus.size(); a++) {
            if (consensus[a] >= 0 && next_consensus[a] < 0) {
                next_consensus[a] = 1-consensus[a];
            }
        }
        ((MethylPatternUniModelOneNoise*)(m_models[c]))->init_params(next_consensus, noise);
        adjust_params(c);
    }

    init_models(data);
}

void MethylPatternMixtureModelOneNoise::adjust_params(const int& c) {
    // if the noise is half we may never change this parameter
    if (0.5 == ((MethylPatternUniModelOneNoise*)(m_models[c]))->get_noise()) {
        ((MethylPatternUniModelOneNoise*)(m_models[c]))->set_noise(0.49);
    }
}





