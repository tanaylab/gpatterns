using namespace std;
/*
 * MethylPatternMixtureModelGainLoss.cpp
 *
 *  Created on: Jul 23, 2015
 *      Author: aviezerl
 */

#include "MethylPatternMixtureModelGainLoss.h"
#include "MethylPatternUniModelGainLoss.h"
#include "UniformModel.h"


MethylPatternMixtureModelGainLoss::MethylPatternMixtureModelGainLoss() {
    // TODO Auto-generated constructor stub

}

MethylPatternMixtureModelGainLoss::~MethylPatternMixtureModelGainLoss() {
    // TODO Auto-generated destructor stub
}

MethylPatternMixtureModelGainLoss::MethylPatternMixtureModelGainLoss(const MethylPatternMixtureModelGainLoss& m): MethylPatternMixtureModel(m) {
    for (unsigned int i=0; i<m.m_models.size()-1; i++) {
        MethylPatternUniModelGainLoss* m_uni = (MethylPatternUniModelGainLoss*)m.m_models[i];
        MethylPatternUniModel* uni = new MethylPatternUniModelGainLoss(*m_uni);
        add_model(uni);
    }
    UniformModel* m_uniform = (UniformModel*)m.m_models[m.m_models.size()-1];
    UniformModel* uniform = new UniformModel(*m_uniform);
    add_model(uniform);
}

MethylPatternMixtureModelGainLoss::MethylPatternMixtureModelGainLoss(int number_of_mixtures, float uniform_prior): MethylPatternMixtureModel(number_of_mixtures, uniform_prior) {}

MethylPatternMixtureModelGainLoss::MethylPatternMixtureModelGainLoss(int number_of_mixtures, float uniform_prior, const int& K): MethylPatternMixtureModel(number_of_mixtures, uniform_prior, K) {
    for (int i=0; i<number_of_mixtures; i++) {
        MethylPatternUniModel* uni = new MethylPatternUniModelGainLoss();
        add_model(uni);
    }
    UniformModel* uniform = new UniformModel();
    add_model(uniform);
}

MethylPatternMixtureModelGainLoss& MethylPatternMixtureModelGainLoss::operator=(const MethylPatternMixtureModelGainLoss& m) {
    if(&m == this) {
        return(*this);
    }
    MethylPatternMixtureModelGainLoss temp(m);
    swap(temp);
    return(*this);
}

MethylPatternMixtureModel* MethylPatternMixtureModelGainLoss::clone() {
    return(new MethylPatternMixtureModelGainLoss(*this));
}

void MethylPatternMixtureModelGainLoss::copy(MethylPatternMixtureModel* m) {
    (*this) = *((MethylPatternMixtureModelGainLoss*)m);
}

void MethylPatternMixtureModelGainLoss::init_models_hamming(const ModelData& data) {
    //how to initialize the data....start by initializing one model with consensus, and
    //then selecting the pattern from the data that is farthest away (hamming distance)
    //from the existing and so on...
    //mixture coefficient is then computed by determining for each sample if it is closer to one or the other
    //and gain and loss is learned according to the center
    m_models[0]->learn(data);
    vector<int> consensus = ((MethylPatternUniModel*)(m_models[0]))->get_consensus_pattern();

    //initializing the uniform model (always the last one in the vector)
    UniformModel* uniform = (UniformModel*)m_models[m_models.size()-1];
    uniform->set_size(consensus.size());

    //initializing all other methylpattern unimodels
//    float gain = ((MethylPatternUniModelGainLoss*)(m_models[0]))->get_gain();
//    float loss = ((MethylPatternUniModelGainLoss*)(m_models[0]))->get_loss();
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
        ((MethylPatternUniModelGainLoss*)(m_models[c]))->init_params(next_consensus, data);

        adjust_params(c);

    }

    init_models(data);
}

void MethylPatternMixtureModelGainLoss::adjust_params(const int& c) {
    // if the gain or loss are half we may never change these parameters
    if (0.5 == ((MethylPatternUniModelGainLoss*)(m_models[c]))->get_gain()) {
        ((MethylPatternUniModelGainLoss*)(m_models[c]))->set_gain(0.49);
    }
    if (0.5 == ((MethylPatternUniModelGainLoss*)(m_models[c]))->get_loss()) {
        ((MethylPatternUniModelGainLoss*)(m_models[c]))->set_loss(0.49);
    }
    if (0 == ((MethylPatternUniModelGainLoss*)(m_models[c]))->get_gain()) {
        ((MethylPatternUniModelGainLoss*)(m_models[c]))->set_gain(0.001);
    }
    if (0 == ((MethylPatternUniModelGainLoss*)(m_models[c]))->get_loss()) {
        ((MethylPatternUniModelGainLoss*)(m_models[c]))->set_loss(0.001);
    }

}



