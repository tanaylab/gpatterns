/*
 * MethylPatternMixtureModelsCompare.h
 *
 *  Created on: Aug 6, 2015
 *      Author: aviezerl
 */

#ifndef METHYLPATTERNMIXTUREMODELSCOMPARE_H_
#define METHYLPATTERNMIXTUREMODELSCOMPARE_H_

#include "MethylPatternMixtureModelScanner.h"

class MethylPatternMixtureModelsCompare {
public:
    MethylPatternMixtureModelsCompare();
    virtual ~MethylPatternMixtureModelsCompare();

    void compare_models(const MethylPatternData& data, MethylPatternMixtureModelScanner* model1, MethylPatternMixtureModelScanner* model2, const int& K, const int& min_k, const int& num_of_sims, ostream& model1_out, ostream& model2_out);
};


#endif /* METHYLPATTERNMIXTUREMODELSCOMPARE_H_ */
