/*
 * SamplingUtil.h
 *
 *  Created on: Nov 12, 2014
 *      Author: aviezerl
 */

#ifndef SAMPLINGUTIL_H_
#define SAMPLINGUTIL_H_


#include <vector>
#include <random>
#include <algorithm>

inline std::vector<int> random_sampling_without_replacement(std::vector<int>& vData, int N) {
    // Random generator
    std::random_device rd;
    std::mt19937 randomGenerator(rd());

    int max = static_cast<int>(vData.size()-1);
    std::vector<int> vResult;

    for (int n=1; n<=N; n++) {
        std::uniform_int_distribution<> uniformDistribution(0,max);
        int index = uniformDistribution(randomGenerator);
        std::swap(vData[index],vData[max]);
        vResult.push_back(vData[max]);
        max--;
    }

    return vResult;
}

template <typename T>
std::vector<size_t> ordered(std::vector<T> const& values, bool decreasing=false) {
    std::vector<size_t> indices(values.size());
    std::iota(begin(indices), end(indices), static_cast<size_t>(0));

    if (decreasing) {
        std::sort(
            begin(indices), end(indices),
            [&](size_t a, size_t b) {
                return values[a] > values[b];
            }
        );
    } else {
        std::sort(
            begin(indices), end(indices),
            [&](size_t a, size_t b) {
                return values[a] < values[b];
            }
        );
    }
    return indices;
}

inline vector<float> p_adjust(const vector<float>& p) {
    int n = p.size();
    vector<float> q(n, -1);
    vector<int> inds(n);

    std::iota(inds.rbegin(), inds.rend(), static_cast<size_t>(0));
    vector<size_t> o = ordered<float>(p, true);
    vector<size_t> ro = ordered<size_t>(o);

    float min = 1.0;
    for (int i=0; i < n; i++) {
        float f = ((float)n / (float)(inds[i]+1)) * p[o[i]];
        if (f < min) {
            q[i] = f;
            min = f;
        } else {
            q[i] = min;
        }
    }

    vector<float> qvals(n);
    for (unsigned int i=0; i < ro.size(); i++){
    	qvals[i] = q[ro[i]];
    }
    return(qvals);
}



//inline double factorial(double nValue) {
//    double result = nValue;
//    double result_next;
//    double pc = nValue;
//    do {
//        result_next = result*(pc-1);
//        result = result_next;
//        pc--;
//    } while(pc>2);
//    return result;
//}

inline int nchoosek(const int& n, const int& g_k)
{
	int k = g_k;
    int res = 1;

    // Since C(n, k) = C(n, n-k)
    if ( k > n - k )
        k = n - k;

    // Calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
    for (int i = 0; i < k; ++i)
    {
        res *= (n - i);
        res /= (i + 1);
    }

    return res;
}
//inline int nchoosek(const int& n, const int& k)   {
//    if(k == 1)
//        return n;
//    int result = (factorial(n))/(factorial(k)*factorial((n - k)));
//    assert(result <= 1000000);
//    return (factorial(n))/(factorial(k)*factorial((n - k)));
//}


#endif /* SAMPLINGUTIL_H_ */
