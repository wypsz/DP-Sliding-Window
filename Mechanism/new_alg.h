//
// Created by 12293 on 2023/9/7.
//

#ifndef DP_SLIDING_WINDOW_NEW_ALG_H
#define DP_SLIDING_WINDOW_NEW_ALG_H
#include <math.h>
#include <vector>
#include <string>
#include <functional>
#include "md5.h"
#include <random>
using namespace std;

class CountMinSketch {
public:
    CountMinSketch(double gamma, double beta, double rho = 0.0,unsigned int seed=0,unsigned int hseed=0);
    void update(const std::string& x, double value = 1.0);
    double query(const std::string x);
    // vector<int> query_vector(vector<string> dic);
    vector<double> get_parameter();

private:
    // d行t列d个hash function 映射到t个bit
    int t;
    int d;
    double sigma;
    double E;
    double budget;
    unsigned int hseed;
    
    vector<vector<double>> C;
    void generateNoises(std::vector<double>& noises) const;
    vector<int> h(const std::string& x) const;
    double E1;
    std::mt19937 gen;
};

class SmoothHistogram {
public:
    SmoothHistogram(double alpha, int w,int step);
    void addItem();
    vector<int> getCheckpoints();
    vector<int> getIndices();
    void Add();
private:
    double alpha;
    int s;
    vector<int> indices;
    vector<int> checkpoints;
    int N;
    int w;
    int step;
};
#endif //DP_SLIDING_WINDOW_NEW_ALG_H
