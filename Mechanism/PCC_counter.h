//
// Created by 12293 on 2023/12/10.
//
#include <iostream>
#include <unordered_map>
#include <cmath>
#include <random>
#include <vector>
#include <math.h>
#include <algorithm>
#include <memory>
#include <string>
#include <chrono>

using namespace std;
#ifndef DP_SLIDING_WINDOW_PCC_COUNTER_H
#define DP_SLIDING_WINDOW_PCC_COUNTER_H

class MG{
public:
    MG(double alpha);
    void processStream(const std::vector<string>& stream);
    std::unordered_map<string, double> getResult();
    double alpha;
    int threshold;
private:
    //double alpha;
    //int T;
    std::unordered_map<string, double> f;
};


class PCC_counter {
public:
    PCC_counter(double epsilon,double alpha, int W,double true_the,const std::vector<std::string>& initialData);

    void processItem(string stream);
    std::unordered_map<std::string, double> query_all();
    //double query(string& element);
    void deleteExpiredBlocks();
    vector<double> generateLaplacianRandom(std::default_random_engine& generator, double scale, int numSamples);
    unordered_map<string,int> convertToFrequencyMap(vector<string>& vec);
private:
    enum BlockState {
        FUTURE,
        UNDER_CONSTRUCTION,
        ACTIVE,
        EXPIRED
    };

    struct Block {
        int startTime;
        int endTime;
        BlockState state;
        unique_ptr<MG> pmg = nullptr;
        vector<string> buffer;
        unordered_map<string,int> fre_true;
    };

    double true_the;
    double epsilon;
    double alpha;
    int W; // total time steps
    int S1; // block size at level 0
    int L; // max level
    vector<vector<Block>> blocks;
    int currentTime;
    void initBlocks(const std::vector<std::string>& initialData);
    std::default_random_engine generator;
    unsigned int initialSeed;
    unsigned int currentSeed;
};


#endif //DP_SLIDING_WINDOW_PCC_COUNTER_H
