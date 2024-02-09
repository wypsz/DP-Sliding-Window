//
// Created by 12293 on 2023/10/12.
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

using namespace  std;
#ifndef DP_SLIDING_WINDOW_PRIVATEHEAVYHITTERMG_H
#define DP_SLIDING_WINDOW_PRIVATEHEAVYHITTERMG_H


class PrivateMisraGries {
public:
    PrivateMisraGries(double epsilon, double lambda);
    void processStream(const std::vector<string>& stream);
    std::unordered_map<string, double> getResult();
    int generateSymmetricGeometricRandomNumber(double alpha);

private:
    double epsilon;
    double lambda;
    int threshold;
    std::unordered_map<string, double> f;
    unsigned int seed;
    std::mt19937 gen;

};



class PCC {
public:
    PCC(double epsilon, double lambda,double true_value, int W,const std::vector<std::string>& initialData);
    void processItem(string stream);
    std::unordered_map<std::string, double> query_all();
    double query(string& element);
    void deleteExpiredBlocks();

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
        unique_ptr<PrivateMisraGries> pmg = nullptr;
        vector<string> buffer;
    };


    double epsilon;
    double lambda;
    int W; // total time steps
    int W0; // block size at level 0
    int l; // max level
    vector<vector<Block>> blocks;
    int currentTime;
    void initBlocks(const std::vector<std::string>& initialData);
};



#endif //DP_SLIDING_WINDOW_PRIVATEHEAVYHITTERMG_H

