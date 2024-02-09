//
// Created by 12293 on 2023/9/7.
//

#include "new_alg.h"
#include "math.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <sstream>
#include <algorithm>
using namespace std;

CountMinSketch::CountMinSketch(double gamma, double beta, double rho,unsigned int seed,unsigned int hseed)
        : t(1.0 / gamma), d(ceil(log(1.0 / beta))), sigma(0.0), E(0.0),gen(seed) {
    this->hseed = hseed;
    if (rho == 0.0){
        for (int i = 0; i < d; ++i) {
            vector<double> table(t, 0.0);
            C.push_back(table);
        }
    }else{
        sigma = sqrt(log(2.0 / beta) / rho);
        E1 = sqrt(2.0 * log(2.0 / beta) / rho) *
            sqrt(log(std::log(2.0 / beta) * 4.0 / gamma) / beta);
        E=0;
        budget = rho;
        for (int i = 0; i < d; ++i) {
            vector<double> noises(t);
            generateNoises(noises);
            C.push_back(noises);
        }
    }
}

void CountMinSketch::update(const string& x, double value) {
    auto hash_values = h(x);
    for (int i = 0; i < d; ++i) {
        int index = hash_values[i];
        C[i][index] += value;
    }
}

double CountMinSketch::query(const std::string x) {
    vector<int> hash_values = h(x);
    double min_value = C[0][hash_values[0]];
    for (int i = 1; i<d; ++i){
        min_value = min(C[i][hash_values[i]], min_value);
    }
    return min_value;
}

void CountMinSketch::generateNoises(std::vector<double>& noises) const {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> distribution(0, sigma);

    for (int i = 0; i < t; ++i) {
        noises[i] = distribution(gen);
    }
}



std::vector<int>CountMinSketch::h(const std::string& x) const {
    std::vector<int> hash_values;
    for (int i = 0; i < d; ++i) {
        MD5 md5(x + std::to_string(i+hseed));
        const byte* digest = md5.getDigest(); 

        std::uint64_t part1 = (static_cast<std::uint64_t>(digest[0]) << 56) |
                                (static_cast<std::uint64_t>(digest[1]) << 48) |
                                (static_cast<std::uint64_t>(digest[2]) << 40) |
                                (static_cast<std::uint64_t>(digest[3]) << 32);
        std::uint64_t part2 = (static_cast<std::uint64_t>(digest[4]) << 24) |
                                (static_cast<std::uint64_t>(digest[5]) << 16) |
                                (static_cast<std::uint64_t>(digest[6]) << 8) |
                                static_cast<std::uint64_t>(digest[7]);
        std::uint64_t combined = part1 | part2;

        int hash_value = combined % t;
        hash_values.push_back(hash_value);
    }
    return hash_values;
}


vector<double> CountMinSketch::get_parameter() {
    cout<<d<<" hash function"<<endl;
    cout<<t<<" bit"<<endl;
    cout<<"E:"<<E<<endl;
    cout<<"sigma:"<<sigma<<endl;
    cout<<"rho:"<<budget<<endl;
    cout<<"E1:"<<E1<<endl;
    vector<double> res;
    res.push_back(d);
    res.push_back(t);
    res.push_back(E);
    res.push_back(sigma);
    res.push_back(budget);
    res.push_back(E1);
    return res;
}



SmoothHistogram::SmoothHistogram(double alpha, int w,int step)
        : alpha(alpha), s(0), indices(), checkpoints(), N(0), w(w),step(step) {}


void SmoothHistogram::addItem() {
    for (int i = 0; i < s; ++i) {
        checkpoints[i]++;
    }

    s++;
    N++;
    indices.push_back(N);
    checkpoints.push_back(1);
    // cout<<"here"<<endl;


    for (int i = 0; i < s - 1; ++i) {
        std::vector<int> tmp_delete;
        for (int j = i + 1; j < s - 1; ++j) {
            double threshold = (1 - alpha) * checkpoints[i];
            if (checkpoints[j] >= threshold) {
                tmp_delete.push_back(j);
            } else {
                continue;
            }
        }
        for (int k = tmp_delete.size() - 1; k >= 0; --k) {
            s--;
            indices.erase(indices.begin() + tmp_delete[k]);
            checkpoints.erase(checkpoints.begin() + tmp_delete[k]);
        }
    }

}

void SmoothHistogram::Add() {
    for (int i=0;i<w;i++){
        addItem();
    }

}


vector<int> SmoothHistogram::getCheckpoints() {
    return checkpoints;
}

vector<int> SmoothHistogram::getIndices() {
    return indices;
}
