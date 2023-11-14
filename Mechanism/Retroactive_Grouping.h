//
// Created by 12293 on 2023/9/7.
//

#ifndef DP_SLIDING_WINDOW_RETROACTIVE_GROUPING_H
#define DP_SLIDING_WINDOW_RETROACTIVE_GROUPING_H
#include "string.h"
#include <vector>
#include <string>
#include <iostream>
#include <random>
#include <cmath>
#include <map>
#include <algorithm>
#include <set>

using namespace std;


class Retroactive_Grouping {
public:
    Retroactive_Grouping(const vector<string>& xw, vector<int> indices, double epsilon1, double epsilon2,vector<string>& dic, int w, int sub_num);
    void publication();
    void monitor(vector<string>& stream);
    vector<double> generateLaplacianRandom(double scale,int numSamples);
    bool bernoulliSample(double probability);
    double query(string& item);
    void insert(string& newItem);
private:
    vector<string> dic;
    int be;
    vector<int> index;
    vector<double> noise_histogram;// 当前时刻的histogram情况
    vector<int> histogram;
    vector<vector<int>> H;  // 存储所有时刻的counter情况
    vector<vector<double>> H_noise;
    vector<double> difference;
    vector<vector<double>> G;
    vector<double> Lamb; // 记录每个item上一次发布时加的噪声尺度Lap(lambda)
    int w;  // 窗口大小
    int m;
    double epsilon1;
    double epsilon2;
    vector<int> indices;
    int sub_num;

    int last_indices;
    vector<string> last_update;

};



#endif //DP_SLIDING_WINDOW_RETROACTIVE_GROUPING_H
