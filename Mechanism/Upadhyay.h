//
// Created by 12293 on 2023/12/4.
//

#ifndef DP_SLIDING_WINDOW_UPADHYAY_H
#define DP_SLIDING_WINDOW_UPADHYAY_H
#include <math.h>
#include <vector>
#include <string>
#include <functional>
#include "md5.h"
#include "new_alg.h"
#include <unordered_map>
#include <random>

using namespace std;


class Private_Heavy{
public:
    Private_Heavy(const vector<string>& xw, int w, int sub_num,int step,double epsilon, double gamma, double beta,double zeta,vector<string> dic,int flag);
    vector<CountMinSketch> ProcessSubWindow(const vector<string>& newItems);
    void ProcessNew(string item);
    // vector<pair<int, int>> windowAggregate();
    //pair<vector<double>,unordered_map<string,double>> Query_no();
    pair<vector<double>,unordered_map<string,double>> Query_all();
    vector<vector<double>> show_parameter();
    vector<double> Query();
    vector<double> generateLaplacianRandom(std::default_random_engine& generator, double scale, int numSamples);
private:
    int w;
    int step;
    int flag;
    double epsilon;
    double gamma;
    double beta;
    double zeta;
    int sub_size;
    int sub_num;
    int currentTime;
    vector<string> dic;
    vector<int> indices;
    vector<int> checkpoints;
    vector<vector<CountMinSketch>> Window_CMs;
    int last_win=0;
    vector<CountMinSketch> last_win_cms;
    std::default_random_engine generator;
    unsigned int initialSeed;
    unsigned int currentSeed;
};

#endif //DP_SLIDING_WINDOW_UPADHYAY_H
