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
using namespace std;

class CountMinSketch {
public:
    CountMinSketch(double gamma, double beta, double rho = 0.0);
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
    vector<vector<double>> C;
    void generateNoises(std::vector<double>& noises) const;
    vector<int> h(const std::string& x) const;
};

class SmoothHistogram {
public:
    SmoothHistogram(double alpha, int w,int step);
    void addItem();
    vector<int> getCheckpoints();
    vector<int> getIndices();
    void Add();
    void cut();
private:
    double alpha;
    int s;
    vector<int> indices;
    vector<int> checkpoints;
    int N;
    int w;
    int step;
};

class Private_counter{
public:
    Private_counter(const vector<string>& xw, int w,int step, int sub_num,double rho, double gamma, double beta,double q,double alpha, double rest_budget);
    vector<CountMinSketch> ProcessSubWindow(const vector<string>& newItems);
    void ProcessNew(string item);
    // vector<pair<int, int>> windowAggregate();
    double Query(string item);
    vector<vector<double>> show_parameter();
private:
     int w; //窗口大小
     int step;
    double rho;
    double gamma;
    double beta;
    double q;
    double alpha;
    int sub_size;
    int sub_num;
    double rest_budget;
    double all_budget;
    // vector<string> dic;
    vector<int> indices;
    vector<int> checkpoints;
    vector<vector<CountMinSketch>> Window_CMs; // 所有的子窗口中CM汇总
    int last_win=0;//统计最后一个子窗口内有多少item
    vector<CountMinSketch> last_win_cms;
};

#endif //DP_SLIDING_WINDOW_NEW_ALG_H
