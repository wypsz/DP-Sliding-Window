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

CountMinSketch::CountMinSketch(double gamma, double beta, double rho)
        : t(1.0 / gamma), d(ceil(log(1.0 / beta))), sigma(0.0), E(0.0) {
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

vector<int> CountMinSketch::h(const std::string& x) const {
    std::vector<int> hash_values;
    for (int i = 0; i < d; ++i) {
        auto val = MD5(x).toStr() + to_string(i);
        std::hash<std::string> hasher;
        auto hash_value = hasher(val) % t;
        hash_values.push_back(hash_value);

    }

    return hash_values;

}

//vector<int> CountMinSketch::query_vector(vector<std::string> dic) {
//    vector<int> freq_s;
//    for (int i=0;i<dic.size();i++){
//        int freq = query(dic[i]);
//        freq_s.push_back(freq);
//    }
//    return freq_s;
//}

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
        // L1-norm是(alpha,alpha)-smooth函数
        // alpha=zeta/2
        // s = O(1/zeta log w log(1/beta)) checkpoint数量
        // indices : 维护s个检查点的indices  X1<X2<...<XS=n X1=1
        // checkpoints 维护相应的L1 norm  对于第一个元素，传入一个frequency vector，L1-norm相当于统计有多少个元素
        // N=1 stream的长度

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
    cut();
    if (checkpoints[checkpoints.size()-1] < step){
        checkpoints.pop_back();
        indices.pop_back();
    }
}

void SmoothHistogram::cut() {
    int size = indices.size();
    for (int i=0;i<size-1;i++){
        for(int j=i+1;j<size;j++){
            while (indices[j]-indices[i]<step){
                indices.erase(indices.begin()+j);
                checkpoints.erase(checkpoints.begin()+j);
                size--;
                if (j >= size){
                    break;
                }
            }
            break;
        }
    }
}

vector<int> SmoothHistogram::getCheckpoints() {
    return checkpoints;
}

vector<int> SmoothHistogram::getIndices() {
    return indices;
}


Private_counter::Private_counter(const vector<string>& xw, int w,int step, int sub_num,double rho, double gamma,
                                 double beta,double q, double alpha, double all_budget) {
    /*
     * xw:初始化时需要有一个窗口的数据灌入，构造第一个窗口相应的结构根号w个子窗口(批量更新，传入的是vector，统计好数字的)
     * w:窗口大小（可开根号）
     * rho:隐私预算
     * gamma:count-min sketch 参数
     * zeta:count-min sketch 参数
     * beta:count-min sketch 参数
     * dic:vector<string>结构，按顺序存放item，便于后续将一个item转为一个vector
     */

    this->w = w;
    this->step = step;
    this->rho = rho;
    this->gamma = gamma;
    this->beta = beta;
    this->q = q;
    this->alpha = alpha;
    this->all_budget = all_budget;
    this->sub_num = sub_num;
    // 划分子窗口
    this->sub_size = w/sub_num;//子窗口大小
//    vector<double> checkpoints_noise;
//    vector<double> checkpoints;
    // 获得子窗口大小后直接调用sh获得checkpoint分割的list  s;
    cout<<"hi"<<endl;
    SmoothHistogram sh(alpha,sub_size,step);
    cout<<"hi"<<endl;
    sh.Add();
    indices = sh.getIndices();
    checkpoints = sh.getCheckpoints();
    for (int i=0;i<indices.size();i++){
        cout<<indices[i]<<"     ";
    }
    cout<<endl;
    cout<<"indices finish"<<endl;
    for (int i = 0; i < checkpoints.size(); ++i) {
        cout<<checkpoints[i]<<"     ";
    }
    cout<<endl;
    cout<<"checkpoints finish"<<endl;
    cout<<"size:"<<indices.size()<<"=="<<checkpoints.size()<<endl;


    vector<vector<string>> sub_windows;
    for (int i = 0; i < sub_num; i++) {
        vector<string> sub_win;
        for (int j=0;j<sub_size;j++){
            sub_win.push_back(xw[i*sub_size+j]);
        }
        sub_windows.push_back(sub_win);
    }
    // 处理每个子窗口
//    for (int i=0;i<sub_num;i=i+sub_size){
//        vector<CountMinSketch> cms = ProcessSubWindow((xw.begin()+i,xw.begin()+i+sub_size),gamma,beta,rho,q);
//    }
    for (int i=0;i<sub_num;i++){
        vector<CountMinSketch> cms = ProcessSubWindow(sub_windows[i]);
        Window_CMs.push_back(cms);
        vector<string>().swap(sub_windows[i]);
        cout<<i<<"th sub_window finish!"<<endl;
//        if (i==0){
//            for (int j=0;j<cms.size();j++){
//                cms[j].get_parameter();
//            }
//        }
    }
}

vector<vector<double>> Private_counter::show_parameter() {
    vector<CountMinSketch> cms = Window_CMs[0];
    vector<vector<double >> res;
    for (int i=0;i<cms.size();i++){
        vector<double> tmp = cms[i].get_parameter();
        res.push_back(tmp);
    }
    return res;

}

vector<CountMinSketch> Private_counter::ProcessSubWindow(const vector<std::string> &newItems) {
    // 对于整个子窗口的元素进行处理
    vector<CountMinSketch> cms;
    int index=0;//定位是第几个checkpoint
    for (int i=0;i<newItems.size();i++){
        // 第一个元素，新建一个cm并插入一个元素
        if (i == 0){
            CountMinSketch cm (gamma,beta,rho);
            rest_budget = all_budget - rho;
            cm.update(newItems[i],1);
            cms.push_back(cm);
            index++;
        }else{
            // 更新已有的所有CM
            for (int j=0;j<index;j++){
                cms[j].update(newItems[i]);
            }
            if (i+1 == indices[index]){
                double budget;
                // 这里是checkpoint，新开sketch
                if (index == indices.size()-1){
                    // 最后一个checkpoint，分配所有的预算给它
                    budget = rest_budget;
                }else{
                    budget = rho*pow(q,index);
                    rest_budget -= budget;
                }
                CountMinSketch cm(gamma,beta,budget);
                cm.update(newItems[i]);
                cms.push_back(cm);
                index++;
            }


        }
    }
     // cout<<cms.size()<<endl;
    return cms;
}


void Private_counter::ProcessNew(std::string item) {
    // 判断是不是第一个元素，如果是需要新建cm
    if (last_win == 0){
        CountMinSketch cm(gamma,beta,rho);
        rest_budget = all_budget - rho;
        cm.update(item);
        last_win ++;
        last_win_cms.push_back(cm);
    }else{
        for (int i=0;i<last_win_cms.size();i++){
            last_win_cms[i].update(item);
        }
        vector<int>::iterator it = find(indices.begin(),indices.end(),last_win+1);
        if ( it != indices.end()){
            //判断这个点是否需要是checkpoint，如果是也需要新建cm
            // 找到这个点是第几个checkpoint，便于计算所分配的预算rhoi
            int ind = distance(indices.begin(),it);
            double budget;
            if (ind == indices.size()-1){
                // 如果是最后一个checkpoint，分配给它所有的预算
                budget = rest_budget;
            }else{
                budget = rho * pow(q,ind);
                rest_budget -= budget;
            }

            CountMinSketch cm(gamma,beta,budget);
            cm.update(item);
            last_win_cms.push_back(cm);
        }
        last_win++;
    }

    // 如若这个元素到达好正好一个窗口满了，设置一下操作。。。
    if (last_win == sub_size){
        // 1. 将这个子窗口的所有cm更新到全局中
        Window_CMs.push_back(last_win_cms);
        cout<<"Window_CMs"<<Window_CMs.size()<<endl;
        // 删除
        last_win = 0;
        last_win_cms.clear();
        Window_CMs.erase(Window_CMs.begin());
        // cout<<"last win"<<endl;
    }
}


double Private_counter::Query(string item) {
    // 1. 找到当前的窗口需要使用那些checkpoints
    double sum=0.0;
    if (last_win == 0 || indices.size()==1){
        // 直接把第一个checkpoints都加起来
        for (int i=0;i<sub_num;i++){
            sum += Window_CMs[i][0].query(item);
        }
        return sum;
    }
    int tmp = last_win+1;
    // int tmp = 12;
    //去第一个窗口中找最相近的点
    auto it = std::lower_bound(indices.begin(), indices.end(), tmp);
    int pos = 0; //最接近的checkpoint的位置（第pos个checkpoint）

    if (*it == tmp) {
        //std::cout << "Found x in arr\n";
        pos = distance(indices.begin(),it);
    } else {
        int idx = it - indices.begin();
        int diff1 = tmp - indices[idx - 1];
        int diff2 = indices[idx] - tmp;
        if (diff1 < diff2) {
            // std::cout << "The element closest to x is " << checkpoints[idx - 1] << std::endl;
            pos = idx-1;
        } else {
            // std::cout << "The element closest to x is " << checkpoints[idx] << std::endl;
            pos = idx;
        }
    }

    // 2. 开始组装每个窗口的结果
    //double sum = 0.0;
    sum += Window_CMs[0][pos].query(item);
    for (int i=1;i<sub_num;i++){
        sum += Window_CMs[i][0].query(item);
    }
    return sum;
}

