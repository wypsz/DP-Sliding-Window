//
// Created by 12293 on 2023/12/22.
//

#include "DPSW.h"
#include "math.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <sstream>
#include <algorithm>
using namespace std;



DPSW::DPSW(const vector<string>& xw, int w,int step, int sub_num,double rho, double gamma,
                                 double beta,double q, double alpha,unsigned int seedGenerator){



    this->w = w;
    this->step = step;
    this->rho = rho;
    this->gamma = gamma;
    this->beta = beta;
    this->q = q;// 0.8
    this->alpha = alpha;
    this->sub_num = sub_num;
    this->rho1 = rho*q;
    this->seedGenerator=seedGenerator;
    // 划分子窗口
    this->sub_size = w / sub_num;//子窗口大小
//    vector<double> checkpoints_noise;
//    vector<double> checkpoints;
    // 获得子窗口大小后直接调用sh获得checkpoint分割的list  s;
    cout << "hi" << endl;
    SmoothHistogram sh(this->alpha, sub_size, this->step);
    cout << "hi" << endl;
    sh.Add();
    indices = sh.getIndices();
    checkpoints = sh.getCheckpoints();
    for (int i = 0; i < indices.size(); i++) {
        cout << indices[i] << "     ";
    }
    cout << endl;
    cout << "indices finish" << endl;
    for (int i = 0; i < checkpoints.size(); ++i) {
        cout << checkpoints[i] << "     ";
    }
    cout << endl;
    cout << "checkpoints finish" << endl;
    cout << "size:" << indices.size() << "==" << checkpoints.size() << endl;

    for (int i = 0; i < sub_num; i++) {
        vector<CountMinSketch> cms = ProcessSubWindow(
                vector<string>(xw.begin() + (i * sub_size), xw.begin() + (i * sub_size + sub_size)));
        Window_CMs.push_back(cms);
        cout << i << "th sub_window finish!" << endl;
    }
    last_finish_index = -1;

    cout<<"seed Generator = "<<seedGenerator<<endl;

}

/*
 * indices: 1,7753,8884
 * checkpoints: 10000,2248,1117
 */


vector<CountMinSketch> DPSW::ProcessSubWindow(const vector<std::string> &newItems) {


    vector<CountMinSketch> cms;
    int index=0;
    cout<<"len indices:"<<indices.size()<<endl;
    for (int i=0;i<newItems.size();i++){
        if (i == 0){
            CountMinSketch cm (gamma,beta,rho1,seedGenerator++);
            //rest_budget = all_budget - rho;
            cm.update(newItems[i],1);
            cms.push_back(cm);
            index++;
        }else{
            // 更新已有的所有CM
            for (int j=0;j<index;j++){
                cms[j].update(newItems[i]);
            }
            if (i+1 == indices[index]){
                double budget=(rho1*pow((1-0.0001-q),index))/2.0;
                // 这里是checkpoint，新开sketch
                CountMinSketch cm(gamma,beta,budget,seedGenerator++);
                cm.update(newItems[i]);
                cms.push_back(cm);
                index++;
            }


        }
    }
    // cout<<cms.size()<<endl;
    return cms;
}

void DPSW::ProcessNew(std::string item) {
    if (last_win == 0){
        CountMinSketch cm(gamma,beta,rho1,seedGenerator++);
        cm.update(item);
        last_win ++;
        last_win_cms.push_back(cm);
        // 新建cms为滑出做准备
        for (int j=indices.size()-1;j>0;j--) {
            double budget=(rho1*pow((1-0.0001-q),j))/2.0;
            CountMinSketch cm_ready(gamma, beta, budget,seedGenerator++);// 正向维护的从小到大，每次访问第一个
            cm_ready.update(item);
            building_CM.push_back(cm_ready);
        }
        cout<<"building_CM size:"<<building_CM.size()<<endl;
        last_finish_index=-1;

    }else{
        for (int i=0;i<last_win_cms.size();i++){
            last_win_cms[i].update(item);
        }
        vector<int>::iterator it = find(indices.begin(),indices.end(),last_win+1);
        if ( it != indices.end()){
            int ind = distance(indices.begin(),it);
            double budget=(rho1*pow((1-0.0001-q),ind))/2.0;
            CountMinSketch cm(gamma,beta,budget,seedGenerator++);
            cm.update(item);
            last_win_cms.push_back(cm);
        }
        last_win++;

        for (int i=0;i<building_CM.size();i++){
            //cout<<"update_build"<<endl;
            building_CM[i].update(item);
        }
        if (last_finish_index == -1 && last_win == checkpoints[building_CM.size()]) {
            cout << "begin to use the first one" << endl;
            last_finish_index = 1;
        }else if (!building_CM.empty() && last_win == checkpoints[building_CM.size()-1]){
                cout<<"delete one"<<endl;
                building_CM.erase(building_CM.begin());
        }

        }

    if (last_win == sub_size){
        Window_CMs.push_back(last_win_cms);
        cout<<"Window_CMs"<<Window_CMs.size()<<endl;
        last_win = 0;
        last_win_cms.clear();
        Window_CMs.erase(Window_CMs.begin());

        last_finish_index = -1;
        if (!building_CM.empty()){
            //cout<<"problem"<<endl;
            for (int i=0;i<building_CM.size();i++){
                building_CM.erase(building_CM.begin());
            }
        }
        // cout<<"last win"<<endl;
    }
}


double DPSW::Query(string item) {
    double sum=0.0;
    if (last_win == 0 || indices.size()==1){
        for (int i=0;i<sub_num;i++){
            sum += Window_CMs[i][0].query(item);
        }
        return sum;
    }


    int tmp = last_win +1;
    auto it = std::lower_bound(indices.begin(), indices.end(), tmp);
    int pos = 0;

    if (it == indices.begin()) {
        pos = -1;
    } else {
        if (it == indices.end() || *it != tmp) {

            --it;
        }
        pos = std::distance(indices.begin(), it);
    }


    //double sum = 0.0;
    sum += Window_CMs[0][pos].query(item);
    for (int i=1;i<sub_num-1;i++){
        sum += Window_CMs[i][0].query(item);
    }
    if (last_finish_index != -1){
        // cout<<"using last"<<endl;
        sum += building_CM[0].query(item);
    }
    return sum;
}


vector<vector<double>> DPSW::show_parameter() {
    vector<CountMinSketch> cms = Window_CMs[0];
    vector<vector<double >> res;
    for (int i=0;i<cms.size();i++){
        vector<double> tmp = cms[i].get_parameter();
        res.push_back(tmp);
    }
    return res;

}