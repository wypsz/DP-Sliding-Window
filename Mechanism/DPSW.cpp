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

    // 对于整个子窗口的元素进行处理
    vector<CountMinSketch> cms;
    int index=0;//定位是第几个checkpoint
    cout<<"len indices:"<<indices.size()<<endl;
    for (int i=0;i<newItems.size();i++){
        // 第一个元素，新建一个cm并插入一个元素
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
                //cout<<"i+1,"<<i+1<<"indices[index]"<<index<<","<<indices[index]<<endl;
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
    // 查看last_finish_CM 是否新到来元素后，cms1达到足够的数量可以更新last_finish_CM
    // 判断是不是第一个元素，如果是需要新建cm
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
            //判断这个点是否需要是checkpoint，如果是也需要新建cm
            // 找到这个点是第几个checkpoint，便于计算所分配的预算rhoi
            int ind = distance(indices.begin(),it);
            double budget=(rho1*pow((1-0.0001-q),ind))/2.0;
            CountMinSketch cm(gamma,beta,budget,seedGenerator++);
            cm.update(item);
            last_win_cms.push_back(cm);
        }
        last_win++;
        // 更新building_cm中的每一个
        for (int i=0;i<building_CM.size();i++){
            //cout<<"update_build"<<endl;
            building_CM[i].update(item);
        }
        // 某个cm是否已经建满, building_CM有x个元素，说明现在
        if (last_finish_index == -1 && last_win == checkpoints[building_CM.size()]) {
            cout << "begin to use the first one" << endl;
            last_finish_index = 1;
        }else if (!building_CM.empty() && last_win == checkpoints[building_CM.size()-1]){
            // 有一个已经满了，将他放在finish中，可以支持查询，并将他从当前vector中删除，相应地更新last
                cout<<"delete one"<<endl;
                // 说明是第二个开始可以支持查询了，先删掉第一个，再更新index
                building_CM.erase(building_CM.begin());// 删掉第一个
        }

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
    // 1. 找到当前的窗口需要使用那些checkpoints
    double sum=0.0;
    if (last_win == 0 || indices.size()==1){
        // 直接把第一个checkpoints都加起来
        for (int i=0;i<sub_num;i++){
            sum += Window_CMs[i][0].query(item);
        }
        return sum;
    }


    // 查找方式1. 对于第一个窗口，找到最接近的，对于最后一个窗口直接使用当前保存的状态。
    
    /*
    int tmp = last_win+1;
    //去第一个窗口中找最相近的点
    auto it = std::lower_bound(indices.begin(), indices.end(), tmp);
    int pos = 0; //最接近的checkpoint的位置（第pos个checkpoint）

    if (*it == tmp) {
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
    */
    // 查找方式2，往前一个找
    int tmp = last_win +1;
    auto it = std::lower_bound(indices.begin(), indices.end(), tmp);
    int pos = 0;

    if (it == indices.begin()) {
        // tmp 小于所有元素，没有上一个位置
        // 根据具体需求处理这种情况
        pos = -1; // 或者任何表示 "没有有效位置" 的值
    } else {
        if (it == indices.end() || *it != tmp) {
            // 如果 tmp 大于所有元素，或者没有在数组中找到完全匹配的元素
            // it 指向的是第一个大于 tmp 的元素或者是 end()
            --it; // 移动到上一个元素
        }
        pos = std::distance(indices.begin(), it);
    }



    // 开始组装每个窗口的结果
    //double sum = 0.0;
    sum += Window_CMs[0][pos].query(item);
    for (int i=1;i<sub_num-1;i++){
        sum += Window_CMs[i][0].query(item);
    }
    // 最后一个窗口用last_finish
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