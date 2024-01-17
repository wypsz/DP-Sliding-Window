//
// Created by 12293 on 2023/12/4.
//

#include "Upadhyay.h"
#include "new_alg.h"
#include "math.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <sstream>
#include <algorithm>
using namespace std;

Private_Heavy::Private_Heavy(const vector<string>& xw, int w, int sub_num,int step,double epsilon,
                             double gamma,double beta, double zeta,vector<string> dic,int flag){
    this->w = w;
    this->epsilon = epsilon;
    this->step = step;
    this->gamma = gamma;
    this->beta = beta;
    this->zeta = zeta;
    this->sub_num = sub_num;
    this->dic = dic;
    this->sub_size = w/sub_num;
    this->flag=flag;
    this->initialSeed = 0;
    this->currentSeed = initialSeed;
    cout<<"hi"<<endl;

    SmoothHistogram sh(beta,sub_size,step);
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
    if (flag == 1){
        vector<string> sub_win;
        for (int i=sub_num*sub_size;i<w;i++){
            sub_win.push_back(xw[i]);
        }
        sub_windows.push_back(sub_win);
        this->sub_num +=1;
    }
    for (int i=0;i<sub_windows.size();i++) {
        vector<CountMinSketch> cms = ProcessSubWindow(sub_windows[i]);
        Window_CMs.push_back(cms);
        vector<string>().swap(sub_windows[i]);
        cout << i << "th sub_window finish!" << endl;
    }

    //cout<<"seed Generator = "<<seedGenerator<<endl;
}

vector<CountMinSketch> Private_Heavy::ProcessSubWindow(const vector<std::string> &newItems) {
    vector<CountMinSketch> cms;
    int index=0;
    cout<<"len indices:"<<indices.size()<<endl;
    for (int i=0;i<newItems.size();i++){
        // 第一个元素，新建一个cm并插入一个元素
        if (i == 0){
            CountMinSketch cm (gamma/exp(1),zeta,0,0);
            cm.update(newItems[i],1);
            cms.push_back(cm);
            index++;
        }else{
            for (int j=0;j<index;j++){
                cms[j].update(newItems[i]);
            }
            if (i+1 == indices[index]){
                //cout<<"i+1,"<<i+1<<"indices[index]"<<index<<","<<indices[index]<<endl;
                CountMinSketch cm(gamma/exp(1),zeta,0);
                cm.update(newItems[i]);
                cms.push_back(cm);
                index++;
                //cout<<"2222"<<endl;
            }


        }
    }
    // cout<<cms.size()<<endl;
    return cms;
}

void Private_Heavy::ProcessNew(std::string item) {
    if (last_win == 0){
        CountMinSketch cm(gamma/exp(1),zeta,0,0);
        cm.update(item);
        last_win ++;
        last_win_cms.push_back(cm);
    }else{
        for (int i=0;i<last_win_cms.size();i++){
            last_win_cms[i].update(item);
        }
        vector<int>::iterator it = find(indices.begin(),indices.end(),last_win+1);
        if ( it != indices.end()){
            int ind = distance(indices.begin(),it);

            CountMinSketch cm(gamma/exp(1),zeta,0,0);
            cm.update(item);
            last_win_cms.push_back(cm);
        }
        last_win++;
    }

    if (last_win == sub_size){
        Window_CMs.push_back(last_win_cms);
        cout<<"Window_CMs"<<Window_CMs.size()<<endl;
        last_win = 0;
        last_win_cms.clear();
        Window_CMs.erase(Window_CMs.begin());
    }

}

vector<double> Private_Heavy::generateLaplacianRandom(std::default_random_engine& generator, double scale, int numSamples) {
    std::uniform_real_distribution<double> uniformDist(0.0, 1.0);
    std::vector<double> randomNumbers;
    randomNumbers.reserve(numSamples);

    for (int i = 0; i < numSamples; ++i) {
        double u = uniformDist(generator);
        if (u < 0.5) {
            randomNumbers.push_back(scale * log(2.0 * u));
        } else {
            randomNumbers.push_back(-scale * log(2.0 * (1.0 - u)));
        }
    }

    return randomNumbers;
}


// pair<vector<double>, unordered_map<string,double>> Private_Heavy::Query_no() {
//     //vector<double> fre_no_noise;
//     //unordered_map<string,double> item_fre_no_noise;
//     vector<double> fre_noise;
//     unordered_map<string,double> item_fre_noise;
//     double L1;
//     pair<vector<double>,unordered_map<string,double>> res;
//     vector<double> noise_vec = generateLaplacianRandom(1.0*sub_size / epsilon, dic.size());
//     for (int k=0;k<dic.size();k++){
//         // 获取每个元素的频率
//         double sum = 0.0;
//         if (flag ==0){

//             // 1. 找到当前的窗口需要使用那些checkpoints
//             //double sum = 0.0;
//             if (last_win == 0 || indices.size() == 1) {
//                 // 直接把第一个checkpoints都加起来
//                 for (int i = 0; i < sub_num; i++) {
//                     double r1 = Window_CMs[i][0].query(dic[k]);
//                     sum += r1;
//                     //cout<<sum<<"     ";
//                 }
//                 //cout<<sum;
//                 fre_noise.push_back(sum + noise_vec[k]);
//                 L1 = w;
//             } else {
//                 int tmp = last_win +1;
//                 auto it = std::lower_bound(indices.begin(), indices.end(), tmp);
//                 int pos = 0;

//                 if (it == indices.begin()) {
//                     // tmp 小于所有元素，没有上一个位置
//                     // 根据具体需求处理这种情况
//                     pos = -1; // 或者任何表示 "没有有效位置" 的值
//                 } else {
//                     if (it == indices.end() || *it != tmp) {
//                         // 如果 tmp 大于所有元素，或者没有在数组中找到完全匹配的元素
//                         // it 指向的是第一个大于 tmp 的元素或者是 end()
//                         --it; // 移动到上一个元素
//                     }
//                     pos = std::distance(indices.begin(), it);
//                 }

//                 L1 = checkpoints[pos] + last_win + sub_size * (sub_num - 2);

//                 // 2. 开始组装每个窗口的结果
//                 //double sum = 0.0;
//                 sum += Window_CMs[0][pos].query(dic[k]);
//                 for (int i = 1; i < sub_num; i++) {
//                     sum += Window_CMs[i][0].query(dic[k]);
//                 }
//                 fre_noise.push_back(sum + noise_vec[k]);

//             }
//         }else{
//             // 需要查找真实的窗口的大小，
//             int tmp1 =(w-last_win);
//             int n = 0 ;
//             while (tmp1>sub_size){
//                 tmp1 -= sub_size;
//                 n++;
//             }
//             int tmp = sub_size - tmp1;
//             auto it = std::lower_bound(indices.begin(), indices.end(), tmp);
//             int pos = 0; //最接近的checkpoint的位置（第pos个checkpoint）

//             if (*it == tmp) {
//                 //std::cout << "Found x in arr\n";
//                 pos = distance(indices.end(), it);
//             } else {
//                 int idx = it - indices.begin();
//                 if (idx != 0){
//                     pos = idx-1;
//                 }else{
//                     pos = idx;
//                 }
//             }
//             L1 = checkpoints[pos] + last_win + sub_size * n;

//             // 2. 开始组装每个窗口的结果
//             //double sum = 0.0;
//             //cout<<"windows_CMs size:"<<Window_CMs.size()<<endl;
//             sum += Window_CMs[0][pos].query(dic[k]);
//             for (int i = 1; i < sub_num; i++) {
//                 sum += Window_CMs[i][0].query(dic[k]);
//             }
//             fre_noise.push_back(sum + noise_vec[k]);

//         }

//         //double scale1 = 2*sub_size * log(w) / epsilon;
//         //cout<<"L1:"<<L1<<"   ";
//         //L1 += generateLaplacianRandom(scale1,1)[0];
//         //cout<<L1<<"   ";
//         double the = (1.0-zeta/2.0)*(gamma) *L1;
//         //cout<<"the:"<<the<<endl;
//         //double the1 = (1.0-zeta/2.0)*gamma *L1-100;
//         if (sum> the){
//             item_fre_noise.insert(pair<string,double> (dic[k],sum));
//         }

//     }

//     //double scale2 = 2 * sub_size / epsilon;
//     //vector<double> noise = generateLaplacianRandom(scale2,item_fre_noise.size());
//     // unordered_map<string,double>::iterator iter;
//     // int ix = 0;
//     // for(iter = item_fre_noise.begin();iter!=item_fre_noise.end();iter++){
//     //     //item_fre_noise.insert(pair<string,double>(iter->first,iter->second+noise[ix]));
//     //     iter->second += noise[ix];
//     //     ix ++;
//     // }
//     return make_pair(fre_noise,item_fre_noise);

// }



pair<vector<double>, unordered_map<string,double>> Private_Heavy::Query_all() {
    //vector<double> fre_no_noise;
    //unordered_map<string,double> item_fre_no_noise;
    generator.seed(currentSeed);
    vector<double> fre_noise;
    unordered_map<string,double> item_fre_noise;
    double L1;
    pair<vector<double>,unordered_map<string,double>> res;
    generator.seed(currentSeed);
    vector<double> noise_vec = generateLaplacianRandom(generator,1.0*sub_size / epsilon, dic.size());
    currentSeed++;
    for (int k=0;k<dic.size();k++){
        double sum = 0.0;
        if (flag ==0){

            //double sum = 0.0;
            if (last_win == 0 || indices.size() == 1) {
                // 直接把第一个checkpoints都加起来
                for (int i = 0; i < sub_num; i++) {
                    double r1 = Window_CMs[i][0].query(dic[k]);
                    sum += r1;
                    //cout<<sum<<"     ";
                }
                //cout<<sum;
                fre_noise.push_back(sum + noise_vec[k]);
                L1 = w;
            } else {
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

                L1 = checkpoints[pos] + last_win + sub_size * (sub_num - 2);

                //double sum = 0.0;
                sum += Window_CMs[0][pos].query(dic[k]);
                for (int i = 1; i < sub_num; i++) {
                    sum += Window_CMs[i][0].query(dic[k]);
                }
                fre_noise.push_back(sum + noise_vec[k]);

            }
        }else{
            int tmp1 =(w-last_win);
            int n = 0 ;
            while (tmp1>sub_size){
                tmp1 -= sub_size;
                n++;
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
            L1 = checkpoints[pos] + last_win + sub_size * n;


            sum += Window_CMs[0][pos].query(dic[k]);
            for (int i = 1; i < sub_num; i++) {
                sum += Window_CMs[i][0].query(dic[k]);
            }
            fre_noise.push_back(sum + noise_vec[k]);

        }

        double scale1 = 2*sub_size * log(w) / epsilon;
        //cout<<"L1:"<<L1<<"   ";
        generator.seed(currentSeed);
        L1 += generateLaplacianRandom(generator,scale1,1)[0];
        currentSeed++;
        //cout<<L1<<"   ";
        double the = (1.0-zeta/2.0)*(gamma) *L1 +0.005*( pow(w,3/4)/(epsilon*1.0))*log(1/beta);
        //cout<<"the:"<<the<<endl;
        //double the1 = (1.0-zeta/2.0)*gamma *L1-100;
        if (sum> the){
            item_fre_noise.insert(pair<string,double> (dic[k],sum));
        }

    }

    double scale2 = 2 * sub_size / epsilon;
    generator.seed(currentSeed);
    vector<double> noise = generateLaplacianRandom(generator,scale2,item_fre_noise.size());
    currentSeed++;
    unordered_map<string,double>::iterator iter;
    int ix = 0;
    for(iter = item_fre_noise.begin();iter!=item_fre_noise.end();iter++){
        iter->second += noise[ix];
        ix ++;
    }
    return make_pair(fre_noise,item_fre_noise);

}

vector<double> Private_Heavy::Query() {

    generator.seed(currentSeed);
    vector<double> fre_noise;
    vector<double> noise_vec = generateLaplacianRandom(generator,2.0*sub_size / epsilon, dic.size());
    currentSeed++;
    for (int k=0;k<dic.size();k++){
        double sum = 0.0;
        if (flag ==0){

            //double sum = 0.0;
            if (last_win == 0 || indices.size() == 1) {
                for (int i = 0; i < sub_num; i++) {
                    double r1 = Window_CMs[i][0].query(dic[k]);
                    sum += r1;
                    //cout<<sum<<"     ";
                }
                //cout<<sum;
                fre_noise.push_back(sum + noise_vec[k]);
            } else {
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


                sum += Window_CMs[0][pos].query(dic[k]);
                for (int i = 1; i < sub_num; i++) {
                    sum += Window_CMs[i][0].query(dic[k]);
                }
                fre_noise.push_back(sum + noise_vec[k]);

            }
        }else{
            int tmp1 =(w-last_win);
            int n = 0 ;
            while (tmp1>sub_size){
                tmp1 -= sub_size;
                n++;
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

            sum += Window_CMs[0][pos].query(dic[k]);
            for (int i = 1; i < sub_num; i++) {
                sum += Window_CMs[i][0].query(dic[k]);
            }
            fre_noise.push_back(sum + noise_vec[k]);

        }


    }
    return fre_noise;

}


vector<vector<double>> Private_Heavy::show_parameter() {
    vector<CountMinSketch> cms = Window_CMs[0];
    vector<vector<double >> res;
    for (int i=0;i<cms.size();i++){
        vector<double> tmp = cms[i].get_parameter();
        res.push_back(tmp);
    }
    return res;

}