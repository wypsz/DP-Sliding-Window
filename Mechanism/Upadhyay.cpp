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
                             double gamma,double beta, double zeta,vector<string> dic,int flag,unsigned int initialSeed){
    this->w = w;
    this->epsilon = epsilon;
    this->step = step;
    this->gamma = gamma;
    this->beta = beta;
    this->zeta = zeta;
    this->sub_num = sub_num;
    this->dic = dic;
    this->sub_size = w/(this->sub_num);
    this->flag=flag;
    this->initialSeed = initialSeed;
    this->currentSeed = initialSeed;
    cout<<"hi"<<endl;

    SmoothHistogram sh(beta,this->sub_size,step);
    cout<<"hi "<<"sub_sum"<<this->sub_num<<" subsize:"<<this->sub_size<<endl;

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
    for (int i = 0; i < this->sub_num; i++) {
        vector<string> sub_win;
        for (int j=0;j<this->sub_size;j++){
            sub_win.push_back(xw[i*this->sub_size+j]);
        }
        sub_windows.push_back(sub_win);
    }
    if (flag == 1){
        vector<string> sub_win;
        for (int i=this->sub_num*this->sub_size;i<w;i++){
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

}

vector<CountMinSketch> Private_Heavy::ProcessSubWindow(const vector<std::string> &newItems) {
    vector<CountMinSketch> cms;
    int index=0;
    cout<<"len indices:"<<indices.size()<<endl;
    for (int i=0;i<newItems.size();i++){
        if (i == 0){
            CountMinSketch cm (gamma/exp(1),zeta,0,0,initialSeed);
            cm.update(newItems[i],1);
            cms.push_back(cm);
            index++;
        }else{
            for (int j=0;j<index;j++){
                cms[j].update(newItems[i]);
            }
            if (i+1 == indices[index]){
                CountMinSketch cm(gamma/exp(1),zeta,0,0,initialSeed);
                cm.update(newItems[i]);
                cms.push_back(cm);
                index++;
            }


        }
    }
    return cms;
}

void Private_Heavy::ProcessNew(std::string item) {
    if (last_win == 0){
        CountMinSketch cm(gamma/exp(1),zeta,0,0,initialSeed);
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

            CountMinSketch cm(gamma/exp(1),zeta,0,0,initialSeed);
            cm.update(item);
            last_win_cms.push_back(cm);
        }
        last_win++;
    }

    if (last_win == this->sub_size){
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



vector<double> Private_Heavy::Query() {

    generator.seed(currentSeed);
    vector<double> fre_noise;
    int sqrtw = static_cast<int>(ceil(sqrt(w)));
    vector<double> noise_vec = generateLaplacianRandom(generator,2.0*sqrtw / (epsilon*1.0), dic.size());
    currentSeed++;
    for (int k=0;k<dic.size();k++){
        double sum = 0.0;
        if (flag ==0){

            if (last_win == 0 || indices.size() == 1) {
                for (int i = 0; i < this->sub_num; i++) {
                    double r1 = Window_CMs[i][0].query(dic[k]);
                    sum += r1;
                }
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
                for (int i = 1; i < this->sub_num; i++) {
                    sum += Window_CMs[i][0].query(dic[k]);
                }
                fre_noise.push_back(sum + noise_vec[k]);

            }
        }else{
            int tmp1 =(w-last_win);
            int n = 0 ;
            while (tmp1>this->sub_size){
                tmp1 -= this->sub_size;
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
            for (int i = 1; i < this->sub_num; i++) {
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