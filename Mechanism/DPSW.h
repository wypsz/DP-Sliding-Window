//
// Created by 12293 on 2023/12/22.
//

#ifndef DP_SLIDING_WINDOW_DPSW_H
#define DP_SLIDING_WINDOW_DPSW_H

# include "new_alg.h"

class DPSW {
public:
    DPSW(const vector<string>& xw, int w,int step, int sub_num,double rho, double gamma, double beta,double q,double alpha,unsigned int seedGenerator,unsigned int hseed);
    vector<CountMinSketch> ProcessSubWindow(const vector<string>& newItems);
    void ProcessNew(string item);
    vector<vector<double>> show_parameter();
    double Query(string item);

private:
    int w;
    int step;
    double rho;
    double gamma;
    double beta;
    double q;
    double alpha;
    double rho1;
    int sub_size;
    int sub_num;
    double all_budget;
    vector<double> privacy_budget;
    vector<int> indices;
    vector<int> checkpoints;
    vector<vector<CountMinSketch>> Window_CMs;
    vector<CountMinSketch> building_CM;
    int last_win=0;
    vector<CountMinSketch> last_win_cms;
    int last_finish_index;
    //int flag;
    unsigned int seedGenerator;
    unsigned int hseed;
    //double rest_budget;
};



#endif //DP_SLIDING_WINDOW_DPSW_H