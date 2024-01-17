//
// Created by 12293 on 2023/12/22.
//

#ifndef DP_SLIDING_WINDOW_DPSW_H
#define DP_SLIDING_WINDOW_DPSW_H

# include "new_alg.h"

class DPSW {
public:
    DPSW(const vector<string>& xw, int w,int step, int sub_num,double rho, double gamma, double beta,double q,double alpha,unsigned int seedGenerator);
    vector<CountMinSketch> ProcessSubWindow(const vector<string>& newItems);
    void ProcessNew(string item);
    vector<vector<double>> show_parameter();
    double Query(string item);

private:
    int w; //窗口大小
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
    double rest_budget;
    // vector<string> dic;
    vector<int> indices;
    vector<int> checkpoints;
    vector<vector<CountMinSketch>> Window_CMs; // 所有的子窗口中CM汇总
    vector<CountMinSketch> building_CM;// 每个子窗口中构建好的最长的CMs
    int last_win=0;//统计最后一个子窗口内有多少item
    vector<CountMinSketch> last_win_cms;
    int last_finish_index;
    int flag;
    unsigned int seedGenerator;
};



#endif //DP_SLIDING_WINDOW_DPSW_H
