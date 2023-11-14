#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <vector>
#include <stdio.h>
#include <cstdlib>
#include <algorithm>
#include <iomanip>
#include "Mechanism/new_alg.h"
#include "Mechanism/Retroactive_Grouping.h"
#include "Mechanism/PCC.h"
using namespace std;

void WriteVectorToCSV(std::ofstream& outFile, const std::vector<double>& vec) {
    for (size_t i = 0; i < vec.size(); ++i) {
        outFile << vec[i];
        if (i != vec.size() - 1) outFile << ","; // 最后一个值后不加逗号
    }
    outFile << std::endl;
}

void WriteVectorToCSVint(std::ofstream& outFile, const std::vector<int>& vec) {
    for (size_t i = 0; i < vec.size(); ++i) {
        outFile << vec[i];
        if (i != vec.size() - 1) outFile << ","; // 最后一个值后不加逗号
    }
    outFile << std::endl;
}

vector<int> find_index_dic(vector<string>& items,vector<string>& dic){
    vector<int> res;
    for (int i=0;i<items.size();i++){
        auto it = find(dic.begin(),dic.end(),items[i]);
        if (it != dic.end()){
            int index = distance(dic.begin(),it);
            res.push_back(index);
        }else{
            cout<<"Element not found in dic."<<endl;
        }
    }
    return res;
}

std::vector<std::vector<std::string>> loadFromFile(const std::string& filename) {
    std::vector<std::vector<std::string>> data;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        return data;
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        std::istringstream ss(line);
        std::string item;
        std::vector<std::string> row;

        while (std::getline(ss, item, ',')) {
            try{
                row.push_back(item);
            }catch (const exception& e){
                std::cerr << "Error parsing number from '" << item << "': " << e.what() << std::endl;

            }

        }

        data.push_back(row);
    }

    file.close();
    return data;
}


std::vector<std::vector<int>> loadFromFile_num(const std::string& filename) {
    std::vector<std::vector<int>> data;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        return data;
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        std::istringstream ss(line);
        std::string item;
        std::vector<int> row;

        while (std::getline(ss, item, ',')) {
            try{
                row.push_back(std::stoi(item)); // Convert string to int
            }catch (const exception& e){
                std::cerr << "Error parsing number from '" << item << "': " << e.what() << std::endl;

            }
        }

        data.push_back(row);
    }

    file.close();
    return data;
}
double MAE_noise_true(vector<int> frequency, vector<double> frequency_noise){
    // 平均绝对误差
    int size = frequency.size();
    if (size!=frequency_noise.size()){
        cout<<size<<","<<frequency_noise.size()<<endl;
        cout<<"error!"<<endl;
    }
    double sum = 0.0;
    for (int i=0;i<size;i++){
        sum += fabs(frequency_noise[i]-frequency[i]);
    }
    if (isnan(sum)){
        for (int i=0;i<size;i++){
            cout<<"fre_noise"<<frequency_noise[i]<<",";
        }
        cout<<endl;
        for (int i=0;i<size;i++){
            cout<<"fre:"<<frequency[i]<<",";
        }
        cout<<endl;
    }
    //cout<<"sum:"<<sum<<endl;
    //cout<<"size:"<<size<<endl;
    return sum/size;
}

double MAE_cm_noise(vector<double> frequency_cm, vector<double> frequency_noise){
    int size = frequency_cm.size();
    if (size!=frequency_noise.size()){
        cout<<"error!"<<endl;
    }
    double sum = 0.0;
    for (int i=0;i<size;i++){
        sum = sum + fabs(frequency_noise[i]-frequency_cm[i]);
    }

    if (isnan(sum)){
        for (int i=0;i<size;i++){
            cout<<"fre_noise"<<frequency_noise[i]<<",";
        }
        cout<<endl;
        for (int i=0;i<size;i++){
            cout<<"fre:"<<frequency_cm[i]<<",";
        }
        cout<<endl;
    }
    //cout<<"double"<<endl;
    return sum/size;
}


double MRE_noise_true(vector<int> frequency, vector<double> frequency_noise){
    int size = frequency.size();
    if (size!=frequency_noise.size()){
        cout<<"error!"<<endl;
    }
    double sum = 0.0;
    for (int i=0;i<size;i++){
        sum += (fabs(frequency_noise[i]-frequency[i]))/frequency[i];
    }

    if (isnan(sum)){
        for (int i=0;i<size;i++){
            cout<<"fre_noise"<<frequency_noise[i]<<",";
        }
        cout<<endl;
        for (int i=0;i<size;i++){
            cout<<"fre:"<<frequency[i]<<",";
        }
        cout<<endl;
    }
    return sum/size;
}

double MRE_cm_noise(vector<double>& frequency, vector<double>& frequency_noise){
    int size = frequency.size();
    if (size!=frequency_noise.size()){
        cout<<"error!"<<endl;
    }
    double sum = 0.0;
    for (int i=0;i<size;i++){
        sum += (fabs(frequency_noise[i]-frequency[i]))/frequency[i];
    }

    if (isnan(sum)){
        for (int i=0;i<size;i++){
            cout<<"fre_noise"<<frequency_noise[i]<<",";
        }
        cout<<endl;
        for (int i=0;i<size;i++){
            cout<<"fre:"<<frequency[i]<<",";
        }
        cout<<endl;
    }
    return sum/size;
}

//vector<int> true_frequency(vector<string>& part){
//
//}


//double ARE_noise_true(vector<int> frequency, vector<double> frequency_noise){
//    int size = frequency.size();
//    if (size!=frequency_noise.size()){
//        cout<<"error!"<<endl;
//    }
//    double sum = 0.0;
//    for (int i=0;i<size;i++){
//        sum += (fabs(frequency_noise[i]-frequency[i]))/frequency[i];
//    }
//    return sum;
//}
//
//
//double ARE_cm_noise(vector<double> frequency, vector<double> frequency_noise){
//    int size = frequency.size();
//    if (size!=frequency_noise.size()){
//        cout<<"error!"<<endl;
//    }
//    double sum = 0.0;
//    for (int i=0;i<size;i++){
//        sum += (fabs(frequency_noise[i]-frequency[i]))/frequency[i];
//    }
//    return sum;
//}
// 10*10w
// checkpoint 中100
int main() {



//(double epsilon, double lambda, int W)
/*
 * 测试PCC
    vector<string> items = {"1","2","3","4","5","6","7"};
    PCC pcc(1.0,0.5,7,items);
    unordered_map<std::string, double> res;

    pcc.processItem("7");
    pcc.processItem("7");

    res = pcc.query_all();
    for (auto  [element,freq]:res){
        cout<<element<<":"<<freq<<endl;
    }
 */

//    for (auto  [element,freq]:res){
//        cout<<element<<":"<<freq<<endl;
//    }

//    vector<vector<string>> High_item = {{"1","1"}};
//    vector<vector<int>> High_freq = {{1,2},{1,2}};
//    vector<vector<string>> Middle_item = {{"1","1"}};
//    vector<vector<int>> Middle_freq = {{1,2},{1,2}};
//    vector<vector<string>> Low_item = {{"1","1"}};
//    vector<vector<int>> Low_freq = {{1,2},{1,2}};
// * new_alg 和 RG一起测试
    //ifstream inFile("C://Users/12293/CLionProjects/DP_Sliding_window/dataset/sx-stackoverflow-c2q.txt");
    ifstream inFile("C://Users/12293/CLionProjects/DP_Sliding_window/dataset/worldCup.txt");
    if (!inFile) {
        cerr << "failed" << endl;
        return 1;
    }

    // 读取txt文件中的数据并对数据进行处理
    string line;
    vector<string> mini_data;
    int x = 0;
    while (getline(inFile, line)) {
        // 对读取到的每行数据进行处理
        stringstream ls;
        string t;
        ls.str(line);
        ls>>t;
        mini_data.push_back(t);
        x++;
//        if (x == 20000000){
//            break;
//        }
       if (x == 30000000){
            break;
        }
//        if (x == 3000){
//            break;
//        }
    }
    // 关闭txt文件
    inFile.close();

    cout<<mini_data.size()<<endl;


    // 测试两个完整窗口，开始准备参数
    // 10个子窗口
//    int w = 1000000;
//    int sub_num = 10;
//    int step = 1000;
//    double all_budget = 1;
//    double rho = 0.799;
//    double gamma = 0.0001;
//    double beta = 0.01;
//    double alpha = 0.3;
//    double q = 0.2;


//    int w = 1000;
//    int sub_num = 10;
//    int step = 10;
//    double rho = 0.8;
//    double gamma = 0.01;
//    double beta = 0.1;
//    double alpha = 0.1;
//    double q = 0.8;


//    int w = 1000.01m = 10;
//    int step = 10000;
//    double all_budget = 1;
//    double rho = 0.799;
//    double gamma = 0.000005;
//    double beta = 0.01;
//    double alpha = 0.9;
//    double q = 0.2;

    int w = 10000000;
    int sub_num = 10;
    int step = 10000;
    double all_budget = 1;
    double rho = 0.6999;
    double gamma = 0.0005;
    double beta = 0.0002;
    double alpha = 0.5;
    double q = 0.3;

    std::ofstream outFile("C://Users/12293/CLionProjects/DP_Sliding_window/res/2000bit/9hash.csv");
    if (!outFile.is_open()) {
        std::cerr << "Unable to open file";
        return 1;
    }
    outFile << "w = "<< w<<" , "<< "sub_num = "<< sub_num<<" , "<< "step = "<< step<<" , "<< "all_budget = "<< all_budget<<" , "<< "rho1 = "<< rho<<" , "<< "gamma = "<< gamma<<" , "<< "beta = "<< beta<<" , "<< "alpha = "<<alpha<<" , "<< "q = "<<q<<endl;

    vector<vector<string>> High_item = loadFromFile("C://Users/12293/CLionProjects/DP_Sliding_window/dataset/High_item.csv");
//    for (int x=0;x<High_item[0].size();x++){
//        cout<<High_item[0][x]<<";";
//    }
    vector<vector<string>> Middle_item = loadFromFile("C://Users/12293/CLionProjects/DP_Sliding_window/dataset/Middle_item.csv");
    vector<vector<string>> Low_item = loadFromFile("C://Users/12293/CLionProjects/DP_Sliding_window/dataset/Low_item.csv");
    vector<vector<int>> High_freq = loadFromFile_num("C://Users/12293/CLionProjects/DP_Sliding_window/dataset/High_freq.csv");
    vector<vector<int>> Middle_freq = loadFromFile_num("C://Users/12293/CLionProjects/DP_Sliding_window/dataset/Middle_freq.csv");
    vector<vector<int>> Low_freq = loadFromFile_num("C://Users/12293/CLionProjects/DP_Sliding_window/dataset/Low_freq.csv");
    // 准备dic
    // 对于RG算法，需要统计所有元素的dic
    vector<string> dic;
    set<string> d_all(mini_data.begin(),mini_data.end());
    dic.assign(d_all.begin(),d_all.end());
    cout<<"all dic size:"<<dic.size()<<endl;

    // 统计第一个窗口中的dic
//    vector<string> dic;
//    set<string> d(mini_data.begin(),mini_data.begin()+w);
//    dic.assign(d.begin(),d.end());
//    cout<<"the first window size:"<<dic.size()<<endl;

    //cout<<mini_data.size();
    //vector<string> tmp(mini_data.begin(),mini_data.begin()+w);
    cout<<"hi"<<endl;
    Private_counter pc (vector<string> (mini_data.begin(),mini_data.begin()+w), w, step, sub_num, rho, gamma, beta, q, alpha,all_budget);
    Private_counter pc_no_noise(vector<string> (mini_data.begin(),mini_data.begin()+w), w, step, sub_num, 0, gamma, beta, q, alpha,0);


    cout<<"Private CS and CS first window initialize finish"<<endl;

    SmoothHistogram sh(alpha,w/sub_num,step);
    sh.Add();
    vector<int> indices = sh.getIndices();
    vector<int>checkpoints = sh.getCheckpoints();
    outFile<<"indices:"<<endl;
    WriteVectorToCSVint(outFile,indices);
    outFile<<"checkpoints:"<<endl;
    WriteVectorToCSVint(outFile,checkpoints);

    vector<vector<double>> parameter = pc.show_parameter();

//    for (int a=0;a<parameter.size();a++){
//        for (int b=0;b<parameter[0].size();b++){
//            outFile << parameter[a][b];
//            if (a!=parameter.size()-1 && b!= parameter[0].size()-1){
//                outFile<<",";
//            }
//        }
//    }
//    outFile<<endl;
    outFile<<"Parameters in sub_window"<<endl;
    for (int a=0;a<parameter.size();a++){
        WriteVectorToCSV(outFile,parameter[a]);
    }
    cout<<"---------------------------------------------------------------------------------"<<endl;
    /*
    cout<<"RG begin"<<endl;

    SmoothHistogram sh(0.99,w/sub_num,step);
    sh.Add();
    vector<int> indices = sh.getIndices();
    vector<int> checkpoints = sh.getCheckpoints();
    indices.push_back(w/sub_num);
    checkpoints.push_back(1);
    cout<<"indices:"<<endl;
    for (int i=0;i<indices.size();i++){
        indices[i] -= 1;
        cout<<indices[i]<<"     ";
    }
    cout<<endl;
    cout<<"checkpoints:"<<endl;
    for (int i=0;i<checkpoints.size();i++){
        checkpoints[i]-= 1;
        cout<<checkpoints[i]<<"     ";
    }
    cout<<endl;
    cout<<checkpoints.size()<<"=="<<indices.size();
    cout<<"hello!"<<endl;
    double epsilon1 = 0.1;
    double epsilon2 = 0.9;
    Retroactive_Grouping RG( vector<string> (mini_data.begin(),mini_data.begin()+w), indices,epsilon1, epsilon2,dic, w,10);
    */

    vector<double> MAE_Noise_True_H;
    vector<double> MAE_Noise_True_M;
    vector<double> MAE_Noise_True_L;
    vector<double> MAE_Noise_True_sum;

    vector<double> MAE_Cm_Noise_H;
    vector<double> MAE_Cm_Noise_M;
    vector<double> MAE_Cm_Noise_L;
    vector<double> MAE_Cm_Noise_sum;

    vector<double> MAE_Cm_True_H;
    vector<double> MAE_Cm_True_M;
    vector<double> MAE_Cm_True_L;
    vector<double> MAE_Cm_True_sum;

    vector<double> MRE_Noise_True_H;
    vector<double> MRE_Noise_True_M;
    vector<double> MRE_Noise_True_L;
    vector<double> MRE_Noise_True_sum;

    vector<double> MRE_Cm_Noise_H;
    vector<double> MRE_Cm_Noise_M;
    vector<double> MRE_Cm_Noise_L;
    vector<double> MRE_Cm_Noise_sum;

    vector<double> MRE_Cm_True_H;
    vector<double> MRE_Cm_True_M;
    vector<double> MRE_Cm_True_L;
    vector<double> MRE_Cm_True_sum;

    vector<double> MAE_Noise_True_RG;
    vector<double> MRE_Noise_True_RG;

    cout<<endl;
    cout<<"--------------------------------------------------------------------------------"<<endl;
    cout<<"begin test window"<<endl;
    //vector<vector<int>> query_list;
    //int now = w;
    //int query_space = 20000;
    //int query_space = 50;
    //int query_space = 50000;
    int query_space = 50000; // 子窗口大小为100w
    int ind = 0; // 标记在indices中的位置
    int index = 0;
    for (int i=w;i<mini_data.size();i=i+query_space) {
        //now = i+query_space;
        // 添加了10000个元素，查询一次
//        vector<int> items;
//        vector<int> frequency;
//        vector<double> frequency_cm;
//
//        vector<double> frequency_noise;
        // RG取消注释下面一行
        //vector<double> noise_fre_rg;
        cout<<"Begin test:"<<endl;
        /*
        for (int k=0;k<100;k++){
            int random = (rand()%(dic.size()-0+1));
             // cout<<random<<endl;
            // 查找元素random的真实频率
            //cout<<"find ture fre"<<now-w<<","<<now<<endl;
            int fre = std::count(mini_data.begin()+now-w, mini_data.begin()+now, dic[random]);
            //cout<<"Bad"<<endl;
            if (fre != 0){
                items.push_back(random);
                frequency.push_back(fre);
                //cout<<"true fre:"<<fre<<endl;
                // 查找元素random不加噪的频率
                double fre_cm = pc_no_noise.Query(dic[random]);
                //cout<<"find cm fre"<<endl;
                frequency_cm.push_back(fre_cm);
                // 查找元素random加噪后的频率
                double fre_noise = pc.Query(dic[random]);
                frequency_noise.push_back(fre_noise);
                //cout<<"find noise fre:"<<fre_noise<<endl;
                // RG的下面两行取消注释
                //double noise_f = RG.query(dic[random]);
                //noise_fre_rg.push_back(noise_f);
            }else{
                k--;
            }
        }
         */
        // 选取const中的第i组高、中、低频元素进行查询
        cout<<"hi:"<<index<<endl;
        vector<string> hh_i = High_item[index];
        vector<int> H_item_index = find_index_dic(hh_i,dic);
        vector<string> h_i;
        for (int x=0;x<H_item_index.size();x++){
            h_i.push_back(dic[H_item_index[x]]);
        }
        auto mm_i = Middle_item[index];
        vector<int> M_item_index = find_index_dic(mm_i,dic);
        vector<string> m_i;
        for (int x=0;x<M_item_index.size();x++){
            m_i.push_back(dic[M_item_index[x]]);
        }

        auto ll_i = Low_item[index];
        vector<int> L_item_index = find_index_dic(ll_i,dic);
        vector<string> l_i;
        for (int x=0;x<L_item_index.size();x++){
            l_i.push_back(dic[L_item_index[x]]);
        }
        // 高中低频元素的真实frequency
        auto H_f = High_freq[index];
        auto M_f = Middle_freq[index];
        auto L_f = Low_freq[index];
        index++;

        // 高中低频元素查询的frequency
        vector<double> freq_noise_h;
        vector<double> freq_cm_h;
        vector<double> freq_noise_m;
        vector<double> freq_cm_m;
        vector<double> freq_noise_l;
        vector<double> freq_cm_l;
        //cout<<"h2"<<endl;
        for (int k=0;k<h_i.size();k++){
            // 查询sketch
            //cout<<"h3"<<endl;
            double a_h = pc_no_noise.Query(h_i[k]);
            double a_m = pc_no_noise.Query(m_i[k]);
            double a_l = pc_no_noise.Query(l_i[k]);
            //cout<<"h4"<<endl;
            freq_cm_h.push_back(a_h);
            freq_cm_m.push_back(a_m);
            freq_cm_l.push_back(a_l);
            // 查询加噪后的结果
            double b_h = pc.Query(h_i[k]);
            double b_m = pc.Query(m_i[k]);
            double b_l = pc.Query(l_i[k]);
            freq_noise_h.push_back(b_h);
            freq_noise_m.push_back(b_m);
            freq_noise_l.push_back(b_l);
        }
        //cout<<"items:"<<items.size()<<endl;
        //query_list.push_back(items);
        cout<<i<<endl;
        // 高频结果
        double mae_noise_true_h = MAE_noise_true(H_f,freq_noise_h);
        double mae_cm_noise_h = MAE_cm_noise(freq_noise_h,freq_cm_h);
        double mae_cm_true_h = MAE_noise_true(H_f,freq_cm_h);
        double mre_noise_true_h = MRE_noise_true(H_f,freq_noise_h);
        double mre_cm_noise_h = MRE_cm_noise(freq_cm_h,freq_noise_h);
        double mre_cm_true_h = MRE_noise_true(H_f,freq_cm_h);
        // 中频结果
        double mae_noise_true_m = MAE_noise_true(M_f,freq_noise_m);
        double mae_cm_noise_m = MAE_cm_noise(freq_noise_m,freq_cm_m);
        double mae_cm_true_m = MAE_noise_true(M_f,freq_cm_m);
        double mre_noise_true_m = MRE_noise_true(M_f,freq_noise_m);
        double mre_cm_noise_m = MRE_cm_noise(freq_cm_m,freq_noise_m);
        double mre_cm_true_m = MRE_noise_true(M_f,freq_cm_m);
        // 低频结果
        double mae_noise_true_l = MAE_noise_true(L_f,freq_noise_l);
        double mae_cm_noise_l = MAE_cm_noise(freq_noise_l,freq_cm_l);
        double mae_cm_true_l = MAE_noise_true(L_f,freq_cm_l);
        double mre_noise_true_l = MRE_noise_true(L_f,freq_noise_l);
        double mre_cm_noise_l = MRE_cm_noise(freq_cm_l,freq_noise_l);
        double mre_cm_true_l = MRE_noise_true(L_f,freq_cm_l);
        // 总的结果，三个值求和再除3
        double mae_noise_true_sum = (mae_noise_true_h+mae_noise_true_m+mae_noise_true_l)/3;
        double mae_cm_noise_sum = (mae_cm_noise_h+mae_cm_noise_m+mae_cm_noise_l)/3;
        double mae_cm_true_sum = (mae_cm_true_h+mae_cm_true_m+mae_cm_true_l)/3;
        double mre_noise_true_sum = (mre_noise_true_h+mre_noise_true_m+mre_noise_true_l)/3;
        double mre_cm_noise_sum = (mre_cm_noise_h+mre_cm_noise_m+mre_cm_noise_l)/3;
        double mre_cm_true_sum = (mre_cm_true_h+mre_cm_true_m+mre_cm_true_l)/3;

        //RG注释下面两行
        //double mae_RG = MAE_noise_true(frequency,noise_fre_rg);
        //double mre_RG = MRE_noise_true(frequency,noise_fre_rg);
        MAE_Noise_True_H.push_back(mae_noise_true_h);
        MAE_Noise_True_M.push_back(mae_noise_true_m);
        MAE_Noise_True_L.push_back(mae_noise_true_l);
        MAE_Noise_True_sum.push_back(mae_noise_true_sum);

        MAE_Cm_Noise_H.push_back(mae_cm_noise_h);
        MAE_Cm_Noise_M.push_back(mae_cm_noise_m);
        MAE_Cm_Noise_L.push_back(mae_cm_noise_l);
        MAE_Cm_Noise_sum.push_back(mae_cm_noise_sum);

        MAE_Cm_True_H.push_back(mae_cm_true_h);
        MAE_Cm_True_M.push_back(mae_cm_true_m);
        MAE_Cm_True_L.push_back(mae_cm_true_l);
        MAE_Cm_True_sum.push_back(mae_cm_true_sum);

        MRE_Noise_True_H.push_back(mre_noise_true_h);
        MRE_Noise_True_M.push_back(mre_noise_true_m);
        MRE_Noise_True_L.push_back(mre_noise_true_l);
        MRE_Noise_True_sum.push_back(mre_noise_true_sum);

        MRE_Cm_Noise_H.push_back(mre_cm_noise_h);
        MRE_Cm_Noise_M.push_back(mre_cm_noise_m);
        MRE_Cm_Noise_L.push_back(mre_cm_noise_l);
        MRE_Cm_Noise_sum.push_back(mre_cm_noise_sum);

        MRE_Cm_True_H.push_back(mre_cm_true_h);
        MRE_Cm_True_M.push_back(mre_cm_true_m);
        MRE_Cm_True_L.push_back(mre_cm_true_l);
        MRE_Cm_True_sum.push_back(mre_cm_true_sum);
        for (int j=0;j<query_space;j++){
            pc.ProcessNew(mini_data[i+j]);
            //cout<<"Process:"<<i+j<<"th item"<<endl;
            // RG取消注释下面一行
            //RG.insert(mini_data[i+j]);
        }
        //RG取消注释下面两行
        //MAE_Noise_True_RG.push_back(mae_RG);
        //MRE_Noise_True_RG.push_back(mre_RG);
//        for(int x = 0;x<noise_fre_rg.size();x++){
//            cout<<"noise RG:"<<noise_fre_rg[x]<<endl;
//            cout<<"true fre:"<<frequency[x]<<endl;
//            cout<<"MRE RG:"<<mre_RG<<endl;
//        }

    }
    cout<<"High : MAE between noise and true:"<<endl;
    for(int i=0;i<MAE_Noise_True_H.size();i++){
        cout<<MAE_Noise_True_H[i]<<",";
    }
    cout<<endl;

    outFile<<"High : MAE between noise and true:"<<endl;
    WriteVectorToCSV(outFile,MAE_Noise_True_H);

    cout<<"Middle : MAE between noise and true:"<<endl;
    for(int i=0;i<MAE_Noise_True_M.size();i++){
        cout<<MAE_Noise_True_M[i]<<",";
    }
    cout<<endl;

    outFile<<"Middle : MAE between noise and true:"<<endl;
    WriteVectorToCSV(outFile,MAE_Noise_True_M);

    cout<<"Low : MAE between noise and true:"<<endl;
    for(int i=0;i<MAE_Noise_True_L.size();i++){
        cout<<MAE_Noise_True_L[i]<<",";
    }
    cout<<endl;
    outFile<<"Low : MAE between noise and true:"<<endl;
    WriteVectorToCSV(outFile,MAE_Noise_True_L);


    cout<<"sum : MAE between noise and true:"<<endl;
    for(int i=0;i<MAE_Noise_True_sum.size();i++){
        cout<<MAE_Noise_True_sum[i]<<",";
    }
    cout<<endl;
    outFile<<"sum : MAE between noise and true:"<<endl;
    WriteVectorToCSV(outFile,MAE_Noise_True_sum);
    cout<<"---------------------------------------------------------------------------------"<<endl;
    cout<<"High:MAE between noise and cm:"<<endl;
    for(int i=0;i<MAE_Cm_Noise_H.size();i++){
        cout<<MAE_Cm_Noise_H[i]<<",";
    }
    cout<<endl;

    outFile<<"High:MAE between noise and cm:"<<endl;
    WriteVectorToCSV(outFile,MAE_Cm_Noise_H);

    cout<<"Middle:MAE between noise and cm:"<<endl;
    for(int i=0;i<MAE_Cm_Noise_M.size();i++){
        cout<<MAE_Cm_Noise_M[i]<<",";
    }
    cout<<endl;
    outFile<<"Middle:MAE between noise and cm:"<<endl;
    WriteVectorToCSV(outFile,MAE_Cm_Noise_M);


    cout<<"Low : MAE between noise and cm:"<<endl;
    for(int i=0;i<MAE_Cm_Noise_L.size();i++){
        cout<<MAE_Cm_Noise_L[i]<<",";
    }
    cout<<endl;
    outFile<<"Low : MAE between noise and cm:"<<endl;
    WriteVectorToCSV(outFile,MAE_Cm_Noise_L);

    cout<<"sum :MAE between noise and cm:"<<endl;
    for(int i=0;i<MAE_Cm_Noise_sum.size();i++){
        cout<<MAE_Cm_Noise_sum[i]<<",";
    }
    cout<<endl;
    outFile<<"sum :MAE between noise and cm:"<<endl;
    WriteVectorToCSV(outFile,MAE_Cm_Noise_sum);
    cout<<"----------------------------------------------------------------"<<endl;
    cout<<"High : MAE between cm and true:"<<endl;
    for(int i=0;i<MAE_Cm_True_H.size();i++){
        cout<<MAE_Cm_True_H[i]<<",";
    }
    cout<<endl;
    outFile<<"High : MAE between cm and true:"<<endl;
    WriteVectorToCSV(outFile,MAE_Cm_True_H);


    cout<<"Middle : MAE between cm and true:"<<endl;
    for(int i=0;i<MAE_Cm_True_M.size();i++){
        cout<<MAE_Cm_True_M[i]<<",";
    }
    cout<<endl;
    outFile<<"Middle : MAE between cm and true:"<<endl;
    WriteVectorToCSV(outFile,MAE_Cm_True_M);

    cout<<"Low : MAE between cm and true:"<<endl;
    for(int i=0;i<MAE_Cm_True_L.size();i++){
        cout<<MAE_Cm_True_L[i]<<",";
    }
    cout<<endl;
    outFile<<"Low : MAE between cm and true:"<<endl;
    WriteVectorToCSV(outFile,MAE_Cm_True_L);


    cout<<"sum : MAE between cm and true:"<<endl;
    for(int i=0;i<MAE_Cm_True_sum.size();i++){
        cout<<MAE_Cm_True_sum[i]<<",";
    }
    cout<<endl;
    outFile<<"sum : MAE between cm and true:"<<endl;
    WriteVectorToCSV(outFile,MAE_Cm_True_sum);


    cout<<"-------------------------------------------------------"<<endl;


    cout<<"High : MRE between noise and true:"<<endl;
    for(int i=0;i<MRE_Noise_True_H.size();i++){
        cout<<MRE_Noise_True_H[i]<<",";
    }
    cout<<endl;
    outFile<<"High : MRE between noise and true:"<<endl;
    WriteVectorToCSV(outFile,MRE_Noise_True_H);


    cout<<"Middle : MRE between noise and true:"<<endl;
    for(int i=0;i<MRE_Noise_True_M.size();i++){
        cout<<MRE_Noise_True_M[i]<<",";
    }
    cout<<endl;
    outFile<<"Middle : MRE between noise and true:"<<endl;
    WriteVectorToCSV(outFile,MRE_Noise_True_M);


    cout<<"Low : MRE between noise and true:"<<endl;
    for(int i=0;i<MRE_Noise_True_L.size();i++){
        cout<<MRE_Noise_True_L[i]<<",";
    }
    cout<<endl;
    outFile<<"Low : MRE between noise and true:"<<endl;
    WriteVectorToCSV(outFile,MRE_Noise_True_L);


    cout<<"sum : MRE between noise and true:"<<endl;
    for(int i=0;i<MRE_Noise_True_sum.size();i++){
        cout<<MRE_Noise_True_sum[i]<<",";
    }
    cout<<endl;
    outFile<<"sum : MRE between noise and true:"<<endl;
    WriteVectorToCSV(outFile,MRE_Noise_True_sum);
    cout<<"----------------------------------------------------------"<<endl;


    cout<<"High : MRE between noise and cm:"<<endl;
    for(int i=0;i<MRE_Cm_Noise_H.size();i++){
        cout<<MRE_Cm_Noise_H[i]<<",";
    }
    cout<<endl;
    outFile<<"High : MRE between noise and cm:"<<endl;
    WriteVectorToCSV(outFile,MRE_Cm_Noise_H);

    cout<<"Middle : MRE between noise and cm:"<<endl;
    for(int i=0;i<MRE_Cm_Noise_M.size();i++){
        cout<<MRE_Cm_Noise_M[i]<<",";
    }
    cout<<endl;
    outFile<<"Middle : MRE between noise and cm:"<<endl;
    WriteVectorToCSV(outFile,MRE_Cm_Noise_M);


    cout<<"Low : MRE between noise and cm:"<<endl;
    for(int i=0;i<MRE_Cm_Noise_L.size();i++){
        cout<<MRE_Cm_Noise_L[i]<<",";
    }
    cout<<endl;
    outFile<<"Low : MRE between noise and cm:"<<endl;
    WriteVectorToCSV(outFile,MRE_Cm_Noise_L);

    cout<<"sum : MRE between noise and cm:"<<endl;
    for(int i=0;i<MRE_Cm_Noise_sum.size();i++){
        cout<<MRE_Cm_Noise_sum[i]<<",";
    }
    cout<<endl;
    outFile<<"sum : MRE between noise and cm:"<<endl;
    WriteVectorToCSV(outFile,MRE_Cm_Noise_sum);
    cout<<"--------------------------------------------------------"<<endl;
    cout<<"High:MRE between cm and true:"<<endl;
    for(int i=0;i<MRE_Cm_True_H.size();i++){
        cout<<MRE_Cm_True_H[i]<<",";
    }
    cout<<endl;
    outFile<<"High:MRE between cm and true:"<<endl;
    WriteVectorToCSV(outFile,MRE_Cm_True_H);


    cout<<"Middle:MRE between cm and true:"<<endl;
    for(int i=0;i<MRE_Cm_True_M.size();i++){
        cout<<MRE_Cm_True_M[i]<<",";
    }
    cout<<endl;
    outFile<<"Middle:MRE between cm and true:"<<endl;
    WriteVectorToCSV(outFile,MRE_Cm_True_M);


    cout<<"Low:MRE between cm and true:"<<endl;
    for(int i=0;i<MRE_Cm_True_L.size();i++){
        cout<<MRE_Cm_True_L[i]<<",";
    }
    cout<<endl;
    outFile<<"Low:MRE between cm and true:"<<endl;
    WriteVectorToCSV(outFile,MRE_Cm_True_L);


    cout<<"sum:MRE between cm and true:"<<endl;
    for(int i=0;i<MRE_Cm_True_sum.size();i++){
        cout<<MRE_Cm_True_sum[i]<<",";
    }
    cout<<endl;
    outFile<<"sum:MRE between cm and true:"<<endl;
    WriteVectorToCSV(outFile,MRE_Cm_True_sum);

    outFile.close();


//    cout<<"query_list"<<endl;
//    for (int i=0;i<query_list.size();i++){
//        for (int j=0;j<query_list[i].size();j++){
//                cout<<query_list[i][j]<<",";
//
//        }
//        cout<<endl;
//    }
    /*
    cout<<"MAE between RG and true:"<<endl;
    for(int i=0;i<MAE_Noise_True_RG.size();i++){
        cout<<MAE_Noise_True_RG[i]<<",";
    }

    cout<<endl;

    cout<<"MRE between RG and true:"<<endl;
    for(int i=0;i<MRE_Noise_True_RG.size();i++) {
        cout << MRE_Noise_True_RG[i] << ",";
    }

*/


/*
    测试RG
    // little data 选取1w2条，分成100个大小为100的子窗口,窗口大小为1w，相当于窗口向后滑动20次
    // 先看看元素的所有情况
    // 去重获得dic
    set<string> d(mini_data.begin(),mini_data.end());
    vector<string> dic;
    dic.assign(d.begin(),d.end());
    cout<<"dic-size"<<dic.size()<<endl;

    SmoothHistogram sh(0.1,100,10);
    sh.Add();
    vector<int> indices = sh.getIndices();
    vector<int> checkpoints = sh.getCheckpoints();
    indices.push_back(100);
    checkpoints.push_back(1);
    cout<<"indices:"<<endl;
    for (int i=0;i<indices.size();i++){
        indices[i] -= 1;
        cout<<indices[i]<<"     ";
    }
    cout<<endl;
    cout<<"checkpoints:"<<endl;
    for (int i=0;i<checkpoints.size();i++){
        checkpoints[i]-= 1;
        cout<<checkpoints[i]<<"     ";
    }
    cout<<endl;
    cout<<checkpoints.size()<<"=="<<indices.size();
    cout<<"hello!"<<endl;
    int w = 1000;
    double epsilon1 = 0.1;
    double epsilon2 = 0.9;
    Retroactive_Grouping RG(mini_data, indices,epsilon1, epsilon2,dic, w,10);
    vector<double> res;
    res = RG.query();

    RG.insert(vector<string>(mini_data.begin()+w,mini_data.begin()+w+16));
    res = RG.query();
    for (int i=0;i<res.size();i++){
        cout<<res[i]<<"           ";
    }

*/


/*
 * 测试New——alg
     //开始测试，准备参数
     cout<<"hello!"<<endl;
     int w = 1000000;
     int sub_num = 10;// 子窗口个数
     Private_counter pc(mini_data,w,100,sub_num,0.8, 0.001,0.01,0.8,0.05);
     Private_counter pc_no_noise(mini_data,w,100,sub_num,0, 0.001,0.01,0.8,0.05);
     // 第一个大窗口查找结果
    // 按照dic中的顺序统计每个item准确的frequency
    cout<<"finish1"<<endl;


    // 生成随机数，随机查询frequency
    vector<int> items;
    vector<int> frequency;
    vector<double> frequency_cm;
    vector<double> frequency_noise;
    for (int i=0;i<100;i++){
        int random = (rand()%(dic.size()-0+1));
        // cout<<random<<endl;
        // 查找元素random的真实频率
        int fre = std::count(mini_data.begin(), mini_data.begin()+w, dic[random]);
        if (fre != 0){
            items.push_back(random);
            frequency.push_back(fre);
            // 查找元素random不加噪的频率
            double fre_cm = pc_no_noise.Query(dic[random]);
            frequency_cm.push_back(fre_cm);
            // 查找元素random加噪后的频率
            double fre_noise = pc.Query(dic[random]);
            frequency_noise.push_back(fre_noise);
        }else{
            i--;
        }
    }
    vector<double> MAE_Noise_True;
    vector<double> MAE_Cm_Noise;
    vector<double> MAE_Cm_True;
    vector<double> MRE_Noise_True;
    vector<double> MRE_Cm_Noise;
    vector<double> MRE_Cm_True;

    // 计算误差MAE,MRE,ARE,AAE
    double mae_noise_true = MAE_noise_true(frequency,frequency_noise);
    double mae_cm_noise = MAE_cm_noise(frequency_noise,frequency_cm);
    double mae_cm_true = MAE_noise_true(frequency,frequency_cm);
    double mre_noise_true = MRE_noise_true(frequency,frequency_noise);
    double mre_cm_noise = MRE_cm_noise(frequency_noise,frequency_cm);
    double mre_cm_true = MRE_noise_true(frequency,frequency_cm);
    cout<<"MAE between noise and true:"<<mae_noise_true<<endl;
    cout<<"MAE between noise and cm:"<<mae_cm_noise<<endl;
    cout<<"MAE between cm and true:"<<mae_cm_true<<endl;
    cout<<"MRE between noise and true:"<<mre_noise_true<<endl;
    cout<<"MRE between noise and cm:"<<mre_cm_noise<<endl;
    cout<<"MRE between cm and true:"<<mre_cm_true<<endl;
    MAE_Noise_True.push_back(mae_noise_true);
    MAE_Cm_Noise.push_back(mae_cm_noise);
    MAE_Cm_True.push_back(mae_cm_true);
    MRE_Noise_True.push_back(mre_noise_true);
    MRE_Cm_Noise.push_back(mre_cm_noise);
    MRE_Cm_True.push_back(mre_cm_true);


    cout<<endl;
    cout<<"--------------------------------------------------------------------------------"<<endl;
    cout<<"begin test window"<<endl;
    int now = w;
    int query_space = 1000;
    for (int i=w;i<w+(w/sub_num);i=i+query_space) {
        for (int j=0;j<query_space;j++){
            pc.ProcessNew(mini_data[i+j]);
            //cout<<"Process:"<<i+j<<"th item"<<endl;
        }
        now = i+query_space;
        // 添加了1000个元素，查询一次
        vector<int> items;
        vector<int> frequency;
        vector<double> frequency_cm;
        vector<double> frequency_noise;
        cout<<"Begin test:"<<endl;
        for (int k=0;k<100;k++){
            int random = (rand()%(dic.size()-0+1));
            // cout<<random<<endl;
            // 查找元素random的真实频率
            //cout<<"find ture fre"<<now-w<<","<<now<<endl;
            int fre = std::count(mini_data.begin()+now-w, mini_data.begin()+now, dic[random]);
            //cout<<"Bad"<<endl;
            if (fre != 0){
                items.push_back(random);
                frequency.push_back(fre);
                // 查找元素random不加噪的频率
                double fre_cm = pc_no_noise.Query(dic[random]);
                //cout<<"find cm fre"<<endl;
                frequency_cm.push_back(fre_cm);
                // 查找元素random加噪后的频率
                double fre_noise = pc.Query(dic[random]);
                frequency_noise.push_back(fre_noise);
                //cout<<"kth query :"<<k<<endl;
            }else{
                k--;
            }
        }
        cout<<i<<endl;
        double mae_noise_true = MAE_noise_true(frequency,frequency_noise);
        double mae_cm_noise = MAE_cm_noise(frequency_noise,frequency_cm);
        double mae_cm_true = MAE_noise_true(frequency,frequency_cm);
        double mre_noise_true = MRE_noise_true(frequency,frequency_noise);
        double mre_cm_noise = MRE_cm_noise(frequency_noise,frequency_cm);
        double mre_cm_true = MRE_noise_true(frequency,frequency_cm);
        MAE_Noise_True.push_back(mae_noise_true);
        MAE_Cm_Noise.push_back(mae_cm_noise);
        MAE_Cm_True.push_back(mae_cm_true);
        MRE_Noise_True.push_back(mre_noise_true);
        MRE_Cm_Noise.push_back(mre_cm_noise);
        MRE_Cm_True.push_back(mre_cm_true);

    }
    cout<<"MAE between noise and true:"<<endl;
    for(int i=0;i<MAE_Noise_True.size();i++){
        cout<<MAE_Noise_True[i]<<",";
    }
    cout<<endl;
    cout<<"MAE between noise and cm:"<<mae_cm_noise<<endl;
    for(int i=0;i<MAE_Cm_Noise.size();i++){
        cout<<MAE_Cm_Noise[i]<<",";
    }
    cout<<endl;

    cout<<"MAE between cm and true:"<<mae_cm_true<<endl;
    for(int i=0;i<MAE_Cm_True.size();i++){
        cout<<MAE_Cm_True[i]<<",";
    }
    cout<<endl;
    cout<<"MRE between noise and true:"<<mre_noise_true<<endl;
    for(int i=0;i<MRE_Noise_True.size();i++){
        cout<<MRE_Noise_True[i]<<",";
    }
    cout<<endl;
    cout<<"MRE between noise and cm:"<<mre_cm_noise<<endl;
    for(int i=0;i<MRE_Cm_Noise.size();i++){
        cout<<MAE_Cm_Noise[i]<<",";
    }
    cout<<endl;
    cout<<"MRE between cm and true:"<<mre_cm_true;
    for(int i=0;i<MRE_Cm_True.size();i++){
        cout<<MAE_Cm_True[i]<<",";
    }
    cout<<endl;

 */



    return 0;
}
