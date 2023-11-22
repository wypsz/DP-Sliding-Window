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
        std::cerr << "Failed to open the file.1" << std::endl;
        return data;
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        std::istringstream ss(line);
        std::string item;
        std::vector<int> row;

        while (std::getline(ss, item, ',')) {
            try {
                float num = std::stof(item); // 将字符串转换为浮点数
                row.push_back(static_cast<int>(num)); // 将浮点数转换为整数
            } catch (const std::exception& e) {
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
    }
        cout<<"error!"<<endl;
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


std::vector<std::string> loadInitialData(const std::string& filename, int maxLines) {
    std::vector<std::string> data;
    std::ifstream inFile(filename);
    std::string line;

    while (getline(inFile, line) && data.size() < maxLines) {
        data.push_back(line);
    }

    inFile.close();
    return data;
}


std::string readSpecificLine(const std::string& filename, int lineNumber) {
    std::ifstream inFile(filename);
    inFile.seekg(std::ios::beg);

    // 忽略直到指定的行
    for(int i = 0; i < lineNumber - 1; ++i) {
        inFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    std::string line;
    getline(inFile, line);
    inFile.close();
    return line;
}


std::vector<std::string> buildDataSet(const std::string& filename,int maxLines) {
    std::ifstream inFile(filename);
    std::set<std::string> dataSet;
    std::string line;
    int l = 0;
    while (getline(inFile, line) && l < maxLines) {
        dataSet.insert(line);
        l++;
    }

    inFile.close();
    std::vector<std::string> dataVector(dataSet.begin(), dataSet.end());
    return dataVector;
}
int main(int argc, char* argv[]) {



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

//    SmoothHistogram shh(0.24,100000,2000);
//    shh.Add();
//    cout<<shh.getCheckpoints().size()<<endl;
//    return 0;

// 统计所有元素的dic
//    int maxLines = 1200000000;
//    vector<string> dic = buildDataSet(filename, maxLines);
//    std::cout << "Number of unique entries: " << dic.size() << std::endl;


    // 测试两个完整窗口，开始准备参数
//    int w = 1000000;
//    int sub_num = 10;
//    int step = 5000;
//    double all_budget = 1;
//    double rho = 0.6999;
//    double gamma = 0.0001;
//    double beta = 0.005;
//    double alpha = 0.5;
//    double q = 0.3;
    cout<<"aaaaaaaaa"<<endl;
//    if (argc < 16) { // 检查是否有足够的参数
//        std::cerr << "Usage: " << argv[0] << " <w> <sub_num> <step> <all_budget> <rho> <gamma> <beta> <alpha> <q>" << std::endl;
//        return 1;
//    }

    // 将字符串参数转换为整数和浮点数
    int w = std::atoi(argv[1]);
    int sub_num = std::atoi(argv[2]);
    int step = std::atoi(argv[3]);
    double all_budget = std::atof(argv[4]);
    double rho = std::atof(argv[5]);
    double gamma = std::atof(argv[6]);
    double beta = std::atof(argv[7]);
    double alpha = std::atof(argv[8]);
    double q = std::atof(argv[9]);
    int n = std::atoi(argv[10]);
    int query_space = std::atoi(argv[11]);

    string file_path1 = argv[12];
    string file_path2 = argv[13];
    string file_path3 = argv[14];

    // 先读一个窗口的数据
//    vector<string> initialData = loadInitialData(filename, 100000000);
    //ofstream outFile2("C://Users/12293/CLionProjects/DP_Sliding_window/res/network/1.txt");
    //ofstream outFile3("C://Users/12293/CLionProjects/DP_Sliding_window/res/network/1_no.txt");

    ofstream outFile(file_path1);
    ofstream outFile2(file_path2);
    ofstream outFile3(file_path3);
    //    if (!outFile2.is_open()) {
//        std::cerr << "Unable to open file";
//        return 1;
//    }
    ifstream inFile("./dataset/sx-stackoverflow-c2q.txt");
    // 读取txt文件中的数据并对数据进行处理
    if (!inFile) {
        cerr << "failed" << endl;
        return 1;
    }
    string l;
    vector<string> mini_data;
    int wh = 0;
    while (getline(inFile, l)) {
        // 对读取到的每行数据进行处理
        stringstream ls;
        string t;
        ls.str(l);
        ls>>t;
        mini_data.push_back(t);
        wh++;

        if (wh == n){
            break;
        }
    }
    // 关闭txt文件
    inFile.close();
    cout<<mini_data.size()<<endl;
    vector<string> dic;
    set<string> d_all(mini_data.begin(),mini_data.end());
    dic.assign(d_all.begin(),d_all.end());
    cout<<"all dic size:"<<dic.size()<<endl;

    //std::ofstream outFile("C://Users/12293/CLionProjects/DP_Sliding_window/res/network/1.csv");

    if (!outFile.is_open()) {
        std::cerr << "Unable to open file";
        return 1;
    }
    outFile << "w = "<< w<<" , "<< "sub_num = "<< sub_num<<" , "<< "step = "<< step<<" , "<< "all_budget = "<< all_budget<<" , "<< "rho1 = "<< rho<<" , "<< "gamma = "<< gamma<<" , "<< "beta = "<< beta<<" , "<< "alpha = "<<alpha<<" , "<< "q = "<<q<<endl;

//    vector<vector<string>> High_item = loadFromFile("C://Users/12293/CLionProjects/DP_Sliding_window/dataset/High_item.csv");
//    vector<vector<string>> Middle_item = loadFromFile("C://Users/12293/CLionProjects/DP_Sliding_window/dataset/Middle_item.csv");
//    vector<vector<string>> Low_item = loadFromFile("C://Users/12293/CLionProjects/DP_Sliding_window/dataset/Low_item.csv");
//    vector<vector<int>> High_freq = loadFromFile_num("C://Users/12293/CLionProjects/DP_Sliding_window/dataset/High_freq.csv");
//    vector<vector<int>> Middle_freq = loadFromFile_num("C://Users/12293/CLionProjects/DP_Sliding_window/dataset/Middle_freq.csv");
//    vector<vector<int>> Low_freq = loadFromFile_num("C://Users/12293/CLionProjects/DP_Sliding_window/dataset/Low_freq.csv");
//

//    vector<vector<string>> Query_set = loadFromFile("C://Users/12293/CLionProjects/DP_Sliding_window/dataset/worldCup_queries1y.txt");
//    vector<vector<int>> Query_freq = loadFromFile_num("C://Users/12293/CLionProjects/DP_Sliding_window/dataset/worldCup_frequencies1y.txt");


    cout<<"hi"<<endl;
    Private_counter pc (vector<string>(mini_data.begin(),mini_data.begin()+w), w, step, sub_num, rho, gamma, beta, q, alpha,all_budget);
    Private_counter pc_no_noise(vector<string>(mini_data.begin(),mini_data.begin()+w), w, step, sub_num, 0, gamma, beta, q, alpha,0);


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


    cout<<endl;
    cout<<"--------------------------------------------------------------------------------"<<endl;
    cout<<"begin test window"<<endl;
    //vector<vector<int>> query_list;
    //int now = w;
    //int query_space = 20000;
    //int query_space = 50;
    //int query_space = 50000;
    //int query_space = 10000; // 子窗口大小为100w
    outFile2<<"Query_space="<<query_space<<endl;
    outFile3<<"Query_space="<<query_space<<endl;
    int ind = 0; // 标记在indices中的位置
    //int index = 0;
    for (int i=w;i<n;i=i+query_space) {

        // 获取真实frequency

        cout << i << endl;


        // 存每个元素的frequency
        for (int g=0;g<dic.size();g++){
            double fr = pc.Query(dic[g]);
            //outFile2<<"("<<dic[g]<<","<<fr<<");";
            outFile2<<fr<<",";
            double frr = pc_no_noise.Query(dic[g]);
            //outFile3<<"("<<dic[g]<<","<<frr<<");";
            outFile3<<frr<<",";
        }
        outFile2<<endl;
        outFile3<<endl;
        for (int j = 0; j < query_space; j++) {
            //string line = readSpecificLine(filename, i + j + 1);
            pc.ProcessNew(mini_data[i+j]);
            pc_no_noise.ProcessNew(mini_data[i+j]);
        }


    }
    outFile.close();
    outFile2.close();
    outFile3.close();

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
