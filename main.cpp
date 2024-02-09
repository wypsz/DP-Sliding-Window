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
#include "Mechanism/PCC_counter.h"
#include "Mechanism/PCC.h"
#include <chrono>
#include <iomanip>
#include <ctime>
#include "Mechanism/Upadhyay.h"
#include "Mechanism/DPSW.h"
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

double calculateLeftSide(double rho, double delta) {
    return rho + 2 * sqrt(rho * log(1.0 / delta));
}


double findRho(double epsilon, double delta, double tol = 1e-9) {
    double low = 0, high = epsilon;
    while (high - low > tol) {
        double mid = low + (high - low) / 2;
        double leftSide = calculateLeftSide(mid, delta);

        if (leftSide > epsilon) {
            high = mid;
        } else {
            // 当左侧结果等于或略低于epsilon时，适当增加low
            low = mid;
        }
    }
    // 检查并确保不超过epsilon
    if (calculateLeftSide((low + high) / 2, delta) > epsilon) {
        return low;
    }
    return (low + high) / 2;
}




int main(int argc, char* argv[]) {


    std::vector<int> differences = {
            22733, 28067, 12465, 29859, 11736, 29160, 7014, 33083, 9352,
            31182, 6055, 38278, 1104, 38519, 2243, 36072, 5972, 23323,
            20705, 19115, 16127, 22909, 19779, 23026, 14169, 19825, 25728,
            14897, 32690, 4913, 30101, 10555, 20042, 36246, 3380, 35898,
            7016, 25399, 17951, 14027, 37962, 1475, 31572, 9484, 29161,
            16853, 18800, 17601, 19938, 29116, 8965, 32113, 11230, 19369,
            29280, 10438, 28130, 9904, 32522, 9315, 27441, 15867, 16937,
            30021, 8238, 27335, 13175, 33342, 3786, 38551, 1689, 32738,
            8009, 33253, 10927, 28009, 7944, 33108, 9606, 25441, 12325,
            29807, 15607, 23861, 11998, 21512, 29696, 6395, 29237, 15998,
            15677, 25893, 24505, 8386, 23506, 29842, 8350, 34267, 6214,
            17564
    };


    // run DPSW
    int w = std::atoi(argv[1]);
    int sub_num = std::atoi(argv[2]);
    int step = std::atoi(argv[3]);
    double epsilon = std::atof(argv[4]);
    double q = std::atof(argv[5]);
    double gamma = std::atof(argv[6]);
    double beta = std::atof(argv[7]);
    double alpha = std::atof(argv[8]);
    int n = std::atoi(argv[9]);
    int hseed = std::atoi(argv[10]);

    string file_path1 = argv[11];
    string file_path2 = argv[12];
    string infile_path = argv[13];
    int datasize = std::atoi(argv[14]);


    ofstream outFile(file_path1);
    ofstream outFile2(file_path2);
    ifstream inFile(infile_path);
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

    if (!outFile2.is_open()) {
        std::cerr << "Unable to open file";
        return 1;
    }
    if (!outFile.is_open()) {
        std::cerr << "Unable to open file";
        return 1;
    }

    double delta = 1.0 / pow(datasize,1.5);
    double rho = findRho(epsilon,delta);



    cout<< "w = "<< w<<" , "<< "sub_num = "<< sub_num<<" , "<< "step = "<< step<<" , "<< "epsilon = "<<epsilon<<" , "<<"rho = "<<rho<<", n = "<<n<<" , "<< "gamma = "<< gamma<<" , "<< "beta = "<< beta<<" , "<< "alpha = "<<alpha<<" , "<< "q = "<<q<<" datasize:"<<datasize<<", dicsize:"<<dic.size()<<endl;

    outFile<< "w = "<< w<<" , "<< "sub_num = "<< sub_num<<" , "<< "step = "<< step<<" , "<< "epsilon = "<<epsilon<<" , "<<"rho = "<<rho<<", n = "<<n<<" , "<< "gamma = "<< gamma<<" , "<< "beta = "<< beta<<" , "<< "alpha = "<<alpha<<" , "<< "q = "<<q<<" datasize:"<<datasize<<", dicsize:"<<dic.size()<<endl;

    if (w != 1000000){
        for (int i=0;i<differences.size();i++){
            double diff = (w*1.0)/(1000000.0)*differences[i];
            differences[i] = floor(diff);
            cout<<differences[i]<<" ";
        }
    }

    cout<<"hi"<<endl;
    auto start = std::chrono::high_resolution_clock::now();
    DPSW pc (vector<string>(mini_data.begin(),mini_data.begin()+w), w, step, sub_num, rho, gamma, beta, q, alpha,hseed,hseed*5);

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    outFile<<"initialize time:"<<elapsed.count()<<endl;

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

    cout<<endl;
    cout<<"--------------------------------------------------------------------------------"<<endl;
    cout<<"begin test window"<<endl;

    std::chrono::duration<double> total_elapsed(0);
    for (int g=0;g<dic.size();g++){
        start = std::chrono::high_resolution_clock::now();
        double fr = pc.Query(dic[g]);
        end = std::chrono::high_resolution_clock::now();
        total_elapsed += end - start;
        outFile2<<fr<<",";
    }

    outFile<<"query "<<dic.size()<<":"<<total_elapsed.count()<<endl;
    outFile2<<endl;
    int ind = w;
    for (int i=0;i<differences.size();i++){
        int process_num = differences[i];
        start = std::chrono::high_resolution_clock::now();
        for (int j = 1; j <= process_num; j++) {
            pc.ProcessNew(mini_data[ind+j]);
        }
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        outFile<<"process 20000:"<<elapsed.count()<<endl;
        ind += process_num;
        cout<<"ind: "<<ind<<endl;
        std::chrono::duration<double> total_elapsed(0);
        for (int g=0;g<dic.size();g++){
            start = std::chrono::high_resolution_clock::now();
            double fr = pc.Query(dic[g]);
            end = std::chrono::high_resolution_clock::now();
            total_elapsed += end - start;
            outFile2<<fr<<",";
        }
        outFile<<"query "<<dic.size()<<":"<<total_elapsed.count()<<endl;
        outFile2<<endl;
    }
    outFile.close();
    outFile2.close();

    /*
    // run U-Sketch
    int w = std::atoi(argv[1]);
    int sub_num = std::atoi(argv[2]);
    int step = std::atoi(argv[3]);
    double epsilon = std::atof(argv[4]);
    double gamma = std::atof(argv[5]);
    double beta = std::atof(argv[6]);
    double zeta = std::atof(argv[7]);
    int n = std::atoi(argv[8]);
    string file_path1 = argv[9];
    string time_path = argv[10];
    string in_path = argv[11];


    ofstream time_file(time_path);
    ofstream outFile_fre(file_path1);



    ifstream inFile(in_path);
    // 读取txt文件中的数据并对数据进行处理
    if (!inFile) {
        cerr << "failed" << endl;
        return 1;
    }
    if (!outFile_fre) {
        cerr << "failed" << endl;
        return 1;
    }
    if (!time_file) {
        cerr << "failed" << endl;
        return 1;
    }
    string l;
    vector<string> mini_data;
    int wh = 0;
    while (getline(inFile, l)) {
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
    time_file<< "w = "<< w<<" , "<< "sub_num = "<< sub_num<<" , "<< "step = "<< step<<" , "<< "epsilon = "<<epsilon<<" , "<<"n = "<<n<<" , "<< "gamma = "<< gamma<<" , "<< "beta = "<< beta<<" , zeta"<<zeta<<"datasize:"<<mini_data.size()<<", dicsize:"<<dic.size()<<endl;
    int flag = 1;
    if (w % sub_num == 0){
        flag = 0;
    }
    cout<<flag<<endl;

    if (w != 1000000){
        for (int i=0;i<differences.size();i++){
            double diff = (w*1.0)/(1000000.0)*differences[i];
            differences[i] = floor(diff);
            cout<<differences[i]<<" ";
        }
    }


    cout<<"hi"<<endl;

    auto start = std::chrono::high_resolution_clock::now();
    Private_Heavy ph(vector<string>(mini_data.begin(),mini_data.begin()+w),w,sub_num,step,epsilon,gamma,beta,zeta,dic,flag,0);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    time_file<<"initialize time:"<<elapsed.count()<<endl;
    cout<<"Private HH and CS first window initialize finish"<<endl;

    SmoothHistogram sh(beta,w/(sub_num/10),step);
    sh.Add();
    vector<int> indices = sh.getIndices();
    vector<int>checkpoints = sh.getCheckpoints();
    time_file<<"indices:"<<endl;
    WriteVectorToCSVint(time_file,indices);
    time_file<<"checkpoints:"<<endl;
    WriteVectorToCSVint(time_file,checkpoints);
    vector<vector<double>> parameter = ph.show_parameter();
    time_file<<"Parameters in sub_window"<<endl;
    for (int a=0;a<parameter.size();a++){
        WriteVectorToCSV(time_file,parameter[a]);
    }
    //outFile_fre<<"Parameters in sub_window"<<endl;
    cout<<"---------------------------------------------------------------------------------"<<endl;


    cout<<endl;
    cout<<"--------------------------------------------------------------------------------"<<endl;
    cout<<"begin test window"<<endl;

    start = std::chrono::high_resolution_clock::now();
    vector<double> res = ph.Query();
    cout<<"error"<<endl;
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    time_file<<"query "<<dic.size()<<":"<<elapsed.count()<<endl;
    for (int z= 0;z<dic.size();z++){
        outFile_fre<<res[z]<<',';
    }
    outFile_fre<<endl;


    int ind = w;
    for (int i=0;i<differences.size();i++) {
        int process_num = differences[i];
        start = std::chrono::high_resolution_clock::now();
        for (int j = 1; j <= process_num; j++){
            ph.ProcessNew(mini_data[ind+j]);
        }
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        time_file<<"process 20000:"<<elapsed.count()<<endl;
        ind += process_num;
        cout<<"ind: "<<ind<<endl;
        start = std::chrono::high_resolution_clock::now();
        vector<double> res = ph.Query();
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        time_file<<"query "<<dic.size()<<":"<<elapsed.count()<<endl;
        for (int z= 0;z<dic.size();z++){
            outFile_fre<<res[z]<<',';
        }
        outFile_fre<<endl;
    }
    time_file.close();
    outFile_fre.close();
*/


/*
    // run PCC

    int w = std::atoi(argv[1]);
    double alpha = std::atof(argv[2]);
    double all_budget = std::atof(argv[3]);
    int n = std::atoi(argv[4]);
    double true_the = std::atof(argv[5]);
    string file_path1 = argv[6];
    string time_path = argv[7];
    string infile_path = argv[8];

    ofstream outFile(file_path1);
    ifstream inFile(infile_path);
    ofstream timeFile(time_path);
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

    inFile.close();
    cout<<"mini_data size:"<<mini_data.size()<<endl;
    vector<string> dic_part;
    set<string> d_part(mini_data.begin(),mini_data.end());
    dic_part.assign(d_part.begin(),d_part.end());
    cout<<"part dic size:"<<dic_part.size()<<endl;
    if (!outFile.is_open()) {
        std::cerr << "Unable to open file";
        return 1;
    }

    auto start = std::chrono::high_resolution_clock::now();
    PCC pcc_counter(all_budget,alpha,true_the,w,vector<string> (mini_data.begin(),mini_data.begin()+w));
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    timeFile<<"w="<<w<<" , dicsize="<<dic_part.size()<<" , datasize:"<<mini_data.size()<<endl;
    timeFile<<"initialize time:"<<elapsed.count()<<endl;

    cout<<"first window initialize finish"<<endl;
    if (w!=1000000){
        for(int i=0;i<differences.size();i++){
            differences[i] = floor(differences[i] *(w/1000000.0));
        }
    }


    unordered_map<std::string, double> res= pcc_counter.query_all();
    cout<<"1:"<<res.size()<<endl;
    for (const auto &pair :res){
        outFile<<pair.first<<":"<<pair.second<<",";
    }
    outFile<<endl;

    int ind = w;
    for (int i=0;i<differences.size();i++) {
        cout<<ind + differences[i]<<endl;
        int process_num = differences[i];
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        start = std::chrono::high_resolution_clock::now();
        for (int j = 1; j <= process_num; j++) {
            pcc_counter.processItem(mini_data[ind+j]);
        }
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        timeFile<<"process 20000:"<<elapsed.count()<<endl;

        start = std::chrono::high_resolution_clock::now();
        unordered_map<std::string, double> res= pcc_counter.query_all();
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        timeFile<<"query time:"<<elapsed.count()<<endl;
        start = std::chrono::high_resolution_clock::now();
        for (const auto &pair :res){
            outFile<<pair.first<<":"<<pair.second<<",";
        }
        outFile<<endl;
    }
    outFile.close();
    timeFile.close();

    */




/*
// run BLMZ-Sketch

// 将字符串参数转换为整数和浮点数
    int w = std::atoi(argv[1]);
    int level = std::atoi(argv[2]);
    double alpha = std::atof(argv[3]);
    double all_budget = std::atof(argv[4]);
    int n = std::atoi(argv[5]);
    string file_path1 = argv[6];
    string time_path = argv[7];
    string infile_path = argv[8];


    ofstream outFile(file_path1);
    ifstream inFile(infile_path);
    ofstream timeFile(time_path);
    if (!inFile) {
        cerr << "failed" << endl;
        return 1;
    }
    string l;
    vector<string> mini_data;
    int wh = 0;


    while (getline(inFile, l)) {
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

    inFile.close();
    cout<<"mini_data size:"<<mini_data.size()<<endl;
    vector<string> dic_part;
    set<string> d_part(mini_data.begin(),mini_data.end());
    dic_part.assign(d_part.begin(),d_part.end());
    cout<<"part dic size:"<<dic_part.size()<<endl;

    if (!outFile.is_open()) {
        std::cerr << "Unable to open file";
        return 1;
    }

    auto start = std::chrono::high_resolution_clock::now();
    PCC_counter pcc_counter(all_budget,level,alpha,w,vector<string> (mini_data.begin(),mini_data.begin()+w));
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    timeFile<<"w="<<w<<" , dicsize="<<dic_part.size()<<" , datasize:"<<mini_data.size()<<endl;
    timeFile<<"initialize time:"<<elapsed.count()<<endl;

    cout<<"first window initialize finish"<<endl;

    unordered_map<std::string, double> res= pcc_counter.query_all();
    for (const auto &pair :res){
        outFile<<pair.first<<":"<<pair.second<<",";
    }
    outFile<<endl;
    if (w!=1000000){
        for(int i=0;i<differences.size();i++){
            differences[i] = floor(differences[i] *(w/1000000.0));
        }
    }

    int ind = w;
    for (int i=0;i<differences.size();i++) {
        int process_num = differences[i];
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        start = std::chrono::high_resolution_clock::now();
        for (int j = 1; j <= process_num; j++) {
            pcc_counter.processItem(mini_data[ind+j]);
        }
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        timeFile<<"process 20000:"<<elapsed.count()<<endl;


        start = std::chrono::high_resolution_clock::now();
        unordered_map<std::string, double> res= pcc_counter.query_all();
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        timeFile<<"query time:"<<elapsed.count()<<endl;
        start = std::chrono::high_resolution_clock::now();
        for (const auto &pair :res){
            outFile<<pair.first<<":"<<pair.second<<",";
        }
        outFile<<endl;
    }
    outFile.close();
    timeFile.close();

    */
    return 0;
}
