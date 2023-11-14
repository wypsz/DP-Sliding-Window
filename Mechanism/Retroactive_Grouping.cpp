//
// Created by 12293 on 2023/9/7.
//

#include "Retroactive_Grouping.h"

/*
 * 考虑发布直方图，直方图中的每个bin对应一个属性值，计数是窗口中落入bin的总数
 *
 * 在每个时间戳使用完整的隐私预算
 *
 * 在第1行中，监控模块决定适当的采样率并计算差值向量d，该差值向量d以平方误差的形式当前时间单位和前一时刻之间所有bin的计数差值。
 * 由于监控模块对原始数据进行操作，因此需要通过差分隐私进行保护。
 * 为此，为其分配隐私预算ε1。我们没有比较整个直方图的 SSE，而是观察到不同箱的计数可能以不同的速率变化，因此单独处理每个箱。
 * d中的每个元素记录了 bin 计数的差异。
 */

Retroactive_Grouping::Retroactive_Grouping(const vector<string>& xw, vector<int> indices, double epsilon1, double epsilon2,vector<string>& dic, int w, int sub_num)
:noise_histogram(dic.size(),0.0),histogram(dic.size(),0),difference(dic.size(),0),G(dic.size()),Lamb(dic.size(),0){
    // 初始化一些变量
    this->epsilon1 = epsilon1;
    this->epsilon2 = epsilon2;
    this->indices = indices;
    this->dic  = dic;
    this->w = w;
    this->sub_num = sub_num;
    /* xw是第一个窗口的元素, indices是Private counter中每个子窗口内checkpoint位置进行处理，在这个算法中每个时刻处理多少条数据，sub_num是子窗口个数
     *
     *  即在这个算法中，有|indices| * sub_sum个时刻！
     * 按照同样也是就近查询，会有一定的偏差，但同样是把每个时刻的histogram结果加起来
     *
     */
    // 开始处理第一个窗口window
    int sub_size = w/sub_num; // 每个子窗口的大小
    this->be = 0;
    index.push_back(be);
    for (int i=0;i<sub_num;i++){
//    while(be<xw.size()){
        for (int j=1;j<indices.size();j++){
            // 每个时刻处理[be , be + indices[j]]的元素
            int a = (i*sub_size)+indices[j-1];
            int b = (i*sub_size)+indices[j];
            vector<string> tmp(xw.begin()+a,xw.begin()+b);
            // 1. 先计算好tmp中元素的频率，即这个时刻的真实frequency
            //cout<<indices[j-1]<<" "<<indices[j]<<" ";
            for (int k = 0;k<tmp.size();k++){
                // 找到tmp[k]在dic中对应的索引，并更新histogram中的计数值
                auto it = find(dic.begin(),dic.end(),tmp[k]);
                if (it != dic.end()){
                    // 元素找到
                    int ind = distance(dic.begin(),it);
                    histogram[ind] ++;
                }
            }
            H.push_back(histogram);

            // 初始化第一个时刻，发布所有item的noise版本
            if (i == 0 && j == 1){
                for ( int k=0;k<dic.size();k++){
                    vector<double> Lap = generateLaplacianRandom(1/epsilon2,1);
                    noise_histogram[k] = histogram[k] + Lap[0];
                    Lamb[k] = Lap[0];
                }
                H_noise.push_back(noise_histogram);
                for (int k=0;k<dic.size();k++){
                    noise_histogram[k] = 0.0;
                }
            }else{
                monitor(tmp);// 更新了difference
                // 开始构建publication
                publication();
            }
            // 全部重置为0
            for (int k = 0;k<histogram.size();k++){
                histogram[k] = 0;
            }

            // histogram和noise histogram都已经计算完毕，并push入H和H_noise中，可以将其设置为0，等待下一批数据来计算。
            be = i*sub_size+indices[j];
            index.push_back(be);
            cout<<"index:"<<be<<endl;
            cout<<"i,j,indices_size:"<<i<<","<<j<<","<<indices.size()<<endl;
        }
    }
}

void Retroactive_Grouping::insert(string& newItem) {
    // 到达了一个元素，先判断这个元素indices是否该更新；
    last_update.push_back(newItem);
    if (last_update.size() == indices[last_indices]-indices[last_indices-1]){// last_indices从1开始
        // ind是这个元素在stream中的位置index，需要先确定好这个元素在哪个子窗口，以及在子窗口的哪个checkpoint中
        for (int i = 0;i<dic.size();i++){
            histogram[i] = 0;
        }
        // 1. 统计好这个时刻真实的frequency
        for (int i=0;i<last_update.size();i++){
            // 找到newItem[i]在dic中对应的索引，并更新histogram中的计数值
            auto it = find(dic.begin(),dic.end(),last_update[i]);
            if(it != dic.end()){
                // 元素找到
                int x = distance(dic.begin(),it);
                histogram[x] ++;
            }
        }
        be += last_update.size();
        index.push_back(be);
        H.push_back(histogram);
        // 更新difference
        monitor(last_update);
        // 开始构建publication
        publication();
        // histogram和noise histogram都已经计算完毕，并push入H和H_noise中，可以将其设置为0，等待下一批数据来计算。
        // 适当删除过期的时刻,H和H_noise

        H_noise.erase(H_noise.begin());
        H.erase(H.begin());

        // 适当更新last_indices和last_update
        // last_indices +1 或置为0
        last_indices ++ ;
        if (last_indices > indices.size()){
            last_indices = 1;
        }
        last_update.clear();

    }

}


double Retroactive_Grouping::query(string& item) {
    /*
    // 查询窗口大小为w内的元素频率
    //q_time是最后的时刻，从后往前找、
    // 从前往后，找到第一个大于 q_begin 的索引
    int q_begin = be+1 - w;
    int pos=0;
    for (int i = 0;i<index.size();i++){
        if (index[i]>=q_begin){
            if (i == 0 && index[0]>q_begin){
                cout<<"error"<<endl;
            }else{
                pos = i;
                break;
            }
            if (q_begin-index[i-1] < index[i]-q_begin){
                pos = i-1;
                break;
            }else{
                pos = i;
                break;
            }
        }
    }
    // 找到窗口最接近的pos，开始组装结果
    vector<double> res(dic.size(),0);
    for(int i=pos;i<index.size()-1;i++){
        for (int j=0;j<dic.size();j++){
            res[j]+=H_noise[i][j];
        }
    }
     */
    // cout<<H_noise.size()<<endl;
    // 找到item在dic中的索引
    auto it = find(dic.begin(),dic.end(),item);


    int pos = distance(dic.begin(),it);
    double res = 0.0;
    for (int i=0;i<H_noise.size();i++){
        res += H_noise[i][pos];
    }
    return res;

}

vector<double> Retroactive_Grouping::generateLaplacianRandom(double scale,int numSamples) {
    std::random_device rd;
    std::default_random_engine generator(rd());
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


bool Retroactive_Grouping::bernoulliSample(double probability) {
    random_device rd;
    default_random_engine generator(rd());
    bernoulli_distribution distribution(probability);

    return distribution(generator);
}

/*
 * 在每个时刻ti，对于每一个bin Bj，估计当前时间单元[t(i-1),t(i)]和前一个时间单元[t(i-2),t(i-1)]中的计数差异。
 * 1. 确定Bj的适当采样率rj，(使用伯努利抽样，对于每个即将发生的时间，以成功概率0<= rj<=1，成功iu作为样本，失败就丢弃这个事件)
 *      总体大小Nj，抽样的样本大小nj，
 */

void Retroactive_Grouping::monitor(vector<std::string>& NewItems) {
    // 计算当前时间单位的计数与前一个时间单位计数之间的差异
    // 此时的noise_histogram和dic是上一个时间段[t(i-2),t(i-1)]的情况
    for(int i=0;i<dic.size();i++){
        // 对于元素Bi 查看情况更新差异d[i]

        // 1. 计算sample rate rj
        float sample_rate = 0.1*epsilon1;
        // 2. 统计出NewItem中的Bi的数量
        int c = count(NewItems.begin(),NewItems.end(),dic[i]);
        // 3. 根据sample rate进行抽样（伯努利）
        int Si_count = 0;
        for (int j = 0;j<c;j++){
            // 伯努利抽样成功的概率为rj
            bool res = bernoulliSample(sample_rate);
            if (res){
                // 伯努利采样为true，采样这个事件
                Si_count++;
            }

        }
        // 4. 开始计算N_j,加噪估计出来的Bi的frequency
        vector<double> Lap = generateLaplacianRandom(1/(log(exp(epsilon1)-1+sample_rate)- log(sample_rate)),1);
        double N_i = (Si_count + Lap[0])/sample_rate;
        // 通过访问H_noise中的最后一个元素，即H[t(i-2),t(i-1)]
        if (H_noise.empty()){
            difference[i] = pow(0-N_i,2);
        }else{
            difference[i] = pow((H_noise[H_noise.size()-1][i] - N_i),2);

        }
    }

}

void Retroactive_Grouping::publication() {
    // 当前时刻应该如何发布？
    for (int i=0;i<dic.size();i++){
        double threshold = 2.0/pow(epsilon2,2);
        if (difference[i] >= threshold){
            // 直接发布真实fre+noise
            vector<double> Lap = generateLaplacianRandom(1/epsilon2,1);
            noise_histogram[i] = histogram[i] + Lap[0];
            Lamb[i] = 1/epsilon2;
            if(!G[i].empty()){
                G[i].clear();
            }
        }else{
            // 将当前时刻的真实fre添加到G中
            G[i].push_back(histogram[i]);
            if(Lamb[i] == 0){
                cout<<"lambda == 0"<<endl;
            }
            m = 2/(Lamb[i] * epsilon2);
            if (G[i].size() >= m){
                vector<double> Lap = generateLaplacianRandom(1/(G[i].size() * epsilon2),1);
                double tmp = 0.0;
                for (int j = 0;j<G[i].size();j++){
                    tmp += G[i][j];
                }
                noise_histogram[i] = tmp + Lap[0];
                Lamb[i] = 1/(G[i].size() * epsilon2);
                G[i].clear();
            } else{
                // 被动发布
                if (!H_noise.empty()){
                    noise_histogram[i] = H_noise[H_noise.size()-1][i];
                }else{
                    cout<<"H_noise is empty"<<endl;
                    noise_histogram[i] = 0.0;
                }
            }
        }
    }
    // 将整理好的当前时刻的counter情况push入H中方便后续汇总计算
    // cout<<"finish publicate"<<endl;
    H_noise.push_back(noise_histogram);
    for(int k=0;k<noise_histogram.size();k++){
        noise_histogram[k] = 0.0;
    }



}