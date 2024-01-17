//
// Created by 12293 on 2023/12/10.
//


#include <cassert>
#include "PCC_counter.h"

MG::MG(double alpha) {
    this->alpha = alpha;
    this->threshold = ceil((2/ alpha));
}

void MG::processStream(const std::vector<string> &stream) {
    for(int t = 0; t < stream.size(); t++) {
        string sigma = stream[t];
        f[sigma]++;

        // 查找此时f(x)>0有多少个
        if(f.size() > threshold) {
            for(auto it = f.begin(); it != f.end(); ) {
                if(it->second > 0) it->second--;
                if(it->second == 0) {
                    it = f.erase(it);  // 删除计数为0的元素
                } else {
                    ++it;
                }
            }
        }
    }
}

unordered_map<string,double> MG::getResult() {
    return f;
}


unordered_map<string,int> PCC_counter::convertToFrequencyMap(vector<std::string> &vec) {
    std::unordered_map<std::string, int> freqMap;
    for (const auto& str : vec) {
        ++freqMap[str];
    }
    return freqMap;
}

void PCC_counter::initBlocks(const std::vector<std::string>& initialData) {
    int currentLevelBlockSize = S1;

    // 计算有多少叶子节点
//    int leaf_num = W/W0;
//    // 排满子树应该有多少叶子节点
//    int new_num = pow(2, ceil(log2(leaf_num)));
//    // 此时新的W值为
//    int new_w = new_num * W0;
    // 1. 按照初始化好第一个窗口中block的情况构建好跟窗口相关的所有block
    for (int i = 0; i < L; i++) {
        // Calculate the number of blocks for this level
        int numOfBlocks = ceil(W*1.0 / currentLevelBlockSize);
        cout<<i<<" level "<<"numOfBlocks:"<<numOfBlocks<<"block size:"<<currentLevelBlockSize<<endl;


        // Create a vector for the current level
        vector<Block> currentLevel;

        // Initialize blocks for this level
        for (int j = 0; j < numOfBlocks; j++) {
            Block block;
            block.startTime = j * currentLevelBlockSize + 1;
            block.endTime = block.startTime + currentLevelBlockSize - 1;
            block.state = FUTURE;
            currentLevel.push_back(move(block));
        }

        // Add the current level to the blocks
        blocks.push_back(move(currentLevel));

        // Double the block size for the next level
        currentLevelBlockSize *= 2;
    }

    // 2. 遍历blocks，根据提前算好的每个block的startTime和endTime，获取initialData中对应的元素，用PMG处理这部分元素，并更改状态
    // for(auto &level : blocks) {
    for (int i=0 ; i<L; i++){

        //double epsilon_i = epsilon / pow(2,(l-i+1));
        double alpha_i = 1.0 / (pow(2,i) * (L+1));
        //double si = pow(2,i-2) * this->alpha* static_cast<int>(sqrt(W)) /(100 * log2(W));
        for(auto &block : blocks[i]) {
            if(block.endTime <= initialData.size()) {
                // This block is entirely within the initialData
                block.buffer.insert(block.buffer.end(),
                                    initialData.begin() + block.startTime - 1,
                                    initialData.begin() + block.endTime);
                block.pmg = make_unique<MG>(alpha_i);
                block.pmg->processStream(block.buffer);
                block.state = ACTIVE;
                block.fre_true = convertToFrequencyMap(block.buffer);
                block.buffer.clear();  // Clear the buffer after processing
            } else if(block.startTime <= initialData.size()) {
                // This block is partially filled with initialData
                block.buffer.insert(block.buffer.end(),
                                    initialData.begin() + block.startTime - 1,
                                    initialData.end());
                block.pmg = make_unique<MG>(alpha_i);
                block.state = UNDER_CONSTRUCTION;
            }
        }
    }
    currentTime = W;

}



void PCC_counter::processItem(string stream) {
    currentTime++;

    // Block* underConstructionLeaf = nullptr;
    // for (auto &block : blocks[0]) {  // 只检查叶子节点
    //     if (block.state == UNDER_CONSTRUCTION) {
    //         underConstructionLeaf = &block;
    //         break;
    //     }
    // }
    Block* underConstructionLeaf = nullptr;
    if (!blocks[0].empty() && blocks[0].back().state == UNDER_CONSTRUCTION) {
        assert(blocks[0].back().pmg!= nullptr);
        underConstructionLeaf = &blocks[0].back();
    }

    if (!underConstructionLeaf) {
        // 如果没有UNDER_CONSTRUCTION的叶子节点
        // 直接为叶子节点创建一个新的block
        Block newLeaf;
        newLeaf.state = UNDER_CONSTRUCTION;
        newLeaf.startTime = currentTime;
        newLeaf.endTime = currentTime+S1-1;
        newLeaf.buffer.push_back(stream);

        if (newLeaf.endTime - newLeaf.startTime + 1 == 1) {  // 如果叶子节点只需要一个元素
            cout<<"1 item"<<endl;
            newLeaf.state = ACTIVE;
            newLeaf.pmg = make_unique<MG>(alpha);
            newLeaf.pmg->processStream(newLeaf.buffer);
            newLeaf.fre_true = convertToFrequencyMap(newLeaf.buffer);
            newLeaf.buffer.clear();
        } else {
            if (!newLeaf.pmg) {
//                double epsilon_i = epsilon / pow(2,(l-0+1));
                double alpha_i = 1.0 / (pow(2,0) * (L+1));
                //double si = pow(2,1-2) * this->alpha* static_cast<int>(sqrt(W)) /(100 * log2(W));
                newLeaf.pmg = make_unique<MG>(alpha_i);
            }
        }

        blocks[0].push_back(move(newLeaf));

        // 对于更高的层次，根据最后一个block的状态进行操作
        for (int i = 1; i < L; i++) {

//            double epsilon_i = epsilon / pow(2,(l-i+1));
            double alpha_i = 1.0 / (pow(2,i) * (L+1));
            //double si = pow(2,i-2) * this->alpha* static_cast<int>(sqrt(W)) /(100 * log2(W));


            Block& lastBlock = blocks[i].back();
            if (lastBlock.state == UNDER_CONSTRUCTION) {
                lastBlock.buffer.push_back(stream);
                if (lastBlock.buffer.size() == S1*pow(2, i)) {
                    lastBlock.state = ACTIVE;
                    if (!lastBlock.pmg){
                        lastBlock.pmg = make_unique<MG>(alpha_i);
                        assert(lastBlock.pmg!= nullptr);
                    }
                    lastBlock.pmg->processStream(lastBlock.buffer);
                    lastBlock.fre_true = convertToFrequencyMap(lastBlock.buffer);

                    lastBlock.buffer.clear();
                }
            } else if (lastBlock.state == ACTIVE) {
                Block newBlock;
                newBlock.state = UNDER_CONSTRUCTION;
                newBlock.startTime = currentTime;
                newBlock.endTime = currentTime + S1*pow(2, i) - 1;
                newBlock.buffer.push_back(stream);
                if (!newBlock.pmg) {
                    newBlock.pmg = make_unique<MG>(alpha_i);
                    assert(newBlock.pmg!= nullptr);
                }
                blocks[i].push_back(move(newBlock));
            }
        }
    } else {
        underConstructionLeaf->buffer.push_back(stream);
        if (underConstructionLeaf->buffer.size() == S1*pow(2, 0)) {
            underConstructionLeaf->state = ACTIVE;
            //assert(underConstructionLeaf->pmg != nullptr);
            underConstructionLeaf->pmg->processStream(underConstructionLeaf->buffer);
            underConstructionLeaf->fre_true = convertToFrequencyMap(underConstructionLeaf->buffer);

            underConstructionLeaf->buffer.clear();

            // 删除过期的
        }
        for (int i = 1; i < L; i++) {

//            double epsilon_i = epsilon / pow(2,(l-i+1));
            double alpha_i = 1.0 / (pow(2,i) * (L+1));
            //double si = pow(2,1-2) * this->alpha* static_cast<int>(sqrt(W)) /(100 * log2(W));


            Block& lastBlock = blocks[i].back();
            if (lastBlock.state == UNDER_CONSTRUCTION) {
                lastBlock.buffer.push_back(stream);
                if (lastBlock.buffer.size() == S1*pow(2, i)) {
                    lastBlock.state = ACTIVE;
                    if (!lastBlock.pmg){
                        lastBlock.pmg = make_unique<MG>(alpha_i);
                        assert(lastBlock.pmg!= nullptr);
                    }
                    lastBlock.pmg->processStream(lastBlock.buffer);
                    lastBlock.fre_true = convertToFrequencyMap(lastBlock.buffer);

                    lastBlock.buffer.clear();
                }
            } else if (lastBlock.state == ACTIVE) {
                Block newBlock;
                newBlock.state = UNDER_CONSTRUCTION;
                newBlock.startTime = currentTime;
                newBlock.endTime = currentTime + S1*pow(2, i) - 1;
                newBlock.buffer.push_back(stream);
                if (!newBlock.pmg) {
                    newBlock.pmg = make_unique<MG>(alpha_i);
                    assert(newBlock.pmg!= nullptr);
                }
                blocks[i].push_back(move(newBlock));
            }
        }

    }
    if (currentTime % S1 == 0){
        deleteExpiredBlocks();
    }

}

std::unordered_map<string, double> PCC_counter::query_all() {
    std::vector<std::pair<int, int>> uncoveredTimeRanges;
    // 这里要改！！！！！！！ 应该是离当前时间最近的上一个，
    int mod = currentTime % S1;
    cout<<"current:"<<currentTime-mod<<endl;
    uncoveredTimeRanges.push_back({currentTime -mod - W + 1, currentTime - mod});

    // 结果map
    std::unordered_map<string, double> resultMap;

    for (int i = blocks.size()-1 ; i >= 0 && !uncoveredTimeRanges.empty(); --i) {
        auto& currentLayer = blocks[i];
        std::vector<std::pair<int, int>> newUncoveredTimeRanges;

        for (Block & block : currentLayer) {
            std::vector<std::pair<int, int>> tempUncovered;  //int windowEndTime = it->second;

                for (auto& range : uncoveredTimeRanges ){
                    if (block.endTime<=range.second && block.startTime >= range.first && block.state == ACTIVE){
                        auto blockResult = block.pmg->getResult();
                        //cout<<"size:"<<blockResult.size()<<endl;
                        //cout<<"result size1:"<<resultMap.size()<<endl;
                        // 获得lambda_i,求epsilon_i
                        double alpha_i = block.pmg->alpha;
                        double epsilon_i = (2/alpha_i) / (epsilon);
                        //double epsilon_i = epsilon/L
                        cout<<"alpha_i:"<<alpha_i<<"epsilon_i:"<<epsilon_i<<endl;
                        int x = 0;
                        // 获取buffer
                    
                        for (const auto& pair : blockResult) {
                            x++;
                            //double noise = generateLaplacianRandom(alpha*(static_cast<int>(sqrt(W)) /(16 * log2(W) *log2(W))),1)[0];
                            // 通过fre_true获得元素中真实频率
                            int fre = block.fre_true[pair.first];
                            //cout<<"true fre:"<<fre<<endl;
                            generator.seed(currentSeed);
                            double noise = generateLaplacianRandom(generator,epsilon_i,1)[0];
                            currentSeed++;
                            //cout<<"noise:"<<noise<<endl;
                            //resultMap[pair.first] = resultMap[pair.first]+fre;
                            resultMap[pair.first] = resultMap[pair.first]+fre+noise;
                            //cout<<"resultMap fre:"<<resultMap[pair.first]<<endl;
                        }
                        if (block.startTime > range.first){
                            tempUncovered.push_back({range.first,block.startTime-1});
                        }
                        if (block.endTime < range.second){
                            tempUncovered.push_back({block.endTime+1,range.second});
                        }
                    } else{
                        tempUncovered.push_back(range);
                    }
                    }
                uncoveredTimeRanges = tempUncovered;

                }
        }

        //
        //uncoveredTimeRanges.insert(uncoveredTimeRanges.end(), newUncoveredTimeRanges.begin(), newUncoveredTimeRanges.end());
    // 遍历resultMap,再次筛选出heavy hitter
    cout<<"original size:"<<resultMap.size()<<endl;
    cout<<"true_the:"<<true_the<<endl;
    //double the = alpha * W;
    auto temp = resultMap["209"];
    for (auto it = resultMap.begin(); it != resultMap.end();) {
        if (it->second < true_the) {
            // 删除这组key-value对，并更新迭代器
            it = resultMap.erase(it);
        } else {
            // 否则，只移动迭代器
            ++it;
        }
    }


    cout<<"all resultMap size:"<<resultMap.size()<<endl;
    return resultMap;
}



vector<double> PCC_counter::generateLaplacianRandom(std::default_random_engine& generator, double scale, int numSamples) {
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

void PCC_counter::deleteExpiredBlocks() {
    for (int i=0;i<L;i++){
        for (auto it = blocks[i].begin(); it != blocks[i].end();) {
            if (it->endTime < currentTime - W+1) {  // 当前时间超过窗口W的结束时间，表示block已经过期
                it = blocks[i].erase(it);
            } else {
                ++it;
            }
        }
    }
}

PCC_counter::PCC_counter(double epsilon, double alpha, int W,double true_the, const vector<std::string> &initialData) {
    this->alpha = alpha;
    this->W =W;
    this->currentTime = 0;
    this->epsilon = epsilon;
    this->true_the = true_the*W;
    this->initialSeed = 0;
    this->currentSeed = initialSeed;
    //this-> L = ceil(log2(100*epsilon* static_cast<int>(sqrt(W))/alpha)  + log2(log2(W))) +2;
    this->L=log2(4/ alpha);
    //this-> S1 = pow(2,(1-2)) * alpha * static_cast<int>(sqrt(W)) /(100 * log (W));
    this->S1 = ceil(alpha * W / 4);

    initBlocks(initialData);
}
