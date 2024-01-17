//
// Created by 12293 on 2023/10/12.
//

#include "PCC.h"
#include <climits>
#include <cassert>


PrivateMisraGries::PrivateMisraGries(double epsilon, double lambda) {
    this->epsilon = epsilon;
    this->lambda = lambda;
    this->threshold = ceil(2/lambda);
}

void PrivateMisraGries::processStream(const std::vector<string>& stream) {
    gen.seed(seed);
    for(int t = 0; t < stream.size(); t++) {
        string sigma = stream[t];
        f[sigma]++;

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

    int sensitivity = threshold + 1;



    for(auto& [key, value] : f) {
        int r = generateSymmetricGeometricRandomNumber(exp(epsilon / (sensitivity + 1)));
        if (r!=0){
            cout<<r<<endl;
        }
        f[key] = std::max(f[key] + r, 0.0);
        seed++;
        //cout<<f[key]<<endl;
    }

    if(f.size() > threshold) {
        std::vector<std::pair<string, int>> counts;

        for (const auto& kv : f) {
            counts.push_back(kv);
        }

        std::sort(counts.begin(), counts.end(), [](const auto& a, const auto& b) {
            return a.second < b.second;
        });

        while(f.size() > threshold) {
            int minValue = counts.front().second;
            auto it = counts.begin();
            while (it != counts.end() && it->second == minValue) {
                f.erase(it->first);
                it = counts.erase(it);
            }
        }
    }
}

std::unordered_map<string, double> PrivateMisraGries::getResult() {
    return f;
}

// int PrivateMisraGries::generateSymmetricGeometricRandomNumber(double alpha){
//     static random_device rd;
//     static mt19937 gen(rd());

//     geometric_distribution<> d((alpha-1) / alpha);

//     uniform_int_distribution<> dis(0,1);
//     int geom_random = d(gen);

//     int sign = dis (gen)*2 -1;
//     return geom_random * sign;
// }


int PrivateMisraGries::generateSymmetricGeometricRandomNumber(double alpha) {
    std::geometric_distribution<> d((alpha - 1) / alpha);
    std::uniform_int_distribution<> dis(0, 1);
    int geom_random = d(gen);
    int sign = dis(gen) * 2 - 1;
    return geom_random * sign;
}

PCC::PCC(double epsilon, double lambda,double true_value, int W,const std::vector<std::string>& initialData) : epsilon(epsilon), lambda(lambda), W(W),currentTime(0) {
    this->W0 = ceil(lambda * W / 4);
    cout<<lambda<<" "<<W<<" "<<W0<<endl;
    this->l = log(4 / lambda);
    initBlocks(initialData);

}


void PCC::processItem(string stream) {
    currentTime++;

    Block* underConstructionLeaf = nullptr;
    if (!blocks[0].empty() && blocks[0].back().state == UNDER_CONSTRUCTION) {
        assert(blocks[0].back().pmg!= nullptr);
        underConstructionLeaf = &blocks[0].back();
    }

    if (!underConstructionLeaf) {
        Block newLeaf;
        newLeaf.state = UNDER_CONSTRUCTION;
        newLeaf.startTime = currentTime;
        newLeaf.endTime = currentTime+W0-1;
        newLeaf.buffer.push_back(stream);

        if (newLeaf.endTime - newLeaf.startTime + 1 == 1) {
            newLeaf.state = ACTIVE;
            newLeaf.pmg = make_unique<PrivateMisraGries>(epsilon, lambda);
            newLeaf.pmg->processStream(newLeaf.buffer);
            newLeaf.buffer.clear();
        } else {
            if (!newLeaf.pmg) {
                double epsilon_i = epsilon / pow(2,(l-0+1));
                double lambda_i = 1.0 / (pow(2,0) * (l+1));

                newLeaf.pmg = make_unique<PrivateMisraGries>(epsilon_i, lambda_i);
            }
        }

        blocks[0].push_back(move(newLeaf));

        for (int i = 1; i < l; i++) {

            double epsilon_i = epsilon / pow(2,(l-i+1));
            double lambda_i = 1.0 / (pow(2,i) * (l+1));

            Block& lastBlock = blocks[i].back();
            if (lastBlock.state == UNDER_CONSTRUCTION) {
                lastBlock.buffer.push_back(stream);
                if (lastBlock.buffer.size() == W0*pow(2, i)) {
                    lastBlock.state = ACTIVE;
                    if (!lastBlock.pmg){
                        lastBlock.pmg = make_unique<PrivateMisraGries>(epsilon_i, lambda_i);
                    }
                    lastBlock.pmg->processStream(lastBlock.buffer);
                    lastBlock.buffer.clear();
                }
            } else if (lastBlock.state == ACTIVE) {
                Block newBlock;
                newBlock.state = UNDER_CONSTRUCTION;
                newBlock.startTime = currentTime;
                newBlock.endTime = currentTime + W0*pow(2, i) - 1;
                newBlock.buffer.push_back(stream);
                if (!newBlock.pmg) {
                    newBlock.pmg = make_unique<PrivateMisraGries>(epsilon_i, lambda_i);
                }
                blocks[i].push_back(move(newBlock));
            }
        }
    } else {
        underConstructionLeaf->buffer.push_back(stream);
        if (underConstructionLeaf->buffer.size() == W0*pow(2, 0)) {
            underConstructionLeaf->state = ACTIVE;
            underConstructionLeaf->pmg->processStream(underConstructionLeaf->buffer);
            underConstructionLeaf->buffer.clear();
        }
        for (int i = 1; i < l; i++) {

            double epsilon_i = epsilon / pow(2,(l-i+1));
            double lambda_i = 1.0 / (pow(2,i) * (l+1));

            Block& lastBlock = blocks[i].back();
            if (lastBlock.state == UNDER_CONSTRUCTION) {
                lastBlock.buffer.push_back(stream);
                if (lastBlock.buffer.size() == W0*pow(2, i)) {
                    lastBlock.state = ACTIVE;
                    if (!lastBlock.pmg){
                        lastBlock.pmg = make_unique<PrivateMisraGries>(epsilon_i, lambda_i);

                    }
                    lastBlock.pmg->processStream(lastBlock.buffer);
                    lastBlock.buffer.clear();
                }
            } else if (lastBlock.state == ACTIVE) {
                Block newBlock;
                newBlock.state = UNDER_CONSTRUCTION;
                newBlock.startTime = currentTime;
                newBlock.endTime = currentTime + W0*pow(2, i) - 1;
                newBlock.buffer.push_back(stream);
                if (!newBlock.pmg) {
                    newBlock.pmg = make_unique<PrivateMisraGries>(epsilon_i, lambda_i);
                }
                blocks[i].push_back(move(newBlock));
            }
        }

    }

    if (currentTime % W0 == 0){
        deleteExpiredBlocks();
    }
}



std::unordered_map<string, double> PCC::query_all() {
    std::vector<std::pair<int, int>> uncoveredTimeRanges;

    int mod = currentTime % W0;
    cout<<"current:"<<currentTime-mod<<endl;
    uncoveredTimeRanges.push_back({currentTime -mod - W + 1, currentTime - mod});

    // 结果map
    std::unordered_map<string, double> resultMap;

    for (int i = blocks.size() - 1; i >= 0 && !uncoveredTimeRanges.empty(); --i) {
        auto& currentLayer = blocks[i];
        std::vector<std::pair<int, int>> newUncoveredTimeRanges;

        for (Block& block : currentLayer) {
            std::vector<std::pair<int, int>> tempUncovered;  //int windowEndTime = it->second;

            for (auto& range : uncoveredTimeRanges ){
                if (block.endTime<=range.second && block.startTime >= range.first && block.state == ACTIVE){
                    auto blockResult = block.pmg->getResult();
                    cout<<"size:"<<blockResult.size()<<endl;
                    //cout<<"result size1:"<<resultMap.size()<<endl;
                    for (const auto& pair : blockResult) {
                        resultMap[pair.first] += pair.second;
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
    cout<<"all resultMap size:"<<resultMap.size()<<endl;
    return resultMap;
}


double PCC::query(string& element){
    unordered_map<string, double> resMap = query_all();
    double res = 0;
    auto it = resMap.find(element);
    if (it!=resMap.end()){
        res = it->second;
    }else{
        res = -100;
        cout<<"can't find"<<endl;
    }
    return res;
}


void PCC::deleteExpiredBlocks() {
    for (int i=0;i<l;i++){
        for (auto it = blocks[i].begin(); it != blocks[i].end();) {
            if (it->startTime <= currentTime - W) {  // 当前时间超过窗口W的结束时间，表示block已经过期
                it = blocks[i].erase(it);
            } else {
                ++it;
            }
        }
    }
}


void PCC::initBlocks(const std::vector<std::string>& initialData) {
    int currentLevelBlockSize = W0;

    for (int i = 0; i < l; i++) {
        // Calculate the number of blocks for this level
        int numOfBlocks = ceil(W*1.0 / currentLevelBlockSize);
        cout<<i<<" level "<<numOfBlocks<<endl;

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


        blocks.push_back(move(currentLevel));

        // Double the block size for the next level
        currentLevelBlockSize *= 2;
    }

    for (int i=0 ; i<l; i++){

        double epsilon_i = epsilon / pow(2,(l-i+1));
        double lambda_i = 1.0 / (pow(2,i) * (l+1));

        for(auto &block : blocks[i]) {
            if(block.endTime <= initialData.size()) {
                // This block is entirely within the initialData
                block.buffer.insert(block.buffer.end(),
                                    initialData.begin() + block.startTime - 1,
                                    initialData.begin() + block.endTime);
                block.pmg = make_unique<PrivateMisraGries>(epsilon_i,lambda_i);
                block.pmg->processStream(block.buffer);
                block.state = ACTIVE;
                block.buffer.clear();  // Clear the buffer after processing
            } else if(block.startTime <= initialData.size()) {
                // This block is partially filled with initialData
                block.buffer.insert(block.buffer.end(),
                                    initialData.begin() + block.startTime - 1,
                                    initialData.end());
                block.pmg = make_unique<PrivateMisraGries>(epsilon_i, lambda_i);
                block.pmg->processStream(block.buffer);
                block.state = UNDER_CONSTRUCTION;
            }
        }
    }
    currentTime = W;
}

