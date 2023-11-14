//
// Created by 12293 on 2023/10/12.
//

#include "PCC.h"


PrivateMisraGries::PrivateMisraGries(double epsilon, double lambda, int T)
        : epsilon(epsilon), lambda(lambda), T(T), threshold(5) {}

void PrivateMisraGries::processStream(const std::vector<string>& stream) {
    for(int t = 0; t < stream.size(); t++) {
        string sigma = stream[t];
        f[sigma]++;

        // 查找此时f(x)>0有多少个
        int countOverThreshold = 0;
        for(auto& [key, value] : f) {
            if(value > threshold) countOverThreshold++;
        }

        if(countOverThreshold > threshold) {
            for(auto& [key, value] : f) {
                if(value > 0) value--;
            }
        }
    }

//    int sensitivity = threshold + 1;
//    std::geometric_distribution<int> distribution(std::exp(-epsilon / (sensitivity + 1)));
//    std::default_random_engine generator;
//
//    for(auto& [key, value] : f) {
//        int r = distribution(generator);
//        f[key] = std::max(f[key] + r, 0.0);
//    }

    int nonZeroCount = 0;
    string smallestNonZeroItem;
    int smallestNonZeroValue = INT_MAX;

    for(auto& [key, value] : f) {
        if(value > 0) {
            nonZeroCount++;
            if(value < smallestNonZeroValue) {
                smallestNonZeroValue = value;
                smallestNonZeroItem = key;
            }
        }
    }

    if(nonZeroCount > threshold) {
        f[smallestNonZeroItem] = 0;
    }
}

std::unordered_map<string, double> PrivateMisraGries::getResult() {
    return f;
}





PCC::PCC(double epsilon, double lambda, int W,const std::vector<std::string>& initialData) : epsilon(epsilon), lambda(lambda), W(W),currentTime(0) {
    this->W0 = ceil(lambda * W / 4);
    this->l = log2(4 / lambda);
    initBlocks(initialData);
//    for (int i=0;i<initialData.size();i++){
//        processItem(initialData[i]);
//    }
}


//void PCC::processItem(string stream) {
//    currentTime++;
//    for (int i = 0; i <l; i++) {
//        for (auto &block: blocks[i]) {
//            if (block.startTime > currentTime) {
//                block.state = FUTURE;
//            } else if (block.startTime <= currentTime && block.endTime > currentTime) {
//                block.state = UNDER_CONSTRUCTION;
//                if (!block.pmg) { // Only initiate once
//                    block.pmg = make_unique<PrivateMisraGries>(epsilon, lambda, block.endTime - block.startTime + 1);
//                }
//                block.buffer.push_back(stream);
//            }else if (currentTime - W + 1 <= block.startTime && block.endTime <= currentTime && block.state!= ACTIVE) {
//                block.state = ACTIVE;
//                if (!block.pmg){
//                    block.pmg = make_unique<PrivateMisraGries>(epsilon, lambda, block.endTime - block.startTime + 1);
//                    block.pmg->processStream({stream});
//                }
//                if (block.buffer.size() > 0) {
//                    block.pmg->processStream(block.buffer);
//                    block.buffer.clear(); // Empty the buffer
//                }
//                // Now, you can retrieve results if you want
//                // auto result = block.pmg->getResult();
//            }else if (block.startTime < currentTime - W + 1) {
//                    block.state = EXPIRED;
//                }
//            }
//        }
//    }

void PCC::processItem(string stream) {
    currentTime++;

    Block* underConstructionLeaf = nullptr;
    for (auto &block : blocks[0]) {  // 只检查叶子节点
        if (block.state == UNDER_CONSTRUCTION) {
            underConstructionLeaf = &block;
            break;
        }
    }

    if (!underConstructionLeaf) {  // 如果没有UNDER_CONSTRUCTION的叶子节点
        // 直接为叶子节点创建一个新的block
        Block newLeaf;
        newLeaf.state = UNDER_CONSTRUCTION;
        newLeaf.startTime = currentTime;
        newLeaf.endTime = currentTime;
        newLeaf.buffer.push_back(stream);

        if (newLeaf.endTime - newLeaf.startTime + 1 == 1) {  // 如果叶子节点只需要一个元素
            newLeaf.state = ACTIVE;
            newLeaf.pmg = make_unique<PrivateMisraGries>(epsilon, lambda, 1);
            newLeaf.pmg->processStream(newLeaf.buffer);
            newLeaf.buffer.clear();
        } else {
            if (!newLeaf.pmg) {
                newLeaf.pmg = make_unique<PrivateMisraGries>(epsilon, lambda, newLeaf.endTime - newLeaf.startTime + 1);
            }
        }

        blocks[0].push_back(move(newLeaf));

        // 对于更高的层次，根据最后一个block的状态进行操作
        for (int i = 1; i < l; i++) {
            Block& lastBlock = blocks[i].back();
            if (lastBlock.state == UNDER_CONSTRUCTION) {
                lastBlock.buffer.push_back(stream);
                if (lastBlock.buffer.size() == W0*pow(2, i)) {
                    lastBlock.state = ACTIVE;
                    if (!lastBlock.pmg){
                        lastBlock.pmg = make_unique<PrivateMisraGries>(epsilon, lambda, 1);
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
                    newBlock.pmg = make_unique<PrivateMisraGries>(epsilon, lambda, newBlock.endTime - newBlock.startTime + 1);
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
            Block& lastBlock = blocks[i].back();
            if (lastBlock.state == UNDER_CONSTRUCTION) {
                lastBlock.buffer.push_back(stream);
                if (lastBlock.buffer.size() == W0*pow(2, i)) {
                    lastBlock.state = ACTIVE;
                    if (!lastBlock.pmg){
                        lastBlock.pmg = make_unique<PrivateMisraGries>(epsilon, lambda, 1);

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
                    newBlock.pmg = make_unique<PrivateMisraGries>(epsilon, lambda, newBlock.endTime - newBlock.startTime + 1);
                }
                blocks[i].push_back(move(newBlock));
            }
        }

    }


    for (int i=0;i<l;i++){
        for (auto it = blocks[i].begin(); it != blocks[i].end();) {
            if (it->startTime <= currentTime - W) {  // 当前时间超过窗口W的结束时间，表示block已经过期
                it = blocks[i].erase(it);
            } else {
                ++it;
            }
        }
    }

    deleteExpiredBlocks();
}



std::unordered_map<string, double> PCC::query_all() {
    std::vector<std::pair<int, int>> uncoveredTimeRanges;
    uncoveredTimeRanges.push_back({currentTime - W + 1, currentTime});

    // 结果map
    std::unordered_map<string, double> resultMap;

    for (int i = blocks.size() - 1; i >= 0 && !uncoveredTimeRanges.empty(); --i) {
        auto& currentLayer = blocks[i];
        std::vector<std::pair<int, int>> newUncoveredTimeRanges;

        for (Block& block : currentLayer) {
            for (auto it = uncoveredTimeRanges.begin(); it != uncoveredTimeRanges.end(); ) {
                int windowStartTime = it->first;
                int windowEndTime = it->second;

                if (block.endTime >= windowStartTime && block.startTime <= windowEndTime && block.state == ACTIVE) {
                    auto blockResult = block.pmg->getResult();
                    for (const auto& pair : blockResult) {
                        resultMap[pair.first] += pair.second;
                    }

                    if (block.startTime > windowStartTime) {
                        newUncoveredTimeRanges.push_back({windowStartTime, block.startTime - 1});
                    }
                    if (block.endTime < windowEndTime) {
                        newUncoveredTimeRanges.push_back({block.endTime + 1, windowEndTime});
                    }

                    // 使用返回的迭代器来避免迭代器失效的问题
                    it = uncoveredTimeRanges.erase(it);
                } else {
                    ++it;
                }
            }
        }

        //
        uncoveredTimeRanges.insert(uncoveredTimeRanges.end(), newUncoveredTimeRanges.begin(), newUncoveredTimeRanges.end());
    }

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

    // 计算有多少叶子节点
//    int leaf_num = W/W0;
//    // 排满子树应该有多少叶子节点
//    int new_num = pow(2, ceil(log2(leaf_num)));
//    // 此时新的W值为
//    int new_w = new_num * W0;
    // 1. 按照初始化好第一个窗口中block的情况构建好跟窗口相关的所有block
    for (int i = 0; i < l; i++) {
        // Calculate the number of blocks for this level
        int numOfBlocks = ceil(W*1.0 / currentLevelBlockSize);

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
    for(auto &level : blocks) {
        for(auto &block : level) {
            if(block.endTime <= initialData.size()) {
                // This block is entirely within the initialData
                block.buffer.insert(block.buffer.end(),
                                    initialData.begin() + block.startTime - 1,
                                    initialData.begin() + block.endTime);
                block.pmg = make_unique<PrivateMisraGries>(epsilon, lambda, block.endTime - block.startTime + 1);
                block.pmg->processStream(block.buffer);
                block.state = ACTIVE;
                block.buffer.clear();  // Clear the buffer after processing
            } else if(block.startTime <= initialData.size()) {
                // This block is partially filled with initialData
                block.buffer.insert(block.buffer.end(),
                                    initialData.begin() + block.startTime - 1,
                                    initialData.end());
                if (block.endTime<=initialData.size()){
                    block.pmg = make_unique<PrivateMisraGries>(epsilon, lambda, block.buffer.size());
                    block.pmg->processStream(block.buffer);
                }
                block.state = UNDER_CONSTRUCTION;
            }
        }
    }
    currentTime = W;
}

