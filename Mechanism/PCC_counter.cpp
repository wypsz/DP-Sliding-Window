//
// Created by 12293 on 2023/12/10.
//


#include <cassert>
#include "PCC_counter.h"

MG::MG(double alpha) {
    this->alpha = alpha;
    this->threshold = ceil((1/ alpha));
}

void MG::processStream(const std::vector<string> &stream) {
    for(int t = 0; t < stream.size(); t++) {
        string sigma = stream[t];
        f[sigma]++;

        if(f.size() > threshold) {
            for(auto it = f.begin(); it != f.end(); ) {
                if(it->second > 0) it->second--;
                if(it->second == 0) {
                    it = f.erase(it);
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

    for (int i = 0; i <= L; i++) {
        int numOfBlocks = ceil(W*1.0 / currentLevelBlockSize);
        cout<<i<<" level "<<"numOfBlocks:"<<numOfBlocks<<"block size:"<<currentLevelBlockSize<<endl;


        vector<Block> currentLevel;

        for (int j = 0; j < numOfBlocks; j++) {
            Block block;
            block.startTime = j * currentLevelBlockSize + 1;
            block.endTime = block.startTime + currentLevelBlockSize - 1;
            block.state = FUTURE;
            currentLevel.push_back(move(block));
        }

        blocks.push_back(move(currentLevel));

        currentLevelBlockSize *= 2;
    }

    for (int i=0 ; i<=L; i++){

        double alpha_i = alpha / pow(2,(i+1));
        for(auto &block : blocks[i]) {
            if(block.endTime <= initialData.size()) {
                block.buffer.insert(block.buffer.end(),
                                    initialData.begin() + block.startTime - 1,
                                    initialData.begin() + block.endTime);
                block.pmg = make_unique<MG>(alpha_i);
                block.pmg->processStream(block.buffer);
                block.state = ACTIVE;
                block.fre_true = convertToFrequencyMap(block.buffer);
                block.buffer.clear();
            } else if(block.startTime <= initialData.size()) {
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

    Block* underConstructionLeaf = nullptr;
    if (!blocks[0].empty() && blocks[0].back().state == UNDER_CONSTRUCTION) {
        assert(blocks[0].back().pmg!= nullptr);
        underConstructionLeaf = &blocks[0].back();
    }

    if (!underConstructionLeaf) {
        Block newLeaf;
        newLeaf.state = UNDER_CONSTRUCTION;
        newLeaf.startTime = currentTime;
        newLeaf.endTime = currentTime+S1-1;
        newLeaf.buffer.push_back(stream);

        if (newLeaf.endTime - newLeaf.startTime + 1 == 1) {
            cout<<"1 item"<<endl;
            newLeaf.state = ACTIVE;
            newLeaf.pmg = make_unique<MG>(alpha);
            newLeaf.pmg->processStream(newLeaf.buffer);
            newLeaf.fre_true = convertToFrequencyMap(newLeaf.buffer);
            newLeaf.buffer.clear();
        } else {
            if (!newLeaf.pmg) {
                double alpha_i = alpha / pow(2,1);

                newLeaf.pmg = make_unique<MG>(alpha_i);
            }
        }

        blocks[0].push_back(move(newLeaf));

        for (int i = 1; i <= L; i++) {

            double alpha_i = alpha / pow(2,(i+1));

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

        }
        for (int i = 1; i <= L; i++) {

            double alpha_i = alpha / pow(2,(i+1));


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
                    double alpha_i = block.pmg->alpha;
                    double epsilon_i = (2/alpha_i) / (epsilon/((L+1)*1.0));
                    cout<<"alpha_i:"<<alpha_i<<"epsilon_i:"<<epsilon_i<<endl;
                    int x = 0;

                    for (const auto& pair : blockResult) {
                        x++;
                        int fre = block.fre_true[pair.first];
                        generator.seed(currentSeed);
                        double noise = generateLaplacianRandom(generator,epsilon_i,1)[0];
                        currentSeed++;
                        resultMap[pair.first] = resultMap[pair.first]+fre+noise;
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
    cout<<"original size:"<<resultMap.size()<<endl;


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
    for (int i=0;i<=L;i++){
        for (auto it = blocks[i].begin(); it != blocks[i].end();) {
            if (it->endTime < currentTime - W+1) {
                it = blocks[i].erase(it);
            } else {
                ++it;
            }
        }
    }
}

PCC_counter::PCC_counter(double epsilon,int level, double alpha, int W, const vector<std::string> &initialData) {
    this->alpha = alpha;
    this->W =W;
    this->currentTime = 0;
    this->epsilon = epsilon;
    this->initialSeed = 0;
    this->currentSeed = initialSeed;
    this->L = level;
    this->S1 = ceil(W / (pow(2,(level+1))));
    initBlocks(initialData);
}
