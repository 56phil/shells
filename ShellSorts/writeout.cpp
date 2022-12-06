//
//  writeout.cpp
//  writeout
//
//  Created by Phil Huffman on 12/3/21.
//

#include "writeout.hpp"
#include "core.hpp"

void randomWrite(const std::string fn, const ul nExperiments ) {
    int minInt(std::numeric_limits<int>::min()), maxInt(std::numeric_limits<int>::min());
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(minInt, maxInt);
    std::ofstream ofs;
     
    ofs.open (fn, std::ofstream::out);
    
    for (int i(0); i < nExperiments; ++i) {
        int n(distribution(generator));
        ofs << n << '\n';
    }
    
    ofs << std::endl;
    
    ofs.close();
}

int getRandyBe() {
    std::bernoulli_distribution dist(0.5);
    std::random_device rd;
    return dist(rd);
}

int getRandyBi() {
    std::binomial_distribution<> dist(1000, 0.5);
    std::random_device rd;
    return dist(rd);
}

int getRandyN() {
    std::normal_distribution<> dist(10.0, 5.0);
    std::random_device rd;
    return dist(rd);
}

int getRandyP() {
    std::poisson_distribution<> dist(1000.0);
    std::random_device rd;
    return dist(rd);
}

int getRandyU() {
    std::uniform_int_distribution<> dist(rMin, rMax);
    std::random_device rd;
    return dist(rd);
}

void randomFill(ul n, vi &v, std::string distroName) {
    int rMin(std::numeric_limits<int>::min()), rMax(std::numeric_limits<int>::max());
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    
    if (distroName == "Normal") {
        std::generate_n(v.begin(), n, getRandyN);
    } else if(distroName == "Poisson") {
        std::generate_n(v.begin(), n, getRandyP);
    } else if(distroName == "Bernoulli") {
        std::generate_n(v.begin(), n, getRandyBe);
    } else if(distroName == "Binomial") {
        std::generate_n(v.begin(), n, getRandyBi);
    } else if(distroName == "Uniform") {
        std::generate_n(v.begin(), n, getRandyU);
    } else if(distroName == "Uniform - Sorted") {
        std::generate_n(v.begin(), n, getRandyU);
        std::sort(v.begin(), v.end());
    } else if(distroName == "Uniform - Sorted & Reversed") {
        std::generate_n(v.begin(), n, getRandyU);
        std::sort(v.begin(), v.end());
        std::reverse(v.begin(), v.end());
    } else {
        std::cerr << "Unknown distribution requested (" << distroName << "). Using uniform." << std::endl;
        std::generate_n(v.begin(), n, getRandyU);
    }
}
