//
//  writeout.cpp
//  writeout
//
//  Created by Phil Huffman on 12/3/21.
//

#include "writeout.hpp"
#include "core.hpp"

//void randomWrite(const std::string fn, const ul nExperiments ) {
//    int minInt(std::numeric_limits<int>::min()), maxInt(std::numeric_limits<int>::min());
//    std::default_random_engine generator;
//    std::uniform_int_distribution<int> distribution(minInt, maxInt);
//    std::ofstream ofs;
//
//    ofs.open (fn, std::ofstream::out);
//
//    for (int i(0); i < nExperiments; ++i) {
//        int n(distribution(generator));
//        ofs << n << '\n';
//    }
//
//    ofs << std::endl;
//
//    ofs.close();
//}

void getRandyBe(vi &v, ul n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::bernoulli_distribution dist(0.5);
    while (n--) {
        v.push_back(dist(rd));
    }
}

void getRandyBi(vi &v, ul n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::binomial_distribution<int> dist(1000, 0.5);
    while (n--) {
        v.push_back(dist(rd));
    }
}

void getRandyN(vi &v, ul n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(10.0, 5.0);
    while (n--) {
        v.push_back(dist(rd));
    }
}

void getRandyP(vi &v, ul n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::poisson_distribution<int> dist(1000.0);
    while (n--) {
        v.push_back(dist(rd));
    }
}

void getRandyU(vi &v, ul n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(rMin, rMax);
    while (n--) {
        v.push_back(dist(gen));
    }
}

void randomFill(ul n, vi &v, std::string distroName) {
    v.clear();
    if (distroName == "Normal") {
        getRandyN(v,n);
    } else if(distroName == "Poisson") {
        getRandyP(v,n);
    } else if(distroName == "Bernoulli") {
        getRandyBe(v,n);
    } else if(distroName == "Binomial") {
        getRandyBi(v,n);
    } else if(distroName == "Uniform") {
        getRandyU(v,n);
    } else if(distroName == "Uniform - Sorted") {
        getRandyU(v,n);
        std::sort(v.begin(), v.end());
    } else if(distroName == "Uniform - Sorted & Reversed") {
        getRandyU(v,n);
        std::sort(v.begin(), v.end());
        std::reverse(v.begin(), v.end());
    } else {
        std::cerr << "Unknown distribution requested (" << distroName << "). Using uniform." << std::endl;
        getRandyU(v,n);
    }
}
