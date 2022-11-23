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

void randomFill(ul n, vi &v, std::string distroName) {
    int rMin(std::numeric_limits<int>::min()), rMax(std::numeric_limits<int>::max());
    std::default_random_engine generator;
    std::normal_distribution<int> distU(rMin, rMax);
    std::uniform_int_distribution<int> distN(rMin,rMax);
    std::bernoulli_distribution distBe(0.5);
    std::binomial_distribution<int> distBi(rMin, rMax);
    std::lognormal_distribution<int> distLo(rMin, rMax);
    std::chi_squared_distribution<int> distCh(2.8);
    std::cauchy_distribution<int> distCa(rMin, rMax);
    std::fisher_f_distribution<int> distFi(rMin, rMax);
    std::student_t_distribution<int> distSt(42.0);
    
    if (distroName == "Normal") {
        while (n--) {
            int r(distU(generator));
            v.push_back(r);
        }
    } else if(distroName == "Bernoulli") {
        while (n--) {
            int r(distBe(generator));
            v.push_back(r);
        }
    } else if(distroName == "Student T") {
        while (n--) {
            int r(distSt(generator));
            v.push_back(r);
        }
    } else if(distroName == "Fisher F") {
        while (n--) {
            int r(distFi(generator));
            v.push_back(r);
        }
    } else if(distroName == "Cauchy") {
        while (n--) {
            int r(distCa(generator));
            v.push_back(r);
        }
    } else if(distroName == "Chi Squared") {
        while (n--) {
            int r(distCh(generator));
            v.push_back(r);
        }
    } else if(distroName == "Lognormal") {
        while (n--) {
            int r(distBi(generator));
            v.push_back(r);
        }
    } else if(distroName == "Binomial") {
        while (n--) {
            int r(distBi(generator));
            v.push_back(r);
        }
    } else {
        while (n--) {
            int r(distU(generator));
            v.push_back(r);
        }
    }
}
