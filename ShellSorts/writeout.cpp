//
//  writeout.cpp
//  writeout
//
//  Created by Phil Huffman on 12/3/21.
//

#include "writeout.hpp"
#include "core.hpp"

void randomWrite(const std::string fn, const ul nExperiments ) {
    int minInt(std::numeric_limits<double>::min()), maxInt(std::numeric_limits<double>::min());
    std::default_random_engine generator;
    std::uniform_int_distribution<double> distribution(minInt, maxInt);
    std::ofstream ofs;
     
    ofs.open (fn, std::ofstream::out);
    
    for (int i(0); i < nExperiments; ++i) {
        int n(distribution(generator));
        ofs << n << '\n';
    }
    
    ofs << std::endl;
    
    ofs.close();
}

void uniFill(vd &v, ul n) {
    srand(time(NULL));
    while (n--) {
        v.push_back(rand());
    }
}

void randomFill(ul n, vd &v, std::string distroName) {
    int rMin(std::numeric_limits<double>::min()), rMax(std::numeric_limits<double>::max());
    std::default_random_engine generator;
    std::normal_distribution<double> distU(rMin, rMax);
    std::uniform_int_distribution<double> distN(rMin,rMax);
    std::bernoulli_distribution distBe(0.5);
    std::binomial_distribution<double> distBi(rMin, rMax);
    std::lognormal_distribution<double> distLo(rMin, rMax);
    std::chi_squared_distribution<double> distCh(2.8);
    std::cauchy_distribution<double> distCa(rMin, rMax);
    std::fisher_f_distribution<double> distFi(rMin, rMax);
    std::student_t_distribution<double> distSt(42.0);
    
    if (distroName == "Normal") {
        while (n--) {
            int r(distN(generator));
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
            int r(distLo(generator));
            v.push_back(r);
        }
    } else if(distroName == "Binomial") {
        while (n--) {
            int r(distBi(generator));
            v.push_back(r);
        }
    } else if(distroName == "Uniform") {
        uniFill(v, n);
    } else if(distroName == "Uniform - Sorted") {
        uniFill(v, n);
        std::sort(v.begin(), v.end());
    } else if(distroName == "Uniform - Sorted & Reversed") {
        uniFill(v, n);
        std::sort(v.begin(), v.end());
        std::reverse(v.begin(), v.end());
    } else {
        std::cerr << "Unknown distribution requested. Using uniform." << std::endl;
        uniFill(v, n);
    }
}
