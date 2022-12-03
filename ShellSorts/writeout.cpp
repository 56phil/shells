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

void uniFill(vi &v, ul n, std::uniform_int_distribution<int> distrib, std::mt19937 gen) {
    while (n--) {
        v.push_back(distrib(gen));
    }
}

void randomFill(ul n, vi &v, std::string distroName) {
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> distrib(1, 6);

    std::random_device r;
 
    std::default_random_engine e1(r());
 
    std::seed_seq seed2{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937 e2(seed2);

    
    int rMin(std::numeric_limits<int>::min()), rMax(std::numeric_limits<int>::max());
    std::default_random_engine generator;
    std::normal_distribution<int> distN(0.0, 100.0);
    std::uniform_int_distribution<int> distU(rMin,rMax);
    std::bernoulli_distribution distBe(0.5);
    std::binomial_distribution<int> distBi(4, 0.5);
    std::lognormal_distribution<int> distLo(0.5, 0.25);
    std::chi_squared_distribution<int> distCh(2.8);
    std::cauchy_distribution<int> distCa(rMin, 0.25);
    std::fisher_f_distribution<int> distFi(rMin, rMax);
    std::student_t_distribution<int> distSt(42.0);
    
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
    } else if(distroName == "Binomial") {
        while (n--) {
            int r(distBi(generator));
            v.push_back(r);
        }
    } else if(distroName == "Uniform") {
        uniFill(v, n, distrib, gen);
    } else if(distroName == "Uniform - Sorted") {
        uniFill(v, n, distrib, gen);
        std::sort(v.begin(), v.end());
    } else if(distroName == "Uniform - Sorted & Reversed") {
        uniFill(v, n, distrib, gen);
        std::sort(v.begin(), v.end());
        std::reverse(v.begin(), v.end());
    } else {
        std::cerr << "Unknown distribution requested. Using uniform." << std::endl;
        uniFill(v, n, distrib, gen);
    }
}
