//
//  writeout.cpp
//  writeout
//
//  Created by Phil Huffman on 12/3/21.
//

#include "writeout.hpp"

void randomWrite(const std::string fn, const long nExperiments ) {
    
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(std::numeric_limits<int>::min(),
                                                    std::numeric_limits<int>::max());
    std::ofstream ofs;
     
    ofs.open (fn, std::ofstream::out);
    
    for (int i(0); i < nExperiments; ++i) {
        int n(distribution(generator));
        ofs << n << '\n';
    }
    
    ofs << std::endl;
    
    ofs.close();
}

void randomFill(int n, vi &v) {
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(std::numeric_limits<int>::min(),
                                                    std::numeric_limits<int>::max());
    while (n--) {
        int r(distribution(generator));
        v.push_back(r);
    }
}
