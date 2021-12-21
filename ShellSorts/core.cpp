//
//  core.cpp
//  core
//
//  Created by Phil Huffman on 12/21/21.
//

#include "core.hpp"
struct sortMetrics {
    enum errorState {
        aOK = 0,
        order = 1,
    };
    long time;
    long sampleSize;
};

struct sortStruct {
    std::string name;
    std::function<vi(int)> gapFn;
    vi gaps;
    std::vector<sortMetrics> runData;
};

struct my_numpunct : std::numpunct<char> {
    std::string do_grouping() const {return "\03";}
};

