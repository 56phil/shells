//
//  core.hpp
//  core
//
//  Created by Phil Huffman on 12/21/21.
//

#ifndef core_hpp
#define core_hpp

#include <chrono>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <locale>
#include <random>
#include <set>
#include <sstream>
#include <vector>

typedef std::vector<int> vi;

#include "formattime.hpp"
#include "formatMicroSeconds.hpp"
#include "verifyarray.hpp"
#include "writeout.hpp"
#include "shell.hpp"

using namespace std::chrono;

struct my_numpunct : std::numpunct<char> {
    std::string do_grouping() const {return "\03";}
};

struct sortMetrics {
    long time;
    long sampleSize;
    enum errorState {
        ok = 0,
        outOfOrder = 0x00000001,
        unknown = 0x80000000,
    } status;
};

struct gapStruct {
    vi gaps;
    std::string name;
    std::function<void(vi &, int)> gapFn;
    std::vector<sortMetrics> runData;
};

void errorFunction(vi, vi);
void makeFile(std::vector<gapStruct>);
void setup();

void shell1959(vi &, int);
void frank1960(vi &, int);
void hibbard1963(vi &, int);
void papernov1965(vi &, int);
void pratt1971(vi &, int);
void kunth1973(vi &, int);
void sedgewick1982(vi &, int);
void sedgewick1985(vi &, int);
void gonnet1991(vi &, int);
void empirical2001(vi &, int);


#endif /* core_hpp */
