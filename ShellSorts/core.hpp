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

typedef unsigned long long ull;
typedef std::vector<int> vi;
typedef std::vector<ull> vull;

#include "formatTime.hpp"
#include "formatMicroSeconds.hpp"
#include "verifyArray.hpp"
#include "writeout.hpp"
#include "shell.hpp"

using namespace std::chrono;

struct my_numpunct : std::numpunct<char> {
    std::string do_grouping() const {return "\03";}
};

struct sortMetrics {
    ull time;
    ull sampleSize;
};
typedef std::vector<sortMetrics> vsm;

struct gapStruct {
    vull gaps;
    std::string name;
    enum errorState {
        ok = 0,
        outOfOrder = 1,
        deactivated = 1 << 1,
        unknown = 1 << 31,
    } status;
    std::function<void(vull &, ull)> gapFn;
    vsm runData;
};
typedef std::vector<gapStruct> vgs;

void setup();

void shell1959(vull &, ull);
void frank1960(vull &, ull);
void hibbard1963(vull &, ull);
void papernov1965(vull &, ull);
void pratt1971(vull &, ull);
void kunth1973(vull &, ull);
void sedgewick1982(vull &, ull);
void sedgewick1985(vull &, ull);
void sedgewick1986(vull &, ull);
void gonnet1991(vull &, ull);
void tokuda1992(vull &, ull);
void empirical2001(vull &, ull);
void huffman2022(vull &, ull);


#endif /* core_hpp */
