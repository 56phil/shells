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
typedef std::vector<ull> vul;

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
    ull time;
    ull sampleSize;
};
typedef std::vector<sortMetrics> vsm;

struct gapStruct {
    vul gaps;
    std::string name;
    enum errorState {
        ok = 0,
        outOfOrder = 1,
        deactivated = 1 << 1,
        unknown = 1 << 31,
    } status;
    std::function<void(vul &, ull)> gapFn;
    vsm runData;
};
typedef std::vector<gapStruct> vgs;

void setup();

void shell1959(vul &, ull);
void frank1960(vul &, ull);
void hibbard1963(vul &, ull);
void papernov1965(vul &, ull);
void pratt1971(vul &, ull);
void kunth1973(vul &, ull);
void sedgewick1982(vul &, ull);
void sedgewick1985(vul &, ull);
void sedgewick1986(vul &, ull);
void gonnet1991(vul &, ull);
void tokuda1992(vul &, ull);
void empirical2001(vul &, ull);
void huffman2022(vul &, ull);


#endif /* core_hpp */
