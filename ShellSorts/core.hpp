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
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <vector>

typedef unsigned long long ul;
typedef std::vector<int> vi;
typedef std::vector<ul> vul;

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
    ul time;
    ul sampleSize;
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
    std::function<void(vul &, ul)> gapFn;
    vsm runData;
};
typedef std::vector<gapStruct> vgs;

void setup();

void shell1959(vul &, ul);
void frank1960(vul &, ul);
void hibbard1963(vul &, ul);
void papernov1965(vul &, ul);
void pratt1971(vul &, ul);
void kunth1973(vul &, ul);
void sedgewick1982(vul &, ul);
void sedgewick1985(vul &, ul);
void sedgewick1986(vul &, ul);
void gonnet1991(vul &, ul);
void tokuda1992(vul &, ul);
void empirical2001(vul &, ul);
void huffman2022(vul &, ul);


#endif /* core_hpp */
