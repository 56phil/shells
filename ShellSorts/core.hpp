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
typedef std::vector<u_long> vul;

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
    u_long time;
    u_long sampleSize;
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
    std::function<void(vul &, u_long)> gapFn;
    vsm runData;
};
typedef std::vector<gapStruct> vgs;

static void getGaps(vgs &, u_long);
static void cullAlgorithms(vgs &, u_long, u_long);
static void cullAlgorithms(vgs &, u_long, u_long);
static void errorFunction(vi, vi);
static void makeFile(vgs &);
void setup();

void shell1959(vul &, u_long);
void frank1960(vul &, u_long);
void hibbard1963(vul &, u_long);
void papernov1965(vul &, u_long);
void pratt1971(vul &, u_long);
void kunth1973(vul &, u_long);
void sedgewick1982(vul &, u_long);
void sedgewick1986(vul &, u_long);
void gonnet1991(vul &, u_long);
void tokuda1992(vul &, u_long);
void empirical2001(vul &, u_long);
void huffman2022(vul &, u_long);


#endif /* core_hpp */
