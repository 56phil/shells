//
//  core.hpp
//  core
//
//  Created by Phil Huffman on 12/21/21.
//

#ifndef core_hpp
#define core_hpp

#include <algorithm>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <locale>
#include <list>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <vector>

typedef unsigned long ul;
typedef unsigned long long ull;
typedef std::vector<int> vi;
typedef std::vector<double> vd;
typedef std::vector<ul> vul;
typedef std::vector<long> vl;
typedef std::map<std::string, vi> msvi;
typedef std::vector<std::string> vs;

#include "formatMicroSeconds.hpp"
#include "formatTime.hpp"
//#include "median.hpp"
#include "shell.hpp"
#include "verifyArray.hpp"
#include "writeout.hpp"

using namespace std::chrono;

const ul MEDIAN_TrialSize(3);  // keep this number ODD!
const ul MAX_SampleSize(250000000);
const ul MIN_ActiveGapStructs(5);
const long MAX_DistroLines(3250);
const int MAX_Warnings(5);
const int MAX_Passes(1);
const bool CULL_SlowerGapSequences(false);
const bool FULL_Run(true);
const int rMin(std::numeric_limits<int>::min());
const int rMax(std::numeric_limits<int>::max());
const std::string FN_Base("/Users/prh/Keepers/code/xCode/shells/results/");

struct my_numpunct : std::numpunct<char> {
    std::string do_grouping() const {return "\03";}
};

struct sortMetrics {
    long time;
    ul sampleSize;
};
typedef std::vector<sortMetrics> vsm;
typedef std::map<std::string, vsm> msm;

struct gapStruct {
    int warnings;
    vul gaps;
    std::string name;
    enum errorState {
        ok = 0,
        outOfOrder = 1,
        unknown = 1 << 31,
    } status;
    std::function<void(vul &, ul)> gapFn;
    msm runData;
};
//typedef std::vector<gapStruct> vgs;
typedef std::list<gapStruct> lgs;

void setup();
void shell(vul &, ul);
void frank(vul &, ul);
void hibbard(vul &, ul);
void papernov(vul &, ul);
void pratt(vul &, ul);
void pratt_A(vul &, ul);
void kunth(vul &, ul);
void sedgewick82(vul &, ul);
void sedgewick1985(vul &, ul);
void sedgewick86(vul &, ul);
void gonnet(vul &, ul);
void tokuda(vul &, ul);
void ciura(vul &, ul);
void a(vul &, ul);
void b(vul &, ul);
void c(vul &, ul);
void d(vul &, ul);
void e(vul &, ul);

#endif /* core_hpp */
