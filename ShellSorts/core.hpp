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

typedef unsigned long ul;
typedef unsigned long long ull;
typedef std::vector<int> vi;
typedef std::vector<double> vd;
typedef std::vector<ul> vul;
typedef std::map<std::string, vd> msvd;
typedef std::vector<std::string> vs;

#include "formatTime.hpp"
#include "formatMicroSeconds.hpp"
#include "verifyArray.hpp"
#include "writeout.hpp"
#include "shell.hpp"

using namespace std::chrono;

#define MAX_OUTER_LOOP 3
#define MAX_DISTRO_LINES 2500
#define MAX_SAMPLE_SIZE 99999999

struct my_numpunct : std::numpunct<char> {
    std::string do_grouping() const {return "\03";}
};

struct sortMetrics {
    ul time;
    ul sampleSize;
};
typedef std::vector<sortMetrics> vsm;
typedef std::map<std::string, vsm> msm;

struct gapStruct {
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
typedef std::vector<gapStruct> vgs;

void setup();

void shell(vul &, ul);
void frank(vul &, ul);
void hibbard(vul &, ul);
void papernov(vul &, ul);
void pratt(vul &, ul);
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

#endif /* core_hpp */
