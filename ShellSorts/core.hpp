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

using namespace std::chrono;

const bool FULL_Run(true);
const bool WARN_Lagards(false);
const double dMax(std::numeric_limits<double>::max());
const double dMin(std::numeric_limits<double>::min());
const int MAX_DistroLines(250);
const int MAX_Passes(1);
const int MAX_Warnings(5);
const int iMax(std::numeric_limits<int>::max());
const int iMin(std::numeric_limits<int>::min());
const int GAPPER_Length(29);
const int DISTRO_Length(27);
const int FORMATTED_MicroSecondLength(13);
const int MICROSECOND_length(12);
const ul MAX_SampleSize(30000000);
const ul MIN_SampleSize(1000);
const ul MEDIAN_TrialSize(3);
const ul MIN_ActiveGapStructs(5);
const ul ulMax(std::numeric_limits<ul>::max());
const vs DISTRO_NAMES({
    "Bernoulli",
    "Binomial",
    "Normal",
    "Poisson",
    "Uniform",
    "Uniform - Sorted",
    "Uniform - Sorted & Reversed"
});
const std::string FN_Base("/Users/prh/Keepers/code/xCode/shells/results/");

struct topGapper {
    std::string distro;
    std::string gapper;
    ul time;
};
typedef topGapper tg;

struct my_numpunct : std::numpunct<char> {
    std::string do_grouping() const {return "\03";}
};

struct distroData {
    ul time;
};
typedef std::map<std::string,distroData> m_s_dd;

struct runData {
    vul gaps;
    m_s_dd dData;
};
typedef std::map<ul,runData> m_ul_rd;

struct gapStruct {
    int warnings;
    enum errorState {
        ok = 0,
        outOfOrder = 1,
        unknown = 1 << 31,
    } status;
    std::function<void(vul &, ul)> gapFn;
    m_ul_rd results;
};
typedef std::map<std::string,gapStruct> m_s_gs;

typedef std::map<ul,vul> m_ul_vul;

struct distroStruct {
    m_ul_vul originals;
};
typedef std::map<std::string,distroStruct> m_s_ds;

void init();
void shell(vul &, ul);
void frank(vul &, ul);
void hibbard(vul &, ul);
void papernov(vul &, ul);
void pratt(vul &, ul);
void pratt_A(vul &, ul);
void knuth(vul &, ul);
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
void shellSort(vul &, vul &);
void randomFill(ul, vul &, std::string);

#endif /* core_hpp */
