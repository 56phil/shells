//
//  core.cpp
//
//  Created by Phil Huffman on 12/21/21.
//

#include "core.hpp"

template<typename T>
T median(std::vector<T> a) {
    if (a.empty())
        return -1;
    
    auto it(a.begin() + a.size() / 2);
    std::sort(a.begin(), a.end());
    if ((a.size() & 1) == 0) {
        return ((*it + *(it+1)) / 2);
    }
    return *it;
}

template<typename T>
T average(std::vector<T> a) {
    if (a.empty())
        return -1;
    T sum(0);
    for (auto tmp : a) {
        sum += tmp;
    }
    return sum / a.size();
}

std::string gaps2string(vul gaps, std::string delimiter) {
    std::stringstream sst;
    for (auto gap : gaps) {
        sst << delimiter << gap;
    }
    std::string oString(sst.str());
    return oString;
}

static void minFillGaps(vul &gaps, ul vSize) {
    gaps.clear();
    gaps.push_back(vSize >> 1);
    gaps.push_back(vSize >> 2);
    gaps.push_back(1);
}

static void inspectGaps(vul &gaps, const ul vSize) {
    if (gaps.empty()) {
        std::cerr << " \tGap error. Empty.\n";
        minFillGaps(gaps, vSize);
    }

    if (gaps.front() < gaps.back())
        std::reverse(gaps.begin(), gaps.end());

    if (gaps.back() != 1)
        gaps.push_back(1);

    while (gaps.front() >= vSize)
        gaps.erase(gaps.begin());

    for (auto it(gaps.begin() + 1); it != gaps.end();  it++) {
        if (*(it - 1) <= *it) {
            std::cerr << " \tGap Error. Sequence.\n";
            minFillGaps(gaps, vSize);
        }
    }

    gaps.shrink_to_fit();
}

static void writeDistros(m_s_ds dMap, vul sizes) {
//    for (auto &d0 : dMap) {
    if (dMap.empty()) {
        std::cerr << formatTime(true, true) << " \tdMap error. Empty.";
        exit(2);
    }
    
    auto maxPossibleLines(dMap.begin()->second.originals.begin()->second.size());   // size of the first (smallest) sample vector
    auto max_dl(MAX_DistroLines < maxPossibleLines ? MAX_DistroLines : maxPossibleLines);
    std::string fnBase("/Users/prh/Keepers/code/xCode/shells/results/");
    fnBase += formatTime(true, true);
    fnBase += "-distros.csv";
    std::fstream fst;
    fst.open(fnBase, std::ios::out);
    for(long indx(-1); indx < max_dl; indx++) {
        for (auto d0 : dMap) {
                if (indx < 0) {
                    fst << ',' << d0.first;
                } else {
                    fst << ',' << d0.second.originals[sizes[0]][indx];
                }
        }
        fst << '\n';
    }
    fst << std::endl;
    fst.close();
}

static void prepG1(m_s_gs &gMap, vul sizes) {
    for (auto &g0 : gMap) {
        for (auto dName : DISTRO_NAMES) {
            for (auto size : sizes) {
                g0.second.results[size].dData[dName].time = 0;
            }
        }
    }
}

static void makeGapSequences(m_s_gs &gMap) {
    for (auto &g0 : gMap) {
        for (auto &g1 : g0.second.results) {
            g1.second.gaps.clear();
            g0.second.gapFn(g1.second.gaps, g1.first);
            inspectGaps(g1.second.gaps, g1.first);
        }
    }
}

static void make_gMap(m_s_gs &gMap, vul sizes) {
    gMap["Shell 1959"].gapFn = shell;
    gMap["Frank & Lazarus 1960"].gapFn = frank;
    gMap["Hibbard 1963"].gapFn = hibbard;
    gMap["Papernov & Stasevich 1965"].gapFn = papernov;
//    gMap["Pratt 1971"].gapFn = pratt;
    gMap["Knuth 1973"].gapFn = knuth;
    gMap["Sedgewick 1982"].gapFn = sedgewick82;
    gMap["Sedgewick 1986"].gapFn = sedgewick86;
    gMap["Gonnet & Baeza-Yates 1991"].gapFn = gonnet;
    gMap["Tokuda 1992"].gapFn = tokuda;
    gMap["Ciura 2001"].gapFn = ciura;
    gMap["a 2022"].gapFn = a;
    gMap["b 2022"].gapFn = b;
    gMap["c 2022"].gapFn = c;
    gMap["d 2022"].gapFn = d;
    gMap["e 2022"].gapFn = e;
    
    prepG1(gMap, sizes);
    
    for (auto &g0 : gMap) {
        g0.second.status = gapStruct::ok;
        g0.second.warnings = 0;
        for (auto &g1 : g0.second.results) {
            g1.second.gaps.clear();
        }
    }
    
    makeGapSequences(gMap);
    std::cerr << formatTime(true, true) << " \tgMap constructed.\n";
}

static void make_dMap(m_s_ds &dMap, const vul sizes) {
    for (auto dName : DISTRO_NAMES) {
        for (auto size : sizes) {
            vul t_vul;
            randomFill(size, t_vul, dName);
            dMap[dName].originals[size] = t_vul;
        }
        std::cerr << formatTime(true, true) << " \tdMap constructed " << dName << ".\n";
    }
    std::cerr << formatTime(true, true) << " \tdMap constructed.\n";
}

static void errorFunction(vul &wc, vul &cc) {
    int n(0), w(28), maxLines(24);
    
    auto itw(wc.begin());
    auto itc(cc.begin());
    
    std::cout << std::right << std::setw(4) << "n"
    << std::right << std::setw(w) << "expected"
    << std::right << std::setw(w) << "result" << '\n';
    while (itw != wc.end() && itc != cc.end() && n < maxLines)
        std::cout << std::right << std::setw(4) << ++n
        << std::right << std::setw(w) << *itc++
        << std::right << std::setw(w) << *itw++ << '\n';
}

static void writeTimes(m_s_gs gMap, m_s_ds dMap) {
    if (FULL_Run) {
        std::fstream fst;
        std::string fileName(FN_Base + formatTime(true, true));
        fileName += "-Times.csv";
        fst.open(fileName, std::ios::out);
        fst << "Size,Distribution,Name,Time\n";
        for (auto g0 : gMap) {
            for (auto g1 : g0.second.results) {
                for (auto d0 : dMap) {
                    for (auto d1 : d0.second.originals) {
                        fst << d1.first << ',' << d0.first << ',' << g0.first << ',' << g1.second.dData[d0.first].time << '\n';
                    }
                }
            }
        }
        fst << std::endl;
        fst.close();
    }
}

static void writeGaps(m_s_gs gMap) {
    std::vector<std::string> alsoRans;
    std::fstream gst;
    std::string fnBase(FN_Base + formatTime(true, true));
    fnBase += "-Gaps.csv";
    gst.open(fnBase, std::ios::out);
    gst << "Gapper,Size,First,Second,..." << '\n';
    for (auto g0 : gMap) {
        for (auto g1 : g0.second.results) {
            gst << g0.first << ',' << g1.first << gaps2string(g1.second.gaps, ",") << '\n';
        }
    }
    gst << std::endl;
    gst.close();
}

static void findFastest(m_s_gs gMap, vul sizes) {
    typedef std::map<ul, std::string> q;
    std::map<std::string,q> ranked;  // map of fastest for size within fastest for distro
    std::cout << "\n;;;;;;;;;;;;;;;;;;;;;\n\n";
    
    for (auto dName : DISTRO_NAMES) {
        auto bestGapperDistroName(gMap.begin()->first);
        auto bestGapperDistro(gMap.begin()->second);
        for (auto size : sizes) {
            auto bestGapperSizeName(gMap.begin()->first);
            auto bestGapperSize(gMap.begin()->second);
            for (auto g0 : gMap) {  // evaluate each gapper for this distro
                for (auto g1 : g0.second.results) {
                    if (g1.second.dData[dName].time < bestGapperDistro.results[size].dData[dName].time) {
                        bestGapperDistro = g0.second;
                        bestGapperDistroName = g0.first;
                    }
                    if (g1.second.dData[dName].time < bestGapperSize.results[size].dData[dName].time) {
                        bestGapperSize = g0.second;
                        bestGapperSizeName = g0.first;
                    }
                }
            }
            std::cout << formatTime(true, true) << " Fastest for size: " << bestGapperSizeName << " \t" << size <<'\n';
        }
        std::cout << formatTime(true, true) << " Fastest for distro: " << bestGapperDistroName << " \t" << dName << '\n';
    }
}

static void summerize(m_s_gs gMap) {
    std::cout << "\n=======================\n\n";
    for (auto g0 : gMap) {  // gapper
        std::cout << "Gapper: " << g0.first;
        for (auto g1 : g0.second.results) {     // sample size
            std::cout << "\nn: "<< std::right << std::setw(11) << g1.first
            << " Gaps: " << gaps2string(g1.second.gaps, " ") << '\n';
            for (auto g2 : g1.second.dData) {   // distro
                std::cout << formatTime(true, true)
                << " Distribution: " << std::setw(27) << g2.first
                << " \tLapsed Time: " << std::right << std::setw(11) << formatMicroSeconds(g2.second.time) << '\n';
            }
        }
        std::cout << '\n';
    }
}

static void eoj(m_s_gs gMap, m_s_ds dMap, vul sizes) {
    summerize(gMap);
//    findFastest(gMap, sizes); this function needs a total rewrite
    writeDistros(dMap, sizes);
    writeGaps(gMap);
    writeTimes(gMap, dMap);
    std::cout << std::endl;
    std::cerr << std::endl;
}

static void warnLagards(const std::pair<const std::string, distroStruct> &d0, m_s_gs &gMap, const vul &gTimes) {
    const auto avgTime(average(gTimes));
    const auto lim(avgTime + (avgTime >> 2));
    for (auto g0 : gMap) {
        for (auto g1 : g0.second.results) {
            if (g1.second.dData[d0.first].time > lim) {
                g0.second.warnings++;
                std::cerr << formatTime(true, true) << " \tWarned " << g0.first << " (" << g0.second.warnings << " of " << MAX_Warnings << ")\n";
            }
        }
    }
}

static void testDur(long long dur, std::pair<const std::string, gapStruct> &g0, tg &tgx, std::string dName) {
    if (dur < tgx.time) {
        tgx.distro = dName;
        tgx.gapper = g0.first;
        tgx.time = dur;
        assert(tgx.gapper.size() > 1);
    }
}

static void printIt(const std::pair<const std::string, distroStruct> &d0, long dur, std::pair<const std::string, gapStruct> &g0, std::pair<const unsigned long, runData> &g1, vul &gTimes,
                    tg &tgs) {
    g1.second.dData[d0.first].time = dur;
    gTimes.push_back(dur);
    testDur(dur, g0, tgs, d0.first);
    std::cout << formatTime(false, true)
    << std::right << std::setw(GAPPER_Length) << g0.first
    << std::right << std::setw(MICROSECOND_length) << dur << "µs "
    << std::right << std::setw(FORMATTED_MicroSecondLength) << formatMicroSeconds(dur)
    << '\n';
}

static void doSort(const std::pair<const std::string, distroStruct> &d0, const std::pair<const unsigned long,
                   std::vector<unsigned long>> &d1, std::pair<const std::string, gapStruct> &g0, vul &gTimes, ul size,
                   tg &tgs) {
    for (auto &g1 : g0.second.results ) {
        if (g1.first == d1.first) {
            auto checkCopy(d1.second);
            std::sort(checkCopy.begin(), checkCopy.end());
            vl times(MEDIAN_TrialSize);
            times.clear();
            auto workCopy(d1.second);
            while (times.size() < MEDIAN_TrialSize) {
                auto start = high_resolution_clock::now();
                shellSort(workCopy, g0.second.results[d1.first].gaps);
                auto stop = high_resolution_clock::now();
                auto durT = duration_cast<microseconds>(stop - start).count();
                if (workCopy == checkCopy) {
                    times.push_back(durT);
                } else {
                    std::cerr << formatTime(true, true) << d0.first << " \t" << g0.first << " \t" << g1.first << " \tSort error\n";
                    g0.second.status = gapStruct::outOfOrder;
                    g1.second.dData[d0.first].time = ulMax;
                    g0.second.warnings = MAX_Warnings;
                    errorFunction(workCopy, checkCopy);
                    break;
                }
                workCopy = d1.second;
            }
            long dur(median(times));
            if (dur > 0)
                printIt(d0, dur, g0, g1, gTimes, tgs);
        }
    }
}

static void work(m_s_gs &gMap, m_s_ds &dMap, vul sizes) {
    for (auto d0 : dMap) { // each distro
        tg tgs;
        tgs.time = ulMax;
        tgs.gapper = "";
        for (auto d1 : d0.second.originals) {   // each sample
            std::cout << "\n\n" << formatTime(true, true) << " \tNew size: " << d1.first << " distro: " << d0.first << '\n';
            vul gTimes;
            for (auto &g0 : gMap) {     // each gapper
                if (g0.second.warnings < MAX_Warnings) {    // skip slow sequences
                    doSort(d0, d1, g0, gTimes, d1.first, tgs);
                } else {
                    std::cerr << formatTime(true, true) << " \tSkipping " << g0.first << " too slow\n";
                }
            }
            std::cout << formatTime(true, true) << " \tBest gapper for size of " << d1.first
            << " in distro " << d0.first << " is " << tgs.gapper << '\n';
            tgs.time = ulMax;
            tgs.gapper = "";
            if (WARN_Lagards)
                warnLagards(d0, gMap, gTimes);
        }
    }
}

void init() {
    std::cerr << formatTime(true, true) << " \tInitialization begins.\n";
    vul sizes({1234567, 123456789, 42});
    for (auto &size : sizes) {
        size = size > MAX_SampleSize ? MAX_SampleSize : size < MIN_SampleSize ? MIN_SampleSize : size;
    }
    std::sort(sizes.begin(), sizes.end());
    
    m_s_ds dMap;
    m_s_gs gMap;
    
    make_gMap(gMap, sizes);
    make_dMap(dMap, sizes);
    std::cerr << formatTime(true, true) << " \tInitialization ends.\n";
    
    if (FULL_Run)
        work(gMap, dMap, sizes);
    
    eoj(gMap, dMap, sizes);
}

void shell(vul &gaps, ul vSize) {
    ul gap(vSize);
    while (gap > 1) {
        gap >>= 1;
        gaps.push_back(gap);
    }
}

void frank(vul &gaps, ul vSize) {
    ul gap(vSize >> 1);
    while (gap) {
        gaps.push_back(gap | 1);
        gap >>= 1;
    }
}

void hibbard(vul &gaps, ul vSize) {
    ul gap(1);
    while (gap  < vSize) {
        gaps.push_back(gap);
        gap <<= 1;
        gap |= 1;
    }
}

void papernov(vul &gaps, ul vSize) {
    ul n(1);
    gaps.push_back(n);
    while (gaps.back() < vSize)
        gaps.push_back((2 << n++) + 1);
    gaps.pop_back();
}

bool is3smooth(ul n) {
    while (n % 3 == 0)
        n /= 3;
    while (n % 2 == 0)
        n /= 2;
    return n == 1;
}

void pratt(vul &gaps, ul vSize) {
    gaps.clear();
    for (ul n(1); n < vSize; n++)
        if (is3smooth(n))
            gaps.push_back(n);

}

void pratt_A(vul &gaps, ul vSize) {
    gaps.clear();
    for (ul n(1); n < vSize / 2; n += 4)
        if (is3smooth(n))
            gaps.push_back(n);
    
}

void knuth(vul &gaps, ul vSize) {
    ul k(1), lim(vSize / 3);
    do {
        k *= 3;
        gaps.push_back((k - 1) >> 1);
    } while (gaps.back() < lim);
    gaps.pop_back();
}

bool mySeq(ul a, ul b) {return a > b;}

ul pw2(ul e) {
    ul rv(1);
    while (e--) {
        rv *= 2;
    }
    return rv;
}

void sedgewick82(vul &gaps, ul vSize) {
    gaps.push_back(1);
    for (ul k(1); gaps.back() < vSize; k++) {
        gaps.push_back(pw2(k+1) + 3 * pw2(k-1) + 1);
    }
}

void sedgewick86(vul &gaps, ul vSize) {
     ul k(0);
    gaps.push_back(1);
    do {
        if (++k & 1) {
            gaps.push_back((2 << (k + 3)) - (6 * (1 << (k >> 1)) + 1));
        } else {
            gaps.push_back(9 * (2 << k) - (1 << k) + 1);
        }
    } while ( gaps.back() < vSize);
    gaps.pop_back();
}

void gonnet(vul &gaps, ul vSize) {
    ul k(vSize);
    while (k > 1) {
        k = (5 * k - 1) / 11;
        k = k > 1 ? k : 1;
        gaps.push_back(k);
    }
}

void tokuda(vul &gaps, ul vSize) {
    gaps.push_back(1);
    for (ul i(1); gaps.back() < vSize; i++) {
        double a(pow(2.25, static_cast<double>(i)));
        ul j((9.0 * a - 4.0) / 5.0);
        gaps.push_back(j | 1);
    }
}

void ciura(vul &gaps, ul vSize) {
    vul t {1, 4, 10, 23, 57, 132, 301, 701, 1750};
    gaps = t;
    while (gaps.back() < vSize / 3) {
        gaps.push_back(gaps.back() << 3);
    }
}

void a(vul &gaps, ul vSize) {
    const int lac(7); // ladd candidate
    const int ladd(vSize > lac ? lac : 1), sra_0(2), sra_1(3), sra_2(7), strt_0(2), strt_1(5), strt_2(5);
    gaps.push_back((vSize - (vSize >> strt_0) - (vSize >> strt_1) - (vSize >> strt_2)) | ladd);
    do {
        auto tmp((gaps.back() - (gaps.back() >> sra_0) - (gaps.back() >> sra_1) - (gaps.back() >> sra_2)) | ladd);
        tmp = tmp < gaps.back() ? tmp : tmp >> 1;
        gaps.push_back(tmp);
    } while(gaps.back() > ladd);
}

void b(vul &gaps, ul vSize) {
    const int lac(7); // ladd candidate
    const int ladd(vSize > lac ? lac : 1), sra_0(2), sra_1(3), sra_2(7), strt_0(1), strt_1(28), strt_2(28);
    gaps.push_back((vSize - (vSize >> strt_0) - (vSize >> strt_1) - (vSize >> strt_2)) | ladd);
    do {
        auto tmp((gaps.back() - (gaps.back() >> sra_0) - (gaps.back() >> sra_1) - (gaps.back() >> sra_2)) | ladd);
        tmp = tmp < gaps.back() ? tmp : tmp >> 1;
        gaps.push_back(tmp);
    } while(gaps.back() > ladd);
}

void c(vul &gaps, ul vSize) {
    const int lac(7); // ladd candidate
    const int ladd(vSize > lac ? lac : 1), sra_0(2), sra_1(3), sra_2(7), strt_0(1), strt_1(2), strt_2(3);
    gaps.push_back((vSize - (vSize >> strt_0) - (vSize >> strt_1) - (vSize >> strt_2)) | ladd);
    do {
        auto tmp((gaps.back() - (gaps.back() >> sra_0) - (gaps.back() >> sra_1) - (gaps.back() >> sra_2)) | ladd);
        tmp = tmp < gaps.back() ? tmp : tmp >> 1;
        gaps.push_back(tmp);
    } while(gaps.back() > ladd);
}

void d(vul &gaps, ul vSize) {
    const int lac(7); // ladd candidate
    const int ladd(vSize > lac ? lac : 1), sra_0(2), sra_1(3), sra_2(7), strt_0(2), strt_1(4), strt_2(6);
    gaps.push_back((vSize - (vSize >> strt_0) - (vSize >> strt_1) - (vSize >> strt_2)) | ladd);
    do {
        auto tmp((gaps.back() - (gaps.back() >> sra_0) - (gaps.back() >> sra_1) - (gaps.back() >> sra_2)) | ladd);
        tmp = tmp < gaps.back() ? tmp : tmp >> 1;
        gaps.push_back(tmp);
    } while(gaps.back() > ladd);
}

void e(vul &gaps, ul vSize) {
    const int lac(7); // ladd candidate
    const int ladd(vSize > lac ? lac : 1), sra_0(2), sra_1(3), sra_2(7), strt_0(2), strt_1(3), strt_2(4);
    gaps.push_back((vSize - (vSize >> strt_0) - (vSize >> strt_1) - (vSize >> strt_2)) | ladd);
    do {
        auto tmp((gaps.back() - (gaps.back() >> sra_0) - (gaps.back() >> sra_1) - (gaps.back() >> sra_2)) | ladd);
        tmp = tmp < gaps.back() ? tmp : tmp >> 1;
        gaps.push_back(tmp);
    } while(gaps.back() > ladd);
}

void shellSort(vul &v, vul &gaps) {
//    inspectGaps(gaps, v.size());
    for (auto gap : gaps) {
        for (auto iti(v.begin() + gap); iti != v.end(); iti++) {
            auto tmp(*iti);
            auto itj(iti);
            for (itj = iti; itj >= v.begin() + gap && *(itj - gap) > tmp; itj -= gap)
                *itj = *(itj - gap);
            *itj = tmp;
        }
    }
}

void getRandyBe(vul &v, ul n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::bernoulli_distribution dist(0.5);
    while (n--) {
        v.push_back(dist(rd));
    }
}

void getRandyBi(vul &v, ul n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::binomial_distribution<int> dist(1000, 0.5);
    while (n--) {
        v.push_back(dist(rd));
    }
}

void getRandyN(vul &v, ul n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(10.0, 5.0);
    while (n--) {
        v.push_back(dist(rd));
    }
}

void getRandyP(vul &v, ul n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::poisson_distribution<int> dist(1000.0);
    while (n--) {
        v.push_back(dist(rd));
    }
}

void getRandyU(vul &v, ul n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(iMin, iMax);
    while (n--) {
        v.push_back(dist(gen));
    }
}

void randomFill(ul n, vul &v, std::string distroName) {
    v.clear();
    if (distroName == "Normal") {
        getRandyN(v,n);
    } else if(distroName == "Poisson") {
        getRandyP(v,n);
    } else if(distroName == "Bernoulli") {
        getRandyBe(v,n);
    } else if(distroName == "Binomial") {
        getRandyBi(v,n);
    } else if(distroName == "Uniform") {
        getRandyU(v,n);
    } else if(distroName == "Uniform - Sorted") {
        getRandyU(v,n);
        std::sort(v.begin(), v.end());
    } else if(distroName == "Uniform - Sorted & Reversed") {
        getRandyU(v,n);
        std::sort(v.begin(), v.end());
        std::reverse(v.begin(), v.end());
    } else {
        std::cerr << "Unknown distribution requested (" << distroName << "). Using uniform." << std::endl;
        getRandyU(v,n);
    }
    v.shrink_to_fit();
}
