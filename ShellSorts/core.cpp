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

static void writeDistros(m_s_ds dMap) {
    //    for (auto &d0 : dMap) {
    if (dMap.empty()) {
        std::cerr << formatTime(true, true) << " \tdMap error. Empty.";
        exit(2);
    }
    
    
    long maxPossibleLines(dMap.begin()->second.originals.begin()->first);   // size of the first (smallest) sample vector
    long maxDistroLines(MAX_DistroLines < maxPossibleLines ? MAX_DistroLines : maxPossibleLines);
    std::string fnBase("/Users/prh/Keepers/code/xCode/shells/results/");
    fnBase += formatTime(true, true);
    fnBase += "-distros.csv";
    std::fstream fst;
    fst.open(fnBase, std::ios::out);
    for(long indx(-1); indx < maxDistroLines; indx++) {
        for (auto d0 : dMap) {
            if (indx < 0) {
                fst << ',' << d0.first;
            } else {
                fst << ',' << d0.second.originals.begin()->second.sample.begin()[indx];
            }
        }
        fst << '\n';
    }
    fst << std::endl;
    fst.close();
}

static void writeGaps(m_s_gs gMap) {
    std::vector<std::string> alsoRans;
    std::fstream gst;
    std::string fnBase(FN_Base + formatTime(true, true));
    fnBase += "-Gaps.csv";
    gst.open(fnBase, std::ios::out);
    gst << "Gapper,Size,First,Second,etc." << '\n';
    for (auto g0 : gMap) {
        for (auto g1 : g0.second.results) {
            gst << g0.first << ',' << g1.first << gaps2string(g1.second.gaps, ",") << '\n';
        }
    }
    gst << std::endl;
    gst.close();
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
//    gMap["Shell 1959"].gapFn = shell;
//    gMap["Frank & Lazarus 1960"].gapFn = frank;
//    gMap["Hibbard 1963"].gapFn = hibbard;
//    gMap["Papernov & Stasevich 1965"].gapFn = papernov;
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
    writeGaps(gMap);
    std::cerr << formatTime(true, true) << " \tgMap constructed.\n";
}

static void make_dMap(m_s_ds &dMap, const vul sizes) {
    for (auto dName : DISTRO_NAMES) {
        for (auto size : sizes) {
            dMap[dName].originals[size].results.push_back(tg(ULONG_MAX, ""));
            dMap[dName].originals[size].sample.clear();
            randomFill(size, dMap[dName].originals[size].sample, dName);
        }
        std::cerr << formatTime(true, true) << " \tConstructed " << dName << ".\n";
    }
    writeDistros(dMap);
    std::cerr << formatTime(true, true) << " \tdMap build completed.\n";
}

static void errorFunction(vl &wc, vl &cc) {
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
        fst << "Distro,Gapper,Size,Time\n";
        for (auto d0 : dMap) {
            for (auto d1 : d0.second.originals) {
                for (auto g0 : gMap) {
                    for (auto g1 : g0.second.results) {
                        for (auto t : g1.second.dData) {
                            std::stringstream sst;
//                            sst << d0.first << ',' << d1.first << ',' << g0.first << ',' << t.second.time << '\n';
//                            std::string dmy(sst.str());
                            fst << d0.first << ',' << g0.first << ',' << d1.first << ',' << t.second.time << '\n';
                        }
                    }
                }
            }
        }
        fst << std::endl;
        fst.close();
    }
}

static void summerize(m_s_gs gMap) {
    for (auto g0 : gMap) {  // gapper
        for (auto g1 : g0.second.results) {     // sample size
            std::cout << std::right << std::setw(GAPPER_Length) << "\nGapper: " << g0.first << " \tn: "<< std::right << std::setw(11) << g1.first
            << "\nGaps: " << gaps2string(g1.second.gaps, " ") << '\n';
            for (auto g2 : g1.second.dData) {   // distro
                std::cout << formatTime(true, true)
                << " Distribution: " << std::setw(27) << g2.first
                << " \tLapsed Time: " << std::right << std::setw(FORMATTED_MicroSecondLength) << formatMicroSeconds(g2.second.time) << '\n';
            }
        }
        std::cout << '\n';
    }
}

void listWinners(m_s_ds &dMap) {
    if (FULL_Run) {
        for (auto &d0 : dMap) {
            for (auto &d1 : d0.second.originals) {
                sort(d1.second.results.begin(), d1.second.results.end(), [](tg &lhs, tg &rhs) {
                    return lhs.time < rhs.time;
                });
                std::cout << "Quickest gapper for distro " << d0.first << " with a size of " << d1.first << " is " << d1.second.results.front().gapper << '\n';
            }
            std::cout << '\n';
        }
    }
}

static void eoj(m_s_gs gMap, m_s_ds dMap) {
    summerize(gMap);
    listWinners(dMap);
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

static void printIt(const std::pair<const std::string, distroStruct> &d0, long dur, std::pair<const std::string,gapStruct> &g0,
                    std::pair<const unsigned long, runData> &g1, vul &gTimes) {
    g1.second.dData[d0.first].time = dur;   //???
    gTimes.push_back(dur);
    std::cout << formatTime(false, true)
    << std::right << std::setw(GAPPER_Length) << g0.first
    << std::right << std::setw(MICROSECOND_length) << dur << "µs "
    << std::right << std::setw(FORMATTED_MicroSecondLength) << formatMicroSeconds(dur)
    << '\n';
}

static void doSort(std::pair<const std::string, distroStruct> &d0, std::pair<const unsigned long,originalSample> &d1,
                   std::pair<const std::string, gapStruct> &g0, vul &gTimes, const ul size) {
    for (auto &g1 : g0.second.results ) {
        if (g1.first == d1.first) {
            auto checkCopy(d1.second.sample);
            std::sort(checkCopy.begin(), checkCopy.end());
            vl times(MEDIAN_TrialSize);
            times.clear();
            auto workCopy(d1.second.sample);
            while (times.size() < MEDIAN_TrialSize) {
                auto start = high_resolution_clock::now();
                shellSort(workCopy, g0.second.results[d1.first].gaps);
                auto stop = high_resolution_clock::now();
                auto durT = duration_cast<microseconds>(stop - start).count();
                if (workCopy == checkCopy) {
                    times.push_back(durT);
                    d1.second.results.push_back(tg(durT, g0.first));
                } else {
                    std::cerr << formatTime(true, true) << d0.first << " \t" << g0.first << " \t" << g1.first << " \tSort error\n";
                    g0.second.status = gapStruct::outOfOrder;
                    g1.second.dData[d0.first].time = UL_MAX;
                    g0.second.warnings = MAX_Warnings;
                    errorFunction(workCopy, checkCopy);
                    break;
                }
                workCopy = d1.second.sample;
            }
            long dur(median(times));
            if (dur > 0) {
                printIt(d0, dur, g0, g1, gTimes);
                d1.second.results.push_back(tg(dur, g0.first));
            }
        }
    }
}

static void work(m_s_gs &gMap, m_s_ds &dMap) {
    for (auto &d0 : dMap) { // each distro
        for (auto &d1 : d0.second.originals) {   // each sample
            std::cout << "\n\n" << formatTime(true, true) << " \tNew size: " << d1.first << " distro: " << d0.first << '\n';
            vul gTimes;
            for (auto &g0 : gMap) {     // each gapper
                if (g0.second.warnings < MAX_Warnings) {    // skip slow sequences
                    doSort(d0, d1, g0, gTimes, d1.first);
                } else {
                    std::cerr << formatTime(true, true) << " \tSkipping " << g0.first << " too slow\n";
                }
            }
            sort(d1.second.results.begin(), d1.second.results.end(), [](tg &lhs, tg &rhs) {
                return lhs.time < rhs.time;
             });
            std::cout << formatTime(true, true) << " \tBest gapper for size of " << d1.first
            << " with distro " << d0.first << " is " << d1.second.results.front().gapper << '\n';
            if (WARN_Lagards)
                warnLagards(d0, gMap, gTimes);
        }
    }
    writeTimes(gMap, dMap);
    eoj(gMap, dMap);
}

void init() {
    std::cerr << formatTime(true, true) << " \tInitialization begins.\n";
    vul sizes({250000, 2500000, 25000000000});
    for (auto &size : sizes) {
        size = size > MAX_SampleSize ? MAX_SampleSize : size < MIN_SampleSize ? MIN_SampleSize : size;
    }
    std::sort(sizes.begin(), sizes.end());
    
    m_s_ds dMap;
    m_s_gs gMap;
    
    make_gMap(gMap, sizes);
    make_dMap(dMap, sizes);
    std::cerr << formatTime(true, true) << " \tInitialization complete.\n";
    
    if (FULL_Run) 
        work(gMap, dMap);
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

void shellSort(vl &v, vul &gaps) {
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

void getRandyBe(vl &v, ul n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::bernoulli_distribution dist(0.5);
    while (n--) {
        v.push_back(dist(rd));
    }
}

void getRandyBi(vl &v, ul n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::binomial_distribution<int> dist(1000, 0.5);
    while (n--) {
        v.push_back(dist(rd));
    }
}

void getRandyN(vl &v, ul n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(100.0, 5.0);
    while (n--) {
        v.push_back(dist(rd));
    }
}

void getRandyP(vl &v, ul n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::poisson_distribution<int> dist(1000.0);
    while (n--) {
        v.push_back(dist(rd));
    }
}

void getRandyU(vl &v, ul n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(iMin, iMax);
    while (n--) {
        v.push_back(dist(gen));
    }
}

void randomFill(ul n, vl &v, std::string distroName) {
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

static void gap22(vul &gaps, ul vSize, int bits, vi sri, vi srj) {
    bits = bits < 1 ? 1 : bits > 6 ? 6 : bits;
    auto ladd(1);
    while (--bits > 0) {
        ladd <<= 1;
        ladd |= 1;
    }
    auto tmp(vSize);
    for (auto i : sri) {
        tmp -= vSize >> i;
    }
    gaps.push_back(tmp < vSize ? tmp > 0 ? tmp | 1 : vSize >> 2 : vSize >> 1);
    while (gaps.back() > ladd) {
        long tmp(gaps.back());
        for (auto j : srj) {
            long t0(gaps.back() >> j);
            if (t0 <= ladd)
                break;
            tmp -= gaps.back() >> j;
        }
        if (tmp >= gaps.back())
            break;
        gaps.push_back(tmp | ladd);
    }
    if (gaps.back() != 1)
        gaps.push_back(1);
}

void a(vul &gaps, ul vSize) {
    int bits(1); // ladd candidate
    vi sri({1,2,3});
    vi srj({1,2});
    gap22(gaps, vSize, bits, sri, srj);
}

void b(vul &gaps, ul vSize) {
    int bits(2); // ladd candidate
    vi sri({1,2,3});
    vi srj({1,2});
    gap22(gaps, vSize, bits, sri, srj);
}

void c(vul &gaps, ul vSize) {
    int bits(3); // ladd candidate
    vi sri({1,2,3});
    vi srj({1,2});
    gap22(gaps, vSize, bits, sri, srj);
}

void d(vul &gaps, ul vSize) {
    int bits(4); // ladd candidate
    vi sri({1,2,3});
    vi srj({1,2});
    gap22(gaps, vSize, bits, sri, srj);
}

void e(vul &gaps, ul vSize) {
    int bits(4); // ladd candidate
    vi sri({1,2,3});
    vi srj({1,3,4});
    gap22(gaps, vSize, bits, sri, srj);
}
