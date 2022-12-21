//
//  core.cpp
//
//  Created by Phil Huffman on 12/21/21.
//

#include "core.hpp"

template<typename T>
T median(std::vector<T> a) {
    if (a.empty()) {
        std::cerr << formatTime(true, true) << " \tMedian error. Empty.";
        exit(3);
    }
    
    auto it(a.begin() + a.size() / 2);
    std::sort(a.begin(), a.end());
    if ((a.size() & 1) == 0) {
        return ((*it + *(it+1)) / 2);
    }
    return *it;
}

template<typename T>
T sum(std::vector<T> a) {
    if (a.empty())
        return -1;
    T sum(0);
    for (auto tmp : a) {
        sum += tmp;
    }
    return sum;
}

template<typename T>
T average(std::vector<T> a) {
    if (a.empty())
        return -1;
    return sum(a) / a.size();
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
    
    std::map<std::string, vi> mit;
    for (auto d0 : dMap) {
        mit[d0.first] = d0.second.originals.begin()->second.sample;
    }
    std::cerr << formatTime(true, true) << " \tWrite distros begins.\n";
    long maxPossibleLines(dMap.begin()->second.originals.begin()->first);   // size of the first (smallest) sample vector
    long maxDistroLines(MAX_DistroLines < maxPossibleLines ? MAX_DistroLines : maxPossibleLines);
    std::string fnBase(FN_Base);
    fnBase += formatTime(true, true);
    fnBase += "-distros.csv";
    std::fstream fst;
    fst.open(fnBase, std::ios::out);
    for(long indx(-1); indx < maxDistroLines; indx++) {
        for (auto vit : mit) {
            if (indx < 0) {
                fst << ',' << vit.first;     // header
            } else {
                fst << ',' << vit.second[indx];
            }
        }
        fst << '\n';
    }
    fst << std::endl;
    fst.close();
    std::cerr << formatTime(true, true) << " \tWrite distros ends.\n";
}

static void writeGaps(m_s_gs gMap) {
    std::cerr << formatTime(true, true) << " \tWrite gaps begins.\n";
    std::fstream gst;
    std::string fnBase(FN_Base + formatTime(true, true));
    fnBase += "-Gaps.csv";
    gst.open(fnBase, std::ios::out);
    gst << "Sequence,Size,First,Second,etc." << '\n';
    for (auto g0 : gMap) {
        for (auto g1 : g0.second.results) {
            gst << g0.first << ',' << g1.first << gaps2string(g1.second.gaps, ",") << '\n';
        }
    }
    gst << std::endl;
    gst.close();
    std::cerr << formatTime(true, true) << " \tWrite gaps ends.\n";
}

static void prepG1(m_s_gs &gMap, vul sizes) {
    for (auto &g0 : gMap) {
        for (auto dName : DISTRO_NAMES) {
            for (auto size : sizes) {
                g0.second.results[size].gaps.clear();
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
//    gMap["1959 Shell"].gapFn = shell;
    gMap["1960 Frank & Lazarus"].gapFn = frank;
    gMap["1963 Hibbard"].gapFn = hibbard;
    gMap["1965 Papernov & Stasevich"].gapFn = papernov;
//    gMap["1971 Pratt"].gapFn = pratt;
    gMap["1973 Knuth"].gapFn = knuth;
    gMap["1982 Sedgewick"].gapFn = sedgewick82;
    gMap["1986 Sedgewick"].gapFn = sedgewick86;
    gMap["1991 Gonnet & Baeza-Yates"].gapFn = gonnet;
    gMap["1992 Tokuda"].gapFn = tokuda;
//    gMap["2001 Ciura"].gapFn = ciura;
//    gMap["System"].gapFn = sys;
    gMap["2022 a"].gapFn = a;
    gMap["2022 b"].gapFn = b;
    gMap["2022 c"].gapFn = c;
    gMap["2022 d"].gapFn = d;
    gMap["2022 e"].gapFn = e;
    gMap["2022 f"].gapFn = e;
    gMap["2022 g"].gapFn = e;
    
    prepG1(gMap, sizes);
    
    for (auto &g0 : gMap) {
        g0.second.status = gs::ok;
        g0.second.warnings = 0;
        for (auto &g1 : g0.second.results) {
            g1.second.gaps.clear();
        }
    }
    
    makeGapSequences(gMap);
    writeGaps(gMap);
    std::cerr << formatTime(true, true) << " \tgMap built " << gMap.size() << " gappers running.\n";
}

static void make_dMap(m_s_ds &dMap, const vul sizes) {
    std::cerr // << formatTime(true, true)
    << "\t \t \t \t \t \tBuilding dMap for " << DISTRO_NAMES.size() << " distributions each with " << sizes.size() << " samples."
    << "\n\t \t \t \t \t \tBuilding large samples using complicated distributions will take time."
    << '\n';
    for (auto dName : DISTRO_NAMES) {
        for (auto size : sizes) {
            dMap[dName].originals[size].sample.clear();
            randomFill(size, dMap[dName].originals[size].sample, dName);
        }
        std::cerr << formatTime(true, true) << " \t" << dName << " distribution samples built.\n";
    }
    writeDistros(dMap);
    std::cerr << formatTime(true, true) << " \tdMap built.\n";
}

static void errorFunction(vi &wc, vi &cc) {
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
        fst << "Distro,Size,Sequence,Time\n";
        for (auto d0 : dMap) {
            for (auto d1 : d0.second.originals) {
                sort(d1.second.results.begin(), d1.second.results.end(), [](tg &lhs, tg &rhs) {
                    return lhs.time < rhs.time;
                });
                for (auto result : d1.second.results) {
                    fst << d0.first << ',' << d1.first << ',' << result.gapper << ',' << result.time << '\n';
                }
            }
        }
        fst << std::endl;
        fst.close();
    }
}

static void summerize(m_s_ds &dMap, m_s_gs gMap) {
    for (auto &d0 : dMap) {
        for (auto &d1 : d0.second.originals) {
            auto firstTime(true);
            for (auto result : d1.second.results) {
                if (firstTime) {
                    firstTime = false;
                    std::cout
                    << " \nn: "<< std::right << std::setw(11) << d1.first
                    << " \tDistribution: " << d0.first
                    << '\n';
                }
                std::cout
                << " Lapsed Time: " << std::right << std::setw(FORMATTED_MicroSecondLength) << formatMicroSeconds(result.time, 3)
                << std::right << std::setw(GAPPER_Length) << result.gapper
                << " Gaps: " << gaps2string(gMap[result.gapper].results[d1.first].gaps, " ")
                << '\n';
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
    summerize(dMap, gMap);
    listWinners(dMap);
}

static void doSort(std::pair<const std::string, distroStruct> &d0, std::pair<const unsigned long,originalSample> &d1,
                   std::pair<const std::string, gs> &g0, vul &gTimes, const ul size, vi checkCopy) {
    for (auto &g1 : g0.second.results ) {
        if (g1.first == d1.first) {
            vl times;
            times.clear();
            while (times.size() < MEDIAN_TrialSize) {
                auto workCopy(d1.second.sample);
                auto start = high_resolution_clock::now();
                shellSort(workCopy, g0.second.results[d1.first].gaps);
                auto stop = high_resolution_clock::now();
                auto durT = duration_cast<microseconds>(stop - start).count();
                if (workCopy == checkCopy) {
                    times.push_back(durT);
                } else {
                    std::cerr << formatTime(true, true) << d0.first << " \t" << g0.first << " \t" << g1.first << " \tSort error\n";
                    g0.second.status = gs::outOfOrder;
                    g0.second.warnings = iMax;
                    errorFunction(workCopy, checkCopy);
                    break;
                }
            }
            auto dur(times.size() == 1 ? times.front() : median(times));
            std::cout << formatTime(true, true)
            << std::right << std::setw(MICROSECOND_length) << dur << "µs "
            << std::right << std::setw(FORMATTED_MicroSecondLength) << formatMicroSeconds(dur, 3)
            << " \t" << g0.first
            << '\n';
            
            d1.second.results.push_back(tg(dur, g0.first));
            gTimes.push_back(dur);
        }
    }
}

static void checkForLagards(vtg results, m_s_gs &gMap, const vul &gTimes, ul warnLimit) {
    auto kudo(average(gTimes));
    auto limt(kudo + (kudo >> 2));
    for (auto result : results) {
        if (gMap[result.gapper].warnings < warnLimit) {
            if (result.time > limt) {
                gMap[result.gapper].warnings++;
                std::cerr << formatTime(true, true)
                << " \tWarned " << result.gapper << " (" << gMap[result.gapper].warnings << "/" << warnLimit << ")\n";
            } else if (result.time < kudo && gMap[result.gapper].warnings > 0) {
                gMap[result.gapper].warnings--;
                std::cerr << formatTime(true, true)
                << " \tBlessed " << result.gapper << " (" << gMap[result.gapper].warnings << "/" << warnLimit << ")\n";
            }
        }
    }
}

static void work(m_s_gs &gMap, m_s_ds &dMap) {
    auto warnLimit(SIZES.size() - 1 < MAX_Warnings ? MAX_Warnings : SIZES.size() - 1);
    for (auto &d0 : dMap) { // each distro
        for (auto &d1 : d0.second.originals) {   // each sample
            std::cout << "\n\n" << formatTime(true, true) << " \tn: " << d1.first << " distro: " << d0.first << '\n';
            auto checkCopy(d1.second.sample);
            std::sort(checkCopy.begin(), checkCopy.end());
            vul gTimes;
            for (auto &g0 : gMap) {     // each sequence
                if (g0.second.warnings < warnLimit) {    // skip slow sequences
                    doSort(d0, d1, g0, gTimes, d1.first, checkCopy);
                } else {
                    std::cerr << formatTime(false, true) << " \tSkipping " << g0.first << " (too slow)\n";
                }
            }
            if (WARN_Lagards) {
                checkForLagards(d1.second.results  , gMap, gTimes, warnLimit);
            }
            sort(d1.second.results.begin(), d1.second.results.end(), [](tg &lhs, tg &rhs) {
                return lhs.time < rhs.time;
            });
            std::cout << formatTime(true, true) << " \tBest sequence for size of " << d1.first
            << " with distro " << d0.first << " is " << d1.second.results.front().gapper << '\n';
        }
    }
    writeTimes(gMap, dMap);
    eoj(gMap, dMap);
}

void init() {
    for (int passes(0); passes < MAX_Passes; passes++) {
        auto t0(system_clock::now());
        std::cerr << formatTime(true, true) << " \tInitialization begins.\n";
        auto sizes(SIZES);
        for (auto &size : sizes) {
            size = size > MAX_SampleSize ? MAX_SampleSize : size < MIN_SampleSize ? MIN_SampleSize : size;
        }
        
        // sizes must be unique because it is used as an index in both dMap & gMap.
        std::sort(sizes.begin(), sizes.end());
        sizes.erase(std::unique(sizes.begin(), sizes.end()), sizes.end());  // ensure sizes are unique
     
        m_s_ds dMap;
        m_s_gs gMap;
        
        make_dMap(dMap, sizes);
        make_gMap(gMap, sizes);
        auto t1(system_clock::now());
        auto durT = duration_cast<microseconds>(t1 - t0).count();
        std::cerr << formatTime(true, true) << " \tInitialization complete. Required "
        << formatMicroSeconds(durT,1) << " seconds.\n";
        
        if (FULL_Run)
            work(gMap, dMap);
            auto t2(system_clock::now());
        auto durE = duration_cast<microseconds>(t2 - t1).count();
        std::cerr << formatTime(true, true) << " \tPass complete. Required "
        << formatMicroSeconds(durE,1) << ".\n";
 }
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

void sys(vul &gaps, ul vSize) {
    vul t({1});
    gaps = t;
}

void shellSort(std::vector<int>::iterator start, std::vector<int>::iterator stop, vul &gaps) {
    for (auto gap : gaps) {
        for (auto iti(start + gap); iti != stop; iti++) {
            auto tmp(*iti);
            auto itj(iti);
            for (itj = iti; itj >= start + gap && *(itj - gap) > tmp; itj -= gap)
                *itj = *(itj - gap);
            *itj = tmp;
        }
    }
}

void shellSort(vi &v, vul &gaps) {
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

void getRandyBe(vi &v, ul n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::bernoulli_distribution dist(0.5);
    while (n--) {
        v.push_back(dist(rd));
    }
}

void getRandyBi(vi &v, ul n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::binomial_distribution<int> dist(1000, 0.5);
    while (n--) {
        v.push_back(dist(rd));
    }
}

void getRandyE(vi &v, ul n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::exponential_distribution<double> dist(1000.0);
    while (n--) {
        v.push_back(dist(rd));
    }
}

void getRandyG(vi &v, ul n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::gamma_distribution<double> dist(1000.0, 5.0);
    while (n--) {
        v.push_back(dist(rd));
    }
}

void getRandyN(vi &v, ul n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(100.0, 5.0);
    while (n--) {
        v.push_back(dist(rd));
    }
}

void getRandyP(vi &v, ul n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::poisson_distribution<int> dist(100000.0);
    while (n--) {
        v.push_back(dist(rd));
    }
}

void getRandyU(vi &v, ul n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(iMin, iMax);
    while (n--) {
        v.push_back(dist(gen));
    }
}

void randomFill(ul n, vi &v, std::string distroName) {
    v.clear();
    if (distroName == "Normal") {
        getRandyN(v,n);
    } else if(distroName == "Poisson") {
        getRandyP(v,n);
    } else if(distroName == "Bernoulli") {
        getRandyBe(v,n);
    } else if(distroName == "Binomial") {
        getRandyBi(v,n);
    } else if(distroName == "Exponetial") {
        getRandyE(v,n);
    } else if(distroName == "Gamma") {
        getRandyG(v,n);
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
    int bits(2); // size of mask
    vi sri({3});
    vi srj({1,3,12});
    gap22(gaps, vSize, bits, sri, srj);
}

void b(vul &gaps, ul vSize) {
    int bits(2); // size of mask
    vi sri({3,4});
    vi srj({1,3,12});
    gap22(gaps, vSize, bits, sri, srj);
}

void c(vul &gaps, ul vSize) {
    int bits(2); // size of mask
    vi sri({3,5});
    vi srj({1,3,12});
    gap22(gaps, vSize, bits, sri, srj);
}

void d(vul &gaps, ul vSize) {
    int bits(2); // size of mask
    vi sri({3,6});
    vi srj({1,3,12});
    gap22(gaps, vSize, bits, sri, srj);
}

void e(vul &gaps, ul vSize) {
    int bits(2); // size of mask
    vi sri({3,7});
    vi srj({1,3,12});
    gap22(gaps, vSize, bits, sri, srj);
}

void f(vul &gaps, ul vSize) {
    int bits(2); // size of mask
    vi sri({3,8});
    vi srj({1,3,12});
    gap22(gaps, vSize, bits, sri, srj);
}

void g(vul &gaps, ul vSize) {
    int bits(2); // size of mask
    vi sri({3,9});
    vi srj({1,3,12});
    gap22(gaps, vSize, bits, sri, srj);
}
