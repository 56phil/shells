//
//  core.cpp
//  core
//
//  Created by Phil Huffman on 12/21/21.
//

#include "core.hpp"

static void makeAlgorithmElements(vgs &algorithms) {
    algorithms.clear();
    
    gapStruct shell;
    shell.name = "Shell 1959";
    shell.gapFn = shell1959;
    shell.status = gapStruct::ok;
    shell.runData.clear();
    algorithms.emplace_back(shell);
    
    gapStruct frank;
    frank.name = "Frank & Lazarus 1960";
    frank.gapFn = frank1960;
    frank.runData.clear();
    frank.status = gapStruct::ok;
    algorithms.emplace_back(frank);
    
    gapStruct hibbard;
    hibbard.name = "Hibbard 1963";
    hibbard.gapFn = hibbard1963;
    hibbard.status = gapStruct::ok;
    algorithms.emplace_back(hibbard);
    
    gapStruct papernov;
    papernov.name = "Papernov & Stasevich 1965";
    papernov.gapFn = papernov1965;
    papernov.status = gapStruct::ok;
    algorithms.emplace_back(papernov);
    
    gapStruct pratt;
    pratt.name = "Pratt 1971";
    pratt.gapFn = pratt1971;
    pratt.status = gapStruct::ok;
    algorithms.emplace_back(pratt);
    
    gapStruct kunth;
    kunth.name = "Kunth 1973";
    kunth.gapFn = kunth1973;
    kunth.status = gapStruct::ok;
    algorithms.emplace_back(kunth);
    
    gapStruct sedgewick82;
    sedgewick82.name = "Sedgewick 1982";
    sedgewick82.gapFn = sedgewick1982;
    sedgewick82.status = gapStruct::ok;
    algorithms.emplace_back(sedgewick82);
    
    gapStruct sedgewick86;
    sedgewick86.name = "Sedgewick 1986";
    sedgewick86.gapFn = sedgewick1986;
    sedgewick86.status = gapStruct::ok;
    algorithms.emplace_back(sedgewick86);
    
    gapStruct gonnet;
    gonnet.name = "Gonnet & Baeza-Yates 1991";
    gonnet.gapFn = gonnet1991;
    gonnet.status = gapStruct::ok;
    algorithms.emplace_back(gonnet);
    
    gapStruct tokuda;
    tokuda.name = "Tokuda 1992";
    tokuda.gapFn = tokuda1992;
    tokuda.status = gapStruct::ok;
    algorithms.emplace_back(tokuda);
    
    gapStruct empirical;
    empirical.name = "empirical 2001";
    empirical.gapFn = empirical2001;
    empirical.status = gapStruct::ok;
    algorithms.emplace_back(empirical);
    
    gapStruct huffman;
    huffman.name = "Huffman 2022";
    huffman.gapFn = huffman2022;
    huffman.status = gapStruct::ok;
    algorithms.emplace_back(huffman);
}

static void errorFunction(vi &wc, vi &cc) {
    int n(0), w(28);
    
    auto itw(wc.begin());
    auto itc(cc.begin());
    
    std::cout << std::right << std::setw(4) << "n"
    << std::right << std::setw(w) << "original"
    << std::right << std::setw(w) << "result" << '\n';
    while (itw != wc.end() && itc != cc.end() && n < 33)
        std::cout << std::right << std::setw(4) << ++n
        << std::right << std::setw(w) << *itc++
        << std::right << std::setw(w) << *itw++ << '\n';
}

static void makeFile(vgs &v) {
    std::fstream fst, gst;
    std::string fName("/Users/prh/Keepers/code/xCode/shells/");
    fName += formatTime(true, true);
    std::string fNameResults = fName + "-Results.csv";
    std::string fNameGaps = fName + "-Gaps.csv";
    fst.open(fNameResults, std::ios::out);
    gst.open(fNameGaps, std::ios::out);
    fst << "Algorithm";
    for (auto rd : v[0].runData)
        fst << ',' << rd.sampleSize;
    fst << '\n';
    for (auto s : v) {
        fst << s.name;
        gst << s.name;
        for (auto rd : s.runData)
            fst << ',' << rd.time;
        fst << '\n';
        for (auto n : s.gaps) {
            gst << ',';
            gst << n;
        }
        gst << '\n';
    }
    fst.close();
    gst.close();
}

static void summerize(vgs &algorithms) {
    for (auto a : algorithms) {
        std::cout << '\n' << a.name;
        long n(a.gaps.size());
        for (auto itr(a.gaps.rbegin()); n--; itr++)
            std::cout << "  " << *itr;
        std::cout << '\n';
        for (auto d : a.runData)
            std::cout << std::right << std::setw(12) << d.sampleSize
            << std::right << std::setw(14) << d.time
            << std::setw(18)<< convertMicroSeconds(d.time) << '\n';
        std::cout << '\n';
    }
    std::cout << std::endl;
}

static void cullSlowerAlgorithms(vgs &algorithms, ull averageFuncTime) {
    /*
     Deactivates algorithms that had subpar sort times
     */
    ull lim(averageFuncTime + (averageFuncTime >> 1));
    
    for (auto &a : algorithms) {
        if (a.status == gapStruct::ok && a.runData.back().time > lim) {
            a.status = gapStruct::deactivated;
        }
    }
}

static void getGaps(vgs &algorithms, ull sampleSize) {
    for (auto &a : algorithms) {
        if (a.status == gapStruct::ok) {
            a.gaps.clear();
            a.gapFn(a.gaps, sampleSize);
            if (a.gaps.front() == 1) {
                std::reverse(a.gaps.begin(), a.gaps.end());
            }
            while (a.gaps.front() >= sampleSize) {
                a.gaps.erase(a.gaps.begin());
            }
            a.gaps.shrink_to_fit();
        }
    }
}

static void traverseAlgorithmVector(ull &activeFuncCount, vgs &algorithms, vi &checkCopy, const vi &orginalCopy, ull sampleSize, ull &totalFuncTime, int wdth, vi &workCopy) {
    for (auto &a : algorithms) {
        if (a.status == gapStruct::ok) {
            workCopy = orginalCopy;
            auto start = high_resolution_clock::now();
            shellsort(workCopy, a.gaps);
            auto stop = high_resolution_clock::now();
            long duration = duration_cast<microseconds>(stop - start).count();
            std::cout << formatTime(false, true) << std::right << std::setw(31)
            << a.name << ": " <<std::setw(wdth) <<std::right << duration << " µs"
            << convertMicroSeconds(duration) << "\n";
            sortMetrics sortMetrics;
            a.status = verify(workCopy, checkCopy) ? gapStruct::ok : gapStruct::outOfOrder;
            sortMetrics.time = duration;
            sortMetrics.sampleSize = sampleSize;
            if (a.gapStruct::status == a.gapStruct::outOfOrder)
                errorFunction(workCopy, checkCopy);
            a.runData.emplace_back(sortMetrics);
            totalFuncTime += duration;
            activeFuncCount++;
        }
    }
}

static void runActiveAlgorithms(vgs &algorithms, vi &checkCopy, const vi &orginalCopy, ull sampleSize,  int wdth, vi &workCopy) {
    ull totalFuncTime(0), activeFuncCount(0);
    getGaps(algorithms, sampleSize);
    traverseAlgorithmVector(activeFuncCount, algorithms, checkCopy, orginalCopy, sampleSize, totalFuncTime, wdth, workCopy);
    cullSlowerAlgorithms(algorithms, totalFuncTime / activeFuncCount);
    totalFuncTime = 0;
    activeFuncCount = 0;
}

static void eoj(vgs &algorithms) {
    summerize(algorithms);
    makeFile(algorithms);
}

static void prep4size(vi &checkCopy, vi &orginalCopy, ull sampleSize) {
    randomFill(sampleSize, orginalCopy);
    orginalCopy.shrink_to_fit();
    checkCopy = orginalCopy;
    std::sort(checkCopy.begin(), checkCopy.end());
    std::cout << '\n' << formatTime(true, true) << " n: " << sampleSize << std::endl;
}

static void work(vgs &algorithms) {
    int wdth(11);
    ull  ssMin(100000), ssMax(1000000000000);
    std::cout << "\nStart: " << ssMin << "  Max: " << ssMax << '\n';
    
    vi orginalCopy, workCopy, checkCopy;
    
    for (ull sampleSize(ssMin); sampleSize < ssMax; sampleSize *= 10) {
        prep4size(checkCopy, orginalCopy, sampleSize);
        runActiveAlgorithms(algorithms, checkCopy, orginalCopy, sampleSize, wdth, workCopy);
    }
}

void setup() {
    vgs algorithms;
    makeAlgorithmElements(algorithms);
    for (int i(0); i < 3; i++) {
        work(algorithms);
    }
    eoj(algorithms);
}

void shell1959(vull &gaps, ull vSize) {
    ull gap(vSize);
    while (gap > 1) {
        gap >>= 1;
        gaps.push_back(gap);
    }
}

void frank1960(vull &gaps, ull vSize) {
    ull gap(vSize >> 1);
    while (gap) {
        gaps.push_back(gap | 1);
        gap >>= 1;
    }
}

void hibbard1963(vull &gaps, ull vSize) {
    ull gap(1);
    while (gap  < vSize) {
        gaps.push_back(gap);
        gap <<= 1;
        gap |= 1;
    }
}

void papernov1965(vull &gaps, ull vSize) {
    ull n(1);
    gaps.push_back(n);
    while (gaps.back() < vSize)
        gaps.push_back((2 << n++) + 1);
    gaps.pop_back();
}

bool is3smooth(ull n) {
    while (n % 2 == 0)
        n /= 2;
    while (n % 3 == 0)
        n /= 3;
    return n == 1;
}

bool is5smooth(ull n) {
    while (n % 5 == 0)
        n /= 5;
    return is3smooth(n);
}

void pratt1971(vull &gaps, ull vSize) {
    gaps.clear();
    for (ull n(1); n < vSize; n++)
        if (is3smooth(n))
            gaps.push_back(n);
    
}

void kunth1973(vull &gaps, ull vSize) {
    ull k(1), lim(vSize / 3);
    do {
        k *= 3;
        gaps.push_back((k - 1) >> 1);
    } while (gaps.back() < lim);
    gaps.pop_back();
}

bool mySeq(ull a, ull b) {return a > b;}

void sedgewick1982(vull &gaps, ull vSize) {
    ull k4(4), k2(1), gap(0), lim(vSize - (vSize / 8));
    std::set<ull, std::greater<>> gSet;
    gSet.insert(1);
    do {
        gap = k4 + 3 * k2 + 1;
        gSet.insert(gap);
        k4 *= 4;
        k2 *= 2;
    } while (gap < lim);
    gaps.clear();
    for (ull g : gSet)
        gaps.push_back(g);
}

void sedgewick1986(vull &gaps, ull vSize) {
     ull k(0);
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

void gonnet1991(vull &gaps, ull vSize) {
    ull k(vSize);
    while (k > 1) {
        k = (5 * k - 1) / 11;
        k = k > 1 ? k : 1;
        gaps.push_back(k);
    }
}

void tokuda1992(vull &gaps, ull vSize) {
    gaps.push_back(1);
    for (ull i(1); gaps.back() < vSize; i++) {
        double a(pow(2.25, static_cast<double>(i)));
        ull j((9.0 * a - 4.0) / 5.0);
        gaps.push_back(j | 1);
    }
}

void empirical2001(vull &gaps, ull vSize) {
    vull t {701, 301, 132, 57, 23, 10, 4, 1};
    gaps = t;
    ull nextGap((gaps.front() << 2) | 1);
    while (nextGap < vSize) {
        gaps.insert(gaps.begin(), nextGap);
        nextGap = (nextGap << 3) | 1;
        nextGap = nextGap % 5 ? nextGap : nextGap - 4;
    }
}

void huffman2022(vull &gaps, ull vSize) {
    gaps.push_back(1);
    ull lim(vSize / 3 + (vSize >> 1));
    while (gaps.back() < lim) {
        ull k(2), gap(gaps.back() << 1), t(gaps.back());
        while (t > 0) {
            t = gaps.back() >> k;
            gap += t;
            k += 1;
        }
        gaps.push_back(gap);
    }
}
