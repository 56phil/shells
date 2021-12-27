//
//  core.cpp
//  core
//
//  Created by Phil Huffman on 12/21/21.
//

#include "core.hpp"

static void makeAlgorithmElements(std::vector<gapStruct> &algorithms) {
    algorithms.clear();
    
    gapStruct shell;
    shell.name = "Shell 1959";
    shell.gapFn = shell1959;
    shell.runData.clear();
    algorithms.emplace_back(shell);
    
    gapStruct frank;
    frank.name = "Frank & Lazarus 1960";
    frank.gapFn = frank1960;
    frank.runData.clear();
    algorithms.emplace_back(frank);
    
    gapStruct hibbard;
    hibbard.name = "Hibbard 1963";
    hibbard.gapFn = hibbard1963;
    algorithms.emplace_back(hibbard);
    
    gapStruct papernov;
    papernov.name = "Papernov & Stasevich 1965";
    papernov.gapFn = papernov1965;
    algorithms.emplace_back(papernov);
    
    gapStruct pratt;
    pratt.name = "Pratt 1971";
    pratt.gapFn = pratt1971;
    algorithms.emplace_back(pratt);
    
    gapStruct kunth;
    kunth.name = "Kunth 1973";
    kunth.gapFn = kunth1973;
    algorithms.emplace_back(kunth);
    
    gapStruct sedgewick82;
    sedgewick82.name = "Sedgewick 1982";
    sedgewick82.gapFn = sedgewick1982;
    algorithms.emplace_back(sedgewick82);
    
    gapStruct sedgewick85;
    sedgewick85.name = "Sedgewick 1985";
    sedgewick85.gapFn = sedgewick1985;
    algorithms.emplace_back(sedgewick85);
    
    gapStruct gonnet;
    gonnet.name = "Gonnet & Baeza-Yates 1991";
    gonnet.gapFn = gonnet1991;
    algorithms.emplace_back(gonnet);
    
    gapStruct tokuda;
    tokuda.name = "Tokuda 1992";
    tokuda.gapFn = tokuda1992;
    algorithms.emplace_back(tokuda);
    
    gapStruct empirical;
    empirical.name = "empirical 2001";
    empirical.gapFn = empirical2001;
    algorithms.emplace_back(empirical);
}

static void summerize(std::vector<gapStruct> &algorithms) {
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

static void work(std::vector<gapStruct> &algorithms) {
    int ssMin(512), ssMax(1028 << 16), wdth(17);
    std::cout << "\nStart: " << ssMin << "  Max: " << ssMax << '\n';
    
    vi orginalCopy, workCopy, checkCopy;
    
    for (int sampleSize(ssMin - 1); sampleSize < ssMax; (sampleSize <<= 1) |= 1) {
        randomFill(sampleSize, orginalCopy);
        checkCopy = orginalCopy;
        std::sort(checkCopy.begin(), checkCopy.end());
        std::cout << '\n' << formatTime(true, true) << "        n: " << std::right
        << std::setw(20) << sampleSize << " ----------" << std::endl;
        for (auto &a : algorithms) {
            a.gaps.clear();
            a.gapFn(a.gaps, sampleSize);
            workCopy = orginalCopy;
            auto start = high_resolution_clock::now();
            shellsort(workCopy, a.gaps);
            auto stop = high_resolution_clock::now();
            long duration = duration_cast<microseconds>(stop - start).count();
            std::cout << formatTime(false, true) << std::right << std::setw(31)
            << a.name << ": " <<std::setw(wdth) <<std::right << duration << " µs"
            << convertMicroSeconds(duration) << "\n";
            sortMetrics tm;
            tm.status = verify(workCopy, checkCopy) ? tm.ok : tm.outOfOrder;
            tm.time = duration;
            tm.sampleSize = sampleSize;
            if (tm.status == tm.outOfOrder)
                errorFunction(workCopy, checkCopy);
            a.runData.emplace_back(tm);
        } // function loop
        makeFile(algorithms);
    } // sample size loop
    
    summerize(algorithms);
}

void setup() {
    std::vector<gapStruct> algorithms;
    makeAlgorithmElements(algorithms);
    
    work(algorithms);
}

void errorFunction(vi wc, vi cc) {
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

void makeFile(std::vector<gapStruct> v) {
    std::fstream fst;
    fst.open("/Users/prh/Keepers/code/cpp/sorts/ShellSorts/list.csv", std::ios::out);
    fst << "Algorithm";
    for (auto rd : v[0].runData)
        fst << ',' << rd.sampleSize;
    fst << '\n';
    for (auto s : v) {
        fst << s.name;
        for (auto rd : s.runData)
            fst << ',' << rd.time;
        fst << '\n';
    }
    fst.close();
}

void shell1959(vi &gaps, int vSize) {
    int wSize(static_cast<int>(vSize));
    while (wSize > 1) {
        wSize >>= 1;
        gaps.push_back(static_cast<int>(wSize));
    }
}

void frank1960(vi &gaps, int vSize) {
    int wSize(static_cast<int>(vSize));
    gaps.clear();
    while (wSize) {
        wSize >>= 1;
        gaps.push_back(wSize + 1);
    }
}

void hibbard1963(vi &gaps, int vSize) {
    int wSize(1);
    while (wSize  < vSize) {
        gaps.push_back(wSize);
        wSize <<= 1;
        wSize |= 1;
    }
    std::reverse(gaps.begin(), gaps.end());
}

void papernov1965(vi &gaps, int vSize) {
    int n(1);
    gaps.push_back(n);
    while (gaps.back() < vSize)
        gaps.push_back((2 << n++) + 1);
    gaps.pop_back();
    std::reverse(gaps.begin(), gaps.end());
}

bool is3smooth(int n) {
    while (n % 2 == 0)
        n /= 2;
    while (n % 3 == 0)
        n /= 3;
    return n == 1;
}

void pratt1971(vi &gaps, int vSize) {
    std::set<int, std::greater<>>smooths;
    smooths.clear();
    smooths.insert(1);
    for (int n(1); n < vSize/3; n++)
        if (is3smooth(n))
            smooths.insert(n);
    
    auto its(smooths.begin());
    while (its != smooths.end()) {
            gaps.push_back(*its++);
    }
}

void kunth1973(vi &gaps, int vSize) {
    int k(1), lim(vSize / 3);
    do {
        k *= 3;
        gaps.push_back((k - 1) >> 1);
    } while (gaps.back() < lim);
    gaps.pop_back();
    std::reverse(gaps.begin(), gaps.end());
}

bool mySeq(int a, int b) {return a > b;}

void sedgewick1982(vi &gaps, int vSize) {
    int k4(4), k2(1), gap(0), lim(vSize - (vSize >> 3));
    std::set<int, std::greater<>> gSet;
    gSet.insert(1);
    do {
        gap = k4 + 3 * k2 + 1;
        gSet.insert(gap);
        k4 *= 4;
        k2 *= 2;
    } while (gap < lim);
    gaps.clear();
    for (auto g : gSet)
        gaps.push_back(g);
    if (gaps.front() < gaps.back())
        std::reverse(gaps.begin(), gaps.end());
}

void sedgewick1985(vi &gaps, int vSize) {
    int k(0);
    gaps.push_back(1);
    do {
        if (++k & 1) {
            gaps.push_back((2 << (k + 3)) - (6 * (1 << (k >> 1)) + 1));
        } else {
            gaps.push_back(9 * (2 << k) - (1 << k) + 1);
        }
    } while ( gaps.back() < vSize);
    gaps.pop_back();
    std::reverse(gaps.begin(), gaps.end());
}

void gonnet1991(vi &gaps, int vSize) {
    int k(vSize);
    while (k > 1) {
        k = std::max((5 * k - 1) / 11, 1);
        gaps.push_back(k);
    }
}

void tokuda1992(vi &gaps, int vSize) {
    gaps.push_back(1);
    for (int i(1); gaps.back() < vSize; i++) {
        double a(pow(2.25, static_cast<double>(i)));
        int j((9.0 * a - 4.0) / 5.0);
        gaps.push_back(j | 1);
    }
    std::reverse(gaps.begin(), gaps.end());
}

void empirical2001(vi &gaps, int vSize) {
    vi t {701, 301, 132, 57, 23, 10, 4, 1};
    gaps = t;
    int nextGap((gaps.front() << 2) | 1);
    while (nextGap < vSize) {
        gaps.insert(gaps.begin(), nextGap);
        nextGap = (nextGap << 3) | 1;
        nextGap = nextGap % 5 ? nextGap : nextGap - 4;
    }
}
