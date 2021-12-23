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
    
    gapStruct sedgewick;
    sedgewick.name = "Sedgewick 1982";
    sedgewick.gapFn = sedgewick1982;
    algorithms.emplace_back(sedgewick);
}

void setup() {
    std::vector<gapStruct> algorithms;
    makeAlgorithmElements(algorithms);
    
    int ssMin(8192), ssMax(8192 << 12), wdth(21);
    std::cout << "\nStart: " << ssMin << "  Max: " << ssMax << '\n';
    
    for (auto &a : algorithms)
        a.gapFn(a.gaps, ssMax);
    
    vi orginalCopy, workCopy, checkCopy;
    
    for (int sampleSize(ssMin - 1); sampleSize < ssMax; (sampleSize <<= 2) |= 3) {
        randomFill(sampleSize, orginalCopy);
        checkCopy = orginalCopy;
        std::sort(checkCopy.begin(), checkCopy.end());
        std::cout << '\n' << formatTime(true, true) << "    n: " << std::right
        << std::setw(24) << sampleSize << " ----------" << std::endl;
        for (auto &a : algorithms) {
            workCopy = orginalCopy;
            auto start = high_resolution_clock::now();
            shellsort(workCopy, a.gaps);
            auto stop = high_resolution_clock::now();
            long duration = duration_cast<microseconds>(stop - start).count();
            std::cout << formatTime(false, true) << std::right << std::setw(27)
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
    for (auto rd : v[0].runData) {
        fst << ',' << rd.sampleSize;
    }
    fst << '\n';
    for (auto s : v) {
        fst << s.name;
        for (auto rd : s.runData) {
            fst << ',' << rd.time;
        }
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
    gaps.clear();
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
    gaps.clear();
    gaps.push_back(1);
    while (gaps.back() < vSize)
        gaps.push_back((2 << n++) + 1);
    gaps.pop_back();
    std::reverse(gaps.begin(), gaps.end());
}

bool is3smooth(int n) {
    int f(n);
    while (f % 2 == 0)
        f /= 2;
    while (f % 3 == 0)
        f /= 3;
    return f == 1;
}

void pratt1971(vi &gaps, int vSize) {
    std::set<int>smooths;
    smooths.clear();
    smooths.insert(1);
    int s2(0), s3(0), s23(0);
    int lim(vSize >> 1);
    for (int n(1); s23 < lim; n++) {
        s2 = n << 1;
        s3 = s2 + n;
        s23 = s2 * s3;
        smooths.insert(s2);
        
        if (s3 < vSize)
            smooths.insert(s3);
        
        if (s23 < vSize)
            smooths.insert(s23);
    }
    
    gaps.clear();
    auto rits(smooths.rbegin());
    while (rits != smooths.rend()) {
        if (is3smooth(*rits))
            gaps.push_back(*rits);
        rits++;
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

bool mySeq(int a, int b) {return a < b;}

void sedgewick1982(vi &gaps, int vSize) {
    int k4(4), k2(1), k(1), lim(vSize >> 1);
    gaps.push_back(1);
    do {
        gaps.push_back(k4 + 3 * k2 + 1);
        k++;
        k4 *= 4;
        k2 *= 2;
    } while (gaps.back() < lim);
    gaps.pop_back();
    std::reverse(gaps.begin(), gaps.end());
    for (auto it(gaps.begin() + 1); it != gaps.end(); it++) {
        if (*it > *(it - 1)) {
            std::sort(gaps.begin(), gaps.end(), mySeq);
            break;
        }
    }
}

