//
//  core.cpp
//  core
//
//  Created by Phil Huffman on 12/21/21.
//

#include "core.hpp"
static void fillDistros(std::vector<std::string> &distros) {
    distros.push_back("Bernoulli");
//    distros.push_back("Binomial");
//    distros.push_back("Chi Squared");
//    distros.push_back("Cauchy");
//    distros.push_back("Lognormal");
    distros.push_back("Normal");
//    distros.push_back("Student T");
    distros.push_back("Uniform");
    std::sort(distros.begin(), distros.end());
}

static void makeAlgorithmElements(vgs &algorithms) {
    algorithms.clear();
    
    gapStruct temp;
    temp.name = "Shell 1959";
    temp.gapFn = shell1959;
    temp.status = gapStruct::ok;
    temp.runData.clear();
    algorithms.emplace_back(temp);
    
    temp.name = "Frank & Lazarus 1960";
    temp.gapFn = frank1960;
    temp.runData.clear();
    temp.status = gapStruct::ok;
    algorithms.emplace_back(temp);
    
    temp.name = "Hibbard 1963";
    temp.gapFn = hibbard1963;
    temp.status = gapStruct::ok;
    algorithms.emplace_back(temp);
    
    temp.name = "Papernov & Stasevich 1965";
    temp.gapFn = papernov1965;
    temp.status = gapStruct::ok;
    algorithms.emplace_back(temp);
    
    temp.name = "Pratt 1971";
    temp.gapFn = pratt1971;
    temp.status = gapStruct::ok;
//    algorithms.emplace_back(temp);
    
    temp.name = "Knuth 1973";
    temp.gapFn = kunth1973;
    temp.status = gapStruct::ok;
    algorithms.emplace_back(temp);
    
    temp.name = "Sedgewick 1982";
    temp.gapFn = sedgewick1982;
    temp.status = gapStruct::ok;
    algorithms.emplace_back(temp);
    
    temp.name = "Sedgewick 1986";
    temp.gapFn = sedgewick1986;
    temp.status = gapStruct::ok;
    algorithms.emplace_back(temp);
    
    temp.name = "Gonnet & Baeza-Yates 1991";
    temp.gapFn = gonnet1991;
    temp.status = gapStruct::ok;
    algorithms.emplace_back(temp);
    
    temp.name = "Tokuda 1992";
    temp.gapFn = tokuda1992;
    temp.status = gapStruct::ok;
    algorithms.emplace_back(temp);
    
    temp.name = "empirical 2001";
    temp.gapFn = empirical2001;
    temp.status = gapStruct::ok;
    algorithms.emplace_back(temp);
    
    temp.name = "t1";
    temp.gapFn = huffman_A2022;
    temp.status = gapStruct::ok;
//    algorithms.emplace_back(temp);
    
    temp.name = "t2";
    temp.gapFn = huffman_B2022;
    temp.status = gapStruct::ok;
//    algorithms.emplace_back(temp);
    
    temp.name = "A";
    temp.gapFn = a;
    temp.status = gapStruct::ok;
    algorithms.emplace_back(temp);
    
    temp.name = "B";
    temp.gapFn = b;
    temp.status = gapStruct::ok;
    algorithms.emplace_back(temp);
}

static void errorFunction(vi &wc, vi &cc) {
    int n(0), w(28), maxLines(24);
    
    auto itw(wc.begin());
    auto itc(cc.begin());
    
    std::cout << std::right << std::setw(4) << "n"
    << std::right << std::setw(w) << "original"
    << std::right << std::setw(w) << "result" << '\n';
    while (itw != wc.end() && itc != cc.end() && n < maxLines)
        std::cout << std::right << std::setw(4) << ++n
        << std::right << std::setw(w) << *itc++
        << std::right << std::setw(w) << *itw++ << '\n';
}

static void doTimes(vgs algorithms, std::string fnBase) {
    std::fstream fst;
    fnBase += formatTime(true, true);
    fnBase += "-Results.csv";
    fst.open(fnBase, std::ios::out);
    fst << "Algorithm" << '\n';
    for (auto a : algorithms) {
        fst << a.name;
        for (auto pair : a.runData) {
            fst << ',' << pair.first;
            for (auto pp : pair.second) {
                fst << ',' << pp.time;
            }
            fst << '\n';
        }
    }
    fst << std::endl;
    fst.close();
}

static void doGaps(vgs algorithms, std::string fnBase) {
    std::vector<std::string> alsoRans;
    std::fstream gst;
    fnBase += formatTime(true, true);
    fnBase += "-Gaps.csv";
    gst.open(fnBase, std::ios::out);
    gst << "Algorithm" << '\n';
    for (auto a : algorithms) {
        gst << a.name;
        for (auto gap : a.gaps) {
            gst << ',' << gap;
        }
        gst << '\n';
    }
    gst << std::endl;
    gst.close();
}

static void makeFile(vgs algorithms) {
    std::string fnBase("/Users/prh/Keepers/code/xCode/shells/");
    doTimes(algorithms, fnBase);
    doGaps(algorithms, fnBase);
}

static void summerize(vgs &algorithms) {
    for (auto a : algorithms) {
        std::cout << '\n' << a.name;
        long n(a.gaps.size());
        for (auto itr(a.gaps.rbegin()); n--; itr++)
            std::cout << "  " << *itr;
        std::cout << '\n';
        for (auto pair : a.runData) {
            std::cout << pair.first << std::endl;
            for (auto t : pair.second) {
                std::cout << std::right << std::setw(12) << t.sampleSize
                << std::right << std::setw(14) << t.time
                << std::setw(18)<< convertMicroSeconds(t.time) << '\n';
            }
            std::cout << '\n';
        }
        std::cout << std::endl;
    }
}

static void eraseSlowerGaps(vgs &algorithms, ul averageFuncTime) {
    /*
     Removes elements that had subpar sort times
     */
    ul lim(averageFuncTime + (averageFuncTime >> 2) + (averageFuncTime >> 4)), laggardIndex(0);
    
    vul laggards;
    
    for (auto &a : algorithms) {
        for (auto &pair : a.runData) {
            if (pair.second.back().time > lim) {
                laggards.push_back(laggardIndex);
                std::cerr << "Removed " << a.name << '\n';
            }
        }
        laggardIndex++;
    }
    
    if (!laggards.empty()) {
        std::reverse(laggards.begin(), laggards.end());
        for (auto laggard : laggards) {
            algorithms.erase(algorithms.begin() + laggard);
        }
    }
}

static void getGaps(vgs &algorithms, ul sampleSize) {
    for (auto &a : algorithms) {
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

static void traverseAlgorithmVector(ul &activeFuncCount, vgs &algorithms, vi &checkCopy, const vi &orginalCopy, ul sampleSize, ul &totalFuncTime, int wdth, vi &workCopy, std::string distroName) {
    for (auto &a : algorithms) {
        workCopy.clear();
        workCopy = orginalCopy;
        auto start = high_resolution_clock::now();
        shellsort(workCopy, a.gaps);
        auto stop = high_resolution_clock::now();
        long duration = duration_cast<microseconds>(stop - start).count();
        std::cout << formatTime(false, true) << std::right << std::setw(31)
        << a.name << ": " <<std::setw(wdth) <<std::right << duration << " µs"
        << convertMicroSeconds(duration) << std::endl;
        sortMetrics sortMetrics;
        a.status = verify(workCopy, checkCopy) ? gapStruct::ok : gapStruct::outOfOrder;
        sortMetrics.time = duration;
        sortMetrics.sampleSize = sampleSize;
        if (a.gapStruct::status == a.gapStruct::outOfOrder)
            errorFunction(workCopy, checkCopy);
        a.runData[distroName].emplace_back(sortMetrics);
        totalFuncTime += duration;
        activeFuncCount++;
    }
}

static void runActiveAlgorithms(vgs &algorithms, vi &checkCopy, const vi &orginalCopy, ul sampleSize,  int wdth, vi &workCopy, std::string distroName) {
    ul totalFuncTime(0), activeFuncCount(0);
    traverseAlgorithmVector(activeFuncCount, algorithms, checkCopy, orginalCopy, sampleSize, totalFuncTime, wdth, workCopy, distroName);
//    eraseSlowerGaps(algorithms, totalFuncTime / activeFuncCount);
    totalFuncTime = 0;
    activeFuncCount = 0;
}

static void eoj(vgs &algorithms) {
    summerize(algorithms);
    makeFile(algorithms);
    algorithms.clear();
}

static void prep4size(vi &checkCopy, vi &orginalCopy, ul sampleSize, std::string distro) {
    randomFill(sampleSize, orginalCopy, distro);
    orginalCopy.shrink_to_fit();
    checkCopy = orginalCopy;
    std::sort(checkCopy.begin(), checkCopy.end());
    std::cout << '\n' << formatTime(true, true) << " n: " << sampleSize << " Distribution: " << distro << std::endl;
}

static void work(vgs &algorithms, std::vector<std::string> &distros) {
    int wdth(14);
    ul  ssMin(999999), ssMax(1999999999);
    std::cout << "\nStart: " << ssMin << "  Max: " << ssMax << '\n';
    
    vi orginalCopy, workCopy, checkCopy;
    for (ul sampleSize(ssMin); sampleSize < ssMax; sampleSize *= 23) {
        sampleSize |= 0xf;
        getGaps(algorithms, sampleSize);
        for (auto distro : distros) {
            prep4size(checkCopy, orginalCopy, sampleSize, distro);
            runActiveAlgorithms(algorithms, checkCopy, orginalCopy, sampleSize, wdth, workCopy, distro);
        }
    }
}

void setup() {
    std::vector<std::string> distros;
    vgs algorithms;
    fillDistros(distros);
    for (int i(0); 1 < 5; i++) {
        makeAlgorithmElements(algorithms);
        work(algorithms, distros);
        eoj(algorithms);
    }
}

void shell1959(vul &gaps, ul vSize) {
    ul gap(vSize);
    while (gap > 1) {
        gap >>= 1;
        gaps.push_back(gap);
    }
}

void frank1960(vul &gaps, ul vSize) {
    ul gap(vSize >> 1);
    while (gap) {
        gaps.push_back(gap | 1);
        gap >>= 1;
    }
}

void hibbard1963(vul &gaps, ul vSize) {
    ul gap(1);
    while (gap  < vSize) {
        gaps.push_back(gap);
        gap <<= 1;
        gap |= 1;
    }
}

void papernov1965(vul &gaps, ul vSize) {
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

bool is5smooth(ul n) {
    while (n % 5 == 0)
        n /= 5;
    return is3smooth(n);
}

void pratt1971(vul &gaps, ul vSize) {
    gaps.clear();
    for (ul n(1); n < vSize; n++)
        if (is3smooth(n))
            gaps.push_back(n);
    
}

void kunth1973(vul &gaps, ul vSize) {
    ul k(1), lim(vSize / 3);
    do {
        k *= 3;
        gaps.push_back((k - 1) >> 1);
    } while (gaps.back() < lim);
    gaps.pop_back();
}

bool mySeq(ul a, ul b) {return a > b;}

void sedgewick1982(vul &gaps, ul vSize) {
    ul k4(4), k2(1), gap(0), lim(vSize - (vSize / 8));
    std::set<ul, std::greater<>> gSet;
    gSet.insert(1);
    do {
        gap = k4 + 3 * k2 + 1;
        gSet.insert(gap);
        k4 *= 4;
        k2 *= 2;
    } while (gap < lim);
    gaps.clear();
    for (ul g : gSet)
        gaps.push_back(g);
}

void sedgewick1986(vul &gaps, ul vSize) {
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

void gonnet1991(vul &gaps, ul vSize) {
    ul k(vSize);
    while (k > 1) {
        k = (5 * k - 1) / 11;
        k = k > 1 ? k : 1;
        gaps.push_back(k);
    }
}

void tokuda1992(vul &gaps, ul vSize) {
    gaps.push_back(1);
    for (ul i(1); gaps.back() < vSize; i++) {
        double a(pow(2.25, static_cast<double>(i)));
        ul j((9.0 * a - 4.0) / 5.0);
        gaps.push_back(j | 1);
    }
}

void empirical2001(vul &gaps, ul vSize) {
    vul t {701, 301, 132, 57, 23, 10, 4, 1};
    gaps = t;
    ul nextGap((gaps.front() << 2) | 1);
    while (nextGap < vSize) {
        gaps.insert(gaps.begin(), nextGap);
        nextGap = (nextGap << 3) | 1;
        nextGap = nextGap % 5 ? nextGap : nextGap - 4;
    }
}
 
void b(vul &gaps, ul vSize) {
    ul tmp(0);
    while (vSize) {
        vSize >>= 2;
        tmp += vSize;
    }
    gaps.push_back(tmp);
    while (gaps.back() > 1) {
        gaps.push_back((gaps.back() >> 3) | 1);
    }
}

void a(vul &gaps, ul vSize) {
    gaps.push_back((vSize >> 1) + (vSize >> 3) | 1);
    while (gaps.back() > 1) {
        gaps.push_back((gaps.back() >> 3) | 1);
    }
}

void huffman_B2022(vul &gaps, ul vSize) {
    gaps.push_back(1);
    gaps.push_back(5);
    gaps.push_back(11);
    gaps.push_back(17);
    ul lim(vSize / 3);
    while (gaps.back() < lim) {
        ul k(2), gap(gaps.back() << 1), t(gaps.back());
        while (t > 0) {
            t = gaps.back() >> k++;
            gap += t;
        }
        gaps.push_back(gap);
    }
}

void huffman_A2022(vul &gaps, ul vSize) {
    gaps.push_back(1);
    gaps.push_back(5);
    gaps.push_back(11);
    gaps.push_back(17);
    ul lim(vSize / 3 + (vSize >> 1));
    while (gaps.back() < lim) {
        ul k(2), gap(gaps.back() << 1), t(gaps.back());
        while (t > 0) {
            t = gaps.back() >> k++;
            gap += t;
        }
        gaps.push_back(gap);
    }
}
