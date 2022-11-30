//
//  core.cpp
//  core
//
//  Created by Phil Huffman on 12/21/21.
//

#include "core.hpp"
static void displayDistros(std::vector<std::string> &distros) {
    int maxLines(22), wdth(30);
    msvd dMap;
    
    for (auto d : distros) {
        vd tmp;
        randomFill(maxLines, tmp, d);
        dMap[d] = tmp;
    }
    
    for(int indx(-1); indx < maxLines; indx++) {
        for (auto d : dMap) {
            if (indx < 0) {
                std::cout << std::right << std::setw(wdth) << d.first;
            } else {
                std::cout << std::right << std::setw(wdth) << d.second[indx];
            }
        }
        std::cout << '\n';
    }
    std::cout << std::endl;
}

static void fillDistros(std::vector<std::string> &distros) {
    distros.push_back("Bernoulli");
    distros.push_back("Binomial");
    distros.push_back("Chi Squared");
    distros.push_back("Cauchy");
    distros.push_back("Lognormal");
    distros.push_back("Normal");
    distros.push_back("Student T");
    distros.push_back("Uniform");
    distros.push_back("Uniform - Sorted - Reversed");
    distros.push_back("Uniform - Sorted");
    std::sort(distros.begin(), distros.end());
    
//    displayDistros(distros);
}

static void makeAlgorithmElements(vgs &algorithms) {
    algorithms.clear();
    
    gapStruct temp;
    temp.name = "Shell 1959";
    temp.gapFn = shell;
    temp.status = gapStruct::ok;
//    algorithms.emplace_back(temp);
    
    temp.name = "Frank & Lazarus 1960";
    temp.gapFn = frank;
    temp.status = gapStruct::ok;
//    algorithms.emplace_back(temp);
    
    temp.name = "Hibbard 1963";
    temp.gapFn = hibbard;
    temp.status = gapStruct::ok;
//    algorithms.emplace_back(temp);
    
    temp.name = "Papernov & Stasevich 1965";
    temp.gapFn = papernov;
    temp.status = gapStruct::ok;
//    algorithms.emplace_back(temp);
    
    temp.name = "Pratt 1971";
    temp.gapFn = pratt;
    temp.status = gapStruct::ok;
//    algorithms.emplace_back(temp);
    
    temp.name = "Knuth 1973";
    temp.gapFn = kunth;
    temp.status = gapStruct::ok;
    algorithms.emplace_back(temp);
    
    temp.name = "Sedgewick 1982";
    temp.gapFn = sedgewick82;
    temp.status = gapStruct::ok;
    algorithms.emplace_back(temp);
    
    temp.name = "Sedgewick 1986";
    temp.gapFn = sedgewick86;
    temp.status = gapStruct::ok;
    algorithms.emplace_back(temp);
    
    temp.name = "Gonnet & Baeza-Yates 1991";
    temp.gapFn = gonnet;
    temp.status = gapStruct::ok;
    algorithms.emplace_back(temp);
    
    temp.name = "Tokuda 1992";
    temp.gapFn = tokuda;
    temp.status = gapStruct::ok;
    algorithms.emplace_back(temp);
    
    temp.name = "Ciura 2001";
    temp.gapFn = ciura;
    temp.status = gapStruct::ok;
    algorithms.emplace_back(temp);
    
    temp.name = "a 2022";
    temp.gapFn = a;
    temp.status = gapStruct::ok;
    algorithms.emplace_back(temp);
    
    temp.name = "b 2022";
    temp.gapFn = b;
    temp.status = gapStruct::ok;
    algorithms.emplace_back(temp);
    
    temp.name = "c 2022";
    temp.gapFn = c;
    temp.status = gapStruct::ok;
    algorithms.emplace_back(temp);
}

static void errorFunction(vd &wc, vd &cc) {
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

static void doTimes(vgs algorithms, std::string fileName) {
    std::fstream fst;
    fileName += formatTime(true, true);
    fileName += "-Results.csv";
    fst.open(fileName, std::ios::out);
    fst << "Algorithm,Distribution";
    auto p(algorithms.front().runData);
    auto q(p.begin());
    for (auto r : q->second) {
        fst << ',' << r.sampleSize;
    }
    fst << '\n';
    
    for (auto a : algorithms) {
        for (auto pair : a.runData) {
            fst << a.name;
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
    std::string fnBase("/Users/prh/Keepers/code/xCode/shells/results/");
    doTimes(algorithms, fnBase);
    doGaps(algorithms, fnBase);
}

static void summerize(vgs &algorithms) {
    for (auto a : algorithms) {
        std::cout << '\n' << a.name;
        long n(a.gaps.size());
        for (auto itr(a.gaps.begin()); n--; itr++)
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

static void traverseAlgorithmVector(ul &activeFuncCount, vgs &algorithms, vd &checkCopy, const vd &orginalCopy, ul sampleSize, ul &totalFuncTime, int wdth, vd &workCopy, std::string distroName) {
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

static void runActiveAlgorithms(vgs &algorithms, vd &checkCopy, const vd &orginalCopy, ul sampleSize,  int wdth, vd &workCopy, std::string distroName) {
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

static void prep4size(vd &checkCopy, vd &orginalCopy, ul sampleSize, std::string distro) {
    randomFill(sampleSize, orginalCopy, distro);
    orginalCopy.shrink_to_fit();
    checkCopy = orginalCopy;
    std::sort(checkCopy.begin(), checkCopy.end());
    std::cout << '\n' << formatTime(true, true) << " n: " << sampleSize << " Distribution: " << distro << std::endl;
}

static void work(vgs &algorithms, std::vector<std::string> &distros) {
    int wdth(14);
    
    vi sampleSizes({933333, 944444, 9555555});
    vd orginalCopy, workCopy, checkCopy;
    
    for (auto sampleSize : sampleSizes) {
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
    for (int outerLoopCounter(0); outerLoopCounter < MAX_OUTER_LOOP; outerLoopCounter++) {
        std::cerr << formatTime(true, true) << " Pass " << (outerLoopCounter + 1) << " of " << MAX_OUTER_LOOP << ".\n";
        makeAlgorithmElements(algorithms);
        work(algorithms, distros);
        eoj(algorithms);
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

bool is5smooth(ul n) {
    while (n % 5 == 0)
        n /= 5;
    return is3smooth(n);
}

void pratt(vul &gaps, ul vSize) {
    gaps.clear();
    for (ul n(1); n < vSize; n++)
        if (is3smooth(n))
            gaps.push_back(n);
    
}

void kunth(vul &gaps, ul vSize) {
    ul k(1), lim(vSize / 3);
    do {
        k *= 3;
        gaps.push_back((k - 1) >> 1);
    } while (gaps.back() < lim);
    gaps.pop_back();
}

bool mySeq(ul a, ul b) {return a > b;}

void sedgewick82(vul &gaps, ul vSize) {
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
    gaps.push_back((vSize >> 1) + (vSize >> 2) + (vSize >> 3)| 1);
    while (gaps.back() > 1) {
        gaps.push_back((gaps.back() >> 3) | 1);
    }
}

void b(vul &gaps, ul vSize) {
    gaps.push_back((vSize >> 1) + (vSize >> 2) + (vSize >> 3)| 1);
    while (gaps.back() > 15) {
        gaps.push_back((gaps.back() >> 3) | 15);
    }
    gaps.push_back(1);
}

void c(vul &gaps, ul vSize) {
    gaps.push_back((vSize >> 1) + (vSize >> 2) + (vSize >> 3)| 1);
    while (gaps.back() > 7) {
        gaps.push_back((gaps.back() >> 3) | 7);
    }
    gaps.push_back(1);
}
