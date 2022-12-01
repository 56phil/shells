//
//  core.cpp
//  core
//
//  Created by Phil Huffman on 12/21/21.
//

#include "core.hpp"

static void writeDistros(std::vector<std::string> &distros) {
    std::string fnBase("/Users/prh/Keepers/code/xCode/shells/results/");
    fnBase += formatTime(true, true);
    fnBase += "-distros.csv";
    std::fstream fst;
    fst.open(fnBase, std::ios::out);
    int maxLines(2222);
    msvd dMap;
    
    for (auto d : distros) {
        vd tmp;
        randomFill(maxLines, tmp, d);
        dMap[d] = tmp;
    }
    
    for(int indx(-1); indx < maxLines; indx++) {
        for (auto d : dMap) {
            if (indx < 0) {
                fst << ',' << d.first;
            } else {
                fst << ',' << d.second[indx];
            }
        }
        fst << '\n';
    }
    fst << std::endl;
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
    distros.push_back("Uniform - Sorted & Reversed");
    distros.push_back("Uniform - Sorted");
    std::sort(distros.begin(), distros.end());
    
    writeDistros(distros);
}

static void makeAlgorithmElements(vgs &gapStructs) {
    gapStructs.clear();
    
    gapStruct temp;
    
    temp.name = "Shell 1959";
    temp.gapFn = shell;
    temp.status = gapStruct::ok;
    gapStructs.emplace_back(temp);
    
    temp.name = "Frank & Lazarus 1960";
    temp.gapFn = frank;
    temp.status = gapStruct::ok;
    gapStructs.emplace_back(temp);
    
    temp.name = "Hibbard 1963";
    temp.gapFn = hibbard;
    temp.status = gapStruct::ok;
    gapStructs.emplace_back(temp);
    
    temp.name = "Papernov & Stasevich 1965";
    temp.gapFn = papernov;
    temp.status = gapStruct::ok;
    gapStructs.emplace_back(temp);
    
    temp.name = "Pratt 1971";
    temp.gapFn = pratt;
    temp.status = gapStruct::ok;
    gapStructs.emplace_back(temp);
    
    temp.name = "Knuth 1973";
    temp.gapFn = kunth;
    temp.status = gapStruct::ok;
    gapStructs.emplace_back(temp);
    
    temp.name = "Sedgewick 1982";
    temp.gapFn = sedgewick82;
    temp.status = gapStruct::ok;
    gapStructs.emplace_back(temp);
    
//    temp.name = "Sedgewick 1986";
//    temp.gapFn = sedgewick86;
//    temp.status = gapStruct::ok;
//    gapStructs.emplace_back(temp);
    
//    temp.name = "Gonnet & Baeza-Yates 1991";
//    temp.gapFn = gonnet;
//    temp.status = gapStruct::ok;
//    gapStructs.emplace_back(temp);
    
    temp.name = "Tokuda 1992";
    temp.gapFn = tokuda;
    temp.status = gapStruct::ok;
    gapStructs.emplace_back(temp);
    
//    temp.name = "Ciura 2001";
//    temp.gapFn = ciura;
//    temp.status = gapStruct::ok;
//    gapStructs.emplace_back(temp);
    
    temp.name = "a 2022";
    temp.gapFn = a;
    temp.status = gapStruct::ok;
    gapStructs.emplace_back(temp);
    
    temp.name = "b 2022";
    temp.gapFn = b;
    temp.status = gapStruct::ok;
    gapStructs.emplace_back(temp);
    
    temp.name = "c 2022";
    temp.gapFn = c;
    temp.status = gapStruct::ok;
    gapStructs.emplace_back(temp);
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

static void doTimes(vgs gapStructs, std::string fileName) {
    std::fstream fst;
    fileName += formatTime(true, true);
    fileName += "-Results.csv";
    fst.open(fileName, std::ios::out);
    fst << "Gap Sequence,Distribution";
    auto p(gapStructs.front().runData);
    auto q(p.begin());
    for (auto r : q->second) {
        fst << ',' << r.sampleSize;
    }
    fst << '\n';
    
    for (auto gapStruct : gapStructs) {
        for (auto pair : gapStruct.runData) {
            fst << gapStruct.name;
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

static void doGaps(vgs gapStruct, std::string fnBase) {
    std::vector<std::string> alsoRans;
    std::fstream gst;
    fnBase += formatTime(true, true);
    fnBase += "-Gaps.csv";
    gst.open(fnBase, std::ios::out);
    gst << "Algorithm" << '\n';
    for (auto gaoStruct : gapStruct) {
        gst << gaoStruct.name;
        for (auto gap : gaoStruct.gaps) {
            gst << ',' << gap;
        }
        gst << '\n';
    }
    gst << std::endl;
    gst.close();
}

static void makeFile(vgs gapStructs) {
    std::string fnBase("/Users/prh/Keepers/code/xCode/shells/results/");
    doTimes(gapStructs, fnBase);
    doGaps(gapStructs, fnBase);
}

static void summerize(vgs &gapStructs) {
    for (auto gapStruct : gapStructs) {
        std::cout << '\n' << gapStruct.name;
        long n(gapStruct.gaps.size());
        for (auto itr(gapStruct.gaps.begin()); n--; itr++)
            std::cout << "  " << *itr;
        std::cout << '\n';
        for (auto pair : gapStruct.runData) {
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

static void eraseSlowerGaps(vgs &gapStructs, ul averageFuncTime) {
    /*
     Removes elements that had subpar sort times
     */
    ul lim(averageFuncTime + (averageFuncTime >> 2) + (averageFuncTime >> 4)), laggardIndex(0);
    
    vul laggards;
    
    for (auto &gapstruct : gapStructs) {
        for (auto &pair : gapstruct.runData) {
            if (pair.second.back().time > lim) {
                laggards.push_back(laggardIndex);
                std::cerr << "Removed " << gapstruct.name << '\n';
            }
        }
        laggardIndex++;
    }
    
    if (!laggards.empty()) {
        std::reverse(laggards.begin(), laggards.end());
        for (auto laggard : laggards) {
            gapStructs.erase(gapStructs.begin() + laggard);
        }
    }
}

static void getGaps(vgs &gapStructs, ul sampleSize) {
    for (auto &gapStruct : gapStructs) {
        gapStruct.gaps.clear();
        gapStruct.gapFn(gapStruct.gaps, sampleSize);
        if (gapStruct.gaps.front() < gapStruct.gaps.back()) {
            std::reverse(gapStruct.gaps.begin(), gapStruct.gaps.end());
        }
        while (gapStruct.gaps.front() >= sampleSize) {
            gapStruct.gaps.erase(gapStruct.gaps.begin());
        }
        gapStruct.gaps.shrink_to_fit();
    }
}

static void traverseAlgorithmVector(ul &activeFuncCount, vgs &gapstricts, vd &checkCopy, const vd &orginalCopy, ul sampleSize, ul &totalFuncTime, int wdth, vd &workCopy, std::string distroName) {
    for (auto &gapStruct : gapstricts) {
        workCopy.clear();
        workCopy = orginalCopy;
        auto start = high_resolution_clock::now();
        shellsort(workCopy, gapStruct.gaps);
        auto stop = high_resolution_clock::now();
        long duration = duration_cast<microseconds>(stop - start).count();
        std::cout << formatTime(false, true) << std::right << std::setw(31)
        << gapStruct.name << ": " <<std::setw(wdth) <<std::right << duration << " µs"
        << convertMicroSeconds(duration) << std::endl;
        sortMetrics sortMetrics;
        gapStruct.status = verify(workCopy, checkCopy) ? gapStruct::ok : gapStruct::outOfOrder;
        sortMetrics.time = duration;
        sortMetrics.sampleSize = sampleSize;
        if (gapStruct.gapStruct::status == gapStruct.gapStruct::outOfOrder)
            errorFunction(workCopy, checkCopy);
        gapStruct.runData[distroName].emplace_back(sortMetrics);
        totalFuncTime += duration;
        activeFuncCount++;
    }
}

static void runActiveAlgorithms(vgs &gapStructs, vd &checkCopy, const vd &orginalCopy, ul sampleSize,  int wdth, vd &workCopy, std::string distroName) {
    ul totalFuncTime(0), activeFuncCount(0);
    traverseAlgorithmVector(activeFuncCount, gapStructs, checkCopy, orginalCopy, sampleSize, totalFuncTime, wdth, workCopy, distroName);
//    eraseSlowerGaps(gapStructs, totalFuncTime / activeFuncCount);
    totalFuncTime = 0;
    activeFuncCount = 0;
}

static void eoj(vgs &gapStructs) {
    summerize(gapStructs);
    makeFile(gapStructs);
    gapStructs.clear();
}

static void prep4size(vd &checkCopy, vd &orginalCopy, ul sampleSize, std::string distroName) {
    randomFill(sampleSize, orginalCopy, distroName);
    orginalCopy.shrink_to_fit();
    checkCopy = orginalCopy;
    std::sort(checkCopy.begin(), checkCopy.end());
    std::cout << '\n' << formatTime(true, true) << " n: " << sampleSize << " Distribution: " << distroName << std::endl;
}

static void work(vgs &gapStructs, vs distroNames) {
    int wdth(14);
    
    vi sampleSizes({999999, 1000000, 1000001});
    vd orginalCopy, workCopy, checkCopy;
    
    for (auto sampleSize : sampleSizes) {
        sampleSize = sampleSize > MAX_SAMPLE_SIZE ? MAX_SAMPLE_SIZE : sampleSize;
        getGaps(gapStructs, sampleSize);
        for (auto distroName : distroNames) {
            prep4size(checkCopy, orginalCopy, sampleSize, distroName);
            runActiveAlgorithms(gapStructs, checkCopy, orginalCopy, sampleSize, wdth, workCopy, distroName);
        }
    }
}

void setup() {
    vs distroNames;
    vgs gapStructs;
    fillDistros(distroNames);
    for (int outerLoopCounter(0); outerLoopCounter < MAX_OUTER_LOOP; outerLoopCounter++) {
        std::cerr << formatTime(true, true) << " Pass " << (outerLoopCounter + 1) << " of " << MAX_OUTER_LOOP << ".\n";
        makeAlgorithmElements(gapStructs);
        work(gapStructs, distroNames);
        eoj(gapStructs);
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
    int ladd(7), sra(3), strt(3);
    gaps.push_back((vSize - (vSize >> strt)) | 1);
    while (gaps.back() > ladd) {
        gaps.push_back((gaps.back() >> sra) | ladd);
    }
    gaps.push_back(1);
}

void b(vul &gaps, ul vSize) {
    int ladd(5), sra(3), strt(4);
    gaps.push_back((vSize - (vSize >> strt)) | 1);
    while (gaps.back() > ladd) {
        gaps.push_back((gaps.back() >> sra) | ladd);
    }
    gaps.push_back(1);
}

void c(vul &gaps, ul vSize) {
    int ladd(3), sra(3), strt(5);
    gaps.push_back((vSize - (vSize >> strt)) | 1);
    while (gaps.back() > ladd) {
        gaps.push_back((gaps.back() >> sra) | ladd);
    }
    gaps.push_back(1);
}
