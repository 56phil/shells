//
//  core.cpp
//
//  Created by Phil Huffman on 12/21/21.
//

#include "core.hpp"

static void writeDistros(vs &distros) {
//    int maxLines(2222);
    msvi dMap;
    
    for (auto d : distros) {
        vi tmp (MAX_DistroLines);
        randomFill(MAX_DistroLines, tmp, d);
        if (tmp.empty()) {
            std::cerr << "random fill failed.\n";
            exit(1);
        }
        dMap[d] = tmp;
    }
    
    std::string fnBase("/Users/prh/Keepers/code/xCode/shells/results/");
    fnBase += formatTime(true, true);
    fnBase += "-distros.csv";
    std::fstream fst;
    fst.open(fnBase, std::ios::out);
    for(long indx(-1); indx < MAX_DistroLines; indx++) {
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
    fst.close();
}

static void fillDistros(vs &distros) {
    distros.push_back("Bernoulli");
    distros.push_back("Binomial");
    distros.push_back("Normal");
    distros.push_back("Poisson");
    distros.push_back("Uniform");
    distros.push_back("Uniform - Sorted & Reversed");
    distros.push_back("Uniform - Sorted");
    
    std::sort(distros.begin(), distros.end());
    
    writeDistros(distros);
}

static void makeGapSequenceGenerators(vgs &gapStructs) {
    gapStructs.clear();
    
    gapStruct temp;
    temp.warnings = 0;
    temp.status = gapStruct::ok;

//    temp.name = "Shell 1959";
//    temp.gapFn = shell;
//    gapStructs.emplace_back(temp);
//
//    temp.name = "Frank & Lazarus 1960";
//    temp.gapFn = frank;
//    gapStructs.emplace_back(temp);
//
//    temp.name = "Hibbard 1963";
//    temp.gapFn = hibbard;
//    gapStructs.emplace_back(temp);
//
//    temp.name = "Papernov & Stasevich 1965";
//    temp.gapFn = papernov;
//    gapStructs.emplace_back(temp);
//
//    temp.name = "Pratt 1971";
//    temp.gapFn = pratt;
//    gapStructs.emplace_back(temp);
//
//    temp.name = "Pratt 1971 (modified)";
//    temp.gapFn = pratt_A;
//    gapStructs.emplace_back(temp);
//
//    temp.name = "Knuth 1973";
//    temp.gapFn = kunth;
//    gapStructs.emplace_back(temp);

    temp.name = "Sedgewick 1982";
    temp.gapFn = sedgewick82;
    gapStructs.emplace_back(temp);

//    temp.name = "Sedgewick 1986";
//    temp.gapFn = sedgewick86;
//    gapStructs.emplace_back(temp);
//
//    temp.name = "Gonnet & Baeza-Yates 1991";
//    temp.gapFn = gonnet;
//    gapStructs.emplace_back(temp);
//
//    temp.name = "Tokuda 1992";
//    temp.gapFn = tokuda;
//    gapStructs.emplace_back(temp);
//
//    temp.name = "Ciura 2001";
//    temp.gapFn = ciura;
//    gapStructs.emplace_back(temp);
//
    temp.name = "a 2022";
    temp.gapFn = a;
    gapStructs.emplace_back(temp);
//
//    temp.name = "b 2022";
//    temp.gapFn = a;
//    gapStructs.emplace_back(temp);
//
//    temp.name = "c 2022";
//    temp.gapFn = a;
//    gapStructs.emplace_back(temp);
//
//    temp.name = "d 2022";
//    temp.gapFn = d;
//    gapStructs.emplace_back(temp);
//
//    temp.name = "e 2022";
//    temp.gapFn = e;
//    gapStructs.emplace_back(temp);
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
                << std::setw(18)<< formatMicroSeconds(t.time) << '\n';
            }
            std::cout << '\n';
        }
        std::cout << std::endl;
    }
}

static void eraseSlowerGaps(vgs &gapStructs, ull averageFuncTime, std::string distroName) {
    /*
     Removes elements that had subpar sort times
     */
    ul laggardIndex(0);
    vul laggards;
    ull lim(averageFuncTime + (averageFuncTime >> 3));
    
    if (cullSlowerGapSequences && gapStructs.size() > MIN_ActiveGapStructs) {
        for (auto &gapstruct : gapStructs) {
            if (gapstruct.runData[distroName].back().time > lim) {
                gapstruct.warnings++;
                if (gapstruct.warnings < MAX_Warnings) {
                    std::cerr << "Warned " << gapstruct.name << "  (" << gapstruct.warnings << " of " << MAX_Warnings << ")\n";
                } else {
                    laggards.push_back(laggardIndex);
                    std::cerr << "Removed " << gapstruct.name << '\n';
                    if (gapStructs.size() - laggards.size() <= MIN_ActiveGapStructs) {
                        break;
                    }
                }
            }
            laggardIndex++;
        }
    
        if (!laggards.empty()) {
            std::reverse(laggards.begin(), laggards.end());
            std::cerr << "Another " << laggards.size() << " bite" << (laggards.size() == 1 ? "s " : " ") << "the dust!\n";
            for (auto laggard : laggards) {
                gapStructs.erase(gapStructs.begin() + laggard);
            }
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

std::string gaps2string(vul gapVect) {
    std::stringstream sst;
    for (auto gap : gapVect) {
        sst << " " << gap;
    }
    std::string oString(sst.str());
    return oString;
}

ul median(vul vl) {
    if ((vl.size() & 1) == 0 && vl.back() != 0) {
        vl.push_back(vl.front() / vl.back()); // a cheap trick to make sure size is odd
    }
    std::sort(vl.begin(), vl.end());
    return  vl[vl.size() >> 1];
}

static void traverseVGS(ull &activeFuncCount, vgs &gapstructs, vi &chkCpy, const vi &orgCpy, ul sampleSize, ull &totalFuncTime, int wdth, vi &wrkCpy, std::string dName) {
    vul times(MEDIAN_TrialSize);
    for (auto &gapStruct : gapstructs) {
        times.clear();
        for (int ex(0); ex < MEDIAN_TrialSize; ex++) {
            wrkCpy.clear();
            wrkCpy = orgCpy;
            auto start = high_resolution_clock::now();
            shellsort(wrkCpy, gapStruct.gaps);
            auto stop = high_resolution_clock::now();
            ul dur = duration_cast<microseconds>(stop - start).count();
            times.push_back(dur);
        }
        auto duration(median(times));
        std::cout << formatTime(false, true) << std::right << std::setw(31)
        << gapStruct.name << ": " <<std::setw(wdth) <<std::right << duration << " µs"
        << formatMicroSeconds(duration)
//        << '\n' << gaps2string(gapStruct.gaps)
        << '\n';
        sortMetrics sortMetrics;
        gapStruct.status = verify(wrkCpy, chkCpy) ? gapStruct::ok : gapStruct::outOfOrder;
        sortMetrics.time = duration;
        sortMetrics.sampleSize = sampleSize;
        if (gapStruct.gapStruct::status == gapStruct.gapStruct::outOfOrder)
            errorFunction(wrkCpy, chkCpy);
        gapStruct.runData[dName].emplace_back(sortMetrics);
        totalFuncTime += duration;
        activeFuncCount++;
    }
}

static void runActiveAlgorithms(vgs &gapStructs, vi &checkCopy, const vi &orginalCopy, ul sampleSize,  int wdth, vi &workCopy, std::string distroName) {
    ull totalFuncTime(0), activeFuncCount(0);
    traverseVGS(activeFuncCount, gapStructs, checkCopy, orginalCopy, sampleSize, totalFuncTime, wdth, workCopy, distroName);
    eraseSlowerGaps(gapStructs, totalFuncTime / activeFuncCount, distroName);
    totalFuncTime = 0;
    activeFuncCount = 0;
}

static void eoj(vgs &gapStructs) {
    summerize(gapStructs);
    makeFile(gapStructs);
    gapStructs.clear();
}

static void prep4size(vi &checkCopy, vi &orginalCopy, ul sampleSize, std::string distroName) {
    randomFill(sampleSize, orginalCopy, distroName);
    orginalCopy.shrink_to_fit();
    checkCopy = orginalCopy;
    std::sort(checkCopy.begin(), checkCopy.end());
    std::cout << '\n' << formatTime(true, true) << " n: " << sampleSize << " Distribution: " << distroName << std::endl;
}

static void work(vgs &gapStructs, vs distroNames) {
    int wdth(14);
    
    vul sampleSizes({1999999, 0xffffffffff});
    std::sort(sampleSizes.begin(), sampleSizes.end());
    
    for (auto sampleSize : sampleSizes) {
        sampleSize = sampleSize > MAX_SampleSize ? MAX_SampleSize : sampleSize;
        vi originalCopy(sampleSize), workCopy(sampleSize), checkCopy(sampleSize);
        getGaps(gapStructs, sampleSize);
        for (auto distroName : distroNames) {
            prep4size(checkCopy, originalCopy, sampleSize, distroName);
            runActiveAlgorithms(gapStructs, checkCopy, originalCopy, sampleSize, wdth, workCopy, distroName);
        }
    }
}

void setup() {
    vs distroNames;
    fillDistros(distroNames);
    for (int outerLoopCounter(0); outerLoopCounter < MAX_Passes; outerLoopCounter++) {
        std::cerr << formatTime(true, true) << " Pass " << (outerLoopCounter + 1) << " of " << MAX_Passes << ".\n";
        vgs gapStructs;
        makeGapSequenceGenerators(gapStructs);
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
    ul k4(4), k2(1), gap(0), lim(vSize - (vSize >> 3));
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

static void base_22(vul &gaps, int ladd, int sra_0, int sra_1, int sra_2, int strt_0, int strt_1, int strt_2, ul vSize) {
    gaps.push_back(((vSize >> strt_0) + (vSize >> strt_1) - (vSize >> strt_2)) | ladd);
    while (gaps.back() > ladd) {
        gaps.push_back(((gaps.back() >> sra_0) + (gaps.back() >> sra_1) - (gaps.back() >> sra_2)) | ladd);
    }
    gaps.push_back(1);
}

void a(vul &gaps, ul vSize) {
    const int ladd(7), sra_0(2), sra_1(7), sra_2(10), strt_0(1), strt_1(2), strt_2(1);
    base_22(gaps, ladd, sra_0, sra_1, sra_2, strt_0, strt_1, strt_2, vSize);
}

void b(vul &gaps, ul vSize) {
    const int ladd(7), sra_0(2), sra_1(7), sra_2(10), strt_0(1), strt_1(2), strt_2(1);
    base_22(gaps, ladd, sra_0, sra_1, sra_2, strt_0, strt_1, strt_2, vSize);
}

void c(vul &gaps, ul vSize) {
    const int ladd(7), sra_0(2), sra_1(7), sra_2(10), strt_0(1), strt_1(2), strt_2(1);
    base_22(gaps, ladd, sra_0, sra_1, sra_2, strt_0, strt_1, strt_2, vSize);
}

void d(vul &gaps, ul vSize) {
    const int ladd(7), sra_0(2), sra_1(7), sra_2(10), strt_0(1), strt_1(2), strt_2(1);
    base_22(gaps, ladd, sra_0, sra_1, sra_2, strt_0, strt_1, strt_2, vSize);
}

void e(vul &gaps, ul vSize) {
    const int ladd(7), sra_0(2), sra_1(7), sra_2(10), strt_0(1), strt_1(2), strt_2(1);
    base_22(gaps, ladd, sra_0, sra_1, sra_2, strt_0, strt_1, strt_2, vSize);
}
