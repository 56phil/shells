//
//  formatMicroSeconds.cpp
//  formatMicroSeconds
//
//  Created by Phil Huffman on 12/21/21.
//

#include "formatMicroSeconds.hpp"

std::string formatMicroSeconds(const ul tms, int p) {
    const ul ku(1000000);
    const ul fractional(tms % ku);
    ul seconds(tms / ku);
    ul minutes(seconds / 60);
    seconds -= minutes * 60;
    const ul hours(minutes / 60);
    minutes -= hours * 60;
    
    std::stringstream sst;
    sst << fractional;
    std::string frac(sst.str());
    while (frac.size() > p) {
        frac.pop_back();
    }

    sst << std::setfill(' ');
    if (hours > 0)
        sst << (hours > 99 ? "  " : hours > 9 ? "   " : "    " ) << hours << ":";
    else
        sst << "    ";
    
    sst << std::setfill('0');
    if (minutes > 0 || hours > 0)
        sst << (hours > 0 ? "" : minutes > 9 ? "  " : "   ") << std::fixed
        << std::setw(minutes > 9 || hours > 0 ? 2 : 1)
        << minutes << ":";
    
    sst << (hours == 0 && minutes == 0 && seconds < 10 ? "      " : hours == 0 && minutes == 0 ? "     " : "")
    << std::fixed
    << std::setw(hours == 0 && minutes == 0 && seconds < 10 ? 1 : 2) << seconds;
    if (p)
        sst << "." << frac;
    
    return sst.str();
}
