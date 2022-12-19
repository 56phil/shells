//
//  formatMicroSeconds.cpp
//  formatMicroSeconds
//
//  Created by Phil Huffman on 12/21/21.
//

#include "formatMicroSeconds.hpp"

std::string formatMicroSeconds(const ul tms, int p) {
    const double kd(1000000);
    double seconds(tms / kd);
    ul minutes(seconds / 60);
    ul secs(seconds - minutes * 60);
    seconds -= minutes * 60;
    const ul hours(minutes / 60);
    minutes -= hours * 60;
    
    std::stringstream sst;
    sst.precision(p);
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
    
    sst << (hours == 0 && minutes == 0 && secs < 10 ? "      " : hours == 0 && minutes == 0 ? "     " : "")
    << std::fixed
    << std::setw((hours > 0 || minutes > 0) && secs < 10 ? 1 : 0)
    << ((hours > 0 || minutes > 0) && secs < 10 ? "0" : "")
    << seconds;
    
    return sst.str();
}
