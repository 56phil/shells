//
//  formattime.cpp
//  formattime
//
//  Created by Phil Huffman on 12/8/21.
//

#include "formattime.hpp"

std::string formatTime (bool doDate=false, bool doTime=true) {
    time_t rawtime;
    struct tm *timeinfo;
    char buffer[80];
    
    time (&rawtime);
    timeinfo = localtime (&rawtime);
    
    std::stringstream sst;
    
    if (doDate) {
        strftime (buffer,80,"%Y-%m-%d", timeinfo);
        sst << buffer;
        if (doTime)
            sst << "T";
    }
    if (doTime) {
        strftime (buffer,80,"%H:%M:%S", timeinfo);
        sst << buffer;
    }
    std::string rv(sst.str());

    return rv;
}
