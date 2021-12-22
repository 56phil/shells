//
//  main.cpp
//  ShellSorts
//
//  Created by Phil Huffman on 12/21/21.
//

#include "core.hpp"

int main(int argc, const char **argv) {
    int rtnCde(0);
    
    std::locale loc (std::cout.getloc(), new my_numpunct);
    std::cout.imbue(loc);

    setup();
    
    return rtnCde;
}
