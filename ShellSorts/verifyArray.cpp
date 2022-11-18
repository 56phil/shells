//
//  verifyarray.cpp
//  verifyarray
//
//  Created by Phil Huffman on 12/3/21.
//

#include "verifyArray.hpp"

bool verify(vi wc, vi cc) {
    auto itw(wc.begin()), itc(cc.begin());
    while (itw != wc.end() && itc != cc.end())
        if (*itw++ != *itc++)
            return false;
    return true;
}
