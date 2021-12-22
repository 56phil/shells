//
//  shell.cpp
//  shell
//
//  Created by Phil Huffman on 12/4/21.
//

#include "shell.hpp"

void shellsort(vi &v, const vi gaps) {
    long lim(v.size() / 3);
    for (auto gap : gaps) {
        if (gap > lim)
            continue;
        for (auto iti(v.begin() + gap); iti != v.end(); iti++) {
            int tmp = *iti;
            auto itj(iti);
            for (itj = iti; itj >= v.begin() + gap && *(itj - gap) > tmp; itj -= gap)
                *itj = *(itj - gap);
            *itj = tmp;
        }
    }
}
