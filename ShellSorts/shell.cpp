//
//  shell.cpp
//  shell
//
//  Created by Phil Huffman on 12/4/21.
//

#include "shell.hpp"

void shell(lv &v) {
    for (long gap(v.size() / 3 | 1); gap; gap >>= 1) {
        for (auto iti(v.begin() + gap); iti != v.end(); iti++) {
            int tmp = *iti;
            auto itj(iti);
            for (itj = iti; itj >= v.begin() + gap && *(itj - gap) > tmp; itj -= gap)
                *itj = *(itj - gap);
            *itj = tmp;
        }
    }
}
