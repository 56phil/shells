//
//  shell.cpp
//  shell
//
//  Created by Phil Huffman on 12/4/21.
//

#include "shell.hpp"

void shellsort(vi &v, const vull &gaps) {
    int lim(static_cast<int>(v.size() - 1));
    for (auto gap : gaps) {
        if (gap < lim) {
            for (auto iti(v.begin() + gap); iti != v.end(); iti++) {
                int tmp = *iti;
                auto itj(iti);
                for (itj = iti; itj >= v.begin() + gap && *(itj - gap) > tmp; itj -= gap)
                    *itj = *(itj - gap);
                *itj = tmp;
            }
        }
    }
}
