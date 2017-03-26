// SuffixArray_x64.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
//    Advanced Algorithm Assignment I: Sequence Alignment              |
//    Author: Sijie Chen (2016310721), Department of Automation, THU   |
//    Mar. 20, 2017   20:00 pm                                           |

// ********************************************************************
// You may find this program partially coincides with the program in
//     https://github.com/chansigit/ComputationalMolecularBiology/ \
//     blob/master/Attachment1-Alignment%20Code/seq_align/main.cpp
// That is my own Github repository, please do not regard it as plagiarism
//
// The FASTA reading module uses the following code for reference
//         http://rosettacode.org/wiki/FASTA_format#C.2B.2B
// That code obeys GNU Free Documentation License 1.2
//
// Copyright © 2017 Sijie Chen <chensj16-AT-mails-tsinghua-edu-cn>
// This work is free. You can redistribute it and/or modify it under the
// terms of the WTFPL, Version 2, as published by Sam Hocevar.
// See http://www.wtfpl.net/ for more details.
// ********************************************************************
#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cassert>
#include <algorithm>
#include <vector>
#include <iomanip>

namespace SuffixArray {
    const int
        MAXN = 158000000;

    typedef unsigned long long hash;

    const hash BASE = 137;

    int N;
    //char * S;
    std::string S;
    //int sa[MAXN];
    //hash h[MAXN], hPow[MAXN];
    int *sa;
    hash *h, *hPow;

#define getHash(lo, size) (h[lo] - h[(lo) + (size)] * hPow[size])

    inline bool sufCmp(int i, int j) {
        int lo = 1, hi = std::min(N - i, N - j);
        while (lo <= hi) {
            int mid = (lo + hi) >> 1;
            if (getHash(i, mid) == getHash(j, mid))
                lo = mid + 1;
            else
                hi = mid - 1;
        }
        return S[i + hi] < S[j + hi];
    }

    void buildSA() {
        N = (S.size());
        hPow[0] = 1;
        for (int i = 1; i <= N; ++i)
            hPow[i] = hPow[i - 1] * BASE;
        h[N] = 0;
        for (int i = N - 1; i >= 0; --i){
            h[i] = h[i + 1] * BASE + S[i], sa[i] = i;
        }
        
        std::stable_sort(sa, sa + N, sufCmp);
        
    }

} // end namespace HashSuffixArray


using namespace std;
int main() {
    std::clog << "Filename:" << std::endl;
    std::string seqFileName;
    std::cin >> seqFileName;
    std::ifstream infile(seqFileName, std::ios_base::in);
    std::vector<std::string> seqset, idset;
    std::string line, name, content;
    while (std::getline(infile, line).good()) {
        if (line.empty() || line[0] == '>') { // Identifier marker
            if (!name.empty()) { // Print out what we read from the last entry
                transform(content.begin(), content.end(), content.begin(), ::toupper);
                idset.push_back(name);
                seqset.push_back(content);
                name.clear();
            }
            if (!line.empty()) {
                name = line.substr(1);
            }
            content.clear();
        } else if (!name.empty()) {
            if (line.find(' ') != std::string::npos) { // Invalid sequence--no spaces allowed
                name.clear();
                content.clear();
            } else {
                content += line;
            }
        }
    }
    if (!name.empty()) { // Print out what we read from the last entry
        transform(content.begin(), content.end(), content.begin(), ::toupper);
        idset.push_back(name);
        seqset.push_back(content);
    }
    string s = "";
    for (auto reads : seqset) {
        s = s + reads + "$";
    }

    
    SuffixArray::sa = new int[SuffixArray::MAXN];
    SuffixArray::h = (SuffixArray::hash*) malloc(SuffixArray::MAXN * sizeof(SuffixArray::hash));
    SuffixArray::hPow = (SuffixArray::hash*) malloc(SuffixArray::MAXN * sizeof(SuffixArray::hash));

    SuffixArray::S = s;
    SuffixArray::buildSA();

    freopen("SuffixArray.txt", "w", stdout);

    std::cout << "The suffix array is listed as follow:" << std::endl;
    for (int cnt = 0, i = 0; i<s.size(); ++i) {
        std::cout << std::setw(8) << SuffixArray::sa[i] << " ";
        if (cnt == 8) {
            std::cout << std::endl;
            cnt = 0;
        } else {
            ++cnt;
        }
    }
}
