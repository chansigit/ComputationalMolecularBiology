//    Advanced Algorithm Assignment I: Sequence Alignment              |
//    Author: Sijie Chen (2016310721), Department of Automation, THU   |
//    Apr 10, 2017   16:03 pm                                           |

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
#include <list>
#define rep(a,b,c) for(int a=b;a<=c;++a)
#define LL long long
using namespace std;


inline int encode(char c){
    if (c=='a'||c=='A') return 0;
    if (c=='g'||c=='G') return 1;
    if (c=='c'||c=='C') return 2;
    if (c=='t'||c=='T')
        return 3;
    else
        return 4;

}

int hashFunction(string a, int base, int BIGPRIME, int hashSize){
    unsigned int k=a.size();
    unsigned int p=encode(a[k-1]);
    for (int i=k-2; i>=0; --i){
        p=encode(a[i])+base*p;
    }
    return size_t(1+(LL)p % (LL)BIGPRIME % (LL)hashSize);
}

#define HashTableSize  316040895
int BigPrime = 15727379237;
//vector<list<int>> HashTable(HashTableSize);
int hashArray[HashTableSize];
int kmerLen=100;

int main(int argc, char* argv[]){

    ifstream infile("chrX.fna",ios_base::in);
    vector<string> seqset, idset;
    std::string line, name, content;
    while( std::getline( infile, line ).good() ){
        if( line.empty() || line[0] == '>' ){ // Identifier marker
            if( !name.empty() ){ // Print out what we read from the last entry
                transform(content.begin(), content.end(), content.begin(), ::toupper);
                idset.push_back(name);
                seqset.push_back(content);
                name.clear();
            }
            if( !line.empty() ){
                name = line.substr(1);
            }
            content.clear();
        } else if( !name.empty() ){
            if( line.find(' ') != std::string::npos ){ // Invalid sequence--no spaces allowed
                name.clear();
                content.clear();
            } else {
                content += line;
            }
        }
    }
    if( !name.empty() ){ // Print out what we read from the last entry
        transform(content.begin(), content.end(), content.begin(), ::toupper);
        idset.push_back(name);
        seqset.push_back(content);
    }

    FILE* fp= fopen("hashTable.txt","w");
    int seqCnt=seqset.size();

    string genome = seqset[0];
    int genomeLen=seqset[0].size();
    cout<<genomeLen<<endl;
    for (int pos=0;pos<genomeLen-kmerLen;++pos){
        string kmer = genome.substr(pos,kmerLen);
        int hashVal=hashFunction(kmer,17, BigPrime, HashTableSize);
        if (pos%10000==0){
            printf("%d finished. hash= %d\n", pos,  hashVal);
        }
        int idx= hashVal;
        while (hashArray[idx]!=0) idx=(idx+1)%HashTableSize;
        if (idx==hashVal && hashArray[idx]!=0){
            cout<<"Hash Table Full!"<<endl;
            exit(1);
        }
        hashArray[idx]=pos;

        //HashTable[hashVal].push_back(pos);
        //fprintf(fp,"hashTable[%d]=%d", hashVal, pos);

    }
    return 0;
}
