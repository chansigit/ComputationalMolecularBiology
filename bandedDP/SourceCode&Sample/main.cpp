//    Advanced Algorithm Assignment I: Sequence Alignment              |
//    Author: Sijie Chen (2016310721), Department of Automation, THU   |
//    Mar. 12, 2017   15:00 pm                                           |

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
#define rep(a,b,c) for(int a=b;a<=c;++a)
#define INF 0x3f3f3f3f
#define MATCH 0
#define MISMATCH 1
#define INDEL 1
using namespace std;


int** AllocateDpMatrix(int row,int col, int defaultVal=0){
    int** elem = new int* [row];
    for (int i=0; i<row; i++)
        elem[i] = new int[col];
    if (elem!=NULL){
        rep(i,0,row-1)rep(j,0,col-1) elem[i][j]=defaultVal;
    }
    return elem;
}
int ReleaseDpMatrix(int**& elem, int row){
    for (int i=0;i<row;++i){
        delete[] elem[i];
        elem[i]=NULL;
    }
    delete[] elem;
    elem=NULL;
    return 0;
}


int score(string str1, string str2,int i, int j){
    char a= (i==0)? '_': str1[i-1];
    char b= (j==0)? '_': str2[j-1];
    //cout<<"a="<<a<<"  b="<<b<<endl;
    return a==b ? MATCH:MISMATCH;
}
void boundAlignment(const string& A, const string& B,
					int**& dp,
					string& alignA, string& alignB, int& finalScore, int k){
    // Calculate dp matrix
    int lenA=A.size(), lenB=B.size();
    // Initialization
    rep(i,0,lenA)rep(j,0,lenB)
        dp[i][j]=INF;


    rep(i,0,lenA) dp[i][0]=i*INDEL;
    rep(j,0,lenB) dp[0][j]=j*INDEL;

    for (int i=1;i<=lenA;++i){
        int lowerBound=max(1,    i-(k-abs(lenA-lenB))/2  );
        int upperBound=min(lenB, lenA+ i-1+(k-abs(lenA-lenB))/2  );
        //cout<<lowerBound<<"  "<<upperBound<<endl;

        for (int j=lowerBound;j<=upperBound;++j){
            int matchScore=dp[i-1][j-1] + score(A,B,i,j);
            int yGapScore=dp[i-1][j]+INDEL; //Ai matches to B's gap
            int xGapScore=dp[i][j-1]+INDEL; //Bj matches to A's gap
            dp[i][j]= std::min(matchScore,std::min(xGapScore,yGapScore) );
        }

    }

    finalScore = dp[lenA][lenB];
    if (abs(finalScore-INF) < INF*0.95){
        alignA=alignB="No Solution!!!";
        return;
    }
    // Reconstructing alignment sequence
    // Traceback from the optimum alignment position
    int i=lenA, j=lenB;
    while (i!=0 || j!=0){
        if (dp[i][j] == dp[i-1][j-1]+ score(A,B,i,j)){
            alignA = A[--i] + alignA;
            alignB = B[--j] + alignB;
        }else if (dp[i][j] == dp[i-1][j]+INDEL){ //Ai matches to B's gap
            alignA = A[--i] + alignA;
            alignB = "-" + alignB;
        }else if (dp[i][j] == dp[i][j-1]+INDEL){ //Bj matches to A's gap
            alignA = "-" + alignA;
            alignB = B[--j] + alignB;
        }else{
            ;
        }
    }
}

void printAlignment(string align1, string align2){
    string cacheA,cacheB,cacheS;
    int n=align1.size();
    int counter=0;
    for (int i=0;i<n;++i){
        counter++;
        cacheA=cacheA+align1[i];
        cacheB=cacheB+align2[i];
        if (align1[i]==align2[i])
            cacheS=cacheS+"|";
        else if (align1[i]=='-' || align2[i]=='-')
            cacheS=cacheS+" ";
        else
            cacheS=cacheS+".";
        if (counter%70==0){
            counter=0;
            //flush
            cout<<cacheA<<endl<<cacheS<<endl<<cacheB<<endl<<endl<<endl;
            cacheA=cacheB=cacheS="";
        }
    }
}

int main(int argc, char* argv[]){
    ifstream infile("HW1FILE1.txt",ios_base::in);
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

    freopen("AlignmentResult.txt","w",stdout);
    int seqCnt=seqset.size();
    for (int i=0;i<seqCnt;++i){
        for (int j=i+1;j<seqCnt;++j){
            cout<<"===================================================================="<<endl;
            string seqA=seqset[i];
            string seqB=seqset[j];
            string idA=idset[i];
            string idB=idset[j];
            if (seqA.size()>seqB.size()){
                swap(seqA,seqB);
                swap(idA, idB);
            }
            string gAlignA, gAlignB;
            int finalScore=0;
            int**  dp1= AllocateDpMatrix(seqA.size()+1, seqB.size()+1, 0); //New a DP matrix
            boundAlignment(seqA,seqB,dp1,gAlignA,gAlignB,finalScore,400);
            cout<<"["<<idA<<"]"<<endl<<"v.s."<<endl<<"["<<idB<<"]"<<endl;
            cout<<"Alignment Distance="<<finalScore<<"    Alignment Length="<<gAlignB.length()<<endl;
            printAlignment(gAlignA,gAlignB);
            ReleaseDpMatrix(dp1, seqA.size()+1);                        //Delete the DP matrix
        }
    }

    return 0;
}

