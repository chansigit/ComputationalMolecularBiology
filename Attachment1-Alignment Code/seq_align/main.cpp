//    Computational Molecular Biology Assignment III: Sequence Alignment |
//    Author: Sijie Chen (2016310721), Department of Automation, THU     |
//    Oct. 31, 2016   11:02 am                                           |
#include <iostream>
#include <cstdio>
#include <string>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cassert>
#include <algorithm>
#define rep(a,b,c) for(int a=b;a<=c;++a)
#define INF 0x3f3f3f3f
//#define DEBUG_MODEL
using namespace std;
//************************************************************************
int prtable[256]={0};
void prtableInit(){
    prtable['A']=0;  prtable['R']=1;  prtable['N']=2;  prtable['D']=3;
    prtable['C']=4;  prtable['Q']=5;  prtable['E']=6;  prtable['G']=7;
    prtable['H']=8;  prtable['I']=9;  prtable['L']=10; prtable['K']=11;
    prtable['M']=12; prtable['F']=13; prtable['P']=14; prtable['S']=15;
    prtable['T']=16; prtable['W']=17; prtable['Y']=18; prtable['V']=19;
    prtable['B']=20; prtable['Z']=21; prtable['X']=22; prtable['*']=23;
}
int elemEncode(char elem,string type="protein"){
    int code=-1;
    if (type=="protein"){
        code = prtable[toupper(elem)];
    }
    return code;
}
void unittest_elemEncode(){
    string sarray="A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *";
    for (size_t i=0;i!=sarray.size();++i){
        if (sarray[i]!=' ')
            cout<<sarray[i]<<" "<<elemEncode(sarray[i])<<endl;
    }
}
//************************************************************************
int scoreMatrix[24][24]={0};
void scoreMatrixInit(string type="blosum62"){
    if (type=="blosum62"){
        int blosum62[24][24]={
            { 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4},
            {-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4},
            {-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4},
            {-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4},
            { 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4},
            {-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4},
            {-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
            { 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4},
            {-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4},
            {-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4},
            {-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4},
            {-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4},
            {-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4},
            {-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4},
            {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4},
            { 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4},
            { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4},
            {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4},
            {-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4},
            { 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4},
            {-2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4},
            {-1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
            { 0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4},
            {-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1}
        };
        memcpy(scoreMatrix, blosum62, sizeof(blosum62));
    }
}
int score(char x, char y, string type="protein"){
    int FirstIdx  = elemEncode(x, type);
    int SecondIdx = elemEncode(y, type);
    return scoreMatrix[FirstIdx][SecondIdx];
}
void unittest_score(){
    string sarray="ARNDCQEGHILKMFPSTWYVBZX*";
    for (size_t i=0;i!=sarray.size();++i){
        for (size_t j=0;j<=i;++j){
            if (sarray[i]!=' ' && sarray[j]!=' ')
                cout<<score(sarray[i],sarray[j])<<" ";
        }
        cout<<endl;
    }
}
//************************************************************************
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

void unittest_heapMatrix(){
    int** a=AllocateDpMatrix(3555,5551, 888);
    assert(a!=NULL);
    cout<<a<<endl;
    rep(i,0,2){
        rep(j,0,4) cout<<a[i][j]<<" ";
        cout<<endl;
    }
    ReleaseDpMatrix(a, 3); //# seems a memory leak issue
    cout<<a<<endl;
    assert(a==0);
}
//************************************************************************
void globalAlignment(const string& A, const string& B,
                    int gapPenalty, int**& dp,
                    string& alignA, string& alignB, int& finalScore){
    // Calculate dp matrix
    int lenA=A.size(), lenB=B.size();
    rep(i,0,lenA) dp[i][0]=gapPenalty*i;
    rep(j,0,lenB) dp[0][j]=gapPenalty*j;
    rep(i,1,lenA) rep(j,1,lenB){
        int matchScore =dp[i-1][j-1]+score(A[i-1], B[j-1], "protein");
        int deleteScore=dp[i-1][j]+gapPenalty; //Ai matches to B's gap
        int insertScore=dp[i][j-1]+gapPenalty; //Bj matches to A's gap
        dp[i][j]=std::max(std::max(matchScore,deleteScore), insertScore);
    }
    finalScore = dp[A.size()][B.size()];
    // Reconstructing alignment sequence
    int i=A.size(), j=B.size();
    alignA=""; alignB="";
    while (i>0 || j>0){
        if (i>0 && j>0 && dp[i][j]==dp[i-1][j-1]+score(A[i-1],B[j-1],"protein") ){
            alignA = A[--i] + alignA;
            alignB = B[--j] + alignB;
        }else if (i>0 && dp[i][j]==dp[i-1][j]+gapPenalty){
            alignA = A[--i] + alignA; //Ai matches to B's gap
            alignB = "-" + alignB;
        }else{
            alignA = "-" + alignA; //Bj matches to A's gap
            alignB = B[--j] + alignB;
        }
    }
}


void localAlignment(const string& A, const string& B,
                     int gapPanelty, int**& dp,
                     string& alignA, string& alignB, int& finalScore){
    // Calculate dp matrix
    int lenA=A.size(), lenB=B.size();
    rep(i,0,lenA) dp[i][0]=0;
    rep(j,0,lenB) dp[0][j]=0;
    int optScore=-INF, optI=0,optJ=0;
    rep(i,1,lenA)rep(j,1,lenB){
        int matchScore=dp[i-1][j-1]+score(A[i-1], B[j-1],"protein");
        int yGapScore=dp[i-1][j]+gapPanelty; //Ai matches to B's gap
        int xGapScore=dp[i][j-1]+gapPanelty; //Bj matches to A's gap
        dp[i][j]=std::max(0,
                          std::max(xGapScore,
                                   std::max(matchScore,yGapScore) ));
        if (dp[i][j]>optScore){
            optI=i;
            optJ=j;
            optScore=dp[i][j];
        }
        #ifdef DEBUG_MODE
        if(dp[i][j]==matchScore)
            printf("dp[%d][%d]=dp[%d][%d]+%c:%c=%d\n",i,j,i-1,j-1,A[i-1], B[j-1],score(A[i-1], B[j-1],"protein"));
        if(dp[i][j]==yGapScore)
            printf("dp[%d][%d]=dp[%d][%d]-8\n",i,j,i-1,j);
        if(dp[i][j]==xGapScore)
            printf("dp[%d][%d]=dp[%d][%d]-8\n",i,j,i,j-1);
        #endif // DEBUG_MODE
    }
    finalScore = dp[optI][optJ];
    // Reconstructing alignment sequence
    // Traceback from the optimum alignment position
    int i=optI, j=optJ;
    while (dp[i][j]!=0){
        if (dp[i][j] == dp[i-1][j-1]+score(A[i-1], B[j-1],"protein")){
            alignA = A[--i] + alignA;
            alignB = B[--j] + alignB;
        }else if (dp[i][j] == dp[i-1][j]+gapPanelty){ //Ai matches to B's gap
            alignA = A[--i] + alignA;
            alignB = "-" + alignB;
        }else if (dp[i][j] == dp[i][j-1]+gapPanelty){ //Bj matches to A's gap
            alignA = "-" + alignA;
            alignB = B[--j] + alignB;
        }else{
            ;
        }
    }
}
//************************************************************************
int main(){
    prtableInit();
    scoreMatrixInit();
    string a;
    string b;
    string response;
    bool termination=false;
    cout<<"*----------------------------------------------------------------------"<<endl;
    cout<<"* Computational Molecular Biology Assignment III: Sequence Alignment"<<endl;
    cout<<"* Author: Sijie Chen (2016310721)"<<endl;
    cout<<"*----------------------------------------------------------------------"<<endl;
    while (1){
        // Data Acquisition
        cout   <<"Enter Sequence 1:\n>>>"<<flush;  cin>>a;
        cout   <<"Enter Sequence 2:\n>>>"<<flush;  cin>>b;
        // Global Alignment
        cout   <<"Global Alignment Result:"<<endl;
        int**  dp1= AllocateDpMatrix(a.size()+1, b.size()+1, 0); //New a DP matrix
        string gAlignA, gAlignB;
        int    finalScoreG= 0;
        globalAlignment(a,b,-8,dp1,gAlignA,gAlignB,finalScoreG); //Call Global Alignment
        cout   <<finalScoreG <<endl;                             //Final Global Alignment Score
        cout   <<gAlignA     <<endl;                             //Global Alignment Sequence on A
        cout   <<gAlignB     <<endl;                             //Global Alignment Sequence on B
        ReleaseDpMatrix(dp1, a.size()+1);                        //Delete the DP matrix

        // Local Alignment
        int**  dp2=AllocateDpMatrix(a.size()+1, b.size()+1, 0);  //New a DP matrix
        cout   << "Local Alignment Result:" << endl;
        string lAlignA, lAlignB;
        int finalScoreL=0;
        localAlignment(a,b,-8,dp2,lAlignA,lAlignB,finalScoreL);  //Call Local Alignment
        cout   <<finalScoreL <<endl;                             //Final Local Alignment Score
        cout   <<lAlignA     <<endl;                             //Local Alignment Sequence on A
        cout   <<lAlignB     <<endl;                             //Local Alignment Sequence on B
        ReleaseDpMatrix(dp2, a.size()+1);                        //Delete the DP matrix

        // Hints for continue
        while (1){
            cout << "Continue?(y/n)\n>>>" <<flush;
            cin  >> response;
            if (response=="n"){
                termination=true;
                break;
            }
            if (response=="y") break;
        }
        if (termination){
            cout <<"bye!" <<endl;
            break;
        }
    }
    return 0;
}
