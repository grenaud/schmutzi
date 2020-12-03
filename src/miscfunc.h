/*
 * miscfunc
 * Date: Jun-08-2014 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef miscfunc_h
#define miscfunc_h

#include <stdlib.h>
#include <vector>
#include <set>
#include <string>
#include <gzstream.h>

#include "libgab.h"

using namespace std;

typedef struct { 
    double s[12];
 } substitutionRates;


//  model->obs
//  0  A->A 
//  1  A->C 
//  2  A->G 
//  3  A->T 
//  4  C->A 
//  5  C->C 
//  6  C->G 
//  7  C->T 
//  8  G->A 
//  9  G->C 
//  10 G->G 
//  11 G->T 
//  12 T->A 
//  13 T->C 
//  14 T->G 
//  15 T->T 

typedef struct { 
    double s[16];
} probSubstition;

typedef struct { 
    double p[4][4];
} diNucleotideProb;


//frequency of A,C,G,T
typedef struct { 
    long double f[4];
 } alleleFrequency;


//To store consensus information
typedef struct { 
    long double perror[4];
    long double phred[4];
    long double perrorC[4];
    long double phredC[4];

    char ref;
    char consensus; //for the endogenous
} PHREDgeno;


typedef struct { 

    char ref;
    char base; //for the endogenous
    int  pos;
    double q;

    double aprob;
    double cprob;
    double gprob;
    double tprob;

} logRecord;

//To store a single read
typedef struct { 
    char base;
    int  qual;
    int  mapq;       
    int dist5p;
    int dist3p;
    bool isReversed;
} singleRead;

//To store contamination likelihood
typedef struct { 
    string filename;    
    long double contaminationRate;
    long double logLike;
    //    long double logLike;
} contaminationEstimate;


typedef struct { 
    vector<singleRead> readsVec;
    long double mapqAvg;
    int cov;
    char refBase;
    int posAlign;
    bool skipPosition;
} positionInformation;


typedef struct { 
    long double likeBaseNoindel[4];
    long double likeBaseNoindelCont[4][4];
    long double likeBaseNoindelNoBoundary[4];           //when we do not consider bases at the ends of reads
    long double likeBaseNoindelContNoBoundary[4][4];    //when we do not consider bases at the ends of reads

    int  covPerBase[4];
    int  covPerBaseNoBoundary[4];

    long double mapqAvg;
    
    int numDel;
    long double llikDeletion;
    long double llikNoDeletion;

    long double llikDeletionBoth;
    long double llikDeletionCont;
    long double llikDeletionEndo;
    long double llikDeletionNone;

    set<string> allInserts;

    vector<string> insertionRight;
    map< pair<string,string> , long double> insertion2loglikeEndoCont; //key is (endo ins,cont ins) to log likelihood

    map<string,int> insertion2count;
    map<string,long double> insertion2loglike;
    //map<string,long double> insertion2loglikeCont;

    int cov;
} singlePosInfo;


void readNucSubstitionFreq(const string filename,vector<probSubstition> & subVec);
void readIlluminaError(const string errFile,probSubstition & illuminaErrorsProb);
void readMTConsensus(const string consensusFile,map<int, PHREDgeno> & pos2phredgeno,int & sizeGenome,vector<int> & posOfIndels);
void readMTAlleleFreq(const string freqFile,	map<int, alleleFrequency> & pos2allelefreq);



// Returns log10( pow(10,x)+pow(10,y) ), but does so without causing
// overflow or loss of precision.
/* template <typename T> */
/* inline T oplusInit(T x,T y ){ */
/*     if( x == 0 ){ //no initialized, a log = 0 should not exist */
/* 	return y; */
/*     } */

/*     return x > y  */
/*         ? x + log1p( pow( 10, y-x ) ) / log(10) */
/*         : y + log1p( pow( 10, x-y ) ) / log(10) ; */
/* } */


#endif
