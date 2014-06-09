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
#include <string>
#include <gzstream.h>

#include "utils.h"

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


//frequency of A,C,G,T
typedef struct { 
    double f[4];
 } alleleFrequency;


//To store consensus information
typedef struct { 
    double perror[4];
    double phred[4];
    char consensus;
 } PHREDgeno;

//To store a single read
typedef struct { 
    char base;
    int  qual;
    int  mapq;       
    int dist5p;
    int dist3p;
    bool isReversed;
} singleRead;

void readNucSubstitionFreq(const string filename,vector<probSubstition> & subVec);
void readIlluminaError(const string errFile,probSubstition & illuminaErrorsProb);

#endif
