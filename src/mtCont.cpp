/*
 * mtCont
 * Date: Mar-26-2014 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

// #define DEBUG1
// #define DEBUG2
// #define DEBUGPOS 714 //position to debug
// #define DEBUGPOS 608 //position to debug

//#define DEBUGPOS 11298
//#define DEBUGPOS 2830

//#define DEBUGPOSLOGLIKE
//#define DEBUGPOSPRINTSINGLEPOS
//#define DEBUGSUBPATTERN //to debug sub. pattern
//#define DEBUGPOSEACHREAD //to debug each read
//#define DEBUGCONTPOS  //print pos that are potential contaminants
//#define DEBUGMTPOSSKIP
//#define DEBUGCONTRATE 0.01

#define MAXCOV 5000            // beyond that coverage, we stop computing
#define MINPRIOR 0.1         // below that prior for contamination at a given base, we give up
#define MAXMAPPINGQUAL 257     // maximal mapping quality, should be sufficient as mapping qualities are encoded using 8 bits
#define MINDISTFROMINDEL 5     // ignore positions that are within X bp of an indel


//TODO 
// mapping TODO
// probEndoVec seems to have been used twice


// CODE ORGANIZATION
//
// main()
//   call initScores() to initialize probability scores
//   read arguments, contaminant profiles are stored in queueFilesToprocess
//   read error and deamination profiles
//   read the endogenous consensus profile
//   read the fasta reference
//   read BAM file with MyPileupVisitor where MyPileupVisitor::Visit will be called on each position which populates infoPPos
//   creates threads using pthread_create using mainContaminationThread()
//   print to outlog and exit
//
// mainContaminationThread()
//   gets a contamination profile to use via locking the mutex "mutexQueue"
//   if no data is left to read, die
//   otherwise,  get an allele frequency file freqFileNameToUse  from queueFilesToprocess
//   call readMTAlleleFreq() on the freqFileNameToUse which will populate freqFromFile
//   Begin pre-computations that are static for each contamination rate
//        For every endogenous and contaminant pair
//            Compute prior on the pair
//            Compute the probability that the reads match them for the endogenous and contaminant probEndoVec


#include <api/BamConstants.h>
#include <api/BamMultiReader.h>
#include <utils/bamtools_fasta.h>
#include <utils/bamtools_options.h>
#include <utils/bamtools_pileup_engine.h>
#include <utils/bamtools_utilities.h>
#include <gzstream.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>
#include <queue>
#include <algorithm>

#include "utils.h"
#include "miscfunc.h"


using namespace BamTools;
using namespace std;


//! Chunk of code to check if a certain thread call failed
/*!
  This block is calls by the pthread

*/				
#define checkResults(string, val) {             \
 if (val) {                                     \
     cerr<<"Failed with "<<val<<" at "<<string<<endl;	\
   exit(1);                                     \
 }                                              \
}
 

typedef struct{
    long double logLike;
    string fname;
    long double contEst;
    long double contEstLow;
    long double contEstHigh;
} contRecord;

struct argsptcont {
    bool doubleDeam;
};

bool compContRecord (contRecord i, contRecord j) { return (i.logLike<j.logLike); }

char        offsetQual=33;

long double likeMatch[MAXMAPPINGQUAL];
long double likeMismatch[MAXMAPPINGQUAL];
long double likeMatchMQ[MAXMAPPINGQUAL][MAXMAPPINGQUAL];
long double likeMismatchMQ[MAXMAPPINGQUAL][MAXMAPPINGQUAL];

long double likeMatchProb[MAXMAPPINGQUAL];
long double likeMismatchProb[MAXMAPPINGQUAL];
long double likeMatchProbMQ[MAXMAPPINGQUAL][MAXMAPPINGQUAL];
long double likeMismatchProbMQ[MAXMAPPINGQUAL][MAXMAPPINGQUAL];

probSubstition illuminaErrorsProb;
vector<probSubstition> sub5p;
vector<probSubstition> sub3p;
probSubstition defaultSubMatch;

long double probMapping[MAXMAPPINGQUAL];
long double probMismapping[MAXMAPPINGQUAL];

string dnaAlphabet="ACGT";
long double stepContEst = 0.01;
long double topCont = 1.0;

#define MIN(a,b) (((a)<(b))?(a):(b))


queue<string>  queueFilesToprocess;
vector< vector<contaminationEstimate> > outputToPrint;

pthread_mutex_t  mutexQueue   = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t  mutexCounter = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t  mutexRank    = PTHREAD_MUTEX_INITIALIZER;


//GLOBALLY accessed
int sizeGenome               =0;
vector<positionInformation>  infoPPos;
map<int, PHREDgeno>          pos2phredgeno;
vector<int>                  posOfIndels;
map<unsigned int, int>       threadID2Rank;



vector<bool>  * definedSite;      //= new vector<bool>(sizeGenome+1,false); // if there is data
vector<bool>  * skipPositions;    //= new vector<bool>(sizeGenome+1,false); // if the site has such a low prior that it can be overlooked



//! Method to find positions to use where there is potential for contamination
/*!
  

*/				
void findPosToSkip(bool printPosToskip,bool verbose){
    queue<string>  queueFilesToRead = queueFilesToprocess;
    vector< map<int, alleleFrequency> > freqsFromFile;
    cerr<<"Finding positions to skip"<<endl;

    while(!queueFilesToRead.empty()){    
 	string freqFileNameToread = queueFilesToRead.front();
 	queueFilesToRead.pop();
	// cerr<<"file="<<freqFileNameToread<<"\t"<<queueFilesToproc.empty()<<endl;	

	map<int, alleleFrequency> freqToAdd;
	readMTAlleleFreq(freqFileNameToread,   freqToAdd);
	// cerr<<"file2="<<freqFileNameToread<<endl;	

	freqsFromFile.push_back(freqToAdd);
	// cerr<<"file3="<<freqFileNameToread<<"\t"<<queueFilesToprocess.empty()<<endl;	
    }


    vector<bool> definedFreqPos  = vector<bool>(sizeGenome+1,true); 

    //detect undefined positions
    for(unsigned int fileFreq=0;fileFreq<freqsFromFile.size();fileFreq++){ //   each frequency file found
	for(int i=0;i<sizeGenome;i++){
	    
	    double freqSum=0.0;
	    for(unsigned int nuc=0;nuc<4;nuc++){ //     b = endogenous

		freqSum+=freqsFromFile[fileFreq][ infoPPos[i].posAlign ].f[nuc];

	    }    
	    
	    if(freqSum<0.99){
		definedFreqPos[i]=false;
	    }
	}
    }

#ifdef DEBUGMTPOSSKIP
	for(int i=0;i<sizeGenome;i++){
	    if(!definedFreqPos[i])
		cout<<"empty\t"<<i<<"\t"<<infoPPos[i].posAlign<<"\t"<<definedFreqPos[i]<<endl;
	}
#endif

    // cout<<sizeGenome<<endl;
    
    for(int i=0;i<sizeGenome;i++){
	//for(int i=261;i<=262;i++){
	//	cout<<"Thread #"<<rankThread <<" test2 "<<i<<endl;

	//cout<<"pre i"<<i<<"\t"<<endl;
	// cout<<infoPPos[i].cov<<endl;
	// cout<<(pos2phredgeno.find( infoPPos[i].posAlign ) == pos2phredgeno.end())<<endl;
	if( (infoPPos[i].cov == 0)                                              || //no coverage
	    (pos2phredgeno.find( infoPPos[i].posAlign ) == pos2phredgeno.end()) || //not found in endogenous consensus, could be due to deletion
	    !definedFreqPos[i]                                                   ){  //frequency files disagree
	    definedSite->at(i) = false;
	    continue;
	}
	definedSite->at(i)     = true;

	
	//cout<<"Thread #"<<rankThread <<" test3 "<<i<<endl;
	bool hasPriorAboveThreshold=false;
	//computing prior
	//                 p(nuc1) is the prob. of endogenous                  *  p(contaminant)
	for(unsigned int fileFreq=0;fileFreq<freqsFromFile.size();fileFreq++){ //   each frequency file found
	    for(unsigned int nuc1=0;nuc1<4;nuc1++){ //     b = endogenous
		for(unsigned int nuc2=0;nuc2<4;nuc2++){//  c = contaminant
		    long double priortemp = (1-pos2phredgeno[ infoPPos[i].posAlign ].perror[nuc1]) * freqsFromFile[fileFreq][ infoPPos[i].posAlign ].f[nuc2];
		    //		    priorDiNuc->p[nuc1][nuc2] = priortemp;
		
#ifdef DEBUGPOS
		    if(i==DEBUGPOS){
			cout<<"file#"<<fileFreq<<"\tMIN i="<<i<<"\te="<<dnaAlphabet[nuc1]<<"\tc="<<dnaAlphabet[nuc2]<<"\tprior="<<priortemp<<"\tposal="<<infoPPos[i].posAlign<<"\tproblog="<<(1-pos2phredgeno[ infoPPos[i].posAlign ].perror[nuc1])<<"\tfreq="<< freqsFromFile[fileFreq][ infoPPos[i].posAlign ].f[nuc2]<<endl;
		    }
#endif
		    if( (nuc1 != nuc2) && //){
			(priortemp > MINPRIOR) ){ //this is to speed up computation and only look at sites that are likely to be contaminated

			// #ifdef DEBUGPOS
			// 		    cout<<"MIN i"<<i<<"\t"<<nuc1<<"\t"<<nuc2<<"\t"<<priortemp<<"\t"<<infoPPos[i].posAlign<<"\t"<<(1-pos2phredgeno[ infoPPos[i].posAlign ].perror[nuc1])<<"\t"<< freqFromFile[ infoPPos[i].posAlign ].f[nuc2]<<endl;
			// #endif

#ifdef DEBUGCONTPOS
			cout<<"MIN i="<<i<<"\te="<<nuc1<<"\tc="<<nuc2<<"\tprior="<<priortemp<<"\tposal="<<infoPPos[i].posAlign<<"\tproblog="<<(1-pos2phredgeno[ infoPPos[i].posAlign ].perror[nuc1])<<"\tfreq="<< freqsFromFile[fileFreq][ infoPPos[i].posAlign ].f[nuc2]<<endl;
#endif
			hasPriorAboveThreshold=true;
		    
		    }

		}
	    }
	}

	//	cout<<"Thread #"<<rankThread <<" test4 "<<i<<endl;
#ifdef DEBUGCONTPOS
	if(hasPriorAboveThreshold)
	    cout<<endl;
#endif

	//priorDiNucVec[ i ] = priorDiNuc;
	// if(i==3105){ exit(1); }
	skipPositions->at(i) = (!hasPriorAboveThreshold);

	if(printPosToskip && !skipPositions->at(i))
	    cerr<<"Keeping position "<<(i+1)<<endl;

	if(verbose)
	    cerr<<(skipPositions->at(i)?"Skipping":"Keeping")<<" position "<<(i+1)<<endl;
	    //cout<<i<<"\t"<<skipPositions->at(i)<<endl;
	// if(i%1000==0)
	//     skipPositions->at(i) = (false); //TO REMOVE

	//cout<<i<<"\t"<<skipPositions->at(i)<<endl;
    }//end for each pos in the genome

    cerr<<"done"<<endl;


}


//! Method called for each thread to compute contamination likelihood
/*!
  

*/				
void *mainContaminationThread(void * argc){


    int   rc;
    // int   stackIndex;
    string freqFileNameToUse;
    int rankThread=0;

    rc = pthread_mutex_lock(&mutexRank);
    checkResults("pthread_mutex_lock()\n", rc);

    threadID2Rank[*(int *)pthread_self()]  = threadID2Rank.size()+1;
    rankThread = threadID2Rank[*(int *)pthread_self()];

    
    rc = pthread_mutex_unlock(&mutexRank);
    checkResults("pthread_mutex_unlock()\n", rc);

    bool doubleDeam = ((struct argsptcont*)argc)->doubleDeam;


 checkqueue:    
    // stackIndex=-1;
    //check stack

    
    rc = pthread_mutex_lock(&mutexQueue);
    checkResults("pthread_mutex_lock()\n", rc);


    bool foundData=false;
    
    cerr<<"Thread #"<<rankThread <<" started and is requesting data"<<endl;

    //cerr<<"Thread #"<<(unsigned int)pthread_self() <<" started "<<endl;

	
    // cout<<"Thread "<<(unsigned int)pthread_self()<<" taking mutex queue "<<endl;
    if(!queueFilesToprocess.empty()){    
 	foundData=true;
 	freqFileNameToUse = queueFilesToprocess.front();
 	queueFilesToprocess.pop();
 	cerr<<"Thread #"<<rankThread<<" is reading "<<freqFileNameToUse<<endl;
    }

    
  

    if(!foundData){
 	//if(doneReading){

	rc = pthread_mutex_unlock(&mutexQueue);
	checkResults("pthread_mutex_unlock()\n", rc);

	cerr<<"Thread #"<<rankThread<<" is done"<<endl;
	return NULL;	
 	// }else{
 	//     sleep(1);
 	//     goto checkqueue;	   
 	// }
    }else{

    //release stack
	rc = pthread_mutex_unlock(&mutexQueue);
	checkResults("pthread_mutex_unlock()\n", rc);
    }
    //////////////////////////////////////////////////////////////
    //                BEGIN COMPUTATION                         //
    //////////////////////////////////////////////////////////////



    vector<contaminationEstimate> toAddToLog;
    long double contaminationRate=0.0;

    map<int, alleleFrequency> freqFromFile;
    readMTAlleleFreq(freqFileNameToUse,	freqFromFile);


    //////////////////////
    //pre-computations ///
    //////////////////////

    cerr<<"Thread #"<<rankThread <<" started  pre-computations"<<endl;

    // map<int, diNucleotideProb> priorDiNucVec;
    // map<int, vector<diNucleotideProb> > probConsVec;
    // map<int, vector<diNucleotideProb> > probContVec;


    vector<         diNucleotideProb * >      priorDiNucVec;
    vector< vector< diNucleotideProb * > * >  probEndoVec;
    vector< vector< diNucleotideProb * > * >  probContVec;
    vector< vector< bool               > * >  strandVec;

    priorDiNucVec.resize(sizeGenome+1);
    probEndoVec.resize(sizeGenome+1);
    probContVec.resize(sizeGenome+1);
    strandVec.resize(sizeGenome+1);

    // vector<bool>  * definedSite      = new vector<bool>(sizeGenome+1,false); // if there is data
    // vector<bool>  * skipPositions    = new vector<bool>(sizeGenome+1,false); // if the site has such a low prior that it can be overlooked


    for(int i=0;i<sizeGenome;i++){
	//for(int i=261;i<=262;i++){
	//	cout<<"Thread #"<<rankThread <<" test2 "<<i<<endl;

	if(    !definedSite->at(i) ){
	    continue;
	}
	if(    skipPositions->at(i) ){
	    continue;
	}
	
	diNucleotideProb * priorDiNuc = 0;
	try{
	    priorDiNuc = new diNucleotideProb;
	}catch( char * str ) {
	    cout << "Exception raised: " << str << endl;
	}
	//cout<<"Thread #"<<rankThread <<" test3 "<<i<<endl;
	//bool hasPriorAboveThreshold=false;
	//computing prior
	//                 p(nuc1) is the prob. of endogenous                  *  p(contaminant)
	for(unsigned int nuc1=0;nuc1<4;nuc1++){ //     b = endogenous
	    for(unsigned int nuc2=0;nuc2<4;nuc2++){//  c = contaminant
		long double priortemp = (1-pos2phredgeno[ infoPPos[i].posAlign ].perror[nuc1]) * freqFromFile[ infoPPos[i].posAlign ].f[nuc2];
		priorDiNuc->p[nuc1][nuc2] = priortemp;
		
		// if(i==3105){
		//cout<<"i"<<i<<"\t"<<nuc1<<"\t"<<nuc2<<"\t"<<priortemp<<"\t"<<(1-pos2phredgeno[ infoPPos[i].posAlign ].perror[nuc1])<<"\t"<<freqFromFile[ infoPPos[i].posAlign ].f[nuc2]<<"\t"<<infoPPos[i].posAlign<<endl;			    			
		// }
#ifdef DEBUGPOS
		if(i==DEBUGPOS){
		    cout<<"MIN i="<<i<<"\te="<<dnaAlphabet[nuc1]<<"\tc="<<dnaAlphabet[nuc2]<<"\tprior="<<priortemp<<"\tposal="<<infoPPos[i].posAlign<<"\tproblog="<<(1-pos2phredgeno[ infoPPos[i].posAlign ].perror[nuc1])<<"\tfreq="<< freqFromFile[ infoPPos[i].posAlign ].f[nuc2]<<endl;
		}
#endif



	    }//end nuc2
	}//end nuc1

	//	cout<<"Thread #"<<rankThread <<" test4 "<<i<<endl;
// #ifdef DEBUGCONTPOS
// 	if(hasPriorAboveThreshold)
// 	    cout<<endl;
// #endif

	priorDiNucVec[ i ] = priorDiNuc;
	// if(i==3105){ exit(1); }
	// skipPositions->at(i) = (!hasPriorAboveThreshold);



	vector<diNucleotideProb *> * probEndoVecToAdd = new vector<diNucleotideProb *>();
	vector<diNucleotideProb *> * probContVecToAdd = new vector<diNucleotideProb *>();
	vector<bool              > * strandVecToAdd   = new vector<bool>();



	// continue;
	for(unsigned int k=0;k<infoPPos[i].readsVec.size();k++){ //for every read at that position
	    //cout<<"Thread #"<<rankThread <<" test5\t"<<i<<"\t"<<k<<endl;
	    diNucleotideProb * probEndoDinuc=0;
	    diNucleotideProb * probContDinuc=0;
	    try {
		probEndoDinuc = new diNucleotideProb;
		probContDinuc = new diNucleotideProb;
	    }catch( char * str ) {
		cout << "Exception raised: " << str << endl;
	    }

	    //cout<<"Thread #"<<rankThread <<" test6\t"<<i<<"\t"<<k<<endl;
	    int baseIndex = baseResolved2int(infoPPos[i].readsVec[k].base);
	    int qual      = infoPPos[i].readsVec[k].qual;
	    int dist5p    = infoPPos[i].readsVec[k].dist5p;
	    int dist3p    = infoPPos[i].readsVec[k].dist3p;

	    //int mapq      = infoPPos[i].readsVec[k].mapq;

	    probSubstition * probSubMatchToUseEndo = &defaultSubMatch ;
	    probSubstition * probSubMatchToUseCont = &defaultSubMatch ; //leave it that way, we only allow deamination for the endogenous only
		    

	    //consider deamination to be possible for the endogenous only
	    if(dist5p <= (int(sub5p.size()) -1)){
		probSubMatchToUseEndo = &sub5p[ dist5p ];			
		if(doubleDeam){ //copy G->A from last 3' position
		    probSubMatchToUseEndo->s[  8 ] = sub3p[ (int(sub3p.size()) -1) ].s[  8 ];
		    probSubMatchToUseEndo->s[  9 ] = sub3p[ (int(sub3p.size()) -1) ].s[  9 ];
		    probSubMatchToUseEndo->s[ 10 ] = sub3p[ (int(sub3p.size()) -1) ].s[ 10 ];
		    probSubMatchToUseEndo->s[ 11 ] = sub3p[ (int(sub3p.size()) -1) ].s[ 11 ];
		}
	    }
		    
	    if(dist3p <= (int(sub3p.size()) -1)){
		probSubMatchToUseEndo = &sub3p[ dist3p ];
		if(doubleDeam){ //copy C->T from last 5' position
		    probSubMatchToUseEndo->s[ 4 ] = sub5p[ (int(sub5p.size()) -1) ].s[ 4 ];
		    probSubMatchToUseEndo->s[ 5 ] = sub5p[ (int(sub5p.size()) -1) ].s[ 5 ];
		    probSubMatchToUseEndo->s[ 6 ] = sub5p[ (int(sub5p.size()) -1) ].s[ 6 ];
		    probSubMatchToUseEndo->s[ 7 ] = sub5p[ (int(sub5p.size()) -1) ].s[ 7 ];
		}
	    }
		    
	    //we have substitution probabilities for both... take the closest
	    if(dist5p <= (int(sub5p.size()) -1) &&
	       dist3p <= (int(sub3p.size()) -1) ){
			
		if(dist5p < dist3p){
		    probSubMatchToUseEndo = &sub5p[ dist5p ];
		    if(doubleDeam){ //copy G->A from last 3' position
			probSubMatchToUseEndo->s[  8 ] = sub3p[ (int(sub3p.size()) -1) ].s[  8 ];
			probSubMatchToUseEndo->s[  9 ] = sub3p[ (int(sub3p.size()) -1) ].s[  9 ];
			probSubMatchToUseEndo->s[ 10 ] = sub3p[ (int(sub3p.size()) -1) ].s[ 10 ];
			probSubMatchToUseEndo->s[ 11 ] = sub3p[ (int(sub3p.size()) -1) ].s[ 11 ];
		    }
		}else{
		    probSubMatchToUseEndo = &sub3p[ dist3p ];
		    if(doubleDeam){ //copy C->T from last 5' position
			probSubMatchToUseEndo->s[ 4 ] = sub5p[ (int(sub5p.size()) -1) ].s[ 4 ];
			probSubMatchToUseEndo->s[ 5 ] = sub5p[ (int(sub5p.size()) -1) ].s[ 5 ];
			probSubMatchToUseEndo->s[ 6 ] = sub5p[ (int(sub5p.size()) -1) ].s[ 6 ];
			probSubMatchToUseEndo->s[ 7 ] = sub5p[ (int(sub5p.size()) -1) ].s[ 7 ];
		    }
		}
			
	    }
	    
#ifdef DEBUGSUBPATTERN
	    if(i==DEBUGPOS){

		for(int sub=0;sub<16;sub++){ //     b = endogenous
		    cerr<<"read#"<<k<<"\t"<<dist5p<<"\t"<<dist3p<<"\tsub. damage P["<<sub<<"] = "<<probSubMatchToUseEndo->s[sub]<<endl;
		}

	    }	    
#endif

	    //iterate over each possible endogenous and contaminant base
	    for(int nuc1=0;nuc1<4;nuc1++){ //     b = endogenous
		for(int nuc2=0;nuc2<4;nuc2++){//  c = contaminant
		    //skip when contaminant is the endogenous
		    // if(nuc1 == nuc2)
		    // 	continue;


		    ////////////////////
		    // Endogenous base /
		    ////////////////////

		    //                  model*4  + obs
		    int dinucIndexEndo = nuc1*4+baseIndex;

		    if(infoPPos[i].readsVec[k].isReversed)
			dinucIndexEndo = (3-nuc1)*4+(3-baseIndex);
#ifdef DEBUGPOSEACHREAD

		    if(i==DEBUGPOS){
			if(infoPPos[i].readsVec[k].isReversed){
			    dinucIndexEndo = (3-nuc1)*4+(3-baseIndex);
			    cerr<<"read#"<<k<<"\t"<<dnaAlphabet[nuc1]<<"\t"<<dnaAlphabet[baseIndex]<<"\tidxDiNuc:"<<dinucIndexEndo<<"\t"<<probSubMatchToUseEndo->s[dinucIndexEndo]<<"\t"<<infoPPos[i].readsVec[k].isReversed<<endl;
			}else{
			    cerr<<"read#"<<k<<"\t"<<dnaAlphabet[nuc1]<<"\t"<<dnaAlphabet[baseIndex]<<"\tidxDiNuc:"<<dinucIndexEndo<<"\t"<<probSubMatchToUseEndo->s[dinucIndexEndo]<<"\t"<<infoPPos[i].readsVec[k].isReversed<<endl;
			}

		    }
		    // if(i==DEBUGPOS){
		    // 	cout<<dist5p<<"\t"<<dist3p<<endl;			
		    // }
#endif
		    
		    //                        (1-e)           *  p(sub|1-e)                                  + (e) *  p(sub|1-e)
		    long double probEndo=likeMatchProb[qual]  * (probSubMatchToUseEndo->s[dinucIndexEndo] )  + (1.0 - likeMatchProb[qual])*(illuminaErrorsProb.s[dinucIndexEndo]);
			
		    ///////////////////
		    //Contaminant base/
		    ///////////////////

		    //                  model*4  + obs
		    int dinucIndexCont = nuc2*4+baseIndex;
		    if(infoPPos[i].readsVec[k].isReversed)
			dinucIndexCont = (3-nuc2)*4+(3-baseIndex);


		    //                        (1-e)           *  p(sub|1-e)                                  + (e) *  p(sub|1-e)
		    long double probCont=likeMatchProb[qual]  * (probSubMatchToUseCont->s[dinucIndexCont] )  + (1.0 - likeMatchProb[qual])*(illuminaErrorsProb.s[dinucIndexCont]);


		    probEndoDinuc->p[nuc1][nuc2] = probEndo;
		    probContDinuc->p[nuc1][nuc2] = probCont;
		  
		}
	    }//end for each di-nucleotide
		
	    
	    probEndoVecToAdd->push_back(probEndoDinuc);
	    probContVecToAdd->push_back(probContDinuc);
	    strandVecToAdd->push_back(infoPPos[i].readsVec[k].isReversed);

	    
	} //end for each read at that position
	
	//	cout<<"Thread #"<<rankThread <<" test3\t"<<probEndoVec.size()<<"\t"<<probContVec.size()<<endl;
	//cout<<"adding vector at pos "<<i<<endl;
	probEndoVec[i] = probEndoVecToAdd;
	probContVec[i] = probContVecToAdd;
	strandVec[i]   = strandVecToAdd;

	//	cout<<"Thread #"<<rankThread <<" test4"<<endl;

    } //end for each position in the genome
    cerr<<"Thread #"<<rankThread <<" is done with  pre-computations"<<endl;

#ifdef DEBUGCONTRATE 
    for(contaminationRate=DEBUGCONTRATE;contaminationRate<DEBUGCONTRATE+stepContEst;contaminationRate+=stepContEst){
#else        
    for(contaminationRate=0.0;contaminationRate<topCont;contaminationRate+=stepContEst){
#endif
	long double logLike=0.0;
    

#ifdef DEBUGPOSPRINTSINGLEPOS
 	for(int i=DEBUGPOS;i<=DEBUGPOS;i++){
#else
        for(int i=0;i<sizeGenome;i++){
#endif
	    if(    !definedSite->at(i) ){//site is not defined
		continue;
	    }
	    if(    skipPositions->at(i) ){//position has a tiny prior on contamination and can be safely skipped
		continue;
	    }
	    
	    if((infoPPos[i].cov == 0) || //no coverage
	       (pos2phredgeno.find( infoPPos[i].posAlign ) == pos2phredgeno.end()) ){  //not found in endogenous, could be due to deletion 
		// skipPositions->at(i) ){ //position has a tiny prior on contamination and can be safely skipped
		continue;
	    }

#ifdef DEBUGPOSEACHREAD
	    cerr<<"after i="<<i<<"\t"<<infoPPos[i].posAlign<<endl;
#endif

	    // continue;
	    for(unsigned int k=0;k<infoPPos[i].readsVec.size();k++){ //for every read at that position
	    //for(unsigned int k=0;k<1;k++){ //for every read at that position
		long double mappedProb  =0.0; //initialized
		long double misappedProb=0.25;
	   
		int mapq      = infoPPos[i].readsVec[k].mapq;

		//iterate over each possible endogenous and contaminant base
		for(int nuc1=0;nuc1<4;nuc1++){ //     b = endogenous
		    for(int nuc2=0;nuc2<4;nuc2++){//  c = contaminant
			//skip when contaminant is the endogenous
			// if(nuc1 == nuc2)
			//     continue;

			
		  	long double probForDinuc = priorDiNucVec[i]->p[nuc1][nuc2]*
			    ( ( (1.0-contaminationRate) * probEndoVec[i]->at(k)->p[nuc1][nuc2]) +
			      ( (contaminationRate)     * probContVec[i]->at(k)->p[nuc1][nuc2]   ) ) ;
		  	// long double probForDinuc = priorDiNucVec[i]->p[nuc1][nuc2]*
			//     ( ( (1.0-contaminationRate) * probEndoVec[i]->at(k)->p[nuc1][nuc2]) +
			//       ( (contaminationRate)     * probContVec[i]->at(k)->p[nuc1][nuc2]   ) ) ;
		     
			mappedProb+=probForDinuc;

#ifdef DEBUGPOSEACHREAD
			if(i==DEBUGPOS){
			    if(0){
				cout<<endl<<"k="<< k <<"\ti="<<i<<"\te="<<dnaAlphabet[nuc1]<<"\tc="<<dnaAlphabet[nuc2]<<endl;			    
				cout<<"c               "<<contaminationRate<<endl;
				cout<<"base read       "<<infoPPos[i].readsVec[k].base<<endl;
				cout<<"base qual       "<<infoPPos[i].readsVec[k].qual<<endl;
				cout<<"prior           "<<priorDiNucVec[i]->p[nuc1][nuc2]<<endl;
				cout<<"probCont        "<<probContVec[i]->at(k)->p[nuc1][nuc2]<<endl;
				cout<<"probEndo        "<<probEndoVec[i]->at(k)->p[nuc1][nuc2]<<endl;
				cout<<"probForDinuc    "<<probForDinuc<<endl;
			    }else{
				//if(priorDiNucVec[i].p[nuc1][nuc2]>0.9){
				cerr.precision(15);
				cerr<<"read#"<< k <<"\ti="<<i<<"\te="<<dnaAlphabet[nuc1]<<"\tc="<<dnaAlphabet[nuc2]<<"\tbase="<<infoPPos[i].readsVec[k].base<<"\tqual=" <<infoPPos[i].readsVec[k].qual<<"\tprior="<<priorDiNucVec[i]->p[nuc1][nuc2]<<"\tp(e)="<<probEndoVec[i]->at(k)->p[nuc1][nuc2]<<"\tp(c)="<<probContVec[i]->at(k)->p[nuc1][nuc2]<<"\tp="<<probForDinuc<<"\tp(ma)="<<mappedProb<<"\t"<<mapq<<"\te="<<dnaAlphabet[nuc1]<<"\tc="<<dnaAlphabet[nuc2]<<"\t"<<strandVec[i]->at(k)<<endl;
				//}
			    }
			}
#endif

		    }
		}//end for each di-nucleotide

		
		//        m = prob. of mismapping
		//                  (1-m)     * p(data|mapped) +           m         * p(data|mismapped)	
		long double pread = (probMapping[mapq]*(mappedProb) + probMismapping[mapq]*(misappedProb) );
		// cout<<k<<"\t"<<pread<<endl;
		    
		//cout<<pread<<endl;
		// cout<<probContVec[i][k].p[nuc1][nuc2]<<endl;
		// cout<<probEndoVec[i][k].p[nuc1][nuc2]<<endl;
		
		//exit(1);
		
		logLike += log  (   pread   );

#ifdef DEBUGPOSLOGLIKE
		if(i==DEBUGPOS){
		    cout<<freqFileNameToUse<<"\tread#"<<k<<"\tpread\t"<<pread<<"\tloglike="<< logLike <<"\tcont="<<contaminationRate<<endl;		
		    //cout<<"-------------"<<endl;
		}
#endif
	    }//end for each read at that position

	    // cout<<"logLike "<<i<<"\t"<<logLike<<endl;
#ifdef DEBUGPOS
	    // if(i==DEBUGPOS){
	    // 	exit(1);
	    // }
#endif
	} //end for each position in the genome
	// exit(1);

	//cout.precision(15);
#ifdef DEBUGPOS
	cerr<<"freqFile "<<freqFileNameToUse<<"\t"<<contaminationRate<<"\t"<<logLike<<endl;
#endif
	contaminationEstimate cse;
	cse.filename           = freqFileNameToUse ;
	cse.contaminationRate  = contaminationRate ;
	cse.logLike            = logLike ;

	//toAddToLog+=(""+freqFileNameToUse+"\t"+stringify(contaminationRate)+"\t"+stringify(logLike)+"\n");
	toAddToLog.push_back(cse);
	
    }//end each contaminationRate
    // }//end for each cont freq file
	


    cerr<<"Thread #"<<rankThread <<" is done with computations"<<endl;

    //////////////////////////////////////////////////////////////
    //                END   COMPUTATION                         //
    //////////////////////////////////////////////////////////////
	
    // cout<<"done deleting"<<endl;
    for(unsigned i=0;i<priorDiNucVec.size();i++){

	//cout<<"i1 delete"<<i<<endl;
	if(i<definedSite->size() &&
	   !definedSite->at(i) ){
	    continue;
	}

	if(i<skipPositions->size() && 
	   skipPositions->at(i) ){
	    continue;
	}

	///cout<<"i2 delete"<<i<<endl;
	delete priorDiNucVec[i];

	for(unsigned j=0;j<probEndoVec[i]->size();j++){
	    // cout<<"j "<<j<<endl;
	    delete probEndoVec[i]->at(j);
	    delete probContVec[i]->at(j);
	}


	delete probEndoVec[i];
	delete probContVec[i];	    
    }
    

    //COUNTERS
    rc = pthread_mutex_lock(&mutexCounter);
    checkResults("pthread_mutex_lock()\n", rc);

    outputToPrint.push_back(toAddToLog);
    

    rc = pthread_mutex_unlock(&mutexCounter);
    checkResults("pthread_mutex_unlock()\n", rc);

    cerr<<"Thread #"<<rankThread <<" is re-starting"<<endl;

    goto checkqueue;	   


    

    
    cerr<<"Thread "<<rankThread<<" ended "<<endl;
    return NULL;

}






//unsigned int posFound;
//vector<positionInfo> * positionsInfoFound;
// map<string, map<unsigned int,contaminationInfo> > contFreq;

class MyPileupVisitor : public PileupVisitor {
  
    public:
    MyPileupVisitor(const RefVector& references, Fasta * fastaReference,vector<positionInformation> * infoPPos,int sizeGenome,bool ignoreMQ)
            : PileupVisitor()
	    , m_references(references)
	    , m_fastaReference(fastaReference)
	    , m_infoPPos(infoPPos)
	    , sizeGenome(sizeGenome)
	    , ignoreMQ(ignoreMQ)
        { 

	}
        ~MyPileupVisitor(void) { }
  
    // PileupVisitor interface implementation
    public:
	// prints coverage results ( tab-delimited )
        void Visit(const PileupPosition& pileupData) {

	    char referenceBase = 'N';
	    
	    unsigned int posAlign = pileupData.Position+1;
	    int posVector=int(pileupData.Position)%sizeGenome;
	    // cout<<posVector<<endl;
	    // cout<<endl<<"pos = "<<posAlign<<"\t"<<posVector;

	    if( (posAlign%100) == 0){
		//cerr<<"pos  = "<<posAlign<<endl;
	    }
	    //for some reason, we have to do -1 on the .Position
	    if ( !m_fastaReference->GetBase(pileupData.RefId, posAlign-1, referenceBase ) ) {
		cerr << "bamtools convert ERROR: pileup conversion - could not read reference base from FASTA file" << endl;
		exit(1);
	    }

	    
	    //Like in the endogenous calling, there are 3 possibilities, except here we only use 3) :
	    //1) There is a insertion in the sample (or deletion in the reference)
	    //2) There is a deletion in the sample (or insertion in the reference)
	    //3) There is potentially a variation of a single nucleotide
	    	    

	    


	    //insertion in the reads/deletion in the reference, skipping
	    for(unsigned int i=0;i<pileupData.PileupAlignments.size();i++){		
		if( !pileupData.PileupAlignments[i].IsCurrentDeletion &&
		    pileupData.PileupAlignments[i].IsNextInsertion &&
		    (pileupData.PileupAlignments[i].InsertionLength>0)){
		    continue;
		}
	    }


	    //deletion in the reads/insertion in the reference, skipping
	    unsigned int numReads=0;
	    for(unsigned int i=0;i<pileupData.PileupAlignments.size();i++){		
		if( pileupData.PileupAlignments[i].IsCurrentDeletion &&
		    pileupData.PileupAlignments[i].IsNextInsertion &&
		    (pileupData.PileupAlignments[i].InsertionLength == 0)){
		    continue;
		    //cout<<"del"<<endl;
		    //m_infoPPos->at(posVector).numDel++;
		}
	    }
	    

	    // typedef struct { 
	    //     vector<singleRead> readsVec;
	    //     int cov;
	    // } positionInformation;


	    for(unsigned int i=0;i<pileupData.PileupAlignments.size();i++){		
		// cerr<<pileupData.PileupAlignments[i].Alignment.Name<<endl;	    		    
		//skip deletion in the reads/insertion in the reference
		if( pileupData.PileupAlignments[i].IsCurrentDeletion &&
		    pileupData.PileupAlignments[i].IsNextInsertion ){
		    continue;
		}
		   
		char b   = pileupData.PileupAlignments[i].Alignment.QueryBases[pileupData.PileupAlignments[i].PositionInAlignment];
		char q   = pileupData.PileupAlignments[i].Alignment.Qualities[pileupData.PileupAlignments[i].PositionInAlignment]-offsetQual;
		int  m   = int(pileupData.PileupAlignments[i].Alignment.MapQuality);
		
		if(b == 'N') //uninformative anyway...
		    continue;

		singleRead sr;
		sr.base=b;
		sr.qual=int(q);
		sr.mapq=int(m);

		sr.dist5p=-1;
		sr.dist3p=-1;
		sr.isReversed= pileupData.PileupAlignments[i].Alignment.IsReverseStrand();
		
		if( sr.isReversed  ){
		    sr.dist5p = pileupData.PileupAlignments[i].Alignment.QueryBases.size() - pileupData.PileupAlignments[i].PositionInAlignment-1;
		    sr.dist3p = pileupData.PileupAlignments[i].PositionInAlignment;
		}else{
		    sr.dist5p = pileupData.PileupAlignments[i].PositionInAlignment;
		    sr.dist3p = pileupData.PileupAlignments[i].Alignment.QueryBases.size() - pileupData.PileupAlignments[i].PositionInAlignment-1;
		}

		m_infoPPos->at(posVector).cov++;		    

		//Add mapq
		m_infoPPos->at(posVector).mapqAvg += pow(10.0, (long double)(pileupData.PileupAlignments[i].Alignment.MapQuality)/ (long double)(-10.0) );
		m_infoPPos->at(posVector).readsVec.push_back(sr);
		
		// numReads++;
		if(numReads >= MAXCOV){
		    break;
		}

	    }//end each read
	    
	    // m_infoPPos->at(posVector).mapqAvg  = m_infoPPos->at(posVector).mapqAvg/double(m_infoPPos->at(posVector).cov);
	    // m_infoPPos->at(posVector).mapqAvg  = -10.0*( log( m_infoPPos->at(posVector).mapqAvg )/log(10.0) );
	    // m_infoPPos->at(posVector).refBase  = referenceBase;
	    //cout<<"pv = "<<posVector<<"\tpa ="<<posAlign<<endl;
	    m_infoPPos->at(posVector).posAlign = int(pileupData.Position+1)%sizeGenome; //posVector;//posAlign;
	    // m_infoPPos->at(posVector).skipPosition = false;

        }
        
    private:
    RefVector m_references;
    Fasta * m_fastaReference;
    vector<positionInformation> * m_infoPPos;
    int sizeGenome;
    bool ignoreMQ;
    //        ostream*  m_out;
};










int main (int argc, char *argv[]) {





    int numberOfThreads=1;
    // string outSeq  = "/dev/stdout";
    string outLog  = "/dev/stdout";
    // string nameMT  = "MT";
    //string fileFreq = "alleleFreqMT/1000g/freqHumans.dat";

    // ofstream outSeqFP ;
    ofstream outLogFP;
    // int minQual=0;
    bool ignoreMQ=false;
    string line;
    

    /////////////////////////////////////////
    //        INITIALIZE SCORES            //
    /////////////////////////////////////////


    for(int i=0;i<2;i++){
        likeMatch[i]        = log1p(    -pow(10.0,2.0/-10.0) )    /log(10);         
        likeMismatch[i]     = log  (     pow(10.0,2.0/-10.0)/3.0 )/log(10);

	likeMatchProb[i]           = 1.0-pow(10.0,2.0/-10.0) ;
        likeMismatchProb[i]        =     pow(10.0,2.0/-10.0)/3.0 ;
    }


    //Computing for quality scores 2 and up
    for(int i=2;i<MAXMAPPINGQUAL;i++){
        likeMatch[i]        = log1p(    -pow(10.0,i/-10.0) )     /log(10);          
        likeMismatch[i]     = log  (     pow(10.0,i/-10.0)/3.0  )/log(10);

        likeMatchProb[i]           = 1.0-pow(10.0,i/-10.0);
        likeMismatchProb[i]        =     pow(10.0,i/-10.0)/3.0;
    }


    //Adding mismapping probability
    for(int m=0;m<MAXMAPPINGQUAL;m++){

	//m = prob of mismapping
	long double incorrectMappingProb   =     pow(10.0,m/-10.0); //m
	long double correctMappingProb     = 1.0-pow(10.0,m/-10.0); //1-m


	probMapping[m]    = correctMappingProb;    //1-m
	probMismapping[m] = incorrectMappingProb;  //m

#ifdef DEBUG1
	cerr<<"m\t"<<m<<"\t"<<incorrectMappingProb<<"\t"<<correctMappingProb<<endl;
#endif
	
    	for(int i=0;i<2;i++){
    	    likeMatchMQ[m][i]           = log(  correctMappingProb*(1.0-pow(10.0,2.0/-10.0)    ) + incorrectMappingProb/4.0   )/log(10);         
    	    likeMismatchMQ[m][i]        = log(  correctMappingProb*(    pow(10.0,2.0/-10.0)/3.0) + incorrectMappingProb/4.0   )/log(10);
    	    likeMatchProbMQ[m][i]       = correctMappingProb*(1.0-pow(10.0,2.0/-10.0)    ) + incorrectMappingProb/4.0;
    	    likeMismatchProbMQ[m][i]    = correctMappingProb*(    pow(10.0,2.0/-10.0)/3.0) + incorrectMappingProb/4.0;
    	}


    	//Computing for quality scores 2 and up
    	for(int i=2;i<MAXMAPPINGQUAL;i++){
	    //  (1-m)(1-e) + m/4  = 1-m-e+me +m/4  = 1+3m/4-e+me
    	    likeMatchMQ[m][i]         = log(  correctMappingProb*(1.0-pow(10.0,i/-10.0)    ) + incorrectMappingProb/4.0    )/log(10);    
	    //  (1-m)(e/3) + m/4  = e/3 -me/3 + m/4
    	    likeMismatchMQ[m][i]      = log(  correctMappingProb*(    pow(10.0,i/-10.0)/3.0) + incorrectMappingProb/4.0    )/log(10);    
	    
    	    likeMatchProbMQ[m][i]           = correctMappingProb*(1.0-pow(10.0,i/-10.0)    ) + incorrectMappingProb/4.0;
    	    likeMismatchProbMQ[m][i]        = correctMappingProb*(    pow(10.0,i/-10.0)/3.0) + incorrectMappingProb/4.0;
    	}


#ifdef DEBUG1
    	for(int i=0;i<MAXMAPPINGQUAL;i++){
	    cerr<<"m\t"<<m<<"\t"<<i<<"\t"<<likeMatchMQ[m][i]<<"\t"<<likeMismatchMQ[m][i]<<"\t"<<likeMatchProbMQ[m][i]<<"\t"<<likeMismatchProbMQ[m][i]<<endl;
	}
#endif

    }


    /////////////////////////////////////////
    //        PARSING ARGUMENTS            //
    /////////////////////////////////////////

    //    return 1;
    string errFile    = getCWD(argv[0])+"../share/schmutzi/illuminaProf/error.prof";
    string deam5pfreq = getCWD(argv[0])+"../share/schmutzi/deaminationProfile/none.prof";
    string deam3pfreq = getCWD(argv[0])+"../share/schmutzi/deaminationProfile/none.prof";
    bool verbose      = false;
    bool summary      = false;
    bool printPosToskip = false;
    bool doubleDeam         = false;
    const string usage=string("\t"+string(argv[0])+
			      " [options]  [endogenous consensus file] [reference fasta file] [bam file] [cont file1] [cont file2] ..."+"\n\n"+
			      
			      "\n\tOutput options:\n"+	
			      "\t\t"+"-s" +"\t\t\t\t"+"Print summary, print the MAP value for each contaminant by decreasing order of likelihood"+"\n"			      "\t\t"+"-v" +"\t\t\t\t"+"Verbose output"+"\n"+
			      "\t\t"+"-o  [output log]" +"\t\t"+"Output log (default: stdout)"+"\n"+


			      //"\t\t"+"-log  [log file]" +"\t\t"+"Output log  (default: stderr)"+"\n"+
			      // "\t\t"+"-name [name]" +"\t\t\t"  +"Name  (default "+nameMT+") "+"\n"+
			      // "\t\t"+"-qual [minimum quality]" +"\t\t"  +"Filter bases with quality less than this  (default "+stringify(minQual)+") "+"\n"+
			       // "\n\tOutput options:\n"+				      
			      "\n\tDeamination options:\n"+				      

			      "\t\t"+"-deam5p [.prof file]" +"\t\t"+"5p deamination frequency (default: "+deam5pfreq+")"+"\n"+
			      "\t\t"+"-deam3p [.prof file]" +"\t\t"+"3p deamination frequency (default: "+deam3pfreq+")"+"\n"+
			      "\t\t"+"-double" +"\t\t"+"Consider residual C->T even close to the 3' endand residual G->A close to 5' (default: "+booleanAsString(doubleDeam)+")"+"\n"+
			      // // "\t\t"+"-deam5" +"\t\t"+"5p deamination frequency"+"\n"+
			      // // "\t\t"+"-deam3" +"\t\t"+"3p deamination frequency"+"\n"+
			      "\n\tComputation options:\n"+	
			      "\t\t"+"-maxc" +"\t\t\t\t"+"Maximum contamination value to check (default: "+stringify(topCont)+")"+"\n"+

			      "\t\t"+"-step" +"\t\t\t\t"+"Step for reporting the contamination estimate (default: "+stringify(stepContEst)+")"+"\n"+
			      "\t\t"+"-err" +"\t\t\t\t"+"Illumina error profile (default: "+errFile+")"+"\n"+
			      "\t\t"+"-nomq" +"\t\t\t\t"+"Ignore mapping quality (default: "+booleanAsString(ignoreMQ)+")"+"\n"+
			      "\t\t"+"--phred64" +"\t\t\t"+"Use PHRED 64 as the offset for QC scores (default : PHRED33)"+"\n"+

			      "\n\tMisc. options:\n"+	
			      "\t\t"+"-t" +"\t\t\t\t"+"Number of cores to use (default: "+stringify(numberOfThreads)+")"+"\n"+
			      "\t\t"+"-p" +"\t\t\t\t"+"Print considered positions  (default: "+booleanAsString(printPosToskip)+")"+"\n"+
			      
			      // "\t\t"+"-cont [file1,file2,...]" +"\t\t"+"Contamination allele frequency (default : "+fileFreq+")"+"\n"+ 
			      // "\n\tReference options:\n"+	
			      // "\t\t"+"-l [length]" +"\t\t\t"+"Actual length of the genome used for"+"\n"+
			      // "\t\t"+"  " +"\t\t\t\t"+"the reference as been wrapped around"+"\n"+
			      // "\t\t"+"  " +"\t\t\t\t"+"by default, the length of the genome will be used "+"\n"+
			      
			      "");
			      
    if( (argc== 1) ||
	(argc== 2 && string(argv[1]) == "-h") ||
	(argc== 2 && string(argv[1]) == "-help") ||
	(argc== 2 && string(argv[1]) == "--help") ){
	cout<<"Usage:"<<endl;
	cout<<""<<endl;
	cout<<usage<<endl;
	return 1;
    }


    int lastOpt=1;
    for(int i=1;i<(argc-3);i++){ //all but the last 3 args

	if(string(argv[i]) == "--phred64"  ){
	    offsetQual=64;
	    continue;
	}

	if(string(argv[i]) == "-p" ){
	    printPosToskip=true;
	    continue;
	}

	if(string(argv[i]) == "-s" ){
	    summary=true;
	    continue;
	}
	
	if(string(argv[i]) == "-v" ){
	    verbose=true;
	    continue;
	}

	if(string(argv[i])[0] != '-'  ){
	    //cout<<"end"<<i<<endl;
	    lastOpt=i;
	    break;
	}

	if(string(argv[i]) == "-err"  ){
	    errFile=string(argv[i+1]);
	    i++;
	    continue;
	}

	if(string(argv[i]) == "-deam5p"  ){
	    deam5pfreq=string(argv[i+1]);
	    i++;
	    continue;
	}

	if(string(argv[i]) == "-maxc"  ){
	    topCont=destringify<long double>(argv[i+1]);
	    i++;
	    continue;
	}

	if(string(argv[i]) == "-step"  ){
	    stepContEst=destringify<long double>(argv[i+1]);
	    i++;
	    continue;
	}


	if(string(argv[i]) == "-deam3p"  ){
	    deam3pfreq=string(argv[i+1]);
	    i++;
	    continue;
	}

	if(string(argv[i]) == "-double"  ){
	    doubleDeam=true;
	    continue;
	}

	if(string(argv[i]) == "-t"  ){
	    numberOfThreads =destringify<int>(argv[i+1]);
	    i++;
	    continue;
	}

	if( string(argv[i]) == "-o"  ){
	    outLog=string(argv[i+1]);
	    i++;
	    continue;
	}

	if(string(argv[i]) == "-nomq" ){
	    ignoreMQ=true;
	    continue;
	}

	cerr<<"Wrong option "<<string(argv[i])<<endl;
	return 1;

    }


    if(numberOfThreads<=0){
	cerr<<"Number of cores must be a positive integer "<<endl;
	return 1;
    }
    

    string consensusFile    = string(argv[lastOpt+0]); //consensus file
    string fastaFile        = string(argv[lastOpt+1]); //fasta file
    string bamfiletopen     = string(argv[lastOpt+2]); //bam file

    for(int i=(lastOpt+3);i<(argc);i++){ //all but the last 3 args
	//cout<<"cont "<<string(argv[i])<<endl;
	queueFilesToprocess.push(string(argv[i]));	
    }






    outLogFP.open(outLog.c_str());

    if (!outLogFP.is_open()){
	cerr << "Unable to write to qual file "<<outLog<<endl;
	return 1;
    }






    ////////////////////////////////////////////////////////////////////////
    //
    // BEGIN READING ERROR PROFILE
    //
    ////////////////////////////////////////////////////////////////////////

    readIlluminaError(errFile,illuminaErrorsProb);

    ////////////////////////////////////////////////////////////////////////
    //
    // END  READING ERROR PROFILE
    //
    ////////////////////////////////////////////////////////////////////////











    ////////////////////////////////////////////////////////////////////////
    //
    // BEGIN DEAMINATION PROFILE
    //
    ////////////////////////////////////////////////////////////////////////
    readNucSubstitionFreq(deam5pfreq,sub5p);
    readNucSubstitionFreq(deam3pfreq,sub3p);
   

    int defaultSubMatchIndex=0;

    for(int nuc1=0;nuc1<4;nuc1++){
    	for(int nuc2=0;nuc2<4;nuc2++){	    
    	    if(nuc1==nuc2)
		defaultSubMatch.s[ defaultSubMatchIndex++ ] = 1.0;
	    else
		defaultSubMatch.s[ defaultSubMatchIndex++ ] = 0.0;
    	}
    	
    }


    
    ////////////////////////////////////////////////////////////////////////
    //
    // END  DEAMINATION PROFILE
    //
    ////////////////////////////////////////////////////////////////////////








    ////////////////////////////////////////////////////////////////////////
    //
    // BEGIN ENDOGENOUS PROFILE
    //
    ////////////////////////////////////////////////////////////////////////

    readMTConsensus(consensusFile,pos2phredgeno,sizeGenome,posOfIndels);


    definedSite    = new vector<bool>(sizeGenome+1,false); // if there is data
    skipPositions  = new vector<bool>(sizeGenome+1,false); // if the site has such a low prior that it can be overlooked


#ifdef MINDISTFROMINDEL //delete records in pos2phredgeno found near indels

    // map<int, PHREDgeno>          pos2phredgeno;
    // vector<int>                  posOfIndels;    
    for(unsigned int idxindel=0;idxindel<posOfIndels.size();idxindel++){
	//cout<<posOfIndels[idxindel]<<endl;
	for(int indelPos=(int(posOfIndels[idxindel])-MINDISTFROMINDEL);
	    indelPos<=(int(posOfIndels[idxindel])+MINDISTFROMINDEL);
	    indelPos++){

	    map<int,PHREDgeno>::iterator itposgeno = pos2phredgeno.find(indelPos);

		if(itposgeno != pos2phredgeno.end()){ //present, erase it
		    //cout<<"Removing "<<indelPos<<endl;
		    pos2phredgeno.erase(itposgeno);
		}
	}
    }
    
#endif


    ////////////////////////////////////////////////////////////////////////
    //
    // END ENDOGENOUS PROFILE
    //
    ////////////////////////////////////////////////////////////////////////

    for(int i=0;i<=sizeGenome;i++){
	 positionInformation toadd;

	 toadd.cov            = 0;
	 toadd.mapqAvg        = 0.0;
	 
	 infoPPos.push_back(toadd);
     }





     Fasta fastaReference;
     if ( !fastaReference.Open(fastaFile , fastaFile+".fai") ){ 
	 cerr << "ERROR: failed to open fasta file " <<fastaFile<<" and index " << fastaFile<<".fai"<<endl;
	 return false;
     }


    cerr<<"Begin reading BAM file ..."<<endl;

    BamReader reader;

    if ( !reader.Open(bamfiletopen) ) {
	cerr << "Could not open input BAM files " <<bamfiletopen<< endl;
    	return 1;
    }



    //  // retrieve reference data
    const RefVector  references = reader.GetReferenceData();


    MyPileupVisitor* cv = new MyPileupVisitor(references,&fastaReference,&infoPPos,sizeGenome,ignoreMQ);
    PileupEngine pileup;
    pileup.AddVisitor(cv);


    BamAlignment al;
    unsigned int numReads=0;
    while ( reader.GetNextAlignment(al) ) {
	// cerr<<al.Name<<endl;
	numReads++;
	if(numReads !=0 && (numReads%100000)==0){
	    cerr<<"Processed "<<thousandSeparator(numReads)<<" reads"<<endl;
	}

	if(al.IsMapped() && !al.IsFailedQC()){
	    // cerr<<al.Name<<endl;
	    pileup.AddAlignment(al);
	}
	
    }
    //    cerr<<"done"<<endl;

    


    //clean up
    pileup.Flush();

    reader.Close();
    fastaReference.Close();



    cerr<<"... done reading BAM file"<<endl;
    










    findPosToSkip(printPosToskip,verbose);




    //init mutex
    // pthread_mutex_t  mutexQueue   = PTHREAD_MUTEX_INITIALIZER;
    // pthread_mutex_t  mutexCounter = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_init(&mutexQueue,   NULL);
    pthread_mutex_init(&mutexCounter, NULL);
    pthread_mutex_init(&mutexRank ,   NULL);

    pthread_t             thread[numberOfThreads];
    int                   rc=0;

    struct argsptcont * argptc = (struct argsptcont *)malloc(sizeof(struct argsptcont));
    argptc->doubleDeam = doubleDeam;
    
    // doneReading=true;    

    for(int i=0;i<numberOfThreads;i++){
	rc = pthread_create(&thread[i], NULL, mainContaminationThread, (void*)argptc);
	checkResults("pthread_create()\n", rc);
    }


    //threads are running here

    //waiting for threads to finish
    for (int i=0; i <numberOfThreads; ++i) {
	rc = pthread_join(thread[i], NULL);
	checkResults("pthread_join()\n", rc);
    }

    pthread_mutex_destroy(&mutexRank);
    pthread_mutex_destroy(&mutexQueue);
    pthread_mutex_destroy(&mutexCounter);

    vector< contRecord > vecPairfilenameLike;


    for(unsigned int i=0;i<outputToPrint.size();i++){
	contRecord toaddCt;
	toaddCt.fname="";
	// long double maxLike;
	// string fname = "";
	for(unsigned int j=0;j<outputToPrint[i].size();j++){

	    if(j==0){
		toaddCt.logLike      = outputToPrint[i][j].logLike;
		toaddCt.contEst      = outputToPrint[i][j].contaminationRate;
		toaddCt.fname        = outputToPrint[i][j].filename;
	    }

	    outLogFP<<outputToPrint[i][j].filename<<"\t"<<outputToPrint[i][j].contaminationRate<<"\t"<<outputToPrint[i][j].logLike<<endl;
	    
	    
	    //s+=pow(2.0,outputToPrint[i][j].logLike);
	    if(toaddCt.logLike < outputToPrint[i][j].logLike){
		toaddCt.logLike = outputToPrint[i][j].logLike;
		toaddCt.contEst = outputToPrint[i][j].contaminationRate;
	    }
	}
	
	vecPairfilenameLike.push_back( toaddCt );
    }


    
    if(summary){
	sort(vecPairfilenameLike.begin(),vecPairfilenameLike.end(),compContRecord);
	// outLogFP<<"#"<<vecPairfilenameLike.size()<<endl;
	cerr<<"################################################"<<endl<<"Sources of contamination (sorted bottom to top by log likelihood)"<<endl<<"################################################"<<endl;

	for(unsigned int i=0;i<vecPairfilenameLike.size();i++){
	    cerr<<vecPairfilenameLike[i].fname<<"\t"<<vecPairfilenameLike[i].contEst<<"\t"<<vecPairfilenameLike[i].logLike<<endl;
	}
	outLogFP<<"################################################"<<endl<<"Sources of contamination (sorted bottom to top by log likelihood)"<<endl<<"################################################"<<endl;

	for(unsigned int i=0;i<vecPairfilenameLike.size();i++){
	    outLogFP<<vecPairfilenameLike[i].fname<<"\t"<<vecPairfilenameLike[i].contEst<<"\t"<<vecPairfilenameLike[i].logLike<<endl;
	}

	outLogFP.close();
    }

    delete cv;

    pthread_exit(NULL);

    delete definedSite;
    delete skipPositions;


    return 0;
}

