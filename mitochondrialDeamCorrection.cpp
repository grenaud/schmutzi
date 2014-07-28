/*
 * mitonchondrialDeam
 * Date: Mar-26-2014 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */



// #include <boost/math/special_functions/beta.hpp>

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
#include <math.h>

#include "utils.h"
#include "miscfunc.h"
#include "ReconsReferenceBAM.h"

using namespace BamTools;
using namespace std;

#define MAXCOV 5000
const long double PI  = atanl(1.0L)*4;   


//! Computes log of probability density function
/*!
  Computes the log base 10 of pdf(x)
  
  \param mu The mean (location)
  \param sigma The variance (scale)
  \param x The value for which we want the pdf
  \return The values of log base 10 of (pdf(x))
*/
long double logcomppdf(long double mu,long double sigma,long double x){
    if(x==0){
	x=1.0;
    }

    double two = 2.0;   
    double exponent = log(x) - mu;
    exponent *= -exponent;
    exponent /= two * sigma * sigma;
    
    long double result = exponent/logl(10);
    result -= logl(sigma * sqrtl(two * PI) * x)/logl(10);

    return result;
}


// #define DEBUG1
// #define DEBUG2


// typedef struct { 
//     double likeBaseNoindel[4];
//     int  covPerBase[4];
//     double mapqAvg;
    
//     int numDel;
//     vector<string> insertionRight;
//     int cov;
// } positionInformation;


char offsetQual=33;
double likeMatch[64];
double likeMismatch[64];
double likeMatchMQ[64][64];
double likeMismatchMQ[64][64];

double likeMatchProb[64];
double likeMismatchProb[64];
double likeMatchProbMQ[64][64];
double likeMismatchProbMQ[64][64];

double probMapping[64];
double probMismapping[64];


double probLengthEndo[1000];

probSubstition illuminaErrorsProb;
vector<probSubstition> sub5p;
vector<probSubstition> sub3p;
probSubstition defaultSubMatch;

string dnaAlphabet="ACGT";
map<int, PHREDgeno> pos2phredgeno;

map<string, double > read2endoProb; //map seq id to probability that the read is endogenous using a deamination model
double read2endoProbInit=false;

#define MIN(a,b) (((a)<(b))?(a):(b))


template <typename T>
inline string arrayToStringInt(const T toPrint[] ,const int size,const string separator=","){
    if(size == 0){
    	return "";
    }
    string toReturn="";
    for(int i=0;i<(size-1);i++){
    	toReturn+=(stringify(int(toPrint[i]))+separator);
    }
    toReturn+=(stringify(int(toPrint[ size -1 ])));
    return toReturn;
}


inline void transformRef(char * refeBase,char * readBase){
    if( (*refeBase) == 'M'){
	(*refeBase)=(*readBase);
    }
    
}


inline bool hasIinfirstOrLastTwoBases(const string & reconstructedReference){
    if(reconstructedReference.length() <= 4){
	cerr<<"ERROR read has length less than 4 bp"<<endl;
	exit(1);
    }

    for(unsigned int j=0;j<2;j++){
	if(reconstructedReference[j] == 'I')
	    return true;
    }


    for(unsigned int j=(reconstructedReference.length()-2);
	j<(reconstructedReference.length());
	j++){
	if(reconstructedReference[j] == 'I')
	    return true;
    }

    return false;
}

inline bool deletionsNextToTwo(const BamAlignment  * al){
    vector<int> lengthOfNonDels;
    vector<CigarOp> cigarData=al->CigarData;
    bool foundDel=false;
    for(unsigned int i=0;i<cigarData.size();i++){
        if(cigarData[i].Type == 'D'){
	    foundDel=true;
	}else{
	    lengthOfNonDels.push_back(cigarData[i].Length);
	}
    }

    if(foundDel){
	if(lengthOfNonDels[0]<=2)
	    return true;
	if(lengthOfNonDels[ lengthOfNonDels.size() -1 ]<=2)
	    return true;
    }
    
    return false;
}

 

//checks for an 'R' or 'S' for soft clip
inline bool hasBadCharacter(const string & reconstructedReference){

    for(unsigned int j=0;j<(reconstructedReference.length());j++){
	if(reconstructedReference[j] == 'R'  || 
	   reconstructedReference[j] == 'S' ){
	    return true;
	}
    }
    return false;
}

inline bool skipAlign(const string & reconstructedReference,const BamAlignment  * al,unsigned int * skipped){
    if(hasBadCharacter(reconstructedReference)){
	(*skipped)++;
	return true;
    }
	

    if(hasIinfirstOrLastTwoBases(reconstructedReference)){
	(*skipped)++;
	return true;
    }

    if(deletionsNextToTwo(al)){
	(*skipped)++;
	return true;
    }

    return false;
}


// typedef struct { 
//     char ref;
//     char alt;
//      unsigned int pos;
//     double contFreq;
//     unsigned char cov;
//     unsigned char refCov;
//     unsigned char altCov;

//     //    unsigned int ;

//     bool      refOrAlt[MAXCOV]; //ref = false, alt=true
//     unsigned char mapq[MAXCOV];
//     unsigned char quals[MAXCOV];;
//     unsigned char dist5p[MAXCOV];;
//     unsigned char dist3p[MAXCOV];;   
//  } positionInfo;


void  printLogAndGenome(int sizeGenome,
			const vector<singlePosInfo> & infoPPos,
			const string outSeq,
			const string outLog, 
			const string genomeRef,
			const int minQual,
			const string nameMT){
    
    ofstream outSeqFP ;
    ofstream outLogFP;
    outSeqFP.open(outSeq.c_str());

    if (!outSeqFP.is_open()){
	cerr << "Unable to write to seq file "<<outSeq<<endl;
	exit(1);
    }

    outLogFP.open(outLog.c_str());

    if (!outLogFP.is_open()){
	cerr << "Unable to write to qual file "<<outLog<<endl;
	exit(1);
    }


    string genomeToPrint="";
    stringstream logToPrint;


    logToPrint<<"pos\trefBase\tbase\tqual\tavgmapq\tcov\tsupp\tpa\tpc\tpg\tpt\n";
    //genomeRef
    for(int i=0;i<sizeGenome;i++){

	 //if half of the reads support an deletion in reads (ins in reference) 
	 //move to next position since deletion
	 if( (double(infoPPos[i].numDel)/double(infoPPos[i].cov)) >= 0.5){
	     //min quality ? based on what ?

	     logToPrint<<(i+1)<<"\t"<<genomeRef[i]<<"\t"<<"D"<<"\t"<<-10.0*(log(1.0-(double(infoPPos[i].numDel)/double(infoPPos[i].cov)))/log(10.0))<<"\t"<<infoPPos[i].mapqAvg<<"\t"<<infoPPos[i].cov<<"\t"<<infoPPos[i].numDel<<"\t0.0\t0.0\t0.0\t0.0"<<endl;

	     PHREDgeno toadd;

	     toadd.ref       = genomeRef[i];
	     toadd.consensus = 'D';
	     
	     pos2phredgeno[   (i+1)   ] = toadd;

	     continue;

	     
	 }


	 
	 if(infoPPos[i].cov == 0){
	     genomeToPrint+="N";
	     logToPrint<<(i+1)<<"\t"<<genomeRef[i]<<"\tN\t0\t0\t0\t0\t0.0\t0.0\t0.0\t0.0"<<endl;
	     continue;
	 }
	     

	 //cout<<i<<"\t"<<infoPPos[i].cov<<"\t"<<vectorToString(infoPPos[i].insertionRight)<<"\t"<<infoPPos[i].numDel<<"\t";
	 double bestLike=-INFINITY;
	 int    bestNuc=-1;
	 for(unsigned int nuc=0;nuc<4;nuc++){
	     //cout<<infoPPos[i].likeBaseNoindel[nuc]<<"\t";
	     if(infoPPos[i].likeBaseNoindel[nuc] > bestLike){
		 bestLike=infoPPos[i].likeBaseNoindel[nuc];
		 bestNuc=nuc;
	     }
	 }

	 vector<int> bestNucs;

	 for(unsigned int nuc=0;nuc<4;nuc++){
	     if(infoPPos[i].likeBaseNoindel[nuc] == bestLike){
		 bestNucs.push_back(nuc);
	     }
	 }
	 
	 if(bestNucs.size() > 1){ // multiple equally likely nuc, use coverage to call best one
	     // cerr<<"size "<<bestNucs.size()<<endl;
	     // return 1;
	     int bestCov=-1;
	     int bestCovN=bestNucs[0];
	     
	     for(unsigned int bc=0;bc<bestNucs.size();bc++){
		 if(infoPPos[i].covPerBase[ bestNucs[bc] ] > bestCov){
		     bestCov  =infoPPos[i].covPerBase[ bestNucs[bc] ];
		     bestCovN =                        bestNucs[bc];		     
		 }
	     }
	     
	     bestNuc = bestCovN;
	 }
	 
	 double sumLogLikeAll          = 0.0;

	 double sumLogLikeOnlyBest     = 0.0;
	 double sumLogLikeAllButBest   = 0.0;
	 bool sumLogLikeAllB           = true;
	 bool sumLogLikeOnlyBestB      = true;
	 bool sumLogLikeAllButBestB    = true;
	 
	 // int nuc=0;
	 // sumLogLikeAll              =  infoPPos[i].likeBaseNoindel[nuc];  //oplus= log10( pow(10,x)+pow(10,y) )
	 // if(nuc==bestNuc)
	 //     sumLogLikeOnlyBest     =  infoPPos[i].likeBaseNoindel[nuc];
	 // else
	 //     sumLogLikeAllButBest   =  infoPPos[i].likeBaseNoindel[nuc];

	 double sumLogForNucs[4];
	 bool sumLogForNucsB[4];
	 for(int nuc=0;nuc<4;nuc++){
	     sumLogForNucsB[nuc]=true;
	 }

	 for(int nuc=0;nuc<4;nuc++){
	     //cout<<"nuc\t"<<nuc<<"\t"<<bestNuc<<"\t"<<infoPPos[i].likeBaseNoindel[nuc]<<"\t"<<pow(10.0,infoPPos[i].likeBaseNoindel[nuc])<<endl;
	     
	     for(int nuc2=0;nuc2<4;nuc2++){
		 if(nuc!=nuc2){
		     if(sumLogForNucsB[nuc]){
			 sumLogForNucs[nuc]         = infoPPos[i].likeBaseNoindel[nuc2]; 
			 sumLogForNucsB[nuc]        = false;
		     }else{
			 sumLogForNucs[nuc]         = oplus(sumLogForNucs[nuc], infoPPos[i].likeBaseNoindel[nuc2]);
		     }
		 }
	     }

	     if(sumLogLikeAllB){
		 sumLogLikeAll              =  infoPPos[i].likeBaseNoindel[nuc];  //oplus= log10( pow(10,x)+pow(10,y) )
		 sumLogLikeAllB=false;
	     }else{
		 sumLogLikeAll              =  oplus(sumLogLikeAll,        infoPPos[i].likeBaseNoindel[nuc]);  //oplus= log10( pow(10,x)+pow(10,y) )
	     }

	     if(nuc==bestNuc){

		 if(sumLogLikeOnlyBestB){
		     sumLogLikeOnlyBest     = infoPPos[i].likeBaseNoindel[nuc];
		     sumLogLikeOnlyBestB    = false;
		 }else{
		     sumLogLikeOnlyBest     = oplus(sumLogLikeOnlyBest,infoPPos[i].likeBaseNoindel[nuc]);
		 }

	     }else{
		 //sumLogLikeAllButBest   = oplus(sumLogLikeAllButBest, infoPPos[i].likeBaseNoindel[nuc]);
		 if(sumLogLikeAllButBestB){
		     sumLogLikeAllButBest   = infoPPos[i].likeBaseNoindel[nuc];
		     sumLogLikeAllButBestB  = false;		     
		 }else{
		     sumLogLikeAllButBest   = oplus(sumLogLikeAllButBest, infoPPos[i].likeBaseNoindel[nuc]);
		 }
	     }		 	    
	 }//end for nuc
	 

	 if( (-10.0*(sumLogLikeAllButBest-sumLogLikeAll)) >= minQual){
	     genomeToPrint+=dnaAlphabet[bestNuc];
	 }else{
	     genomeToPrint+="N";
	 }



	 logToPrint<<(i+1)<<"\t"<<genomeRef[i]<<"\t"<<dnaAlphabet[bestNuc]<<"\t"<<(-10*(sumLogLikeAllButBest-sumLogLikeAll)+0.0)<<"\t"<<infoPPos[i].mapqAvg<<"\t"<<infoPPos[i].cov<<"\t"<<infoPPos[i].covPerBase[bestNuc];//<<endl;
	 PHREDgeno toadd;

	 for(int nuc=0;nuc<4;nuc++){
	     logToPrint<<"\t"<<(-10*(sumLogForNucs[nuc]-sumLogLikeAll));
	     toadd.phred[nuc]  =     (-10*(sumLogForNucs[nuc]-sumLogLikeAll));
	     toadd.perror[nuc] = pow(10.0,(sumLogForNucs[nuc]-sumLogLikeAll));
	 }
	 logToPrint<<endl;
	 toadd.ref       = genomeRef[i];
	 toadd.consensus = dnaAlphabet[bestNuc];

	 pos2phredgeno[   (i+1)   ] = toadd;
	 // cout<<(i+1)<<"\t"<<toadd.ref<<"\t"<<toadd.consensus<<endl;


	 if(infoPPos[i].insertionRight.size() != 0){

	     map<string,int> insert2count;
	     for(unsigned int k=0;k<infoPPos[i].insertionRight.size();k++){		 
		 insert2count[ infoPPos[i].insertionRight[k] ]++;
	     }
	     
	     
	     int     mostCommonInsCount=-1;
	     string  mostCommonIns="";


	     for (map<string,int>::iterator itdel=insert2count.begin(); 
		  itdel!=insert2count.end(); ++itdel){
		 //cout<<itdel->first<<"\t"<<itdel->second<<endl;
		 if( itdel->second > mostCommonInsCount){
		     mostCommonInsCount = itdel->second;
		     mostCommonIns      = itdel->first;
		 }		     
	     }


	     //if half of the reads support an insertions in reads (deletions in reference) 
	     if( (double(mostCommonInsCount)/double(infoPPos[i].cov)) >= 0.5){
		 //outSeqFP<<mostCommonIns;
		 //return 1;
		 //if( (-10.0*(log(1.0-double(mostCommonInsCount)/double(infoPPos[i].cov))/log(10.0)))>=minQual){
		 //min quality ? based on what ?
		 genomeToPrint+=mostCommonIns;
		 //}else{
		 //do not add insert
		 //}

		 //logToPrint<<string('!',mostCommonIns.size());
		 for(unsigned int k=0;k<(mostCommonIns.size());k++){
		     logToPrint<<(i+1)<<"i\t"<<"-"<<"\t"<<mostCommonIns[k]<<"\t"<<-10.0*(log(1.0-double(mostCommonInsCount)/double(infoPPos[i].cov))/log(10.0))<<"\t"<<infoPPos[i].mapqAvg<<"\t"<<infoPPos[i].cov<<"\t"<<mostCommonInsCount<<"\t0.0\t0.0\t0.0\t0.0"<<endl;

		 }
		 
	     }

	 }
     }


     string genomeToPrintCopy="";
     for(unsigned int i=1;i<(genomeToPrint.size()+1);i++){
	 genomeToPrintCopy+=genomeToPrint[i-1];
	 if(i!=0 && (i%80 == 0))
	     genomeToPrintCopy+="\n";
     }
     outSeqFP<<nameMT+"\n"+genomeToPrintCopy+"\n";
     

     outSeqFP.close() ;
     outLogFP<<logToPrint.str()<<endl;
     outLogFP.close();
}



//unsigned int posFound;
//vector<positionInfo> * positionsInfoFound;
// map<string, map<unsigned int,contaminationInfo> > contFreq;

class MyPileupVisitor : public PileupVisitor {
  
    public:
    MyPileupVisitor(const RefVector& references, Fasta * fastaReference,vector<singlePosInfo> * infoPPos,int sizeGenome,bool ignoreMQ)
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
	    //cout<<"visit"<<endl;
	    char referenceBase = 'N';
	    
	    unsigned int posAlign = pileupData.Position+1;
	    int posVector=int(pileupData.Position)%sizeGenome;
	    

	    // if( (posAlign%100) == 0){
	    // 	cerr<<"pos  = "<<posAlign<<endl;
	    // }
	    //cout<<endl<<"pos = "<<posAlign<<"\t"<<posVector;
	    

	    //for some reason, we have to do -1 on the .Position
	    if ( !m_fastaReference->GetBase(pileupData.RefId, posAlign-1, referenceBase ) ) {
		cerr << "bamtools convert ERROR: pileup conversion - could not read reference base from FASTA file at chr "<<pileupData.RefId<<" position "<<(posAlign-1) << endl;
		exit(1);
	    }


	    //cout<<"\t"<<referenceBase<<"\t"<<endl;


	    for(unsigned int i=0;i<pileupData.PileupAlignments.size();i++){
		m_infoPPos->at(posVector).cov++;

		//Add mapq
		m_infoPPos->at(posVector).mapqAvg += pow(10.0, double(pileupData.PileupAlignments[i].Alignment.MapQuality)/-10.0);

		// cout<<pileupData.PileupAlignments[i].Alignment.Name<<"\t"<<
		//     pileupData.PileupAlignments[i].IsCurrentDeletion<<"\t"<<
		//     pileupData.PileupAlignments[i].IsNextDeletion<<"\t"<<
		//     pileupData.PileupAlignments[i].IsNextInsertion<<"\t"<<
		//     pileupData.PileupAlignments[i].DeletionLength<<"\t"<<
		//     pileupData.PileupAlignments[i].InsertionLength<<"\t"<<
		//     pileupData.PileupAlignments[i].IsSegmentBegin<<"\t"<<
		//     pileupData.PileupAlignments[i].IsSegmentEnd<<"\t"<<
		//     endl;
	    }

	    m_infoPPos->at(posVector).mapqAvg = m_infoPPos->at(posVector).mapqAvg/double(m_infoPPos->at(posVector).cov);
	    m_infoPPos->at(posVector).mapqAvg = -10.0*( log( m_infoPPos->at(posVector).mapqAvg )/log(10.0) );

	    //insertion in the reads/deletion in the reference
	    for(unsigned int i=0;i<pileupData.PileupAlignments.size();i++){		
		if( !pileupData.PileupAlignments[i].IsCurrentDeletion &&
		    pileupData.PileupAlignments[i].IsNextInsertion &&
		    (pileupData.PileupAlignments[i].InsertionLength>0)){
		    string insert="";
		    for(int del=1;del<=pileupData.PileupAlignments[i].InsertionLength;del++){
			insert+=  pileupData.PileupAlignments[i].Alignment.QueryBases[ pileupData.PileupAlignments[i].PositionInAlignment+del ];
		    }
		    //cout<<"ins "<<insert<<endl;
		    m_infoPPos->at(posVector).insertionRight.push_back(insert);
		}
	    }


	    //deletion in the reads/insertion in the reference
	    for(unsigned int i=0;i<pileupData.PileupAlignments.size();i++){		
		if( pileupData.PileupAlignments[i].IsCurrentDeletion &&
		    pileupData.PileupAlignments[i].IsNextInsertion &&
		    (pileupData.PileupAlignments[i].InsertionLength == 0)){
		    //continue;
		    //cout<<"del"<<endl;
		    m_infoPPos->at(posVector).numDel++;
		}
	    }
	    


	    for(unsigned int nuc=0;nuc<4;nuc++){
		unsigned int numReads=0;
		char currentNuc=dnaAlphabet[nuc];
		for(unsigned int i=0;i<pileupData.PileupAlignments.size();i++){		
		    
		    //skip deletion in the reads/insertion in the reference
		    if( pileupData.PileupAlignments[i].IsCurrentDeletion &&
			pileupData.PileupAlignments[i].IsNextInsertion ){
			continue;
		    }
		   		    
		    char b   = pileupData.PileupAlignments[i].Alignment.QueryBases[pileupData.PileupAlignments[i].PositionInAlignment];
		    char q   = pileupData.PileupAlignments[i].Alignment.Qualities[pileupData.PileupAlignments[i].PositionInAlignment]-offsetQual;
		    int  m   = int(pileupData.PileupAlignments[i].Alignment.MapQuality);
		    
		    //skip unresolved
		    if(b == 'N')
			continue;

		    // BEGIN DEAMINATION COMPUTATION
		    //zero base distance to the 5p/3p end
		    int dist5p=-1;
		    int dist3p=-1;

		    if( pileupData.PileupAlignments[i].Alignment.IsReverseStrand() ){
			dist5p = pileupData.PileupAlignments[i].Alignment.QueryBases.size() - pileupData.PileupAlignments[i].PositionInAlignment-1;
			dist3p = pileupData.PileupAlignments[i].PositionInAlignment;
		    }else{
			dist5p = pileupData.PileupAlignments[i].PositionInAlignment;
			dist3p = pileupData.PileupAlignments[i].Alignment.QueryBases.size() - pileupData.PileupAlignments[i].PositionInAlignment-1;
		    }
		    		    
		    probSubstition * probSubMatchToUse = &defaultSubMatch ;
		    
		    // for(int jas=0;jas<16;jas++){
		    // 	cout<<jas<<"\t"<<probSubMatchToUse->s[jas]<<endl;
		    // }
		    // cout<<"sub\t"<<dist5p<<"\t"<<dist3p<<"\t"<<int(sub5p.size())<< "\t"<<int(sub3p.size()) <<endl;

		    if(dist5p <= (int(sub5p.size()) -1)){
			probSubMatchToUse = &sub5p[ dist5p ];			
			 // cout<<"5p sub"<<endl;
		    }
		    
		    if(dist3p <= (int(sub3p.size()) -1)){
			probSubMatchToUse = &sub3p[ dist3p ];
			 // cout<<"3p sub"<<endl;
		    }
		    
		    //we have substitution probabilities for both... take the closest
		    if(dist5p <= (sub5p.size() -1) &&
		       dist3p <= (sub3p.size() -1) ){
			
			if(dist5p < dist3p){
			    probSubMatchToUse = &sub5p[ dist5p ];
			    // cout<<"5p sub"<<endl;
			}else{
			    probSubMatchToUse = &sub3p[ dist3p ];
			    // cout<<"3p sub"<<endl;
			}
			
		    }
		    // END DEAMINATION COMPUTATION



		    if(currentNuc == b){//match
			m_infoPPos->at(posVector).covPerBase[nuc]++;			
		    }else{//mismatch
		    }


		    // for(int jas=0;jas<16;jas++){
		    // 	cout<<jas<<"\t"<<probSubMatchToUse->s[jas]<<endl;
		    // }
		    //exit(1);
		    // b   is the observed
		    // nuc is the model
		    int dinucIndex;//The index depends on the strand
		    if( pileupData.PileupAlignments[i].Alignment.IsReverseStrand() ){
			dinucIndex= (3-nuc)*4+baseResolved2int(complement(b));
		    }else{
			dinucIndex=    nuc *4+baseResolved2int(b);
		    }
		    // = nuc*4+baseResolved2int(b);
		    
		    //                        (1-e)           *  p(sub|1-e)                          + (e) *  p(sub|1-e)
		    double probBase =  likeMatchProb[int(q)]  * (probSubMatchToUse->s[dinucIndex] )  + (1.0 - likeMatchProb[int(q)])*(illuminaErrorsProb.s[dinucIndex]);
		    //m_infoPPos->at(posVector).likeBaseNoindel[nuc] += 
		    double probFinal;
		    
		    //     read2endoProb[ al.Name+"#"+ stringify(al.AlignmentFlag) ] = probEndo;
		    
		    // }
		    // cerr<<"...  done"<<endl;
		    // read2endoProbInit=true;
		    if(read2endoProbInit){

			map<string,double>::iterator itRead2endoProb = read2endoProb.find( pileupData.PileupAlignments[i].Alignment.Name+
											   "#"+ 
											   stringify(pileupData.PileupAlignments[i].Alignment.AlignmentFlag) );

			double probEndogenous=0.5; 
			if( itRead2endoProb == read2endoProb.end() ){ //skipped due to deletions near the end or something			    
			    // cerr<<"Error: cannot find the ID "<<( pileupData.PileupAlignments[i].Alignment.Name+
			    // 					  "#"+ 
			    // 					  stringify(pileupData.PileupAlignments[i].Alignment.AlignmentFlag) )<<endl;
			    // exit(1);

			    probEndogenous = 0.5;                         //assume that a read is equally to belong to either the endo or the contaminant
			}else{
			    probEndogenous = itRead2endoProb->second;
			}

			//redifining probBase
			//            p(endo)                 *
			//cout<<probBase<<"\t";
			probBase = ( probEndogenous )*(probBase)  + ( 1.0-probEndogenous )*(0.25) ;
			//cout<< itRead2endoProb->second<<"\t"<<probBase<<endl<<endl;
		    }



		    if(ignoreMQ){ //ignore MQ
			probFinal = (               probBase                          );
		    }else{
			probFinal = (probMapping[m]*probBase + probMismapping[m]*0.25);
		    }
		    
		    m_infoPPos->at(posVector).likeBaseNoindel[nuc] += log(probFinal)/log(10);;

#ifdef DEBUG2		    
		    cout<<"b_obs="<<b<<" n_model="<<currentNuc<<"\tindex="<<dinucIndex<<" q="<<int(q)<<" m="<<m<<endl;
		    cout<<"p(match)="<<likeMatchProb[int(q)]  <<"\t"
			<<"p(sub|match)="<< (probSubMatchToUse->s[dinucIndex] ) <<"\t"
			<<"p(match)*p(sub|match)"<< (likeMatchProb[int(q)]  * (probSubMatchToUse->s[dinucIndex] )) <<"\t"
			<<"p(mismatch)="<<(1.0 - likeMatchProb[int(q)]) <<"\t"
			<<"p(sub|mismatch)="<<(illuminaErrorsProb.s[dinucIndex])<<"\t"
			<<"p(mismatch)*p(sub|mismatch)"<<( (1.0 - likeMatchProb[int(q)]) *(illuminaErrorsProb.s[dinucIndex]) )
			<<endl;
		    cout<<"final="<<(probFinal)<<"\t"<<log(probFinal)/log(10)<<endl;
#endif


		    
		    // numReads++;
		    if(numReads >= MAXCOV){
			break;
		    }

		}//end each read
#ifdef DEBUG2	
		cout<<posAlign<<"\tllik\t"<<(m_infoPPos->at(posVector).likeBaseNoindel[nuc])<<endl;
		cout<<"-----------------"<<endl;
#endif

	    }//end for each nuc
	    
	    //exit(1);
	    //store
	    // toStore.cov = numReads;

	    // toStore.alt = alt;
	    //  toStore.pos = posAlign;

	    // positionsInfoFound->push_back( toStore );

	    // if(positionsInfoFound->size()%10000 == 0){
	    // 	cerr<<"Found  "<<thousandSeparator(positionsInfoFound->size())<<" positions "<<endl;
	    // }
	    
	    
	    // skiptonextpos:
	    //     return;

	    // cout <<m_references[pileupData.RefId].RefName << "\t" 
	    // 	 <<referenceBase<<"\t"
	    // 	 << pileupData.Position << "\t" 
	    // 	 	 << pileupData.PileupAlignments.size() << endl;
        }
        
    private:
    RefVector m_references;
    Fasta * m_fastaReference;
    vector<singlePosInfo> * m_infoPPos;
    int sizeGenome;
    bool ignoreMQ;
    //        ostream*  m_out;
};

























void iterateOverReads(const string fastaFile,
		      const string bamfiletopen,
		      vector<singlePosInfo>  * infoPPos,
		      const int sizeGenome,
		      const bool ignoreMQ ){
    infoPPos->clear(); //clear previous data

    cerr<<"Reading genome file ..."<<endl;
     for(int i=0;i<sizeGenome;i++){
	 singlePosInfo toadd;
	 toadd.numDel         = 0;

	 toadd.cov            = 0;
	 toadd.mapqAvg        = 0.0;

	 for(unsigned int nuc=0;nuc<4;nuc++){
	     toadd.likeBaseNoindel[nuc] = 0;
	     toadd.covPerBase[nuc]      = 0;
	 }

	 infoPPos->push_back(toadd);
     }
    cerr<<"... done"<<endl;

    // cout<<fastaFile<<"\t"<<bamfiletopen<<"\t"<<infoPPos->size()<<"\t"<<sizeGenome<<"\t"<<ignoreMQ<<endl;


    Fasta fastaReference;
    if ( !fastaReference.Open(fastaFile , fastaFile+".fai") ){ 
	cerr << "ERROR: failed to open fasta file " <<fastaFile<<" and index " << fastaFile<<".fai"<<endl;
	exit(1);
    }
	
	
    BamReader reader;

    if ( !reader.Open(bamfiletopen) ) {
	cerr << "Could not open input BAM files " <<bamfiletopen<< endl;
	exit(1);
    }

    if ( !reader.OpenIndex(bamfiletopen+".bai") ) {
	cerr << "Could not open input index for BAM files" <<bamfiletopen+".bai"<< endl;
	exit(1);
    }
	
    //  // retrieve reference data
    const RefVector  references = reader.GetReferenceData();



    MyPileupVisitor* cv = new MyPileupVisitor(references,&fastaReference,infoPPos,sizeGenome,ignoreMQ);
    PileupEngine pileup;
    pileup.AddVisitor(cv);


    BamAlignment al;
    unsigned int numReads=0;
    cerr<<"Reading BAM file ..."<<endl;

    while ( reader.GetNextAlignment(al) ) {
	//cerr<<al.Name<<endl;
	numReads++;
	if(numReads !=0 && (numReads%100000)==0){
	    cerr<<"Read "<<thousandSeparator(numReads)<<" reads"<<endl;
	}

	if(al.IsMapped() && 
	   !al.IsFailedQC()){
	    // cerr<<al.Name<<endl;
	    pileup.AddAlignment(al);
	}
	    
    }
    cerr<<"...  done"<<endl;
    
    

    

    //clean up
    pileup.Flush();

    reader.Close();
    fastaReference.Close();
    delete cv;



}


void initScores(){

    for(int i=0;i<2;i++){
        likeMatch[i]        = log1p(    -pow(10.0,2.0/-10.0) )    /log(10);         
        likeMismatch[i]     = log  (     pow(10.0,2.0/-10.0)/3.0 )/log(10);

	likeMatchProb[i]           = 1.0-pow(10.0,2.0/-10.0) ;
        likeMismatchProb[i]        =     pow(10.0,2.0/-10.0)/3.0 ;
    }


    //Computing for quality scores 2 and up
    for(int i=2;i<64;i++){
        likeMatch[i]        = log1p(    -pow(10.0,i/-10.0) )     /log(10);          
        likeMismatch[i]     = log  (     pow(10.0,i/-10.0)/3.0  )/log(10);

        likeMatchProb[i]           = 1.0-pow(10.0,i/-10.0);
        likeMismatchProb[i]        =     pow(10.0,i/-10.0)/3.0;
    }


    //Adding mismapping probability
    for(int m=0;m<64;m++){

	double incorrectMappingProb   =     pow(10.0,m/-10.0); //m
	double correctMappingProb     = 1.0-pow(10.0,m/-10.0); //1-m
	
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
    	for(int i=2;i<64;i++){
	    //  (1-m)(1-e) + m/4  = 1-m-e+me +m/4  = 1+3m/4-e+me
    	    likeMatchMQ[m][i]         = log(  correctMappingProb*(1.0-pow(10.0,i/-10.0)    ) + incorrectMappingProb/4.0    )/log(10);    
	    //  (1-m)(e/3) + m/4  = e/3 -me/3 + m/4
    	    likeMismatchMQ[m][i]      = log(  correctMappingProb*(    pow(10.0,i/-10.0)/3.0) + incorrectMappingProb/4.0    )/log(10);    
	    
    	    likeMatchProbMQ[m][i]           = correctMappingProb*(1.0-pow(10.0,i/-10.0)    ) + incorrectMappingProb/4.0;
    	    likeMismatchProbMQ[m][i]        = correctMappingProb*(    pow(10.0,i/-10.0)/3.0) + incorrectMappingProb/4.0;
    	}


#ifdef DEBUG1
    	for(int i=0;i<64;i++){
	    cerr<<"m\t"<<m<<"\t"<<i<<"\t"<<likeMatchMQ[m][i]<<"\t"<<likeMismatchMQ[m][i]<<"\t"<<likeMatchProbMQ[m][i]<<"\t"<<likeMismatchProbMQ[m][i]<<endl;
	}
#endif

    }

}

int main (int argc, char *argv[]) {

    int sizeGenome=0;
    string outSeq  = "/dev/stdout";
    string outLog  = "/dev/stderr";
    string nameMT  = "MT";

  
    int minQual=0;
    bool ignoreMQ=false;


    ////////////////////////////////////
    // BEGIN Initializing scores      //
    ////////////////////////////////////
    initScores();
    ////////////////////////////////////
    //    END Initializing scores     //
    ////////////////////////////////////


    double locatione=0.0;
    double locationc=0.0;
    double scalee   =0.0;
    double scalec   =0.0;



    ////////////////////////////////////
    // BEGIN Parsing arguments        //
    ////////////////////////////////////

    //    return 1;
    string errFile    = getCWD(argv[0])+"illuminaProf/error.prof";
    string deam5pfreq = getCWD(argv[0])+"deaminationProfile/none.prof";
    string deam3pfreq = getCWD(argv[0])+"deaminationProfile/none.prof";
    // substitutionRates freqIlluminaError;
    vector<substitutionRates>    deam5Psub;
    vector<substitutionRates>    deam3Psub;
    bool deamread=false;
    double contaminationPrior=0.1;

    const string usage=string("\nThis program takes an aligned BAM file for a mitonchondria and calls a\nconsensus for the endogenous material\n\n\t"+
			      string(argv[0])+			      
			      " [options]  [reference fasta] [bam file] "+"\n\n"+

			      "\n\tOutput options:\n"+	
			      "\t\t"+"-seq  [fasta file]" +"\t\t"+"Output fasta file (default: stdout)"+"\n"+
			      "\t\t"+"-log  [log file]" +"\t\t"+"Output log  (default: stderr)"+"\n"+
			      "\t\t"+"-name [name]" +"\t\t\t"  +"Name  (default "+nameMT+") "+"\n"+
			      "\t\t"+"-qual [minimum quality]" +"\t\t"  +"Filter bases with quality less than this  (default "+stringify(minQual)+") "+"\n"+
			      // "\t\t"+"-cont" +"\t\t"+"Contamination allele frequency"+"\n"+

			      "\n\tDeamination options:\n"+				      
			      "\t\t"+"-deam5p [.prof file]" +"\t\t"+"5p deamination frequency (default: "+deam5pfreq+")"+"\n"+
			      "\t\t"+"-deam3p [.prof file]" +"\t\t"+"3p deamination frequency (default: "+deam3pfreq+")"+"\n"+
			      "\t\t"+"-deamread" +"\t\t\t"+"Set a prior on reads according to their deamination pattern (default: "+ booleanAsString(deamread) +")"+"\n"+
			      "\t\t"+"-cont [cont prior]"+"\t\t"+"If the -deamread option is specified, this is the contamination prior (default: "+ stringify(contaminationPrior) +")"+"\n"+
			      "\n\tLength options:\n"+				      
			      "\t\t"+"--loce"+  "\t\t\t\t"+"Location for lognormal dist for the endogenous sequences (default none)"+"\n"+
			      "\t\t"+"--scalee"+"\t\t\t"+"Scale for lognormal dist for the endogenous sequences (default none)"+"\n"+
			      "\t\t"+"--locc"+  "\t\t\t\t"+"Location for lognormal dist for the contaminant sequences (default none)"+"\n"+
			      "\t\t"+"--scalec"+"\t\t\t"+"Scale for lognormal dist for the contaminant sequences (default none)"+"\n"+

			      "\n\tComputation options:\n"+	
			      "\t\t"+"-nomq" +"\t\t\t\t"+"Ignore mapping quality (default: "+booleanAsString(ignoreMQ)+")"+"\n"+
			      "\t\t"+"-err" +"\t\t\t\t"+"Illumina error profile (default: "+errFile+")"+"\n"+
			      
			      "\n\tReference options:\n"+	
			      "\t\t"+"-l [length]" +"\t\t\t"+"Actual length of the genome used for"+"\n"+
			      "\t\t"+"  " +"\t\t\t\t"+"the reference as been wrapped around"+"\n"+
			      "\t\t"+"  " +"\t\t\t\t"+"by default, the length of the genome will be used "+"\n"+
			      
			      "");
			      
    if( (argc== 1) ||
	(argc== 2 && string(argv[1]) == "-h") ||
	(argc== 2 && string(argv[1]) == "-help") ||
	(argc== 2 && string(argv[1]) == "--help") ){
	cout<<usage<<endl;
	return 1;
    }

    string bamfiletopen = string(argv[argc-1]);//bam file
    string fastaFile    = string(argv[argc-2]);//fasta file
    bool specifiedContPrior=false;
    bool specifiedLoce   = false;
    bool specifiedLocc   = false;
    bool specifiedScalee = false;
    bool specifiedScalec = false;
    bool useLengthPrior  = false;

    for(int i=1;i<(argc-2);i++){ //all but the last 3 args

	
	if(strcmp(argv[i],"--loce") == 0 ){
	    locatione =destringify<double>(argv[i+1]);
	    i++;
	    specifiedLoce=true;
	    continue;
	}

	if(strcmp(argv[i],"--scalee") == 0 ){
	    scalee =destringify<double>(argv[i+1]);
	    i++;
	    specifiedScalee=true;
	    continue;
	}

	if(strcmp(argv[i],"--locc") == 0 ){
	    locationc =destringify<double>(argv[i+1]);
	    i++;
	    specifiedLocc=true;
	    continue;
	}

	if(strcmp(argv[i],"--scalec") == 0 ){
	    scalec =destringify<double>(argv[i+1]);
	    i++;
	    specifiedScalec=true;
	    continue;
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

	if(string(argv[i]) == "-deam3p"  ){
	    deam3pfreq=string(argv[i+1]);
	    i++;
	    continue;
	}

	if(string(argv[i]) == "-deamread"  ){
	    deamread=true;
	    continue;
	}


	if(strcmp(argv[i],"-cont") == 0 ){
	    contaminationPrior=destringify<double>(argv[i+1]);
	    specifiedContPrior=true;
	    i++;
	    continue;
	}

	// if(strcmp(argv[i],"-deam5") == 0 ){
	//     deam5File=string(argv[i+1]);
	//     i++;
	//     continue;
	// }

	if(strcmp(argv[i],"-nomq") == 0 ){
	    ignoreMQ=true;
	    continue;
	}

	if(strcmp(argv[i],"-seq") == 0 ){
	    outSeq=string(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-log") == 0 ){
	    outLog=string(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-qual") == 0 ){
	    minQual=destringify<int>(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-name") == 0 ){
	    nameMT=string(argv[i+1]);
	    i++;
	    continue;
	}



	if(strcmp(argv[i],"-l") == 0 ){
	    sizeGenome=atoi(argv[i+1]);
	    i++;
	    continue;
	}


	cerr<<"Error: unknown option "<<string(argv[i])<<endl;
	return 1;
    }


    if(specifiedLoce     ||
       specifiedLocc     || 
       specifiedScalee   ||
       specifiedScalec ){

	if( !(specifiedLoce     &&
	      specifiedLocc     && 
	      specifiedScalee   &&
	      specifiedScalec ) ){
	    cerr<<"Error: need to specify the location and scale for both the endogenous and contaminant"<<endl;
	    return 1;
	}

	for(int lengthMolecule=0;lengthMolecule<1000;lengthMolecule++){
	    double pdfEndo = pow(10.0,double(logcomppdf(locatione,scalee ,double(lengthMolecule))));
	    double pdfCont = pow(10.0,double(logcomppdf(locationc,scalec ,double(lengthMolecule))));

	    //double pdfCont = logcomppdf  
	    //cout<<lengthMolecule<<"\t"<<pdfEndo<<"\t"<<pdfCont<<"\t"<<(pdfEndo/(pdfEndo+pdfCont))<<endl;
	    //dobuel 
	    probLengthEndo[lengthMolecule] = (pdfEndo/(pdfEndo+pdfCont));
	    
	}
	useLengthPrior=true;

    }

    if(specifiedContPrior ){
	if( !(deamread || useLengthPrior) ){
	    cerr<<"Error: cannot specify -cont if you do not specify -deamread or use the length distribution priors, exiting"<<endl;
	    return 1;	
	}
    }


    ////////////////////////////////////
    // END Parsing arguments        //
    ////////////////////////////////////











    ////////////////////////////////////////////////////////////////////////
    //
    // BEGIN READING ERROR PROFILE
    //
    ////////////////////////////////////////////////////////////////////////

    readIlluminaError(errFile,illuminaErrorsProb);

    // for(int nuc1=0;nuc1<4;nuc1++){
    // 	for(int nuc2=0;nuc2<4;nuc2++){
	    
    // 	    cout<<illuminaErrorsProb.s[nuc2+nuc1*4]<<"\t";
    // 	}
    // 	cout<<endl;
    // }
    // return 1;
    
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

    // cout<<sub5p.size()<<endl;
    // cout<<sub3p.size()<<endl;

    // return 1;
    

    int defaultSubMatchIndex=0;

    for(int nuc1=0;nuc1<4;nuc1++){
    	for(int nuc2=0;nuc2<4;nuc2++){	    
    	    if(nuc1==nuc2)
		defaultSubMatch.s[ defaultSubMatchIndex++ ] = 1.0;
	    else
		defaultSubMatch.s[ defaultSubMatchIndex++ ] = 0.0;
    	}
    	
    }


    // for(int nuc1=0;nuc1<4;nuc1++){
    // 	for(int nuc2=0;nuc2<4;nuc2++){	    
    // 	    cout<<sub5p[0].s[nuc2+nuc1*4]<<"\t";
    // 	}
    // 	cout<<endl;
    // }

    // for(int nuc1=0;nuc1<4;nuc1++){
    // 	for(int nuc2=0;nuc2<4;nuc2++){	    
    // 	    cout<<sub3p[0].s[nuc2+nuc1*4]<<"\t";
    // 	}
    // 	cout<<endl;
    // }

    // return 1;
    
    ////////////////////////////////////////////////////////////////////////
    //
    // END  DEAMINATION PROFILE
    //
    ////////////////////////////////////////////////////////////////////////










    nameMT = ">"+nameMT;














     
    vector<singlePosInfo> infoPPos;

    if(isDos(fastaFile) || isMac(fastaFile) ){
	cerr << "File  "<<fastaFile<<" must be unix formatted, exiting"<<endl;
	return 1;
    }


    string line;

    bool firstFastaLine=true;
    string genomeRef="";

    igzstream fastaRefFile;
    fastaRefFile.open(fastaFile.c_str(),ios::in);
    if (fastaRefFile.good()){

	while ( getline (fastaRefFile,line)){
	
	    if(firstFastaLine){
		if(line[0]!='>'){
		    cerr << "File  "<<fastaFile<<" does not appear to be in fasta format"<<endl;
		    return 1;
		}
		firstFastaLine=false;
		continue;		
	    }

	    genomeRef+=line;
	}
	fastaRefFile.close();

    }else{
	cerr << "Cannot open fasta file  "<<fastaFile<<""<<endl;
	return 1;
    }

    // cerr<<genomeRef<<endl;

    // return 1;
    if(sizeGenome == 0){
	sizeGenome=genomeRef.size();//use the one from the fasta
    }

     // double ** likelihoodPBasePPos = new double * [genome.size()];
     // for(unsigned int i = 0; i< genome.size(); i++) {
     // 	likelihoodPBasePPos[i] = new double  [4];
     // }
     // //return 1;


// void iterateOverReads(const string fastaFile,
// 		      const string bamfiletopen,
// 		      vector<singlePosInfo>  infoPPos,
// 		      const int sizeGenome,
// 		      const bool ignoreMQ ){

    iterateOverReads(fastaFile,
		     bamfiletopen,
		     &infoPPos,
		     sizeGenome,
		     ignoreMQ );


     // Fasta fastaReference;
     // if ( !fastaReference.Open(fastaFile , fastaFile+".fai") ){ 
     // 	 cerr << "ERROR: failed to open fasta file " <<fastaFile<<" and index " << fastaFile<<".fai"<<endl;
     // 	 return false;
     // }



    // BamReader reader;

    // if ( !reader.Open(bamfiletopen) ) {
    // 	cerr << "Could not open input BAM files " <<bamfiletopen<< endl;
    // 	return 1;
    // }



    //  if ( !reader.OpenIndex(bamfiletopen+".bai") ) {
    // 	 cerr << "Could not open input index for BAM files" <<bamfiletopen+".bai"<< endl;
    // 	 return 1;
    //  }
     
    // //  // retrieve reference data
    // const RefVector  references = reader.GetReferenceData();


     // char referenceBase = 'N';
     // if(fastaReference.GetBase(20, 9411241, referenceBase ) ) {
	 
     // }
     // cout<<referenceBase<<endl;
     // if(fastaReference.GetBase(20, 9411242, referenceBase ) ) {
	 
     // }
     // cout<<referenceBase<<endl;
     // if(fastaReference.GetBase(20, 9411243, referenceBase ) ) {
	 
     // }
     // cout<<referenceBase<<endl;
     // if(fastaReference.GetBase(20, 9411244, referenceBase ) ) {
	 
     // }
     // cout<<referenceBase<<endl;

     // return 1;



     // vector<unsigned int> * coverageCounter=  new vector<unsigned int>(250+1,0);

     // positionsInfoFound = new vector<positionInfo>();
     // try { 
     // 	 positionsInfoFound->reserve(10000000);
     // } catch (bad_alloc const&) {
     // 	 cerr<<"Cannot allocate sufficient memory, try another system"<<endl;
     // 	 exit(1);
     // }
     //     posFound =  0;

     // cerr<<"done"<<endl;
     // sleep(100);
     // return 1;

    // MyPileupVisitor* cv = new MyPileupVisitor(references,&fastaReference,&infoPPos,sizeGenome,ignoreMQ);
    // PileupEngine pileup;
    // pileup.AddVisitor(cv);


    // BamAlignment al;
    // unsigned int numReads=0;
    // cerr<<"Reading BAM file ..."<<endl;

    // while ( reader.GetNextAlignment(al) ) {
    // 	//cerr<<al.Name<<endl;
    // 	numReads++;
    // 	if(numReads !=0 && (numReads%100000)==0){
    // 	    cerr<<"Read "<<thousandSeparator(numReads)<<" reads"<<endl;
    // 	}

    // 	if(al.IsMapped() && !al.IsFailedQC()){
    // 	    // cerr<<al.Name<<endl;
    // 	    pileup.AddAlignment(al);
    // 	}
	
    // }
    // cerr<<"...  done"<<endl;
    
    

    

    // //clean up
    // pileup.Flush();

    // reader.Close();
    // fastaReference.Close();
    // delete cv;


    //printLogAndGenome(sizeGenome, infoPPos,outSeq,outLog);
    printLogAndGenome(sizeGenome, infoPPos,outSeq,outLog, genomeRef,minQual,nameMT);

    // delete coverageCounter;














    if(deamread || useLengthPrior){
	//iterate over the reads once again and compute prob of deam
	
	cerr<<"Reading BAM to set priors for each read ..."<<endl;
	BamReader reader;
	BamAlignment al;
	if ( !reader.Open(bamfiletopen) ) {
	    cerr << "Could not open input BAM files " <<bamfiletopen<< endl;
	    return 1;
	}
	unsigned int skipped      =0;

	while ( reader.GetNextAlignment(al) ) { //for each read
	    //cout<<al.Name<<endl;
	    //those reads are ignored later anyway..
	    if(!al.IsMapped()){   continue; }
	    if( al.IsFailedQC()){ continue; }

	    pair< string,vector<int> > reconstructedReference = reconstructRefWithPos(&al);
	    
	    if( skipAlign(reconstructedReference.first,&al,&skipped) ){ 
		// cout<<al.QueryBases<<endl<<reconstructedReference.first<<endl<<endl;
		continue; 
	    }

	    if(al.QueryBases.size() != reconstructedReference.first.size()){
		cerr<<"Query bases line is not the same size as the reconstructed reference for read "<<al.Name<<endl;
		return 1;
	    }

	    if(al.QueryBases.size() != reconstructedReference.second.size()){
		cerr<<"Query bases line is not the same size as the reconstructed positions on the reference for read "<<al.Name<<endl;
		return 1;
	    }

	    // for(unsigned int i=0;i<reconstructedReference.first.size();i++){
	    // 	cout<<reconstructedReference.first[i]<<"\t"<<reconstructedReference.second[i]<<endl;
	    // }


	    double deamLogLike=0.0;
	    double nullLogLike=0.0;
	    if(deamread){
		for(unsigned int i=0;i<al.QueryBases.size();i++){//for each base
		
		    char refeBase =toupper(reconstructedReference.first[i]);
		    char readBase =toupper(         al.QueryBases[i]);
		    char q        = al.Qualities[i]-offsetQual;


		    int pos       = reconstructedReference.second[i]+1;
		    // cout<<i<<"\t"<<reconstructedReference.second[i]<<endl;
		    transformRef(&refeBase,&readBase);

		    if(refeBase == 'I'   ){
			continue;
		    }

		    if(refeBase == 'N' ){
			continue;
		    }

		    if(readBase == 'N' ){
			continue;
		    }

		    if(pos2phredgeno[ pos ].consensus == 'D'){//skip deletions
			continue;
		    }



		    
		    if(pos2phredgeno[ pos ].ref != refeBase){
			cerr<<"Query reference base is not the same for read "<<al.Name<<endl;

			cout<<pos<<"\t"<<al.QueryBases[i]<<"\t"<<reconstructedReference.first[i]<<"\t"<<refeBase<<"\t"<<readBase<<"\tR="<<pos2phredgeno[ pos ].ref<<"\tC="<<pos2phredgeno[ pos ].consensus<<"\t"<<al.Name<<endl;
			for(unsigned int j=0;j<al.QueryBases.size();j++){
			    cout<<j<<"\t"<<reconstructedReference.first[j]<<"\t"<<reconstructedReference.second[j]<<endl;
			}

			return 1;
		    }
		
		    //deam model
		

		
		    int dist5p=-1;
		    int dist3p=-1;

		    if( al.IsReverseStrand() ){
			dist5p = int(al.QueryBases.size()) - int(i)-1;
			dist3p = int(i);
		    }else{
			dist5p = int(i);
			dist3p = int(al.QueryBases.size()) - int(i)-1;
		    }
		    		    
		    probSubstition * probSubMatchDeam = &defaultSubMatch ;
		    probSubstition * probSubMatchNull = &defaultSubMatch ;

		    if(dist5p <= (int(sub5p.size()) -1)){
			probSubMatchDeam = &sub5p[ dist5p ];			
		    }
		    
		    if(dist3p <= (int(sub3p.size()) -1)){
			probSubMatchDeam = &sub3p[ dist3p ];
		    }
		    
		    //we have substitution probabilities for both... take the closest
		    if(dist5p <= (sub5p.size() -1) &&
		       dist3p <= (sub3p.size() -1) ){
			
			if(dist5p < dist3p){
			    probSubMatchDeam = &sub5p[ dist5p ];			
			}else{
			    probSubMatchDeam = &sub3p[ dist3p ];
			}
			
		    }


		    // cout<<pos<<"\t"<<al.QueryBases[i]<<"\t"<<reconstructedReference.first[i]<<"\t"<<refeBase<<"\t"<<readBase<<"\tR="<<pos2phredgeno[ pos ].ref<<"\tC="<<pos2phredgeno[ pos ].consensus<<"\t"<<al.Name<<"\t"<<dist5p<<"\t"<<dist3p<<endl;
		
		    // b   is the observed
		    // nuc is the model
		
		    int obsReadInt;
		    if( al.IsReverseStrand() ){
			obsReadInt = baseResolved2int(complement(readBase));
		    }else{
			obsReadInt = baseResolved2int(readBase);
		    }

		    double probBaseDeam = 0.0;
		    double probBaseNull = 0.0;

		    for(unsigned int nuc=0;nuc<4;nuc++){


			int dinucIndex;
			if( al.IsReverseStrand() ){
			    dinucIndex= (3-nuc)*4+obsReadInt;
			}else{
			    dinucIndex=     nuc*4+obsReadInt;
			}
			probBaseDeam +=
			    (1-pos2phredgeno[ pos ].perror[nuc])
			    *
			    //      (1-e)           *  p(sub|1-e)                         + (e)                          *  p(sub|1-e)
			    (likeMatchProb[int(q)]  * (probSubMatchDeam->s[dinucIndex] )  + (1.0 - likeMatchProb[int(q)])*(illuminaErrorsProb.s[dinucIndex]));


			probBaseNull +=
			    (1-pos2phredgeno[ pos ].perror[nuc])
			    *
			    //      (1-e)           *  p(sub|1-e)                         + (e)                          *  p(sub|1-e)
			    (likeMatchProb[int(q)]  * (probSubMatchNull->s[dinucIndex] )  + (1.0 - likeMatchProb[int(q)])*(illuminaErrorsProb.s[dinucIndex])); 

			
			// if(probBaseDeam!=probBaseNull)
			// if(al.Name=="SN7001204_0274_AH72EDADXX_R_PEdi_VS_A4961_2:1:1215:18712:56271c"){
			// 	cout<<i<<"\t"<<dnaAlphabet[nuc]<<"\t"<<obsReadInt<<"\t"<<pos2phredgeno[ pos ].perror[nuc]<<"\t5p"<<dist5p <<"\t3p"<< dist3p<<"\tdi="<<dinucIndex<<"\tsub="<< probSubMatchDeam->s[dinucIndex]<<"\tsubN="<< probSubMatchNull->s[dinucIndex]<<"\t"<<readBase<<"\t"<<probBaseDeam<<"\t"<<probBaseNull<<"\t"<<probBaseDeam/(probBaseDeam+probBaseNull)<<"\t"<<probBaseNull/(probBaseDeam+probBaseNull)<<endl;
			// }


		    }//end for each possible base

		    deamLogLike+=log(probBaseDeam);
		    nullLogLike+=log(probBaseNull);

	
		    //m_infoPPos->at(posVector).likeBaseNoindel[nuc] += 
		    // double probFinal;
		
		    // if(ignoreMQ){ //ignore MQ
		    //     probFinal = (               probBase                          );
		    // }else{
		    //     probFinal = (probMapping[m]*probBase + probMismapping[m]*0.25);
		    // }
		    

		    //null model

		}//end for all bases
	    }else{//if we use deamination priors
		deamLogLike=0.5;
		nullLogLike=0.5;
	    }
	    



	    double probDeamUnscaled = 
                                exp(deamLogLike)
		                    /
		(    exp(deamLogLike) + exp(nullLogLike) );
	    
	    //probLengthEndo[max(al.QueryBases.size(),1000)];
	    double probEndoProd;
	    double probContProd;
	    double probLengthEndoForRead;
	    
	    if(useLengthPrior){
		probLengthEndoForRead = probLengthEndo[min(int(al.QueryBases.size()),999)];
	    }else{
		probLengthEndoForRead = 0.5;
	    }

	    //if(useLengthPrior){
	    probEndoProd   =      probDeamUnscaled  *      probLengthEndoForRead;
	    probContProd   = (1.0-probDeamUnscaled) * (1.0-probLengthEndoForRead);
	    // }else{
	    // 	probEndoProd   = 0.5;
	    // 	probContProd   = 0.5;

	    // }
	    double probEndo         = 
                         		(1-contaminationPrior)*probEndoProd
		                                          /
		(  ((1-contaminationPrior)*probEndoProd) + contaminationPrior*probContProd );

	    //double probEndo         = boost::math::ibeta(contaminationPrior,1.0-contaminationPrior,probDeamUnscaled);
	    //cout<<al.Name<<"\t"<<int(al.QueryBases.size())<<"\t"<<probEndo<<"\t"<<probDeamUnscaled<<"\t"<<probLengthEndoForRead<<"\t"<<probEndoProd<<"\t"<<probContProd<<endl;


	    //  if(al.Name=="SN7001204_0274_AH72EDADXX_R_PEdi_VS_A4961_2:1:1215:18712:56271c"){
	    // 	cout<<deamLogLike<<"\t"<<nullLogLike<<endl;
	    // 	cout<<al.Name<<"\t"<<probEndo<<endl;
	    // 	return 1;
	    // }

	     read2endoProb[ al.Name+"#"+ stringify(al.AlignmentFlag) ] = probEndo;

	} //for each read
	cerr<<"...  done"<<endl;
	read2endoProbInit=true;
	








	//iterate over the reads once again
	
	//iterateOverReads();
	// cerr<<"Reading BAM to find to compute probability of deamination ..."<<endl;

	iterateOverReads(fastaFile,
			 bamfiletopen,
			 &infoPPos,
			 sizeGenome,
			 ignoreMQ );
	

	


	//printLogAndGenome(sizeGenome, infoPPos,outSeq,outLog);
	printLogAndGenome(sizeGenome, infoPPos,outSeq,outLog, genomeRef,minQual,nameMT);






	// // fastaReference;
	// if ( !fastaReference.Open(fastaFile , fastaFile+".fai") ){ 
	//     cerr << "ERROR: failed to open fasta file " <<fastaFile<<" and index " << fastaFile<<".fai"<<endl;
	//     return false;
	// }

	// if ( !reader.Open(bamfiletopen) ) {
	//     cerr << "Could not open input BAM files " <<bamfiletopen<< endl;
	//     return 1;
	// }
	
	// if ( !reader.OpenIndex(bamfiletopen+".bai") ) {
	//     cerr << "Could not open input index for BAM files" <<bamfiletopen+".bai"<< endl;
	//     return 1;
	// }
     
	// //  // retrieve reference data
	// const RefVector  	references2 = reader.GetReferenceData();

	//  cv = new MyPileupVisitor(references2,&fastaReference,&infoPPos,sizeGenome,ignoreMQ);
	//  //pileup;
	//  pileup.AddVisitor(cv);


	// // BamAlignment al;
	// //	unsigned int 
	// numReads=0;
	// cerr<<"Reading BAM file ..."<<endl;

	// while ( reader.GetNextAlignment(al) ) {
	//     //cerr<<al.Name<<endl;
	//     numReads++;
	//     if(numReads !=0 && (numReads%100000)==0){
	// 	cerr<<"Read "<<thousandSeparator(numReads)<<" reads"<<endl;
	//     }

	//     if(al.IsMapped() && !al.IsFailedQC()){
	// 	// cerr<<al.Name<<endl;
	// 	pileup.AddAlignment(al);
	//     }
	
	// }
	// cerr<<"...  done"<<endl;

	// //clean up
	// pileup.Flush();

	// reader.Close();
	// fastaReference.Close();
	// delete cv;

    }


    return 0;
}

