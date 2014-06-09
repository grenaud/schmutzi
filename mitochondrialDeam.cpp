/*
 * mitonchondrialDeam
 * Date: Mar-26-2014 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */




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

#include "utils.h"
#include "miscfunc.h"

using namespace BamTools;
using namespace std;

#define MAXCOV 5000

// #define DEBUG1
// #define DEBUG2


typedef struct { 
    double likeBaseNoindel[4];
    int  covPerBase[4];
    double mapqAvg;
    
    int numDel;
    vector<string> insertionRight;
    int cov;
} positionInformation;


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

probSubstition illuminaErrorsProb;
vector<probSubstition> sub5p;
vector<probSubstition> sub3p;
probSubstition defaultSubMatch;

string dnaAlphabet="ACGT";

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
		    int dinucIndex = nuc*4+baseResolved2int(b);
		    
		    //                        (1-e)           *  p(sub|1-e)                          + (e) *  p(sub|1-e)
		    double probBase =  likeMatchProb[int(q)]  * (probSubMatchToUse->s[dinucIndex] )  + (1.0 - likeMatchProb[int(q)])*(illuminaErrorsProb.s[dinucIndex]);
		    //m_infoPPos->at(posVector).likeBaseNoindel[nuc] += 
		    double probFinal;

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
    vector<positionInformation> * m_infoPPos;
    int sizeGenome;
    bool ignoreMQ;
    //        ostream*  m_out;
};
























int main (int argc, char *argv[]) {

    int sizeGenome=0;
    string outSeq  = "/dev/stdout";
    string outLog  = "/dev/stderr";
    string nameMT  = "MT";

    ofstream outSeqFP ;
    ofstream outLogFP;
    int minQual=0;
    bool ignoreMQ=false;


    ////////////////////////////////////
    // BEGIN Initializing scores      //
    ////////////////////////////////////
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

    ////////////////////////////////////
    //    END Initializing scores     //
    ////////////////////////////////////






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

    const string usage=string("\t"+string(argv[0])+
			      " [options]  [reference fasta] [bam file] "+"\n\n"+
			       "\n\tOutput options:\n"+	
			      "\t\t"+"-seq  [fasta file]" +"\t\t"+"Output fasta file (default: stdout)"+"\n"+
			      "\t\t"+"-log  [log file]" +"\t\t"+"Output log  (default: stderr)"+"\n"+
			      "\t\t"+"-name [name]" +"\t\t\t"  +"Name  (default "+nameMT+") "+"\n"+
			      "\t\t"+"-qual [minimum quality]" +"\t\t"  +"Filter bases with quality less than this  (default "+stringify(minQual)+") "+"\n"+
			      // "\t\t"+"-cont" +"\t\t"+"Contamination allele frequency"+"\n"+
			      "\t\t"+"-deam5p [.prof file]" +"\t\t"+"5p deamination frequency (default: "+deam5pfreq+")"+"\n"+
			      "\t\t"+"-deam3p [.prof file]" +"\t\t"+"3p deamination frequency (default: "+deam3pfreq+")"+"\n"+
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
	cout<<"Usage:"<<endl;
	cout<<""<<endl;
	cout<<usage<<endl;
	return 1;
    }

    string bamfiletopen = string(argv[argc-1]);//bam file
    string fastaFile    = string(argv[argc-2]);//fasta file

    for(int i=1;i<(argc-2);i++){ //all but the last 3 args
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


	// if(strcmp(argv[i],"-cont") == 0 ){
	//     contFile=string(argv[i+1]);
	//     i++;
	//     continue;
	// }

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



    outSeqFP.open(outSeq.c_str());

    if (!outSeqFP.is_open()){
	cerr << "Unable to write to seq file "<<outSeq<<endl;
	return 1;
    }

    outLogFP.open(outLog.c_str());

    if (!outLogFP.is_open()){
	cerr << "Unable to write to qual file "<<outLog<<endl;
	return 1;
    }














     
    vector<positionInformation> infoPPos;

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

     for(int i=0;i<sizeGenome;i++){
	 positionInformation toadd;
	 toadd.numDel         = 0;

	 toadd.cov            = 0;
	 toadd.mapqAvg        = 0.0;

	 for(unsigned int nuc=0;nuc<4;nuc++){
	     toadd.likeBaseNoindel[nuc] = 0;
	     toadd.covPerBase[nuc]      = 0;
	 }

	 infoPPos.push_back(toadd);
     }


     // double ** likelihoodPBasePPos = new double * [genome.size()];
     // for(unsigned int i = 0; i< genome.size(); i++) {
     // 	likelihoodPBasePPos[i] = new double  [4];
     // }
     // //return 1;

     Fasta fastaReference;
     if ( !fastaReference.Open(fastaFile , fastaFile+".fai") ){ 
	 cerr << "ERROR: failed to open fasta file " <<fastaFile<<" and index " << fastaFile<<".fai"<<endl;
	 return false;
     }



    BamReader reader;

    if ( !reader.Open(bamfiletopen) ) {
	cerr << "Could not open input BAM files " <<bamfiletopen<< endl;
    	return 1;
    }



    //  if ( !reader.OpenIndex(bamfiletopen+".bai") ) {
    // 	 cerr << "Could not open input index for BAM files" <<bamfiletopen+".bai"<< endl;
    // 	 return 1;
    //  }
     
    //  // retrieve reference data
    const RefVector  references = reader.GetReferenceData();


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

    MyPileupVisitor* cv = new MyPileupVisitor(references,&fastaReference,&infoPPos,sizeGenome,ignoreMQ);
    PileupEngine pileup;
    pileup.AddVisitor(cv);


    BamAlignment al;
    unsigned int numReads=0;
    while ( reader.GetNextAlignment(al) ) {
	//cerr<<al.Name<<endl;
	numReads++;

	if(al.IsMapped() && !al.IsFailedQC()){
	    // cerr<<al.Name<<endl;
	    pileup.AddAlignment(al);
	}
	
    }
    //    cerr<<"done"<<endl;

    

    
    // for(unsigned int i=0;i<positionsInfoFound->size();i++){

    // 	char alt = positionsInfoFound->at(i).alt;
    // 	if(alt == 'N')
    // 	    alt = returnRandomNuc(positionsInfoFound->at(i).ref,dnaEvolCumu);

    // 	cout<<"#"<<"\t"<<positionsInfoFound->at(i).ref<<"\t"
    // 	    <<positionsInfoFound->at(i).pos<<"\t"
    // 	    <<alt<<"\t"
    // 	    <<(int)positionsInfoFound->at(i).cov<<"\t"
    // 	    <<positionsInfoFound->at(i).contFreq<<"\t"<<
    // 	    arrayToStringInt(positionsInfoFound->at(i).mapq,(int)positionsInfoFound->at(i).cov)<<"\t"<<
    // 	    arrayToString(positionsInfoFound->at(i).refOrAlt,(int)positionsInfoFound->at(i).cov)<<"\t"<<
    // 	    arrayToStringInt(positionsInfoFound->at(i).quals,(int)positionsInfoFound->at(i).cov)<<"\t"<<
    // 	    arrayToStringInt(positionsInfoFound->at(i).dist5p,(int)positionsInfoFound->at(i).cov)<<"\t"<<
    // 	    arrayToStringInt(positionsInfoFound->at(i).dist3p,(int)positionsInfoFound->at(i).cov)<<endl;


    // 	//heterozygous
	

    // 	//<<"\t"<<(int)positionsInfoFound->at(i).pos<<endl;

    //  	 // cout<<i<<"\t"<<coverageCounter->at(i)<<endl;
    // 	 // total     += coverageCounter->at(i);
    // 	 // totalBase += coverageCounter->at(i)*i;
    //  }
    
    // unsigned int total=0;
    // unsigned int totalBase=0;

    //  for(unsigned int i=0;i<MAXCOV;i++){
    //  	 cout<<i<<"\t"<<coverageCounter->at(i)<<endl;
    // 	 total     += coverageCounter->at(i);
    // 	 totalBase += coverageCounter->at(i)*i;
    //  }
    //  cerr<<"avg:"<<double(totalBase)/double(total)<<"\tnum reads:"<<numReads<<endl;

    //clean up
    pileup.Flush();

    reader.Close();
    fastaReference.Close();



    string genomeToPrint="";
    outLogFP<<"pos\trefBase\tbase\tqual\tavgmapq\tcov\tsupp\tpa\tpc\tpg\tpt\n";
    //genomeRef
    for(int i=0;i<sizeGenome;i++){

	 //if half of the reads support an deletion in reads (ins in reference) 
	 //move to next position since deletion
	 if( (double(infoPPos[i].numDel)/double(infoPPos[i].cov)) >= 0.5){
	     //min quality ? based on what ?

	     outLogFP<<(i+1)<<"\t"<<genomeRef[i]<<"\t"<<"D"<<"\t"<<-10.0*(log(1.0-(double(infoPPos[i].numDel)/double(infoPPos[i].cov)))/log(10.0))<<"\t"<<infoPPos[i].mapqAvg<<"\t"<<infoPPos[i].cov<<"\t"<<infoPPos[i].numDel<<"\t0.0\t0.0\t0.0\t0.0"<<endl;
	     continue;
	 }


	 
	 if(infoPPos[i].cov == 0){
	     genomeToPrint+="N";
	     outLogFP<<(i+1)<<"\t"<<genomeRef[i]<<"\tN\t0\t0\t0\t0\t0.0\t0.0\t0.0\t0.0"<<endl;
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
	 

	 //cout<<dnaAlphabet[bestNuc]<<"\t"<<sumLogLikeOnlyBest<<"\t"<<sumLogLikeAllButBest<<"\t"<<sumLogLikeAll<<"\tnum:"<<pow(10.0,sumLogLikeOnlyBest)<<"\tnum:"<<pow(10.0,sumLogLikeAllButBest)<<"\tall:"<<pow(10.0,sumLogLikeAll)<<"\t"<<-10*(sumLogLikeAllButBest-sumLogLikeAll)<<endl;
	 //cout<<dnaAlphabet[bestNuc]<<"\t"<<(sumLogLikeOnlyBest/sumLogLikeAll)<<"\t"<<-10.0*(log(1.0-(sumLogLikeOnlyBest/sumLogLikeAll))/log(10.0))<<endl;
	 //outSeqFP<<dnaAlphabet[bestNuc];
	 if( (-10.0*(sumLogLikeAllButBest-sumLogLikeAll)) >= minQual){
	     genomeToPrint+=dnaAlphabet[bestNuc];
	 }else{
	     genomeToPrint+="N";
	 }

	 outLogFP<<(i+1)<<"\t"<<genomeRef[i]<<"\t"<<dnaAlphabet[bestNuc]<<"\t"<<(-10*(sumLogLikeAllButBest-sumLogLikeAll)+0.0)<<"\t"<<infoPPos[i].mapqAvg<<"\t"<<infoPPos[i].cov<<"\t"<<infoPPos[i].covPerBase[bestNuc];//<<endl;
	 for(int nuc=0;nuc<4;nuc++){
	     outLogFP<<"\t"<<(-10*(sumLogForNucs[nuc]-sumLogLikeAll)+0.0);
	 }
	 outLogFP<<endl;

	 //outLogFP<< char(-10.0*(log(sumLikeAllButBest/sumLikeAll)/log(10.0)) +offsetQual);
	 //insertions in the reads (deletion in reference)
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
	     // cout<<mostCommonIns<<endl;
	     // cout<<mostCommonInsCount<<endl;
	     // cout<<(double(mostCommonInsCount)/double(infoPPos[i].cov))<<endl;


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

		 //outLogFP<<string('!',mostCommonIns.size());
		 for(unsigned int k=0;k<(mostCommonIns.size());k++){
		     outLogFP<<(i+1)<<"i\t"<<"-"<<"\t"<<mostCommonIns[k]<<"\t"<<-10.0*(log(1.0-double(mostCommonInsCount)/double(infoPPos[i].cov))/log(10.0))<<"\t"<<infoPPos[i].mapqAvg<<"\t"<<infoPPos[i].cov<<"\t"<<mostCommonInsCount<<"\t0.0\t0.0\t0.0\t0.0"<<endl;

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
     outLogFP.close();

    delete cv;
    // delete coverageCounter;



    return 0;
}

