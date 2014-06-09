/*
 * mtCont
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


using namespace BamTools;
using namespace std;

#define MAXCOV 5000

// #define DEBUG1
// #define DEBUG2

char   offsetQual=33;
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

string dnaAlphabet="ACGT";

#define MIN(a,b) (((a)<(b))?(a):(b))


// // Returns log10( pow(10,x)+pow(10,y) ), but does so without causing
// // overflow or loss of precision.
// double oplus( double x, double y )
// {
//     return x > y 
//         ? x + log1p( pow( 10, y-x ) ) / log(10)
//         : y + log1p( pow( 10, x-y ) ) / log(10) ;
// }

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



typedef struct { 
    char base;
    int  qual;
    int  mapq;       	 
} singleRead;

typedef struct { 
    vector<singleRead> readsVec;
    double mapqAvg;
    int cov;
    char refBase;
    int posAlign;
} positionInformation;


typedef struct { 
    double f[4];
 } alleleFrequency;


typedef struct { 
    double perror[4];
    double phred[4];
    char consensus;
 } PHREDgeno;



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
	    //cout<<endl<<"pos = "<<posAlign<<"\t"<<posVector;

	    if( (posAlign%100) == 0){
		//cerr<<"pos  = "<<posAlign<<endl;
	    }
	    //for some reason, we have to do -1 on the .Position
	    if ( !m_fastaReference->GetBase(pileupData.RefId, posAlign-1, referenceBase ) ) {
		cerr << "bamtools convert ERROR: pileup conversion - could not read reference base from FASTA file" << endl;
		exit(1);
	    }


	    // cout<<""<<referenceBase<<"\t"<<endl;
	    


	    //insertion in the reads/deletion in the reference, skipping
	    for(unsigned int i=0;i<pileupData.PileupAlignments.size();i++){		
		if( !pileupData.PileupAlignments[i].IsCurrentDeletion &&
		    pileupData.PileupAlignments[i].IsNextInsertion &&
		    (pileupData.PileupAlignments[i].InsertionLength>0)){
		    continue;
		}
	    }


	    //deletion in the reads/insertion in the reference
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
		m_infoPPos->at(posVector).cov++;		    

		//Add mapq
		m_infoPPos->at(posVector).mapqAvg += pow(10.0, double(pileupData.PileupAlignments[i].Alignment.MapQuality)/-10.0);
		m_infoPPos->at(posVector).readsVec.push_back(sr);

		// numReads++;
		if(numReads >= MAXCOV){
		    break;
		}

	    }//end each read
	    
	    m_infoPPos->at(posVector).mapqAvg  = m_infoPPos->at(posVector).mapqAvg/double(m_infoPPos->at(posVector).cov);
	    m_infoPPos->at(posVector).mapqAvg  = -10.0*( log( m_infoPPos->at(posVector).mapqAvg )/log(10.0) );
	    m_infoPPos->at(posVector).refBase  = referenceBase;
	    m_infoPPos->at(posVector).posAlign = posAlign;

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
    // string outSeq  = "/dev/stdout";
    string outLog  = "/dev/stderr";
    // string nameMT  = "MT";
    string fileFreq = "alleleFreqMT/1000g/freqHumans.dat";

    // ofstream outSeqFP ;
    ofstream outLogFP;
    // int minQual=0;
    bool ignoreMQ=false;
    string line;

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

	//m = prob of mismapping
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

    //    return 1;

    const string usage=string("\t"+string(argv[0])+
			      " [options]  [consensus file] [reference fasta file] [bam file] "+"\n\n"+
			       "\n\tOutput options:\n"+	
			      
			      // "\t\t"+"-seq  [fasta file]" +"\t\t"+"Output fasta file (default: stdout)"+"\n"+
			      // "\t\t"+"-log  [log file]" +"\t\t"+"Output log  (default: stderr)"+"\n"+
			      // "\t\t"+"-name [name]" +"\t\t\t"  +"Name  (default "+nameMT+") "+"\n"+
			      // "\t\t"+"-qual [minimum quality]" +"\t\t"  +"Filter bases with quality less than this  (default "+stringify(minQual)+") "+"\n"+
			      
			      // // "\t\t"+"-deam5" +"\t\t"+"5p deamination frequency"+"\n"+
			      // // "\t\t"+"-deam3" +"\t\t"+"3p deamination frequency"+"\n"+
			      "\n\tComputation options:\n"+	
			       "\t\t"+"-nomq" +"\t\t\t"+"Ignore mapping quality (default: "+booleanAsString(ignoreMQ)+")"+"\n"+
			     "\t\t"+"-cont [file1,file2,...]" +"\t\t"+"Contamination allele frequency (default : "+fileFreq+")"+"\n"+ 
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

    string bamfiletopen     = string(argv[argc-1]);//bam file
    string fastaFile        = string(argv[argc-2]);//fasta file
    string consensusFile    = string(argv[argc-3]);//fasta file

    for(int i=1;i<(argc-3);i++){ //all but the last 3 args
	
	if(strcmp(argv[i],"-cont") == 0 ){
	    fileFreq=string(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-nomq") == 0 ){
	    ignoreMQ=true;
	    continue;
	}

	// if(strcmp(argv[i],"-deam5") == 0 ){
	//     deam5File=string(argv[i+1]);
	//     i++;
	//     continue;
	// }


	// if(strcmp(argv[i],"-seq") == 0 ){
	//     outSeq=string(argv[i+1]);
	//     i++;
	//     continue;
	// }

	// if(strcmp(argv[i],"-log") == 0 ){
	//     outLog=string(argv[i+1]);
	//     i++;
	//     continue;
	// }

	// if(strcmp(argv[i],"-qual") == 0 ){
	//     minQual=destringify<int>(argv[i+1]);
	//     i++;
	//     continue;
	// }

	// if(strcmp(argv[i],"-name") == 0 ){
	//     nameMT=string(argv[i+1]);
	//     i++;
	//     continue;
	// }



	// if(strcmp(argv[i],"-l") == 0 ){
	//     sizeGenome=atoi(argv[i+1]);
	//     i++;
	//     continue;
	// }
	cerr<<"Wrong option "<<string(argv[i])<<endl;
	return 1;


    }




    // outSeqFP.open(outSeq.c_str());

    // if (!outSeqFP.is_open()){
    // 	cerr << "Unable to write to seq file "<<outSeq<<endl;
    // 	return 1;
    // }

    outLogFP.open(outLog.c_str());

    if (!outLogFP.is_open()){
	cerr << "Unable to write to qual file "<<outLog<<endl;
	return 1;
    }










    vector<string> allFreqFiles  = allTokens(fileFreq,',');
    vector< map<int, alleleFrequency>  > vecAlleleFreq ;

    for(unsigned int ffreqf=0;ffreqf<allFreqFiles.size();ffreqf++){
	// cerr<<"ffreqf "<<ffreqf<<endl;    
	map<int, alleleFrequency> pos2allelefreq;
	string genomeRef="";

	igzstream freqAlleleFile;
	freqAlleleFile.open(allFreqFiles[ffreqf].c_str());
	if (freqAlleleFile.good()){

	    while ( getline (freqAlleleFile,line)){

		vector<string> fields = allTokens(line,'\t');
		alleleFrequency freqToadd;
	    
		if(fields.size() != 5){
		    cerr << "line "<<line<<"  in file  "<<freqAlleleFile<<" does not have 5 fields"<<endl;
		    return 1;
		}
	    
		for(int nuc=0;nuc<4;nuc++){
		    freqToadd.f[nuc]=destringify<double>(fields[nuc+1]);
		}

		pos2allelefreq[ destringify<int>( fields[0])  ] = freqToadd;
	    	    
	    }
	    freqAlleleFile.close();

	}else{
	    cerr << "Cannot open allele frequency file  "<<fileFreq<<""<<endl;
	    return 1;
	}

	// cerr<<"ffreqfe "<<ffreqf<<endl;    

	vecAlleleFreq.push_back(pos2allelefreq);
    }
    
    // cerr<<"ffreqfe "<<endl;    






    // typedef struct { 
    //     double perror[4];
    //     double phred[4];
    //  } PHREDgeno;
    map<int, PHREDgeno> pos2phredgeno;

    igzstream consensusFD;
    consensusFD.open(consensusFile.c_str());
    if (consensusFD.good()){
	getline (consensusFD,line);

	while ( getline (consensusFD,line)){

	    vector<string> fields = allTokens(line,'\t');
	    PHREDgeno toadd;
	    // cerr<<line<<endl;


	    if(fields.size() != 11){
		cerr << "line "<<line<<"  in file  "<<consensusFile<<" does not have 11 fields"<<endl;
		return 1;
	    }
	    
	    toadd.consensus = fields[2][0];
	    for(int nuc=0;nuc<4;nuc++){
		
		toadd.phred[nuc]  = destringify<double>(fields[nuc+7]);		
		toadd.perror[nuc] = pow(10.0,toadd.phred[nuc]/(-10.0));
		
	    }

	    pos2phredgeno[     destringify<int>( fields[0])   ] = toadd;
	    sizeGenome =  max( destringify<int>( fields[0]), sizeGenome);
	}
	consensusFD.close();

    }else{
	cerr << "Cannot open consensus file  "<<consensusFile<<""<<endl;
	return 1;
    }



    // cerr<<"cons "<<endl;    




    // typedef struct { 
    //     vector<singleRead> readsVec;
    //     double mapqAvg;
    //     int cov;
    // } positionInformation;
    vector<positionInformation> infoPPos;
    for(int i=0;i<=sizeGenome;i++){
	 positionInformation toadd;

	 toadd.cov            = 0;
	 toadd.mapqAvg        = 0.0;
	 
	 infoPPos.push_back(toadd);
     }
    //return 1;




    // for(int nuc=0;nuc<4;nuc++){
    // 	cout<<nuc<<"\t"<<alleFreqVec[1].f[nuc]<<endl;
    // }
    // return 1;

    // if(sizeGenome == 0){
    // 	sizeGenome=genomeRef.size();//use the one from the fasta
    // }

    //  for(int i=0;i<sizeGenome;i++){
    // 	 positionInformation toadd;
    // 	 toadd.numDel         = 0;

    // 	 toadd.cov            = 0;
    // 	 toadd.mapqAvg        = 0.0;

    // 	 for(unsigned int nuc=0;nuc<4;nuc++){
    // 	     toadd.likeBaseNoindel[nuc] = 0;
    // 	     toadd.covPerBase[nuc]      = 0;
    // 	 }

    // 	 infoPPos.push_back(toadd);
    //  }


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


    // cerr<<"Begin reading BAM file"<<endl;

    MyPileupVisitor* cv = new MyPileupVisitor(references,&fastaReference,&infoPPos,sizeGenome,ignoreMQ);
    PileupEngine pileup;
    pileup.AddVisitor(cv);


    BamAlignment al;
    unsigned int numReads=0;
    while ( reader.GetNextAlignment(al) ) {
	// cerr<<al.Name<<endl;
	numReads++;
	if(numReads !=0 && (numReads%100000)==0){
	    cerr<<"Read "<<thousandSeparator(numReads)<<" reads"<<endl;
	}

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



    // cerr<<"... done reading BAM file"<<endl;
    

    // string genomeToPrint="";
    // outLogFP<<"pos\trefBase\tbase\tqual\tavgmapq\tcov\tsupp\n";
    // //genomeRef

    //for each file in cont
    for(unsigned int ffreqf=0;ffreqf<allFreqFiles.size();ffreqf++){
	// vector<string> allFreqFiles     = allTokens(fileFreq,',');
	// vector< map<int, alleleFrequency>  > vecAlleleFreq ;
	// pos2allelefreq
	double contaminationRate=0.0;


	for(contaminationRate=0.0;contaminationRate<0.5;contaminationRate+=0.005){
	    double logLike=0.0;
    
	    for(int i=0;i<sizeGenome;i++){
		//for(int i=261;i<=262;i++){

#ifdef DEBUG2
		cout<<"cov "<<infoPPos[i].cov<<endl;
#endif


		if(infoPPos[i].cov == 0){
		    continue;
		}

#ifdef DEBUG2
		cout<<"\npos="<<infoPPos[i].posAlign<<endl;
		cout<<"freq "<<vecAlleleFreq[ffreqf][ infoPPos[i].posAlign ].f[0]<<"\t"<<vecAlleleFreq[ffreqf][ infoPPos[i].posAlign ].f[1]<<"\t"<<vecAlleleFreq[ffreqf][ infoPPos[i].posAlign ].f[2]<<"\t"<<vecAlleleFreq[ffreqf][ infoPPos[i].posAlign ].f[3]<<endl;
		cout<<"cons "<<pos2phredgeno[ infoPPos[i].posAlign ].consensus<<"\t"<<pos2phredgeno[ infoPPos[i].posAlign ].perror[0]<<"\t"<<pos2phredgeno[ infoPPos[i].posAlign ].perror[1]<<"\t"<<pos2phredgeno[ infoPPos[i].posAlign ].perror[2]<<"\t"<<pos2phredgeno[ infoPPos[i].posAlign ].perror[3]<<endl;
#endif

		double priorDiNuc[4][4];
	 
		//computing prior
		//                 p(nuc1) is the prob. of consensus                  *  p(contaminant)
		for(unsigned int nuc1=0;nuc1<4;nuc1++){ //     b = consensus
		    for(unsigned int nuc2=0;nuc2<4;nuc2++){//  c = contaminant
			double priortemp = (1-pos2phredgeno[ infoPPos[i].posAlign ].perror[nuc1]) * vecAlleleFreq[ffreqf][ infoPPos[i].posAlign ].f[nuc2];
			priorDiNuc[nuc1][nuc2] = priortemp;		 
		    }
		}

#ifdef DEBUG2
		for(unsigned int nuc1=0;nuc1<4;nuc1++){ //     b = consensus
		    for(unsigned int nuc2=0;nuc2<4;nuc2++){//  c = contaminant
			cout<<"cons="<<dnaAlphabet[nuc1]<<"\tcont="<<dnaAlphabet[nuc2]<<"\tprior="<<priorDiNuc[nuc1][nuc2]<<endl;		 
		    }
		}
#endif

		// continue;
		for(unsigned int k=0;k<infoPPos[i].readsVec.size();k++){ //for every read at that position
		    double mappedProb  =0.0; //initialized
		    double misappedProb=0.25;
	   
		    int baseIndex = base2int(infoPPos[i].readsVec[k].base)-1;
		    int qual      = infoPPos[i].readsVec[k].qual;
#ifdef DEBUG2
		    bool potentialCont=false;
#endif
		    //iterate over each possible consensus and contaminant base
		    for(int nuc1=0;nuc1<4;nuc1++){ //     b = consensus
			for(int nuc2=0;nuc2<4;nuc2++){//  c = contaminant
			    //skip when contaminant is the consensus
			    if(nuc1 == nuc2)
				continue;

			    double probNonCont;		     
			    if(baseIndex == nuc1){ //matches the consensus
				probNonCont  = likeMatchProb[qual];
			    }else{
				probNonCont  = likeMismatchProb[qual];
			    }

			    double probCont;
			    if(baseIndex == nuc2){ //matches the contaminant
				probCont  = likeMatchProb[qual];
			    }else{
				probCont  = likeMismatchProb[qual];
			    }

		     
			    double probForDinuc =  priorDiNuc[nuc1][nuc2]*
				( ( (1.0-contaminationRate) * probNonCont) +
				  ( (contaminationRate)     * probCont   ) ) ;
		     
			    mappedProb+=probForDinuc;

#ifdef DEBUG2
			    cout<<"b="<<dnaAlphabet[nuc1]<<"\tc="<<dnaAlphabet[nuc2]<<"\tp(b,c)="<<priorDiNuc[nuc1][nuc2]<<"\tp(noncont)="<<probNonCont<<"\tp(cont)="<<probCont<<"\tp(read)="<<probForDinuc<<"\tbase="<<dnaAlphabet[baseIndex]<<"\tqual="<<qual<<"\tp(map)="<<mappedProb<<endl;

			    if(priorDiNuc[nuc1][nuc2] >0.9){
				potentialCont=true;
				// 	 //return 1;
			    }
#endif

			}
		    }//end for each di-nucleotide


		    double pread = (probMapping[infoPPos[i].readsVec[k].mapq]*(mappedProb) + probMismapping[infoPPos[i].readsVec[k].mapq]*(misappedProb) );
#ifdef DEBUG2     
		    if(potentialCont){
			cout<<"potential "<<endl;
		    }

		    cout<<"c="<<contaminationRate<<" mapq = "<<infoPPos[i].readsVec[k].mapq<<"\tprio(map)="<<probMapping[infoPPos[i].readsVec[k].mapq]<<"\tp(map)="<<mappedProb<<"\tprio(nmap)="<<probMismapping[infoPPos[i].readsVec[k].mapq]<<"\tp(nmap)="<<misappedProb<<"\tpread="<<pread<<endl;
#endif


#ifdef DEBUG2     
		    cout<<"logLike before="<<logLike<<"\tadd="<<log  (   pread   )<<endl;
#endif
		

		    // }
		    // m = prob. of mismapping
		    //           (1-m)   * p(data|mapped)                             +      m  * p(data|mismapped)		
		    logLike += log  (   pread   );

#ifdef DEBUG2     
		    cout<<"logLike "<<logLike<<endl<<endl;
#endif

		} //end for each read at that position
	    
	    } //end for each position in the genome

	    //cout.precision(15);
	    cout<<allFreqFiles[ffreqf]<<"\t"<<contaminationRate<<"\t"<<logLike<<endl;

	}//end each contaminationRate
    }//end for each cont freq file

    outLogFP.close();

    delete cv;



    return 0;
}

