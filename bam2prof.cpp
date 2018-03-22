#include <iostream>
#include <vector>
#include <set>
#include <ctype.h>
#include <stdlib.h>

#include <api/BamMultiReader.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <api/BamAux.h>
#include <api/SamSequenceDictionary.h>

#include "utils.h"
#include "ReconsReferenceBAM.h"

using namespace std;
using namespace BamTools;

const int offset=33;
int numberOfCycles;

#define MAXLENGTH 1000


#define MAX(a,b) (((a)>(b))?(a):(b))


// vector<unsigned int>    matches;
// vector<unsigned int> mismatches;

vector< vector<unsigned int> > typesOfDimer5p; //5' deam rates
vector< vector<unsigned int> > typesOfDimer3p; //3' deam rates
vector< vector<unsigned int> > typesOfDimer5pDouble; //5' deam rates when the 3' is deaminated according to a double str.
vector< vector<unsigned int> > typesOfDimer3pDouble; //3' deam rates when the 5' is deaminated according to a double str.
vector< vector<unsigned int> > typesOfDimer5pSingle; //5' deam rates when the 3' is deaminated according to a single str.
vector< vector<unsigned int> > typesOfDimer3pSingle; //3' deam rates when the 5' is deaminated according to a single str.


//increases the counters mismatches and typesOfMismatches of a given BamAlignment object
inline void increaseCounters(const BamAlignment & al,string & reconstructedReference, const int & minQualBase){ // ,int firstCycleRead,int increment

    char refeBase;
    char readBase;
    int  qualBase;
    //    int cycleToUse=firstCycleRead;
    // cout<<"name "<<al.Name<<endl;
    // cout<<"firstCycleRead "<<firstCycleRead<<endl;
    // cout<<"increment      "<<increment<<endl;

    //Checking if the 5' is deaminated
    bool isDeam5pS=false; //C->T 5'
    bool isDeam3pS=false; //C->T 3'
    bool isDeam5pD=false; //C->T 5'
    bool isDeam3pD=false; //G->A 3'

    int i;


    
    i=0; //5p for forward str, 3p for reverse
    refeBase=toupper(reconstructedReference[i]);
    readBase=toupper(         al.QueryBases[i]);
    qualBase=int(             al.Qualities[i])-offset;
    
    if(qualBase < minQualBase)
	goto eval3pdeam;

    if(refeBase == 'S' ||refeBase == 'I'){ //don't care about soft clipped or indels
	goto eval3pdeam;
    }
    
    if(refeBase == 'M'){//match
	refeBase =  readBase;
    }

    if( isResolvedDNA(refeBase)  && 
	isResolvedDNA(readBase) ){
	if(al.IsReverseStrand()){ //need to take the complement
	    refeBase=complement(refeBase);
	    readBase=complement(readBase);
	}
	
	if(refeBase == 'C' &&
	   readBase == 'T' ){ //C->T

	    if(al.IsReverseStrand()){ //3'
		isDeam3pS=true;		
	    }else{                    //5'
		isDeam5pS=true;
		isDeam5pD=true;
	    }
	}


	if(refeBase == 'G' &&
	   readBase == 'A' ){ //G->A

	    if(al.IsReverseStrand()){ //3'
		isDeam3pD=true;		
	    }else{                    //5'
	    }
	}
	   
    }


 eval3pdeam:
    i=int(al.QueryBases.size())-1; //3p for forward str, 5p for reverse
    refeBase=toupper(reconstructedReference[i]);
    readBase=toupper(         al.QueryBases[i]);
    qualBase=int(              al.Qualities[i])-offset;
    
    if(qualBase < minQualBase)
	goto iterateLoop;
    
    if(refeBase == 'S' ||refeBase == 'I'){ //don't care about soft clipped or indels
	goto iterateLoop;
    }
    
    if(refeBase == 'M'){//match
	refeBase =  readBase;
    }

    if( isResolvedDNA(refeBase)  && 
	isResolvedDNA(readBase) ){
	if(al.IsReverseStrand()){ //need to take the complement
	    refeBase=complement(refeBase);
	    readBase=complement(readBase);
	}
	
	if(refeBase == 'C' &&
	   readBase == 'T' ){ //C->T

	    if(al.IsReverseStrand()){ //5'
		isDeam5pS=true;
		isDeam5pD=true;		
	    }else{                    //3'
		isDeam3pS=true;
	    }
	}

	if(refeBase == 'G' &&
	   readBase == 'A' ){ //G->A

	    if(al.IsReverseStrand()){ //5'
	    }else{                    //3'
		isDeam3pD=true;
	    }
	}
	   


    }

 iterateLoop:


    for(i=0;i<int(al.QueryBases.size());i++){
	//cout<<i<<endl;
	refeBase=toupper(reconstructedReference[i]);
	readBase=toupper(         al.QueryBases[i]);
	qualBase=int(              al.Qualities[i])-offset;
	//cout<<i<<"\t"<<qualBase<<"\t"<<minQualBase<<endl;
	if(qualBase < minQualBase)
	    continue;
	//cout<<"-"<<i<<"\t"<<qualBase<<"\t"<<minQualBase<<endl;		     
	if(refeBase == 'S' ||refeBase == 'I'){ //don't care about soft clipped or indels
	    continue;
	}

	if(refeBase == 'M'){//match
	    refeBase =  readBase;
	}

	if( isResolvedDNA(refeBase)  && 
	    isResolvedDNA(readBase) ){
	    int dist5p=i;
	    int dist3p=int(al.QueryBases.size())-1-i;
	    
	    if(al.IsReverseStrand()){ //need to take the complement
		refeBase=complement(refeBase);
		readBase=complement(readBase);
		dist5p=int(al.QueryBases.size())-1-i;
		dist3p=i;
	    }

	    if(dist5p > MAXLENGTH ||
	       dist3p > MAXLENGTH ){
		cerr<<"Molecule found "<<al.Name<<" with length greater than limit"<<endl;
		exit(1);;
	    }
	       

	    //mismatches[cycleToUse]++;
	    typesOfDimer5p[dist5p][twoBases2index(refeBase,readBase)]++;
	    typesOfDimer3p[dist3p][twoBases2index(refeBase,readBase)]++;

	    if(isDeam5pS){
		typesOfDimer3pSingle[dist3p][twoBases2index(refeBase,readBase)]++;
	    }

	    if(isDeam3pS){
		typesOfDimer5pSingle[dist5p][twoBases2index(refeBase,readBase)]++;
	    }


	    if(isDeam5pD){
		typesOfDimer3pDouble[dist3p][twoBases2index(refeBase,readBase)]++;
	    }

	    if(isDeam3pD){
		typesOfDimer5pDouble[dist5p][twoBases2index(refeBase,readBase)]++;
	    }


	}
    }
}


double dbl2log(const double d,bool phred){
    double t= -10.0*(log(d)/log(10.0));
    // if(d == 0){
    // 	t = 
    // }
    if(phred)
	return t;
    else 
	return d;
}

int main (int argc, char *argv[]) {

    string file5p="/dev/stdout";
    string file3p="/dev/stdout";
    bool endo=false;

    bool allStr   =true;
    bool singleStr=false;
    bool doubleStr=false;
    bool singAnddoubleStr=false;

    int lengthMaxToPrint = 5;
    int minQualBase      = 0;

    bool dpFormat=false;
    bool hFormat=false;
    double errorToRemove=0.0;
    bool phred=false;

    string usage=string(""+string(argv[0])+" <options>  [in BAM file]"+
			"\nThis program reads a BAM file and produces a deamination profile for the\n"+
			"5' and 3' ends\n"+

			// "\nreads and the puts the rest into another bam file.\n"+
			// "\nTip: if you do not need one of them, use /dev/null as your output\n"+

			"\n\n\tOther options:\n"+
			"\t\t"+"-minq\t\t\tRequire the base to have at least this quality to be considered (Default: "+stringify( minQualBase )+")\n"+
			"\t\t"+"-endo\t\t\tRequire the 5' end to be deaminated to compute the 3' end and vice-versa (Default: "+stringify( endo )+")\n"+
			"\t\t"+"-length\t[length]\tDo not consider bases beyond this length  (Default: "+stringify(lengthMaxToPrint)+" ) \n"+
			"\t\t"+"-err\t[error rate]\tSubstract [error rate] from the rates to account for sequencing errors  (Default: "+stringify(errorToRemove)+" ) \n"+
			"\t\t"+"-log\t\t\tPrint substitutions on a PHRED logarithmic scale  (Default: "+stringify(phred)+" ) \n"+
			

			"\n\n\tYou can specify either one of the two:\n"+
			"\t\t"+"-single\t\t\tUse the deamination profile of a single strand library  (Default: "+booleanAsString( singleStr )+")\n"+
			"\t\t"+"-double\t\t\tUse the deamination profile of a double strand library  (Default: "+booleanAsString( doubleStr )+")\n"+
			"\n\tor specify this option:\n"+
			"\t\t"+"-both\t\t\tReport both C->T and G->A regardless of stand  (Default: "+booleanAsString( singAnddoubleStr )+")\n"+

			"\n\n\tOutput options:\n"+
			"\t\t"+"-5p\t[output file]\tOutput profile for the 5' end (Default: "+stringify(file5p)+")\n"+
			"\t\t"+"-3p\t[output file]\tOutput profile for the 3' end (Default: "+stringify(file3p)+")\n"+
			"\t\t"+"-dp\t\t\tOutput in damage-patterns format (Default: "+booleanAsString(dpFormat)+")\n"+
			"\t\t"+"-h\t\t\tMore human readible output (Default: "+booleanAsString(hFormat)+")\n"+
		       
			"\n");

    if(argc == 1 ||
       (argc == 2 && (string(argv[0]) == "-h" || string(argv[0]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }

    for(int i=1;i<(argc-1);i++){ //all but the last 3 args


        if(string(argv[i]) == "-dp"  ){
            dpFormat=true;
            continue;
        }

        if(string(argv[i]) == "-log"  ){
            phred=true;
            continue;
        }

        if(string(argv[i]) == "-h"  ){
            hFormat=true;
            continue;
        }

        if(string(argv[i]) == "-minq"  ){
            minQualBase=destringify<int>(argv[i+1]);
            i++;
            continue;
        }


        if(string(argv[i]) == "-length"  ){
            lengthMaxToPrint=destringify<int>(argv[i+1]);
            i++;
            continue;
        }

        if(string(argv[i]) == "-err"  ){
            errorToRemove=destringify<double>(argv[i+1]);
            i++;
            continue;
        }

        if(string(argv[i]) == "-5p" ){
	    file5p = string(argv[i+1]);
	    i++;
            continue;
        }

        if(string(argv[i]) == "-3p" ){
	    file3p = string(argv[i+1]);
	    i++;
            continue;
        }

        if(string(argv[i]) == "-endo" ){
	    endo   = true;
            continue;
        }

        if(string(argv[i]) == "-both" ){
	    //doubleStr=true;

	    allStr           = false;
	    singleStr        = false;
	    doubleStr        = false;
	    singAnddoubleStr = true;
            continue;
        }


        if(string(argv[i]) == "-single" ){

	    allStr    = false;
	    singleStr = true;
	    doubleStr = false;

            continue;
        }

        if(string(argv[i]) == "-double" ){
	    //doubleStr=true;

	    allStr    = false;
	    singleStr = false;
	    doubleStr = true;

            continue;
        }


	cerr<<"Error: unknown option "<<string(argv[i])<<endl;
	return 1;
    }

    
    if(phred && hFormat){
	cerr<<"Error: cannot specify both -log and -h"<<endl;
	return 1;
    }

    if(dpFormat && hFormat){
	cerr<<"Error: cannot specify both -dp and -h"<<endl;
	return 1;
    }

    if(endo){
	if(singAnddoubleStr){
	    cerr<<"Error: cannot use -singAnddoubleStr with -endo"<<endl;
	    return 1;
	}
	
	if( !singleStr &&
	    !doubleStr ){
	    cerr<<"Error: you have to provide the type of protocol used (single or double) when using endogenous"<<endl;
	    return 1;
	}

    }

    typesOfDimer5p       = vector< vector<unsigned int> >();
    typesOfDimer3p       = vector< vector<unsigned int> >();
    typesOfDimer5pDouble = vector< vector<unsigned int> >();
    typesOfDimer3pDouble = vector< vector<unsigned int> >();
    typesOfDimer5pSingle = vector< vector<unsigned int> >();
    typesOfDimer3pSingle = vector< vector<unsigned int> >();

    for(int l=0;l<MAXLENGTH;l++){
	//for(int i=0;i<16;i++){
	typesOfDimer5p.push_back( vector<unsigned int> ( 16,0 ) );
	typesOfDimer3p.push_back( vector<unsigned int> ( 16,0 ) );
	typesOfDimer5pDouble.push_back( vector<unsigned int> ( 16,0 ) );
	typesOfDimer3pDouble.push_back( vector<unsigned int> ( 16,0 ) );
	typesOfDimer5pSingle.push_back( vector<unsigned int> ( 16,0 ) );
	typesOfDimer3pSingle.push_back( vector<unsigned int> ( 16,0 ) );

	//}
    }

    string bamfiletopen = string( argv[ argc-1 ] );
    // string deambam      = string( argv[ argc-2 ] );
    // string nondeambam   = string( argv[ argc-1 ] );


    BamReader reader;
    
    if ( !reader.Open(bamfiletopen) ) {
    	cerr << "Could not open input BAM file"<< bamfiletopen << endl;
    	return 1;
    }








    //iterating over the alignments for these regions
    BamAlignment al;
    // bool pairedEnd=false;
    // bool firstRead=true;

    while ( reader.GetNextAlignment(al) ) {

	//cout<<"Read "<<al.Name<<" is wrong, cannot have a mixture of paired and unpaired read for this program"<<endl;
	if( al.IsPaired()  ){
	    // cerr<<"Read "<<al.Name<<" is wrong, cannot have a mixture of paired and unpaired read for this program"<<endl;
	    // return 1;
	    continue;
	}

	//skip unmapped
	if(!al.IsMapped()){
	    continue;
	}



	string reconstructedReference = reconstructRef(&al);
	if(al.Qualities.size() != reconstructedReference.size()){
	    cerr<<"Quality line is not the same size as the reconstructed reference"<<endl;
	    return 1;
	}


	
	increaseCounters(al,reconstructedReference,minQualBase); //start cycle numberOfCycles-1

	

      
    }//end while  each read
	


    reader.Close();

    ofstream file5pFP;
    file5pFP.open(file5p.c_str());

    if (!file5pFP.is_open()){
	cerr << "Unable to write to 5p file "<<file5p<<endl;
	exit(1);
    }




    //cout<<"cycle\tmatches\tmismatches\tmismatches%\tA>C\tA>C%\tA>G\tA>G%\tA>T\tA>T%\tC>A\tC>A%\tC>G\tC>G%\tC>T\tC>T%\tG>A\tG>A%\tG>C\tG>C%\tG>T\tG>T%\tT>A\tT>A%\tT>C\tT>C%\tT>G\tT>G%"<<endl;
    if(dpFormat)
	file5pFP<<"\t";
    if(hFormat)
	file5pFP<<"pos\t";
    file5pFP<<"A>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G"<<endl;
  

    vector< vector<unsigned int> > * typesOfDimer5pToUse;

    if(endo){
	if(doubleStr)
	    typesOfDimer5pToUse = &typesOfDimer5pDouble;
	else
	    typesOfDimer5pToUse = &typesOfDimer5pSingle;
    }else{
	typesOfDimer5pToUse     = &typesOfDimer5p;
    }

    for(int l=0;l<lengthMaxToPrint;l++){
	if(dpFormat)
	    file5pFP<<l<<"\t";

	if(hFormat)
	    file5pFP<<printIntAsWhitePaddedString(l,int(log10(lengthMaxToPrint))+1)<<"\t";
	
	for(int n1=0;n1<4;n1++){   
	    int totalObs=0;
	    for(int n2=0;n2<4;n2++){   
		totalObs+=(*typesOfDimer5pToUse)[l][4*n1+n2];
	    }

	    for(int n2=0;n2<4;n2++){   
		if(n1==n2)
		    continue;
		if(allStr){
		    if(dpFormat)
			file5pFP<<dbl2log(MAX(0.0,double( (*typesOfDimer5pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred)<<" [0..0]";
		    else
			if(hFormat)
			    file5pFP<<printDoubleAsWhitePaddedString( MAX(0.0,double( (*typesOfDimer5pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove) ,1,5);
			else
			    file5pFP<<dbl2log( MAX(0.0,double( (*typesOfDimer5pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred);  
		}else{ 
		    if(singAnddoubleStr){
			if(         (n1==1 && n2==3) || (n1==2 && n2==0 )  ) { 
			    if(dpFormat)
				file5pFP<<dbl2log( MAX(0.0,double((*typesOfDimer5pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred)<<" [0..0]";
			    else
				if(hFormat)
				    file5pFP<<printDoubleAsWhitePaddedString( MAX(0.0,double((*typesOfDimer5pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),1,5);
				else
				    file5pFP<<dbl2log( MAX(0.0,double((*typesOfDimer5pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred); 
			} else { 
			    if(dpFormat)
				file5pFP<<(phred?"-Inf":"0.0")<<" [0..0]";
			    else
				if(hFormat)
				    file5pFP<<printDoubleAsWhitePaddedString( 0.0,1,5);
				else
				    file5pFP<<(phred?"-Inf":"0.0")<<"";
			}

		    }else{
			if(doubleStr){
			    //          C        T
			    if(         n1==1 && n2==3  ) { 
				if(dpFormat)
				    file5pFP<<dbl2log( MAX(0.0,double((*typesOfDimer5pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred)<<" [0..0]";
				else
				    if(hFormat)
					file5pFP<<printDoubleAsWhitePaddedString( MAX(0.0,double((*typesOfDimer5pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove) ,1,5); 
				    else
					file5pFP<<dbl2log( MAX(0.0,double((*typesOfDimer5pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred); 
			    } else { 
				if(dpFormat)
				    file5pFP<<(phred?"-Inf":"0.0")<<" [0..0]";
				else
				    if(hFormat)
					file5pFP<<printDoubleAsWhitePaddedString( 0.0,1,5);
				    else					
					file5pFP<<(phred?"-Inf":"0.0"); 
			    }
			}else{ 
			    if(singleStr){
				//      C        T
				if(     n1==1 && n2==3  ) { 
				    if(dpFormat)
					file5pFP<<dbl2log( MAX(0.0,double((*typesOfDimer5pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred)<<" [0..0]"; 
				    else
					if(hFormat)
					    file5pFP<<printDoubleAsWhitePaddedString(MAX(0.0,double((*typesOfDimer5pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),1,5);
					else
					    file5pFP<<dbl2log( MAX(0.0,double((*typesOfDimer5pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred); 
					    
				} else { 
				    if(dpFormat)
					file5pFP<<(phred?"-Inf":"0.0")<<" [0..0]";
				    else
					if(hFormat)
					    file5pFP<<printDoubleAsWhitePaddedString(0.0,1,5);
					else										    
					    file5pFP<<(phred?"-Inf":"0.0"); 
				}
			    }
			}
		    }
		}

		
		if(!(n1 ==3 && n2 == 2 ))
		    file5pFP<<"\t";
	    }


	}
	file5pFP<<endl;
    }


    file5pFP.close();

    ofstream file3pFP;
    if(file3p == "/dev/stdout"){
	file3pFP.open(file3p.c_str(), ofstream::out | ofstream::app);
    }else{
	file3pFP.open(file3p.c_str());
    }

    if (!file3pFP.is_open()){
	cerr << "Unable to write to 3p file "<<file3p<<endl;
	exit(1);
    }

    if(dpFormat)
	file3pFP<<"\t";
    if(hFormat)
	file3pFP<<"pos\t";

    file3pFP<<"A>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G"<<endl;


    vector< vector<unsigned int> > * typesOfDimer3pToUse;

    if(endo){
	if(doubleStr)
	    typesOfDimer3pToUse = &typesOfDimer3pDouble;
	else
	    typesOfDimer3pToUse = &typesOfDimer3pSingle;
    }else{
	typesOfDimer3pToUse = &typesOfDimer3p;
    }



    for(int le=0;le<lengthMaxToPrint;le++){

	int l=le;
	if(dpFormat){
	    l=lengthMaxToPrint-1-le;
	    if(l==0)
		file3pFP<<""<<l<<"\t";	    
	    else
		file3pFP<<"-"<<l<<"\t";	    
	}

	if(hFormat){
	    //l=lengthMaxToPrint-1-le;
	    // if(l==0)
	    // 	file3pFP<<""<<printIntAsWhitePaddedString(l,int(log10(lengthMaxToPrint)))<<"\t";	    
	    // else
	    // 	file3pFP<<"-"<<printIntAsWhitePaddedString(l,int(log10(lengthMaxToPrint)))<<"\t";	    
	    file3pFP<<""<<printIntAsWhitePaddedString(l,int(log10(lengthMaxToPrint))+1)<<"\t";	    
	}

	for(int n1=0;n1<4;n1++){   
	    int totalObs=0;
	    for(int n2=0;n2<4;n2++){   
		totalObs+=(*typesOfDimer3pToUse)[l][4*n1+n2];
	    }

	    for(int n2=0;n2<4;n2++){   
		if(n1==n2)
		    continue;
		if(allStr){
		    if(dpFormat)
			file3pFP<<dbl2log( MAX(0.0,double( (*typesOfDimer3pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred)<<" [0..0]";
		    else
			if(hFormat)
			    file3pFP<<printDoubleAsWhitePaddedString( MAX(0.0,double( (*typesOfDimer3pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove) ,1,5);
			else
			    file3pFP<<dbl2log( MAX(0.0,double( (*typesOfDimer3pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred);
		}else{ 
		    if(singAnddoubleStr){			
			if(   (n1==1 && n2==3) || (n1==2 && n2==0 )  ) { 
			    if(dpFormat)
				file3pFP<<dbl2log( MAX(0.0,double((*typesOfDimer3pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred)<<" [0..0]"; 
			    else
				if(hFormat)
				    file3pFP<<printDoubleAsWhitePaddedString( MAX(0.0,double((*typesOfDimer3pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove) ,1,5);
				else
				    file3pFP<<dbl2log( MAX(0.0,double((*typesOfDimer3pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred); 
			} else { 
			    if(dpFormat)
				file3pFP<<(phred?"-Inf":"0.0")<<" [0..0]"; 
			    else
				if(hFormat)
				    file3pFP<<printDoubleAsWhitePaddedString( 0.0 ,1,5);
				else
				    file3pFP<<(phred?"-Inf":"0.0"); 
			}

		    }else{
			if(doubleStr){
			    //          G        A
			    if(         n1==2 && n2==0  ) { 
				if(dpFormat)
				    file3pFP<<dbl2log( MAX(0.0,double((*typesOfDimer3pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred)<<" [0..0]"; 
				else
				    if(hFormat)
					file3pFP<<printDoubleAsWhitePaddedString( MAX(0.0,double((*typesOfDimer3pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove) ,1,5);
				    else					
					file3pFP<<dbl2log( MAX(0.0,double((*typesOfDimer3pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred); 
			    } else { 				
				if(dpFormat)
				    file3pFP<<(phred?"-Inf":"0.0")<<" [0..0]"; 
				else
				    if(hFormat)
					file3pFP<<printDoubleAsWhitePaddedString(  0.0 ,1,5);
				    else	
					file3pFP<<(phred?"-Inf":"0.0"); 
			    }
			}else{ 
			    if(singleStr){
				//      C        T
				if(     n1==1 && n2==3  ) { 
				    if(dpFormat)
					file3pFP<<dbl2log( MAX(0.0,double((*typesOfDimer3pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred)<<" [0..0]"; 
				    else
					if(hFormat)
					    file3pFP<<printDoubleAsWhitePaddedString( MAX(0.0,double((*typesOfDimer3pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove) ,1,5);
					else					
					    file3pFP<<dbl2log( MAX(0.0,double((*typesOfDimer3pToUse)[l][4*n1+n2])/double(totalObs)-errorToRemove),phred); 
				} else { 				    
				    if(dpFormat)
					file3pFP<<(phred?"-Inf":"0.0")<<" [0..0]"; 
				    else
					if(hFormat)
					    file3pFP<<printDoubleAsWhitePaddedString(  0.0 , 1,5);
					else	
					    file3pFP<<(phred?"-Inf":"0.0"); 

				}
				
			    }
			}
		    }
		}

		
		if(!(n1 ==3 && n2 == 2 ))
		    file3pFP<<"\t";
	    }


	}
	file3pFP<<endl;
    }



    file3pFP.close();
   
    return 0;
}

