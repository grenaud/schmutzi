#include <iostream>
#include <vector>
#include <set>
#include <ctype.h>
#include <stdlib.h>

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"
#include "api/SamSequenceDictionary.h"

#include "utils.h"
#include "ReconsReferenceBAM.h"

using namespace std;
using namespace BamTools;

const int offset=33;


bool isTransition(char reference,char read){

    if(  (reference == 'A' && read == 'G') ||
	 (reference == 'G' && read == 'A') ||
	 (reference == 'C' && read == 'T') ||
	 (reference == 'T' && read == 'C') ){
	return true;
    }
	 
    return false;
}



bool isTransversions(char reference,char read){
    if(  (reference == 'A' && read == 'C') ||
	 (reference == 'A' && read == 'T') ||
	 
	 (reference == 'G' && read == 'C') ||
	 (reference == 'G' && read == 'T') ||
	 
	 (reference == 'C' && read == 'A') ||
	 (reference == 'C' && read == 'G') ||

	 (reference == 'T' && read == 'A') ||
	 (reference == 'T' && read == 'G') ){
	return true;
    }

    return false;	
}





int main (int argc, char *argv[]) {

    int  minBaseQuality = 0;

    string usage=string(""+string(argv[0])+"  [in BAM file]"+
			"\nThis program computes the joint freq of deamination\n"+
			// "\nreads and the puts the rest into another bam file.\n"+
			// "\nTip: if you do not need one of them, use /dev/null as your output\n"+
			"arguments:\n"+
			"\t"+"--bq  [base qual]   : Minimum base quality to flag a deaminated site (Default: "+stringify(minBaseQuality)+")\n"+
			"\n");

    if(argc == 1 ||
       argc < 2  ||
       (argc == 2 && (string(argv[0]) == "-h" || string(argv[0]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }


    for(int i=1;i<(argc-2);i++){ 

	
        if(string(argv[i]) == "--bq"){
	    minBaseQuality=destringify<int>(argv[i+1]);
            i++;
            continue;
	}

    }

    string bamfiletopen = string( argv[ argc-1 ] );


    BamReader reader;
    
    if ( !reader.Open(bamfiletopen) ) {
    	cerr << "Could not open input BAM file"<< bamfiletopen << endl;
    	return 1;
    }

    vector<RefData>  testRefData=reader.GetReferenceData();
    const SamHeader header = reader.GetHeader();
    const RefVector references = reader.GetReferenceData();

    // BamWriter writerDeam;
    // if ( !writerDeam.Open(deambam,      header, references) ) {
    // 	cerr << "Could not open output BAM file" << endl;
    // 	return 1;
    // }

    // BamWriter writerNoDeam;
    // if ( !writerNoDeam.Open(nondeambam, header, references) ) {
    // 	cerr << "Could not open output BAM file" << endl;
    // 	return 1;
    // }







    //iterating over the alignments for these regions
    BamAlignment al;
    int i;
    unsigned int d5pgiven3p  =0;
    unsigned int nod5pgiven3p=0;
    unsigned int d5not3p     =0;
    unsigned int nod5not3p   =0;
    unsigned int d5     =0;
    unsigned int nod5   =0;

    unsigned int d3pgiven5p  =0;
    unsigned int nod3pgiven5p=0;
    unsigned int d3not5p     =0;
    unsigned int nod3not5p   =0;
   unsigned int d3     =0;
    unsigned int nod3   =0;

    while ( reader.GetNextAlignment(al) ) {


	//skip unmapped
	if(!al.IsMapped())
	    continue;

	if(al.IsPaired() ){  
	    continue;
	    return 1;
	}


	string reconstructedReference = reconstructRef(&al);
	char refeBase;
	char readBase;
	bool isDeaminated5p=false;
	bool isDeaminated3p=false;
	bool count5p=false;
	bool count3p=false;

	if(al.Qualities.size() != reconstructedReference.size()){
	    cerr<<"Quality line is not the same size as the reconstructed reference"<<endl;
	    return 1;
	}



	if(al.IsReverseStrand()){
	    //continue;
	    //first base next to 3'
	    i = 0 ;
	    refeBase=toupper(reconstructedReference[i]);
	    readBase=toupper(         al.QueryBases[i]);
	    if(  readBase  == 'A' && refeBase  == 'G' && int(al.Qualities[i]-offset) >= minBaseQuality){  isDeaminated3p=true; count3p=true;}
	    if(  readBase  == 'G' && refeBase  == 'M' && int(al.Qualities[i]-offset) >= minBaseQuality){  isDeaminated3p=false;count3p=true; }


	    //last  base next to 5'
	    i = (al.QueryBases.length()-1) ;
	    refeBase=toupper(reconstructedReference[i]);
	    readBase=toupper(         al.QueryBases[i]);
	    if( readBase  == 'A' && refeBase  == 'G'  && int(al.Qualities[i]-offset) >= minBaseQuality){  isDeaminated5p=true;  count5p=true; }
	    if( readBase  == 'G' && refeBase  == 'M'  && int(al.Qualities[i]-offset) >= minBaseQuality){  isDeaminated5p=false; count5p=true; }


	}else{

		
	    //first base next to 5'
	    i = 0;
	    refeBase=toupper(reconstructedReference[i]);
	    readBase=toupper(         al.QueryBases[i]);
	    if( readBase  == 'T' && refeBase  == 'C'  && int(al.Qualities[i]-offset) >= minBaseQuality){  isDeaminated5p=true;  count5p=true; }
	    if( readBase  == 'C' && refeBase  == 'M'  && int(al.Qualities[i]-offset) >= minBaseQuality){  isDeaminated5p=false; count5p=true; }


	    //last base next to 3'
	    i = (al.QueryBases.length()-1);
	    refeBase=toupper(reconstructedReference[i]);
	    readBase=toupper(         al.QueryBases[i]);
	    if( readBase  == 'T' && refeBase  == 'C'  && int(al.Qualities[i]-offset) >= minBaseQuality){  isDeaminated3p=true;  count3p=true; }	
	    if( readBase  == 'C' && refeBase  == 'M'  && int(al.Qualities[i]-offset) >= minBaseQuality){  isDeaminated3p=false; count3p=true; }	
	    
	}
		  
	// if(!(count5p && count3p))
	//     continue;
	if( count5p ){

	    if(  isDeaminated5p)
		d5++;
	    if( !isDeaminated5p)
		nod5++;
	}
	if( count3p ){

	    if(  isDeaminated3p)
		d3++;
	    if( !isDeaminated3p)
		nod3++;
	}

	if( count3p && 	count5p ){

	    if(  isDeaminated3p) {

		if(  isDeaminated5p)
		    d5pgiven3p++;
		if( !isDeaminated5p)
		    nod5pgiven3p++;

	    }else{

		if(  isDeaminated5p)
		    d5not3p++;
		if( !isDeaminated5p)
		    nod5not3p++;

	    }
	// }

	// if(count5p){
	    if( isDeaminated5p) {

		if(  isDeaminated3p)
		    d3pgiven5p++;
		if( !isDeaminated3p)
		    nod3pgiven5p++;

	    }else{

		if(  isDeaminated3p)
		    d3not5p++;
		if( !isDeaminated3p)
		    nod3not5p++;

	    }
	}
	// if(   isDeaminated5p &&  isDeaminated3p ) b53d++;
	// if(   isDeaminated5p &&  isDeaminated3p ) o5d++;
	// 	if(  !isDeaminated5p &&  isDeaminated3p ) o3d++;
	// 	if(  !isDeaminated5p && !isDeaminated3p ) nod++;


    
    }//end for each read

    //    cout<<b53d<<"\t"<<o5d<<"\t"<<o3d<<"\t"<<nod<<endl;
    cout<<d5pgiven3p<<"\t"<<nod5pgiven3p<<"\t"<<d5not3p<<"\t"<<nod5not3p<<"\t"<<d5<<"\t"<<nod5<<"\t"<<double(d5pgiven3p)/double( d5pgiven3p+nod5pgiven3p)<<"\t"<<double(d5not3p)/double( d5not3p+nod5not3p)<<"\t"<<double(d5)/double( d5+nod5)<<endl;
    
    cout<<d3pgiven5p<<"\t"<<nod3pgiven5p<<"\t"<<d3not5p<<"\t"<<nod3not5p<<"\t"<<d3<<"\t"<<nod3<<"\t"<<double(d3pgiven5p)/double( d3pgiven5p+nod3pgiven5p)<<"\t"<<double(d3not5p)/double( d3not5p+nod3not5p)<<"\t"<<double(d3)/double( d3+nod3)<<endl;

    






    reader.Close();



   
    return 0;
}

