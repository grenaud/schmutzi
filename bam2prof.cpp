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






vector<unsigned int>    matches;
vector<unsigned int> mismatches;
vector< vector<unsigned int> > typesOfMismatches;


//increases the counters mismatches and typesOfMismatches of a given BamAlignment object
inline void increaseCounters(BamAlignment & al,string & reconstructedReference,int firstCycleRead,int increment){

    char refeBase;
    char readBase;
    int cycleToUse=firstCycleRead;
    // cout<<"name "<<al.Name<<endl;
    // cout<<"firstCycleRead "<<firstCycleRead<<endl;
    // cout<<"increment      "<<increment<<endl;

    for(int i=0;i<numberOfCycles;i++,cycleToUse+=increment){
        // cout<<"i = "<<i<<" cyc "<<cycleToUse<<endl;
	refeBase=toupper(reconstructedReference[i]);
	readBase=toupper(         al.QueryBases[i]);
		     
	//match
	if(refeBase == 'M'){
	    matches[cycleToUse]++;
	    continue;
	}

	if(refeBase == 'S' ||refeBase == 'I'){ //don't care about soft clipped or indels
	    continue;
	}
		    
	//mismatch
	if( isResolvedDNA(refeBase)  && 
	    isResolvedDNA(readBase) ){
	    if(al.IsReverseStrand()){ //need to take the complement
		refeBase=complement(refeBase);
		readBase=complement(readBase);
	    }
	    if(readBase == refeBase){
		cerr<<"Internal error in reconstruction of read "<<al.Name<<", contact developer"<<endl;
		exit(1);;
	    }						
	    mismatches[cycleToUse]++;
	    typesOfMismatches[dimer2index(refeBase,readBase)][cycleToUse]++;
	    continue;
	}		    		     
    }
}


int main (int argc, char *argv[]) {

    string file5p="/dev/stdout";
    string file3p="/dev/stdout";

    string usage=string(""+string(argv[0])+"  [in BAM file]"+
			"\nThis program reads a BAM file and produces a deamination profile for the\n"+
			"5' and 3' ends\n"+

			// "\nreads and the puts the rest into another bam file.\n"+
			// "\nTip: if you do not need one of them, use /dev/null as your output\n"+
			"\n\n\tOutput options:\n"+
			 "\t\t"+"-5p\t[output file]\tOutput profile for the 5' end (Default: "+stringify(file5p)+")\n"+
			 "\t\t"+"-3p\t[output file]\tOutput profile for the 3' end (Default: "+stringify(file3p)+")\n"+
		       
			"\n");

    if(argc == 1 ||
       (argc == 2 && (string(argv[0]) == "-h" || string(argv[0]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }

    for(int i=1;i<(argc-1);i++){ //all but the last 3 args

        if(strcmp(argv[i],"-5p") == 0 ){
            continue;
        }
        if(strcmp(argv[i],"-5p") == 0 ){
            continue;
        }



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
    bool pairedEnd=false;
    bool firstRead=true;

    while ( reader.GetNextAlignment(al) ) {

	if(firstRead){ //reads are either all paired end or single end, I don't allow a mix
	    numberOfCycles=al.QueryBases.size();
	    // cout<<"numberOfCycles "<<numberOfCycles<<endl;
	    if(al.IsPaired() ){  
		pairedEnd=true;		
		matches    = vector<unsigned int> (2*numberOfCycles,0);
		mismatches = vector<unsigned int> (2*numberOfCycles,0);
		typesOfMismatches = vector< vector<unsigned int> >();
		for(int i=0;i<12;i++)
		    typesOfMismatches.push_back( vector<unsigned int> (2*numberOfCycles,0) );
	    }else{
		matches    = vector<unsigned int> (  numberOfCycles,0);
		mismatches = vector<unsigned int> (  numberOfCycles,0);
		typesOfMismatches = vector< vector<unsigned int> >();
		for(int i=0;i<12;i++)
		    typesOfMismatches.push_back( vector<unsigned int> (  numberOfCycles,0) );
	    }
	    firstRead=false;
	}

	if( (  pairedEnd && !al.IsPaired()) ||  
	    ( !pairedEnd &&  al.IsPaired())   ){
	    cerr<<"Read "<<al.Name<<" is wrong, cannot have a mixture of paired and unpaired read for this program"<<endl;
	    return 1;
	}


	//skip unmapped
	if(!al.IsMapped())
	    continue;

	if(numberOfCycles!=int(al.QueryBases.size())){
	    cerr<<"The length of read "<<al.Name<<" is wrong, should be "<<numberOfCycles<<"bp"<<endl;
	    return 1;
	}


	string reconstructedReference = reconstructRef(&al);
	if(al.Qualities.size() != reconstructedReference.size()){
	    cerr<<"Quality line is not the same size as the reconstructed reference"<<endl;
	    return 1;
	}


	if( pairedEnd ){
	    if( al.IsFirstMate() ){ //start cycle 0

		if( al.IsReverseStrand()  ){ 
		    increaseCounters(al,reconstructedReference,numberOfCycles-1,-1); //start cycle numberOfCycles-1
		}else{
		    increaseCounters(al,reconstructedReference,0               , 1); //start cycle 0
		}
       
	    }else{

		if( al.IsSecondMate()  ){ 
		    if( al.IsReverseStrand()  ){ 
			increaseCounters(al,reconstructedReference,2*numberOfCycles-1,-1); //start cycle 2*numberOfCycles-1
		    }else{
			increaseCounters(al,reconstructedReference,numberOfCycles    , 1); //start cycle numberOfCycles
		    }
		}else{
	    	    cerr<<"Reads "<<al.Name<<" must be either first or second mate"<<endl;
	    	    return 1;
	    	}
	    }	
	}else{ //single end 

	    if( al.IsReverseStrand()  ){ 
		increaseCounters(al,reconstructedReference,numberOfCycles-1,-1); //start cycle numberOfCycles-1
	    }else{
		increaseCounters(al,reconstructedReference,0               , 1); //start cycle 0
	    }
	    
	}


    }//end while  each read
	


    reader.Close();

    cout<<"cycle\tmatches\tmismatches\tmismatches%\tA>C\tA>C%\tA>G\tA>G%\tA>T\tA>T%\tC>A\tC>A%\tC>G\tC>G%\tC>T\tC>T%\tG>A\tG>A%\tG>C\tG>C%\tG>T\tG>T%\tT>A\tT>A%\tT>C\tT>C%\tT>G\tT>G%"<<endl;



    for(unsigned int i=0;i<matches.size();i++){
cout<<(i+1);
	if( (matches[i]+mismatches[i]!=0) )
	    cout<<"\t"<<matches[i]<<"\t"<<mismatches[i]<<"\t"<< 100.0*(double(mismatches[i])/double(matches[i]+mismatches[i])) ;
	else
	    cout<<"\t"<<matches[i]<<"\t"<<mismatches[i]<<"\tNA";

	for(int j=0;j<12;j++){   
	    cout<<"\t"<<typesOfMismatches[j][i];
	    if( (matches[i]+mismatches[i]!=0) )
		cout<<"\t"<<100.0*double(typesOfMismatches[j][i])/double(matches[i]+mismatches[i]);
	    else
		cout<<"\tNA";
	}

	cout<<endl;
    }




   
    return 0;
}

