/*
 * miscfunc
 * Date: Jun-08-2014 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */


#include "miscfunc.h"


void readIlluminaError(const string errFile,probSubstition & illuminaErrorsProb){

    igzstream errFileSt;

    errFileSt.open(errFile.c_str(), ios::in);

    //    unsigned int counterCont=0;
    if (errFileSt.good()){
	vector<string> fields;
	string line;
	//header
	if ( !getline (errFileSt,line)){
	    cerr << "Unable to open file "<<errFile<<endl;
	    exit(1);
	}
	fields = allTokens(line,'\t');
	
	if(fields.size() != 12){
	    cerr << "line from error profile does not have 12 fields "<<line<<endl;
	    exit(1);
	}

	//raw sums
	if ( !getline (errFileSt,line)){
	    cerr << "Unable to open file "<<errFile<<endl;
	    exit(1);
	}
	
	fields = allTokens(line,'\t');

	if(fields.size() != 12){
	    cerr << "line from error profile does not have 12 fields "<<line<<endl;
	    exit(1);
	}

	//probs
	if ( !getline (errFileSt,line)){
	    cerr << "Unable to open file "<<errFile<<endl;
	    exit(1);
	}
	
	fields = allTokens(line,'\t');

	if(fields.size() != 12){
	    cerr << "line from error profile does not have 12 fields "<<line<<endl;
	    exit(1);
	}
	substitutionRates tempFreq;	
	    
	
	for(unsigned int k=0;k<=9;k+=3){	

	    for(unsigned int t=0;t<=2;t++){	
		tempFreq.s[k+t]=destringify<double>(fields[k+t]);
		//cerr<<freqIlluminaError.s[k+t]<<endl;
	    }

	}


	int indexFirstArray =0;
	int indexSecondArray=0;

	for(int nuc1=0;nuc1<4;nuc1++){
	    for(int nuc2=0;nuc2<4;nuc2++){
		if(nuc1==nuc2) // prob of error is 0 if both nucleotides are identical
		    illuminaErrorsProb.s[indexFirstArray++]=0.0;
		else //           rely on the substitution frequency
		    illuminaErrorsProb.s[indexFirstArray++]=tempFreq.s[indexSecondArray++];
	    }
	}
	
	             	              
	errFileSt.close();
    }else{
	cerr << "Unable to open file "<<errFile<<endl;
	exit(1);
    }



}

void readNucSubstitionFreq(const string filename,vector<probSubstition> & subVec){
    igzstream subFP;

    subFP.open(filename.c_str(), ios::in);

    //    unsigned int counterCont=0;
    if (subFP.good()){
	vector<string> fields;
	string line;

	//header
	if ( !getline (subFP,line)){
	    cerr << "Unable to open file "<<filename<<endl;
	    exit(1);
	}
	fields = allTokens(line,'\t');
	
	if(fields.size() != 12){
	    cerr << "line from error profile does not have 12 fields "<<line<<endl;
	    exit(1);
	}


	//probs
	while ( getline (subFP,line)){
	    
	    fields = allTokens(line,'\t');

	    if(fields.size() != 12){
		cerr << "line from error profile does not have 12 fields "<<line<<endl;
		exit(1);
	    }

	    substitutionRates tempFreq;	
	    probSubstition toaddSub;


	    for(unsigned int k=0;k<=9;k+=3){	

		for(unsigned int t=0;t<=2;t++){	
		    tempFreq.s[k+t]=destringify<double>(fields[k+t]);
		}

	    }


	    int indexFirstArray =0;
	    int indexSecondArray=0;

	    for(int nuc1=0;nuc1<4;nuc1++){
		double sumMismatchProb=0.0;
		int indexInArrayMatch=1;
		for(int nuc2=0;nuc2<4;nuc2++){
		    if(nuc1==nuc2){ // prob of error is 0 if both nucleotides are identical
			indexInArrayMatch                       = indexFirstArray;
			toaddSub.s[indexFirstArray++] = 1.0;		    
		    }else{ //           rely on the substitution frequency
			sumMismatchProb                         += tempFreq.s[indexSecondArray];
			toaddSub.s[indexFirstArray++]            = tempFreq.s[indexSecondArray++];
		    }
		}

		toaddSub.s[indexInArrayMatch] = 1.0 - sumMismatchProb;
	    }
	    
	    subVec.push_back( toaddSub );
	}	             	              
	subFP.close();
    }else{
	cerr << "Unable to open file "<<filename<<endl;
	exit(1);
    }



}
