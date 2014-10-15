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
			toaddSub.s[indexFirstArray++]           = 1.0;		    
		    }else{ //           rely on the substitution frequency
			sumMismatchProb                         += tempFreq.s[indexSecondArray];
			toaddSub.s[indexFirstArray++]            = tempFreq.s[indexSecondArray++];
		    }
		}

		toaddSub.s[indexInArrayMatch] = 1.0 - sumMismatchProb;
	    }

	    // for(int nuc1=0;nuc1<4;nuc1++){
	    // 	for(int nuc2=0;nuc2<4;nuc2++){
	    // 	    cout<<(nuc1*4+nuc2)<<"\t"<<toaddSub.s[nuc1*4+nuc2]<<endl;
	    // 	}	       
	    // }
	    
	    // exit(1);

	    subVec.push_back( toaddSub );
	}	             	              
	subFP.close();
    }else{
	cerr << "Unable to open file "<<filename<<endl;
	exit(1);
    }



}


void readMTConsensus(const string consensusFile,
		     map<int, PHREDgeno> & pos2phredgeno,
		     int & sizeGenome,
		     vector<int> & posOfIndels){

    string line;
    igzstream consensusFD;
    consensusFD.open(consensusFile.c_str());
    if (consensusFD.good()){
	getline (consensusFD,line);

	while ( getline (consensusFD,line)){
	    if (line.empty())
		continue;

	    vector<string> fields = allTokens(line,'\t');
	    PHREDgeno toadd;
	    // cerr<<line<<endl;


	    if(fields.size() != 11){
		cerr << "line "<<line<<"  in file  "<<consensusFile<<" does not have 11 fields"<<endl;
		exit(1);
	    }
	    

	    if(fields[0][fields[0].size()-1] == 'i'){ //skip insertion
		posOfIndels.push_back( destringify<int>( fields[0]) );
		continue;
	    }

	    if(fields[2] == "D"){ //skip deletions
		posOfIndels.push_back( destringify<int>( fields[0]) );
		continue;
	    }	    

	    toadd.consensus = fields[2][0];
	    for(int nuc=0;nuc<4;nuc++){		
		toadd.phred[nuc]  = destringify<double>(fields[nuc+7]);		
		toadd.perror[nuc] = pow(10.0,toadd.phred[nuc]/(-10.0));		
	    }

	    pos2phredgeno[     destringify<int>( fields[0])   ] = toadd;
	    sizeGenome =  max( destringify<int>( fields[0]), sizeGenome);
	    // cout<<destringify<int>( fields[0])<<endl;
	    
	}
	consensusFD.close();

    }else{
	cerr << "Cannot open consensus file  "<<consensusFile<<""<<endl;
	exit(1);
    }


}

void readMTAlleleFreq(const string freqFile,	map<int, alleleFrequency> & pos2allelefreq){
    // map<int, alleleFrequency> pos2allelefreq;

    string line;
    igzstream freqAlleleFile;
    freqAlleleFile.open(freqFile.c_str());
    if (freqAlleleFile.good()){

	while ( getline (freqAlleleFile,line)){

	    vector<string> fields = allTokens(line,'\t');
	    alleleFrequency freqToadd;
	    
	    if(fields.size() != 5){
		cerr << "line "<<line<<"  in file  "<<freqFile<<" does not have 5 fields"<<endl;
		exit(1);
	    }
	   

	    for(int nuc=0;nuc<4;nuc++){
		freqToadd.f[nuc]=destringify<double>(fields[nuc+1]);
	    }

	    pos2allelefreq[ destringify<int>( fields[0])  ] = freqToadd;
	    	    
	}
	freqAlleleFile.close();

    }else{
	cerr << "Cannot open allele frequency file  "<<freqFile<<""<<endl;
	exit(1);
    }

}
