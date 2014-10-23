#include <gzstream.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>


#include "miscfunc.h"
#include "utils.h"


using namespace std;




int main (int argc, char *argv[]) {

    double qcutoff = 0;
    const string usage=string("\t"+string(argv[0])+
			      " This program takes two log files and prints the positions where they differ"+
			      " [options]  [log file endo]  [log file cont]  "+"\n\n"+			      
			      "\n\tOptions:\n"+	
			      "\t\t"+"-q  [quality cutoff]" +"\t\t"+"Only consider sites with quality greater than this cutoff (default: 0)"+"\n"+
			      // "\t\t"+"-err  [error rate]" +"\t\t"+"Output log (default "+stringify(errorRate)+")"+
			      
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

    for(int i=1;i<(argc-2);i++){ //all but the last arg

	if(string(argv[i]) == "-q" ){
	    qcutoff = destringify<double>(string(argv[i+1]));
	    i++;
	    continue;
	}
       
	cerr<<"Unknown option "<<argv[i] <<" exiting"<<endl;
	return 1;
    }


    string logFileE     = string(argv[argc-2]);
    string logFileC     = string(argv[argc-1]);
     


    string lineE;
    string lineC;

    if(isDos(logFileE) || isMac(logFileE) ){
	cerr << "File  "<<logFileE<<" must be unix formatted, exiting"<<endl;
	return 1;
    }

    if(isDos(logFileC) || isMac(logFileC) ){
	cerr << "File  "<<logFileC<<" must be unix formatted, exiting"<<endl;
	return 1;
    }


    igzstream logFDE;
    igzstream logFDC;

    logFDE.open(logFileE.c_str(),ios::in);
    logFDC.open(logFileC.c_str(),ios::in);

    if (!logFDE.good()){
	cerr << "Cannot open log file  "<<logFileE<<""<<endl;
	return 1;
    }

    if (!logFDC.good()){
	cerr << "Cannot open log file  "<<logFileC<<""<<endl;
	return 1;
    }

    getline (logFDE,lineE); //header

    vector<logRecord> logRecordsE;
    vector<logRecord> logRecordsC;

    while ( getline (logFDE,lineE)){
	if(lineE.empty())
	    continue;
	vector<string> fields= allTokens(lineE,'\t');

	if(strEndsWith(fields[0],"i") ){ //skip insertions
	    continue;
	}

	logRecord toadd;

	toadd.pos     =    destringify<int>( fields[0]);
	toadd.ref     =   destringify<char>( fields[1]);
	toadd.base    =   destringify<char>( fields[2]);
	  
	toadd.q       = destringify<double>( fields[4] );
	toadd.aprob   = destringify<double>( fields[7] );
	toadd.cprob   = destringify<double>( fields[8] );
	toadd.gprob   = destringify<double>( fields[9] );
	toadd.tprob   = destringify<double>( fields[10] );

	logRecordsE.push_back(toadd);
    }
    logFDE.close();

    getline (logFDC,lineC); //header

    while ( getline (logFDC,lineC)){
	if(lineC.empty())
	    continue;
	vector<string> fields= allTokens(lineC,'\t');

	if(strEndsWith(fields[0],"i") ){ //skip insertions
	    continue;
	}

	logRecord toadd;

	toadd.pos     =    destringify<int>( fields[0]);
	toadd.ref     =   destringify<char>( fields[1]);
	toadd.base    =   destringify<char>( fields[2]);
	  
	toadd.q       = destringify<double>( fields[3] );
	// toadd.aprob   = destringify<double>( fields[7] );
	// toadd.cprob   = destringify<double>( fields[8] );
	// toadd.gprob   = destringify<double>( fields[9] );
	// toadd.tprob   = destringify<double>( fields[10] );

	logRecordsC.push_back(toadd);
    }
    logFDC.close();



    for(unsigned int pos=0;pos<logRecordsE.size();pos++){
	if(logRecordsE[pos].pos != logRecordsC[pos].pos ){
	    cerr<<"Positions between the endogenous and contaminant file differs at pos="<<pos<<"\t"<<logRecordsE[pos].pos<<"\t"<<logRecordsC[pos].pos<<endl;
	    return 1;
	}
	if(logRecordsE[pos].q<qcutoff)
	    continue;
	if(logRecordsC[pos].q<qcutoff)
	    continue;

	if(logRecordsE[pos].base ==  logRecordsC[pos].base )
	    continue;

	if(logRecordsE[pos].base ==  'D' ||
	   logRecordsC[pos].base ==  'D' )
	    continue;

	cout<<logRecordsC[pos].pos<<"\t"<<logRecordsC[pos].base<<"\tcont"<<endl;
	cout<<logRecordsE[pos].pos<<"\t"<<logRecordsE[pos].base<<"\tendo"<<endl;

	
    }


    cerr<<"program "<<argv[0]<<" finished successfully"<<endl;

    return 0;
}

