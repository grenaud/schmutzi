#include <gzstream.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>


#include "utils.h"


using namespace std;




int main (int argc, char *argv[]) {

    double qcutoff = 0;
    const string usage=string("\t"+string(argv[0])+
			      " This program transform a log file into a contamination profile "+
			      " [options]  [log file]  "+"\n\n"+			      
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

    for(int i=1;i<(argc-1);i++){ //all but the last arg

	if(string(argv[i]) == "-q" ){
	    qcutoff = destringify<double>(string(argv[i+1]));
	    i++;
	    continue;
	}
       
	cerr<<"Unknown option "<<argv[i] <<" exiting"<<endl;
	return 1;
    }


    string logFile     = string(argv[argc-1]);
     


    string line;
    if(isDos(logFile) || isMac(logFile) ){
	cerr << "File  "<<logFile<<" must be unix formatted, exiting"<<endl;
	return 1;
    }


    igzstream logFD;

    logFD.open(logFile.c_str(),ios::in);
    if (logFD.good()){
	getline (logFD,line); //header
	while ( getline (logFD,line)){
	    if(line.empty())
		continue;
	    vector<string> fields= allTokens(line,'\t');
	    string  pos    =                      fields[0];
	    string  alt    =                      fields[2];

	    if(strEndsWith(pos,"i") || //skip indel
	       alt  == "D" ){
		continue;
	    }
	    

	    double q   = destringify<double>( fields[4] );
	    if(q>=qcutoff){
		double aprob   = destringify<double>( fields[ 7] );
		double cprob   = destringify<double>( fields[ 8] );
		double gprob   = destringify<double>( fields[ 9] );
		double tprob   = destringify<double>( fields[10] );
		double sum = aprob+cprob+gprob+tprob;
		if(sum == 0 )
		    cout<<pos<<"\t"<<abs(0.25)<<"\t"<<abs(0.25)<<"\t"<<abs(0.25)<<"\t"<<abs(0.25)<<endl;
		else
		    cout<<pos<<"\t"<<abs(aprob/sum)<<"\t"<<abs(cprob/sum)<<"\t"<<abs(gprob/sum)<<"\t"<<abs(tprob/sum)<<endl;
	    }
	}
	logFD.close();

    }else{
	cerr << "Cannot open log file  "<<logFile<<""<<endl;
	return 1;
    }




    cerr<<"program "<<argv[0]<<" finished successfully"<<endl;

    return 0;
}

