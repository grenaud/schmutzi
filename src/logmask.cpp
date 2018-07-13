#include <gzstream.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>


#include "utils.h"


using namespace std;

typedef struct{
    int start;
    int end;   
} region;


int main (int argc, char *argv[]) {



    const string usage=string("\t"+string(argv[0])+
			      " This program transform a log file by masking the positions in a bed file  "+
			      " [options]  [log file] [bed file] "+"\n\n"+			      
			      "\n\tOptions:\n"+	

			      // "\t\t"+"-q"+"\t[quality cutoff]" +"\t\t"+"Only consider sites with quality greater than this cutoff (default: 0)"+"\n"+
			      // "\t\t"+"-name"+"\t[name]" +"\t\t\t\t"+"Name of the sequence (default "+stringify(name)+")\n"+ 
			      // "\t\t"+"-indel"+"\t[quality cutoff]" +"\t\t"+"Filter indels according to this threshold (default: 0)"+"\n"+
			      // "\t\t"+"\t"+"\t\t" +"\t\t"+"Warning: This means that for indels, you will produce what the reference has"+"\n"+
			      // "\t\t"+"\t"+"\t\t" +"\t\t"+"         hence adding a potential reference bias"+"\n"+

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

	// if(string(argv[i]) == "-q" ){
	//     qcutoff = destringify<double>(string(argv[i+1]));
	//     i++;
	//     continue;
	// }

	// if(string(argv[i]) == "-indel" ){
	//     qcutoffIndel = destringify<double>(string(argv[i+1]));
	//     i++;
	//     continue;
	// }

	// if(string(argv[i]) == "-name" ){
	//     name = string(argv[i+1]);
	//     i++;
	//     continue;
	// }
       
	cerr<<"Unknown option "<<argv[i] <<" exiting"<<endl;
	return 1;
    }


    string bedFile     = string(argv[argc-1]);
    string logFile     = string(argv[argc-2]);

    // if(name[0] != '>')
    // 	name=">"+name; 

    string line;
    if(isDos(logFile) || isMac(logFile) ){
	cerr << "File  "<<logFile<<" must be unix formatted, exiting"<<endl;
	return 1;
    }

    vector<region> regions;
    igzstream bedFD;
    bedFD.open(bedFile.c_str(),ios::in);
    if (bedFD.good()){

	while ( getline (bedFD,line)){
	    if(line.empty())
		continue;
	    vector<string> fields = allTokens(line,'\t');
	    //string  pos    =                      fields[0];
	    //cout<<line<<endl;
	    region toadd;
	    toadd.start    =    destringify<int>(fields[1])+1;
	    toadd.end      =    destringify<int>(fields[2]);
	    if(toadd.start > toadd.end){
		cerr << "Error at line "<<line<<", start is greater than end"<<endl;
		return 1;
	    }
	    regions.push_back(toadd);
	}
	bedFD.close();

    }else{
	cerr << "Cannot open bed file  "<<bedFile<<""<<endl;
	return 1;
    }

    for(unsigned int i=1;i<regions.size();i++){
	if(! ( (regions[i-1].start < regions[i].start ) &&
	       (regions[i-1].end   < regions[i].end   ) ) ){
	    cerr << "Error with line: "<<regions[i].start<<" "<<regions[i].end<<", file is not sorted"<<endl;
	    return 1;	    
	}
    }

    igzstream logFD;
    string sequence = "";

    logFD.open(logFile.c_str(),ios::in);
    if (logFD.good()){
	getline (logFD,line); //header
	cout<<line<<endl;
	while ( getline (logFD,line)){
	    if(line.empty())
		continue;
	    vector<string> fields = allTokens(line,'\t');
	    int pos    =destringify<int>(                      fields[0]);
	    string  ref    =                      fields[1];
	    string  alt    =                      fields[2];
	    bool inMasked=false;
	    for(unsigned int i=0;i<regions.size();i++){
		if( pos >= regions[i].start  &&
		    pos <= regions[i].end    ){
		    inMasked=true;
		}		    

	    }
	    if(inMasked){
		cout<<pos<<"\t"<<ref<<"\t"<<alt<<"\t0.0\t"<<fields[4]<<"\t"<<fields[5]<<"\t"<<fields[6]<<"\t"<<"0.0\t0.0\t0.0\t0.0"<<endl;
	    }else{
		cout<<line<<endl;
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

