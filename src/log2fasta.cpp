#include <gzstream.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>


#include "libgab.h"


using namespace std;




int main (int argc, char *argv[]) {

    double qcutoff      = 0;
    double qcutoffIndel = 0;

    string name = "MT";
    const string usage=string("\t"+string(argv[0])+
			      " This program transform a log file into a contamination profile "+
			      " [options]  [log file]  "+"\n\n"+			      
			      "\n\tOptions:\n"+	

			      "\t\t"+"-q"+"\t[quality cutoff]" +"\t\t"+"Only consider sites with quality greater than this cutoff (default: 0)"+"\n"+
			      "\t\t"+"-name"+"\t[name]" +"\t\t\t\t"+"Name of the sequence (default "+stringify(name)+")\n"+ 
			      "\t\t"+"-indel"+"\t[quality cutoff]" +"\t\t"+"Filter indels according to this threshold (default: 0)"+"\n"+
			      "\t\t"+"\t"+"\t\t" +"\t\t"+"Warning: This means that for indels, you will produce what the reference has"+"\n"+
			      "\t\t"+"\t"+"\t\t" +"\t\t"+"         hence adding a potential reference bias"+"\n"+

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

	if(string(argv[i]) == "-indel" ){
	    qcutoffIndel = destringify<double>(string(argv[i+1]));
	    i++;
	    continue;
	}

	if(string(argv[i]) == "-name" ){
	    name = string(argv[i+1]);
	    i++;
	    continue;
	}
       
	cerr<<"Unknown option "<<argv[i] <<" exiting"<<endl;
	return 1;
    }


    string logFile     = string(argv[argc-1]);

    if(name[0] != '>')
	name=">"+name; 

    string line;
    if(isDos(logFile) || isMac(logFile) ){
	cerr << "File  "<<logFile<<" must be unix formatted, exiting"<<endl;
	return 1;
    }


    igzstream logFD;
    string sequence = "";

    logFD.open(logFile.c_str(),ios::in);
    if (logFD.good()){
	getline (logFD,line); //header
	while ( getline (logFD,line)){
	    if(line.empty())
		continue;
	    //cout<<line<<"\t"<<sequence<<endl;
	    vector<string> fields = allTokens(line,'\t');
	    string  pos    =                      fields[0];
	    string  ref    =                      fields[1];
	    string  alt    =                      fields[2];
	    double  q;
	    if(fields[3]=="inf")
		    q =  numeric_limits<double>::infinity();
		else
		    q =  destringify<double>(fields[3]);

	    
	    if(alt  == "D" ){ //deletion
		if(q<qcutoffIndel){//skip those below threshold
		    sequence+="N";
		}else{
		    sequence+=""; //do nothing, deletion is above QC threshold
		}       
		continue;
	    }



		
	    //beyond this point, only lines with sufficient QC

	    if(strEndsWith(pos,"i")){ // indel
		//cout<<line<<"\t"<<sequence<<endl;
		if(q<qcutoffIndel)//skip those below threshold		 
		    continue;
		
		sequence+=alt;
		//cout<<line<<"\t"<<sequence<<endl;
		continue;
	    }

	    //not deletion or insertion

	    if(q<qcutoff)//skip those below threshold
		sequence+="N";       	
	    else
		sequence+=alt;
	    
	}
	logFD.close();

    }else{
	cerr << "Cannot open log file  "<<logFile<<""<<endl;
	return 1;
    }


    cerr<<"program "<<argv[0]<<" finished successfully"<<endl;

    cout<<name<<endl;
    for(unsigned int i=0;i<sequence.size();i++){
	cout<<sequence[i];
	if( (i!=0) && ( (i%80)  == 79) ){
	    cout<<endl;
	}
    }
    cout<<endl;
    return 0;
}

