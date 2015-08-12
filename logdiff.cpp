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

    double qcutoff      = 0;
    double qcutoffIndel = 0;

    string name = "MT";
    // bool verbose=false;

    const string usage=string("\t"+string(argv[0])+
			      " This program finds lines in multiple logs where they differ "+
			      " [options]  [log file 1] [log file 2] .. "+"\n\n"+			       
			      "\n\tOptions:\n"+	

			      //"\t\t"+"-v"+"\t\t" +"\t\t"+"Print the consensus and the input sequences"+"\n"+
			      // "\t\t"+"-name"+"\t[name]" +"\t\t\t\t"+"Name of the sequence (default "+stringify(name)+")\n"+ 
			      "\t\t"+"-q"+"\t[quality cutoff]" +"\t\t"+"Only consider sites with quality greater than this cutoff (default: 0)"+"\n"+
			      "\t\t"+"-indel"+"\t[quality cutoff]" +"\t\t"+"Filter indels according to this threshold (default: 0)"+"\n"+
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
    int lastOpt=1;

    for(int i=1;i<(argc-1);i++){ //all but the last arg
	if(string(argv[i])[0] != '-'  ){
            //cout<<"end"<<i<<endl;
            lastOpt=i;
            break;
        }

	// if(string(argv[i]) == "-v" ){
	//     verbose=true;
	//     continue;
	// }
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

	// if(string(argv[i]) == "-name" ){
	//     name = string(argv[i+1]);
	//     i++;
	//     continue;
	// }
       
	cerr<<"Unknown option "<<argv[i] <<" exiting"<<endl;
	return 1;
    }

    vector<string> logFiles;
    for(int i=lastOpt;i<(argc);i++){ 
	string logFile     = string(argv[i]);
	//cout<<logFile<<endl;
	logFiles.push_back( logFile );
	if(isDos(logFile) || isMac(logFile) ){
	    cerr << "File  "<<logFile<<" must be unix formatted, exiting"<<endl;
	    return 1;
	}

    // logFD.open(logFile.c_str(),ios::in);
    // if (logFD.good()){
    }

    string line;
    vector< igzstream * > logFDs;
    // logFDs.reserve(logFiles.size());
    for(unsigned int i=0;i<(logFiles.size());i++){ 
    	string logFile     = logFiles[i];
	igzstream * toadd = new igzstream();
    	toadd->open(logFile.c_str(),ios::in);

    	if (toadd->good()){

	    getline (*toadd,line);//skip header
    	    //toadd->close();
	    
    	}else{
    	    cerr << "Cannot open log file  "<<logFile<<""<<endl;
    	    return 1;
    	}

	logFDs.push_back(toadd);
    }

    unsigned int currentPos=1;

    vector<string> lines     (logFDs.size(),"");
    vector<bool> insertion   (logFDs.size(),false);

    //bool breakLoop;
    while(true){
	bool allDone=true;
	for(unsigned int i=0;i<(logFDs.size());i++){ 

	    if(!insertion[i]){
		getline ( *(logFDs[i])  ,line);//skip header
		if(line.empty()){
		    //break;
		    continue;
		}else{
		    allDone=allDone&&false;
		}
		lines[i]=line;
	    }else{
		//read nothing if it had an insertion
	    }

	    //cout<<i<<"\t"<<currentPos<<"\t"<<lines[i]<<endl;
	    vector<string> fields = allTokens(lines[i],'\t');
	    string  pos    =                      fields[0];

	    string  ref    =                      fields[1];
	    string  alt    =                      fields[2];
	 
	    if(strEndsWith(pos,"i")){ // indel
		//deal with indel
		insertion[i]=true;
		continue;
	    }else{
		insertion[i]=false;
	    }

	    unsigned int  posUI    = destringify<unsigned int>(pos);

	    
	    if(currentPos!=posUI){
		cerr << "Error with position "<<currentPos<<" at line:  "<<lines[i]<<" in file "<<logFiles[i]<<endl;
		return 1;
		
	    }    
	}

	// cout<<"cp"<<currentPos<<endl;
	//if(verbose){
	//cout<<"----"<<endl;
	vector<string> fields0 = allTokens(lines[0],'\t');
	
	string  alt       =                      fields0[2];
	long double qualB0 =  numeric_limits<long double>::infinity(); //destringify<long double>(fields0[2] );
	long double qualI0 =  numeric_limits<long double>::infinity();
	if(strEndsWith(fields0[0],"i") || alt == "D" ){
	    qualI0=destringify<long double>( fields0[3] );
	}else{
	    qualB0=destringify<long double>( fields0[3] );
	}
	bool diff=false;
	for(unsigned int i=1;i<(lines.size());i++){ 
	    vector<string> fields = allTokens(lines[i],'\t');
	    long double qualB =  numeric_limits<long double>::infinity(); //destringify<long double>(fields0[2] );
	    long double qualI =  numeric_limits<long double>::infinity();

	    if(strEndsWith(fields[0],"i") || alt == "D"){
		qualI=destringify<long double>( fields[3] );
	    }else{
		qualB=destringify<long double>( fields[3] );
	    }

	    if(alt    != fields[2]){
		if( qualB0 > qcutoff &&
		    qualB  > qcutoff &&
		    qualI0 > qcutoffIndel &&
		    qualI  > qcutoffIndel ){		    
		    diff=true;
		}
	    }	    
	}
	if(diff){
	    cout<<vectorToString(lines,"\n")<<endl;
	}
	//     cout<<"----"<<endl;
	// }
	// cout<<vectorToString(insertion)<<endl;
	bool noneHaveIndel=true;
	for(unsigned int i=0;i<(logFDs.size());i++){ 
	    //if(!insertion[i]){
	    noneHaveIndel=noneHaveIndel&&(!insertion[i]);
	}


	//cout<<noneHaveIndel<<endl;
	if(!noneHaveIndel){




	    for(unsigned int i=0;i<(logFDs.size());i++){ 
		insertion[i] = !insertion[i];
	    }


	    
	}
	if(noneHaveIndel){
	    currentPos++;
	}
	if(allDone){
	    break;
	}
    }


     cerr<<"program "<<argv[0]<<" finished successfully"<<endl;
    return 0;
}

