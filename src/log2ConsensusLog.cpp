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

    // double qcutoff      = 0;
    // double qcutoffIndel = 0;

    string name = "MT";
    bool verbose=false;

    const string usage=string("\t"+string(argv[0])+
			      " This program transform multiple logs into a single one "+
			      " [options]  [log file 1] [log file 2] .. "+"\n\n"+			       
			      "\n\tOptions:\n"+	

			      "\t\t"+"-v"+"\t\t" +"\t\t"+"Print the consensus and the input sequences"+"\n"+
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
    int lastOpt=1;

    for(int i=1;i<(argc-1);i++){ //all but the last arg
	if(string(argv[i])[0] != '-'  ){
            //cout<<"end"<<i<<endl;
            lastOpt=i;
            break;
        }

	if(string(argv[i]) == "-v" ){
	    verbose=true;
	    continue;
	}

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
    cout<<"pos\trefBase\tbase\tqual\tavgmapq\tcov\tsupp\tpa\tpc\tpg\tpt"<<endl;

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
	if(verbose){
	    cout<<"----"<<endl;
	    cout<<vectorToString(lines,"\n")<<endl;
	    cout<<"----"<<endl;
	}
	// cout<<vectorToString(insertion)<<endl;
	bool noneHaveIndel=true;
	for(unsigned int i=0;i<(logFDs.size());i++){ 
	    //if(!insertion[i]){
	    noneHaveIndel=noneHaveIndel&&(!insertion[i]);
	}


	//cout<<noneHaveIndel<<endl;
	if(noneHaveIndel){
	    long double likeNotA=0.0;
	    long double likeNotC=0.0;
	    long double likeNotG=0.0;
	    long double likeNotT=0.0;
	    long double likeNotDel=0.0;

	    unsigned int sumAllBases=0;
	    unsigned int sumA       =0;
	    unsigned int sumC       =0;
	    unsigned int sumG       =0;
	    unsigned int sumT       =0;
	    unsigned int sumBaseDeletion=0;
	    unsigned int sumWithDeletion=0;
	    
	    string ref;

	    for(unsigned int i=0;i<(logFDs.size());i++){ 
		vector<string> fields = allTokens(lines[i],'\t');
		ref                   =          fields[1];
		string  alt           =          fields[2];
		unsigned int sumBases    =  destringify<unsigned int>( fields[5] );
		unsigned int sumBasesAlt =  destringify<unsigned int>( fields[6] );

		sumAllBases+=sumBases;

		if(alt == "D"){
		    sumWithDeletion++;
		    sumBaseDeletion+=sumBasesAlt;
		    likeNotDel += destringify<long double>( fields[ 3] );
		}else{
		    if(alt == "A"){ sumA+=sumBasesAlt;  }
		    if(alt == "C"){ sumC+=sumBasesAlt;  }
		    if(alt == "G"){ sumG+=sumBasesAlt;  }
		    if(alt == "T"){ sumT+=sumBasesAlt;  }
		    likeNotA+=destringify<long double>( fields[ 7] ); 
		    likeNotC+=destringify<long double>( fields[ 8] ); 
		    likeNotG+=destringify<long double>( fields[ 9] ); 
		    likeNotT+=destringify<long double>( fields[10] );
		}

	    }
	    // cout<<likeNotA<<endl;
	    // cout<<likeNotC<<endl;
	    // cout<<likeNotG<<endl;
	    // cout<<likeNotT<<endl;
	    // cout<<"----"<<endl;
	    // cout<<exp(-likeNotA)<<endl;
	    // cout<<exp(-likeNotC)<<endl;
	    // cout<<exp(-likeNotG)<<endl;
	    // cout<<exp(-likeNotT)<<endl;
	    // cout<<"----"<<endl;
	    // cout<<(double(1)-exp(-likeNotA))<<endl;
	    // cout<<(double(1)-exp(-likeNotC))<<endl;
	    // cout<<(double(1)-exp(-likeNotG))<<endl;
	    // cout<<(double(1)-exp(-likeNotT))<<endl;

	    // cout<<"----"<<endl;
	    // long double t = exp(-likeNotA)  + exp(-likeNotC)  + exp(-likeNotG) + exp(-likeNotT) ;
	    // cout<<"ratio "<< ( ((t-exp(-likeNotG))/t)==1 )<<endl;
	    // cout<<t<<endl;
	    // cout<<(t-exp(-likeNotG))<<endl;
	    long double maxLikeBase = -1.0*numeric_limits<long double>::infinity();
	    maxLikeBase = max(maxLikeBase,likeNotA);
	    maxLikeBase = max(maxLikeBase,likeNotC);
	    maxLikeBase = max(maxLikeBase,likeNotG);
	    maxLikeBase = max(maxLikeBase,likeNotT);
	    //cout<<"final\t"<<likeNotDel<<"\t"<<likeNotA<<"\t"<<likeNotC<<"\t"<<likeNotG<<"\t"<<likeNotT<<"\t"<<maxLikeBase<<endl;
	    if(likeNotDel>maxLikeBase){//call deletion
		//call deletion
		cout<<currentPos<<"\t"<<ref<<"\tD"<<"\t"<<log(exp(likeNotDel)/exp( oplus(maxLikeBase,likeNotDel) ) )<<"\t250\t"<<sumAllBases<<"\t"<<sumBaseDeletion<<"\t0.0\t0.0\t0.0\t0.0"<<endl;
	    }else{
		cout<<currentPos<<"\t"<<ref<<"\t";
		long double sumAllNucProb=likeNotA;
		sumAllNucProb=oplus(sumAllNucProb,likeNotC);
		sumAllNucProb=oplus(sumAllNucProb,likeNotG);
		sumAllNucProb=oplus(sumAllNucProb,likeNotT);
		
		long double likeNotASum=0;
		likeNotASum=oplusInit(likeNotASum,likeNotC);
		likeNotASum=oplusInit(likeNotASum,likeNotG);
		likeNotASum=oplusInit(likeNotASum,likeNotT);

		long double likeNotCSum=0;
		likeNotCSum=oplusInit(likeNotCSum,likeNotA);
		likeNotCSum=oplusInit(likeNotCSum,likeNotG);
		likeNotCSum=oplusInit(likeNotCSum,likeNotT);

		long double likeNotGSum=0;
		likeNotGSum=oplusInit(likeNotGSum,likeNotA);
		likeNotGSum=oplusInit(likeNotGSum,likeNotC);
		likeNotGSum=oplusInit(likeNotGSum,likeNotT);

		long double likeNotTSum=0;
		likeNotTSum=oplusInit(likeNotTSum,likeNotA);
		likeNotTSum=oplusInit(likeNotTSum,likeNotC);
		likeNotTSum=oplusInit(likeNotTSum,likeNotG);
		
		// cout<<endl<<"S"<<maxLikeBase<<"\t"<<sumAllNucProb<<endl;
		// cout<<"D"<<log(double(1)-(sumAllNucProb-maxLikeBase))<<endl;
		// cout<<likeNotA<<endl;
		// cout<<likeNotC<<endl;
		// cout<<likeNotG<<endl;
		// cout<<likeNotT<<endl;
		// cout<<"----"<<endl;
	    
		// cout<<likeNotASum<<endl;
		// cout<<likeNotCSum<<endl;
		// cout<<likeNotGSum<<endl;
		// cout<<likeNotTSum<<endl;

		// cout<<sumAllNucProb-likeNotASum<<endl;
		// cout<<sumAllNucProb-likeNotCSum<<endl;
		// cout<<sumAllNucProb-likeNotGSum<<endl;
		// cout<<sumAllNucProb-likeNotTSum<<endl;

		if(likeNotA == maxLikeBase){
		    long double logLikeForPred=(sumAllNucProb-likeNotASum);
		    if(isnan(logLikeForPred)){
			logLikeForPred=maxLikeBase;
		    }
		    cout<<"A"<<"\t"<<logLikeForPred<<"\t250\t"<<sumAllBases<<"\t"<<sumA<<"\t"<<likeNotA<<"\t"<<likeNotC<<"\t"<<likeNotG<<"\t"<<likeNotT;
		}

		if(likeNotC == maxLikeBase){
		    long double logLikeForPred=(sumAllNucProb-likeNotCSum);
		    if(isnan(logLikeForPred)){
			logLikeForPred=maxLikeBase;
		    }

		    cout<<"C"<<"\t"<<logLikeForPred<<"\t250\t"<<sumAllBases<<"\t"<<sumC<<"\t"<<likeNotA<<"\t"<<likeNotC<<"\t"<<likeNotG<<"\t"<<likeNotT;
		}

		if(likeNotG == maxLikeBase){
		    long double logLikeForPred=(sumAllNucProb-likeNotGSum);
		    if(isnan(logLikeForPred)){
			logLikeForPred=maxLikeBase;
		    }

		    cout<<"G"<<"\t"<<logLikeForPred<<"\t250\t"<<sumAllBases<<"\t"<<sumG<<"\t"<<likeNotA<<"\t"<<likeNotC<<"\t"<<likeNotG<<"\t"<<likeNotT;
		}

		if(likeNotT == maxLikeBase){
		    long double logLikeForPred=(sumAllNucProb-likeNotTSum);
		    if(isnan(logLikeForPred)){
			logLikeForPred=maxLikeBase;
		    }

		    cout<<"T"<<"\t"<<logLikeForPred<<"\t250\t"<<sumAllBases<<"\t"<<sumT<<"\t"<<likeNotA<<"\t"<<likeNotC<<"\t"<<likeNotG<<"\t"<<likeNotT;
		}

		cout<<endl;

	    }
	    
	    
	}else{//have insertion
	    vector<string> insertionsVec;
	    long double likeNotIns=0.0;

	    long double likeNotA=0.0;
	    long double likeNotC=0.0;
	    long double likeNotG=0.0;
	    long double likeNotT=0.0;
	    long double likeNotDel=0.0;

	    unsigned int sumAllBases=0;
	    unsigned int sumA       =0;
	    unsigned int sumC       =0;
	    unsigned int sumG       =0;
	    unsigned int sumT       =0;
	    unsigned int sumBaseDeletion=0;
	    unsigned int sumBaseInsertion=0;
	    
	    
	    unsigned int sumWithDeletion=0;
	    unsigned int sumWithInsertions=0;
	    unsigned int sumWithBases =0;

	    string ref;

	    for(unsigned int i=0;i<(logFDs.size());i++){ 
		vector<string> fields = allTokens(lines[i],'\t');
		string  pos    =                      fields[0];
		if(strEndsWith(pos,"i")){ // indel
		    ref                   =          fields[1];
		    string  alt           =          fields[2];
		    insertionsVec.push_back(alt);
		    unsigned int sumBases    =  destringify<unsigned int>( fields[5] );
		    unsigned int sumBasesAlt =  destringify<unsigned int>( fields[6] );
		    likeNotIns += destringify<long double>( fields[ 3] );
		    sumBaseInsertion += sumBasesAlt;
		    sumAllBases      += sumBases;
		    sumWithInsertions++;
		}else{
		    ref                   =          fields[1];
		    string  alt           =          fields[2];
		    unsigned int sumBases    =  destringify<unsigned int>( fields[5] );
		    unsigned int sumBasesAlt =  destringify<unsigned int>( fields[6] );

		    sumAllBases+=sumBases;

		    if(alt == "D"){
			sumWithDeletion++;
			sumBaseDeletion+=sumBasesAlt;
			likeNotDel += destringify<long double>( fields[ 3] );
		    }else{
			sumWithBases++;
			if(alt == "A"){ sumA+=sumBasesAlt; likeNotA+=destringify<long double>( fields[ 7] );   }
			if(alt == "C"){ sumC+=sumBasesAlt; likeNotC+=destringify<long double>( fields[ 8] );   }
			if(alt == "G"){ sumG+=sumBasesAlt; likeNotG+=destringify<long double>( fields[ 9] );   }
			if(alt == "T"){ sumT+=sumBasesAlt; likeNotT+=destringify<long double>( fields[10] );   }
		    }
		}
	    }

	    long double maxLikeBase =  -1.0*numeric_limits<long double>::infinity();
	    if(sumWithBases >0){
		maxLikeBase = max(maxLikeBase,likeNotA);
		maxLikeBase = max(maxLikeBase,likeNotC);
		maxLikeBase = max(maxLikeBase,likeNotG);
		maxLikeBase = max(maxLikeBase,likeNotT);
	    }

	    //if(likeNotIns>maxLikeBase){//call insertion
	    //cout<<sumWithInsertions<<"\t"<<logFDs.size()<<endl;
	    if( (double(sumWithInsertions)/double(logFDs.size()) ) >=0.5){
		if(insertionsVec.empty()){
		    cerr<<"ERROR insert #1"<<endl;
		    return 1;
		}
		string altInsert= insertionsVec [0];
		for(unsigned int i=1;i<(insertionsVec.size());i++){ 
		    if(altInsert != insertionsVec [i]){
			cerr<<"ERROR insert #2"<<endl;
			return 1;
		    }
		}
		
		// cout<<likeNotIns<<endl;
		// cout<<maxLikeBase<<endl;
		// cout<<oplus(maxLikeBase,likeNotIns)<<endl;
		// cout<<( exp(likeNotIns)/exp( oplus(maxLikeBase,likeNotIns) ))<<endl;
		cout<<(currentPos-1)<<"i\t-"<<"\t"<<altInsert<<"\t"<<likeNotIns<<"\t250\t"<<sumAllBases<<"\t"<<sumBaseInsertion<<"\t0.0\t0.0\t0.0\t0.0"<<endl;
;


	    }
	    
	    // return 1;






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
// if(name[0] != '>')
// 	name=">"+name; 

    // string line;



    // string sequence = "";

    // logFD.open(logFile.c_str(),ios::in);
    // if (logFD.good()){
    // 	getline (logFD,line); //header
    // 	while ( getline (logFD,line)){
    // 	    if(line.empty())
    // 		continue;
    // 	    vector<string> fields = allTokens(line,'\t');
    // 	    string  pos    =                      fields[0];
    // 	    string  ref    =                      fields[1];
    // 	    string  alt    =                      fields[2];
    // 	    double  q;
    // 	    if(fields[3]=="inf")
    // 		    q =  numeric_limits<double>::infinity();
    // 		else
    // 		    q =  destringify<double>(fields[3]);

	    
    // 	    if(alt  == "D" ){ //deletion
    // 		if(q<qcutoffIndel){//skip those below threshold
    // 		    sequence+="N";
    // 		}else{
    // 		    sequence+=""; //do nothing, deletion is above QC threshold
    // 		}       
    // 		continue;
    // 	    }



		
    // 	    //beyond this point, only lines with sufficient QC

    // 	    if(strEndsWith(pos,"i")){ // indel
    // 		//cout<<line<<"\t"<<sequence<<endl;
    // 		if(q<qcutoffIndel)//skip those below threshold		 
    // 		    continue;
		
    // 		sequence+=alt;
    // 		//cout<<line<<"\t"<<sequence<<endl;
    // 		continue;
    // 	    }

    // 	    //not deletion or insertion

    // 	    if(q<qcutoff)//skip those below threshold
    // 		sequence+="N";       	
    // 	    else
    // 		sequence+=alt;
	    
    // 	}
    // 	logFD.close();

    // }else{
    // 	cerr << "Cannot open log file  "<<logFile<<""<<endl;
    // 	return 1;
    // }


     cerr<<"program "<<argv[0]<<" finished successfully"<<endl;

    // cout<<name<<endl;
    // for(unsigned int i=0;i<sequence.size();i++){
    // 	cout<<sequence[i];
    // 	if( (i!=0) && ( (i%80)  == 79) ){
    // 	    cout<<endl;
    // 	}
    // }
    // cout<<endl;
    return 0;
}

