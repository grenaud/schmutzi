/*
 * parseMSA [multiple sequence align] [ref ID]
 * Date: Jun-06-2014 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <set>
#include <algorithm>
#include <string>
#include <fstream>
#include <gzstream.h>
#include "FastQParser.h"
#include <algorithm>   
#include "utils.h"

typedef struct { 
    string id;
    double div;
} iddiv ;

using namespace std;

// bool compdiv (iddiv i,iddiv j) { return (i.div<j.div); }

int main (int argc, char *argv[]) {


    const string usage=string("\nThis program takes an multiple sequence alignment (msa) and outputs the allele frequency for each sequence in the freqs/ directory which must be created in the current working directory.\n\n\t"+
			      string(argv[0])+                        
                              "  [msa]  [name of human reference]"+"\n\n");


    if( (argc!= 3) ||
        (argc== 2 && string(argv[1]) == "-h") ||
        (argc== 2 && string(argv[1]) == "-help") ||
        (argc== 2 && string(argv[1]) == "--help") ){
        cout<<usage<<endl;
        return 1;
    }


    
    

    string IDREF = string(argv[2]);
    FastQParser fqp (string(argv[1]),true );

    vector<bool> isOurSeq;
    vector<string> ids;
    vector<string> seq;


    bool found=false;
    string seqREF;
    while(fqp.hasData()){
	FastQObj * test	=fqp.getData();
	
	string strSeq = *(test->getSeq()); 
	//cout<<strSeq<<endl;
	transform(strSeq.begin(), strSeq.end(),strSeq.begin(), ::toupper);
       

	if(*(test->getID())  == IDREF){
	    found=true;
	    isOurSeq.push_back(true);
	    seqREF = strSeq;
	}else{
	    isOurSeq.push_back(false);
	}

	ids.push_back( *(test->getID()) );	
	seq.push_back( strSeq );

    }

    if(!found){
	cerr<<"ID "<<IDREF<<" was not found"<<endl;
	return 1;       
    }

    // cout<<vectorToString(isOurSeq)<<endl;

    // return 1;

    set<string> idsused;

    vector<iddiv> iddivs;
    string outputs[4] = {
	"1.0\t0.0\t0.0\t0.0",
	"0.0\t1.0\t0.0\t0.0",
	"0.0\t0.0\t1.0\t0.0",
	"0.0\t0.0\t0.0\t1.0"
    };
       
    for(unsigned int i=0;i<ids.size();i++){
	// cout<<i<<endl;
	// cout<<isOurSeq[i]<<endl;

	if(!isOurSeq[i]){ //not ref
	    //cout<<ids[i]<<endl;
	    vector<string> tokens=allTokens(ids[i],' ' );
	    string idtouse;
	    int idIdx=1;
	    while(true){
		idtouse=vectorToString(tokens,"-").substr(1)+stringify(idIdx);
		// cout<<idtouse<<endl;
		if( idsused.find(idtouse) == idsused.end() ){
		    idsused.insert(idtouse);
		    break;		   
		}else{
		    idIdx++;
		}		  
	    }
	   
	    
	    cerr<<"Processing "<<idtouse<<endl;
	    replace( idtouse.begin(), idtouse.end(), '|', '_');//to avoid the | from GenBank
	    string filenameout = "freqs/"+idtouse + ".freq";


	    ofstream outProf;


	    outProf.open(filenameout.c_str());
	    if( !outProf.is_open() ){
		cerr<<"ERROR: cannot write to "<<filenameout<<" does the freqs/ directory exists?"<<endl;
		return 1;
	    }

	    // int totalL=0;
	    // int match=0;

	    if(seqREF.size() != seq[i].size()){
		cerr<<"Size between ref and target differ"<<endl;
		return 1;
	    }

	    int indexRef=1;

	    for(unsigned int j=0;j<seqREF.size();j++){
		if(seqREF[j] == '-' ){
		    //useless and do not increase counter
		    continue;
		}

		if(seqREF[j] == 'N' ){
		    //useless and do not increase counter, put a bunch of zeros even though the other might be defined.
		    outProf<<indexRef<<"\t0.0\t0.0\t0.0\t0.0"<<endl;
		    indexRef++;
		    continue;
		}


		
		if(seq[i][j] == 'N' ||
		   seq[i][j] == '-' ||
		   isAmbiguousUIPAC(seq[i][j])){
		    outProf<<indexRef<<"\t0.0\t0.0\t0.0\t0.0"<<endl;
		}else{	
		    //cout<<indexRef<<"\t"<< idtouse<<"\t"<<(seqREF[j])<<"\t"<<(seq[i][j])<<endl;
		    outProf<<indexRef<<"\t"<<outputs[ baseResolved2int(seq[i][j]) ]<<endl;
		}

		indexRef++;
   	      
	    }

	    outProf.close();
	    //cout<<ids[i]<<"\t"<<(double(match)/double(totalL))<<endl;
   	  
	}

    }
   
    // sort (iddivs.begin(), iddivs.end(), compdiv);

   
    // for(unsigned int i=0;i<iddivs.size();i++){
    //     cout<<iddivs[i].id<<"\t"<<iddivs[i].div<<endl;
    // }
       
    cerr<<"Finished gracefully"<<endl;

    return 0;
}

