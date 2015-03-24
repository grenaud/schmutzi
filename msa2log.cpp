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

    double qual=100;
    const string usage=string("\nThis program takes an pairwise sequence alignment with the mt reference and outputs the per-position log of non-reference sequence \n\t"+
			      string(argv[0])+                        
                              "  [msa]  [name of reference sequence]\n\n"+
			      "\t\t"+"-q [qual]" +"\t\t\t"+"Quality on a PHRED scale for individual bases (Default: "+stringify(qual)+" )"+"\n"+
			      "\n\n");
    

    if( (argc < 2) ||
        (argc== 2 && string(argv[1]) == "-h") ||
        (argc== 2 && string(argv[1]) == "-help") ||
        (argc== 2 && string(argv[1]) == "--help") ){
        cout<<usage<<endl;
        return 1;
    }


    
    
    // int lastOpt=1;
    for(int i=1;i<(argc-2);i++){ //all but the last 3 args

        if(string(argv[i]) == "-q" ){
            qual = destringify<double>(argv[i+1]);
	    i++;
            continue;
        }

        // if(string(argv[i])[0] != '-'  ){
        //     //cout<<"end"<<i<<endl;
        //     lastOpt=i;
        //     break;
        // }



        cerr<<"Wrong option "<<string(argv[i])<<endl;
        return 1;

    }



    FastQParser fqp (string(argv[argc-2]),true );
    string IDREF   = string(argv[argc-1]);

    if(!strBeginsWith(IDREF,">")){
	IDREF = ">"+IDREF;
    }

    //vector<bool> isOurSeq;
    vector<string> ids;
    vector<string> seq;
    int indexRef=-1;
    int indexTarget=-1;

    bool found=false;
    string seqREF;


    while(fqp.hasData()){
	FastQObj * test	=fqp.getData();
	
	string strSeq = *(test->getSeq()); 
	//cout<<strSeq<<endl;
	transform(strSeq.begin(), strSeq.end(),strSeq.begin(), ::toupper);
       
	//cout<<*(test->getID())<<endl;
	if(*(test->getID())  == IDREF){
	    found=true;
	    seqREF = strSeq;

	    indexRef    = int(ids.size());
	}else{
	    indexTarget = int(ids.size());
	}

	ids.push_back( *(test->getID()) );	
	seq.push_back( strSeq );

    }

    if(!found){
	cerr<<"ID "<<IDREF<<" was not found"<<endl;
	return 1;       
    }


    if(ids.size() != 2){
	cerr<<"Cannot have more than 2 sequences in the pairwise alignment"<<endl;
	return 1;       
    }

    // cout<<vectorToString(isOurSeq)<<endl;

    // return 1;

    // set<string> idsused;

    // vector<iddiv> iddivs;
    // string outputs[4] = {
    // 	"1.0\t0.0\t0.0\t0.0",
    // 	"0.0\t1.0\t0.0\t0.0",
    // 	"0.0\t0.0\t1.0\t0.0",
    // 	"0.0\t0.0\t0.0\t1.0"
    // };


    if(seq[indexRef].size() != seq[indexTarget].size()){
	cerr<<"Size between ref and target differ"<<endl;
	return 1;
    }

    int indexOnTheReference=1;
    cout<<"pos\trefBase\tbase\tqual\tavgmapq\tcov\tsupp\tpa\tpc\tpg\tpt"<<endl;

    for(unsigned int j=0;j<seq[indexRef].size();j++){
	if(seq[indexRef][j] == '-' ){
	    cout<<(indexOnTheReference-1)<<"i\t"<<seq[indexRef][j]<<"\t"<<seq[indexTarget][j]<<"\t"<<qual<<"\t"<<250<<"\t"<<100<<"\t"<<100<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<endl;
	    //useless and do not increase counter
	    continue;
	}

	if(seq[indexRef][j] == 'N' ){
	    //useless and do not increase counter, put a bunch of zeros even though the other might be defined.
	    //cout<<indexRef<<"\t0.0\t0.0\t0.0\t0.0"<<endl;
	    cout<<indexOnTheReference<<"\t"<<seq[indexRef][j]<<"\t"<<"D"<<"\t"<<0<<"\t"<<250<<"\t"<<100<<"\t"<<100<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<endl;
	    indexOnTheReference++;
	    continue;
	}

	if(seq[indexTarget][j] == '-' ){
	    cout<<indexOnTheReference<<"\t"<<seq[indexRef][j]<<"\t"<<"D"<<"\t"<<qual<<"\t"<<250<<"\t"<<100<<"\t"<<100<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<endl;
	    indexOnTheReference++;
	    continue;
	}
		
	if(                 seq[indexTarget][j] == 'N' ||
	   isAmbiguousUIPAC(seq[indexTarget][j])){
	    cout<<indexOnTheReference<<"\t"<<seq[indexRef][j]<<"\t"<<seq[indexTarget][j]<<"\t"<<0<<"\t"<<250<<"\t"<<100<<"\t"<<100<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<endl;
	 
	}else{	
	    //cout<<indexRef<<"\t"<< idtouse<<"\t"<<(seqREF[j])<<"\t"<<(seq[i][j])<<endl;
	    //cout<<indexRef<<"\t"<<outputs[ baseResolved2int(seq[i][j]) ]<<endl;
	    cout<<indexOnTheReference<<"\t"<<seq[indexRef][j]<<"\t"<<seq[indexTarget][j]<<"\t"<<qual<<"\t"<<250<<"\t"<<100<<"\t"<<100;
	    if(seq[indexTarget][j] == 'A')
		cout<<"\t"<<qual<<"\t"<<0<<"\t"<<0<<"\t"<<0<<endl;
	    if(seq[indexTarget][j] == 'C')
		cout<<"\t"<<0<<"\t"<<qual<<"\t"<<0<<"\t"<<0<<endl;
	    if(seq[indexTarget][j] == 'G')
		cout<<"\t"<<0<<"\t"<<0<<"\t"<<qual<<"\t"<<0<<endl;
	    if(seq[indexTarget][j] == 'T')
		cout<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<qual<<endl;

	    
	}

	indexOnTheReference++;
   	      
    }

    // outProf.close();

       
    // for(unsigned int i=0;i<ids.size();i++){
    // 	// cout<<i<<endl;
    // 	// cout<<isOurSeq[i]<<endl;

    // 	if(!isOurSeq[i]){ //not ref
    // 	    //cout<<ids[i]<<endl;
    // 	    vector<string> tokens=allTokens(ids[i],' ' );
    // 	    string idtouse;
    // 	    int idIdx=1;
    // 	    while(true){
    // 		idtouse=vectorToString(tokens,"-").substr(1)+stringify(idIdx);
    // 		// cout<<idtouse<<endl;
    // 		if( idsused.find(idtouse) == idsused.end() ){
    // 		    idsused.insert(idtouse);
    // 		    break;		   
    // 		}else{
    // 		    idIdx++;
    // 		}		  
    // 	    }
	   
	    
    // 	    cerr<<"Processing "<<idtouse<<endl;
    // 	    replace( idtouse.begin(), idtouse.end(), '|', '_');//to avoid the | from GenBank
    // 	    string filenameout = "freqs/"+idtouse + ".freq";


    // 	    ofstream outProf;


    // 	    outProf.open(filenameout.c_str());
    // 	    // int totalL=0;
    // 	    // int match=0;

    // 	    //cout<<ids[i]<<"\t"<<(double(match)/double(totalL))<<endl;
   	  
    // 	}

    // }
   
    // // sort (iddivs.begin(), iddivs.end(), compdiv);

   
    // // for(unsigned int i=0;i<iddivs.size();i++){
    // //     cout<<iddivs[i].id<<"\t"<<iddivs[i].div<<endl;
    // // }
       
    cerr<<"Finished gracefully"<<endl;

    return 0;
}

