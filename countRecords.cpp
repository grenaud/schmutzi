/*
 * failQualPair
 * Date: Oct-10-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */


#include <iostream>

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"
#include "utils.h"

using namespace std;
using namespace BamTools;


int main (int argc, char *argv[]) {
    bool onlyMapped=false;

    if( (argc== 1) ||
	(argc== 2 && string(argv[1]) == "-h") ||
	(argc== 2 && string(argv[1]) == "-help") ||
	(argc== 2 && string(argv[1]) == "--help") ){
	cout<<"Usage:"<<argv[0]<<"<options>   [in bam]"<<endl<<"this program returns the number of sequences on the stdout"<<endl;
	cerr<<"Options:"<<endl;
	cerr<<"\t-m\t\t\t\tUse mapped reads only"<<endl;
	return 1;
    }
     

    for(int i=1;i<(argc-1);i++){ //all but the last arg

	if(string(argv[i]) == "-m" ){
	    onlyMapped=true;
	    continue;
	}
      
	cerr<<"Unknown option "<<argv[i] <<" exiting"<<endl;
	return 1;
    }


     string bamfiletopen = string(argv[argc-1]);
     // cout<<bamfiletopen<<endl;
     BamReader reader;
     // cout<<"ok"<<endl;
     if ( !reader.Open(bamfiletopen) ) {
	 cerr << "Could not open input BAM files." << endl;
	 return 1;
     }

     unsigned int total=0;
     BamAlignment al;
     // cout<<"ok"<<endl;
     while ( reader.GetNextAlignment(al) ) {
	 total++;
	 // // cout<<al.Name<<endl;
	 // if(!al.IsMapped())
	 //     continue;

	 // if(al.IsPaired()){
	 //     if( al.IsFirstMate()  ){
	 // 	 if(al.InsertSize > 0)
	 // 	     cout<<al.InsertSize<<endl;
	 // 	 else
	 // 	     cout<<-1.0*al.InsertSize<<endl;
	 //     }
	 // }else{
	 //     cout<<al.Length<<endl;
	 // }		   
     } //while al

     reader.Close();

     cout<<total<<endl;

     return 0;
}

