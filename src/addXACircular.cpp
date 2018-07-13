#include <iostream>

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"
#include "utils.h"

using namespace std;
using namespace BamTools;


int main (int argc, char *argv[]) {

    if( (argc== 1) ||
	(argc== 2 && string(argv[1]) == "-h") ||
	(argc== 2 && string(argv[1]) == "-help") ||
	(argc== 2 && string(argv[1]) == "--help") ){
	cout<<"Usage:"<<argv[0]<<"<options>   [in bam] [ref original length] [extension length]"<<endl;
	cout<<"This program returns the same BAM file except with the XA flag for circular references"<<endl;
	cout<<"Options:"<<endl;
	//cout<<"\t-m\t\t\t\tUse mapped reads only"<<endl;
	return 1;
    }
     

    for(int i=1;i<(argc-3);i++){ //all but the last 3 args

	// if(string(argv[i]) == "-m" ){
	//     onlyMapped=true;
	//     continue;
	// }
      
	cerr<<"Unknown option "<<argv[i] <<" exiting"<<endl;
	return 1;
    }


     string bamfiletopen  = string(          argv[argc-3]);
     int  origLength      = destringify<int>(argv[argc-2]);
     int  extLength       = destringify<int>(argv[argc-1]);
     
     string outputFilename = "/dev/stdout";

     BamReader reader;
     
     if ( !reader.Open(bamfiletopen) ) {
	 cerr << "Could not open input BAM files." << endl;
	 return 1;
     }

    const SamHeader header     = reader.GetHeader();
    const RefVector references = reader.GetReferenceData();

    BamWriter writer;
    if ( !writer.Open(outputFilename, header, references) ) {
	cerr << "Could not open output BAM file" << endl;
	return 1;
    }
    
    BamAlignment al;
    string nameTAG="XA";
    while ( reader.GetNextAlignment(al) ) {
	
	if(al.HasTag(nameTAG)) {
	    cerr << "ERROR: Read "<<al.Name<<" already has XA tags" << endl;
	    return 1;	    
	}
	
	writer.SaveAlignment(al);

		    
     } //while al

     reader.Close();
     writer.Close();

     return 0;
}

