#include <iostream>
#include <string>
#include <cstdlib>
#include <sys/time.h>

#include <api/BamConstants.h>
#include <api/BamMultiReader.h>
#include <utils/bamtools_fasta.h>
#include <utils/bamtools_options.h>
#include <utils/bamtools_pileup_engine.h>
#include <utils/bamtools_utilities.h>

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"

#include "libgab.h"

using namespace std;
using namespace BamTools;

#define MIN(a,b) (((a)<(b))?(a):(b))


class coverageComputeVisitor : public PileupVisitor {
  
public:
    coverageComputeVisitor(const RefVector& references)
	: PileupVisitor()
	, m_references(references)
    { 
	totalBases=0;
	totalSites=0;	
    }
    ~coverageComputeVisitor(void) {};

    // PileupVisitor interface implementation


    void Visit(const PileupPosition& pileupData) {   

	for(unsigned int i=0;i<pileupData.PileupAlignments.size();i++){
	    if( pileupData.PileupAlignments[i].IsCurrentDeletion &&
		pileupData.PileupAlignments[i].IsNextInsertion ){
		continue;
	    }
	    
	    totalBases++;
	}
	
	
	totalSites++;
    }
    
    unsigned int getTotalBases() const{
	return totalBases;
    }
    
    unsigned int getTotalSites() const{
	return totalSites;
    }

private:
    RefVector m_references;
    //Fasta * m_fastaReference;
    unsigned int totalBases;
    unsigned int totalSites;
    
};//end coverageComputeVisitor




int main (int argc, char *argv[]) {

    if( (argc== 1) ||
	(argc== 2 && string(argv[1]) == "-h") ||
	(argc== 2 && string(argv[1]) == "-help") ||
	(argc== 2 && string(argv[1]) == "--help") ){
	cerr<<"Usage:"<<endl;
	cerr<<""<<endl;
	cerr<<"subSampleBAM [COV] in.bam out.bam"<<endl;
	cerr<<""<<endl;
	cerr<<"\twhere [COV] is the desired coverage"<<endl;
	cerr<<endl;
	return 1;
    }
    
    double  coverage      = atof  (argv[1]);
    string bamfiletoopen   = string(argv[2]);
    string outputFilename = string(argv[3]);

    //computing coverage

    

    cerr<<"Computing coverage..";
    BamReader reader_;
    if ( !reader_.Open(bamfiletoopen) ) {
	cerr << "Could not open input BAM file:" << bamfiletoopen <<endl;
    	exit(1);
    }
    
    // reader_.LocateIndex();

    // if(!reader_.HasIndex()){
    // 	cerr << "The BAM file: " << bamfiletoopen <<" does not have an index"<<endl;
    // 	exit(1);
    // }

    // retrieve reference data
    const RefVector  references_ = reader_.GetReferenceData();
    //const int        refID      = reader_.GetReferenceID( currentChunk->rangeGen.getChrName() );




    // BamRegion bregion (refID, 
    // 		       int(currentChunk->rangeGen.getStartCoord()), 
    // 		       refID, 
    // 		       int(currentChunk->rangeGen.getEndCoord() )  );

    // bool setRegionRes = reader.SetRegion( bregion   );

    // if( refID==-1 ||
    //    !setRegionRes){	
    // 	cerr << "Coverage computation: could not set region "<<currentChunk->rangeGen<<" for BAM file:" << bamFileToOpen <<" "<<currentChunk->rangeGen.getStartCoord()<<" "<< currentChunk->rangeGen.getEndCoord()<< "\trefID:"<<refID<<"\tset region fail?:"<<booleanAsString(setRegionRes)<<endl;
    // 	exit(1);
    // }


    coverageComputeVisitor* cv = new coverageComputeVisitor(references_);// ,
							    // currentChunk->rangeGen.getStartCoord(),
							    // currentChunk->rangeGen.getEndCoord()
    //);
    PileupEngine pileup;
    pileup.AddVisitor(cv);

    BamAlignment al_;
    //unsigned int numReads=0;
    while ( reader_.GetNextAlignment(al_) ) {
        pileup.AddAlignment(al_);
    }

    
    //clean up
    pileup.Flush();
    reader_.Close();
    //fastaReference.Close();

    unsigned int totalBasesL=cv->getTotalBases();
    unsigned int totalSitesL=cv->getTotalSites();

    double coverageObserved=double(totalBasesL)/double(totalSitesL);
    delete cv;
    //end computing coverage

    double fraction = MIN( 1.0 , (coverage/coverageObserved) );
    cerr<<"..done";    

    cerr<<"Coverage found "<<coverageObserved<<" target "<<coverage<<" subsampling fraction "<<fraction<<endl;









    BamReader reader;

    if( coverage <= 0.0 ){
    	cerr << "The coverage must be above 0" << endl;
    	return 1;
    }

    if ( !reader.Open(bamfiletoopen) ) {
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
    map<string,bool> ids2add;

    BamAlignment al;
    while ( reader.GetNextAlignment(al) ) {

	if(al.IsPaired()){

	    if(ids2add.find(al.Name) == ids2add.end()){ //check if al.Name has been not found
		
		if(randomProb() < fraction){//if first mate is selected, the second has to be as well
		    writer.SaveAlignment(al);
		    ids2add[al.Name] = true;
		}else{
		    ids2add[al.Name] = false;
		}

	    }else{//found
		
		if(ids2add[al.Name] ){
		    writer.SaveAlignment(al);
		}
		ids2add.erase(al.Name);

	    }



	}else{

	    if( randomProb() < fraction){
		writer.SaveAlignment(al);
	    }	    

	}
    }
    reader.Close();
    writer.Close();
   
    return 0;
}

