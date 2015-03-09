/*
 * pileup
 * Date: Mar-26-2014 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#define SIZEGENOME 17000 //should be enough

#include "api/BamWriter.h"
#include "api/BamReader.h"

#include <api/BamConstants.h>
#include <api/BamMultiReader.h>
#include <utils/bamtools_fasta.h>
#include <utils/bamtools_options.h>
#include <utils/bamtools_pileup_engine.h>
#include <utils/bamtools_utilities.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <gzstream.h>
// #include "GenomicRange.h"

#include "utils.h"

using namespace BamTools;
using namespace std;

#define MAXCOV 5000
#define MIN(a,b) (((a)<(b))?(a):(b))

char   offsetQual=33;

unsigned int coverageTotal;

class myPileVisitor : public PileupVisitor {
  
    public:
    myPileVisitor(const RefVector& references,
		  const vector< vector< pair<char,string> > > * pos2hap,
		  vector< vector< pair<char,int> > > * pos2hapCounter,
		  map<string,string>  * readname2hap,
		  map<string,vector<int> >  * readname2posmask,
		  const int  minQual
		  )
		  //,int coordToFind, 
		  // BamWriter * writerA,
		  // BamWriter * writerC,
		  // BamWriter * writerG,
		  // BamWriter * writerT
	//)
            : PileupVisitor()
	    , m_references(references)
	    , m_pos2hap(pos2hap)
	    , m_pos2hapCounter(pos2hapCounter)
	    , m_readname2hap(readname2hap)
	    , m_readname2posmask(readname2posmask)
	    , m_minQual( minQual)
	    // , m_coverageCounter(coverageCounter)
	    // , m_coordToFind(coordToFind)

	    //,m_writerA(writerA)

	    // ,m_writerC(writerC)
	    // ,m_writerG(writerG)
	    // ,m_writerT(writerT)

	      //, m_out(out)
        { 

	}
        ~myPileVisitor(void) { }
  
    // PileupVisitor interface implementation
    public:
	// prints coverage results ( tab-delimited )
        void Visit(const PileupPosition& pileupData) {
	    //char referenceBase = 'N';
    
	    // if ( !m_fastaReference->GetBase(pileupData.RefId, pileupData.Position, referenceBase ) ) {
	    // 	cerr << "bamtools convert ERROR: pileup conversion - could not read reference base from FASTA file" << endl;
	    // 	return;
	    // }




	    // if(referenceBase!= 'N'){
	    //m_coverageCounter->at(MIN(pileupData.PileupAlignments.size(),MAXCOV))++;
	    // cout <<m_references[pileupData.RefId].RefName << "\t" 
	    //      <<referenceBase<<"\t"
	    //      << pileupData.Position << endl; ;
	    
	    for(unsigned int i=0;i<pileupData.PileupAlignments.size();i++){
		
		//cout<<"p "<<pileupData.PileupAlignments[i].Alignment.Name<<endl;
		//pileupData.Position is zero based
		//cout<<pileupData.Position<<endl;

		// add counter
		for(unsigned int p=0;p<m_pos2hap->at( pileupData.Position+1 ).size();p++){//for every seg nucleotide at position  ( pileupData.Position+1 )
		    
		    if(pileupData.PileupAlignments[i].Alignment.QueryBases[pileupData.PileupAlignments[i].PositionInAlignment] 
		       == 
		       (*m_pos2hap)[ pileupData.Position+1 ][p].first ){ //if matches the base

			if( int( pileupData.PileupAlignments[i].Alignment.Qualities[i]-offsetQual ) > m_minQual){
			    (*m_readname2hap)[      pileupData.PileupAlignments[i].Alignment.Name ] = (*m_pos2hap)[ pileupData.Position+1 ][p].second;
			    (*m_readname2posmask)[  pileupData.PileupAlignments[i].Alignment.Name ].push_back( pileupData.PileupAlignments[i].PositionInAlignment );
			    (*m_pos2hapCounter)[ pileupData.Position+1 ][p].second++;
			}

		    }

		}

	    }

        }
        
    private:
    RefVector m_references;
    // int m_coordToFind;
    // // Fasta * m_fastaReference;
    // // vector<unsigned int> * m_coverageCounter;
    // //        ostream*  m_out;
    // BamWriter * m_writerA;
    // BamWriter * m_writerC;
    // BamWriter * m_writerG;
    // BamWriter * m_writerT;

    const vector< vector< pair<char,string> > >  * m_pos2hap;
    vector< vector< pair<char,int> > >           * m_pos2hapCounter;
    map<string,string>                           * m_readname2hap;
    map<string,vector<int> >                     * m_readname2posmask;
        int m_minQual; 

};



int main (int argc, char *argv[]) {


    //if(argc != 4){
    int minQual=0;

    const string usage= "Usage:\n\t"+string(argv[0])+"(options) [haplogroup diag. pos.]  [bam file in] [bam suffix out]\n\n"+"This program takes a set of diagnostic positions and splits a bam file according to them.\nPlease note that the bases overlapping the diagnostic positions will have their qualities decreased.\n\nOptions:\n"+
	"\t-q\t[qc score]\tMinimum base quality (Default: "+stringify(minQual)+"\n"+

	"\n";

    // 	return 1;
    // }

                                                                                                    
    if( (argc== 1) ||
        (argc== 2 && string(argv[1]) == "-h") || 
        (argc== 2 && string(argv[1]) == "-help") ||
        (argc== 2 && string(argv[1]) == "--help") ){
        cout<<usage<<endl;
        return 1;  
    }    

    //[outbam prefix] bamfile coord
    vector< vector< pair<char,string> > > pos2hap;
    vector< vector< pair<char,int   > > > pos2hapCounter;

    map<string,string>        readname2hap;
    map<string,vector<int> >  readname2posmask;

    for(unsigned int p=0;p<SIZEGENOME;p++){
	vector< pair<char,string> > toadd;

	pos2hap.push_back(toadd);
	vector< pair<char,int> > toadd2;

	pos2hapCounter.push_back(toadd2);
    }
    

    for(int i=1;i<(argc-3);i++){ //all but the last 3 args                                                                                                                                            
        if(string(argv[i]) == "-q"  ){ 
            minQual=destringify<int>( argv[i+1] );
            continue;
        } 
              
    }

    string haploGroupFile =                    string(argv[argc-3]);
    string bamfiletopen   =                    string(argv[argc-2]);
    string bamFileOUTSuf  =                    string(argv[argc-1]);     
    // cout<<coordToFind<<endl;
    // cout<<bamFileOUT<<endl;
    // cout<<bamfiletopen<<endl;
    //int    coordToFind;
    //string bamFileOUT;
    // return 1;
    //string fastaFile    = string(argv[argc-2]);//fasta file
    // string bedFile      = string(argv[argc-3]);//fasta file
    
    string linePos;
    igzstream myFilePos;
    
    myFilePos.open(haploGroupFile.c_str(), ios::in);
    //cout<<haploGroupFile<<endl;
    if (myFilePos.good()){
	while ( getline (myFilePos,linePos)){
	    //cout<<linePos<<endl;
	    vector<string> fields = allTokens(linePos,'\t');
	    int pos    = destringify<int>(fields[0]);
 	    pair<char,string> baseHap;
 	    pair<char,int>    baseHapCounter;

	    baseHap.first          = destringify<char>(fields[1]);
	    baseHap.second         =                  fields[2]; //name
	    baseHapCounter.first   = destringify<char>(fields[1]);
	    baseHapCounter.second  =                  0;         //counter

	    for(unsigned int p=0;p<pos2hap[pos].size();p++){ //check if duplicated record
		if(pos2hap[pos][p].first == baseHap.first){
		    cerr<<"Error with line"<<linePos<<", duplicated allele in the  haplogroup diagnostic positions file"<<endl;
		    return 1;
		}
	    }
	    //cout<<pos<<endl;
	    pos2hap[pos].push_back(baseHap);
	    pos2hapCounter[pos].push_back(baseHapCounter);

	    
	}
	myFilePos.close();
    }else{
	cerr << "Unable to open file "<<haploGroupFile<<endl;
	return 1;
    }


    BamReader reader;

     if ( !reader.Open(bamfiletopen) ) {
	 cerr << "Could not open input BAM files : " << bamfiletopen <<endl;
	 return 1;
     }



     const SamHeader header      = reader.GetHeader();
     const RefVector references  = reader.GetReferenceData();

     myPileVisitor* cv = new myPileVisitor(references,
					   &pos2hap,
					   &pos2hapCounter,
					   &readname2hap,
					   &readname2posmask,
					   minQual);
     //coordToFind,&writerA,&writerC,&writerG,&writerT);
     PileupEngine pileup;
     pileup.AddVisitor(cv);


    BamAlignment al;
    //unsigned int numReads=0;
    while ( reader.GetNextAlignment(al) ) {
	// cout<<"1 "<<al.Name<<endl;
        pileup.AddAlignment(al);
    }


    //clean up
    pileup.Flush();
    reader.Close();
    // fastaReference.Close();
    // cout<<"test"<<endl;
    delete cv;

    BamReader reader2;


    map<string, vector<BamAlignment> > hap2al;

    if ( !reader2.Open(bamfiletopen) ) {
    	cerr << "Could not open input BAM files" << endl;
    	return 1;
    }



    //BamAlignment al;
    //read everything again
    while ( reader2.GetNextAlignment(al) ) {
	//cout<<"2 "<<al.Name<<endl;
	if(readname2hap.find(al.Name) != readname2hap.end()){//found
	    //decreasing the qualities of the base overlapping a diag. pos
	    for(unsigned int i=0;
		i<readname2posmask[al.Name].size();
		i++){
		al.Qualities[ readname2posmask[al.Name][i] ] = '#';//2 on the PHRED scale
	    }

	    hap2al[ readname2hap[al.Name] ] . push_back(al);
	}
    }
    reader2.Close();


    vector<string> allHaps = allKeysMap(hap2al);

    for(unsigned int i=0;
	i<allHaps.size();
	i++){
	BamWriter writer;
	cerr << "writing to   "<<bamFileOUTSuf+"_"+allHaps[i]+".bam"<<endl;
	
	if ( !writer.Open(bamFileOUTSuf+"_"+allHaps[i]+".bam",header,references) ) {
	    cerr << "Could not open output BAM file "<<bamFileOUTSuf+"_"+allHaps[i]+".bam"<<endl;
	    return 1;
	}
	
	for(unsigned int j=0;j<hap2al[ allHaps[i] ].size();j++){
	    writer.SaveAlignment( hap2al[ allHaps[i] ][j] );
	}
	    
	
	writer.Close();


    }



    for(unsigned int pos = 0;pos<pos2hapCounter.size();pos++){

	if(pos2hapCounter[pos].size() == 0 )
	    continue;
	cout<<pos;//
	int total=0;
	for(unsigned int i = 0;i<pos2hapCounter[pos].size();i++){
	    total+=pos2hapCounter[pos][i].second;
	}

	for(unsigned int i = 0;i<pos2hap[pos].size();i++){
	    cout<<"\t"<<pos2hap[pos][i].first<<"\t"<<pos2hap[pos][i].second<<"\t"<<pos2hapCounter[pos][i].second<<"\t"<<(double(pos2hapCounter[pos][i].second)/double(total));
	}
	cout<<endl;
    }


    cerr<<"Program finished successfully"<<endl;

    return 0;
}

