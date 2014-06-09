/*
 * damage2profile
 * Date: Jun-08-2014 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <gzstream.h>
#include <algorithm> 

#include "utils.h"

using namespace std;

int main (int argc, char *argv[]) {


    bool singleS=false;
    bool doubleS=false;
    int length = 5;

    const string usage=string("\t"+string(argv[0])+
			      " [options]  [file produced by damage-patterns] "+"\n\n"+

			      "\t\t"+"-length [length]" +"\t"  +"Do not consider bases beyond this length  (default "+stringify(length)+") "+"\n"+

			      "\t\t"+"-single" +"\t\t\t"+"Use a single strand profile (default: "+booleanAsString(singleS)+")"+"\n"+
			      "\t\t"+"-double" +"\t\t\t"+"Use a double strand profile (default: "+booleanAsString(doubleS)+")"+"\n"+

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


    for(int i=1;i<(argc-1);i++){ 

	if(string(argv[i]) == "-length"  ){
	    length=destringify<int>(argv[i+1]);
	    i++;
	    continue;
	}

	if(string(argv[i]) == "-single" ){
	    singleS=true;
	    continue;
	}

	if(string(argv[i]) == "-double" ){
	    doubleS=true;
	    continue;
	}

	cerr<<"Wrong option "<<string(argv[i])<<endl;
	return 1;

    }



    if(doubleS && singleS){
	cerr<<"Cannot specify both single and double stranded protocols"<<endl;
	return 1;
    }
	       
   string line;
   igzstream myFile;
   string filename = string(argv[argc-1]);

   myFile.open(filename.c_str(), ios::in);
   bool firstLine=true;
   bool is5p=true;
   string dnaAlphabet="ACGT";
   vector<string> header;
   for(int nuc1=0;nuc1<4;nuc1++){
       for(int nuc2=0;nuc2<4;nuc2++){
	   if(nuc1==nuc2){ // prob of error is 0 if both nucleotides are identical
	   }else{ //           rely on the substitution frequency
	       //cout<<dnaAlphabet[nuc1]<<endl;
	       string s1= stringify(dnaAlphabet[nuc1]);
	       string s2= stringify(dnaAlphabet[nuc2]);

	       string toadd = s1+">"+s2;
	       header.push_back(toadd);
	   }
       }

   }
   cout<<vectorToString(header,"\t")<<endl;


   vector<string> stacksStrings;
   if (myFile.good()){
       getline (myFile,line);//header
       while ( getline (myFile,line)){
	   //cout<<line<<endl;
	   vector<string> fields = allTokens(line,'\t');
	   if(firstLine){
	       is5p=(fields[0]=="0");
	       firstLine=false;
	   }

	   int pos = destringify<int>(fields[0]);
	   
	   if(is5p){
	       if(pos>=length)
		   break;	      
	       vector<double> toprint;
	       for(unsigned int i=1;i<fields.size();i++){
		   
		   if(i==6){ //C->T  regardless of protocol
		       vector<string> fields2 = allTokens(fields[i],' ');
		       toprint.push_back(destringify<double>(fields2[0]));
		   }else{
		       toprint.push_back(0.0);
		   }		   
	       }
	       //cout<<vectorToString(toprint,"\t")<<endl;
	       stacksStrings.push_back( vectorToString(toprint,"\t") );
	   }else{
	       //	       cout<<(-1.0*pos)<<endl;
	       if( (-1.0*pos)>=length)
		   continue;
	       

	       vector<double> toprint;
	       for(unsigned int i=1;i<fields.size();i++){
		   
		   if(singleS){ //single strand
		       if(i==6){ //C->T 
			   vector<string> fields2 = allTokens(fields[i],' ');
			   toprint.push_back(destringify<double>(fields2[0]));
		       }else{
			   toprint.push_back(0.0);
		       }
		   }else{ //double
		       if(i==7){ //G->A
			   vector<string> fields2 = allTokens(fields[i],' ');
			   toprint.push_back(destringify<double>(fields2[0]));
		       }else{
			   toprint.push_back(0.0);
		       }

		   }	   
	       }
	       //cout<<vectorToString(toprint,"\t")<<endl;
	       stacksStrings.push_back( vectorToString(toprint,"\t") );

	   }
       }

       myFile.close();
   }else{
       cerr << "Unable to open file "<<filename<<endl;
       return 1;
   }

   if(!is5p)
       reverse(stacksStrings.begin(),stacksStrings.end());


   cout<<vectorToString(stacksStrings,"\n")<<endl;
   return 0;
}

