#ifndef ANN_H   
#define ANN_H
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<unordered_map>
#include<set>
#include<map>
#include "utilities.h"
using namespace std;

class ann_data{
	public:
		unordered_map<string,set<string> > gene_transcript_dict;

		ann_data(){}
		ann_data(const ann_data& a1){
			gene_transcript_dict=a1.gene_transcript_dict;
		}

		void load_annotation_db(string input_address){
			input.open(input_address);
			check_file_open_status(input,input_address);
			string gene;
			int cter=1;
			while(getline(input,line)){
				if(line[0]!='#'){
					if(cter%100000==0){
						cout<<cter<<" lines read."<<endl;
					}
					cter++;
					line_vec=read_char_delim_str(line,'\t');
					gene=line_vec[6];
					if(!gene_transcript_dict.count(gene)){
						gene_transcript_dict[gene]=empty_set;
					}
					gene_transcript_dict[gene].insert(line_vec[10]);
				}
			}
			input.close();
			cout<<"Annotation database loaded."<<endl<<endl;
		}
	private:
		ifstream input;
		set<string> empty_set;
		string line;
		vector<string> line_vec;
};

#endif