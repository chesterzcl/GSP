#ifndef POP_H   
#define POP_H
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<unordered_map>
#include<set>
#include "utilities.h"
using namespace std;

class pop_data{
	public:
		unordered_map<string,set<string> > pop_dict;
		unordered_map<string,string> sample_dict;
		unordered_map<string,set<int> > pop_col_dict;
		unordered_map<string,int> sample_col_dict;
		unordered_map<int,string> col_sample_dict;
		unordered_map<int,double> col_depth_dict;

		pop_data(){}
		pop_data(const pop_data &p1){
			pop_dict=p1.pop_dict;
			sample_dict=p1.sample_dict;
			pop_col_dict=p1.pop_col_dict;
			sample_col_dict=p1.sample_col_dict;
			col_sample_dict=p1.col_sample_dict;
			col_depth_dict=p1.col_depth_dict;
		}

		void load_sample_depth_sumstats(string input_address){
			input.open(input_address);
			check_file_open_status(input,input_address);
			cout<<"Loading sample depth data..."<<endl;	
			col_depth_dict.clear();
			while(getline(input,line)){
				line_vec=read_char_delim_str(line,'\t');
				if(!sample_col_dict.count(line_vec[0])){
					cout<<"Can't find sample: "<<line_vec[0]<<" in the dataset."<<endl;
					cout<<"Shut down the engine."<<endl;
					exit(0);
				}
				col_depth_dict[sample_col_dict[line_vec[0]]]=stod(line_vec[1])-2*stod(line_vec[2]);
			}
			cout<<"Sample read depth load complete."<<endl<<endl;

		}

		void load_population_data(string input_address){
			input.open(input_address);
			check_file_open_status(input,input_address);
			cout<<"Loading population data..."<<endl;
			pop_dict.clear();
			sample_dict.clear();
			pop_col_dict.clear();
			sample_col_dict.clear();
			col_sample_dict.clear();
			while(getline(input,line)){
				line_vec=read_char_delim_str(line,'\t');
				if(!pop_dict.count(line_vec[1])){
					pop_dict[line_vec[1]]=empty_str_set;
				}
				pop_dict[line_vec[1]].insert(line_vec[0]);
				sample_dict[line_vec[0]]=line_vec[1];
			}
			cout<<"Population data loaded for "<<sample_dict.size()<<" samples from "<<pop_dict.size()<<" populations."<<endl<<endl;
			input.close();
		}

		void read_vcf_address(string input_address){
			vcf_add=input_address;
		}

		void index_data(string input_address){
			input.open(input_address);
			check_file_open_status(input,input_address);
			cout<<"Start indexing population info..."<<endl;
			int cter=9;
			sample_col_dict.clear();
			pop_col_dict.clear();
			col_sample_dict.clear();
			while(getline(input,line)){
				line_vec=read_char_delim_str(line,'\t');
				sample_col_dict[line_vec[0]]=cter;
				col_sample_dict[cter]=line_vec[0];
				if(!pop_col_dict.count(line_vec[1])){
					pop_col_dict[line_vec[1]]=empty_int_set;
				}
				pop_col_dict[line_vec[1]].insert(cter);
				cter++;
			}
			input.close();
			cout<<"Indexing complete."<<endl<<endl;
		}

		void index_data(){
			input.open(vcf_add);
			check_file_open_status(input,vcf_add);
			cout<<"Start indexing population info..."<<endl;
			string temp_str;
			sample_col_dict.clear();
			pop_col_dict.clear();
			col_sample_dict.clear();
			while(getline(input,line)){
				if(line[0]!='#'){
					break;
				}
				temp_str=line;
			}
			input.close();
			line_vec=read_char_delim_str(temp_str,'\t');
			for (int i = 9; i < line_vec.size(); ++i){
				sample_col_dict[line_vec[i]]=i;
				col_sample_dict[i]=line_vec[i];
				if(!pop_col_dict.count(line_vec[i])){
					pop_col_dict[line_vec[i]]=empty_int_set;
				}
				pop_col_dict[line_vec[1]].insert(i);				
			}
			cout<<"Indexing complete."<<endl<<endl;
		}

	private:
		ifstream input;
		set<string> empty_str_set;
		set<int> empty_int_set;		
		string line,vcf_add;
		vector<string> line_vec;
};

#endif