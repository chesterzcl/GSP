#ifndef VAR_H   
#define VAR_H
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<unordered_map>
#include<set>
#include<map>
#include "utilities.h"
using namespace std;

class var_list{
	public:
		vector<vector<string> > var_mat;
		vector<vector<string> > var_pop_mat;
		vector<string> var_gt_mat;

		var_list(){}
		var_list(const var_list &v1){
			var_mat=v1.var_mat;
			var_pop_mat=v1.var_pop_mat;
			var_gt_mat=v1.var_gt_mat;
		}

		void load_variant_data(string input_address){
			input.open(input_address);
			vector<string> temp_vec(2);
			check_file_open_status(input,input_address);
			cout<<"Loading variant list..."<<endl;
			string chr_str,pos_str;
			var_mat.clear();
			var_pop_mat.clear();
			var_gt_mat.clear();
			line_cter=1;
			while(getline(input,line)){
				if(line[0]!='#'){
					if(line_cter%100000==0){
						cout<<line_cter<<" variants scanned."<<endl;
					}
					chr_str=find_str_after_nth_char(line,0,'\t');
					pos_str=find_str_after_nth_char(line,1,'\t');
					temp_vec[0]=chr_str;
					temp_vec[1]=pos_str;
					var_mat.push_back(temp_vec);
					line_cter++;
				}
			}
			cout<<"Variant list loading complete. A total of "<<line_cter-1<<" variants loaded."<<endl<<endl;
			input.close();
		}

		void sort_var_pop_mat(){
			int cur_idx,temp_int;
			vector<string> cur_vec;
			for (int i = 1; i < var_pop_mat.size(); ++i){
				cur_vec=var_pop_mat[i];
				cur_idx=stoi(cur_vec.back());
				// cur_vec.pop_back();
				temp_int=i-1;
				while(cur_idx<stoi(var_pop_mat[temp_int].back())&&temp_int>=0){
					var_pop_mat[temp_int+1]=var_pop_mat[temp_int];
					--temp_int;
				}
				var_pop_mat[temp_int+1]=cur_vec;
			}
			for (int i = 0; i < var_pop_mat.size(); ++i){
				var_pop_mat[i].pop_back();
			}
		}	

		void update_var_list(){
			var_mat.clear();
			var_mat=var_pop_mat;
			var_gt_mat.clear();
			var_pop_mat.clear();
			cout<<"A total of "<<var_mat.size()<<" Candidate variants loaded."<<endl;
		}

		void output_gt_data(string output_address){
			output.open(output_address);
			check_file_open_status(output,output_address);
			for (int i = 0; i < var_gt_mat.size(); ++i){
				output<<var_gt_mat[i]<<endl;
			}
			output.close();
		}

		void merge_sort_variant_list(vector<string> input_vec,string opstring){
			string header_str,new_line;
			vector<string> var_vec;
			map<int,map<int,string> > chr_coord_dict;
			map<int,string > empty_coord_dict;
			map<string,int> header_col_dict;
			for (int i = 0; i < input_vec.size(); ++i){
				input.open(input_vec[i]);
				if(!input){
					cout<<"File "<<input_vec[i]<<" can not be opened and will be skipped."<<endl;
				}else{
					while(getline(input,line)){
						if(line[0]=='#'){
							header_col_dict=read_header_line(line);
							header_str="";
							var_vec=read_char_delim_str(line,'\t');
							for (int i = 0; i < 5; ++i){
								header_str+=var_vec[i]+'\t';
							}
							header_str+=var_vec[5];
							for (map<string,int>::iterator i = header_col_dict.begin(); i != header_col_dict.end(); ++i){
								header_str+='\t'+i->first;
							}
						}else{
							if(line.length()!=0){
								new_line="";
								var_vec=read_char_delim_str(line,'\t');
								if(!chr_coord_dict.count(stoi(var_vec[0].substr(7,2)))){
									chr_coord_dict[stoi(var_vec[0].substr(7,2))]=empty_coord_dict;
								}
								if(!chr_coord_dict[stoi(var_vec[0].substr(7,2))].count(stoi(var_vec[1]))){
									for (int i = 0; i < 5; ++i){
										new_line+=var_vec[i]+'\t';
									}
									new_line+=var_vec[5];
									for (map<string,int>::iterator i = header_col_dict.begin(); i != header_col_dict.end(); ++i){
										new_line+='\t'+var_vec[i->second];
									}
									chr_coord_dict[stoi(var_vec[0].substr(7,2))][stoi(var_vec[1])]=new_line;
								}									
							}
						}
					}
				}
				input.close();
			}
			output.open(opstring);
			output<<header_str<<endl;
			for (map<int,map<int,string> >::iterator i = chr_coord_dict.begin(); i != chr_coord_dict.end(); ++i){
				for (map<int,string>::iterator j = i->second.begin(); j != i->second.end(); ++j){
					output<<j->second<<endl;
				}
			}
			output.close();			
		}

		map<string,int> read_header_line(string header){
			vector<string> header_vec;;
			map<string,int> header_col_dict;
			header_vec=read_char_delim_str(header,'\t');
			for (int i = 6; i < header_vec.size(); ++i){
				header_col_dict[header_vec[i]]=i;
			}
			return header_col_dict;
		}

	private:
		ifstream input;
		ofstream output;
		int line_cter;
		string line;
		vector<string> line_vec;
};

#endif