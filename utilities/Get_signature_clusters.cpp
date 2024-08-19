#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<set>
#include<unordered_map>
#include<map>
using namespace std;

vector<string> read_char_delim_str(string& line,char delim){
	vector<string> output;
	string temp="";
	for (int i = 0; i < line.size(); ++i){
		if (line[i]==delim){
			output.push_back(temp);
			temp="";
		}else{
			temp+=line[i];
		}
	}
	output.push_back(temp);
	return output;
}

string find_str_after_nth_char(string& str,int n,char delim){
	int cter=0,idx;
	string op_str;
	if(n>0){
		for (int i = 0; i < str.length(); ++i){
			if(str[i]==delim){
				cter++;
			}
			if(cter==n){
				idx=i+1;
				while(str[idx]!=delim&&idx<str.length()){
					op_str+=str[idx];
					idx++;
				}
				break;
			}
		}		
	}else{
		for (int i = 0; i < str.length(); ++i){
			if(str[i]==delim){
				break;
			}
			op_str+=str[i];
		}
	}
	return op_str;
}

double str_to_freq(string& str){
	int nom=stoi(find_str_after_nth_char(str,0,'/'));
	int denom=stoi(find_str_after_nth_char(str,1,'/'));
	if(denom<3){
		return -1;
	}else{
		return (double)nom/(double)denom;
	}
}

vector<string> combine_two_queries(vector<string> v1,vector<string> v2){
	vector<string> temp_vec;
	for(int i=0;i<6;++i){
		temp_vec.push_back(v1[i]);
	}
	for(int i=6;i<v1.size();++i){
		temp_vec.push_back(to_string(stoi(v1[i])+stoi(v2[i])));
	}
	return temp_vec;
}


int main(int argc, char const *argv[]){
	//file param1 param2
	string ip_file=argv[1];
	string param_ref_num=argv[2];
	string param_min_sig_num=argv[3];
	string param_max_dist=argv[4];
	int min_eff_sample=stoi(param_ref_num);
	ifstream input;
	ofstream output;
	input.open(ip_file);
	//container
	vector<string> line_vec;
	string line;
	map<string,vector<string>> group_sig_dict;
	map<string,int> group_idx_dict;
	map<int,string> idx_group_dict;
	int cter=0;
	//Sort by sample group
	while(getline(input,line)){
		if(line[0]!='#'){
			line_vec=read_char_delim_str(line,'\t');
			vector<string> ann_vec=read_char_delim_str(line_vec[5],'|');
			string tar_group=ann_vec.back();
			set<int> tar_set;
			tar_set.insert(group_idx_dict[tar_group]);
			//check validity
			int valid_num_cnt=0;
			for (int i = 6; i < line_vec.size()&&(!tar_set.count(i)); ++i)
			{
				int denom=stoi(find_str_after_nth_char(line_vec[i],1,'/'));
				valid_num_cnt+=denom;
			}

			if(valid_num_cnt>=min_eff_sample){
				group_sig_dict[tar_group].push_back(line);
			}
		}else{
			line_vec=read_char_delim_str(line,'\t');
			for (int i = 6; i < line_vec.size(); ++i)
			{
				group_idx_dict[line_vec[i]]=i;
				idx_group_dict[i]=line_vec[i];
			}
		}
	}

	//Aggregate for enrichment info
	//param
	int min_sig_num=stoi(param_min_sig_num);
	int max_dist=stoi(param_max_dist);
	for (map<string,vector<string>>::iterator i = group_sig_dict.begin(); i != group_sig_dict.end(); ++i)
	{
		string tar_group=i->first;
		vector<string> segment_vec;
		string prev_chr="";
		int prev_pos=-1;
		for (int j = 0; j != i->second.size(); ++j)
		{
			string cur_chr=find_str_after_nth_char(i->second[j],0,'\t');
			int cur_pos=stoi(find_str_after_nth_char(i->second[j],1,'\t'));

			if((cur_chr!=prev_chr)||(cur_pos-prev_pos)>max_dist){
				if(segment_vec.size()>=min_sig_num){
					string seg_chr=find_str_after_nth_char(segment_vec[0],0,'\t');
					string seg_pos_start=find_str_after_nth_char(segment_vec[0],1,'\t');
					string seg_pos_end=find_str_after_nth_char(segment_vec[segment_vec.size()-1],1,'\t');

					output.open(tar_group+"_"+seg_chr+"_"+seg_pos_start+"_"+seg_pos_end+".txt");
					for (int k = 0; k != segment_vec.size(); ++k)
					{
						output<<segment_vec[k]<<endl;
					}
					output.close();

					//reset
					segment_vec.clear();
				}
			}

			segment_vec.push_back(i->second[j]);
			prev_chr=cur_chr;
			prev_pos=cur_pos;

		}
		if(segment_vec.size()>=min_sig_num){
			string seg_chr=find_str_after_nth_char(segment_vec[0],0,'\t');
			string seg_pos_start=find_str_after_nth_char(segment_vec[0],1,'\t');
			string seg_pos_end=find_str_after_nth_char(segment_vec[segment_vec.size()-1],1,'\t');

			output.open(tar_group+"_"+seg_chr+"_"+seg_pos_start+"_"+seg_pos_end+".txt");
			for (int k = 0; k != segment_vec.size(); ++k)
			{
				output<<segment_vec[k]<<endl;
			}
			output.close();

			//reset
			segment_vec.clear();
		}

	}
	return 0;
}

















