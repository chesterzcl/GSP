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

int main(int argc, char const *argv[])
{
	ifstream input;
	ofstream output;
	vector<string> line_vec,header_vec;
	string line;
	map<string,int> pop_col_dict;
	int min_eff_sample=300;
	string dir=argv[1];
	string input_file_name=argv[2];
	string filtered_file_name=input_file_name+"_filtered_"+to_string(min_eff_sample);
	input.open(dir+"/"+input_file_name+".txt");
	output.open(dir+"/"+filtered_file_name+".txt");
	output<<"#Chromosome\tPosition\tGenomic region\tVariant type\tAnnotation\tSignature carrying population"<<endl;
	int cter=0;
	while(getline(input,line)){
		cter++;
		if(cter%100000==0){
			cout<<cter<<" line analyzed."<<endl;
		}
		if(line[0]!='#'){
			line_vec=read_char_delim_str(line,'\t');
			int eff_sample=0;
			int total_control=0,pos_control=0;
			string tar_pop=find_str_after_nth_char(line_vec[5],3,'|');
			int tar_col=pop_col_dict[tar_pop];
			string op_str="";
			int breed_cnt=0;
			for(int i=6;i!=line_vec.size();++i){
				int denom=stoi(find_str_after_nth_char(line_vec[i],1,'/'));
				if(denom>=3){
					eff_sample+=denom;
					int num=stoi(find_str_after_nth_char(line_vec[i],0,'/'));
					double freq=(double) num/denom;
					if(i==tar_col){
						breed_cnt+=1;
						op_str+='\t'+header_vec[i]+':'+to_string(num)+'/'+to_string(denom);
					}else{
						total_control+=denom;
						pos_control+=num;
					}
				}
			}

			if(total_control>=min_eff_sample){
				output<<stoi(line_vec[0].substr(7,2))-4<<'\t'<<line_vec[1]<<'\t'<<line_vec[2]<<'\t'<<line_vec[3]<<'\t'<<line_vec[4]<<'\t'<<line_vec[5]<<op_str/*<<'\t'<<to_string(pos_control)+'/'+to_string(total_control)*/<<endl;		
			}
		}else{
			header_vec=read_char_delim_str(line,'\t');
			for (int i = 6; i < header_vec.size(); ++i){
				pop_col_dict[header_vec[i]]=i;
			}
		}
	}
	return 0;
}














