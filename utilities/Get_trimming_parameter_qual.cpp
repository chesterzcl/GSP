#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<set>
#include<cmath>
#include<cstdlib>
#include<unordered_map>
#include<map>
#include<sstream>
#include<iomanip>
#include<unordered_set>
#include<algorithm>

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
	string input_file=argv[1];
	string line;
	vector<string> line_vec;
	ifstream input;
	ofstream output;
	bool found=false;
	double total_qual_score=0;
	double total_base_num=0;
	string prob="Per base sequence quality";
	map<int,double> qual_dict;
	input.open(input_file);
	while(getline(input,line)){
			//Start scanning
			if(found){
				if(line[0]!='#'&&line[0]!='>'){
					line_vec=read_char_delim_str(line,'\t');
					vector<string> pos_vec=read_char_delim_str(line_vec[0],'-');
					if(pos_vec.size()==1){
						qual_dict[stoi(pos_vec[0])]=stod(line_vec[5]);
					}else{
						for (int i = stoi(pos_vec[0]); i <= stoi(pos_vec[1]); ++i){
							qual_dict[i]=stod(line_vec[5]);
						}
					}
				}
			}
			//Scan for begining
			if(line.substr(0,2)==">>"){
				string line_suff=line.substr(2);
				line_vec=read_char_delim_str(line_suff,'\t');
				//Exit
				if(line_vec[0]=="END_MODULE"&&found){
					break;
				}
				//Mark the begining
				if(line_vec[0]==prob){
					found=true;
				}
			}
	}
	int start=10;
	int end=qual_dict.size()-10;
	vector<int> qual_vec;
	for (map<int,double>::iterator i = qual_dict.begin(); i !=qual_dict.end() ; ++i)
	{
		if(start<=i->first&&i->first<=end){
			qual_vec.push_back(i->second);
		}	
	}
	sort(qual_vec.begin(),qual_vec.end());
	cout<<qual_vec[qual_vec.size()/4]<<endl;
	return 0;
}


