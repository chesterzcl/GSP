#pragma once
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<set>
#include<unordered_map>
using namespace std;

void check_file_open_status(ifstream& input,string input_address){
	if(!input){
		cout<<"Failed to open the input file: "<< input_address<<endl<<endl;
		abort();
	}
	cout<<"File: "<<input_address<<" successfully opened."<<endl;
}

void check_file_open_status(ofstream& output,string output_address){
	if(!output){
		cout<<"Failed to open the output file: "<< output_address<<endl<<endl;
		abort();
	}
	cout<<"File: "<<output_address<<" successfully opened."<<endl;
}

vector<string> read_char_delim_str(string line,char delim){
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

int check_genotype(string geno_field){
	if(geno_field[0]=='1'){
		return 1;
	}else{
		return 0;
	}
}

string generate_line_header(vector<string> line_vec){
	vector<string> ann_vec=read_char_delim_str(line_vec[7],'|');
	string lheader_str=line_vec[0]+'\t'+line_vec[1]+'\t'+ann_vec[3]+'\t'+ann_vec[1]+'\t'+ann_vec[13]+'\t'+ann_vec[10];
	return lheader_str;
}

string find_str_after_nth_char(string str,int n,char delim){
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

class data_printer{
	public:
		void print(vector<string> a){
			for (int i = 0; i < a.size(); ++i){
				cout<<a[i]<<'\t';
			}
			cout<<endl;
		}

		void print(vector<vector<string> > a){
			for (int i = 0; i < a.size(); ++i){
				print(a[i]);
			}
		}

		void print(set<int> a){
			for (set<int>::iterator i = a.begin(); i != a.end(); ++i){
				cout<<*i<<endl;
			}
		}
		void print(set<string> a){
			for (set<string>::iterator i = a.begin(); i != a.end(); ++i){
				cout<<*i<<endl;
			}
		}
		void print(unordered_map<string, string> a){
			for (unordered_map<string, string>::iterator i = a.begin(); i != a.end(); ++i){
				cout<<i->first<<'\t'<<i->second<<endl;
			}
		}
		void print(unordered_map<string, int> a){
			for (unordered_map<string, int>::iterator i = a.begin(); i != a.end(); ++i){
				cout<<i->first<<'\t'<<i->second<<endl;
			}
		}
		void print(unordered_map<string, set<string> > a){
			for (unordered_map<string, set<string> >::iterator i = a.begin(); i != a.end(); ++i){
				cout<<i->first<<": "<<endl;
				print(i->second);
			}
		}
		void print(unordered_map<string, set<int> > a){
			for (unordered_map<string, set<int> >::iterator i = a.begin(); i != a.end(); ++i){
				cout<<i->first<<": "<<endl;
				print(i->second);
			}
		}
		void print(unordered_map<string, double> a){
			for (unordered_map<string, double>::iterator i = a.begin(); i != a.end(); ++i){
				cout<<i->first<<'\t'<<i->second<<endl;
			}
		}
		void print(unordered_map<int, int> a){
			for (unordered_map<int, int>::iterator i = a.begin(); i != a.end(); ++i){
				cout<<i->first<<'\t'<<i->second<<endl;
			}
		}
		void print(unordered_map<int, string> a){
			for (unordered_map<int, string>::iterator i = a.begin(); i != a.end(); ++i){
				cout<<i->first<<'\t'<<i->second<<endl;
			}
		}
		void print(unordered_map<int, double> a){
			for (unordered_map<int, double>::iterator i = a.begin(); i != a.end(); ++i){
				cout<<i->first<<'\t'<<i->second<<endl;
			}
		}
};



